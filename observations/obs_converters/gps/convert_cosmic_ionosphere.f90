! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

!> An observation sequence converter for ionosphere profiles 
!> from the CDAAC web site http://cosmic-io.cosmic.ucar.edu/cdaac
!> The usual mission is 'cosmic2013', the file type is 'ionPrf'.
!>
!> The strategy here is to read the raw vertical locations and interpolate
!> to a fixed set of vertical locations. If there are any 'bad' observations
!> below the lowest desired level - they are ignored. If there are any 'bad'
!> observations ABOVE the lowest desired level - the entire profile is IGNORED.
!>
!> modified from GPS atmospheric radio/occultation version
!> which was created June 2008 by Ryan Torn, NCAR/MMM
!> nsc 11 mar 2016  
!>
!> Incorporated multiple methods for setting the observation error variance
!>   Feb. 2010  I-TE LEE, NCAR/HAO & NCU
!>             Modify to read the ionPrf data file of COSMIC
!>             and write the data to the DART format.  
!>   Oct. 2010  I-TE LEE, NCAR/HAO & NCU
!>             Remove the code for lower atmospheric observations
!>   Dec. 2010  I-TE LEE, NCAR/HAO & NCU
!>             Add new error valur for log scale testing
!>   Jan. 2011  I-TE LEE, NCAR/HAO & NCU
!>             Modify for real GOX observations
!>   Jan. 2011  I-TE LEE, NACR/HAO & NCU
!>             Add subroutine to calculate the observation variance 
!> the ancestor of this routine modified the observation time with the desired
!> analysis time - this is not done - windowing is done in a post-processing step.
!>
!> The normal usage pattern is to create a list of input files for all satellite
!> pairs for a given day and create a SINGLE observation sequence file. That file
!> can then be post-processed any way you like - but it is a lot more manageable
!> than reprocessing a thousand (literally) netCDF files.

program convert_cosmic_ionosphere

use          types_mod, only : r8, MISSING_R8
use   time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                               increment_time, get_time, set_date, operator(-),  &
                               print_date
use      utilities_mod, only : initialize_utilities, find_namelist_in_file,      &
                               check_namelist_read, nmlfileunit, do_nml_file,    &
                               get_next_filename, error_handler, E_ERR, E_MSG,   &
                               find_textfile_dims, do_nml_term,                  &
                               to_upper, open_file, finalize_utilities
use       location_mod, only : VERTISHEIGHT, set_location
use   obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,        &
                               static_init_obs_sequence, init_obs, destroy_obs,  &
                               write_obs_seq, init_obs_sequence, get_num_obs,    &
                               insert_obs_in_seq, destroy_obs_sequence,          &
                               set_copy_meta_data, set_qc_meta_data, set_qc,     & 
                               set_obs_values, set_obs_def, insert_obs_in_seq,   &
                               get_num_copies, get_num_qc
use        obs_def_mod, only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                               set_obs_def_error_variance, set_obs_def_location, &
                               set_obs_def_key
use  obs_utilities_mod, only : add_obs_to_seq
use       obs_kind_mod, only : COSMIC_ELECTRON_DENSITY

use netcdf_utilities_mod, only : nc_check, nc_get_variable

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"
character(len=*), parameter :: routine  = 'convert_cosmic_ionosphere'

integer, parameter :: METHOD_CONSTANT      = 1   !! = 'constant'
integer, parameter :: METHOD_SCALED        = 2   !! = 'scaled'
integer, parameter :: METHOD_LOOKUP        = 3   !! = 'lookup'
integer, parameter :: METHOD_SCALED_LOOKUP = 4   !! = 'scaled_lookup'
integer            :: method

character (len=512) :: string1, string2, string3
character (len=256) :: next_infile

integer :: ncid, varid, nlevels, k, nfiles, num_new_obs, oday, osec, &
           iyear, imonth, iday, ihour, imin, isec, zloc, obs_num, &
           io, iunit, nobs, filenum, dummy, numrejected
logical :: file_exist, first_obs, from_list = .false.
real(r8) :: hght_miss, elecd_miss, oerr,   & 
            qc, lato, lono, hghto, wght, & 
            obsval, obs_val(1), qc_val(1)

real(r8), allocatable :: lat(:), lon(:), hght(:), elecd(:)

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

integer :: num_copies = 1  ! number of copies in sequence
integer :: num_qc     = 1  ! number of QC entries
integer :: existing_num_copies  ! number of copies in existing sequence
integer :: existing_num_qc      ! number of QC entries in existing sequence

character(len=512) :: hght_units  ! should be long enough

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

integer, parameter :: NMAXLEVELS = 200   !  max number of desired observation levels

character(len=256) :: input_file               = ''
character(len=256) :: input_file_list          = 'input_file_list.txt'
character(len=256) :: output_file              = 'obs_seq.out'
character(len=256) :: observation_error_file   = 'none'
character(len=128) :: observation_error_method = 'scaled_lookup'
logical            :: locations_only           = .false.
real(r8) :: obs_error_factor       =  1.0_r8
real(r8) :: obs_levels(NMAXLEVELS) = -1.0_r8   ! in kilometers
integer  :: verbose                = 0   ! higher is more output

namelist /convert_cosmic_ionosphere_nml/ &
             input_file,                 &
             input_file_list,            &
             output_file,                &
             observation_error_file,     &
             observation_error_method,   &
             obs_error_factor,           &
             locations_only,             &
             obs_levels,                 &
             verbose

! initialize some values
obs_num   = 1
qc        = 0.0_r8
first_obs = .true.
call set_calendar_type(GREGORIAN)

!  read the necessary parameters from input.nml
call initialize_utilities()

call find_namelist_in_file("input.nml", "convert_cosmic_ionosphere_nml", iunit)
read(iunit, nml = convert_cosmic_ionosphere_nml, iostat = io)
call check_namelist_read(iunit, io, "convert_cosmic_ionosphere_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=convert_cosmic_ionosphere_nml)
if (do_nml_term()) write(     *     , nml=convert_cosmic_ionosphere_nml)

! This is basically the OSSE case - we only want locations
if (locations_only) then
   num_copies = 0
   num_qc     = 0
endif

! it is faster to simply check an integer rather than a character string
! each observation will need to check the method
method = check_error_method()

! if more than MAXLEVELS levels are specified in the namelist,
! a namelist read error is generated.
nlevels = 0
LEVELLOOP: do k = 1, NMAXLEVELS
  if ( obs_levels(k) == -1.0_r8 )  exit LEVELLOOP
  nlevels = k
enddo LEVELLOOP
do k = 2, nlevels
  if ( obs_levels(k-1) >= obs_levels(k) ) then
    call error_handler(E_ERR, routine, 'Observation levels must increase', &
                       source, revision, revdate)
  endif
enddo

! cannot have both a single filename and a list;
if (input_file /= '' .and. input_file_list /= '') then
  call error_handler(E_ERR, routine,  &
                     'One of input_file or input_file_list must be NULL', &
                     source, revision, revdate)
endif
if (input_file_list /= '') from_list = .true.

! need to know a reasonable max number of obs that could be added here.
if (from_list) then
   call find_textfile_dims(input_file_list, nfiles, dummy)
   num_new_obs = nlevels * nfiles
else
   num_new_obs = nlevels
endif

! either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=output_file, exist=file_exist)
if ( file_exist ) then

   write(string1,*) "found existing obs_seq file, appending to ", trim(output_file)
   write(string2,*) "adding up to a maximum of ", num_new_obs, " new observations"
   call error_handler(E_MSG, routine, string1, &
                      source, revision, revdate, text2=string2)
   call read_obs_seq(output_file, 0, 0, num_new_obs, obs_seq)

   ! check to see if existing file is compatible with locations_only setting
   existing_num_copies = get_num_copies(obs_seq)
   existing_num_qc     = get_num_qc(obs_seq)

   if (existing_num_copies /= num_copies .or.  existing_num_qc /= num_qc) then
      write(string1,*)'incompatible existing observation sequence file'
      write(string2,'(A,i4,A,i4)')'expected ',num_copies, &
                           ' obs copies got ',existing_num_copies
      write(string3,'(A,i4,A,i4)')'expected ',num_qc, &
                           ' QC  copies got ',existing_num_qc
      call error_handler(E_ERR, routine, string1, &
                  source, revision, revdate, text2=string2, text3=string3)
   endif

else

   write(string1,*) "no existing obs_seq file, creating ", trim(output_file)
   write(string2,*) "with up to a maximum of ", num_new_obs, " observations"
   call error_handler(E_MSG, routine, string1, &
                      source, revision, revdate, text2=string2)

   call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
   if (.not. locations_only) then
      call set_copy_meta_data(obs_seq, 1, 'observation')
      call set_qc_meta_data(obs_seq, 1, 'COSMIC QC')
   endif

endif

! main loop that does either a single file or a list of files.
! the data is currently distributed as a single profile per file
! with hundreds (even thousands) of profiles per day.

filenum = 1
numrejected = 0

fileloop: do      ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then 
      next_infile = get_next_filename(input_file_list, filenum)
   else
      next_infile = input_file
      if (filenum > 1) next_infile = ''
   endif
   if (next_infile == '') exit fileloop

   if ( verbose > 0 ) write(*,*)' '  ! improves readability 
  
   !  open the next occultation profile
   call nc_check( nf90_open(next_infile, nf90_nowrite, ncid), 'open ',next_infile)

   ! process the profile
   call nc_check( nf90_get_att(ncid,nf90_global,'year',  iyear) ,'get_att year',   next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'month', imonth),'get_att month',  next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'day',   iday)  ,'get_att day',    next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'hour',  ihour) ,'get_att hour',   next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'minute',imin)  ,'get_att minute', next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'second',isec)  ,'get_att second', next_infile)
   
   time_obs = set_date(iyear, imonth, iday, ihour, imin, isec)
   call get_time(time_obs,  osec, oday)
   
   call nc_check( nf90_inq_dimid(ncid, "MSL_alt", varid),        'inq dimid MSL_alt', next_infile)
   call nc_check( nf90_inquire_dimension(ncid, varid, len=nobs), 'inq dim   MSL_alt', next_infile)
   
   allocate(  lat(nobs))
   allocate(  lon(nobs))
   allocate( hght(nobs))
   allocate(elecd(nobs))
   
   ! read the latitude array
   call nc_check( nf90_inq_varid(ncid, "GEO_lat", varid) ,'inq varid GEO_lat', next_infile)
   call nc_check( nf90_get_var(ncid, varid, lat)         ,'get var   GEO_lat', next_infile)
   
   ! read the latitude array
   call nc_check( nf90_inq_varid(ncid, "GEO_lon", varid) ,'inq varid GEO_lon', next_infile)
   call nc_check( nf90_get_var(ncid, varid, lon)         ,'get var   GEO_lon', next_infile)
   
   ! read the altitude array  - could check that units are 'km'
   call nc_check( nf90_inq_varid(ncid, "MSL_alt", varid) ,'inq varid MSL_alt', next_infile)
   call nc_check( nf90_get_var(ncid, varid, hght)        ,'get_var   MSL_alt', next_infile)
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', hght_miss) ,'get_att _FillValue MSL_alt', next_infile)
   call nc_check( nf90_get_att(ncid, varid, 'units', hght_units) ,'get_att units MSL_alt', next_infile)
   
   ! read the electron density
   call nc_check( nf90_inq_varid(ncid, "ELEC_dens", varid) ,'inq varid ELEC_dens', next_infile)
   call nc_check( nf90_get_var(ncid, varid, elecd)         ,'get var   ELEC_dens', next_infile)
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', elecd_miss) ,'get_att _FillValue ELEC_dens', next_infile)
   
   call nc_check( nf90_close(ncid) , 'close file', next_infile)

   ! may want to figure out where the problem observations come from
   if ( verbose > 1 ) then
      call summarize( hght, 'MSL_alt',   next_infile)
      call summarize(elecd, 'ELEC_dens', next_infile)
   endif

   ! There are several ways that bad values are present.
   ! If the electron density is flagged bad, 
   ! then interp_height_wght() will fail and we can cycle.

   where (  hght ==  hght_miss   ) elecd = MISSING_R8
   where ( elecd == elecd_miss   ) elecd = MISSING_R8
   where ( elecd <  0.0_r8       ) elecd = MISSING_R8
   where ( elecd >  1000000.0_r8 ) elecd = MISSING_R8

   ! only need the densities above the minimum desired level
   ! discard entire profile if there are negative densities above lowest level
    
   if ( bad_profile_check(obs_levels(1), hght, elecd, next_infile) ) then
      deallocate( lat, lon, hght, elecd )
      filenum     = filenum     + 1
      numrejected = numrejected + 1
      cycle fileloop
   endif
   
   LEVELS: do k = 1, nlevels
   
     call interp_height_wght(hght, obs_levels(k), nobs, zloc, wght)

     if ( zloc < 1 ) cycle LEVELS
   
     ! lon(zloc) and lon(zloc+1) range from -180 to +180
     ! call a subroutine to handle the wrap point, and convert to [0,360].
     lono   = compute_lon_wrap(lon(zloc), lon(zloc+1), wght)
     lato   = wght *   lat(zloc) + (1.0_r8 - wght) *   lat(zloc+1)
     hghto  = wght *  hght(zloc) + (1.0_r8 - wght) *  hght(zloc+1)
     obsval = wght * elecd(zloc) + (1.0_r8 - wght) * elecd(zloc+1)
 
     ! if you need a different observation error variance, modify it by extending
     ! the electron_density_error() function. Please do not just hardcode something here.
     oerr  = electron_density_error(lono, lato, hghto, ihour, imin, method, &
                                     obs_error_factor, obsval)
   
     call set_obs_def_location(obs_def,set_location(lono,lato,hghto*1000.0_r8,VERTISHEIGHT))
     call set_obs_def_type_of_obs(obs_def, COSMIC_ELECTRON_DENSITY)
     call set_obs_def_time(obs_def, set_time(osec, oday))
     call set_obs_def_error_variance(obs_def, oerr * oerr)
     call set_obs_def_key(obs_def, obs_num)
     call set_obs_def(obs, obs_def)
  
     if (.not. locations_only) then
        obs_val(1) = obsval
        qc_val(1)  = qc
        call set_obs_values(obs, obs_val)
        call set_qc(obs, qc_val)
     endif

     call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

     obs_num = obs_num + 1

   enddo LEVELS

  ! clean up and loop if there is another input file
  deallocate( lat, lon, hght, elecd)

  filenum = filenum + 1

enddo fileloop

! done with main loop.  if we added any new obs to the sequence, write it out.
if (obs_num > 1) call write_obs_seq(obs_seq, output_file)

! cleanup memory
call destroy_obs_sequence(obs_seq)
call destroy_obs(obs)   ! do not destroy prev_obs, which is same as obs

write(string1, *) 'processed ', filenum-1, ' total profiles'
call error_handler(E_MSG, routine, string1, source, revision, revdate)

if (numrejected > 0) then
   write(string1, *) numrejected, ' profiles rejected for bad incoming data.'
   call error_handler(E_MSG, routine, string1, source, revision, revdate)
endif

call finalize_utilities()

! END OF MAIN ROUTINE

contains

! local subroutines/functions follow


!-----------------------------------------------------------------------
!> finds the vertical levels closest to the desired height level and
!> the weights to perform vertical interpolation.
!>
!> created June 2008, Ryan Torn NCAR/MMM
!> modified August 2017, Tim Hoar NCAR - changed the order of the
!> argument declarations so compilers can exploit knowing the array lengths

subroutine interp_height_wght(hght, level, iz, zgrid, wght)

integer,  intent(in)  :: iz       !! number of vertical levels
real(r8), intent(in)  :: hght(iz) !! height levels in column
real(r8), intent(in)  :: level    !! height level to interpolate to
integer,  intent(out) :: zgrid    !! index of lowest level for interpolation
real(r8), intent(out) :: wght     !! weight to give to the lower level in interpolation

integer :: k, klev, kbot, ktop, kinc, kleva

if ( hght(1) > hght(iz) ) then
  kbot = iz  ;   ktop = 1   ;  kinc = -1   ;  kleva = 0
else
  kbot = 1   ;   ktop = iz  ;  kinc = 1    ;  kleva = 1
endif

if ( (hght(kbot) <= level) .AND. (hght(ktop) >= level) ) then

  do k = kbot, ktop, kinc  !  search for appropriate level
    if ( hght(k) > level ) then
      klev = k - kleva
      exit
    endif
  enddo

  ! compute the weights
  zgrid = klev
  wght  = (level-hght(klev+1)) / (hght(klev) - hght(klev+1))

else

  zgrid = -1
  wght  = 0.0_r8

endif

end subroutine interp_height_wght


!-----------------------------------------------------------------------
!> interpolate between 2 longitude values taking into account wrapping at -180 degrees
!> returns interpolated longitude between 0 and 360 degrees.

function compute_lon_wrap(lon1, lon2, weight)

real(r8), intent(in) :: lon1   !! longitude in degrees between -180 and +180
real(r8), intent(in) :: lon2   !! longitude in degrees between -180 and +180
real(r8), intent(in) :: weight !! interpolation weight between lon1 and lon2 (0 to 1)
real(r8) :: compute_lon_wrap

! if the longitudes are the same sign (both negative or both positive)
! then do the interpolation with the original values.  if the signs
! are different then we need to decide if they are crossing 0 (where we
! still use the original values) or if they are crossing the -180/180 line
! and we have to wrap the negative value.
!
! to decide between the 0 and 180 cases, take the positive value and subtract 
! the negative value (which adds it on) and see if the sum is > 180.  if not, 
! we're at the 0 crossing and we do nothing.  if yes, then we add 360 to the 
! negative value and interpolate between two positive values.   in either case
! once we have the result, if it's < 0 add 360 so the longitude returned is
! between 0 and 360 in longitude.
!
! this does not try to do anything special if the profile is tracking directly
! over one of the poles.  this is because at the exact poles all longitudes are
! identical, so being off in the longitude in any direction won't be a large
! difference in real distance.
!
! created nancy collins NCAR/IMAGe

real(r8) :: lon1a, lon2a, lono

! r/w temporaries in case we have to change the value.  
lon1a = lon1
lon2a = lon2

! if different signs and crossing the -180/180 boundary, 
! add 360 to the negative value.
if (lon1 <= 0.0_r8 .and. lon2 >= 0.0_r8) then
   if (lon2 - lon1 > 180.0_r8) lon1a = lon1a + 360.0_r8
else if (lon1 >= 0.0_r8 .and. lon2 <= 0.0_r8) then
   if (lon1 - lon2 > 180.0_r8) lon2a = lon2a + 360.0_r8
endif

! linear interpolation, and make return value between 0 and 360.
lono  = weight * lon1a  + (1.0_r8 - weight) * lon2a
if (lono < 0.0_r8) lono = lono + 360.0_r8

compute_lon_wrap = lono

end function compute_lon_wrap


!-----------------------------------------------------------------------
!> computes the observation error given a variety of algorithms
!>
!> METHOD_CONSTANT          percent = factor
!> METHOD_SCALED            percent = factor * obsval
!> METHOD_LOOKUP            percent = f(solar local time, magnetic equator latitude ...)
!> METHOD_SCALED_LOOKUP     percent = factor * f(solar local time, ...) * obsval

function electron_density_error(lon, lat, hght, ihour, imin, method, factor, obsval)

real(r8), intent(in)  :: lon     !! input real value geometric height [km]
real(r8), intent(in)  :: lat     !! latitude in degrees
real(r8), intent(in)  :: hght    !! height of observation
integer,  intent(in)  :: ihour   !! hour of day ... UTC
integer,  intent(in)  :: imin    !! minute of day ... UTC
integer,  intent(in)  :: method  !! integer describing algorithm
real(r8), intent(in)  :: factor  !! multiplier to scale all the errors up or down for testing
real(r8), intent(in)  :: obsval  !! observation value
real(r8)              :: electron_density_error !! output electron density error variance

real(r8) :: percent

if (method == METHOD_CONSTANT) then
   electron_density_error = factor
   return
endif

if (method == METHOD_SCALED) then
   electron_density_error = factor * obsval
   return
endif

! The lookup table is read once in ionprf_obserr_percent().
percent = ionprf_obserr_percent(lon, lat, hght, ihour, imin)

if (method == METHOD_LOOKUP) then
   electron_density_error = percent
   return
endif

if (method == METHOD_SCALED_LOOKUP) then
   electron_density_error = factor * percent * obsval
   return
endif

end function electron_density_error



!-----------------------------------------------------------------------
!> The lookup table is specified by 'observation_error_file'
!> computes the observation error percentage for an electron density observation.
!> These numbers are taken from a Liu's and Yue's paper.
!>
!>    created by I-TE LEE NCAR/HAO & NCU, 01/26/2010

function ionprf_obserr_percent(lone, late, hghte, houre, mine)

real(r8), intent(in) ::  lone  !! longitude of electron density observation
real(r8), intent(in) ::  late  !! latitude of electron density observation
real(r8), intent(in) :: hghte  !! height of electron density observation (in km)
integer,  intent(in) :: houre  !! UTC hour of day
integer,  intent(in) ::  mine  !! UTC minute of hour

character(len=*), parameter :: routine = 'ionprf_obserr_percent'

!add F3/C ionprf observation error /add by ITL 2011.01.31
integer, parameter :: NUMZ   = 15
integer, parameter :: NUMLAT = 37
integer, parameter :: NUMTIME = 25
real(r8), save :: percent(NUMZ,NUMLAT,NUMTIME)
logical,  save :: read_file = .true.

real(r8) :: mag_eq(73), slt, ionprf_obserr_percent, mq, nlat
real(r8) :: err_top1, err_bottom1, err_obs1, err_top2, err_bottom2, err_obs2, err_obs
integer  :: altc, latc, ltc, lonc
integer  :: ios, ncid

!Calculate the magnetic equator latitude (refer to IGRF output in 2008)
data mag_eq/ 11.37_r8, 11.24_r8, 11.06_r8, 10.83_r8, 10.50_r8, 10.03_r8, &
              9.43_r8,  8.75_r8,  8.11_r8,  7.61_r8,  7.32_r8,  7.22_r8, &
              7.27_r8,  7.42_r8,  7.61_r8,  7.82_r8,  8.01_r8,  8.16_r8, &
              8.23_r8,  8.21_r8,  8.08_r8,  7.89_r8,  7.75_r8,  7.64_r8, &
              7.64_r8,  7.72_r8,  7.81_r8,  7.85_r8,  7.77_r8,  7.55_r8, &
              7.18_r8,  6.64_r8,  5.93_r8,  5.07_r8,  4.12_r8,  3.16_r8, &
              2.26_r8,  1.47_r8,  0.80_r8,  0.23_r8, -0.28_r8, -0.78_r8, &
             -1.30_r8, -1.86_r8, -2.46_r8, -3.08_r8, -3.69_r8, -4.31_r8, &
             -4.97_r8, -5.71_r8, -6.56_r8, -7.52_r8, -8.58_r8, -9.69_r8, &
             -10.78_r8, -11.71_r8, -12.34_r8, -12.46_r8, -11.90_r8, -10.61_r8, &
             -8.63_r8, -6.08_r8, -3.17_r8, -0.14_r8,  2.73_r8,  5.23_r8, &
              7.28_r8,  8.85_r8,  9.98_r8, 10.75_r8, 11.20_r8, 11.38_r8, &
             11.38_r8/

! read the lookup table once.
! BTW: if the variable is the wrong shape, nc_get_variable() fails
! with no information about the expected vs actual shape

if (read_file) then
   ios = nf90_open(observation_error_file, NF90_NOWRITE, ncid)
   call nc_check(ios, routine, 'opening file "'//trim(observation_error_file)//'"')
   call nc_get_variable(ncid,'percent', percent, routine, observation_error_file)
   call nc_check(nf90_close(ncid), routine, observation_error_file)

   read_file = .false.
endif

!Convert longitude to solar local time (simple method: longitude difference)
slt = INT(NINT( (lone/15 + houre + mine/60) *10))
if    ( slt > 240) then
    slt = slt - 240
elseif( slt < 0) then
    slt = slt + 240
endif

lonc = INT(lone/5+1)
mq   = mag_eq(lonc) + (mag_eq(lonc+1) - mag_eq(lonc)) / 5 * (lone+5-5*lonc)
nlat = late-mq

if(nlat > 90) then
    nlat = 90
elseif(nlat < -90) then
    nlat = -90
endif

!Linear interpolation the error percentage
latc = INT((nlat+90)/5+1)
altc = INT((hghte-100)/50+1)   ! table is (100 + i*50)  in km-space
 ltc = INT(NINT(slt/10)+1)

if (verbose > 2) then
   ! TJH for debug purposes only
   write(*,*)'hghte  is ',hghte, ';',   'nlat is ',nlat,  ';',   'slt is ',slt,  ';'
   write(*,*)'altc   is ',altc,  ';',   'latc is ',latc,  ';',   'ltc is ',ltc,  ';'
   write(*,*)'altc+1 is ',altc+1,';', 'latc+1 is ',latc+1,';', 'ltc+1 is ',ltc+1,';'
endif

! err_top    = percent(altc+1,latc,ltc)+(percent(altc+1,latc+1,ltc)-percent(altc,latc,ltc))/5*(nlat+95-latc*5)
! err_bottom = percent(altc,latc,ltc)+(percent(altc,latc+1,ltc)-percent(altc+1,latc,ltc))/5*(nlat+95-latc*5)
! err_obs    = err_bottom+ (err_top-err_bottom)/50*(hghte-50-altc*50)

!!3-Dimensional Linear interpolation of the observation error percentage
err_top1    = percent(altc+1,latc,ltc)+(percent(altc+1,latc+1,ltc)-percent(altc,latc,ltc))/5*(nlat+95-latc*5)
err_bottom1 = percent(altc,latc,ltc)+(percent(altc,latc+1,ltc)-percent(altc+1,latc,ltc))/5*(nlat+95-latc*5)
err_obs1    = err_bottom1 + (err_top1-err_bottom1)/50*(hghte-50-altc*50)

err_top2    = percent(altc+1,latc,ltc+1)+(percent(altc+1,latc+1,ltc+1)-percent(altc,latc,ltc+1))/5*(nlat+95-latc*5)
err_bottom2 = percent(altc,latc,ltc+1)+(percent(altc,latc+1,ltc+1)-percent(altc+1,latc,ltc+1))/5*(nlat+95-latc*5)
err_obs2    = err_bottom2 + (err_top2-err_bottom2)/50*(hghte-50-altc*50)

err_obs     = 10.0_r8 + err_obs1 + (err_obs2-err_obs1)/10 * (slt+10-ltc*10)

ionprf_obserr_percent = err_obs

end function ionprf_obserr_percent


!-----------------------------------------------------------------------
!> ensures that the observation error method is supported
!> observation_error_method   is a namelist item specifying the method
!>
!> METHOD_CONSTANT          = 1 = 'constant'
!> METHOD_SCALED            = 2 = 'scaled'
!> METHOD_LOOKUP            = 3 = 'lookup'
!> METHOD_SCALED_LOOKUP     = 4 = 'scaled_lookup'

integer function check_error_method()

character(len=128) :: mystring

mystring = observation_error_method
call to_upper(mystring)

select case(trim(mystring))
   case ('CONSTANT')
      check_error_method = METHOD_CONSTANT
   case ('SCALED')
      check_error_method = METHOD_SCALED
   case ('LOOKUP')
      check_error_method = METHOD_LOOKUP
   case ('SCALED_LOOKUP')
      check_error_method = METHOD_SCALED_LOOKUP
   case default
      write(string1,*)'unknown method "'//trim(observation_error_method)//'"' 
      write(string2,*)'valid values are "constant", "scaled", "lookup", &
                          & or "scaled_lookup"'
      write(string3,*)'this is specified by "observation_error_method"'
      call error_handler(E_ERR,'check_error_method',string1, &
                 source, revision, revdate, text2=string2, text3=string3)
end select

end function check_error_method


!-----------------------------------------------------------------------
!> If there are negative values above the lowest desired level - the profile
!> is discarded.

function bad_profile_check(lowest_height, profile_heights, density, filename) 

real(r8),         intent(in) :: lowest_height
real(r8),         intent(in) :: profile_heights(:)
real(r8),         intent(in) :: density(:)
character(len=*), intent(in) :: filename
logical              :: bad_profile_check

integer :: nobs, lowest_index
real(r8) :: weight ! do not need value for this application

bad_profile_check = .false. ! assume profile is good until proven otherwise

! Need to find the index of the lowest layer 

nobs = size(profile_heights)
call interp_height_wght(profile_heights, lowest_height, nobs, lowest_index, weight)

if (lowest_index < 0 .or. lowest_index == nobs) then
   write(string1,*)'unable to determine an observation below ',lowest_height
   write(string2,*)'"interp_height_wght" returned ',lowest_index
   write(string3,*)'expected something between 1 and ',nobs-1
   call error_handler(E_ERR,'bad_profile_check',string1, &
             source, revision, revdate, text2=string2, text3=string3)
endif

if (verbose > 1) then
   write(*,*)'lowest height ',lowest_height
   write(*,*)'profile    is ',profile_heights(lowest_index  ), density(lowest_index  )
   write(*,*)'profile+1  is ',profile_heights(lowest_index+1), density(lowest_index+1)
endif

! If there are any negative electron densities above this layer it is a bad profile

if ( any(density(lowest_index:nobs) < 0.0_r8) ) then
   bad_profile_check = .true.
   if (verbose > 0 ) &
      call error_handler(E_MSG, 'bad_profile_check', 'rejecting '//trim(filename))
endif

end function bad_profile_check



subroutine summarize(var,varname,filename)

real(r8),         intent(in) :: var(:)
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: filename

integer :: ilast, ilength

! search for last '/' and trim filename
ilength = len_trim(filename)
ilast   = index(filename, '/', back=.true.) + 1

write(*,*)filename(ilast:ilength),' ',trim(varname),' range is ',minval(var),maxval(var)

end subroutine summarize


end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

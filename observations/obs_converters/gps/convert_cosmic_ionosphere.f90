! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> a version of the netcdf -> dart obs_seq converter for ionosphere profiles
!> the file type from the CDAAC data site is 'ionPrf'.
!>
!> i don't know where the obs errors are going to come from. 

program convert_cosmic_ionosphere

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_cosmic_ionosphere - program that reads a COSMIC GPS ionosphere
!                               profile and writes the data to a DART 
!                               observation sequence file. 
!
!     copied from GPS atmospheric radio/occultation version
!     which was created June 2008 by Ryan Torn, NCAR/MMM
!     nsc 11 mar 2016  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use          types_mod, only : r8
use   time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                               increment_time, get_time, set_date, operator(-),  &
                               print_date
use      utilities_mod, only : initialize_utilities, find_namelist_in_file,    &
                               check_namelist_read, nmlfileunit, do_nml_file,   &
                               get_next_filename, error_handler, E_ERR, E_MSG, &
                               nc_check, find_textfile_dims, do_nml_term, finalize_utilities
use       location_mod, only : VERTISHEIGHT, set_location
use   obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,       &
                               static_init_obs_sequence, init_obs, destroy_obs, &
                               write_obs_seq, init_obs_sequence, get_num_obs,   &
                               insert_obs_in_seq, destroy_obs_sequence,         &
                               set_copy_meta_data, set_qc_meta_data, set_qc,    & 
                               set_obs_values, set_obs_def, insert_obs_in_seq
use   obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                               set_obs_def_error_variance, set_obs_def_location, &
                               set_obs_def_key
use       obs_kind_mod, only : COSMIC_ELECTRON_DENSITY
use  obs_utilities_mod, only : add_obs_to_seq

use           netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=512) :: msgstring, msgstring2, next_infile
character (len=NF90_MAX_NAME)  :: name
character (len=6)   :: subset
integer :: ncid, varid, nlevels, k, nfiles, num_new_obs, oday, osec, &
           iyear, imonth, iday, ihour, imin, isec, glat, glon, zloc, obs_num, &
           io, iunit, nobs, filenum, dummy, numrejected
character (len=1) :: badqc
logical :: file_exist, first_obs, from_list = .false.
real(r8) :: hght_miss, refr_miss, azim_miss, elecd_miss, oerr,   & 
            qc, lato, lono, hghto, refro, azimo, wght, nx, ny,   & 
            nz, rfict, obsval, phs, obs_val(1), qc_val(1)

real(r8), allocatable :: lat(:), lon(:), hght(:), elecd(:)

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

integer, parameter :: NMAXLEVELS = 200   !  max number of observation levels

real(r8) :: obs_levels(NMAXLEVELS) = -1.0_r8   ! in kilometers
character(len=256) :: ion_netcdf_file     = 'cosmic_ion_input.nc'
character(len=256) :: ion_netcdf_filelist = 'cosmic_ion_input_list'
character(len=256) :: ion_out_file        = 'obs_seq.iondens'

namelist /convert_cosmic_ion_nml/ obs_levels, ion_netcdf_file,         &
                                  ion_netcdf_filelist, ion_out_file

! initialize some values
obs_num = 1
qc = 0.0_r8
first_obs = .true.
call set_calendar_type(GREGORIAN)

!  read the necessary parameters from input.nml
call initialize_utilities()

call find_namelist_in_file("input.nml", "convert_cosmic_ion_nml", iunit)
read(iunit, nml = convert_cosmic_ion_nml, iostat = io)
call check_namelist_read(iunit, io, "convert_cosmic_ion_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=convert_cosmic_ion_nml)
if (do_nml_term()) write(     *     , nml=convert_cosmic_ion_nml)

!  count observation levels, make sure observation levels increase from 0
nlevels = 0
LEVELLOOP: do k = 1, NMAXLEVELS
  if ( obs_levels(k) == -1.0_r8 )  exit LEVELLOOP
  nlevels = k
end do LEVELLOOP
do k = 2, nlevels
  if ( obs_levels(k-1) >= obs_levels(k) ) then
    call error_handler(E_ERR, 'convert_cosmic_ionosphere',       &
                       'Observation levels should increase',  &
                       source, revision, revdate)
  end if
end do

! cannot have both a single filename and a list; the namelist must
! shut one off.
if (ion_netcdf_file /= '' .and. ion_netcdf_filelist /= '') then
  call error_handler(E_ERR, 'convert_cosmic_ionosphere',                     &
                     'One of ion_netcdf_file or filelist must be NULL', &
                     source, revision, revdate)
endif
if (ion_netcdf_filelist /= '') from_list = .true.

! need to know a reasonable max number of obs that could be added here.
if (from_list) then
   call find_textfile_dims(ion_netcdf_filelist, nfiles, dummy)
   num_new_obs = nlevels * nfiles
else
   num_new_obs = nlevels
endif

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
inquire(file=ion_out_file, exist=file_exist)
if ( file_exist ) then

   write(msgstring, *) "found existing obs_seq file, appending to ", trim(ion_out_file)
   write(msgstring2, *) "adding up to a maximum of ", num_new_obs, " new observations"
   call error_handler(E_MSG, 'convert_cosmic_ionosphere', msgstring, &
                      source, revision, revdate, text2=msgstring2)
   call read_obs_seq(ion_out_file, 0, 0, num_new_obs, obs_seq)

else

   write(msgstring,  *) "no previously existing obs_seq file, creating ", trim(ion_out_file)
   write(msgstring2, *) "with up to a maximum of ", num_new_obs, " observations"
   call error_handler(E_MSG, 'convert_cosmic_ionosphere', msgstring, &
                      source, revision, revdate, text2=msgstring2)

   call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
   call set_copy_meta_data(obs_seq, 1, 'observation')
   call set_qc_meta_data(obs_seq, 1, 'COSMIC QC')

end if


! main loop that does either a single file or a list of files.
! the data is currently distributed as a single profile per file.

filenum = 1
numrejected = 0

fileloop: do      ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then 
      next_infile = get_next_filename(ion_netcdf_filelist, filenum)
   else
      next_infile = ion_netcdf_file
      if (filenum > 1) next_infile = ''
   endif
   if (next_infile == '') exit fileloop
  
   !  open the next occultation profile
   call nc_check( nf90_open(next_infile, nf90_nowrite, ncid), 'file open', next_infile)

   ! process the profile
   call nc_check( nf90_get_att(ncid,nf90_global,'year',  iyear) ,'get_att year',   next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'month', imonth),'get_att month',  next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'day',   iday)  ,'get_att day',    next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'hour',  ihour) ,'get_att hour',   next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'minute',imin)  ,'get_att minute', next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'second',isec)  ,'get_att second', next_infile)
   
   time_obs = set_date(iyear, imonth, iday, ihour, imin, isec)
   call get_time(time_obs,  osec, oday)
   
   call nc_check( nf90_inq_dimid(ncid, "MSL_alt", varid),          'inq dimid MSL_alt', next_infile)
   call nc_check( nf90_inquire_dimension(ncid, varid, name, nobs), 'inq dim MSL_alt',   next_infile)
   
   allocate(  lat(nobs)) ;  allocate(  lon(nobs))
   allocate( hght(nobs)) ;  allocate(elecd(nobs))
   
   ! read the latitude array
   call nc_check( nf90_inq_varid(ncid, "GEO_lat", varid) ,'inq varid GEO_lat', next_infile)
   call nc_check( nf90_get_var(ncid, varid, lat)         ,'get var   GEO_lat', next_infile)
   
   ! read the latitude array
   call nc_check( nf90_inq_varid(ncid, "GEO_lon", varid) ,'inq varid GEO_lon', next_infile)
   call nc_check( nf90_get_var(ncid, varid, lon)         ,'get var   GEO_lon', next_infile)
   
   ! read the altitude array
   call nc_check( nf90_inq_varid(ncid, "MSL_alt", varid) ,'inq varid MSL_alt', next_infile)
   call nc_check( nf90_get_var(ncid, varid, hght)        ,'get_var   MSL_alt', next_infile)
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', hght_miss) ,'get_att _FillValue MSL_alt', next_infile)
   
   ! read the electron density
   call nc_check( nf90_inq_varid(ncid, "ELEC_dens", varid) ,'inq varid ELEC_dens', next_infile)
   call nc_check( nf90_get_var(ncid, varid, elecd)         ,'get var   ELEC_dens', next_infile)
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', elecd_miss) ,'get_att _FillValue elecd', next_infile)
   
   call nc_check( nf90_close(ncid) , 'close file', next_infile)
   
   !> @todo debug only
   ! if(any(hght(:) < -998)) then
   !    write(*,*) 'height array contains missing values:'
   !    do k=1, nobs
   !     if (hght(k) < -999) write(*,*) k, hght(k)
   !    enddo
   ! endif
   ! if(any(elecd(:) < -998)) then
   !    write(*,*) 'electron density array contains missing values:'
   !    do k=1, nobs
   !     if (elecd(k) < -999) write(*,*) k, elecd(k)
   !    enddo
   ! endif
   ! write(*,*) 'min/max elevation in kilometers: ', minval(hght), maxval(hght)

   obsloop2: do k = 1, nlevels
   
     call interp_height_wght(hght, obs_levels(k), nobs, zloc, wght)
     if ( zloc < 1 )  then
       ! write(*,*) 'level ', k, ' unable to find ', obs_levels(k), ' between 1 ', hght(1), ' and ', nobs, hght(nobs)
        cycle obsloop2
     endif
   
     ! lon(zloc) and lon(zloc+1) range from -180 to +180
     ! call a subroutine to handle the wrap point.
     lono = compute_lon_wrap(lon(zloc), lon(zloc+1), wght)
     lato  = wght * lat(zloc)  + (1.0_r8 - wght) * lat(zloc+1)
     hghto = wght * hght(zloc) + (1.0_r8 - wght) * hght(zloc+1)
     obsval = wght * elecd(zloc) + (1.0_r8 - wght) * elecd(zloc+1)
   
     oerr  = electron_density_error(hghto, lato, is_it_global=.true., factor=1.0_r8)
   
     call set_obs_def_location(obs_def,set_location(lono,lato,hghto*1000.0_r8,VERTISHEIGHT))
     call set_obs_def_type_of_obs(obs_def, COSMIC_ELECTRON_DENSITY)
     call set_obs_def_time(obs_def, set_time(osec, oday))
     call set_obs_def_error_variance(obs_def, oerr * oerr)
     call set_obs_def_key(obs_def, obs_num)
     call set_obs_def(obs, obs_def)
   
     obs_val(1) = obsval
     call set_obs_values(obs, obs_val)
     qc_val(1)  = qc
     call set_qc(obs, qc_val)

     !> @todo remove debug when working
     !print *, 'k, val: ', k, obs_val(1)

     call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

     obs_num = obs_num+1

   end do obsloop2

  ! clean up and loop if there is another input file
  deallocate( lat, lon, hght, elecd)

  filenum = filenum + 1

end do fileloop

! done with main loop.  if we added any new obs to the sequence, write it out.
if (obs_num > 1) call write_obs_seq(obs_seq, ion_out_file)

! cleanup memory
call destroy_obs_sequence(obs_seq)
call destroy_obs(obs)   ! do not destroy prev_obs, which is same as obs

write(msgstring, *) 'processed ', filenum-1, ' total profiles'
call error_handler(E_MSG, 'convert_cosmic_ion_obs', msgstring, source, revision, revdate)

if (numrejected > 0) then
   write(msgstring,  *) numrejected, ' profiles rejected for bad incoming quality control'
   call error_handler(E_MSG, 'convert_cosmic_ion_obs', msgstring, source, revision, revdate)
endif

call finalize_utilities()

! END OF MAIN ROUTINE

contains

! local subroutines/functions follow


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   interp_height_wght - subroutine that finds the vertical levels
!                        closest to the desired height level and the
!                        weights to give to these levels to perform
!                        vertical interpolation.
!
!    hght  - height levels in column
!    level - height level to interpolate to
!    zgrid - index of lowest level for interpolation
!    wght  - weight to give to the lower level in interpolation
!    iz    - number of vertical levels
!
!     created June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interp_height_wght(hght, level, iz, zgrid, wght)

use        types_mod, only : r8

implicit none

integer,  intent(in)  :: iz
real(r8), intent(in)  :: hght(iz), level
integer,  intent(out) :: zgrid
real(r8), intent(out) :: wght

integer :: k, klev, kbot, ktop, kinc, kleva

if ( hght(1) > hght(iz) ) then
  kbot = iz  ;   ktop = 1   ;  kinc = -1   ;  kleva = 0
else
  kbot = 1   ;   ktop = iz  ;  kinc = 1    ;  kleva = 1
end if

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

return
end subroutine interp_height_wght


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   compute_lon_wrap - interpolate between 2 longitude values taking
!                      into account the wrap at -180 degrees
!
!    lon1, lon2 - longitude in degrees between -180 and +180
!    weight - interpolation weight between lon1 and lon2 (0 to 1)
! 
!    returns interpolated longitude between 0 and 360 degrees.
!
! if the longitudes are the same sign (both negative or both positive)
! then do the interpolation with the original values.  if the signs
! are different then we need to decide if they are crossing 0 (where we
! still use the original values) or if they are crossing the -180/180 line
! and we have to wrap the negative value.

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
!     created nancy collins NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function compute_lon_wrap(lon1, lon2, weight)

real(r8), intent(in) :: lon1, lon2, weight
real(r8) :: compute_lon_wrap

real(r8) :: lon1a, lon2a, lono

! r/w temporaries in case we have to change the value.  
lon1a = lon1
lon2a = lon2


! if different signs and crossing the -180/180 boundary, add 360
! to the negative value.
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

! error options:


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   electron_density_error - function that computes the observation error
!
!    Input:
!    H   -- input real value geometric height [km]
!    lat -- latitude in degree 
!    is_it_global -- is your forecast model regional or global? (T or F)
!    factor -- quick way to scale all the errors up or down for testing
!
!    Output:
!    electron_density_error -- output electron density error
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function electron_density_error(H, lat, is_it_global, factor)
 real(r8), intent(in)  :: H
 real(r8), intent(in)  :: lat
 logical,  intent(in)  :: is_it_global
 real(r8), intent(in)  :: factor
 real(r8)              :: electron_density_error

 real(r8) :: zkm, rerr
 integer  :: kk
 
 zkm = H

 !> @todo fix the error routine.  this is returning early.  the existing code
 !> is a remainder from the lower atmosphere profile info for refractivity.
 call error_handler(E_MSG, 'the error routine needs to be adapted for electron density', &
    'the existing code is for lower atmos radio occultation obs', source, revision, revdate)
 electron_density_error = 1.0_r8
 return

 if(is_it_global) then    ! for global

 if((lat >= 20.0_r8) .or. (lat <= -20.0_r8)) then
    rerr = -1.321_r8 + 0.341_r8 * zkm - 0.005_r8 * zkm ** 2
 else
    if(zkm > 10.0_r8) then
       rerr = 2.013_r8 - 0.060_r8 * zkm + 0.0045_r8 * zkm ** 2
    else
       rerr = -1.18_r8 + 0.058_r8 * zkm + 0.025_r8 * zkm ** 2
    endif
 endif

 else     ! for regional 

 if((lat >= 20.0_r8) .or. (lat <= -20.0_r8)) then
    if(zkm > 10.0_r8) then
       rerr = -1.321_r8 + 0.341_r8 * zkm - 0.005_r8 * zkm ** 2
    else
       rerr = -1.2_r8 + 0.065_r8 * zkm + 0.021_r8 * zkm ** 2
    endif
 endif

 endif    ! if(is_it_global) then

 if (factor /= 1.0_r8) then
    electron_density_error = 1./(abs(exp(rerr))*factor)
 else
    electron_density_error = 1./abs(exp(rerr))
 endif

 return
end function electron_density_error

! end error options

! next function is currently unused in this converter:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    compute_geopotential_height 
!    subroutine converts geometric height to geopotential height
!
!    Input:
!    H   -- input real value geometric height [m]
!    lat -- latitude in degree 
!
!    Output:
!    Z -- output real value geopotential height [m]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function compute_geopotential_height(H, lat)
 real(r8), intent(in)  :: H 
 real(r8), intent(in)  :: lat
 real(r8)              :: compute_geopotential_height

! -----------------------------------------------------------------------*/
   real(r8) :: pi2, latr
   real(r8) :: semi_major_axis, semi_minor_axis, grav_polar, grav_equator
   real(r8) :: earth_omega, grav_constant, flattening, somigliana
   real(r8) :: grav_ratio, sin2, termg, termr, grav, eccentricity

!  Parameters below from WGS-84 model software inside GPS receivers.
   parameter(semi_major_axis = 6378.1370d3)    ! (m)
   parameter(semi_minor_axis = 6356.7523142d3) ! (m)
   parameter(grav_polar = 9.8321849378)        ! (m/s2)
   parameter(grav_equator = 9.7803253359)      ! (m/s2)
   parameter(earth_omega = 7.292115d-5)        ! (rad/s)
   parameter(grav_constant = 3.986004418d14)   ! (m3/s2)
   parameter(grav = 9.80665d0)                 ! (m/s2) WMO std g at 45 deg lat
   parameter(eccentricity = 0.081819d0)        ! unitless
   parameter(pi2 = 3.14159265358979d0/180.d0)

!  Derived geophysical constants
   parameter(flattening = (semi_major_axis-semi_minor_axis) / semi_major_axis)

   parameter(somigliana = (semi_minor_axis/semi_major_axis)*(grav_polar/grav_equator)-1.d0)

   parameter(grav_ratio = (earth_omega*earth_omega * &
                semi_major_axis*semi_major_axis * semi_minor_axis)/grav_constant)

!  Sanity Check
   if(lat.gt.90 .or. lat.lt.-90) then
      print*, 'compute_geopotential_height: Latitude is not between -90 and 90 degrees: ',lat
      return
   endif

   latr = lat * (pi2)        ! in radians
   sin2  = sin(latr) * sin(latr)
   termg = grav_equator * ( (1.d0+somigliana*sin2) / &
           sqrt(1.d0-eccentricity*eccentricity*sin2) )
   termr = semi_major_axis / (1.d0 + flattening + grav_ratio - 2.d0*flattening*sin2)

   compute_geopotential_height = (termg/grav)*((termr*H)/(termr+H))
   !print '(A,3f15.5)','compute_geopotential_height:',&
   !                   H,compute_geopotential_height,H-compute_geopotential_height

end function compute_geopotential_height

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

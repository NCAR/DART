! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program convert_sat_chl

!	float chlor_a(lat, lon) ;
!		chlor_a:long_name = "Chlorophyll Concentration, OCI Algorithm" ;
!		chlor_a:units = "mg m^-3" ;
!		chlor_a:standard_name = "mass_concentration_of_chlorophyll_in_sea_water" ;
!		chlor_a:_FillValue = -32767.f ;
!		chlor_a:valid_min = 0.001f ;
!		chlor_a:valid_max = 100.f ;
!		chlor_a:reference = "Hu, C., Lee Z., and Franz, B.A. (2012). Chlorophyll-a algorithms for oligotrophic oceans: A novel approach based on three-band reflectance difference, JGR" ;
!		chlor_a:display_scale = "log" ;
!		chlor_a:display_min = 0.01f ;
!		chlor_a:display_max = 20.f ;


use types_mod,            only : r8, digits12
use time_manager_mod,     only : time_type, set_calendar_type, GREGORIAN, &
                                 set_time, get_time, print_time, &
                                 set_date, get_date, print_date, &
                                 operator(+), operator(-)
use utilities_mod,        only : initialize_utilities, find_namelist_in_file, &
                                 check_namelist_read, nmlfileunit, &
                                 error_handler, E_ERR, E_MSG, &
                                 finalize_utilities, do_nml_file, do_nml_term
use location_mod,         only : VERTISSURFACE, set_location
use obs_sequence_mod,     only : obs_type, obs_sequence_type, init_obs, &
                                 static_init_obs_sequence, init_obs_sequence, &
                                 set_copy_meta_data, set_qc_meta_data, &
                                 get_num_obs, write_obs_seq, destroy_obs_sequence
use obs_utilities_mod,    only : create_3d_obs, add_obs_to_seq
use obs_kind_mod,         only : QTY_SURFACE_CHLOROPHYLL, OCEAN_COLOR
use netcdf_utilities_mod, only : nc_check, nc_open_file_readonly, nc_close_file, &
                                 nc_get_variable, nc_get_attribute_from_variable, &
                                 nc_get_dimension_size, nc_get_variable_size 

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'convert_sat_chl'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

real(r8), parameter :: qc = 0.0_r8

character(len=256) :: output_file
character(len=512) :: string1, string2

integer :: ncid, varid, io, iunit
integer :: oday, osec, iday, isec
integer :: year, month, day, hour, minutes, seconds
integer :: num_new_obs, nmissing
integer :: i, j, nlat, nlon, ndays
integer :: itime

logical :: first_obs

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, prev_time
type(time_type)         :: base_time, delta_time

real(digits12), allocatable :: time(:)
real(r8), allocatable       :: lat(:), lon(:)
real(r8), allocatable       :: chl(:,:)
real(r8)                    :: missing_value
real(r8)                    :: log_chl

! Observation uncertainty: MAny studies indicate that 
! obs error sd is 35% [e.g., Moore et al. 2009]
! Discussed this with Ariane Verdy from UCSD, 
! chl = log(data) to go from lognormal to normal pdf
! if in the lognormal space sd is 0.5 the data value, 
! then in the trasnformed space obs_err_sd is simply 0.35 
real(r8) :: chl_error_sd = 0.35_r8

! If we need to mask other sea areas
! For instance, we are only interested in the Red Sea 
! and not the Persian/Arabian Gulf
real(r8) :: lon_mask = 47.5_r8
real(r8) :: lat_mask = 26.5_r8

!------------------------------------------------------------------------
!  Declare namelist parameters

character(len=256) :: file_in         = 'chl_in.nc'
character(len=256) :: file_out        = 'obs_seq.chl'
real(r8)           :: chl_thresh      = 0.03_r8
integer            :: subsample_intv  = 1
logical            :: special_mask    = .false.
logical            :: debug           = .false.

namelist /convert_sat_chl_nml/ file_in,         & 
                               file_out,        &
                               chl_thresh,      &
                               subsample_intv,  &
                               special_mask,    &
                               debug

!------------------------------------------------------------------------
! start of executable code
!------------------------------------------------------------------------

! Read and Record the necessary parameters from input.nml
call initialize_utilities()

call find_namelist_in_file('input.nml', 'convert_sat_chl_nml', iunit)
read(iunit, nml = convert_sat_chl_nml, iostat = io)

if (do_nml_file()) write(nmlfileunit, nml=convert_sat_chl_nml)
if (do_nml_term()) write(     *     , nml=convert_sat_chl_nml)

call set_calendar_type(GREGORIAN)

ncid  = nc_open_file_readonly(file_in) 
ndays = nc_get_dimension_size(ncid, 'time', source)
nlat  = nc_get_dimension_size(ncid, 'lat' , source)
nlon  = nc_get_dimension_size(ncid, 'lon' , source)

allocate(time(ndays), lat(nlat), lon(nlon))

call nc_get_variable(ncid, 'time', time, source) 
call nc_get_variable(ncid, 'lat' ,  lat, source)
call nc_get_variable(ncid, 'lon' ,  lon, source) 

! ensure longitudes are [0,360] 
where(lon < 0.0_r8) lon = lon + 360.0_r8

base_time = set_base_time(ncid)

num_new_obs = nlon*nlat

allocate(chl(nlon,nlat))

call static_init_obs_sequence()

io = nf90_inq_varid(ncid, 'chlor_a', varid) 
call nc_check(io, source, context='getting chl variable ID', ncid=ncid)
call nc_get_attribute_from_variable(ncid, 'chlor_a', '_FillValue', missing_value, source)

TIMELOOP: do itime = 1,ndays

   ! convert to integer days and seconds, and add on to reference time.
   iday = time(itime)
   isec = (time(itime) - iday) * 86400
   delta_time = set_time(isec, iday)
   obs_time = base_time + delta_time
   call get_time(obs_time,  osec, oday)

   call print_time(obs_time, str='obs time is ')
   call print_date(obs_time, str='obs date is ')

   call get_date(obs_time, year, month, day, hour, minutes, seconds)

   seconds = seconds + (hour*60 + minutes)*60

   write(string1,'(i4,''-'',i2.2,''-'',i2.2,''-'',i5.5)') year,month,day,seconds
   write(output_file,'(A)') trim(file_out)//'.'//trim(string1)
   write(*,*)'output file is ',trim(output_file)

   io = nf90_get_var(ncid, varid, chl, start=(/1,1,itime/), count=(/nlon,nlat,1/))
   call nc_check(io, source, context='get_var chl', ncid=ncid)

   first_obs = .true.
   nmissing  = 0
   call init_obs(     obs, num_copies, num_qc)
   call init_obs(prev_obs, num_copies, num_qc)
   call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
   call set_copy_meta_data(obs_seq, 1, 'CHL observation')
   call set_qc_meta_data(  obs_seq, 1, 'CHL QC')

   latloop: do j = 1, nlat, subsample_intv
      lonloop: do i = 1, nlon, subsample_intv

         if (special_mask) then 
            ! Upper right corner of the Red Sea domain (Persian/Arabian Gulf)
            if (lon(i) >= lon_mask .and. lat(j) >= lat_mask .and. chl(i, j) /= missing_value) then  
               ! Make sure it's masked
               chl(i,j) = missing_value
            endif
         endif

         if (chl(i,j) == missing_value) then
            nmissing = nmissing + 1
            cycle lonloop
         endif

         ! log-transform of the data
         ! Assimilate the logarithm of the data
         log_chl = log10(max(chl_thresh, chl(i, j))) 

         call create_3d_obs(lat(j), lon(i), 0.0_r8, VERTISSURFACE, log_chl, &
                            OCEAN_COLOR, chl_error_sd, oday, osec, qc, obs)
         call add_obs_to_seq(obs_seq, obs, obs_time, prev_obs, prev_time, first_obs)

      enddo lonloop
   enddo latloop

   ! if we added any obs to the sequence, write it out to a file now.
   if ( get_num_obs(obs_seq) > 0 ) then
      if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
      if (debug) print *, '                  skipping = ', nmissing
      call write_obs_seq(obs_seq, output_file)
   else
      write(string1,*)'no observations for output file'
      write(string2,*)'"'//trim(output_file)//'"'
      call error_handler(E_MSG, source, string1, text2=string2)
   endif

   call destroy_obs_sequence(obs_seq)

end do TIMELOOP

call nc_close_file(ncid, source)
 
call error_handler(E_MSG, source, 'Finished successfully.')
call finalize_utilities()

contains


!> 	time:units = "days since 1601-01-01 00:00:00" ;
!
function set_base_time(ncid)

integer, intent(in) :: ncid
type(time_type)     :: set_base_time

character(len=256) :: timeunits
integer :: io
integer :: year, month, day, hour, minute, second

call nc_get_attribute_from_variable(ncid,'time','units',timeunits,'set_base_time')

read(timeunits,100,iostat=io) year, month, day, hour, minute, second

if (io /= 0) then 
   write(string1,*)'unable to read time base'
   call error_handler(E_ERR, 'set_base_time', string1, &
              source, revision, revdate, text2=timeunits)
endif
            
100 format(11x,i4,5(1x,i2))

set_base_time = set_date(year, month, day, hour, minute, second)

if (debug) then
   write(*,*)'time units is ',trim(timeunits)
   call print_time(set_base_time, str='obs time is ')
   call print_date(set_base_time, str='obs date is ')
endif

end function set_base_time

end program convert_sat_chl

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

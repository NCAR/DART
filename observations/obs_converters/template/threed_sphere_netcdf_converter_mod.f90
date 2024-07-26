! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

program threed_netcdf_converter_mod

!> title = "Generalized 3D Sphere NetCDF Observation Converter"
!> institution = "NCAR" ;
!> source = "NCAR/DAReS" ;
!> comment = "Generalized converter for 3D sphere data from NetCDF files" ;
!> references = "http://www.image.ucar.edu/DAReS/DART/DART_download" ;
!> dataset_title = "Generalized 3D Sphere NetCDF Data" ;

use         types_mod, only : r8, digits12, MISSING_R8

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_time, get_time, print_time, &
                              set_date, get_date, print_date, &
                              operator(+), operator(-)
                              
use     utilities_mod, only : initialize_utilities, find_namelist_in_file, &
                              check_namelist_read, nmlfileunit, &
                              error_handler, E_ERR, E_MSG, &
                              finalize_utilities, do_nml_file, do_nml_term
                              
use      location_mod, only : get_location, location_type, set_location, VERTISSURFACE, VERTISHEIGHT, VERTISUNDEF

use  obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs, &
                              static_init_obs_sequence, init_obs_sequence, &
                              set_copy_meta_data, set_qc_meta_data, &
                              get_num_obs, write_obs_seq, destroy_obs_sequence
                              
use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

! set the kind of observation here
use obs_kind_mod, only : YOUR_OBS_KIND, get_kind_index

use netcdf_utilities_mod, only : nc_check, nc_open_file_readonly, nc_close_file, &
                                 nc_get_variable, nc_get_attribute_from_variable, &
                                 nc_get_dimension_size, nc_get_variable_size
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = & '$URL$'
character(len=*), parameter :: revision = '$Revision$'
character(len=*), parameter :: revdate  = '$Date$'
character(len=*), parameter :: routine  = 'threed_converter_mod'

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

real(r8), parameter :: qc = 0.0_r8        ! default QC value

! variables for file handling and observations
character(len=256) :: output_file
character(len=512) :: string1, string2

integer :: ncid, varid, io, iunit
integer :: oday, osec, iday, isec
integer :: year, month, day, hour, minutes, seconds
integer :: num_new_obs, nmissing
integer :: i, j, k, nlat, nlon, ndays
integer :: itime

logical :: first_obs

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, prev_time
type(time_type)         :: base_time, delta_time

real(digits12), allocatable :: time(:)
real(r8), allocatable :: lat(:), lon(:), level(:)
! declare name of variable here
real(r8), allocatable :: your_variable(:,:,:)
real(r8) :: missing_value

!------------------------------------------------------------------------
!  Declare namelist parameters

! your_variable_error_std - instrument and representativeness error (std)
real(r8)           :: your_variable_error_std    = 0.3_r8
character(len=256) :: input_file       = 'input.nc'
character(len=256) :: output_file_base = 'obs_seq'
logical            :: debug            = .false.
integer            :: subsample_intv   = 1

namelist /netcdf_to_obs_nml/ input_file, output_file_base, &
                              debug, subsample_intv, your_variable_error_std

!------------------------------------------------------------------------
! Start of executable code
!------------------------------------------------------------------------

! read and record standard parameters from input namelist
call initialize_utilities('netcdf_to_obs', .true., .true.)
call find_namelist_in_file('input.nml', 'netcdf_to_obs_nml', iunit)
read(iunit, nml = netcdf_to_obs_nml, iostat = io)

if (do_nml_file()) write(nmlfileunit, nml=netcdf_to_obs_nml)
if (do_nml_term()) write(     *     , nml=netcdf_to_obs_nml)

call set_calendar_type(GREGORIAN)

! open the NetCDF file and read dimension sizes
ncid  = nc_open_file_readonly(input_file, routine) 
ndays = nc_get_dimension_size(ncid, 'time', routine)
nlat  = nc_get_dimension_size(ncid, 'lat', routine)
nlon  = nc_get_dimension_size(ncid, 'lon', routine)
nlev  = nc_get_dimension_size(ncid, 'level', routine) ! if levels are present

! allocate arrays for time, latitude, and longitude
allocate(time(ndays), lat(nlat), lon(nlon), level(nlev))

! read time, latitude, and longitude variables from NetCDF file
call nc_get_variable(ncid, 'time', time, routine) 
call nc_get_variable(ncid, 'lat',  lat,  routine) 
call nc_get_variable(ncid, 'lon',  lon,  routine) 
call nc_get_variable(ncid, 'level',  level,  routine) 

! ensure longitudes are [0,360] 
where(lon < 0.0_r8) lon = lon + 360.0_r8

! set the base time from the NetCDF file attributes
base_time = set_base_time(ncid)

num_new_obs = nlon * nlat * nlev

! allocate array for the variable of interest
allocate(your_variable(nlev, nlat, nlon))

! initialize observation sequence
call static_init_obs_sequence()

! get variable ID and missing value attribute from NetCDF file
io = nf90_inq_varid(ncid, 'YOUR_VARIABLE', varid) 
call nc_check(io, routine, context='getting YOUR_VARIABLE variable ID',ncid=ncid)
call nc_get_attribute_from_variable(ncid, 'YOUR_VARIABLE', 'missing_value', missing_value, routine)

! loop over each time step in the NetCDF file
TIMELOOP: do itime = 1,ndays

   ! time is stored in the file 2 ways: as real(double) seconds since 1981/1/1,
   ! and as 4 and 2 digit strings for year/mon/day/hr/min
   ! both of these are variables, not attributes

   ! convert to integer days and seconds, and add on to reference time.
   iday = time(itime)
   isec = (time(itime) - iday) * 86400
   delta_time = set_time(isec, iday)
   obs_time = base_time + delta_time
   call get_time(obs_time,  osec, oday)
   
   call get_date(obs_time, year, month, day, hour, minutes, seconds)

   seconds = seconds + (hour*60 + minutes)*60
   
   ! generate output file name based on date and time
   write(string1,'(i4,''-'',i2.2,''-'',i2.2,''-'',i5.5)') year,month,day,seconds
   write(output_file,'(A)') trim(output_file_base)//'.'//trim(string1)
   write(*,*)'output file is ',trim(output_file)

   ! read the variable of interest for the current time step
   io = nf90_get_var(ncid, varid, your_variable(:,:,itime), start=(/1,1,itime/), count=(/nlon,nlat,1/))
   call nc_check(io, routine, context='get_var YOUR_VARIABLE', ncid=ncid)

   first_obs = .true.
   nmissing = 0
   call init_obs(obs, num_copies, num_qc)
   call init_obs(prev_obs, num_copies, num_qc)
   call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
   call set_copy_meta_data(obs_seq, 1, 'YOUR_VARIABLE observation')
   call set_qc_meta_data(obs_seq, 1, 'YOUR_VARIABLE QC')

   ! loop over each spatial point in the 3D grid
   obslooplat: do j = 1, nlat, subsample_intv
   obslooplon: do i = 1, nlon, subsample_intv
   obslooplevel: do k = 1, nlev, subsample_intv
   
     ! skip missing values
     if (your_variable(k,j,i) == missing_value) then
        nmissing = nmissing + 1
        cycle obslooplon
     endif

     ! set observation location
     call set_location(location, lat(j), lon(i), level(k))
     ! create and add observation to sequence
     call create_3d_obs(location, your_variable(k,j,i), VERTISUNDEF, obs_time, obs_seq, first_obs)

   enddo obslooplevel
   enddo obslooplon
   enddo obslooplat

   ! if we added obs to the sequence, write it out to a file
   if ( get_num_obs(obs_seq) > 0 ) then
      if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
      if (debug) print *, '                  skipping = ', nmissing
      call write_obs_seq(obs_seq, output_file)
   else
      write(string1,*)'no observations for output file'
      write(string2,*)'"'//trim(output_file)//'"'
      call error_handler(E_MSG, routine, string1, text2=string2)
   endif

   ! destroy current observation sequence to prepare for the next time step
   call destroy_obs_sequence(obs_seq)

end do TIMELOOP

! close the NetCDF file
call nc_close_file(ncid, routine)
call error_handler(E_MSG, routine, 'Finished successfully.')
call finalize_utilities()

contains

! function to set the base time from the NetCDF file attributes
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

! subroutine to create and add a 3D observation to the sequence
subroutine create_3d_obs(location, value, vert_coord, obs_time, obs_seq, first_obs)
   type(location_type), intent(in) :: location
   real(r8), intent(in) :: value
   integer, intent(in) :: vert_coord
   type(time_type), intent(in) :: obs_time
   type(obs_sequence_type), intent(inout) :: obs_seq
   logical, intent(inout) :: first_obs

   type(obs_def_type) :: obs_def
   type(obs_type) :: obs
   integer :: kind

   ! get the kind index for the variable of interest
   call get_kind_index('YOUR_VARIABLE_KIND', kind)
   call obs_def%initialize()
   call obs_def%set(location, kind, value)

   ! set observation time
   call obs_def%set_time(obs_time)

   ! add observation to sequence
   call add_obs_to_seq(obs_seq, obs, obs_time, first_obs)
   first_obs = .false.

end subroutine create_3d_obs

end program threed_netcdf_converter_mod
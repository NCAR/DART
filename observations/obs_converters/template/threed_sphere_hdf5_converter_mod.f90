! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

program threed_hdf5_converter_mod

!> title = "Generalized 3D Sphere HDF5 Observation Converter"
!> institution = "NCAR" ;
!> source = "NCAR/DAReS" ;
!> comment = "Generalized converter for 3D sphere data from HDF5 files" ;
!> references = "http://www.image.ucar.edu/DAReS/DART/DART_download" ;
!> dataset_title = "Generalized 3D Sphere HDF5 Data" ;

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

use hdf5

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = & '$URL$'
character(len=*), parameter :: revision = '$Revision$'
character(len=*), parameter :: revdate  = '$Date$'
character(len=*), parameter :: routine  = 'threed_hdf5_converter_mod'

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

real(r8), parameter :: qc = 0.0_r8        ! default QC value

! variables for file handling and observations
character(len=256) :: output_file
character(len=512) :: string1, string2

integer :: file_id, dataset_id, dataspace_id, memspace_id, varid, io, iunit
integer :: oday, osec, iday, isec
integer :: year, month, day, hour, minutes, seconds
integer :: num_new_obs, nmissing
integer :: i, j, k, nlat, nlon, ndays, nlev
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
character(len=256) :: input_file       = 'input.h5'
character(len=256) :: output_file_base = 'obs_seq'
logical            :: debug            = .false.
integer            :: subsample_intv   = 1

namelist /hdf5_to_obs_nml/ input_file, output_file_base, &
                              debug, subsample_intv, your_variable_error_std

!------------------------------------------------------------------------
! Start of executable code
!------------------------------------------------------------------------

! read and record standard parameters from input namelist
call initialize_utilities('hdf5_to_obs', .true., .true.)
call find_namelist_in_file('input.nml', 'hdf5_to_obs_nml', iunit)
read(iunit, nml = hdf5_to_obs_nml, iostat = io)

if (do_nml_file()) write(nmlfileunit, nml=hdf5_to_obs_nml)
if (do_nml_term()) write(     *     , nml=hdf5_to_obs_nml)

call set_calendar_type(GREGORIAN)

! open the HDF5 file and read dimension sizes
call h5open_f() ! initialize HDF5 library
call h5fopen_f(input_file, H5F_ACC_RDONLY_F, file_id, io)
call h5check(io, 'opening HDF5 file')

call h5dopen_f(file_id, "time", dataset_id, io)
call h5check(io, 'opening time dataset')
call h5sget_simple_extent_dims_f(dataset_id, ndays, io)
call h5check(io, 'getting time dataset dimensions')
call h5dclose_f(dataset_id, io)
call h5check(io, 'closing time dataset')

call h5dopen_f(file_id, "lat", dataset_id, io)
call h5check(io, 'opening lat dataset')
call h5sget_simple_extent_dims_f(dataset_id, nlat, io)
call h5check(io, 'getting lat dataset dimensions')
call h5dclose_f(dataset_id, io)
call h5check(io, 'closing lat dataset')

call h5dopen_f(file_id, "lon", dataset_id, io)
call h5check(io, 'opening lon dataset')
call h5sget_simple_extent_dims_f(dataset_id, nlon, io)
call h5check(io, 'getting lon dataset dimensions')
call h5dclose_f(dataset_id, io)
call h5check(io, 'closing lon dataset')

call h5dopen_f(file_id, "level", dataset_id, io)
call h5check(io, 'opening level dataset')
call h5sget_simple_extent_dims_f(dataset_id, nlev, io)
call h5check(io, 'getting level dataset dimensions')
call h5dclose_f(dataset_id, io)
call h5check(io, 'closing level dataset')

! allocate arrays for time, latitude, and longitude
allocate(time(ndays), lat(nlat), lon(nlon), level(nlev))

! read time, latitude, and longitude variables from HDF5 file
call h5dopen_f(file_id, "time", dataset_id, io)
call h5check(io, 'opening time dataset for reading')
call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, time, ndays, io)
call h5check(io, 'reading time dataset')
call h5dclose_f(dataset_id, io)
call h5check(io, 'closing time dataset')

call h5dopen_f(file_id, "lat", dataset_id, io)
call h5check(io, 'opening lat dataset for reading')
call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, lat, nlat, io)
call h5check(io, 'reading lat dataset')
call h5dclose_f(dataset_id, io)
call h5check(io, 'closing lat dataset')

call h5dopen_f(file_id, "lon", dataset_id, io)
call h5check(io, 'opening lon dataset for reading')
call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, lon, nlon, io)
call h5check(io, 'reading lon dataset')
call h5dclose_f(dataset_id, io)
call h5check(io, 'closing lon dataset')

call h5dopen_f(file_id, "level", dataset_id, io)
call h5check(io, 'opening level dataset for reading')
call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, level, nlev, io)
call h5check(io, 'reading level dataset')
call h5dclose_f(dataset_id, io)
call h5check(io, 'closing level dataset')

! ensure longitudes are [0,360]
where(lon < 0.0_r8) lon = lon + 360.0_r8

! set the base time from the HDF5 file attributes
base_time = set_base_time(file_id)

num_new_obs = nlon * nlat * nlev

! allocate array for the variable of interest
allocate(your_variable(nlev, nlat, nlon))

! initialize observation sequence
call static_init_obs_sequence()

! get variable ID and missing value attribute from HDF5 file
call h5dopen_f(file_id, "YOUR_VARIABLE", dataset_id, io)
call h5check(io, 'opening YOUR_VARIABLE dataset')
call h5aopen_f(dataset_id, "missing_value", varid, io)
call h5check(io, 'opening missing_value attribute')
call h5aread_f(varid, H5T_NATIVE_DOUBLE, missing_value, io)
call h5check(io, 'reading missing_value attribute')
call h5aclose_f(varid, io)
call h5check(io, 'closing missing_value attribute')
call h5dclose_f(dataset_id, io)
call h5check(io, 'closing YOUR_VARIABLE dataset')

! loop over each time step in the HDF5 file
TIMELOOP: do itime = 1, ndays

   ! convert to integer days and seconds, and add on to reference time.
   iday = time(itime)
   isec = (time(itime) - iday) * 86400
   delta_time = set_time(isec, iday)
   obs_time = base_time + delta_time
   call get_time(obs_time,  osec, oday)
   
   call get_date(obs_time, year, month, day, hour, minutes, seconds)

   seconds = seconds + (hour*60 + minutes)*60
   
   ! generate output file name based on date and time
   write(string1,'(i4,''-'',i2.2,''-'',i2.2,''-'',i5.5)') year, month, day, seconds
   write(output_file,'(A)') trim(output_file_base)//'.'//trim(string1)
   write(*,*)'output file is ',trim(output_file)

   ! read the variable of interest for the current time step
   call h5dopen_f(file_id, "YOUR_VARIABLE", dataset_id, io)
   call h5check(io, 'opening YOUR_VARIABLE dataset for reading')
   call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, your_variable(:,:,itime), (/nlev, nlat, nlon/), io)
   call h5check(io, 'reading YOUR_VARIABLE dataset')
   call h5dclose_f(dataset_id, io)
   call h5check(io, 'closing YOUR_VARIABLE dataset')

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

! close the HDF5 file
call h5fclose_f(file_id, io)
call h5check(io, 'closing HDF5 file')
call h5close_f() ! finalize HDF5 library
call error_handler(E_MSG, routine, 'Finished successfully.')
call finalize_utilities()

contains

! function to set the base time from the HDF5 file attributes
function set_base_time(file_id)

integer, intent(in) :: file_id
type(time_type)     :: set_base_time

character(len=256) :: timeunits
integer :: io, attr_id
integer :: year, month, day, hour, minute, second

call h5aopen_f(file_id, "time_units", attr_id, io)
call h5check(io, 'opening time_units attribute')
call h5aread_f(attr_id, H5T_NATIVE_CHARACTER, timeunits, io)
call h5check(io, 'reading time_units attribute')
call h5aclose_f(attr_id, io)
call h5check(io, 'closing time_units attribute')

read(timeunits,100,iostat=io) year, month, day, hour, minute, second

if (io /= 0) then 
   write(string1,*)'unable to read time base'
   call error_handler(E_ERR, 'set_base_time', string1, &
              source, revision, revdate, text2=timeunits)
endif
            
100 format(11x,i4,5(1x,i2))

set_base_time = set_date(year, month, day, hour, minute, second)

if (debug) then
   write(*,*)'time units are ',trim(timeunits)
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

! ----------------------------------------------------------------------------

subroutine h5check(hdferr, message)

implicit none

integer,                    intent(in) :: hdferr
character(len=*), optional, intent(in) :: message

character(len=*), parameter :: routine = 'h5check'

character(len=512) :: string1

if (hdferr .lt. 0) then
   if (present(message)) then
      write(string1,*) trim(message),', hdf error code: ',hdferr
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   else
      write(string1,*) 'HDF Error code:',hdferr
      call error_handler(E_ERR, routine, string1, source, revision, revdate)
   end if
end if

end subroutine h5check

! ----------------------------------------------------------------------------

end program threed_hdf5_converter_mod

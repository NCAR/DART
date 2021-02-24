! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program gitm_to_netcdf

!----------------------------------------------------------------------
! purpose: interface between the GITM model and DART
!
! method: Read gitm "restart" files of model state
!         This version assumes the individual blocks of gitm data
!         (one per gitm mpi task) have been combined into a single
!         data file.  see the other converter programs for alternatives.
!         Reform fields into a NetCDF format file.
!
! USAGE:  TBD
!----------------------------------------------------------------------

use        types_mod, only : r8, digits12

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, close_file

use time_manager_mod, only : time_type, set_date, get_date, set_time, get_time, &
                             operator(-), set_calendar_type

use netcdf_utilities_mod    ! all for now

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'gitm_to_netcdf.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------


character(len=256) :: gitm_to_netcdf_2d_input_file   = '../testdata2/2DTEC_t110311_204500.bin'
character(len=256) :: gitm_to_netcdf_3d_input_file   = '../testdata2/3DUSR_t110311_204500.bin'
character(len=256) :: gitm_to_netcdf_2d_output_file  = 'gitm_2d_netcdf.nc'
character(len=256) :: gitm_to_netcdf_3d_output_file  = 'gitm_3d_netcdf.nc'

namelist /gitm_to_netcdf_nml/ gitm_to_netcdf_2d_input_file,  &
                              gitm_to_netcdf_3d_input_file,  &
                              gitm_to_netcdf_2d_output_file, &
                              gitm_to_netcdf_3d_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: iunit

! educated guess only !!  my interpretation from dumping a file.
!
! an initial 8-byte record, GIT version?  (3.13  in real format)
! 12 bytes, dimensions?  (a set of 3 integers: 120, 100, 1)
! 4 bytes, a field count? (a single int: 5)
! 5 repeats of 40 bytes, field names?  (ascii characters: Longitude, etc)
! 28 bytes, a time record? (a set of 7 ints: 2011, 3, 11, 20, 45, 0, 0)
!           (what's the last zero? milliseconds?)
!
! n reps of either:
!    96,000 bytes (one 2d variable data? 120 x 100 x      8 byte reals for 2d)
! 5,571,072 bytes (one 3d variable data? 124 x 104 x 54 x 8 byte reals for 3d)


integer :: lon_size
integer :: lat_size
integer :: alt_size

integer :: i, ncid
integer :: nvars_2D, nvars_3D

integer, parameter :: NAMELEN = 40

type var_2D
  character(len=NAMELEN) :: varname
  real(r8), allocatable :: state(:,:)
end type var_2D

type(var_2D) :: gitm_2Dvars(100)

type var_3D
  character(len=NAMELEN) :: varname
  real(r8), allocatable :: state(:,:,:)
end type var_3D

type(var_3D) :: gitm_3Dvars(100)

real(r8) :: gitm_version
real(digits12) :: real_time
integer :: t_year, t_month, t_day, t_hour, t_min, t_sec, t_msec

integer :: ndays, nsecs
type(time_type) :: gitm_time, epoch_time, delta_time

!======================================================================

call initialize_utilities(progname='gitm_to_netcdf')

call set_calendar_type('GREGORIAN')

iunit =  open_file(gitm_to_netcdf_2d_input_file, action='read', form='unformatted')
call read_header(iunit, twod=.true.)

do i=1, nvars_2D
   allocate(gitm_2Dvars(i)%state(lon_size, lat_size))  ! lon/lat in radians, alt in meters
   read(iunit) gitm_2Dvars(i)%state
   print *, gitm_2Dvars(i)%state(1,1)
enddo

call close_file(iunit)

! create 2D netcdf file
ncid = nc_create_file(gitm_to_netcdf_2d_output_file)
print *, 'ncid = ', ncid

! create netcdf dims
call nc_define_dimension(ncid, "lon", lon_size)
call nc_define_dimension(ncid, "lat", lat_size)

! create netcdf vars
do i=1, nvars_2D
print *, 'going to define: ', i, trim(gitm_2Dvars(i)%varname)
   call nc_define_real_variable(ncid, trim(gitm_2Dvars(i)%varname), (/ "lon", "lat" /) )
enddo

! add attributes

! change into data mode
call nc_end_define_mode(ncid)

! fill in coordinate vars?

! fill netcdf vars
do i=1, nvars_2D
print *, 'going to add data for: ', i, trim(gitm_2Dvars(i)%varname)
   call nc_put_variable(ncid, trim(gitm_2Dvars(i)%varname), gitm_2Dvars(i)%state)
enddo


! close netcdf file
call nc_close_file(ncid)


iunit =  open_file(gitm_to_netcdf_3d_input_file, action='read', form='unformatted')
call read_header(iunit, twod=.false.)

do i=1, nvars_3D
   allocate(gitm_3Dvars(i)%state(lon_size, lat_size, alt_size))  ! lon/lat in radians, alt in meters
   read(iunit) gitm_3Dvars(i)%state
   print *, gitm_3Dvars(i)%state(1,1,1)
enddo

call close_file(iunit)

! create netcdf file
ncid = nc_create_file(gitm_to_netcdf_3d_output_file)
print *, 'ncid = ', ncid

! create netcdf dims
call nc_define_dimension(ncid, "lon", lon_size)
call nc_define_dimension(ncid, "lat", lat_size)
call nc_define_dimension(ncid, "alt", alt_size)

! create netcdf vars
do i=1, nvars_3D
print *, 'going to define: ', i, trim(gitm_3Dvars(i)%varname)
   call nc_define_real_variable(ncid, trim(gitm_3Dvars(i)%varname), (/ "lon", "lat", "alt" /) )
enddo

call nc_define_real_scalar(ncid, 'time')
call nc_add_attribute_to_variable(ncid, 'time', 'calendar', 'GREGORIAN')
call nc_add_attribute_to_variable(ncid, 'time', 'units',    &
                                 'days since 1601-01-01 00:00:00')

! add attributes to other vars

! change into data mode
call nc_end_define_mode(ncid)

! fill in coordinate vars?

! fill netcdf vars
do i=1, nvars_3D
print *, 'going to add data for: ', i, trim(gitm_3Dvars(i)%varname)
   call nc_put_variable(ncid, trim(gitm_3Dvars(i)%varname), gitm_3Dvars(i)%state)
enddo

! add time var and attributes
epoch_time = set_date(1601, 1, 1, 0, 0, 0)
gitm_time = set_date(t_year, t_month, t_day, t_hour, t_min, t_sec)
delta_time = gitm_time - epoch_time

call get_time(delta_time, nsecs, ndays)
real_time = ndays + (real(nsecs, r8) / 86400.0_r8)

call nc_put_variable(ncid, 'time', real_time)

! close netcdf file
call nc_close_file(ncid)



! end - close the log, etc
call finalize_utilities()

contains

subroutine read_header(iunit, twod)
 integer, intent(in) :: iunit
 logical, intent(in) :: twod

integer :: nvars, i, j
character(len=NAMELEN) :: tempvar

read(iunit) gitm_version
print *, 'version = ', gitm_version

read(iunit) lon_size, lat_size, alt_size
print *, 'sizes = ', lon_size, lat_size, alt_size

read(iunit) nvars
print *, 'number of vars = ', nvars

if (twod) then
   nvars_2D = nvars 
else 
   nvars_3D = nvars 
endif

do i=1, nvars
   if (twod) then
      read(iunit) tempvar
      tempvar = adjustl(tempvar)
      do j=1, len_trim(tempvar)
         if (tempvar(j:j) == ' ') tempvar(j:j) = '_'
      enddo
      gitm_2Dvars(i)%varname = tempvar
      print *, 'var num, name = ', i, trim(gitm_2Dvars(i)%varname)
   else
      read(iunit) tempvar
      tempvar = adjustl(tempvar)
      do j=1, len_trim(tempvar)
         if (tempvar(j:j) == ' ') tempvar(j:j) = '_'
      enddo
      gitm_3Dvars(i)%varname = tempvar
      print *, 'var num, name = ', i, trim(gitm_3Dvars(i)%varname)
   endif
enddo

read(iunit) t_year, t_month, t_day, t_hour, t_min, t_sec, t_msec
print *, 'time: ', t_year, t_month, t_day, t_hour, t_min, t_sec, t_msec

end subroutine read_header


end program gitm_to_netcdf


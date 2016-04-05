! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod_check.f90 6739 2014-01-15 20:44:54Z hkershaw $

program model_mod_check

!----------------------------------------------------------------------
! purpose: test routines.  this version for models with oned locations.
!----------------------------------------------------------------------

use             types_mod, only : r8, i8, missing_r8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  nc_check, E_MSG, open_file, close_file

use     mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities

use          location_mod, only : location_type, set_location, write_location,  &
                                  get_dist, get_location

use          obs_kind_mod, only : get_obs_kind_name, get_raw_obs_kind_index,    &
                                  KIND_DRY_LAND

use      obs_sequence_mod, only : static_init_obs_sequence

use       assim_model_mod, only : aoutput_diagnostics, init_diag_output,        &
                                  finalize_diag_output, static_init_assim_model,&
                                  netcdf_file_type

use      time_manager_mod, only : time_type, set_calendar_type, GREGORIAN,      &
                                  set_time, print_time, print_date, operator(-)

use  ensemble_manager_mod, only : init_ensemble_manager, ensemble_type

use   state_vector_io_mod, only : state_vector_io_init, setup_read_write,   &
                                  turn_read_copy_on, turn_write_copy_on,    &
                                  filter_read_restart_direct,               &
                                  filter_write_restart_direct

use   state_structure_mod, only : static_init_state_type

use            filter_mod, only : filter_set_initial_time

use      io_filenames_mod, only : io_filenames_init, set_filenames,         &
                                  get_input_file, get_output_file

use             model_mod, only : static_init_model, get_model_size, &
                                  get_state_meta_data,       &
                                  model_interpolate

use  distributed_state_mod, only : create_state_window, free_state_window, &
                                   create_mean_window, free_mean_window

use  test_interpolate_mod, only : test_interpolate_single, test_interpolate_range

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_model_mod_check/models/template/model_mod_check.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 6739 $"
character(len=128), parameter :: revdate  = "$Date: 2014-01-15 13:44:54 -0700 (Wed, 15 Jan 2014) $"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

integer(i8)                   :: x_ind   = -1
integer                       :: num_ens = 1
real(r8), dimension(3)        :: loc_of_interest = -1.0_r8
character(len=metadatalength) :: kind_of_interest = 'ANY'
character(len=metadatalength) :: interp_test_vertcoord = 'VERTISHEIGHT'
logical                       :: verbose = .FALSE.
integer                       :: test1thru = 1
real(r8)               :: interp_test_dlat  = 10.0
real(r8)               :: interp_test_dlon  = 10.0
real(r8)               :: interp_test_dvert = 10.0
real(r8), dimension(2) :: interp_test_latrange  = (/   0.0,  120.0 /)
real(r8), dimension(2) :: interp_test_lonrange  = (/   0.0,  120.0 /)
real(r8), dimension(2) :: interp_test_vertrange = (/   0.0,  100.0 /)
real(r8)               :: interp_test_dx = missing_r8
real(r8)               :: interp_test_dy = missing_r8
real(r8)               :: interp_test_dz = missing_r8
real(r8), dimension(2) :: interp_test_xrange = (/ missing_r8, missing_r8 /)
real(r8), dimension(2) :: interp_test_yrange = (/ missing_r8, missing_r8 /)
real(r8), dimension(2) :: interp_test_zrange = (/ missing_r8, missing_r8 /)

namelist /model_mod_check_nml/ x_ind, num_ens,                         &
                               loc_of_interest, kind_of_interest,      &
                               interp_test_dlat, interp_test_lonrange, &
                               interp_test_dlon, interp_test_latrange, &
                               interp_test_dvert, interp_test_vertrange, &
                               interp_test_dx, interp_test_xrange,     &
                               interp_test_dy, interp_test_yrange,     &
                               interp_test_dz, interp_test_zrange,     &
                               interp_test_vertcoord,                  &
                               verbose, test1thru

! io variables
integer :: iunit, io
integer, allocatable :: ios_out(:)

! model state variables
type(ensemble_type) :: ens_handle

type(time_type)       :: time1, model_time
integer               :: mykindindex
integer(i8)           :: model_size
real(r8), allocatable :: interp_vals(:)

! misc. variables
integer :: imem, num_failed
logical :: cartesian = .false.

! error handler strings
character(len=512) :: string1

!----------------------------------------------------------------------
! This portion checks the geometry information.
!----------------------------------------------------------------------

call initialize_modules_used()

call find_namelist_in_file("input.nml", "model_mod_check_nml", iunit)
read(iunit, nml = model_mod_check_nml, iostat = io)
call check_namelist_read(iunit, io, "model_mod_check_nml")

if ( interp_test_dx  /= missing_r8 .or. &
     interp_test_dy  /= missing_r8 .or. &
     interp_test_dz  /= missing_r8 ) then

   ! if the user defines cartesian coordinates just 
   ! overwrite values for the test_interpolation calls.

   interp_test_dlon  = interp_test_dx
   interp_test_dlat  = interp_test_dy
   interp_test_dvert = interp_test_dz
   
   interp_test_lonrange  = interp_test_xrange 
   interp_test_latrange  = interp_test_yrange 
   interp_test_vertrange = interp_test_zrange 
   
   cartesian = .true.
endif

write(*,*)
write(*,*) '******************** RUNNING TEST 1 *********************'
write(*,*) ' Reading the namelist and running static_init_model       '
write(*,*) '  -- calling get_model_size()'
write(*,*) '*********************************************************'
write(*,*)

call static_init_assim_model()

model_size = get_model_size()

write(*,*)
write(*,'(''state vector has length'',i10)') model_size
write(*,*)

write(*,*) '******************** FINISHED TEST 1 *********************'


if ( test1thru == 1 ) call exit(0)


write(*,*)
write(*,*) '******************** RUNNING TEST 2 *********************'
write(*,*) ' Read and write trivial restart file.                    '
write(*,*) '*********************************************************'
write(*,*)

call set_calendar_type(GREGORIAN)

model_time  = set_time(21600, 149446)   ! 06Z 4 March 2010

! Set up the ensemble storage and read in the restart file
call init_ensemble_manager(ens_handle, num_ens, model_size)

! Reading netcdf restart file:
call setup_read_write(num_ens)
call set_filenames(ens_handle, num_ens, "no_inf1", "no_inf2")

!----------------------------------------------------------------------
! Open a test netcdf initial conditions file.
!----------------------------------------------------------------------

do imem = 1, num_ens
   call turn_read_copy_on(imem)
   write(*,*) 'reading  ', trim( get_input_file(imem, domain=1) )
enddo

call filter_read_restart_direct(ens_handle, time1, num_extras=0, use_time_from_file=.true.)
model_time = time1

do imem = 1, num_ens
   call turn_write_copy_on(imem)
   write(*,*) 'writting ', trim( get_output_file(imem, domain=1, isprior=.true.) )
enddo
call filter_write_restart_direct(ens_handle, num_extras=0, isprior=.true.)

write(*,*) 
call print_date( model_time,' model_mod_check:model date')
call print_time( model_time,' model_mod_check:model time')

write(*,*)
write(*,*) '******************** FINISHED TEST 2 *********************'
write(*,*)

if ( test1thru == 2 ) call exit(0)

!----------------------------------------------------------------------
! Check the meta data
!----------------------------------------------------------------------

write(*,*) '******************** RUNNING TEST 3 *********************'
write(*,*) ' Testing get_state_meta_data                             '
write(*,*) '*********************************************************'

if ( x_ind > 0 .and. x_ind <= model_size ) then
   call check_meta_data( x_ind )
else
   write(*,*) "x_ind = ", x_ind, " not in valid range of model 0-", model_size
endif

write(*,*)
write(*,*) '******************** FINISHED TEST 3 *********************'
write(*,*)

if ( test1thru == 3 ) call exit(0)

!----------------------------------------------------------------------
! Check the interpolation - print initially to STDOUT
!----------------------------------------------------------------------

write(*,*)
write(*,*) '******************** RUNNING TEST 4 *********************'
write(*,*) ' Testing loc_of_interest for model_interpolate           '
write(*,*) '*********************************************************'
write(*,*)

call create_state_window(ens_handle)

mykindindex = get_raw_obs_kind_index(kind_of_interest)

allocate(interp_vals(num_ens), ios_out(num_ens))

write(*,*) "interpolating at ", loc_of_interest

num_failed = test_interpolate_single( ens_handle,            &
                                      num_ens,               &
                                      interp_test_vertcoord, &
                                      loc_of_interest(1),    &
                                      loc_of_interest(2),    &
                                      loc_of_interest(3),    &
                                      mykindindex,           &
                                      interp_vals,           &
                                      ios_out )

write(*,*)
write(*,*) '******************** FINISHED TEST 4 *********************'
write(*,*)

if ( test1thru == 4 ) call exit(0)

write(*,*)
write(*,*) '******************** RUNNING TEST 5 *********************'
write(*,*) ' Testing range of data for model_interpolate             '
write(*,*) '*********************************************************'
write(*,*)

num_failed = test_interpolate_range( ens_handle,            &
                                     num_ens,               &
                                     interp_test_dlon,      &
                                     interp_test_dlat,      &
                                     interp_test_dvert,     &
                                     interp_test_vertcoord, &
                                     interp_test_lonrange,  &
                                     interp_test_latrange,  &
                                     interp_test_vertrange, &
                                     mykindindex,           &
                                     verbose )

write(*,*)
write(*,*) '******************** FINISHED TEST 5 *********************'
write(*,*)

! finalize model_mod_check
call error_handler(E_MSG,'full_model_mod_check','Finished successfully.',source,revision,revdate)
call finalize_mpi_utilities()

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------
subroutine check_meta_data( iloc )

integer(i8), intent(in) :: iloc
type(location_type) :: loc
integer             :: var_type

write(*,*)
write(*,*)'Checking metadata routines.'
write(*,*)

call create_mean_window(ens_handle, 1, .false.)
call get_state_meta_data(ens_handle, iloc, loc, var_type)
call free_mean_window

call write_location(42, loc, fform='formatted', charstring=string1)
write(*,*)
write(*,*)' indx ',iloc,' is type ',var_type
write(*,*)'   ', trim(string1)

end subroutine check_meta_data

!----------------------------------------------------------------------

subroutine initialize_modules_used()

! Standard initialization (mpi not needed to use ensemble manager
! since we are enforcing that this run as a single task).
call initialize_mpi_utilities('model_mod_check')

! Initialize modules used that require it
call register_module(source,revision,revdate)

! Initialize modules used that require it
call static_init_obs_sequence()

call static_init_state_type()
call state_vector_io_init()
call io_filenames_init()

end subroutine initialize_modules_used

!----------------------------------------------------------------------


end program model_mod_check

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/CM1/models/cm1/model_mod_check.f90 $
! $Id: model_mod_check.f90 8729 2015-10-01 22:31:54Z hendric $
! $Revision: 8729 $
! $Date: 2015-10-01 16:31:54 -0600 (Thu, 01 Oct 2015) $

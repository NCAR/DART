! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!----------------------------------------------------------------------
!> purpose: test model_mod routines.  this version works for models
!> with any location type.  depends on a location-specific module
!> for test_interpolate_single and test_interpolate_range.
!----------------------------------------------------------------------

program model_mod_check

use             types_mod, only : r8, i8, missing_r8, metadatalength, MAX_NUM_DOMS

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  nc_check, E_MSG, open_file, close_file, do_output

use     mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities

use          location_mod, only : location_type, write_location

use          obs_kind_mod, only : get_index_for_quantity, get_name_for_quantity

use      obs_sequence_mod, only : static_init_obs_sequence

use       assim_model_mod, only : static_init_assim_model

use      time_manager_mod, only : time_type, set_time, print_time, print_date, operator(-)

use  ensemble_manager_mod, only : init_ensemble_manager, ensemble_type

use   state_vector_io_mod, only : state_vector_io_init, read_state, write_state

use   state_structure_mod, only : get_num_domains

use      io_filenames_mod, only : io_filenames_init, file_info_type,       &
                                  stage_metadata_type, get_stage_metadata, &
                                  get_restart_filename,                    &
                                  set_file_metadata, file_info_dump,       &
                                  set_io_copy_flag, READ_COPY, WRITE_COPY

use distributed_state_mod, only : create_state_window, free_state_window

use             model_mod, only : static_init_model, get_model_size,       &
                                  get_state_meta_data, model_interpolate

use  test_interpolate_mod, only : test_interpolate_single, test_interpolate_range

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

logical                       :: single_file = .false.
integer                       :: num_ens = 1
character(len=256)            :: input_state_files(MAX_NUM_DOMS)  = 'null'
character(len=256)            :: output_state_files(MAX_NUM_DOMS) = 'null'

integer(i8)                   :: x_ind   = -1
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

namelist /model_mod_check_nml/ x_ind, num_ens,                             &
                               loc_of_interest,    kind_of_interest,       &
                               interp_test_dlat,   interp_test_lonrange,   &
                               interp_test_dlon,   interp_test_latrange,   &
                               interp_test_dvert,  interp_test_vertrange,  &
                               interp_test_dx,     interp_test_xrange,     &
                               interp_test_dy,     interp_test_yrange,     &
                               interp_test_dz,     interp_test_zrange,     &
                               verbose, test1thru, interp_test_vertcoord,  &
                               single_file, input_state_files, output_state_files

! io variables
integer                   :: iunit, io
integer, allocatable      :: ios_out(:)
type(file_info_type)      :: file_info_input, file_info_output
type(stage_metadata_type) :: input_restart_files, output_restart_files
logical :: read_time_from_file = .true.

! model state variables
type(ensemble_type)   :: ens_handle

type(time_type)       :: model_time
integer               :: mykindindex
integer(i8)           :: model_size
real(r8), allocatable :: interp_vals(:)

! misc. variables
integer :: idom, imem, num_failed, num_domains
logical :: cartesian = .false.

! error handler strings
character(len=512)              :: my_base, my_desc, my_location
character(len=256), allocatable :: file_array_input(:,:), file_array_output(:,:)

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

call print_test_message('RUNNING TEST 1', &
                        'Reading the namelist and running static_init_model', &
                        'calling get_model_size()')

call static_init_assim_model()

model_size = get_model_size()

if ( do_output() ) then 
   write(*,'(A)') '-------------------------------------------------------------'
   write(*,'(''state vector has length of '',i10)') model_size
   write(*,'(A)') '-------------------------------------------------------------'
endif

call print_test_message('FINISHED TEST 1')
write(*,'(A)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
write(*,'(A)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

if ( test1thru == 1 ) call exit(0)


call print_test_message('RUNNING TEST 2', &
                        'Read and write trivial restart file')

! Set up the ensemble storage and read in the restart file
call init_ensemble_manager(ens_handle, num_ens, model_size)

! Allocate space for file arrays.  contains a matrix of files (num_ens x num_domains)
! If perturbing from a single instance the number of input files does not have to
! be ens_size but rather a single file (or multiple files if more than one domain)

num_domains = get_num_domains()
allocate(file_array_input(num_ens, num_domains), file_array_output(num_ens, num_domains))

file_array_input  = RESHAPE(input_state_files,  (/num_ens,  num_domains/))
file_array_output = RESHAPE(output_state_files, (/num_ens,  num_domains/))

! Initialize input file info
call io_filenames_init(file_info_input,             &
                       ncopies      = num_ens,      &
                       cycling      = single_file,  &
                       single_file  = single_file,  & 
                       restart_files = file_array_input)

do imem = 1, num_ens
   write(my_base,'(A,I2)') 'inens_',    imem
   write(my_desc,'(A,I2)') 'input ens', imem
   call set_file_metadata(file_info_input,                          &
                          cnum     = imem,                          &
                          fnames   = (/input_state_files(imem)/),  &
                          basename = my_base,                       &
                          desc     = my_desc)
   
   call set_io_copy_flag(file_info_input, &
                         cnum    = imem,     &
                         io_flag = READ_COPY)
enddo

! Initialize output file info
call io_filenames_init(file_info_output,           &
                       ncopies      = num_ens,     &
                       cycling      = single_file, &
                       single_file  = single_file, &
                       restart_files = file_array_output)
      
do imem = 1, num_ens
   write(my_base,'(A,I2)') 'outens_',    imem
   write(my_desc,'(A,I2)') 'output ens', imem
   call set_file_metadata(file_info_output,                          & 
                          cnum     = imem,                           &
                          fnames   = (/output_state_files(imem)/),  &
                          basename = my_base,                        &
                          desc     = my_desc)
   
   call set_io_copy_flag(file_info_output,    &
                         cnum    = imem,      &
                         io_flag = WRITE_COPY)  
enddo

!call file_info_dump(file_info_input, 'mmc')

!----------------------------------------------------------------------
! Open a test netcdf initial conditions file.
!----------------------------------------------------------------------
input_restart_files = get_stage_metadata(file_info_input)

do idom = 1, num_domains
do imem = 1, num_ens
   if ( do_output() ) then
      write(*,'(A)')  ''
      write(*,'(A)')  '-------------------------------------------------------------'
      write(*,'(2A)') '- Reading File : ', trim(get_restart_filename(input_restart_files, imem, domain=idom))
      write(*,'(A)')  '-------------------------------------------------------------'
   endif
enddo
enddo

call read_state(ens_handle, file_info_input, read_time_from_file, model_time)

output_restart_files = get_stage_metadata(file_info_output)
do idom = 1, num_domains
do imem = 1, num_ens
   if ( do_output() ) then
      write(*,'(A)')  ''
      write(*,'(A)')  '-------------------------------------------------------------'
      write(*,'(2A)') '- Writing File : ', trim(get_restart_filename(output_restart_files, imem, domain=idom))
      write(*,'(A)')  '-------------------------------------------------------------'
   endif
enddo
enddo

call write_state(ens_handle, file_info_output)

! print date does not work when a model does not have a calendar
!write(*,'(A)') '-- printing model date --------------------------------------'
!call print_date( model_date,' model_mod_check:model date')
write(*,'(A)') '-- printing model time --------------------------------------'
call print_time( model_time,' model_mod_check:model time')
write(*,'(A)') '-------------------------------------------------------------'
write(*,'(A)') ''

call print_test_message('FINISHED TEST 2')
write(*,'(A)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
write(*,'(A)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

if ( test1thru == 2 ) call exit(0)

!----------------------------------------------------------------------
! Check the meta data
!----------------------------------------------------------------------

call print_test_message('RUNNING TEST 3', &
                        'Testing get_state_meta_data')

if ( x_ind >= 1 .and. x_ind <= model_size ) then
   call check_meta_data( x_ind )
else
   if ( do_output() )  then
      write(*,'(A)') '-------------------------------------------------------------'
      write(*,'(A)') "x_ind = ", x_ind, " is not in valid range 1-", model_size
      write(*,'(A)') '-------------------------------------------------------------'
   endif
endif

call print_test_message('FINISHED TEST 3')
write(*,'(A)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
write(*,'(A)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

if ( test1thru == 3 ) call exit(0)

!----------------------------------------------------------------------
! Check the interpolation - print initially to STDOUT
!----------------------------------------------------------------------

call print_test_message('RUNNING TEST 4', &
                        'Testing loc_of_interest for model_interpolate')

call create_state_window(ens_handle)

mykindindex = get_index_for_quantity(kind_of_interest)

allocate(interp_vals(num_ens), ios_out(num_ens))

num_failed = test_interpolate_single( ens_handle,            &
                                      num_ens,               &
                                      interp_test_vertcoord, &
                                      loc_of_interest(1),    &
                                      loc_of_interest(2),    &
                                      loc_of_interest(3),    &
                                      mykindindex,           &
                                      interp_vals,           &
                                      ios_out )

call print_test_message('FINISHED TEST 4')
write(*,'(A)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
write(*,'(A)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

if ( test1thru == 4 ) call exit(0)

call print_test_message('RUNNING TEST 5', &
                        'Testing range of data for model_interpolate')
write(*,'(A)') ''

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

call print_test_message('FINISHED TEST 5')
write(*,'(A)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
write(*,'(A)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

! finalize model_mod_check

write(*,'(A)') ''
write(*,'(A)') '-------------------------------------------------------------'
write(*,'(A)') '- model_mod_check Finished successfully                      '
write(*,'(A)') '-------------------------------------------------------------'

call finalize_mpi_utilities()

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------
subroutine check_meta_data( iloc )

integer(i8), intent(in) :: iloc

type(location_type) :: loc
integer             :: var_type

call get_state_meta_data(iloc, loc, var_type)

if ( do_output() ) then
   write(*,'(A)')  ''
   write(*,'(A)') '-------------------------------------------------------------'
   write(*,'("index = ",I10,", has variable type : ", I4," ", A)') &
      iloc, var_type, trim(get_name_for_quantity(var_type))
   write(*,'(A)') '-------------------------------------------------------------'
endif 

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

call state_vector_io_init()

end subroutine initialize_modules_used

!----------------------------------------------------------------------

subroutine print_test_message(test_msg, msg1, msg2, msg3)

character(len=*), intent(in) :: test_msg
character(len=*), intent(in), optional :: msg1
character(len=*), intent(in), optional :: msg2
character(len=*), intent(in), optional :: msg3

character(len=15) :: my_msg
character(len=64) :: msg_string
character(len=64) :: msg_close

my_msg = test_msg

if ( do_output() ) then
   write(*,'(A)') ''
   write(*,'(A)') ''
   write(msg_string,'(3A)') '******************** ', my_msg,    ' *************************'
   write(msg_close ,'(A)' ) '**************************************************************'

                        write(*,'(A)') trim(msg_string)
   if ( present(msg1) ) write(*,'(2A)') ' -- ', trim(msg1)
   if ( present(msg2) ) write(*,'(2A)') ' -- ', trim(msg2)
   if ( present(msg3) ) write(*,'(2A)') ' -- ', trim(msg3)
   if ( present(msg1) ) write(*,'(A)') trim(msg_close)

endif

end subroutine print_test_message


!----------------------------------------------------------------------


end program model_mod_check

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

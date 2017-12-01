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

use      time_manager_mod, only : time_type, set_time, print_time, print_date, operator(-), &
                                  get_calendar_type, NO_CALENDAR

use  ensemble_manager_mod, only : init_ensemble_manager, ensemble_type

use   state_vector_io_mod, only : state_vector_io_init, read_state, write_state

use   state_structure_mod, only : get_num_domains, get_model_variable_indices

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

integer, parameter :: MAX_TESTS = 15

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

logical                       :: single_file = .false.
integer                       :: num_ens = 1
character(len=256)            :: input_state_files(MAX_NUM_DOMS)  = 'null'
character(len=256)            :: output_state_files(MAX_NUM_DOMS) = 'null'
character(len=256)            :: all_metadata_file = 'metadata.txt'
integer(i8)                   :: x_ind   = -1
real(r8), dimension(3)        :: loc_of_interest = -1.0_r8
character(len=metadatalength) :: kind_of_interest = 'ANY'
character(len=metadatalength) :: interp_test_vertcoord = 'VERTISHEIGHT'
logical                       :: verbose = .FALSE.
integer                       :: test1thru = MAX_TESTS
integer                       :: run_tests(MAX_TESTS) = -1
real(r8)               :: interp_test_dlat  = 10.0_r8
real(r8)               :: interp_test_dlon  = 10.0_r8
real(r8)               :: interp_test_dvert = 10.0_r8
real(r8), dimension(2) :: interp_test_latrange  = (/   0.0_r8,  120.0_r8 /)
real(r8), dimension(2) :: interp_test_lonrange  = (/   0.0_r8,  120.0_r8 /)
real(r8), dimension(2) :: interp_test_vertrange = (/   0.0_r8,  100.0_r8 /)
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
                               verbose, test1thru, run_tests, interp_test_vertcoord,  &
                               single_file, input_state_files, output_state_files, &
                               all_metadata_file

! io variables
integer                   :: iunit, io
integer, allocatable      :: ios_out(:)
type(file_info_type)      :: file_info_input, file_info_output
type(stage_metadata_type) :: input_restart_files, output_restart_files
logical :: read_time_from_file = .true.
logical :: tests_to_run(MAX_TESTS) = .false.

! model state variables
type(ensemble_type)   :: ens_handle

type(time_type)       :: model_time
integer               :: mykindindex
integer(i8)           :: model_size
real(r8), allocatable :: interp_vals(:)

! misc. variables
integer :: idom, imem, num_passed, num_failed, num_domains
logical :: cartesian = .false.

! message strings
character(len=512) :: my_base, my_desc
character(len=512) :: string1, string2, string3

character(len=256), allocatable :: file_array_input(:,:), file_array_output(:,:)

!======================================================================
! start of executable code
!======================================================================

call initialize_modules_used()

call find_namelist_in_file("input.nml", "model_mod_check_nml", iunit)
read(iunit, nml = model_mod_check_nml, iostat = io)
call check_namelist_read(iunit, io, "model_mod_check_nml")

call setup_run_array()
call setup_interp_grid()
call static_init_assim_model() ! must be done for all tests

!----------------------------------------------------------------------
! initialization code, model size
!----------------------------------------------------------------------

if (tests_to_run(1)) then

   call print_test_message('TEST 1', &
                           'Reading the namelist and running static_init_model', &
                           'calling get_model_size()', starting=.true.)

   model_size = get_model_size()

   write(string1, '(A,I10)') 'state vector has length of ', model_size
   call print_info_message('TEST 1', string1)

   call print_test_message('TEST 1', ending=.true.)

endif

!----------------------------------------------------------------------
! read/write restart files
!----------------------------------------------------------------------

if (tests_to_run(2)) then

   call print_test_message('TEST 2', &
                           'Read and write restart file', starting=.true.)

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
                             fnames   = (/file_array_input(imem,:)/),  &
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
                             fnames   = (/file_array_output(imem,:)/),  &
                             basename = my_base,                        &
                             desc     = my_desc)

      call set_io_copy_flag(file_info_output,    &
                            cnum    = imem,      &
                            io_flag = WRITE_COPY)
   enddo

   ! Open a test netcdf initial conditions file.
   input_restart_files = get_stage_metadata(file_info_input)

   do idom = 1, num_domains
      do imem = 1, num_ens
         write(string1, *) '- Reading File : ', trim(get_restart_filename(input_restart_files, imem, domain=idom))
         call print_info_message('TEST 2',string1)
      enddo
   enddo

   call read_state(ens_handle, file_info_input, read_time_from_file, model_time)

   output_restart_files = get_stage_metadata(file_info_output)
   do idom = 1, num_domains
      do imem = 1, num_ens
         write(string1, *) '- Writing File : ', trim(get_restart_filename(output_restart_files, imem, domain=idom))
         call print_info_message('TEST 2',string1)
      enddo
   enddo

   call write_state(ens_handle, file_info_output)

   ! print date does not work when a model does not have a calendar
   if (get_calendar_type() /= NO_CALENDAR) then
      write(*,'(A)') '-- printing model date --------------------------------------'
      call print_date( model_time,' model_mod_check:model date')
   endif

   write(*,'(A)') '-- printing model time --------------------------------------'
   call print_time( model_time,' model_mod_check:model time')
   write(*,'(A)') '-------------------------------------------------------------'
   write(*,'(A)') ''

   call print_test_message('TEST 2', ending=.true.)

endif

!----------------------------------------------------------------------
! Check the metadata
!----------------------------------------------------------------------

if (tests_to_run(3)) then

   call print_test_message('TEST 3', &
                           'Testing get_state_meta_data', starting=.true.)

   if ( x_ind >= 1 .and. x_ind <= model_size ) then
      call check_meta_data( x_ind )
   else
      write(string1, *) "x_ind = ", x_ind, " is not in valid range 1-", model_size
      call print_info_message('TEST 3',string1)
   endif

   call print_test_message('TEST 3', ending=.true.)

endif

!----------------------------------------------------------------------
! Check the interpolation - interpolate a single point
!----------------------------------------------------------------------

!>@todo requires ensemble from test 2

if (tests_to_run(4)) then
   call print_test_message('TEST 4', &
                           'Testing loc_of_interest for model_interpolate', starting=.true.)

   call create_state_window(ens_handle)

   mykindindex = get_index_for_quantity(kind_of_interest)

   allocate(interp_vals(num_ens), ios_out(num_ens))

   num_passed = test_interpolate_single( ens_handle,            &
                                         num_ens,               &
                                         interp_test_vertcoord, &
                                         loc_of_interest(1),    &
                                         loc_of_interest(2),    &
                                         loc_of_interest(3),    &
                                         mykindindex,           &
                                         interp_vals,           &
                                         ios_out )

   ! test_interpolate_single reports individual interpolation failures internally
   if (num_passed == num_ens) then
      write(string1,*)'interpolation successful for all ensemble members.'
      call print_info_message('TEST 4',string1)
   endif

   call print_test_message('TEST 4', ending=.true.)
endif

!----------------------------------------------------------------------
! Check the interpolation with a test grid
!----------------------------------------------------------------------

if (tests_to_run(5)) then
   call print_test_message('TEST 5', &
                           'Testing range of data for model_interpolate', starting=.true.)

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

   ! test_interpolate_range internally reports interpolation metrics.
   write(string1, *)'output values on interpolation grid are in'
   write(string2, *)'check_me_interptest.nc (netcdf) and check_me_interptest.m (matlab)'
   call print_info_message(string1, string2)

   call print_test_message('TEST 5', ending=.true.)
endif

!----------------------------------------------------------------------
! Exhaustive test - print the metadata for every element in the state. 
!----------------------------------------------------------------------

if (tests_to_run(6)) then

   write(string1,*)'Exhaustive test of get_state_meta_data - please be patient.'
   write(string2,*)'There are ',get_model_size(),' items in the state vector.'
   call print_test_message('TEST 6', string1, string2, starting=.true.)

   call check_all_meta_data()

   call print_test_message('TEST 6', ending=.true.)
endif

!----------------------------------------------------------------------
! add more tests here
!----------------------------------------------------------------------

!>@todo possibly add a check to find the i,j,k of the gridcell closest
!> to the location of interest.

! finalize model_mod_check

write(string1,*) '- model_mod_check Finished successfully'
call print_info_message(string1)

call finalize_mpi_utilities()

!======================================================================
contains
!======================================================================

!> initialize modules that need it

subroutine initialize_modules_used()

call initialize_mpi_utilities('model_mod_check')

call register_module(source,revision,revdate)

call static_init_obs_sequence()

call state_vector_io_init()

end subroutine initialize_modules_used

!----------------------------------------------------------------------
!> print the results of get_state_meta_data() at a single location

subroutine check_meta_data( iloc )

integer(i8), intent(in) :: iloc

type(location_type) :: loc
integer             :: ix, iy, iz, dom_id, qty_index, var_type
character(len=256)  :: qty_string

call get_state_meta_data(iloc, loc, var_type)
call get_model_variable_indices(iloc, ix, iy, iz, &
                                   dom_id=dom_id, &
                                   kind_index=qty_index, &
                                   kind_string=qty_string)

write(string1,'(i11,'' is i,j,k'',3(1x,i4),'' and is in domain '',i2)') &
                  iloc, ix, iy, iz, dom_id
write(string2,'(''is quantity '', I4,'', '',A)') var_type, trim(qty_string)
call write_location(0,loc,charstring=string3)

call print_info_message(string1, string2, string3)

end subroutine check_meta_data

!----------------------------------------------------------------------
!> print the labels between the starts of tests

subroutine print_test_message(test_label, msg1, msg2, msg3, starting, ending)

character(len=*), intent(in) :: test_label
character(len=*), intent(in), optional :: msg1
character(len=*), intent(in), optional :: msg2
character(len=*), intent(in), optional :: msg3
logical,          intent(in), optional :: starting
logical,          intent(in), optional :: ending

if (do_output()) &
   call print_message(.true., test_label, msg1, msg2, msg3, starting, ending)

end subroutine print_test_message

!----------------------------------------------------------------------
!> print an informational message

subroutine print_info_message(info_msg, msg1, msg2, msg3, starting, ending)

character(len=*), intent(in) :: info_msg
character(len=*), intent(in), optional :: msg1
character(len=*), intent(in), optional :: msg2
character(len=*), intent(in), optional :: msg3
logical,          intent(in), optional :: starting
logical,          intent(in), optional :: ending

if (do_output()) &
   call print_message(.false., info_msg, msg1, msg2, msg3, starting, ending)

end subroutine print_info_message

!----------------------------------------------------------------------
!> common code for printing

subroutine print_message(is_test_label, msg, msg1, msg2, msg3, starting, ending)

logical,          intent(in) :: is_test_label   ! true is test, false is info
character(len=*), intent(in) :: msg
character(len=*), intent(in), optional :: msg1
character(len=*), intent(in), optional :: msg2
character(len=*), intent(in), optional :: msg3
logical,          intent(in), optional :: starting
logical,          intent(in), optional :: ending

character(len=20) :: test_label
character(len=64) :: msg_label
character(len=64) :: msg_blank
character(len=64) :: msg_close
character(len=64) :: msg_sep1
character(len=64) :: msg_sep2
logical :: is_start, is_end

! if my task isn't writing output, return now.
if (.not. do_output()) return

! setup section - set defaults for optional arguments
! so we don't have to keep testing them.

is_start = set_logical_flag(.false., starting)
is_end   = set_logical_flag(.false., ending)

! is it documented anywhere that this can only be 20 chars long?
! i'm assuming this was set up so the separators would line up.

if (is_test_label) then
   if (is_start) then
      test_label = 'RUNNING    '//trim(msg)
   else if (is_end) then
      test_label = 'FINISHED   '//trim(msg)
   else
      test_label = msg
   endif
   write(msg_label, '(3A)') '***************** ', test_label,     ' ***********************'
endif

write(msg_close ,'(A)' ) '**************************************************************'
write(msg_sep1,  '(A)' ) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
write(msg_sep2,  '(A)' ) '--------------------------------------------------------------'
write(msg_blank, '(A)' ) ''

! ok, here's where the actual writing happens.
! if you want to change the formatting, fool around with
! the order and formatting of these lines and it will affect
! all the output from this program.

if (is_test_label) then
                      write(*,'(A)') trim(msg_blank)
                      write(*,'(A)') trim(msg_blank)

                      write(*,'(A)') trim(msg_label)
   if (present(msg1)) write(*,'(2A)') ' -- ', trim(msg1)
   if (present(msg2)) write(*,'(2A)') ' -- ', trim(msg2)
   if (present(msg3)) write(*,'(2A)') ' -- ', trim(msg3)
   if (present(msg1)) write(*,'(A)') trim(msg_close)

   if (is_end) then
                      write(*,'(A)') trim(msg_sep1)
                      write(*,'(A)') trim(msg_sep1)
   endif
else  ! info message
                      write(*,'(A)') trim(msg_sep2)
                      write(*,'(A)') trim(msg)
   if (present(msg1)) write(*,'(2A)') ' -- ', trim(msg1)
   if (present(msg2)) write(*,'(2A)') ' -- ', trim(msg2)
   if (present(msg3)) write(*,'(2A)') ' -- ', trim(msg3)
                      write(*,'(A)') trim(msg_sep2)
endif

end subroutine print_message

!----------------------------------------------------------------------
!> compute the points to be used when testing interpolation

subroutine setup_interp_grid()

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

end subroutine setup_interp_grid

!----------------------------------------------------------------------
!> configure the tests to run

!> if they used the existing 'test1thru', use that to select tests 1-N.
!> otherwise:
!> run_tests(:) is initialized to all -1.  if the user didn't set
!> test1thru and if run_tests() is still all -1, turn everything on.
!> otherwise if the first entry isn't -1, they did specify
!> something in that namelist entry and we turn on just those tests.

subroutine setup_run_array()

integer :: i

tests_to_run(:) = .false.

! be backwards compatible - set this to -1 to disable.
if (test1thru > 0) then
   if (test1thru > MAX_TESTS) then
      write(string1, *) 'test1thru must be between 1 and ', MAX_TESTS, '; found value ', test1thru
      call error_handler(E_ERR, 'model_mod_check: setup_run_array', string1, &
                         source, revision, revdate)
   endif

   tests_to_run(1:test1thru) = .true.
   return
endif

! Enable specific tests.
! the first -1 ends the list.

do i=1, MAX_TESTS
   if (run_tests(i) <= 0) exit

   if (run_tests(i) > MAX_TESTS) then
      write(string1, *) 'test numbers must be between 1 and ', MAX_TESTS, '; found value ', run_tests(i)
      call error_handler(E_ERR, 'model_mod_check: setup_run_array', string1, &
                         source, revision, revdate)
   endif

   tests_to_run(run_tests(i)) = .true.
enddo

! Make sure they are running something
if (run_tests(1) == -1) then
      write(string1, *) 'No tests selected from the namelist.'
      write(string2, *) 'Either specify "test1thru" to be a positive number - or -'
      write(string3, *) 'specify a list of tests to run in "run_tests".'
      call error_handler(E_ERR, 'model_mod_check: setup_run_array', string1, &
                         source, revision, revdate, text2=string2, text3=string3)
endif

!>@todo enforce or report on unfulfilled dependencies

if ((tests_to_run(4) .or. tests_to_run(5)) .and. .not. tests_to_run(2)) then
   write(string1, *) 'The interpolation tests (Test 4, Test 5) need a model state,'
   write(string2, *) 'so Test 2 must be run.'
   call error_handler(E_MSG, 'model_mod_check: setup_run_array', string1, &
                         text2=string2)
   tests_to_run(2) = .true.
endif

end subroutine setup_run_array

!------------------------------------------------------------------
!> given a default value and an optional value,
!> if the optional value is present return that.
!> otherwise return the default value.

!>@todo this belongs in the utils module, yes?

function set_logical_flag(def_val, user_val)

logical, intent(in)           :: def_val
logical, intent(in), optional :: user_val
logical                       :: set_logical_flag

if (present(user_val)) then
   set_logical_flag = user_val
else
   set_logical_flag = def_val
endif

end function set_logical_flag

!------------------------------------------------------------------
!> report the metadata for every element in the state. This can be
!> very slow but is sometimes useful to determine the exact location
!> for specific indices and to confirm the packing order. 

subroutine check_all_meta_data()

integer(i8)         :: iloc
type(location_type) :: loc
integer             :: ix, iy, iz, dom_id, qty_index, var_type, fid
character(len=256)  :: qty_string, metadata_qty_string

fid = open_file(all_metadata_file)

do iloc = 1,get_model_size()

   call get_model_variable_indices(iloc, ix, iy, iz, &
                                   dom_id=dom_id, &
                                   kind_index=qty_index, &
                                   kind_string=qty_string)

   write(string1,'(i11,1x,''i,j,k'',3(1x,i4),'' domain '',i2)') &
                  iloc, ix, iy, iz, dom_id

   call get_state_meta_data(iloc, loc, var_type)
   metadata_qty_string = trim(get_name_for_quantity(var_type))

   if (trim(qty_string) /= trim(metadata_qty_string) ) then
      write(string2,*)' expected quantity of "'//trim(qty_string)//'"'
      write(string3,*)' got      quantity of "'//trim(metadata_qty_string)//'"'
      call error_handler(E_ERR, 'check_all_meta_data', string1, source, &
                         revision, revdate, text2=string2, text3=string3)
   endif

   call write_location(0,loc,charstring=string2)

   write(string3,'(A,1x,I4,1x,A33,1x,A)') &
         trim(string1), var_type, trim(qty_string), trim(string2)

   if ( do_output()                            ) write(fid,'(A)') trim(string3)
   if ( do_output() .and. mod(iloc,100000) == 0) write( * ,'(A)') trim(string3)

enddo

call close_file(fid)

end subroutine check_all_meta_data


end program model_mod_check

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

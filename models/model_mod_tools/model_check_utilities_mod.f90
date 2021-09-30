! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_check_utilities_mod

!-------------------------------------------------------------------------------
! support routines for interpolation tests
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, metadatalength, missing_r8

use         utilities_mod, only : error_handler, E_MSG, do_output

use     mpi_utilities_mod, only : my_task_id, task_count, &
                                  send_to, receive_from

use          location_mod, only : location_type, &
                                  write_location, &
                                  get_dist, &
                                  query_location, &
                                  has_vertical_choice, &
                                  convert_vertical_state

use          obs_kind_mod, only : get_name_for_quantity, &
                                  get_index_for_quantity

use  ensemble_manager_mod, only : ensemble_type

use             model_mod, only : get_model_size, &
                                  get_state_meta_data, &
                                  model_interpolate

implicit none
private

public :: test_single_interpolation, &
          find_closest_gridpoint, &
          count_error_codes, &
          verify_consistent_istatus

! for messages
character(len=*), parameter :: routine = ''  ! this filename not important in context
character(len=512) :: string1, string2, string3

contains

!-------------------------------------------------------------------------------
!> Do a single interpolation on a given location and kind.
!> Returns the interpolated values and ios_out.
!> Returns the number of ensemble members that passed.

function test_single_interpolation(ens_handle,       &
                                   ens_size,         &
                                   location,         &
                                   quantity_string,  &
                                   interp_vals,      &
                                   ios_out)

type(ensemble_type)   , intent(inout) :: ens_handle
integer               , intent(in)    :: ens_size
type(location_type)   , intent(in)    :: location
character(len=*)      , intent(in)    :: quantity_string
real(r8)              , intent(out)   :: interp_vals(ens_size)
integer               , intent(out)   :: ios_out(ens_size)
integer :: test_single_interpolation

integer :: quantity_index, imem, num_passed
character(len=128) :: my_location

quantity_index = get_index_for_quantity(quantity_string)

if ( do_output() ) then
   call write_location(0, location, charstring=my_location)
   write(*,'(A)') ''
   write(*,'(A)') '-------------------------------------------------------------'
   write(*,'("interpolating at  ",A)') trim(my_location)
   write(*,'("interpolating for ",A)') '"'//trim(quantity_string)//'"'
   write(*,'(A)') '-------------------------------------------------------------'
   write(*,'(A)') ''
endif

call model_interpolate(ens_handle, ens_size, location, quantity_index, interp_vals, ios_out)

num_passed = 0
do imem = 1, ens_size
   if (ios_out(imem) == 0 ) then
      if (do_output()) then
         write(string1,*)'SUCCESS with value    :: ', interp_vals(imem)
         write(*,'(A,I5,A,A)')'member ',imem,',',trim(string1)
         num_passed = num_passed + 1
      endif
   else
      if (do_output()) then
         write(string1,*)'ERROR with error code :: ', ios_out(imem)
         write(*,'(A,I5,A,A)')'member ',imem,',',trim(string1)
      endif
   endif
enddo

if ( do_output() ) write(*,'(A)') ''

test_single_interpolation = num_passed

end function test_single_interpolation


!-------------------------------------------------------------------------------
!> Count the number of different error codes and output the results.
!> This is just a helper function for test_interpolate_range.
!> The input matrix is Nx-by-Nm where:
!> Nx is the number of interpolations and
!> Nm is the ensemble size.

subroutine count_error_codes(error_codes)

integer, intent(in) :: error_codes(:,:)

integer :: i, num_failed

do i =  minval(error_codes(:,:)),  maxval(error_codes(:,:))

   if (i == 0) cycle  ! 0 is a success, not a failure

   num_failed = count(error_codes(:,:) == i)
   if (num_failed > 0 .and. do_output() ) &
       write(*,*) num_failed, " failed with ios_out = ", i

enddo

end subroutine count_error_codes


!-----------------------------------------------------------------------
!> Expensive exhaustive search to find the indices into the
!> state vector of a particular lon/lat/vert. At present, only for a
!> single variable - could be extended to identify the closest location
!> for every variable in each domain. This could help ensure grid
!> staggering is being handled correctly.

subroutine find_closest_gridpoint(location, quantity_string, ens_handle)

type(location_type), intent(in) :: location
character(len=*),    intent(in) :: quantity_string
type(ensemble_type), intent(in) :: ens_handle

! At this point, find_closest_gridpoint() only works for location modules that do not
! have a choice for vertical coordinate systems, so the argument is ignored.
! The vertcoord_string argument is required to make the interface consistent with
! test_interpolate_threed_sphere.f90

type(location_type)   :: loc1
integer               :: i, quantity_index, var_qty, myvars, num_tasks
integer(i8)           :: closest_index, state_index
integer(i8), allocatable :: global_index(:)
real(r8)              :: closest_dist, next_dist
real(r8), allocatable :: global_dist(:)
logical               :: matched
real(r8),   parameter :: FARAWAY = huge(r8)

!>@todo does the following comment means this code doesn't work correctly for multiple 
! domains?  or that it should report on a domain by domain basis?  sorry - i don't know what
! it is trying to say. right now the code iterates the entire state vector and reports 
! the first occurrance of the closest gridpoint, from any domain from any variable of
! the requested quantity.

!>@todo there should be arrays of length state_structure_mod:get_num_variables(domid)
!>      get_num_domains(), get_num_variables() ...


! for single task this is equal to model_size.
! for multiple tasks this is the share of the state
! vector that is on my task.
myvars = ens_handle%my_num_vars

! Trying to support the ability to specify matching a particular QUANTITY.
! With staggered grids, the closest gridpoint might not be of the quantity
! of interest.

quantity_index = get_index_for_quantity(quantity_string)

if (my_task_id() == 0) then
   write(string1,*)'Checking for indices into the state vector that are close to'
   call write_location(0, location, charstring=string2)
   write(string3,*)trim(string2),' for "',trim(quantity_string),'"'
   call error_handler(E_MSG,routine,string1,text2=string3)
   call error_handler(E_MSG,routine,'')
endif

! Change this from previous behavior.  Only save the first
! occurrance of the closest location if more than one has
! the same distance.

closest_dist = FARAWAY
closest_index = -1_r8
matched   = .false.

DISTANCE : do i = 1, myvars

   ! get the global index number for the next item on my task
   ! (same as i if running single task)
   state_index = ens_handle%my_vars(i)

   call get_state_meta_data(state_index, loc1, var_qty)

   if (var_qty .ne. quantity_index) cycle DISTANCE

   next_dist = get_dist(loc1, location)
   if (next_dist < closest_dist) then
      closest_dist = next_dist
      closest_index = state_index
   endif
   matched     = .true.

enddo DISTANCE

if (.not. matched) then
   write(string1,*)'Found no state vector elements of quantity "'//trim(quantity_string)//'"'
   call error_handler(E_MSG, routine, string1)
   return
endif

num_tasks = task_count()

if (num_tasks > 1) then
   
   ! must send distances to task 0 and have it find the closest one
   ! from all tasks.
   if (my_task_id() == 0) then
      allocate(global_dist(0:num_tasks-1))
      allocate(global_index(0:num_tasks-1))

      ! my values, then get the rest from the other tasks
      global_dist(0) = closest_dist
      global_index(0) = state_index
      do i = 1, num_tasks - 1
         call receive_from(i, global_dist(i))
         call receive_from(i, global_index(i))
      enddo

      ! get the index in the distance array that's the smallest
      ! and account for the fact that this array starts at index 0
      ! (mpi tasks are numbered 0 to num_tasks-1, which leads fortran
      ! programming to always be off by one somehow since fortran
      ! arrays start at index 1. sigh.)
      i = minloc(global_dist, 1)
      closest_dist = global_dist(i-1)
      closest_index = global_index(i-1)

      deallocate(global_dist, global_index)
      
   else
      ! send distance and i8 index of my closest item
      call send_to(0, closest_dist)
      call send_to(0, closest_index)
   endif
   
endif

if (my_task_id() == 0) then
   call get_state_meta_data(closest_index, loc1, var_qty)
   call write_location(0, loc1, charstring=string2)
   write(string1,'(A,I12,A)') 'state vector index ', closest_index, ' is closest at '//trim(string2)
   call error_handler(E_MSG, routine, string1)
endif

end subroutine find_closest_gridpoint

!-------------------------------------------------------------------------------

subroutine do_vertical_convert(location, ens_handle)
type(location_type), intent(in) :: location
type(ensemble_type), intent(inout) :: ens_handle

integer :: i, vert_type, myvars, var_qty, istatus
integer(i8) :: state_index
type(location_type) :: loc1
type(location_type), allocatable :: loclist(:)
integer,             allocatable :: qtylist(:)
integer(i8),         allocatable :: indexlist(:)

! if this location type includes more than a single vertical coordinate,
! use the convert routine to convert all the verticals to match the
! 'location of interest'.

if (.not. has_vertical_choice()) return

! type to convert into
vert_type = nint(query_location(location))


! how much of the state vector is on my task?
myvars = ens_handle%my_num_vars

write(string1,*) 'Converting vertical coordinates to match vertical of location of interest'
call write_location(0, location, charstring=string2)
call error_handler(E_MSG,routine,string1,text2=string2,text3=string3)
call error_handler(E_MSG,routine,'')

allocate(loclist(myvars), qtylist(myvars), indexlist(myvars))

! set up the arrays needed for vertical conversion
VERT_CONVERT_SETUP : do i = 1, myvars

   ! get the global index number for the next item on my task
   state_index = ens_handle%my_vars(i)

   call get_state_meta_data(state_index, loc1, var_qty)

   loclist(i) = loc1
   qtylist(i) = var_qty
   indexlist(i) = state_index

enddo VERT_CONVERT_SETUP

! try to convert state locations to match the vertical of the test location
call convert_vertical_state(ens_handle, myvars, loclist, qtylist, indexlist, &
                            vert_type, istatus)


deallocate(loclist, qtylist, indexlist)

if (istatus /= 0) then
   write(string1,*)'unable to convert vertical of state items'
   call error_handler(E_MSG, routine, string1)
endif

end subroutine do_vertical_convert

!-------------------------------------------------------------------------------

subroutine verify_consistent_istatus(ens_size, field, ios_out)
 integer,  intent(in) :: ens_size
 real(r8), intent(in) :: field(ens_size)
 integer,  intent(in) :: ios_out(ens_size)

integer :: i

do i = 1, ens_size
   if (ios_out(i) < 0) then
      write(string1, *) 'ensemble member ', i, &
                        ' inconsistent return: istatus cannot be a negative value.'
      call error_handler(E_MSG, routine, string1)
   endif

   if (ios_out(i) == 0  .and.  field(i) == missing_r8) then
      write(string1, *) 'ensemble member ', i, &
                        ' inconsistent return: istatus = ok but interpolation value = missing data.'
      call error_handler(E_MSG, routine, string1)
   endif

   if (ios_out(i) > 0  .and.  field(i) /= missing_r8) then
      write(string1, *) 'ensemble member ', i, &
                        ' inconsistent return: istatus = error but interpolation value /= missing data.'
      write(string2, *) '  istatus, interp_val = ', ios_out(i), field(i)
      call error_handler(E_MSG, routine, string1, text2=string2)
   endif
enddo

end subroutine verify_consistent_istatus


!-------------------------------------------------------------------------------
! End of model_check_utilities_mod
!-------------------------------------------------------------------------------

end module model_check_utilities_mod


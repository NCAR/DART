! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

module model_check_utilities_mod

!-------------------------------------------------------------------------------
! support routines for interpolation tests
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, metadatalength, missing_r8

use         utilities_mod, only : error_handler, E_MSG, do_output

use          location_mod, only : location_type, &
                                  write_location, &
                                  get_dist, &
                                  set_location

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

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! for messages
character(len=512) :: string1, string2, string3

contains

!-------------------------------------------------------------------------------
!> Do a single interpolation on a given location and kind.
!> Returns the interpolated values and ios_out.
!> Returns the number of ensemble members that passed.

function test_single_interpolation(ens_handle,       &
                                   ens_size,         &
                                   location,         &
                                   vertcoord_string, &
                                   quantity_string,  &
                                   interp_vals,      &
                                   ios_out)

type(ensemble_type)   , intent(inout) :: ens_handle
integer               , intent(in)    :: ens_size
type(location_type)   , intent(in)    :: location
character(len=*)      , intent(in)    :: vertcoord_string
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

subroutine find_closest_gridpoint(loc_of_interest, vertcoord_string, quantity_string)

real(r8),         intent(in) :: loc_of_interest(:)
character(len=*), intent(in) :: vertcoord_string
character(len=*), intent(in) :: quantity_string

! At this point, find_closest_gridpoint() only works for location modules that do not
! have a choice for vertical coordinate systems, so the argument is ignored.
! The vertcoord_string argument is required to make the interface consistent with
! test_interpolate_threed_sphere.f90

character(len=*), parameter :: routine = ''  ! name not important in context

type(location_type)   :: location
type(location_type)   :: loc1
integer(i8)           :: i
integer               :: quantity_index, var_type
real(r8)              :: closest
logical               :: matched
real(r8), allocatable :: thisdist(:)
real(r8),   parameter :: FARAWAY = huge(r8)
character(len=metadatalength) :: myquantity

!>@todo there should be arrays of length state_structure_mod:get_num_variables(domid)
!>      get_num_domains(), get_num_variables() ...

location = set_location(loc_of_interest)

allocate( thisdist(get_model_size()) )
thisdist  = FARAWAY
matched   = .false.

! Trying to support the ability to specify matching a particular QUANTITY.
! With staggered grids, the closest gridpoint might not be of the quantity
! of interest.

quantity_index = get_index_for_quantity(quantity_string)

write(string1,*)'Checking for indices into the state vector that are close to'
call write_location(0, location, charstring=string2)
write(string3,*)trim(string2),' for "',trim(quantity_string),'"'
call error_handler(E_MSG,routine,string1,text2=string3)
call error_handler(E_MSG,routine,'')

! Since there can be/will be multiple variables with
! identical distances, we will just cruise once through
! the array and come back to find all the 'identical' values.

DISTANCE : do i = 1,get_model_size()

   call get_state_meta_data(i, loc1, var_type)

   if (var_type .ne. quantity_index) cycle DISTANCE

   thisdist(i) = get_dist(loc1, location)
   matched     = .true.

enddo DISTANCE

if (.not. matched) then
   write(string1,*)'No state vector elements of type "'//trim(quantity_string)//'"'
   call error_handler(E_MSG, routine, string1)
   deallocate( thisdist )
   return
endif

closest = minval(thisdist)

! Now that we know the distances ... report
! If more than one quantity has the same distance, report all.
! Be aware that if 'approximate_distance' is .true., everything
! in the box has a single location.

REPORT: do i = 1,get_model_size()

   if ( thisdist(i) == closest ) then
      call get_state_meta_data(i, loc1, var_type)
      myquantity = get_name_for_quantity(var_type)

      call write_location(0, loc1, charstring=string2)
      write(string1,'(A,I12,A)')trim(string2)//' is index ',i,' ('//trim(myquantity)//')'
      call error_handler(E_MSG, routine, string1)
   endif

enddo REPORT

deallocate( thisdist )

end subroutine find_closest_gridpoint

!-------------------------------------------------------------------------------

subroutine verify_consistent_istatus(ens_size, field, ios_out)
 integer,  intent(in) :: ens_size
 real(r8), intent(in) :: field(ens_size)
 integer,  intent(in) :: ios_out(ens_size)

character(len=*), parameter :: routine = ''  ! name not important in context
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


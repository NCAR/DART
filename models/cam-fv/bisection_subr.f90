
module bisection

contains

!> given an array of sorted values and a value to find, return the smaller
!> and higher index values, and the fraction across.
!>
!> fraction_across = 0.0 is the 100% the smaller value index, 
!>                   1.0 is the 100% the larger value index.
!>
!> if the array values are inverted (e.g. index 1 is the largest value),
!> set inverted = .true.  the interpretation in the calling code for
!> smaller index, larger index and fraction_across remain the same as the default case.
!>
!> if the fraction_across the enclosing level should be computed using a
!> log scale, set log_scale = .true.
!>
!> my_status values:
!>   0 = good return
!>  -1 = value_to_find is below smallest value
!>   1 = value_to_find is above largest value
!>  96 = cannot use log scale with negative data values
!>  97 = array only has a single value
!>  98 = interval has 0 width or values are inverted
!>  99 = unknown error
!>
!> bad output values use MISSING_I and MISSING_R8
!>
!> this should be in the utilities module.  which should be split into
!> smaller modules because right now it's a dumping ground for every
!> random routine that is useful to more than one module.

integer,  parameter :: r8 = SELECTED_REAL_KIND(12)
integer,  parameter :: MISSING_I = -888888
real(r8), parameter :: MISSING_R8 = -888888.0_r8

subroutine find_enclosing_indices(nitems, data_array, value_to_find,     &
                                  smaller_value_index, larger_value_index, fraction_across, my_status, &
                                  inverted, log_scale)

integer,  intent(in)  :: nitems
real(r8), intent(in)  :: data_array(nitems)
real(r8), intent(in)  :: value_to_find
integer,  intent(out) :: smaller_value_index
integer,  intent(out) :: larger_value_index
real(r8), intent(out) :: fraction_across
integer,  intent(out) :: my_status
logical,  intent(in), optional :: inverted
logical,  intent(in), optional :: log_scale

integer :: i, j, k
logical :: invert, do_log

! set defaults and initialize intent(out) items
! so we can return immediately on error.

invert = .false.
if (present(inverted)) invert = inverted

do_log = .false.
if (present(log_scale)) do_log = log_scale

smaller_value_index = MISSING_I
larger_value_index  = MISSING_I
fraction_across = MISSING_R8
my_status = -99

! exclude malformed call cases
if (nitems <= 1) then
   my_status = 97
   return
endif

! discard out of range values
if ((value_to_find < data_array(1)      .and. .not. invert) .or. &
    (value_to_find < data_array(nitems) .and.       invert)) then
   my_status = -1
   return
endif

if ((value_to_find > data_array(nitems) .and. .not. invert) .or. &
    (value_to_find > data_array(1)      .and.       invert)) then
   my_status = 1
   return
endif

! bisection section (get it?)
i = 1
j = nitems

do
   k=(i+j)/2
   if ((value_to_find < data_array(k) .and. .not. invert) .or. &
       (value_to_find > data_array(k) .and.       invert)) then
      j=k
   else
      i=k
   endif
   if (i+1 >= j) exit
enddo

print *, 'i, j: ', i, i+1
print *, 'data: ', data_array(i), data_array(i+1)

if (.not. invert) then
   smaller_value_index = i
   larger_value_index  = i+1
else
   smaller_value_index = i+1
   larger_value_index  = i
endif

! avoid divide by 0.  the indices have been set but the fraction and status are
! still bad to indicate an error.
if ((data_array(larger_value_index) - data_array(smaller_value_index)) <= 0.0_r8) then
   my_status = 98
   return
endif

! no log computations if any data values are negative
! (do this on 2 lines to avoid testing the data value
!  unless we are planning to take the log.)
if (do_log) then
   if (data_array(smaller_value_index) <= 0.0) then
      my_status = 96
      return
   endif
endif

! compute fraction here
if (.not. do_log) then
   fraction_across = (value_to_find                  - data_array(smaller_value_index)) / &
                     (data_array(larger_value_index) - data_array(smaller_value_index))
else
   fraction_across = (log(value_to_find)                  - log(data_array(smaller_value_index))) / &
                     (log(data_array(larger_value_index)) - log(data_array(smaller_value_index)))

endif

! good return
my_status = 0

end subroutine find_enclosing_indices

end module bisection

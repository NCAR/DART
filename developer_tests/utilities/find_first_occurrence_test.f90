! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod.f90 12563 2018-04-26 21:34:00Z nancy@ucar.edu $

! test of the find_first_occurrence() routine in the utilities module.

program find_first_occurrence_test

use types_mod, only : r8
use utilities_mod, only : error_handler, find_first_occurrence, E_MSG, E_ERR
use sort_mod, only : sort, index_sort
use random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

! make both of these values larger for a real test.
! maybe 50 or 100 for array, 10 or 20 for the samples?
integer, parameter :: MYSIZE = 12
integer, parameter :: NSAMPLES = 2

real(r8), parameter :: BASEVAL = 25.0_r8
real(r8), parameter :: STDDEV = 5.0_r8
type(random_seq_type) :: seq

real(r8) :: array1(MYSIZE), array2(MYSIZE)
real(r8) :: sorted1(MYSIZE), sorted2(MYSIZE)
integer  :: indirect1(MYSIZE)

integer :: seed = 10
integer :: i, this_i, ind_this, istat, j
real(r8) :: thisval, thesevals(MYSIZE), tmp
character(len=512) :: string1

character(len=*), parameter :: routine = 'find_first_occurrence_test'

! calling sequence for routine being tested:
!
!subroutine find_first_occurrence(nitems, data_array, value_to_find, &
!                                 the_index, my_status, inverted, &
!                                 indirect_indicies, the_indirect_index)
!
! should return the largest index number that is less or
! equal to the test value.  last two arguments are optional.

! calling sequence for closely related routine:
!
!subroutine find_enclosing_indices(nitems, data_array, value_to_find,     &
!                                  smaller_index, larger_index, fraction_across, my_status, &
!                                  inverted, log_scale, indirect_indices)
!
! should return the smaller index, larger index, and fraction across.
! last three arguments are optional.



! fill the array with random values and do an initial sanity check
! with none of the special options enabled.

call init_random_seq(seq, seed)

do i=1, MYSIZE
   array1(i) = random_gaussian(seq, BASEVAL, STDDEV)
enddo

! comment this in to add replicated values to the input array
!array1(4) = array1(5)

do i=1, NSAMPLES
   thesevals(i) = random_gaussian(seq, BASEVAL, STDDEV*0.9_r8)
enddo
 
sorted1 = sort(array1)
call index_sort(array1, indirect1, MYSIZE)

! inverted order arrays - cannot do both inverted and indirect!
do i=1, MYSIZE
   array2(i) = array1(MYSIZE - i + 1)
enddo

do i=1, MYSIZE
   sorted2(i) = sorted1(MYSIZE - i + 1)
enddo


! print generated case data

do i=1, MYSIZE
   write(*, '(A50,2(I4,F8.3))') "item, sort, indirect, val: ", i, sorted1(i), indirect1(i), array1(i)
enddo
write(*,*)""

do i=1, MYSIZE
   write(*, '(A50,1(I4,F8.3))') "item, inverted val: ", i, sorted2(i)
enddo
write(*,*)""

do i=1, NSAMPLES
   write(*, '(A50,I4,F8.3)') "test value, data: ", i, thesevals(i)
enddo
write(*,*)""
write(*,*)""


! real start of tests

write(*,*)"direct tests"

! test edge cases
thisval = sorted1(1)
write(string1, *) "lowest value"
call find_first_occurrence(MYSIZE, sorted1, thisval, this_i, istat)
call print_results(thisval, sorted1, this_i, istat, string1)

thisval = sorted1(MYSIZE)
write(string1, *) "highest value"
call find_first_occurrence(MYSIZE, sorted1, thisval, this_i, istat)
call print_results(thisval, sorted1, this_i, istat, string1)

! test outside range
thisval = sorted1(1) - 1.0_r8
write(string1, *) "below lowest value"
call find_first_occurrence(MYSIZE, sorted1, thisval, this_i, istat)
call print_results(thisval, sorted1, this_i,  istat, string1)

thisval = sorted1(MYSIZE) + 1.0_r8
write(string1, *) "above highest value"
call find_first_occurrence(MYSIZE, sorted1, thisval, this_i, istat)
call print_results(thisval, sorted1, this_i, istat, string1)

! tests for equals
do i=1, MYSIZE
   thisval = sorted1(i)
   write(string1, '(A22,I4,F8.3)') "equal, val: ", i, thisval
   call find_first_occurrence(MYSIZE, sorted1, thisval, this_i, istat)
   call print_results(thisval, sorted1, this_i, istat, string1)
enddo

! tests for non-equals
do i=1, NSAMPLES
   write(string1, '(A22,I4,F8.3)') "noneq, val: ", i, thesevals(i)
   call find_first_occurrence(MYSIZE, sorted1, thesevals(i), this_i, istat)
   call print_results(thesevals(i), sorted1, this_i, istat, string1)
enddo

write(*,*)""
write(*,*)"indirect tests"

! test indirect addressing

! test edge cases
thisval = sorted1(1)
write(string1, *) "lowest value"
call find_first_occurrence(MYSIZE, array1, thisval, this_i, istat, &
                           indirect_indices = indirect1, the_indirect_index = ind_this)
call print_results(thisval, array1, this_i, istat, string1, indirect1, indirect_this = ind_this)

thisval = sorted1(MYSIZE)
write(string1, *) "highest value"
call find_first_occurrence(MYSIZE, array1, thisval, this_i, istat, &
                           indirect_indices = indirect1, the_indirect_index = ind_this)
call print_results(thisval, array1, this_i, istat, string1, indirect1, indirect_this = ind_this)

! test outside range
thisval = sorted1(1) - 1.0_r8
write(string1, *) "below lowest value"
call find_first_occurrence(MYSIZE, array1, thisval, this_i, istat, &
                           indirect_indices = indirect1, the_indirect_index = ind_this)
call print_results(thisval, array1, this_i, istat, string1, indirect1, indirect_this = ind_this)

thisval = sorted1(MYSIZE) + 1.0_r8
write(string1, *) "above highest value"
call find_first_occurrence(MYSIZE, array1, thisval, this_i, istat, &
                           indirect_indices = indirect1, the_indirect_index = ind_this)
call print_results(thisval, array1, this_i, istat, string1, indirect1, indirect_this = ind_this)

! tests for equals
do i=1, MYSIZE
   thisval = sorted1(i)
   write(string1, '(A22,I4,F8.3)') "equal, val: ", i, thisval
   call find_first_occurrence(MYSIZE, array1, thisval, this_i, istat, &
                              indirect_indices = indirect1, the_indirect_index = ind_this)
   call print_results(thisval, array1, this_i, istat, string1, indirect1, indirect_this = ind_this)
enddo

! tests for non-equals
do i=1, NSAMPLES
   write(string1, '(A22,I4,F8.3)') "noneq, val: ", i, thesevals(i)
   call find_first_occurrence(MYSIZE, array1, thesevals(i), this_i, istat, &
                              indirect_indices = indirect1, the_indirect_index = ind_this)
   call print_results(thesevals(i), array1, this_i, istat, string1, indirect1, indirect_this = ind_this)
enddo

write(*,*)""
write(*,*)"inverted tests"

! test inverted arrays

! test edge cases
thisval = sorted1(1)
write(string1, *) "lowest value"
call find_first_occurrence(MYSIZE, sorted2, thisval, this_i, istat, &
                            inverted = .true.)
call print_results(thisval, sorted2, this_i, istat, string1, inverted=.true.)

thisval = sorted1(MYSIZE)
write(string1, *) "highest value"
call find_first_occurrence(MYSIZE, sorted2, thisval, this_i, istat, &
                            inverted = .true.)
call print_results(thisval, sorted2, this_i, istat, string1, inverted=.true.)

! test outside range
thisval = sorted1(1) - 1.0_r8
write(string1, *) "below lowest value"
call find_first_occurrence(MYSIZE, sorted2, thisval, this_i, istat, &
                            inverted = .true.)
call print_results(thisval, sorted2, this_i, istat, string1, inverted=.true.)

thisval = sorted1(MYSIZE) + 1.0_r8
write(string1, *) "above highest value"
call find_first_occurrence(MYSIZE, sorted2, thisval, this_i, istat, &
                            inverted = .true.)
call print_results(thisval, sorted2, this_i, istat, string1, inverted=.true.)

! tests for equals
do i=1, MYSIZE
   thisval = sorted1(i)
   write(string1, '(A22,I4,F8.3)') "equal, val: ", i, thisval
   call find_first_occurrence(MYSIZE, sorted2, thisval, this_i, istat, &
                              inverted = .true.)
   call print_results(thisval, sorted2, this_i, istat, string1, inverted=.true.)
enddo

! tests for non-equals
do i=1, NSAMPLES
   write(string1, '(A22,I4,F8.3)') "noneq, val: ", i, thesevals(i)
   call find_first_occurrence(MYSIZE, sorted2, thesevals(i), this_i, istat, inverted=.true.)
   call print_results(thesevals(i), sorted2, this_i, istat, string1, inverted=.true.)
enddo

write(*,*)""
write(*,*)"invalid input tests"

! empty array
write(string1, '(A6,I4)') "single item array"
call find_first_occurrence(0, sorted1(1:1), thesevals(1), this_i, istat)
call print_results(thesevals(1), sorted1, this_i, istat, string1)

! indirect and inverted
write(string1, '(A6,I4)') "indirect and inverted"
call find_first_occurrence(MYSIZE, sorted1, thesevals(1), this_i, istat, &
                           inverted=.true., indirect_indices=indirect1)
call print_results(thesevals(1), sorted1, this_i, istat, string1)

! indirect_this without indirect array
write(string1, '(A6,I4)') "indirect this w/o array"
call find_first_occurrence(MYSIZE, array1, thisval, this_i, istat, &
                           the_indirect_index = ind_this)
call print_results(thisval, array1, this_i, istat, string1, indirect1, indirect_this = ind_this)
call print_results(thesevals(1), array1, this_i, istat, string1)

! inverted (badly sorted) input array
write(string1, '(A6,I4)') "non-monotonic array"
call find_first_occurrence(MYSIZE, array1, array1(MYSIZE/2), this_i, istat)
call print_results(array1(MYSIZE/2), array1, this_i, istat, string1)

! almost sorted input array
write(string1, '(A6,I4)') "mostly sorted array"
j = MYSIZE/3
tmp = sorted1(j)
sorted1(j) = sorted1(j*2)
sorted1(j*2) = tmp
call find_first_occurrence(MYSIZE, sorted1, sorted1(MYSIZE/2), this_i, istat)
call print_results(sorted1(MYSIZE/2), sorted1, this_i, istat, string1)


write(*,*) 'end of test'

contains

subroutine print_results(thisval, thisarray, this_i, istat, label, indirect_a, inverted, indirect_this)
 real(r8), intent(in) :: thisval
 real(r8), intent(in) :: thisarray(:)
 integer,  intent(in) :: this_i
 integer,  intent(in) :: istat
 character(len=*), intent(in) :: label
 integer,  intent(in), optional :: indirect_a(:)
 logical,  intent(in), optional :: inverted
 integer,  intent(in), optional :: indirect_this

integer :: this, next

if (istat /= 0) then
   write(*, '(A56,2I8)') trim(label) // ' index, status = ', this_i, istat
   return
endif

write(*,'(A56,F8.3,I4,2F8.3)') trim(label) // ' val, indx, arrval = ', &
                               thisval, this_i, thisarray(this_i)

! indirect_a has to be present if indirect_this is specified.
if (present(indirect_this)) then
    if (.not. present(indirect_a)) then
       print *, 'bad call to print_results: indirect_this specified but not indirect_a'
       stop
    endif
    this = indirect_a(indirect_this)
    next = indirect_a(min(indirect_this + 1, size(thisarray)))
else if (present(inverted)) then
    this = this_i
    next = max(this_i - 1, 1)
else
    this = this_i
    next = min(this_i + 1, size(thisarray))
endif

if (thisval < thisarray(this) .or. thisval > thisarray(next)) then
   write(*,'(A,3F8.3)') 'unexpected error - val not between the two values, ', &
                         thisval, thisarray(this), thisarray(next)
endif

end subroutine print_results

end program find_first_occurrence_test


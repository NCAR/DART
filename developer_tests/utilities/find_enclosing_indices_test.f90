! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod.f90 12563 2018-04-26 21:34:00Z nancy@ucar.edu $

! test of the find_enclosing_indices() routine in the utilities module.

program find_enclosing_indicies_test

use types_mod, only : r8
use utilities_mod, only : error_handler, find_enclosing_indices, E_MSG, E_ERR
use sort_mod, only : sort, index_sort
use random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

integer, parameter :: MYSIZE = 25
integer, parameter :: NSAMPLES = 4
type(random_seq_type) :: seq

real(r8), parameter :: BASEVAL = 25.0_r8
real(r8), parameter :: STDDEV = 5.0_r8

real(r8) :: array1(MYSIZE), array2(MYSIZE)
real(r8) :: sorted1(MYSIZE), sorted2(MYSIZE)
integer  :: indirect1(MYSIZE), indirect2(MYSIZE)

integer :: i, lower_i, upper_i, istat, j
real(r8) :: fract, thisval, thesevals(MYSIZE), tmp
character(len=512) :: string1

character(len=*), parameter :: routine = 'find_enclosing_indicies_test'

! fill the array with random values and do an initial sanity check
! with none of the special options enabled.

call init_random_seq(seq)

do i=1, MYSIZE
   array1(i) = random_gaussian(seq, BASEVAL, STDDEV)
enddo

do i=1, NSAMPLES
   thesevals(i) = random_gaussian(seq, BASEVAL, STDDEV * 0.9_r8)
enddo
 
sorted1 = sort(array1)
call index_sort(array1, indirect1, MYSIZE)

! inverted order arrays
do i=1, MYSIZE
   array2(i) = array1(MYSIZE - i + 1)
enddo

do i=1, MYSIZE
   sorted2(i) = sorted1(MYSIZE - i + 1)
enddo

do i=1, MYSIZE
   indirect2(i) = indirect1(MYSIZE - i + 1)
enddo



! print generated case data

do i=1, MYSIZE
   write(*, '(A32,2(I4,F8.3))') "original, sorted data, indirect, ", i, sorted1(i), indirect1(i), array1(i)
enddo
write(*,*)""

do i=1, MYSIZE
   write(*, '(A32,I4,2F8.3)') "inverted, sorted data, ", i, array2(i), sorted2(i)
enddo
write(*,*)""

do i=1, NSAMPLES
   write(*, '(A32,I4,F8.3)') "sampled data, item ", i, thesevals(i)
enddo
write(*,*)""
write(*,*)""


! real start of tests

write(*,*)"direct tests"

! test edge cases
thisval = sorted1(1)
write(string1, *) "lowest value"
call find_enclosing_indices(MYSIZE, sorted1, thisval, lower_i, upper_i, fract, istat)
call print_results(thisval, sorted1, lower_i, upper_i, fract, istat, string1)

thisval = sorted1(MYSIZE)
write(string1, *) "highest value"
call find_enclosing_indices(MYSIZE, sorted1, thisval, lower_i, upper_i, fract, istat)
call print_results(thisval, sorted1, lower_i, upper_i, fract, istat, string1)

! test outsize range
thisval = sorted1(1) - 1.0_r8
write(string1, *) "below lowest value"
call find_enclosing_indices(MYSIZE, sorted1, thisval, lower_i, upper_i, fract, istat)
call print_results(thisval, sorted1, lower_i, upper_i, fract, istat, string1)

thisval = sorted1(MYSIZE) + 1.0_r8
write(string1, *) "above highest value"
call find_enclosing_indices(MYSIZE, sorted1, thisval, lower_i, upper_i, fract, istat)
call print_results(thisval, sorted1, lower_i, upper_i, fract, istat, string1)

! basic case
do i=1, NSAMPLES
   write(string1, '(A6,I4)') "item", i
   call find_enclosing_indices(MYSIZE, sorted1, thesevals(i), lower_i, upper_i, fract, istat)
   call print_results(thesevals(i), sorted1, lower_i, upper_i, fract, istat, string1)
enddo

write(*,*)""
write(*,*)"indirect tests"

! test indirect addressing

! test edge cases
thisval = sorted1(1)
write(string1, *) "lowest value"
call find_enclosing_indices(MYSIZE, array1, thisval, lower_i, upper_i, fract, istat, &
                               indirect_indices = indirect1)
call print_results(thisval, array1, lower_i, upper_i, fract, istat, string1)

thisval = sorted1(MYSIZE)
write(string1, *) "highest value"
call find_enclosing_indices(MYSIZE, array1, thisval, lower_i, upper_i, fract, istat, &
                            indirect_indices = indirect1)
call print_results(thisval, array1, lower_i, upper_i, fract, istat, string1)

! test outsize range
thisval = sorted1(1) - 1.0_r8
write(string1, *) "below lowest value"
call find_enclosing_indices(MYSIZE, array1, thisval, lower_i, upper_i, fract, istat, &
                            indirect_indices = indirect1)
call print_results(thisval, array1, lower_i, upper_i, fract, istat, string1)

thisval = sorted1(MYSIZE) + 1.0_r8
write(string1, *) "above highest value"
call find_enclosing_indices(MYSIZE, array1, thisval, lower_i, upper_i, fract, istat, &
                            indirect_indices = indirect1)
call print_results(thisval, array1, lower_i, upper_i, fract, istat, string1)

! basic case
do i=1, NSAMPLES
   write(string1, '(A6,I4)') "item", i
   call find_enclosing_indices(MYSIZE, array1, thesevals(i), lower_i, upper_i, fract, istat, &
                               indirect_indices = indirect1)
   call print_results(thesevals(i), array1, lower_i, upper_i, fract, istat, string1)
enddo

write(*,*)""
write(*,*)"inverted tests"

! test inverted arrays

! test edge cases
thisval = sorted1(1)
write(string1, *) "lowest value"
call find_enclosing_indices(MYSIZE, sorted2, thisval, lower_i, upper_i, fract, istat, &
                            inverted = .true.)
call print_results(thisval, sorted2, lower_i, upper_i, fract, istat, string1)

thisval = sorted1(MYSIZE)
write(string1, *) "highest value"
call find_enclosing_indices(MYSIZE, sorted2, thisval, lower_i, upper_i, fract, istat, &
                            inverted = .true.)
call print_results(thisval, sorted2, lower_i, upper_i, fract, istat, string1)

! test outsize range
thisval = sorted1(1) - 1.0_r8
write(string1, *) "below lowest value"
call find_enclosing_indices(MYSIZE, sorted2, thisval, lower_i, upper_i, fract, istat, &
                            inverted = .true.)
call print_results(thisval, sorted2, lower_i, upper_i, fract, istat, string1)

thisval = sorted1(MYSIZE) + 1.0_r8
write(string1, *) "above highest value"
call find_enclosing_indices(MYSIZE, sorted2, thisval, lower_i, upper_i, fract, istat, &
                            inverted = .true.)
call print_results(thisval, sorted2, lower_i, upper_i, fract, istat, string1)

! basic case
do i=1, NSAMPLES
   write(string1, '(A6,I4)') "item", i
   call find_enclosing_indices(MYSIZE, sorted2, thesevals(i), lower_i, upper_i, fract, istat, inverted=.true.)
   call print_results(thesevals(i), sorted2, lower_i, upper_i, fract, istat, string1)
enddo

write(*,*)""
write(*,*)"inverted indirect tests"

! inverted, indirect arrays - should all fail.

! test edge cases
thisval = sorted1(1)
write(string1, *) "inv lowest value"
call find_enclosing_indices(MYSIZE, array2, thisval, lower_i, upper_i, fract, istat, inverted=.true., &
                            indirect_indices = indirect1)
call print_results(thisval, array2, lower_i, upper_i, fract, istat, string1)

thisval = sorted1(MYSIZE)
write(string1, *) "inv highest value"
call find_enclosing_indices(MYSIZE, array2, thisval, lower_i, upper_i, fract, istat, inverted=.true., &
                            indirect_indices = indirect1)
call print_results(thisval, array2, lower_i, upper_i, fract, istat, string1)

! test outsize range
thisval = sorted1(1) - 1.0_r8
write(string1, *) "inv below lowest value"
call find_enclosing_indices(MYSIZE, array2, thisval, lower_i, upper_i, fract, istat, inverted=.true., &
                            indirect_indices = indirect1)
call print_results(thisval, array2, lower_i, upper_i, fract, istat, string1)

thisval = sorted1(MYSIZE) + 1.0_r8
write(string1, *) "inv above highest value"
call find_enclosing_indices(MYSIZE, array2, thisval, lower_i, upper_i, fract, istat, inverted=.true., &
                            indirect_indices = indirect1)
call print_results(thisval, array2, lower_i, upper_i, fract, istat, string1)

! 1 test indirect addressing

write(string1, '(A6,I4)') "item", 1
call find_enclosing_indices(MYSIZE, array2, thesevals(1), lower_i, upper_i, fract, istat, &
                            indirect_indices = indirect2, inverted=.true.)
call print_results(thesevals(1), array2, lower_i, upper_i, fract, istat, string1)

write(*,*)""
write(*,*)"invalid input tests"

! single array item
write(string1, '(A6,I4)') "single item array"
call find_enclosing_indices(1, sorted1(1:1), thesevals(1), lower_i, upper_i, fract, istat)
call print_results(thesevals(1), sorted1, lower_i, upper_i, fract, istat, string1)

! inverted and indirect
write(string1, '(A6,I4)') "inverted and indirect"
call find_enclosing_indices(1, sorted1(1:1), thesevals(1), lower_i, upper_i, fract, istat, &
                   inverted=.true., indirect_indices = indirect1)
call print_results(thesevals(1), sorted1, lower_i, upper_i, fract, istat, string1)

! inverted (badly sorted) input array
write(string1, '(A6,I4)') "non-monotonic array"
call find_enclosing_indices(MYSIZE, array1, array1(MYSIZE/2), lower_i, upper_i, fract, istat)
call print_results(array1(MYSIZE/2), array1, lower_i, upper_i, fract, istat, string1)

! almost sorted input array
write(string1, '(A6,I4)') "mostly sorted array"
j = MYSIZE/3
tmp = sorted1(j)
sorted1(j) = sorted1(j*2)
sorted1(j*2) = tmp
call find_enclosing_indices(MYSIZE, sorted1, array1(MYSIZE/2), lower_i, upper_i, fract, istat)
call print_results(array1(MYSIZE/2), sorted1, lower_i, upper_i, fract, istat, string1)

write(*,*) 'end of test'

contains

subroutine print_results(thisval, thisarray, lower_i, upper_i, fract, istat, label)
 real(r8), intent(in) :: thisval
 real(r8), intent(in) :: thisarray(:)
 integer,  intent(in) :: lower_i, upper_i
 real(r8), intent(in) :: fract
 integer,  intent(in) :: istat
 character(len=*), intent(in) :: label

real(r8) :: computed_val1, computed_val2

write(*,'(A32,F8.3,A10,I4)') trim(label), thisval, ' status = ', istat

if (istat /= 0) then
   write(*,*) ''
   return
endif

computed_val1 = thisarray(lower_i) + (thisarray(upper_i) - thisarray(lower_i))*fract
computed_val2 = thisarray(lower_i) + (thisarray(upper_i) - thisarray(lower_i))*(1.0 - fract)
write(*, '(A32,2I4,F8.3)') "  lower, upper, fract = ", lower_i, upper_i, fract
write(*, '(A32,4F8.3)') "  low, cval1, up = ", thisarray(lower_i), computed_val1, thisarray(upper_i)
!write(*, '(A32,4F8.3)') "  low, cval1, cval2, up = ", thisarray(lower_i), computed_val1, computed_val2, thisarray(upper_i)

if (computed_val1 /= thisval) write(*, *) "warning! mismatched values: ", thisval, computed_val1

write(*,*) ''

end subroutine print_results

end program find_enclosing_indicies_test


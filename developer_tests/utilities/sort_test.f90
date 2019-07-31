! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod.f90 12563 2018-04-26 21:34:00Z nancy@ucar.edu $

! test of the various sort routines in the utilities module.
!
! things to add: 
!  repeat the sorts (w/ different but reproducible data?) and time them
!  do the sorts on lists of different lengths.
!  output stats

program sort_test

use types_mod, only : r8, digits12
use utilities_mod, only : error_handler, E_MSG, E_ERR
!use sort_mod, only : sort, index_sort, simple_sort
use sort_mod, only : simple_sort, simple_index_sort, sort, index_sort, insertion_sort, index_insertion_sort
use random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian
use mpi_utilities_mod, only : start_mpi_timer, read_mpi_timer

!> usage:
!>  real(digits12) :: base, time_elapsed
!>
!>  call start_mpi_timer(base)
!>  time_elapsed = read_mpi_timer(base)


integer, parameter :: MAXSIZE = 25000
integer, parameter :: SAMPLE = 10
type(random_seq_type) :: seq

real(r8), parameter :: BASEVAL = 25.0_r8
real(r8), parameter :: STDDEV = 5.0_r8

real(r8) :: array1(MAXSIZE), array2(MAXSIZE)
real(r8) :: sorted1(MAXSIZE), sorted2(MAXSIZE)
integer  :: indirect1(MAXSIZE), indirect2(MAXSIZE)

integer  :: i
real(digits12) :: b(10), inc
real(r8) :: fract, thisval, thesevals(MAXSIZE), tmp
character(len=512) :: string1

character(len=*), parameter :: routine = 'sort_test'

! fill the array with random values for sorting.

call init_random_seq(seq)

do i=1, MAXSIZE
   array1(i) = random_gaussian(seq, BASEVAL, STDDEV)
enddo


write(*,*) ""
write(*,*) "random"

call start_mpi_timer(b(1))
sorted1 = sort(array1)
inc = read_mpi_timer(b(1))
call print_time(inc, 'basic sort')

call start_mpi_timer(b(2))
call index_sort(array1, indirect1, MAXSIZE)
inc = read_mpi_timer(b(2))
call print_time(inc, 'index sort')


sorted1 = array1
call start_mpi_timer(b(1))
call insertion_sort(array1)
inc = read_mpi_timer(b(1))
call print_time(inc, 'insertion sort')

call start_mpi_timer(b(2))
call index_insertion_sort(array1, indirect1, MAXSIZE)
inc = read_mpi_timer(b(2))
call print_time(inc, 'index insertion sort')

write(*,*) ""
write(*,*) "inverted order"

! inverted order arrays
do i=1, MAXSIZE
   array2(i) = array1(MAXSIZE - i + 1)
enddo

do i=1, MAXSIZE
   sorted2(i) = sorted1(MAXSIZE - i + 1)
enddo


call start_mpi_timer(b(1))
sorted2 = sort(array2)
inc = read_mpi_timer(b(1))
call print_time(inc, 'basic sort')

call start_mpi_timer(b(2))
call index_sort(sorted2, indirect2, MAXSIZE)
inc = read_mpi_timer(b(2))
call print_time(inc, 'index sort')


call start_mpi_timer(b(1))
call insertion_sort(array2)
inc = read_mpi_timer(b(1))
call print_time(inc, 'insertion sort')

call start_mpi_timer(b(2))
call index_insertion_sort(sorted2, indirect1, MAXSIZE)
inc = read_mpi_timer(b(2))
call print_time(inc, 'index insertion sort')

write(*,*) ""


!! print generated case data
!
!do i=1, SAMPLE
!   write(*, '(A32,2(I4,F8.3))') "original, sorted data, indirect, ", i, array1(i), indirect1(i), sorted1(i)
!enddo
!write(*,*)""
!
!do i=1, SAMPLE
!   write(*, '(A32,I4,2F8.3)') "inverted, sorted data, ", i, array2(i), sorted2(i)
!enddo
!write(*,*)""


! ! inverted input array
! write(string1, '(A6,I4)') "non-monotonic array"

! ! almost sorted input array
! write(string1, '(A6,I4)') "mostly sorted array"
j = MAXSIZE/3
tmp = sorted1(j)
sorted1(j) = sorted1(j*2)
sorted1(j*2) = tmp

write(*,*) 'end of test'

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

subroutine validate(thisarray)
 real(r8), intent(in) :: thisarray(:)

integer :: i, isize

isize = size(thisarray)

do i=2, isize
   if (thisarray(i-1) > thisarray(i)) then
      write(*, *) 'items not in order: ', i-1, thisarray(i-1), &
                  ' not less than or equal to ', i, thisarray(i) 
   endif
enddo

end subroutine validate

!-------------------------------------------------------------------------------

subroutine validate_indx(thisarray, thisindx)
 real(r8), intent(in) :: thisarray(:)
 integer,  intent(in) :: thisindx(:)

integer :: i, isize

isize = size(thisarray)

do i=2, isize
   if (thisarray(thisindx(i-1)) > thisarray(thisindx(i))) then
      write(*, *) 'items not in order: ', i-1, thisindx(i-1), thisarray(thisindx(i-1)), &
                  ' not less than or equal to ', i, thisindx(i), thisarray(thisindx(i))
   endif
enddo

end subroutine validate_indx

!-------------------------------------------------------------------------------

subroutine print_time(timespan, label)
 real(digits12),   intent(in) :: timespan
 character(len=*), intent(in) :: label

write(*, '(A32,F12.2)') trim(label)//' ', timespan

end subroutine print_time

!-------------------------------------------------------------------------------

end program sort_test


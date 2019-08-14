! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! test of the various sort routines in the utilities module.
!
! things to add: 
!  repeat the sorts (w/ different but reproducible data?) and time them
!  to get better statistics

program sort_test

use types_mod,         only : r8, digits12

use utilities_mod,     only : error_handler, E_MSG, E_ERR, initialize_utilities, finalize_utilities

use sort_mod,          only : simple_sort, simple_index_sort, sort, index_sort, insertion_sort, index_insertion_sort

use random_seq_mod,    only: random_seq_type, init_random_seq, random_gaussian

use mpi_utilities_mod, only : start_mpi_timer, read_mpi_timer

!> timer usage example:
!>
!>  real(digits12) :: base, time_elapsed
!>
!>  call start_mpi_timer(base)
!>  ! do stuff here
!>  time_elapsed = read_mpi_timer(base)


! buffer sizes
integer, parameter :: MAXSIZE = 20000
integer, parameter :: SMALLER =  2000
integer, parameter :: TINY    =    80

type(random_seq_type) :: seq

! fill arrays with gaussian distribution of numeric values 
real(r8), parameter :: BASEVAL = 250.0_r8
real(r8), parameter :: STDDEV  = 50.0_r8

real(r8) :: base_real_array(MAXSIZE)
integer  :: base_int_array(MAXSIZE)

! real and integer arrays of 'm' (maxsize), 's' (smaller)
! and 't' (tiny) for sorting and timing.
real(r8) :: mr_array1(MAXSIZE), mr_array2(MAXSIZE)
integer  :: mi_array1(MAXSIZE), mi_array2(MAXSIZE)
real(r8) :: sr_array1(SMALLER), sr_array2(SMALLER)
integer  :: si_array1(SMALLER), si_array2(SMALLER)
real(r8) :: tr_array1(TINY),    tr_array2(TINY)
integer  :: ti_array1(TINY),    ti_array2(TINY)

integer  :: i, itmp
real(r8) :: tmp

character(len=*), parameter :: routine = 'sort_test'

! start of executable code

call initialize_utililities('sort_test')

print *, 'sort test.  array sizes: '
print *, '  large = ', MAXSIZE
print *, 'smaller = ', SMALLER
print *, '   tiny = ', TINY
print *, ''
print *, 'time in milliseconds.'
print *, ''

! fill the base array with random values for sorting.

call init_random_seq(seq)

do i=1, MAXSIZE
   base_real_array(i) = random_gaussian(seq, BASEVAL, STDDEV)
enddo

base_int_array(:) = floor(base_real_array(:))

mr_array1(:) = base_real_array(:)
mi_array1(:) = floor(mr_array1(:))

sr_array1(:) = base_real_array(1:SMALLER)
si_array1(:) = floor(sr_array1(:))

tr_array1(:) = base_real_array(1:TINY)
ti_array1(:) = floor(tr_array1(:))


write(*,*) ""
write(*,*) "random, large"

call do_each_real_type(mr_array1, MAXSIZE)
call do_each_int_type (mi_array1, MAXSIZE)

write(*,*) ""
write(*,*) "random, small"

call do_each_real_type(sr_array1, SMALLER)
call do_each_int_type (si_array1, SMALLER)

write(*,*) ""
write(*,*) "random, tiny"

call do_each_real_type(tr_array1, TINY)
call do_each_int_type (ti_array1, TINY)


! inverted order arrays
do i=1, MAXSIZE
   mr_array2(i) = base_real_array(MAXSIZE - i + 1)
enddo
mi_array2(:) = floor(mr_array2(:))

sr_array2(:) = mr_array2(1:SMALLER)
si_array2(:) = floor(sr_array2(:))

tr_array2(:) = mr_array2(1:TINY)
ti_array2(:) = floor(tr_array2(:))


write(*,*) ""
write(*,*) "inverted order, large"

call do_each_real_type(mr_array2, MAXSIZE)
call do_each_int_type (mi_array2, MAXSIZE)

write(*,*) ""
write(*,*) "inverted order, small"

call do_each_real_type(sr_array2, SMALLER)
call do_each_int_type (si_array2, SMALLER)

write(*,*) ""
write(*,*) "inverted order, tiny"

call do_each_real_type(tr_array2, TINY)
call do_each_int_type (ti_array2, TINY)


! almost sorted order arrays
mr_array2(:) = sort(base_real_array)
mi_array2(:) = floor(mr_array2(:))

call reorder_pair(mr_array2, mi_array2)
call reorder_pair(sr_array2, si_array2)
call reorder_pair(tr_array2, ti_array2)


write(*,*) ""
write(*,*) "almost sorted order, large"

call do_each_real_type(mr_array2, MAXSIZE)
call do_each_int_type (mi_array2, MAXSIZE)

write(*,*) ""
write(*,*) "almost sorted order, small"

call do_each_real_type(sr_array2, SMALLER)
call do_each_int_type (si_array2, SMALLER)

write(*,*) ""
write(*,*) "almost sorted order, tiny"

call do_each_real_type(tr_array2, TINY)
call do_each_int_type (ti_array2, TINY)


write(*,*) 'end of test'

call finalize_utililities('sort_test')

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! here is where templated functions would be nice.

subroutine do_each_real_type(array, asize)
 real(r8), intent(in) :: array(:)
 integer,  intent(in) :: asize
 
real(digits12) :: b, inc
real(r8) :: array1(asize), sorted1(asize)
integer :: indirect1(asize)


call start_mpi_timer(b)
sorted1 = sort(array)
inc = read_mpi_timer(b)
call print_time(inc, 'basic sort (real)')
call validate(sorted1, asize)

call start_mpi_timer(b)
call index_sort(array, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, 'index sort (real)')
call validate_indx(array, indirect1, asize)

array1 = array

call start_mpi_timer(b)
call insertion_sort(array1)
inc = read_mpi_timer(b)
call print_time(inc, 'insertion sort (real)')
call validate(array1, asize)

call start_mpi_timer(b)
call index_insertion_sort(array, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, 'index insertion sort (real)')
call validate_indx(array, indirect1, asize)

array1 = array

call start_mpi_timer(b)
call simple_sort(array1)
inc = read_mpi_timer(b)
call print_time(inc, 'simple sort (real)')
call validate(array1, asize)

call start_mpi_timer(b)
call simple_index_sort(array, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, 'index simple sort (real)')
call validate_indx(array, indirect1, asize)

end subroutine do_each_real_type

!-------------------------------------------------------------------------------

subroutine do_each_int_type(array, asize)
 integer, intent(in) :: array(:)
 integer, intent(in) :: asize
 
real(digits12) :: b, inc
integer :: array1(asize), sorted1(asize)
integer :: indirect1(asize)


call start_mpi_timer(b)
sorted1 = sort(array)
inc = read_mpi_timer(b)
call print_time(inc, 'basic sort (integer)')
call ivalidate(sorted1, asize)

call start_mpi_timer(b)
call index_sort(array, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, 'index sort (integer)')
call ivalidate_indx(array, indirect1, asize)

array1 = array

call start_mpi_timer(b)
call insertion_sort(array1)
inc = read_mpi_timer(b)
call print_time(inc, 'insertion sort (integer)')
call ivalidate(array1, asize)

call start_mpi_timer(b)
call index_insertion_sort(array, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, 'index insertion sort (integer)')
call ivalidate_indx(array, indirect1, asize)

array1 = array

call start_mpi_timer(b)
call simple_sort(array1)
inc = read_mpi_timer(b)
call print_time(inc, 'simple sort (integer)')
call ivalidate(array1, asize)

call start_mpi_timer(b)
call simple_index_sort(array, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, 'index simple sort (integer)')
call ivalidate_indx(array, indirect1, asize)

end subroutine do_each_int_type


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine validate(thisarray, asize)
 real(r8), intent(in) :: thisarray(:)
 integer,  intent(in) :: asize

integer :: i

do i=2, asize
   if (thisarray(i-1) > thisarray(i)) then
      write(*, *) 'items not in order: ', i-1, thisarray(i-1), &
                  ' not less than or equal to ', i, thisarray(i) 
   endif
enddo

end subroutine validate

!-------------------------------------------------------------------------------

subroutine validate_indx(thisarray, thisindx, asize)
 real(r8), intent(in) :: thisarray(:)
 integer,  intent(in) :: thisindx(:)
 integer,  intent(in) :: asize

integer :: i

do i=2, asize
   if (thisarray(thisindx(i-1)) > thisarray(thisindx(i))) then
      write(*, *) 'items not in order: ', i-1, thisindx(i-1), thisarray(thisindx(i-1)), &
                  ' not less than or equal to ', i, thisindx(i), thisarray(thisindx(i))
   endif
enddo

end subroutine validate_indx

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine ivalidate(thisarray, asize)
 integer, intent(in) :: thisarray(:)
 integer, intent(in) :: asize

integer :: i

do i=2, asize
   if (thisarray(i-1) > thisarray(i)) then
      write(*, *) 'items not in order: ', i-1, thisarray(i-1), &
                  ' not less than or equal to ', i, thisarray(i) 
   endif
enddo

end subroutine ivalidate

!-------------------------------------------------------------------------------

subroutine ivalidate_indx(thisarray, thisindx, asize)
 integer, intent(in) :: thisarray(:)
 integer, intent(in) :: thisindx(:)
 integer, intent(in) :: asize

integer :: i

do i=2, asize
   if (thisarray(thisindx(i-1)) > thisarray(thisindx(i))) then
      write(*, *) 'items not in order: ', i-1, thisindx(i-1), thisarray(thisindx(i-1)), &
                  ' not less than or equal to ', i, thisindx(i), thisarray(thisindx(i))
   endif
enddo

end subroutine ivalidate_indx

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! swap the values at about index N*1/3 with the value at about N*2/3

subroutine reorder_pair(rarray, iarray)
 real(r8), intent(inout) :: rarray(:)
 integer,  intent(inout) :: iarray(:)

integer :: j

real(r8) :: rtmp
integer  :: itmp

j = size(rarray) / 3

rtmp        = rarray(j)
rarray(j)   = rarray(j*2)
rarray(j*2) = rtmp

j = size(iarray) / 3

itmp        = iarray(j)
iarray(j)   = iarray(j*2)
iarray(j*2) = itmp

end subroutine reorder_pair

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine print_time(timespan, label)
 real(digits12),   intent(in) :: timespan
 character(len=*), intent(in) :: label

write(*, '(A32,F12.3)') trim(label)//' ', timespan/1000.0_digits12

end subroutine print_time

!-------------------------------------------------------------------------------

end program sort_test


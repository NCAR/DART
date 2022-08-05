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

use sort_mod,          only : sort, index_sort, insertion_sort, index_insertion_sort, index_insertion_sort_orig

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
integer, parameter :: MAXSIZE = 20480
integer, parameter :: SMALLER =  4096
integer, parameter :: SMALLE1 =  2048
integer, parameter :: SMALLE2 =  1024
integer, parameter :: TINY    =    80

type(random_seq_type) :: seq

! fill arrays with gaussian distribution of numeric values 
real(r8), parameter :: BASEVAL = 2500.0_r8
real(r8), parameter :: STDDEV  = 500.0_r8

real(r8) :: base_real_array(MAXSIZE)

! real arrays of 'm' (maxsize), 's' (smaller)
! and 't' (tiny) for sorting and timing.
real(r8) :: mr_array1(MAXSIZE), mr_array2(MAXSIZE)
real(r8) :: sr_array1(SMALLER), sr_array2(SMALLER)
real(r8) :: sr1array1(SMALLE1), sr1array2(SMALLE1)
real(r8) :: sr2array1(SMALLE2), sr2array2(SMALLE2)
real(r8) :: tr_array1(TINY),    tr_array2(TINY)

integer  :: i
character(len=64) :: content_type

character(len=*), parameter :: routine = 'sort_test'

! start of executable code

call initialize_utilities('sort_test')

print *, 'sort test.  array sizes: '
print *, '  large  = ', MAXSIZE
print *, 'smaller  = ', SMALLER
print *, 'smaller1 = ', SMALLE1
print *, 'smaller2 = ', SMALLE2
print *, '   tiny  = ', TINY
print *, ''

! fill the base array with random values for sorting.

call init_random_seq(seq)

do i=1, MAXSIZE
   base_real_array(i) = random_gaussian(seq, BASEVAL, STDDEV)
enddo


mr_array1(:) = base_real_array(:)
sr_array1(:) = base_real_array(1:SMALLER)
sr1array1(:) = base_real_array(1:SMALLE1)
sr2array1(:) = base_real_array(1:SMALLE2)
tr_array1(:) = base_real_array(1:TINY)


content_type = "random"

call run_tests(mr_array1, MAXSIZE, content_type)
call run_tests(sr_array1, SMALLER, content_type)
call run_tests(sr1array1, SMALLE1, content_type)
call run_tests(sr2array1, SMALLE2, content_type)
call run_tests(tr_array1, TINY,    content_type)


content_type = "inverted"

mr_array1 = sort(base_real_array)
do i=1, MAXSIZE
   mr_array2(i) = mr_array1(MAXSIZE - i + 1)
enddo

sr_array2(:) = mr_array2(1:SMALLER)
sr1array2(:) = mr_array2(1:SMALLE1)
sr2array2(:) = mr_array2(1:SMALLE2)
tr_array2(:) = mr_array2(1:TINY)


call run_tests(mr_array2, MAXSIZE, content_type)
call run_tests(sr_array2, SMALLER, content_type)
call run_tests(sr1array2, SMALLE1, content_type)
call run_tests(sr2array2, SMALLE2, content_type)
call run_tests(tr_array2, TINY,    content_type)


content_type = "almost sorted"

! almost sorted order arrays
mr_array2(:) = mr_array1(:)

call reorder_pairs(mr_array2)
call reorder_pairs(sr_array2)
call reorder_pairs(sr1array2)
call reorder_pairs(sr2array2)
call reorder_pairs(tr_array2)


call run_tests(mr_array2, MAXSIZE, content_type)
call run_tests(sr_array2, SMALLER, content_type)
call run_tests(sr1array2, SMALLE1, content_type)
call run_tests(sr2array2, SMALLE2, content_type)
call run_tests(tr_array2, TINY,    content_type)


write(*,*) 'end of test'

call finalize_utilities('sort_test')

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine run_tests(a1, s1, content_type)

integer, intent(in) :: s1
real(r8), intent(in) :: a1(s1)
character(len=*), intent(in) :: content_type

integer :: a2(s1)

call do_each_real_type(a1, s1, content_type)
call do_each_real_type(a1, s1, content_type)
call do_each_real_type(a1, s1, content_type)

a2 = floor(a1)

call do_each_int_type(a2, s1, content_type)
call do_each_int_type(a2, s1, content_type)
call do_each_int_type(a2, s1, content_type)

end subroutine run_tests

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! here is where templated functions would be nice.

subroutine do_each_real_type(array, asize, label)
 real(r8), intent(in) :: array(:)
 integer,  intent(in) :: asize
 character(len=*), intent(in) :: label
 
real(digits12) :: b, inc
real(r8) :: array1(asize), sorted1(asize)
integer :: indirect1(asize), indirect2(asize)


call start_mpi_timer(b)
sorted1 = sort(array)
inc = read_mpi_timer(b)
call print_time(inc, label, 'basic sort, real ', asize)
call validate(sorted1, asize)

call start_mpi_timer(b)
call index_sort(array, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, label, 'index sort, real ', asize)
call validate_indx(array, indirect1, asize)

array1 = array
indirect2 = indirect1

call start_mpi_timer(b)
call insertion_sort(array1)
inc = read_mpi_timer(b)
call print_time(inc, label, 'ins sort, real', asize)
call validate(array1, asize)

! use the sorted index (indirect1) from above with the 
! fully sorted data array (sorted1) for timing the 
! index_insertion_sort.

call start_mpi_timer(b)
call index_insertion_sort(sorted1, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, label, 'index ins sort, real', asize)
call validate_indx(sorted1, indirect1, asize)

indirect1 = indirect2

call start_mpi_timer(b)
call index_insertion_sort_orig(sorted1, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, label, 'index ins sort, real, orig', asize)
call validate_indx(sorted1, indirect1, asize)

call start_mpi_timer(b)
call index_insertion_sort(sorted1, indirect1, asize, .true.)
inc = read_mpi_timer(b)
call print_time(inc, label, 'index ins sort, real, init', asize)
call validate_indx(sorted1, indirect1, asize)

indirect1(1) = -1
call start_mpi_timer(b)
call index_insertion_sort_orig(sorted1, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, label, 'index ins sort, real, orig, init', asize)
call validate_indx(sorted1, indirect1, asize)

print *, ''

end subroutine do_each_real_type

!-------------------------------------------------------------------------------

subroutine do_each_int_type(array, asize, label)
 integer, intent(in) :: array(:)
 integer, intent(in) :: asize
 character(len=*), intent(in) :: label
 
real(digits12) :: b, inc
integer :: array1(asize), sorted1(asize)
integer :: indirect1(asize)


call start_mpi_timer(b)
sorted1 = sort(array)
inc = read_mpi_timer(b)
call print_time(inc, label, 'basic sort, integer', asize)
call ivalidate(sorted1, asize)

call start_mpi_timer(b)
call index_sort(array, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, label, 'index sort, integer', asize)
call ivalidate_indx(array, indirect1, asize)

array1 = array

call start_mpi_timer(b)
call insertion_sort(array1)
inc = read_mpi_timer(b)
call print_time(inc, label, 'ins sort, integer', asize)
call ivalidate(array1, asize)

! use the sorted index (indirect1) from above with the 
! fully sorted data array (sorted1) for timing the 
! index_insertion_sort.

call start_mpi_timer(b)
call index_insertion_sort(sorted1, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, label, 'index ins sort, integer', asize)
call ivalidate_indx(sorted1, indirect1, asize)

print *, ''

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

! swap three pairs of values:
!   at about index N*1/3  with the value at about N*2/3
!   at about index N*1/4  with the value at about N*2/4
!   at about index N*1/10 with the value at about N*2/10

subroutine reorder_pairs(rarray)
 real(r8), intent(inout) :: rarray(:)

integer :: i, j, k, indx(3)

real(r8) :: rtmp

indx(1) = 3
indx(2) = 4
indx(3) = 10

do i=1, 3
   j = size(rarray) / indx(i)
   k = j * 2

   rtmp      = rarray(j)
   rarray(j) = rarray(k)
   rarray(k) = rtmp
enddo

end subroutine reorder_pairs

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine print_time(timespan, label1, label2, asize)
 real(digits12),   intent(in) :: timespan
 character(len=*), intent(in) :: label1, label2
 integer,          intent(in) :: asize

real(digits12) :: mt, t

mt = timespan / 1000.0_digits12
t = timespan

write(*, '(A16,A34,I8,A,2(F22.5,A))') trim(label1)//';', trim(label2)//'; ', asize, ' items; ', mt, ' msec total; ', t/real(asize,digits12), ' usec/item'

end subroutine print_time

!-------------------------------------------------------------------------------

end program sort_test


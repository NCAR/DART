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

use sort_mod,          only : sort, index_sort, insertion_sort, index_insertion_sort

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
integer  :: base_int_array(MAXSIZE)

! real and integer arrays of 'm' (maxsize), 's' (smaller)
! and 't' (tiny) for sorting and timing.
real(r8) :: mr_array1(MAXSIZE), mr_array2(MAXSIZE)
real(r8) :: sr_array1(SMALLER), sr_array2(SMALLER)
real(r8) :: sr1array1(SMALLE1), sr1array2(SMALLE1)
real(r8) :: sr2array1(SMALLE2), sr2array2(SMALLE2)
real(r8) :: tr_array1(TINY),    tr_array2(TINY)
integer  :: mi_array1(MAXSIZE), mi_array2(MAXSIZE)
integer  :: si_array1(SMALLER), si_array2(SMALLER)
integer  :: si1array1(SMALLE1), si1array2(SMALLE1)
integer  :: si2array1(SMALLE2), si2array2(SMALLE2)
integer  :: ti_array1(TINY),    ti_array2(TINY)

integer  :: i, itmp
real(r8) :: tmp
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

base_int_array(:) = floor(base_real_array(:))

mr_array1(:) = base_real_array(:)
mi_array1(:) = floor(mr_array1(:))

sr_array1(:) = base_real_array(1:SMALLER)
si_array1(:) = floor(sr_array1(:))

sr1array1(:) = base_real_array(1:SMALLE1)
si1array1(:) = floor(sr1array1(:))

sr2array1(:) = base_real_array(1:SMALLE2)
si2array1(:) = floor(sr2array1(:))

tr_array1(:) = base_real_array(1:TINY)
ti_array1(:) = floor(tr_array1(:))


content_type = "random"

call do_each_real_type(mr_array1, MAXSIZE, content_type)
call do_each_real_type(sr_array1, SMALLER, content_type)
call do_each_real_type(sr1array1, SMALLE1, content_type)
call do_each_real_type(sr2array1, SMALLE2, content_type)
call do_each_real_type(tr_array1, TINY, content_type)

call do_each_int_type (mi_array1, MAXSIZE, content_type)
call do_each_int_type (si_array1, SMALLER, content_type)
call do_each_int_type (si1array1, SMALLE1, content_type)
call do_each_int_type (si2array1, SMALLE2, content_type)
call do_each_int_type (ti_array1, TINY, content_type)


content_type = "inverted"

mr_array1 = sort(base_real_array)
do i=1, MAXSIZE
   mr_array2(i) = mr_array1(MAXSIZE - i + 1)
enddo

mi_array2(:) = floor(mr_array2(:))

sr_array2(:) = mr_array2(1:SMALLER)
si_array2(:) = floor(sr_array2(:))

sr1array2(:) = mr_array2(1:SMALLE1)
si1array2(:) = floor(sr1array2(:))

sr2array2(:) = mr_array2(1:SMALLE2)
si2array2(:) = floor(sr2array2(:))

tr_array2(:) = mr_array2(1:TINY)
ti_array2(:) = floor(tr_array2(:))


call do_each_real_type(mr_array2, MAXSIZE, content_type)
call do_each_real_type(sr_array2, SMALLER, content_type)
call do_each_real_type(sr1array2, SMALLE1, content_type)
call do_each_real_type(sr2array2, SMALLE2, content_type)
call do_each_real_type(tr_array2, TINY, content_type)

call do_each_int_type (mi_array2, MAXSIZE, content_type)
call do_each_int_type (si_array2, SMALLER, content_type)
call do_each_int_type (si1array2, SMALLE1, content_type)
call do_each_int_type (si2array2, SMALLE2, content_type)
call do_each_int_type (ti_array2, TINY, content_type)


content_type = "almost sorted"

! almost sorted order arrays
mr_array2(:) = mr_array1(:)
mi_array2(:) = floor(mr_array2(:))

call reorder_pair(mr_array2, mi_array2)
call reorder_pair(sr_array2, si_array2)
call reorder_pair(sr1array2, si1array2)
call reorder_pair(sr2array2, si2array2)
call reorder_pair(tr_array2, ti_array2)


call do_each_real_type(mr_array2, MAXSIZE, content_type)
call do_each_real_type(sr_array2, SMALLER, content_type)
call do_each_real_type(sr1array2, SMALLE1, content_type)
call do_each_real_type(sr2array2, SMALLE2, content_type)
call do_each_real_type(tr_array2, TINY, content_type)

call do_each_int_type (mi_array2, MAXSIZE, content_type)
call do_each_int_type (si_array2, SMALLER, content_type)
call do_each_int_type (si1array2, SMALLE1, content_type)
call do_each_int_type (si2array2, SMALLE2, content_type)
call do_each_int_type (ti_array2, TINY, content_type)


write(*,*) 'end of test'

call finalize_utilities('sort_test')

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! here is where templated functions would be nice.

subroutine do_each_real_type(array, asize, label)
 real(r8), intent(in) :: array(:)
 integer,  intent(in) :: asize
 character(len=*), intent(in) :: label
 
real(digits12) :: b, inc
real(r8) :: array1(asize), sorted1(asize)
integer :: indirect1(asize)


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

call start_mpi_timer(b)
call insertion_sort(array1)
inc = read_mpi_timer(b)
call print_time(inc, label, 'insertion sort, real', asize)
call validate(array1, asize)

call start_mpi_timer(b)
call index_insertion_sort(array, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, label, 'index insertion sort, real', asize)
call validate_indx(array, indirect1, asize)

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
call print_time(inc, label, 'insertion sort, integer', asize)
call ivalidate(array1, asize)

call start_mpi_timer(b)
call index_insertion_sort(array, indirect1, asize)
inc = read_mpi_timer(b)
call print_time(inc, label, 'index insertion sort, integer', asize)
call ivalidate_indx(array, indirect1, asize)

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

subroutine print_time(timespan, label1, label2, asize)
 real(digits12),   intent(in) :: timespan
 character(len=*), intent(in) :: label1, label2
 integer,          intent(in) :: asize

real(digits12) :: mt, t

mt = timespan / 1000.0_digits12
t = timespan

write(*, '(A16,A32,I8,A,2(F22.5,A))') trim(label1), trim(label2)//' ', asize, ' items ', mt, ' msec total, ', t/real(asize,digits12), ' usec/item'

end subroutine print_time

!-------------------------------------------------------------------------------

end program sort_test


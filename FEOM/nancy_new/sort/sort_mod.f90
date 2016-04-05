! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module sort_mod

use     types_mod, only : r8
use utilities_mod, only : register_module

implicit none
private

public :: sort, index_sort, &
          slow_sort, slow_index_sort

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

interface sort
   module procedure rsort
   module procedure isort
end interface sort

interface index_sort
   module procedure index_sort_real
   module procedure index_sort_int
   module procedure index_sort_user
end interface


contains

!=======================================================================

subroutine initialize_module()

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!=======================================================================

function slow_sort(x)

! a basic, poorly-scaling sort for real(r8) array data.   returns a sorted
! array.  x() must be allocated or declared to be exactly the intended size; 
! all items in x() are sorted.

real(r8), intent(in) :: x(:)
real(r8)             :: slow_sort(size(x))

real(r8) :: tmp
integer  :: j, k

if ( .not. module_initialized ) call initialize_module

! what will be returned
slow_sort = x

! Do a O(N^2) sort
do j = 1, size(x) - 1
   do k = j + 1, size(x)
      ! Exchange two elements if they are in the wrong order
      if(slow_sort(j) .gt. slow_sort(k)) then
         tmp = slow_sort(k)
         slow_sort(k) = slow_sort(j)
         slow_sort(j) = tmp
      end if
   end do
end do

end function slow_sort

!=======================================================================

subroutine slow_index_sort(dist, index, num)

! slow real(r8) indexed sort. compares the values of the items in the dist()
! array (length 'num') and returns an integer index array with the order
! to traverse the dist() array.  dist() is unchanged, index() doesn't have
! to be initialized before calling this routine.

integer,  intent(in)  :: num
real(r8), intent(in)  :: dist(num)
integer,  intent(out) :: index(num)

integer :: i, j, k, itmp

if ( .not. module_initialized ) call initialize_module

! Initialize the index array to input order
do i = 1, num
   index(i) = i
end do

! Do a O(N^2) sort
do j = 1, num
   do k = 1, num - 1
      ! Exchange two elements if they are in the wrong order
      if(dist(index(k)) > dist(index(k+1))) then
         itmp = index(k)
         index(k) = index(k+1)
         index(k+1) = itmp
      end if
   end do
end do

end subroutine slow_index_sort

!=========================================================================

function rsort(x)

! Uses a heap sort algorithm on real(r8) array x(), returns sorted array.
! x() must be allocated or declared to be exactly the intended size; 
! all items in x() are sorted.

real(r8), intent(in) :: x(:)
real(r8)             :: rsort(size(x))

integer  :: num, level, ind, i, j
real(r8) :: l_val

! Get the size
num = size(x)

! Initial copy over
rsort = x

! Only one element, just send it back
if(num <= 1) return

level = num / 2 + 1
ind = num

! Keep looping until finished
do
   ! Keep going down levels until bottom
   if(level > 1) then
      level = level - 1
      l_val = rsort(level)
   else
      l_val = rsort(ind)
      rsort(ind) = rsort(1)
      ind = ind - 1
      if(ind == 1) then
         rsort(1) = l_val
         return
      end if
   end if

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         if(rsort(j) < rsort(j + 1)) j = j + 1
      end if
      if(l_val < rsort(j)) then
         rsort(i) = rsort(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      end if
      
   end do
   rsort(i) = l_val

end do

end function rsort

!=========================================================================

function isort(x)

! Uses a heap sort algorithm on integer array x(), returns sorted array.
! x() must be allocated or declared to be exactly the intended size; 
! all items in x() are sorted.

integer, intent(in) :: x(:)
integer             :: isort(size(x))

integer :: num, level, ind, i, j
integer :: l_val

! Get the size
num = size(x)

! Initial copy over
isort = x

! Only one element, just send it back
if(num <= 1) return

level = num / 2 + 1
ind = num

! Keep looping until finished
do
   ! Keep going down levels until bottom
   if(level > 1) then
      level = level - 1
      l_val = isort(level)
   else
      l_val = isort(ind)
      isort(ind) = isort(1)
      ind = ind - 1
      if(ind == 1) then
         isort(1) = l_val
         return
      end if
   end if

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         if(isort(j) < isort(j + 1)) j = j + 1
      end if
      if(l_val < isort(j)) then
         isort(i) = isort(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      end if
      
   end do
   isort(i) = l_val

end do

end function isort

!=========================================================================

subroutine index_sort_real(x, index, num)

! Uses a heap sort algorithm on x, returns array of sorted indices.
! x(num) array returned unchanged; index array doesn't need to be 
! initialized before this call.  returns an integer index array,
! usage: do i=1,num;  x(index(i)) is next item in sorted order; enddo

integer,  intent(in)  :: num
real(r8), intent(in)  :: x(num)
integer,  intent(out) :: index(num)

integer  :: ind, i, j, l_val_index, level
real(r8) :: l_val


if ( .not. module_initialized ) call initialize_module

! Initialize the index array to input order
do i = 1, num
   index(i) = i
end do

! Only one element, just send it back
if(num <= 1) return

level = num / 2 + 1
ind = num

! Keep looping until finished
do
   ! Keep going down levels until bottom
   if(level > 1) then
      level = level - 1
      l_val = x(index(level))
      l_val_index = index(level)
   else
      l_val = x(index(ind))
      l_val_index = index(ind)

      index(ind) = index(1)
      ind = ind - 1
      if(ind == 1) then
         index(1) = l_val_index
         return
      end if
   end if

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         if(x(index(j)) < x(index(j + 1))) j = j + 1
      end if
      if(l_val < x(index(j))) then
         index(i) = index(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      end if

   end do
   index(i) = l_val_index

end do

end subroutine index_sort_real

!=========================================================================

subroutine index_sort_int(x, index, num)

! Uses a heap sort algorithm on x, returns array of sorted indices.
! x(num) array returned unchanged; index array doesn't need to be 
! initialized before this call.  returns an integer index array,
! usage: do i=1,num;  x(index(i)) is next item in sorted order; enddo
! This differs from index_sort_real only in the data type of the x() array.

integer, intent(in)  :: num
integer, intent(in)  :: x(num)
integer, intent(out) :: index(num)

integer :: ind, i, j, l_val_index, level
integer :: l_val

if ( .not. module_initialized ) call initialize_module

! Initialize the index array to input order
do i = 1, num
   index(i) = i
end do

! Only one element, just send it back
if(num <= 1) return

level = num / 2 + 1
ind = num

! Keep looping until finished
do
   ! Keep going down levels until bottom
   if(level > 1) then
      level = level - 1
      l_val = x(index(level))
      l_val_index = index(level)
   else
      l_val = x(index(ind))
      l_val_index = index(ind)

      index(ind) = index(1)
      ind = ind - 1
      if(ind == 1) then
         index(1) = l_val_index
         return
      end if
   end if

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         if(x(index(j)) < x(index(j + 1))) j = j + 1
      end if
      if(l_val < x(index(j))) then
         index(i) = index(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      end if

   end do
   index(i) = l_val_index

end do

end subroutine index_sort_int

!=========================================================================

subroutine index_sort_user(index, num, comparefunc)

! Uses a heap sort algorithm on an array of indices.  returns array of 
! sorted indices based on a user-supplied sorting function.  index()
! usage: do i=1,num;  x(index(i)) is next item in sorted order; enddo
! the third argument is a user-supplied function which must take 2 integers
! as arguments and return an integer: -1 if a < b, 0 if a == b, 1 if a > b
! the function only gets the index numbers, the user code should use them 
! to compare the intended items corresponding to those index numbers and 
! return the ordering between them.  

integer,  intent(in)  :: num
integer,  intent(out) :: index(num)
interface
  integer function comparefunc(a, b)
    integer, intent(in) :: a, b
  end function comparefunc
end interface

integer :: ind, i, j, l_val_index, level
integer :: compval

if ( .not. module_initialized ) call initialize_module

! Initialize the index array to input order
do i = 1, num
   index(i) = i
end do

! Only one element, just send it back
if(num <= 1) return

level = num / 2 + 1
ind = num

! Keep looping until finished
do
   ! Keep going down levels until bottom
   if(level > 1) then
      level = level - 1
      l_val_index = index(level)
   else
      l_val_index = index(ind)

      index(ind) = index(1)
      ind = ind - 1
      if(ind == 1) then
        index(1) = l_val_index
        return
      end if
   end if

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         compval = comparefunc(index(j), index(j+1))
         if(compval < 0) j = j + 1
      end if
      compval = comparefunc(l_val_index, index(j))
      if(compval < 0) then
         index(i) = index(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      end if

   end do
   index(i) = l_val_index

end do

end subroutine index_sort_user

!=========================================================================

end module sort_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

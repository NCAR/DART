! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module sort_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use     types_mod, only : r8
use utilities_mod, only : register_module

implicit none
private

public :: slow_sort, slow_index_sort, sort, index_sort

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical, save :: module_initialized = .false.

contains


!=======================================================================

subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module

!=======================================================================

! a silly, inefficient sort for real(r8) array data

function slow_sort(x)

implicit none

real(r8), intent(in) :: x(:)
real(r8) :: slow_sort(size(x))
real(r8) :: tmp
integer j, k

if ( .not. module_initialized ) call initialize_module

! Copy to slow_sort
slow_sort = x

!  DO A SILLY N^2 SORT
do j = 1, size(x) - 1
   do k = j + 1, size(x)
!  EXCHANGE TWO ELEMENTS IF THEY'RE IN THE WRONG ORDER
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

!  real(r8) indexed sort

   implicit none

   integer num, index(num)
   real(r8) dist(num)
   integer i, j, k, itmp

if ( .not. module_initialized ) call initialize_module

!  INITIALIZE THE INDEX ARRAY TO INPUT ORDER
do i = 1, num
   index(i) = i
end do

!  DO A SILLY N^2 SORT
do j = 1, num
   do k = 1, num - 1
!  EXCHANGE TWO ELEMENTS IF THEY RE IN THE WRONG ORDER
      if(dist(index(k)) > dist(index(k+1))) then
         itmp = index(k)
         index(k) = index(k+1)
         index(k+1) = itmp
      endif
   end do
end do

!  TEMPORARY PRINT OUT TO CHECK SORT
!   do 30 j = 1, num
! 30      write(*, *) j, dist(index(j))
!   return
   end subroutine slow_index_sort


!=========================================================================

function sort(x)

! Uses a heap sort alogrithm on x, returns sorted array
implicit none

real(r8), intent(in) :: x(:)
real(r8)             :: sort(size(x))

integer  :: num, level, ind, i, j
real(r8) :: l_val

! Get the size
num = size(x)

! Initial copy over
sort = x

! Only one element, just send it back
if(num <= 1) return

level = num / 2 + 1
ind = num

! Keep looping until finished
do
   ! Keep going down levels until bottom
   if(level > 1) then
      level = level - 1
      l_val = sort(level)
   else
      l_val = sort(ind)
      sort(ind) = sort(1)
      ind = ind - 1
      if(ind == 1) then
         sort(1) = l_val
         return
      endif
   endif

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         if(sort(j) < sort(j + 1)) j = j + 1
      endif
      if(l_val < sort(j)) then
         sort(i) = sort(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      endif
      
   end do
   sort(i) = l_val

end do

end function sort


!=========================================================================


subroutine index_sort(x, index, num)

! Uses a heap sort alogrithm on x, returns array of sorted indices
implicit none

integer,  intent(in)  :: num
real(r8), intent(in)  :: x(num)
integer,  intent(out) :: index(num)

integer  :: evel, ind, i, j, l_val_index, level 
real(r8) :: l_val


if ( .not. module_initialized ) call initialize_module

!  INITIALIZE THE INDEX ARRAY TO INPUT ORDER
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
      endif
   endif

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         if(x(index(j)) < x(index(j + 1))) j = j + 1
      endif
      if(l_val < x(index(j))) then
         index(i) = index(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      endif
      
   end do
   index(i) = l_val_index

end do

end subroutine index_sort


!=========================================================================

end module sort_mod

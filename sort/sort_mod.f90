! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module sort_mod

use     types_mod, only : r8
use utilities_mod, only : register_module

implicit none
private

public :: slow_sort, slow_index_sort, sort, index_sort

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
end interface

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


function rsort(x)

! Uses a heap sort alogrithm on x, returns sorted array
implicit none

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
      endif
   endif

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         if(rsort(j) < rsort(j + 1)) j = j + 1
      endif
      if(l_val < rsort(j)) then
         rsort(i) = rsort(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      endif
      
   end do
   rsort(i) = l_val

end do

end function rsort


!=========================================================================


function isort(x)

! Uses a heap sort alogrithm on x, returns sorted array
implicit none

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
      endif
   endif

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         if(isort(j) < isort(j + 1)) j = j + 1
      endif
      if(l_val < isort(j)) then
         isort(i) = isort(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      endif
      
   end do
   isort(i) = l_val

end do

end function isort


!=========================================================================


subroutine index_sort_real(x, index, num)

! Uses a heap sort alogrithm on x, returns array of sorted indices
implicit none

integer,  intent(in)  :: num
real(r8), intent(in)  :: x(num)
integer,  intent(out) :: index(num)

integer  :: ind, i, j, l_val_index, level
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

end subroutine index_sort_real


!=========================================================================


subroutine index_sort_int(x, index, num)

! Uses a heap sort alogrithm on x (an array of integers)
!  returns array of sorted indices and the sorted array
implicit none

integer,  intent(in)  :: num
integer, intent(in)  :: x(num)
integer,  intent(out) :: index(num)

integer  :: ind, i, j, l_val_index, level
integer :: l_val

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

end subroutine index_sort_int


!=========================================================================

end module sort_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

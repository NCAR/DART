! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> A selection of sorting routines. The simplest version sorts a given array
!> of values and returns a copy of the array with the items in ascending sorted
!> order.  This works on integer and real(r8) arrays.
!> Another version returns an integer index array.  Accessing the original
!> list using these indices will traverse the array in sorted order.
!> The final version of the sort routine is similar to the index sort
!> but a compare routine is supplied by the user, allowing this code
!> to sort a list of any kinds of items.
!>

module sort_mod

use     types_mod, only : r8
use utilities_mod, only : register_module

implicit none
private

public :: sort, index_sort

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

!> single interface to sort real or integer arrays

interface sort
   module procedure rsort
   module procedure isort
end interface sort

!> single interface to return indices in sorted order,
!> leaving the original array in place.

interface index_sort
   module procedure index_sort_real
   module procedure index_sort_int
   module procedure index_sort_user
end interface

contains

!-----------------------------------------------------------------------

subroutine initialize_module()

if (module_initialized) return

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module

!-----------------------------------------------------------------------

!> Uses a heap sort algorithm on real(r8) array x(), returns sorted array.
!> x() must be allocated or declared to be exactly the intended size; 
!> all items in x() are sorted.
!>
!> @param x the (real) array to sort
!>

function rsort(x)

real(r8), intent(in) :: x(:)
real(r8)             :: rsort(size(x))

integer  :: num, level, ind, i, j
real(r8) :: l_val

! Get the array size
num = size(x)

! Copy to output
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

!-------------------------------------------------------------------------

!> Uses a heap sort algorithm on integer array x(), returns sorted array.
!> x() must be allocated or declared to be exactly the intended size; 
!> all items in x() are sorted.
!>
!> @param x the (integer) array to sort
!>

function isort(x)

integer, intent(in) :: x(:)
integer             :: isort(size(x))

integer :: num, level, ind, i, j
integer :: l_val

! Get the array size
num = size(x)

! Copy to output
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

!-------------------------------------------------------------------------

!> Uses a heap sort algorithm on x(), returns integer array of sorted indices.
!> x(num) array returned unchanged; index array doesn't need to be 
!> initialized before this call. 
!> usage: do i=1,num;  x(indx(i)) is next item in sorted order; enddo
!> This differs from index_sort_int only in the data type of the x() array.
!>
!> @param x the (real) array to sort
!> @param indx the array of integer indices to be used to traverse the input array in sorted order
!> @param num the length of x
!>
subroutine index_sort_real(x, indx, num)

integer,  intent(in)  :: num
real(r8), intent(in)  :: x(num)
integer,  intent(out) :: indx(num)

integer  :: ind, i, j, l_val_index, level
real(r8) :: l_val


if ( .not. module_initialized ) call initialize_module

! Initialize the index array to input order
do i = 1, num
   indx(i) = i
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
      l_val = x(indx(level))
      l_val_index = indx(level)
   else
      l_val = x(indx(ind))
      l_val_index = indx(ind)

      indx(ind) = indx(1)
      ind = ind - 1
      if(ind == 1) then
         indx(1) = l_val_index
         return
      endif
   endif

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         if(x(indx(j)) < x(indx(j + 1))) j = j + 1
      endif
      if(l_val < x(indx(j))) then
         indx(i) = indx(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      endif

   end do
   indx(i) = l_val_index

end do

end subroutine index_sort_real

!-------------------------------------------------------------------------

!> Uses a heap sort algorithm on x(), returns integer array of sorted indices.
!> x(num) array returned unchanged; index array doesn't need to be 
!> initialized before this call. 
!> usage: do i=1,num;  x(indx(i)) is next item in sorted order; enddo
!> This differs from index_sort_real only in the data type of the x() array.
!>
!> @param x the (integer) array to sort
!> @param indx the array of integer indices to be used to traverse the input array in sorted order
!> @param num the length of x
!>

subroutine index_sort_int(x, indx, num)

integer,  intent(in)  :: num
integer, intent(in)  :: x(num)
integer, intent(out) :: indx(num)

integer  :: ind, i, j, l_val_index, level
integer :: l_val

if ( .not. module_initialized ) call initialize_module

! Initialize the index array to input order
do i = 1, num
   indx(i) = i
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
      l_val = x(indx(level))
      l_val_index = indx(level)
   else
      l_val = x(indx(ind))
      l_val_index = indx(ind)

      indx(ind) = indx(1)
  ind = ind - 1
    if(ind == 1) then
         indx(1) = l_val_index
    return
    endif
  endif

  i = level
  j = 2 * level

  do while(j <= ind)
    if(j < ind) then
         if(x(indx(j)) < x(indx(j + 1))) j = j + 1
    endif
      if(l_val < x(indx(j))) then
         indx(i) = indx(j)
      i = j
      j = 2 * j
    else
     j = ind + 1
    endif

   end do
   indx(i) = l_val_index

end do

end subroutine index_sort_int

!-------------------------------------------------------------------------

!> Uses a heap sort algorithm to compute a sorted array of indices.  The
!> actual items being sorted are opaque to this routine; the user supplies
!> a sorting function which must have access to the original items.  This
!> routine only manipulates the integer index values which are then returned.
!>
!> The third argument is a user-supplied function which takes 2 integers (a,b) 
!> as arguments and return an integer code: 
!>     -1 if mything(a) <  mything(b)
!>      0 if mything(a) == mything(b)
!>      1 if mything(a) >  mything(b)
!> The function only gets the index numbers; the user code should use them 
!> to compare the items corresponding to those index numbers and return the 
!> ordering between them.  
!> This seems to work best when the comparefunc() is a public routine inside
!> a module so the interface specification is known to the compiler.
!>
!> usage: do i=1,num; mything(indx(i)) is next item in sorted order; enddo
!>
!> @param indx the array of integers
!> @param num the length of x
!> @param comparefunc the name of the comparison function
!>

subroutine index_sort_user(indx, num, comparefunc)

integer,  intent(in)  :: num
integer,  intent(out) :: indx(num)
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
   indx(i) = i
enddo

! Only one element, just send it back
if(num <= 1) return

level = num / 2 + 1
ind = num

! Keep looping until finished
do
   ! Keep going down levels until bottom
   if(level > 1) then
      level = level - 1
      l_val_index = indx(level)
   else
      l_val_index = indx(ind)

      indx(ind) = indx(1)
      ind = ind - 1
      if(ind == 1) then
        indx(1) = l_val_index
        return
      endif
   endif

   i = level
   j = 2 * level

   do while(j <= ind)
      if(j < ind) then
         compval = comparefunc(indx(j), indx(j+1))
         if(compval < 0) j = j + 1
      endif
      compval = comparefunc(l_val_index, indx(j))
      if(compval < 0) then
         indx(i) = indx(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      endif

   enddo
   indx(i) = l_val_index

enddo

end subroutine index_sort_user

!-------------------------------------------------------------------------

end module sort_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

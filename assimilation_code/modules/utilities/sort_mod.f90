! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> A selection of sorting routines. The simplest version sorts a given array
!> of values and returns a copy of the array with the items in ascending sorted
!> order.  This works on integer and real(r8) arrays.
!>
!> Another version returns an integer index array.  Accessing the original
!> list using these indices will traverse the array in sorted order.
!>
!> The final version of the sort routine is similar to the index sort
!> but a compare routine is supplied by the user, allowing this code
!> to sort a list of any kinds of items.
!>
!> @todo FIXME - some routines work in place and some return copies.
!> this affects the calling code in a very visible way.  can we make
!> these consistent?  or add an optional argument that says copy or no?
!> (hard to do because some are functions that return the copy; some
!> are subroutines that work in place.  would need to change the calling
!> code to be all subroutines and pass in either a single array for
!> "in place" or a src and dst array for copy.

module sort_mod

use     types_mod, only : r8, i8

implicit none
private

! what this should be:
!public :: sort, index_sort

! for now, let code indicate what sort they want. 
public :: sort, index_sort, insertion_sort, index_insertion_sort
!public :: simple_sort, simple_index_sort  

logical, save :: module_initialized = .false.

!> single interface to sort real or integer arrays

interface sort
   module procedure rsort
   module procedure isort
end interface sort

interface insertion_sort
   module procedure insertion_sort_real
   module procedure insertion_sort_int
end interface insertion_sort

! the simple sorts are here so we can time them, 
! but in all cases they are slower and a bad choice.

interface simple_sort
   module procedure simple_sort_real
   module procedure simple_sort_int
end interface simple_sort

!> single interface to return indices in sorted order,
!> leaving the original array in place.

interface index_sort
   module procedure index_sort_real_int
   module procedure index_sort_real_i8
   module procedure index_sort_int_int
!  module procedure index_sort_i8_i8
   module procedure index_sort_user
end interface

interface index_insertion_sort
   module procedure index_insertion_sort_real
   module procedure index_insertion_sort_int
end interface index_insertion_sort

interface simple_index_sort
   module procedure simple_index_sort_real
   module procedure simple_index_sort_int
end interface simple_index_sort

contains

!-----------------------------------------------------------------------

subroutine initialize_module()

if (module_initialized) return

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

integer  :: num, level, indx, i, j
real(r8) :: l_val

! Get the array size
num = size(x)

! Copy to output
rsort = x

! Only one element, just send it back
if(num <= 1) return

level = num / 2 + 1
indx = num

! Keep looping until finished
do
   ! Keep going down levels until bottom
   if(level > 1) then
      level = level - 1
      l_val = rsort(level)
   else
      l_val = rsort(indx)
      rsort(indx) = rsort(1)
      indx = indx - 1
      if(indx == 1) then
         rsort(1) = l_val
         return
      endif
   endif

   i = level
   j = 2 * level

   do while(j <= indx)
      if(j < indx) then
         if(rsort(j) < rsort(j + 1)) j = j + 1
      endif
      if(l_val < rsort(j)) then
         rsort(i) = rsort(j)
         i = j
         j = 2 * j
      else
         j = indx + 1
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

integer :: num, level, indx, i, j
integer :: l_val

! Get the array size
num = size(x)

! Copy to output
isort = x

! Only one element, just send it back
if(num <= 1) return

level = num / 2 + 1
indx = num

! Keep looping until finished
do
   ! Keep going down levels until bottom
   if(level > 1) then
      level = level - 1
      l_val = isort(level)
   else
      l_val = isort(indx)
      isort(indx) = isort(1)
      indx = indx - 1
      if(indx == 1) then
         isort(1) = l_val
         return
      endif
   endif

   i = level
   j = 2 * level

   do while(j <= indx)
      if(j < indx) then
         if(isort(j) < isort(j + 1)) j = j + 1
      endif
      if(l_val < isort(j)) then
         isort(i) = isort(j)
         i = j
         j = 2 * j
      else
         j = indx + 1
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
!> This differs from index_sort_int_int only in the data type of the x() array.
!>
!> @param x the (real) array to sort
!> @param indx the array of integer indices to be used to traverse the input array in sorted order
!> @param num the length of x

subroutine index_sort_real_int(x, indx, num)

integer,  intent(in)  :: num
real(r8), intent(in)  :: x(num)
integer,  intent(out) :: indx(num)

integer  :: ind, i, j, l_val_indx, level
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
      l_val_indx = indx(level)
   else
      l_val = x(indx(ind))
      l_val_indx = indx(ind)

      indx(ind) = indx(1)
      ind = ind - 1
      if(ind == 1) then
         indx(1) = l_val_indx
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

   indx(i) = l_val_indx

end do

end subroutine index_sort_real_int

!-------------------------------------------------------------------------

!> Uses a heap sort algorithm on x(), returns long-integer array of sorted indices.
!> x(num) array returned unchanged; index array doesn't need to be 
!> initialized before this call. 
!> usage: do i=1,num;  x(indx(i)) is next item in sorted order; enddo
!> This differs from index_sort_real_int only in the data type of the integer array
!> and the data type of the length of the integer array.
!>
!> @param x the (real) array to sort
!> @param indx the array of long-integer indices to be used to traverse the input array in sorted order
!> @param num the length of x as a long-integer

subroutine index_sort_real_i8(x, indx, num)

integer(i8),  intent(in)  :: num
real(r8),     intent(in)  :: x(num)
integer(i8),  intent(out) :: indx(num)

integer(i8) :: ind, i, j, l_val_indx, level
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
      l_val_indx = indx(level)
   else
      l_val = x(indx(ind))
      l_val_indx = indx(ind)

      indx(ind) = indx(1)
      ind = ind - 1
      if(ind == 1) then
         indx(1) = l_val_indx
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

   indx(i) = l_val_indx

end do

end subroutine index_sort_real_i8

!-------------------------------------------------------------------------

!> Uses a heap sort algorithm on x(), returns integer array of sorted indices.
!> x(num) array returned unchanged; index array doesn't need to be 
!> initialized before this call. 
!> usage: do i=1,num;  x(indx(i)) is next item in sorted order; enddo
!> This differs from index_sort_real_int only in the data type of the x() array.
!>
!> @param x the (integer) array to sort
!> @param indx the array of integer indices to be used to traverse the input array in sorted order
!> @param num the length of x

subroutine index_sort_int_int(x, indx, num)

integer, intent(in)  :: num
integer, intent(in)  :: x(num)
integer, intent(out) :: indx(num)

integer :: ind, i, j, l_val_indx, level
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
      l_val_indx = indx(level)
   else
      l_val = x(indx(ind))
      l_val_indx = indx(ind)

      indx(ind) = indx(1)
      ind = ind - 1
      if(ind == 1) then
         indx(1) = l_val_indx
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

   indx(i) = l_val_indx

end do

end subroutine index_sort_int_int

!-------------------------------------------------------------------------
!> Uses a heap sort algorithm on x(), returns long-integer array of sorted indices.
!> x(num) array returned unchanged; index array doesn't need to be 
!> initialized before this call. 
!> usage: do i=1,num;  x(indx(i)) is next item in sorted order; enddo
!> This differs from index_sort_int_int only in the data type of the x() array.
!>
!> @param x the (long-integer) array to sort
!> @param indx the array of long-integer indices to be used to traverse the input array in sorted order
!> @param num the length of x

subroutine index_sort_i8_i8(x, indx, num)

implicit none

integer(i8), intent(in)  :: num
integer(i8), intent(in)  :: x(num)
integer(i8), intent(out) :: indx(num)

integer(i8) :: ind, i, j, l_val_indx, level
integer(i8) :: l_val

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
      l_val_indx = indx(level)
   else
      l_val = x(indx(ind))
      l_val_indx = indx(ind)
      indx(ind) = indx(1)
      ind = ind - 1
      if(ind == 1) then
         indx(1) = l_val_indx
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

   indx(i) = l_val_indx

end do

end subroutine index_sort_i8_i8

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

integer :: ind, i, j, l_val_indx, level
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
      l_val_indx = indx(level)
   else
      l_val_indx = indx(ind)

      indx(ind) = indx(1)
      ind = ind - 1
      if(ind == 1) then
        indx(1) = l_val_indx
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
      compval = comparefunc(l_val_indx, indx(j))
      if(compval < 0) then
         indx(i) = indx(j)
         i = j
         j = 2 * j
      else
         j = ind + 1
      endif

   enddo
   indx(i) = l_val_indx

enddo

end subroutine index_sort_user

!-------------------------------------------------------------------------
! Insertion sort can be fast for small arrays and for arrays that are 
! already nearly sorted.  works in place.  real(r8) version.

subroutine insertion_sort_real(x)

real(r8), intent(inout) :: x(:)

integer :: i, j, num
real(r8) :: v

num = size(x)
i = 2
do while (i <= num)
  v = x(i)
  j = i - 1 
  do while (j >= 1)
     if (x(j) <= v) exit
     x(j+1) = x(j)
     j = j - 1
  end do
  x(j + 1) = v
  i = i + 1
end do

end subroutine insertion_sort_real

!-------------------------------------------------------------------------
! Insertion sort can be fast for small arrays and for arrays that are 
! already nearly sorted.  works in place.  integer version.

subroutine insertion_sort_int(x)

integer, intent(inout) :: x(:)

integer :: i, j, num, v

num = size(x)
i = 2
do while (i <= num)
  v = x(i)
  j = i - 1 
  do while (j >= 1)
     if (x(j) <= v) exit
     x(j+1) = x(j)
     j = j - 1
  end do
  x(j + 1) = v
  i = i + 1
end do

end subroutine insertion_sort_int

!-------------------------------------------------------------------------
! index version for real(r8) data.

subroutine index_insertion_sort_real(x, indx, num)

integer,  intent(in) :: num
integer, intent(out) :: indx(num)
real(r8), intent(in) :: x(num)

integer :: i, j, v_indx
real(r8) :: v

do i=1, num
   indx(i) = i
enddo

i = 2
do while (i <= num)
   v = x(indx(i))
   v_indx = indx(i)
   j = i - 1 
   do while (j >= 1)
      if (x(indx(j)) <= v) exit
      indx(j+1) = indx(j)
      j = j - 1
   end do
   indx(j + 1) = v_indx
   i = i + 1
end do

end subroutine index_insertion_sort_real

!-------------------------------------------------------------------------
! index version for integer data.

subroutine index_insertion_sort_int(x, indx, num)

integer, intent(in)  :: num
integer, intent(out) :: indx(num)
integer, intent(in)  :: x(num)

integer :: i, j, v_indx
integer :: v

do i=1, num
   indx(i) = i
enddo

i = 2
do while (i <= num)
  v = x(indx(i))
  v_indx = indx(i)
  j = i - 1 
  do while (j >= 1)
     if (x(indx(j)) <= v) exit
     indx(j+1) = indx(j)
     j = j - 1
  end do
  indx(j + 1) = v_indx
  i = i + 1
end do

end subroutine index_insertion_sort_int

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

! a simple sort for real(r8) array data.  do not use - slower in all cases.
! (is this a variant of a selection sort?). works in place.

subroutine simple_sort_real(x)

real(r8), intent(inout) :: x(:)

real(r8) :: tmp
integer :: j, k, num

if ( .not. module_initialized ) call initialize_module

num = size(x)

! O(N^2) sort 
do j = 1, num - 1
   do k = j + 1, num
      ! Exchange two elements if they're in the wrong order
      if(x(j) > x(k)) then
         tmp = x(k)
         x(k) = x(j)
         x(j) = tmp
      end if
   end do
end do

end subroutine simple_sort_real

!-------------------------------------------------------------------------

! a simple sort for integer array data.  do not use - slower in all cases.
! (is this a variant of a selection sort?). works in place.

subroutine simple_sort_int(x)

integer, intent(inout) :: x(:)

integer :: tmp
integer :: j, k, num

if ( .not. module_initialized ) call initialize_module

num = size(x)

! O(N^2) sort 
do j = 1, num - 1
   do k = j + 1, num
      ! Exchange two elements if they're in the wrong order
      if(x(j) > x(k)) then
         tmp = x(k)
         x(k) = x(j)
         x(j) = tmp
      end if
   end do
end do

end subroutine simple_sort_int

!-------------------------------------------------------------------------

!  real(r8) index sort

subroutine simple_index_sort_real(x, indx, num)

integer,  intent(in)  :: num
integer,  intent(out) :: indx(num)
real(r8), intent(in)  :: x(num)

integer :: i, j, k, itmp

if ( .not. module_initialized ) call initialize_module

! Initialize the index array to input order
do i = 1, num
   indx(i) = i
end do

! O(N^2) sort
do j = 1, num
   do k = 1, num - 1
      ! Exchange two elements if they're in the wrong order
      if(x(indx(k)) > x(indx(k+1))) then
         itmp = indx(k)
         indx(k) = indx(k+1)
         indx(k+1) = itmp
      endif
   end do
end do

end subroutine simple_index_sort_real

!-------------------------------------------------------------------------

!  integer  index sort

subroutine simple_index_sort_int(x, indx, num)

integer, intent(in)  :: num
integer, intent(out) :: indx(num)
integer, intent(in)  :: x(num)

integer :: i, j, k, itmp

if ( .not. module_initialized ) call initialize_module

! Initialize the index array to input order
do i = 1, num
   indx(i) = i
end do

! O(N^2) sort
do j = 1, num
   do k = 1, num - 1
      ! Exchange two elements if they're in the wrong order
      if(x(indx(k)) > x(indx(k+1))) then
         itmp = indx(k)
         indx(k) = indx(k+1)
         indx(k+1) = itmp
      endif
   end do
end do

end subroutine simple_index_sort_int

!-------------------------------------------------------------------------


!-------------------------------------------------------------------------

end module sort_mod


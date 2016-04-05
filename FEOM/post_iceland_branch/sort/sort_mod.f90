! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module sort_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

use     types_mod, only : r8
use utilities_mod, only : register_module

implicit none
private

public :: sort, index_sort

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
source   = "$Source$", &
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

function sort(x)

implicit none

real(r8), intent(in) :: x(:)
real(r8) :: sort(size(x))
real(r8) :: tmp
integer j, k

if ( .not. module_initialized ) call initialize_module

! Copy to sort
sort = x

!  DO A SILLY N^2 SORT
do j = 1, size(x) - 1
   do k = j + 1, size(x)
!  EXCHANGE TWO ELEMENTS IF THEY'RE IN THE WRONG ORDER
      if(sort(j) .gt. sort(k)) then
         tmp = sort(k)
         sort(k) = sort(j)
         sort(j) = tmp
      end if
   end do
end do
end function sort

!=======================================================================

   subroutine index_sort(dist, index, num)

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
   end subroutine index_sort

!=========================================================================

end module sort_mod

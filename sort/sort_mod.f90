module sort_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

private
public sort, index_sort

contains

!=======================================================================

! a silly, inefficient sort for double precision array data

function sort(x)

implicit none

double precision, intent(in) :: x(:)
double precision :: sort(size(x))
double precision :: tmp
integer j, k

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

!  double precision indexed sort

   implicit none

   integer num, index(num)
   double precision dist(num)
   integer i, j, k, itmp

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

program main

!  $Source$
!  $Revision$
!  $Date$

implicit none

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer :: i
integer, parameter :: interval = 7

write(*, *) 'set_def.out'
write(*, *) 1
write(*, *) 1800

do i = 1, 1800
   write(*, *) 100.0
   write(*, *) 1 + (i - 1) * 6
end do

end program main

! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program main

! <next three lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$

implicit none

! CVS Generated file description for error handling, do not edit
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

! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program main

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

integer :: i
integer, parameter :: interval = 105

write(*, *) 'set_def.out'
write(*, *) 1
write(*, *) 13440 / interval

do i = 1, 13440, interval
   write(*, *) 100.0
   write(*, *) i
end do

end program main

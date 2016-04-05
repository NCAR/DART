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
integer, parameter :: interval = 7

! Number of observations assumes standard 1800 points
write(*, *) 1800
! No data copies of qc
write(*, *) 0
write(*, *) 0

do i = 1, 1800
! There are more obs
   write(*, *) 0

! Input the index for the next ps variable
   write(*, *) (1 + (i - 1) * 6) * (-1)

! Time set to 0 days 0 seconds for initial sequence creation
   write(*, *) 0, 0

! Error variance is
   write(*, *) 100.0

end do

! Output the file name at the end
write(*, *) 'set_def.out'

end program main

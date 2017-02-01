! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program main

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

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

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

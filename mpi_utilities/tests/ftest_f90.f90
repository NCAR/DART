! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ftest

! very simple fortran program.  used to test compile and run
! of fortran.  if successful, will print a message and exit.

integer :: i, j

   print *, "program start"
   i = 2
   j = 3

   print *, "2 + 3 = ", i + j

   print *, "program end"

end program ftest

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

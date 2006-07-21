! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program ftest

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

! very simple fortran program.  used to test compile and run
! of fortran.  if successful, will print a message and exit.

integer :: i, j

   print *, "program start"
   i = 2
   j = 3

   print *, "2 + 3 = ", i + j

   print *, "program end"

end program ftest


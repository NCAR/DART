! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

       program ftest_f77

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$


C
C  fortran 77 test program, no MPI calls
C
       integer i, j
C
       print *, "program start"
C
       i = 2
       j = 3
C
       print *, "2 + 3 = ", i+j
C
       print *, "program end"
C    
       end


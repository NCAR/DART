! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

namelist /x/ fred, bob

program ftest_mpi

! <next few lines automatically updated by version control software, do not edit>
! $Revision$
! $Date$
! $Id$

! simple MPI fortran program.  use to test running interactively
! with MPI parallel communication libraries.  warning -- this program
! may compile without obvious errors, but at runtime, unless MPI_Init()
! returns 0 as the error code, there is a good chance the compile and
! link phase did not succeed.

! most fortran MPI implementations provide either a fortran 90 module
! which defines the interfaces to the MPI library routines, or an include
! file which defines constants.  try to use the module if it is available.

!use mpi
include "mpif.h"

! integer variables
integer :: ierror, myrank, totalprocs

   print *, "program start"

   ierror = -999
   call MPI_Init(ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Init() did not succeed, error code = ", ierror
      print *, "If error code is -999, the most likely problem is that"
      print *, "the right MPI libraries were not found at compile time."
      stop
   endif

   print *, "MPI initialized successfully"

   write ( *, nml=x)

   myrank = -1
   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Comm_rank() did not succeed, error code = ", ierror
      stop
   endif
   print *, "My MPI rank is: ", myrank

   totalprocs = -1
   call MPI_Comm_size(MPI_COMM_WORLD, totalprocs, ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Comm_size() did not succeed, error code = ", ierror
      stop
   endif
   print *, "Total MPI tasks: ", totalprocs

   ierror = -999
   call MPI_Finalize(ierror)
   if (ierror /= MPI_SUCCESS) then
      print *, "MPI_Finalize() did not succeed, error code = ", ierror
      stop
   endif

   print *, "All MPI calls succeeded, test passed."
   print *, "program end"

end program ftest_mpi


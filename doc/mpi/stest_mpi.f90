! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program fred

! test of named pipes.   proposed way for 2 tasks to communicate
! they are ready to run, and finished running.  the read should
! pause until there is something to read, so the program does not
! have to sleep and loop.
!
! this version has mpi and creates a pipe per mpi task

! <next few lines automatically updated by version control software, do not edit>
! $Revision$
! $Date$
! $Id$

include "mpif.h"


character(len=128) :: junk, pipename
integer :: iam, ierror

print *, "program start"

call MPI_Init(ierror)
if (ierror /= MPI_SUCCESS) stop

call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierror)
if (ierror /= MPI_SUCCESS) stop

print *, "i am here, and i am task ", iam
write(pipename,"(a5,i1)") "pipe.", iam
print *, "pipename = ", trim(pipename)


call system('rm -f '//trim(pipename)//'; mkfifo '//trim(pipename)//'; ls -l '//trim(pipename))
print *, "pipe created"

print *, "starting sleeping process which will write to pipe"
call system('(sleep 30; echo hello > '//trim(pipename)//')&')
print *, "sleeper launched"

print *, "opening pipe back in main program again"
open(unit=9, file=pipename, status="old", action="read", &
     form="formatted")

read(unit=9, fmt=*) junk
print *, "read junk, ready to continue"

close(unit=9)
call system ('rm -f '//trim(pipename))

print *, "pipe gone"

call MPI_Finalize(ierror)

print *, "program end"

end program fred

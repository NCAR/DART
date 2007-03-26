! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program fred

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! second version of test program; this is the main program which creates
! a named pipe (a file which act like a command line pipe).  this program 
! waits to read a line written by the 2b version of the program (which is
! started from here with a system() call.

include "mpif.h"

character(len=128) :: junk, commandline
character(len=7) :: pipename
integer :: iam, ierror

print *, "2a: program start"

call MPI_Init(ierror)
if (ierror /= MPI_SUCCESS) stop

call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierror)
if (ierror /= MPI_SUCCESS) stop

print *, "2a: i am here, and i am task ", iam
write(pipename,"(a5,i2.2)") "pipe.", iam
print *, "2a: pipename = ", pipename

call system('rm -f '//pipename//'; mkfifo '//pipename)
print *, "2a: pipe created"
call system('ls -l '//pipename)

write(commandline,"(a,i2,a)") "echo", iam, " | ./stest2b_mpi &"
print *, "2a: commandline = ", trim(commandline)

print *, "2a: starting second process which will write to pipe"
call system(commandline)
print *, "2a: second process launched"

print *, "2a: opening pipe back in main program again"
!!!open(unit=9, file=pipename, status="old", action="read", &
open(unit=9, file=pipename)

read(unit=9, fmt="(a)") junk
print *, "2a: read junk, ready to continue"

close(unit=9)
call system ('rm -f '//pipename)

print *, "2a: pipe gone"

call MPI_Finalize(ierror)

print *, "2a: program end"

end program fred

! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program fred

! second part of stest program; tests running with named pipes
! (files which act like command line pipes).  this program is expected
! to be launched by stest2a with a system call.  it writes to the pipe.
! the main program will read what this program writes.
!
! this version expects to be able to read the task number from unit 5
! (so it must be started:  echo N | ./stest2b_mpi or the equiv)

! <next few lines automatically updated by version control software, do not edit>
! $Revision$
! $Date$
! $Id$

!!include "mpif.h"

character(len=128) :: junk
character(len=7) :: pipename
integer :: iam, ierror

print *, "2b: program start"

print *, "2b: expecting to read task number from unit 5 here"
read(5, *) iam

print *, "2b: i am here, and i am task ", iam
write(pipename,"(a5,i2.2)") "pipe.", iam
print *, "2b: pipename = ", pipename

print *, "2b: opening pipe to write"
!!!open(unit=9, file=pipename, status="old", action="write", &
open(unit=9, file=pipename)

junk = "hello world"
write(unit=9, fmt="(a)") junk
print *, "2b: wrote junk, ready to continue"

close(unit=9)

print *, "2b: i am done"

print *, "2b: program end"


end program fred

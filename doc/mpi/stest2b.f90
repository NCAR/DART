! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program fred

! second part of stest program; tests running with named pipes
! (files which act like command line pipes).  this program is expected
! to be launched by stest2a with a system call.  it writes to the pipe.
! the main program will read what this program writes.

! <next few lines automatically updated by version control software, do not edit>
! $Revision$
! $Date$
! $Id$

character(len=128) :: junk

print *, "2b: i am here"

print *, "2b: opening pipe in subprogram"
open(unit=9, file="pipe", status="old", action="write", &
     form="formatted")

junk = "hello world"
write(unit=9, fmt=*) junk
print *, "2b: wrote junk, ready to continue"

close(unit=9)

print *, "2b: i am done"

end program fred

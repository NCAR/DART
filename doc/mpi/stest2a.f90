! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program fred

! second version of test program; this is the main program which creates
! a named pipe (a file which act like a command line pipe).  this program 
! waits to read a line written by the 2b version of the program (which is
! started from here with a system() call.

! <next few lines automatically updated by version control software, do not edit>
! $Revision$
! $Date$
! $Id$

character(len=128) :: junk

print *, "2a: i am here"
call system('rm -f pipe; mkfifo pipe; ls -l pipe')
print *, "2a: pipe created"

print *, "2a: starting second process which will write to pipe"
call system('./stest2b &')
print *, "2a: second process launched"

print *, "2a: opening pipe back in main program again"
open(unit=9, file="pipe", status="old", action="read", &
     form="formatted")

read(unit=9, fmt=*) junk
print *, "2a: read junk, ready to continue"

close(unit=9)
call system ('rm -f pipe')

print *, "2a: pipe gone"

end program fred

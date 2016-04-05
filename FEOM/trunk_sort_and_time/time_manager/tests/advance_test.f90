! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> initial stab at test code for the advance_time utility

program advance_test

use     utilities_mod, only : error_handler, E_ERR, E_MSG,  &
                              open_file, close_file
use mpi_utilities_mod, only : shell_execute

! the first argument is the input string to advance_time
! the second argument is what we expect the output to be

! start of code:

call doit("20150505 0", "2015050500")

call doit("20150505 0 -g", "  151334       0")

call doit("20150505 1d -c", "2015-05-06-00000")
call doit("20150505 0 -c", "2015-05-05-00000")
call doit("20150505 -1d -c", "2015-05-04-00000")

call doit("20150505 0 -w",    "2015-05-05_00:00:00")
call doit("20150505 1s -w",   "2015-05-05_00:00:01")
call doit("20150505 1d -w",   "2015-05-06_00:00:00")
call doit("20150505 -1d -w",  "2015-05-04_00:00:00")
call doit("20150505 -12h -w", "2015-05-04_12:00:00")

call doit("20150505 1s", "20150505000001")
call doit("20150505 1", "2015050501")
call doit("20151231 1d", "2016010100")

call doit("2015-05-05-00000 0",   "2015050500")
call doit("2015-05-05-00000 1",   "2015050501")
call doit("2015-05-05-00000 1m",  "201505050001")
call doit("2015-05-05-00000 1s",  "20150505000001")

call doit("2015-05-05-00000 0 -c",   "2015-05-05-00000")
call doit("2015-05-05-00000 1 -c",   "2015-05-05-03600")
call doit("2015-05-05-00000 1m -c",  "2015-05-05-00060")
call doit("2015-05-05-00000 1s -c",  "2015-05-05-00001")

!call doit("20150505 0", "2015050500")
!call doit("20150505 0", "2015050500")

contains

subroutine doit(instr, expected)

character(len=*), intent(in) :: instr
character(len=*), intent(in) :: expected

character(len=128) :: resultstr
character(len=144) :: outmsg
character(len=48) :: outbits(3)
character(len=6) :: label
integer :: retcode, iunit

retcode = shell_execute("echo "//trim(instr)//" | ./advance_time > results")

iunit = open_file("results")
read(iunit, '(a)') resultstr
call close_file(iunit)


if (resultstr /= expected) then
   !> @FIXME: do we make this stop?
   label = "FAIL: "
else
   label = "PASS: "
endif

write(outbits(1), '(A)') 'instr:    '//trim(instr)
write(outbits(2), '(A)') 'expected: '//trim(expected)
write(outbits(3), '(A)') 'got:      '//trim(resultstr)

write(outmsg, '(3A48)') adjustl(outbits(1)), adjustl(outbits(2)), adjustl(outbits(3))
call error_handler(E_MSG, label, trim(outmsg))

end subroutine doit

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

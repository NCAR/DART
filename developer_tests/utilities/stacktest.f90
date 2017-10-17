! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program stacktest

! FIXME: what is the 'right' size?  start with 1 Mb and see where
! it dies.  it seems like 128 Mb to 256 Mb should be sufficient.
integer*8, parameter :: meg = 1000 * 1000
integer*8 :: count = 1 * meg
integer   :: reps = 16

print *, 'This program tests how large an array can be put on the stack.'
print *, 'It will eventually crash.  If it crashes before printing out 256'
print *, 'you may need to use the "limit" command to increase the stacksize.'
print *, ''

! call a subroutine that has large local variables on the stack
! to make sure the stacksize limits aren't too small.  in this
! little test it's easy to see what's going on, but when you run
! filter it tends to die in random places if you've got stack
! limit problems and isn't easy to diagnose.

do n = 1, reps

   print *, ' Testing stacksize with array size of: ', count/meg, ' megabytes'

   call orange_moose(count)

   count = count * 2

enddo

! not expected to reach here

contains


subroutine orange_moose(count)
 integer*8, intent(in) :: count

! large local variable:
integer       :: godzilla(count)
character(32) :: junkbuf
integer       :: bob

! access the array contents to be sure we have enough stack space.
! in some of my tests the program did not crash on either of the
! varieties of assignment statement (disturbing), but then did
! crash with a stack exception when the values were accessed in 
! the write statement.  makes you wonder how much of memory was
! assigned 123, and whether it over-wrote some random memory location.

! should initialize the entire array, but doesn't always trigger a crash
! even if this is too large to fit on the stack.
godzilla(:)     = 123

! try explicitly accessing the first and last elements.  again, not always
! the fatal lines.
godzilla(1)     = 456
godzilla(count) = 789
 
! this seemed to really make it die - not sure why this is different
! than the assignments above -- and it's an integer array so it's not 
! a random illegal bit pattern causing a floating point trap - any bit 
! string is a valid integer value.  
write(junkbuf, "(2I4)") godzilla(1), godzilla(count)

! now take the local var and explicitly pass it as an argument.
! some compilers can be told to put large local variables on the heap.
! but i believe that subroutine arguments have to be on the stack.
bob = deeper(godzilla)
write(junkbuf, "(I4)") bob

! if you get here, it is as good a test as i can devise.
! seems like the stack size is going to be ok.

end subroutine


function deeper(bigone)
 integer, intent(inout) :: bigone(:)
 integer :: deeper

integer :: asize

! again, make sure to read the start and end of the array to be sure
! both parts are accessible.  and make sure that the optimizer
! doesn't skip executing the lines because the results aren't used.

asize = size(bigone)

deeper = bigone(1) + bigone(asize)

end function


end program 

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program stacktest

! FIXME: what is the 'right' size?  this is 16 Mb.
integer*8 :: count = 1000*1000* 16
integer   :: rc

print *, 'If this program crashes, use the "limit" command to raise stacksize.'
print *, 'If successful, a message will print out starting "Test passed:" '

! call a subroutine that has large local variables on the stack,
! and make sure the stacksize limits aren't too small.  in this
! little test it's easy to see what's going on, but when you run
! filter, it tends to die in random places if you've got stack
! limit problems, and isn't easy to diagnose.

call orange_moose()

print *, ''
print *, 'Test passed: stack size appears sufficient to run filter.'

! end main program

contains


subroutine orange_moose()

! large local variable:
integer       :: godzilla(count)
character(32) :: junkbuf

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
write(junkbuf, "(2(I4))") godzilla(1), godzilla(count)

! now take the local var and explicitly pass it as an argument.
! some compilers can be told to put large local variables on the heap.
! but i believe that subroutine arguments have to be on the stack.
call deeper(godzilla)

! if you get here, it is as good a test as i can devise.
! seems like the stack size is going to be ok.

end subroutine


subroutine deeper(bigone)
 integer, intent(inout) :: bigone(:)

integer :: asize, val

! again, make sure to read the start and end of the array to be sure
! both parts are accessible.
asize = size(bigone)

val = bigone(1) 
val = bigone(asize)

end subroutine


end program 

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

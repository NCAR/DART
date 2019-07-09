! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module parse_args_mod

use utilities_mod, only : error_handler, E_ERR

!---------------------------------------------------------------------
! parse a list of blank separated words.  intended to be used to parse
! a line of input read from a terminal/stdin, as opposed to using the
! (non-standard) fortran argc command line arg parsing.
!
! an intended use could be: 
!   % echo "a b c" | program
! or
!   % cat file
!    a b c
!   % program < file
! or
!   % program <<EOF
!   a b c
!   EOF
!
! OR:
!
!   read(unitnum, "(A256)") line
!   call get_args_from_string(line, wordcount, words)
!
! the limit on the number of words and the length of each string
! is determined by the 'words' character array that is passed in.
! it should have already been allocated by the caller.
!
! limitations:  doesn't understand escaped characters (e.g. "\ ")
! and doesn't understand quoted strings (e.g. "this string has spaces")
! it may need to handle quoted strings (with either ' or ") soon.  
!---------------------------------------------------------------------

implicit none
private

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

public :: get_args_from_string


contains

subroutine get_args_from_string(inline, argcount, argwords)
! parse a single string up into blank-separated words

 character(len=*), intent(in)  :: inline
 integer,          intent(out) :: argcount
 character(len=*), intent(out) :: argwords(:)

! in all these offsets, they are relative to 1, left hand char in string:
!  firstoff is offset to next non-blank character starting a word
!  thisoff  is offset to the current character
!  finaloff is offset of the last non-blank character in the string
! inword is a logical which toggles when inside a word or not
! maxw are the max number of words, defined by what the caller passes in
! maxl is the max length of any one word, again defined by the size of the
!  incoming array.

integer :: firstoff, finaloff, thisoff
logical :: inword
integer :: maxw, maxl
integer :: wordlen, i

character(len=len(inline)) :: argline
character(len=128) :: msgstring
character :: endword, thisc


! maxw is max number of arg 'words' allowed
! maxl is the max length of any one 'word'

maxw = size(argwords)
maxl = len(argwords(1))

argwords = ''
argcount = 0

finaloff = len_trim(inline)
if (finaloff <= 0) return

argline = inline

firstoff = 0
thisoff  = 1
inword = .false.
wordlen = 0
endword = ' '

!DEBUG print *, 'line = ', '"'//trim(argline)//'"'

LINE: do
   ! end of input?
   if (thisoff > finaloff) then
      ! if currently in a word, complete it
      if (inword) then
         argcount = argcount + 1
         if (argcount > maxw) exit LINE
         wordlen = thisoff-firstoff-1
!DEBUG print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
         if (wordlen > maxl) exit LINE
         argwords(argcount) = argline(firstoff:firstoff+wordlen)
!DEBUG print *, 'arg ', argcount, ' is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
      endif
      exit LINE
   endif

   ! next character on line
   thisc = argline(thisoff:thisoff)

!DEBUG print *, 'thisoff, finaloff, inword, endword, thisc = ', thisoff, finaloff, inword, '"'//endword//'"', ' ', '"'//thisc//'"'

   ! escaped chars - backslash prevents interpretation of next char
   if (thisc == '\') then
      ! fixme - does this really work?
      !  intent is to skip the backslash and the next char 
      !  so it isn't interpreted.
      do i=thisoff, finaloff-1
         argline(i:i) = argline(i+1:i+1)
      enddo
      argline(finaloff:finaloff) = ' '
      finaloff = finaloff-1
      thisoff = thisoff+1
      cycle LINE
   endif

   ! transition into a word?
   ! start of a space-separated or quoted string.
   if (.not. inword) then 
      if (thisc == '"' .or. thisc == "'") then
         endword = thisc
         inword = .true.
         firstoff = thisoff + 1
      else if (thisc /= ' ') then
         inword = .true.
         firstoff = thisoff
      else if (thisc == ' ') then
         endword = thisc
      endif
      thisoff = thisoff + 1
      cycle LINE
   endif

   ! transition out of a word?
   if (inword) then
      if (thisc == endword) then
         inword = .false.
         argcount = argcount + 1
         if (argcount > maxw) exit LINE
         if (thisc == ' ') endword = thisc
         wordlen = thisoff-firstoff-1
!DEBUG print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
         if (wordlen > maxl) exit LINE
         argwords(argcount) = argline(firstoff:firstoff+wordlen)
!DEBUG print *, 'arg ', argcount, ' is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
      endif
      thisoff = thisoff + 1
      cycle LINE
   endif
  
enddo LINE

if (argcount > maxw) then
   write(msgstring,*) 'more blank-separated args than max number allowed by calling code, ', maxw
   call error_handler(E_ERR,'get_args_from_string',msgstring, source,revision,revdate)
endif

if (wordlen > maxl) then
   write(msgstring,*) 'one or more args longer than max length allowed by calling code, ', maxl
   call error_handler(E_ERR,'get_args_from_string',msgstring, source,revision,revdate)
endif

!DEBUG print *, 'argcount = ', argcount


end subroutine 

end module parse_args_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

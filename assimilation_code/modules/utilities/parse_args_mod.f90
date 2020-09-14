! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! parse a list of blank separated words or name=val pairs from a string.
!
! get_args_from_string() is used to parse a single string into
! multiple blank-separated words.  returns the count of words
! and a character array of each word.
!
! get_name_val_pairs_from_string() is used to parse a line of input
! into two character arrays: the name strings and the value strings.  
!
! these routines can be used to parse text lines into tokens,
! or split up a line of input in lieu of command line args from 
! a terminal/stdin. in earlier years, fortran routines for parsing
! command line args were a non-standard extension.  alternatives
! that are completely portable are:
!
! have the program read a line from stdin (read(*,*) or read(0,*))
!
! run the program thus:
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
! to use from a file, open the file and:
!
!   read(unitnum, "(A256)") line
!   call get_args_from_string(line, wordcount, words)
!
! the limit on the number of words and the length of each string
! is determined by the 'words' character array that is passed in.
! it should have already been allocated by the caller.
!
! handles " " or ' ' quoted strings, and \ escapes the next character.
! for the name=val form, a final & sets a return flag to indicate
! there is a following line that is part of the same context.
!


module parse_args_mod

use utilities_mod, only : error_handler, E_ERR

implicit none
private

public :: get_args_from_string, get_name_val_pairs_from_string


contains

!------------------------------------------------------------------------------
! parse a single string up into blank-separated words

subroutine get_args_from_string(inline, argcount, argwords)

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
character(len=*), parameter :: routine = 'get_args_from_string'


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

!DEBUG print *, 'thisoff, finaloff, inword, endword, thisc = ', thisoff, finaloff, &
!DEBUG           inword, '"'//endword//'"', ' ', '"'//thisc//'"'

   ! escaped chars - backslash prevents interpretation of next char
   if (thisc == '\') then
      ! move the remainder of the string over, overwriting the \ and
      ! skipping the next char.
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
   call error_handler(E_ERR,routine,msgstring)
endif

if (wordlen > maxl) then
   write(msgstring,*) 'one or more args longer than max length allowed by calling code, ', maxl
   call error_handler(E_ERR,routine,msgstring)
endif

!DEBUG print *, 'argcount = ', argcount
!DEBUG do i=1, argcount
!DEBUG    print *, 'arg', i, ' is "'//trim(argwords(i))//'"'
!DEBUG enddo


end subroutine get_args_from_string

!------------------------------------------------------------------------------
! parse a single string up into blank-separated name=value words
! and return an array of names and values, plus a flag indicating
! if the line ended with an &.

subroutine get_name_val_pairs_from_string(inline, argcount, argnames, argvals, continuation)

 character(len=*), intent(in)  :: inline
 integer,          intent(out) :: argcount
 character(len=*), intent(out) :: argnames(:)
 character(len=*), intent(out) :: argvals(:)
 logical,          intent(out) :: continuation

! in all these offsets, they are relative to 1, left hand char in string:
!  firstoff is offset to next non-blank character starting a word
!  thisoff  is offset to the current character
!  finaloff is offset of the last non-blank character in the string
! inpair is a logical which toggles when inside a name=val string
! inname is a logical which is true for the name part of a pair
! inval  is a logical which is true for the value part of a pair

! probably have to allow continuation lines.  & at EOL?
! we don't have next line, so have to indicate to caller that they
! should pass us another line.

! maxw are the max number of words, defined by what the caller passes in
! maxl is the max length of any one word, again defined by the size of the
!  incoming array.

integer :: firstoff, finaloff, thisoff
logical :: inpair, inname, inval
integer :: maxw, maxl
integer :: wordlen, i

character(len=len(inline)) :: argline
character(len=128) :: msgstring
character :: endword, thisc
character(len=*), parameter :: routine = 'get_name_val_pairs_from_string'


! maxw is max number of arg 'words' allowed
! maxl is the max length of any one 'word'

maxw = size(argnames)
maxl = len(argnames(1))

if (size(argvals) /= maxw) then
   write(msgstring,*) 'array size of argnames and argvals must be the same.', size(argvals), ' != ', maxw
   call error_handler(E_ERR,routine,msgstring)
endif
if (len(argvals(1)) /= maxl) then
   write(msgstring,*) 'character length of argnames and argvals must be the same.', len(argvals(1)), ' != ', maxl
   call error_handler(E_ERR,routine,msgstring)
endif

! return vals
argcount = 0
argnames = ''
argvals = ''
continuation = .false.

! empty line?
finaloff = len_trim(inline)
if (finaloff <= 0) return

! something do to.
argline = inline

firstoff = 0
thisoff  = 1
inpair = .false.
inname = .false.
inval = .false.
wordlen = 0
endword = ' '

!DEBUG print *, 'line = ', '"'//trim(argline)//'"'

LINE: do
   ! end of input?
   if (thisoff > finaloff) then
      ! if currently in the value part of a pair, complete it
      if (inval) then
         wordlen = thisoff-firstoff-1
!DEBUG print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
         if (wordlen > maxl) exit LINE
         argvals(argcount) = argline(firstoff:firstoff+wordlen)
!DEBUG print *, 'arg ', argcount, ' is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
      else if (inname) then
         write(msgstring,*) 'name without value found at end of line'
         call error_handler(E_ERR,routine,msgstring)
      endif
      exit LINE
   endif

   ! next character on line
   thisc = argline(thisoff:thisoff)

!DEBUG print *, 'thisoff, finaloff, inpair, inname, inval, endword, thisc = ', &
!DEBUG          thisoff, finaloff, inpair, inname, inval,  '"'//endword1//'" "'//thisc//'"'

   ! escaped chars - backslash prevents interpretation of next char
   ! shift remainder of line 1 char to the left.
   if (thisc == '\') then
      do i=thisoff, finaloff-1
         argline(i:i) = argline(i+1:i+1)
      enddo
      argline(finaloff:finaloff) = ' '
      finaloff = finaloff-1
      thisoff = thisoff+1
      cycle LINE
   endif

   ! transition into a name-value pair?
   ! start of a space-separated or quoted string.
   if (.not. inpair) then 
      if (thisc == '"' .or. thisc == "'") then
         endword = thisc
         inpair = .true.
         inname = .true.
         firstoff = thisoff + 1
      else if (thisc == '&') then
         continuation = .true.
         thisoff = finaloff 
      else if (thisc == ',') then
         continue
      else if (thisc /= ' ') then
         endword = ' '
         inpair = .true.
         inname = .true.
         firstoff = thisoff
      endif
      thisoff = thisoff + 1
      cycle LINE
   endif

   ! transition between name and value, or transition out of value
   if (inpair) then
      if (inname) then
         if (thisc == endword .or. thisc == '=') then
            inname = .false.
            argcount = argcount + 1
            if (argcount > maxw) exit LINE
            if (thisc == ' ') then
               endword = thisc
            endif
            wordlen = thisoff-firstoff-1
!DEBUG print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
            if (wordlen > maxl) exit LINE
            argnames(argcount) = argline(firstoff:firstoff+wordlen)
!DEBUG print *, 'name: arg ', argcount, ' is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
            if (thisc == '=') then
               endword = ' '
               firstoff = thisoff+1
            endif
         endif
      else if (inval) then
         if (thisc == endword .or. thisc == ',') then
            inval = .false.
            inpair = .false.
            wordlen = thisoff-firstoff-1
!DEBUG print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
            if (wordlen > maxl) exit LINE
            argvals(argcount) = argline(firstoff:firstoff+wordlen)
!DEBUG print *, 'vals: arg ', argcount, ' is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
         endif
      else if (.not. inname .and. .not. inval) then
         if (thisc == '"' .or. thisc == "'") then
            endword = thisc
            inval = .true.
            firstoff = thisoff + 1
         else if (thisc == '&') then
            ! error for continue inside pair?
            write(msgstring,*) 'name without value found at end of line'
            call error_handler(E_ERR,routine,msgstring)
         else if (thisc == '=') then
            continue
         else if (thisc /= ' ') then
            endword = ' '
            inval = .true.
            firstoff = thisoff
         endif
      endif

      thisoff = thisoff + 1
      cycle LINE
   endif
  
enddo LINE

if (argcount > maxw) then
   write(msgstring,*) 'more blank-separated args than max number allowed by calling code, ', maxw
   call error_handler(E_ERR,routine,msgstring)
endif

if (wordlen > maxl) then
   write(msgstring,*) 'one or more args longer than max length allowed by calling code, ', maxl
   call error_handler(E_ERR,routine,msgstring)
endif

!DEBUG print *, 'argcount = ', argcount
!DEBUG do i=1, argcount
!DEBUG    print *, 'argname', i, ' is "'//trim(argnames(i))//'"'
!DEBUG    print *, 'argvalu', i, ' is "'//trim(argvals(i))//'"'
!DEBUG enddo


end subroutine get_name_val_pairs_from_string

!------------------------------------------------------------------------------

end module parse_args_mod


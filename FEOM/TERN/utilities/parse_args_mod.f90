! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
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

subroutine get_args_from_string(argline, argcount, argwords)
! parse a single string up into blank-separated words

 character(len=*), intent(in)  :: argline
 integer,          intent(out) :: argcount
 character(len=*), intent(out) :: argwords(:)

! in all these offsets, they are relative to 1, left hand char in string:
!  firstc is next non-blank character starting a word
!  thisc is the current character
!  finalc is the offset of the last non-blank character in the string
! inword is a logical which toggles when inside a word or not
! maxw are the max number of words, defined by what the caller passes in
! maxl is the max length of any one word, again defined by the size of the
!  in coming array.
integer :: firstc, finalc, thisc
logical :: inword
integer :: maxw, maxl
integer :: wordlen

! error handling
character(len=128) :: msgstring


! maxw is max number of arg 'words' allowed
! maxl is the max length of any one 'word'

maxw = size(argwords)
maxl = len(argwords(1))

argwords = ''
argcount = 0

finalc = len_trim(argline)
firstc = 0
thisc  = 1
inword = .false.
wordlen = 0

LINE: do
   ! end of input?
   if (thisc > finalc) then
      ! if currently in a word, complete it
      if (inword) then
         argcount = argcount + 1
         if (argcount > maxw) exit LINE
         wordlen = thisc-firstc+1
         if (wordlen > maxl) exit LINE
         argwords(argcount) = argline(firstc:thisc-1)
      endif
      exit LINE
   endif

   ! transition into a word
   if (.not. inword .and. argline(thisc:thisc) /= ' ') then
      inword = .true.
      firstc = thisc
      thisc = thisc + 1
   endif
 
   ! transition out of a word
   if (inword .and. argline(thisc:thisc) == ' ') then
      inword = .false.
      argcount = argcount + 1
      if (argcount > maxw) exit LINE
      wordlen = thisc-firstc+1
      if (wordlen > maxl) exit LINE
      argwords(argcount) = argline(firstc:thisc-1)
      thisc = thisc + 1
   endif
  
   ! no transition: multiple blanks or consective chars
   if ((.not. inword  .and. argline(thisc:thisc) == ' ') .or. &
       (inword  .and. argline(thisc:thisc) /= ' ')) then
      thisc = thisc + 1
   endif

enddo LINE

if (argcount > maxw) then
   write(msgstring,*) 'more blank-separated args than max allowed by calling code, ', maxw
   call error_handler(E_ERR,'get_args_from_string',msgstring, source,revision,revdate)
endif

if (wordlen > maxl) then
   write(msgstring,*) 'one or more args longer than max length allowed by calling code, ', maxl
   call error_handler(E_ERR,'get_args_from_string',msgstring, source,revision,revdate)
endif

end subroutine 

end module parse_args_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

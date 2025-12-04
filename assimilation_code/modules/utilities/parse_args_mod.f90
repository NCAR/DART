! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

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

use types_mod,     only : r8, i8, i4, MISSING_R8
use utilities_mod, only : error_handler, E_ERR, find_textfile_dims,    &
                          file_exist, open_file, close_file, to_upper, &
                          string_to_real, string_to_integer

implicit none
private

public :: get_args_from_string,           &
          get_name_val_pairs_from_string, &
          get_next_arg,                   &
          csv_get_obs_num,                &
          csv_get_field,                  &
          csv_file_type,                  &
          csv_field_exists,               &
          csv_print_header,               & 
          csv_get_field_index,            & 
          csv_open,                       &
          csv_close

interface csv_get_field
    module procedure csv_get_field_char
    module procedure csv_get_field_int
    module procedure csv_get_field_real
end interface csv_get_field

character(len=*), parameter :: source         = 'parse_args_mod.f90'
character(len=*), parameter :: EMPTY_ENTRY    = '_EMPTY_'            ! Used to fill in missing data in the raw file
integer,          parameter :: MAX_FIELDS_LEN = 15000
integer,          parameter :: MAX_NUM_FIELDS = 1000

! a csv file structure
type csv_file_type
   character(len=256) :: filename = '' 
   integer            :: nrows    = 0  
   integer            :: ncols    = 0  
   character          :: delim    = ','
   character(len=512) :: fields(MAX_NUM_FIELDS)
   logical            :: is_open  = .false.
end type csv_file_type


character(len=512) :: string1, string2

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
character(len=512) :: msgstring, msgstring2
character :: endword, thisc
character(len=*), parameter :: routine = 'get_args_from_string'

logical :: debug = .false. ! true to debug this routine, warning verbose


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

if (debug) print *, 'line = ', '"'//trim(argline)//'"'

NEXTCHAR: do
   ! end of input?
   if (thisoff > finaloff) then
      ! if currently in a word, complete it
      if (inword) then
         argcount = argcount + 1
         if (argcount > maxw) exit NEXTCHAR
         wordlen = thisoff-firstoff-1
if (debug) print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
         if (wordlen > maxl) exit NEXTCHAR
         argwords(argcount) = argline(firstoff:firstoff+wordlen)
if (debug) print *, 'arg ', argcount, ' is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
      endif
      exit NEXTCHAR
   endif

   ! next character on line
   thisc = argline(thisoff:thisoff)

if (debug) print *, 'thisoff, finaloff, inword, endword, thisc = ', thisoff, finaloff, &
                     inword, '"'//endword//'"', ' ', '"'//thisc//'"'

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
      cycle NEXTCHAR
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
      cycle NEXTCHAR
   endif

   ! transition out of a word?
   if (inword) then
      if (thisc == endword) then
         inword = .false.
         argcount = argcount + 1
         if (argcount > maxw) exit NEXTCHAR
         if (thisc == ' ') endword = thisc
         wordlen = thisoff-firstoff-1
if (debug)  print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
         if (wordlen > maxl) exit NEXTCHAR
         argwords(argcount) = argline(firstoff:firstoff+wordlen)
if (debug)  print *, 'arg ', argcount, ' is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
      endif
      thisoff = thisoff + 1
      cycle NEXTCHAR
   endif
  
enddo NEXTCHAR

if (argcount > maxw) then
   write(msgstring,*) 'more blank-separated args than max number allowed by calling code, ', maxw
   call error_handler(E_ERR,routine,msgstring,source)
endif

if (wordlen > maxl) then
   write(msgstring,*) 'one or more args longer than max length allowed by calling code, ', maxl
   call error_handler(E_ERR,routine,msgstring,source)
endif

if (debug) then
   print *, 'argcount = ', argcount
   do i=1, argcount
      print *, 'arg', i, ' is "'//trim(argwords(i))//'"'
   enddo
endif


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
character(len=512) :: msgstring, msgstring2
character :: endword, thisc
character(len=*), parameter :: routine = 'get_name_val_pairs_from_string'

logical :: debug = .false. ! true to debug this routine, warning verbose


! maxw is max number of arg 'words' allowed
! maxl is the max length of any one 'word'

maxw = size(argnames)
maxl = len(argnames(1))

if (size(argvals) /= maxw) then
   write(msgstring,*) 'array size of argnames and argvals must be the same.', size(argvals), ' != ', maxw
   call error_handler(E_ERR,routine,msgstring,source)
endif
if (len(argvals(1)) /= maxl) then
   write(msgstring,*) 'character length of argnames and argvals must be the same.', len(argvals(1)), ' != ', maxl
   call error_handler(E_ERR,routine,msgstring,source)
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

if (debug) print *, 'line = ', '"'//trim(argline)//'"'

NEXTCHAR: do
   ! end of input?
   if (thisoff > finaloff) then
      ! if currently in the value part of a pair, complete it
      if (inval) then
         wordlen = thisoff-firstoff-1
if(debug) print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
         if (wordlen > maxl) exit NEXTCHAR
         argvals(argcount) = argline(firstoff:firstoff+wordlen)
if(debug) print *, 'arg ', argcount, ' is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
      else if (inname) then
         write(msgstring,*) 'name without value found at end of line'
         write(msgstring2,*) 'line = ', '"'//trim(argline)//'"'
         call error_handler(E_ERR,routine,msgstring,source,text2=msgstring2)
      endif
      exit NEXTCHAR
   endif

   ! next character on line
   thisc = argline(thisoff:thisoff)

if(debug) print *, 'thisoff, finaloff, inpair, inname, inval, endword, thisc = ', &
          thisoff, finaloff, inpair, inname, inval,  '"'//endword//'" "'//thisc//'"'

   ! escaped chars - backslash prevents interpretation of next char
   ! shift remainder of line 1 char to the left.
   if (thisc == '\') then
      do i=thisoff, finaloff-1
         argline(i:i) = argline(i+1:i+1)
      enddo
      argline(finaloff:finaloff) = ' '
      finaloff = finaloff-1
      thisoff = thisoff+1
      cycle NEXTCHAR
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
      cycle NEXTCHAR
   endif

   ! transition between name and value, or transition out of value
   if (inpair) then
      if (inname) then
         if (thisc == endword .or. thisc == '=') then
            inname = .false.
            argcount = argcount + 1
            if (argcount > maxw) exit NEXTCHAR
            if (thisc == ' ') then
               endword = thisc
            endif
            wordlen = thisoff-firstoff-1
if(debug) print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
            if (wordlen > maxl) exit NEXTCHAR
            argnames(argcount) = argline(firstoff:firstoff+wordlen)
if(debug) print *, 'name: arg ', argcount, ' is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
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
if(debug) print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
            if (wordlen > maxl) exit NEXTCHAR
            argvals(argcount) = argline(firstoff:firstoff+wordlen)
if(debug) print *, 'vals: arg ', argcount, ' is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
         endif
      else if (.not. inname .and. .not. inval) then
         if (thisc == '"' .or. thisc == "'") then
            endword = thisc
            inval = .true.
            firstoff = thisoff + 1
         else if (thisc == '&') then
            ! error for continue inside pair?
            write(msgstring,*) 'name without value found at end of line'
            call error_handler(E_ERR,routine,msgstring,source)
         else if (thisc == '=') then
            continue
         else if (thisc /= ' ') then
            endword = ' '
            inval = .true.
            firstoff = thisoff
         endif
      endif

      thisoff = thisoff + 1
      cycle NEXTCHAR
   endif
  
enddo NEXTCHAR

if (argcount > maxw) then
   write(msgstring,*) 'more blank-separated args than max number allowed by calling code, ', maxw
   call error_handler(E_ERR,routine,msgstring,source)
endif

if (wordlen > maxl) then
   write(msgstring,*) 'one or more args longer than max length allowed by calling code, ', maxl
   call error_handler(E_ERR,routine,msgstring,source)
endif

if(debug) then
   print *, 'argcount = ', argcount
   do i=1, argcount
      print *, 'argname', i, ' is "'//trim(argnames(i))//'"'
      print *, 'argvalu', i, ' is "'//trim(argvals(i))//'"'
   enddo
endif


end subroutine get_name_val_pairs_from_string

!------------------------------------------------------------------------------
! parse the next blank separated token from a string.
! start parsing at inline(startoff) and return the ending offset

subroutine get_next_arg(inline, startoff, argword, endoff)

 character(len=*), intent(in)  :: inline
 integer,          intent(in)  :: startoff
 character(len=*), intent(out) :: argword
 integer,          intent(out) :: endoff

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
character(len=*), parameter :: routine = 'get_next_arg'

logical :: debug = .false. ! true to debug this routine, warning verbose


! maxw is max number of arg 'words' allowed
! maxl is the max length of any one 'word'

maxw = 1
maxl = len(argword)

argword = ''

finaloff = len_trim(inline)
if (finaloff <= 0) return

argline = inline

firstoff = 0
thisoff = startoff
inword = .false.
wordlen = 0
endword = ' '

if (debug) print *, 'line = ', '"'//trim(argline)//'"'

NEXTCHAR: do
   ! end of input?
   if (thisoff > finaloff) then
      ! if currently in a word, complete it
      if (inword) then
         wordlen = thisoff-firstoff-1
if (debug) print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
         if (wordlen > maxl) exit NEXTCHAR
         argword = argline(firstoff:firstoff+wordlen)
if (debug) print *, '1 arg is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
      endif
      exit NEXTCHAR
   endif

   ! next character on line
   thisc = argline(thisoff:thisoff)

if (debug) print *, 'thisoff, finaloff, inword, endword, thisc = ', thisoff, finaloff, &
                     inword, '"'//endword//'"', ' ', '"'//thisc//'"'

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
      cycle NEXTCHAR
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
      cycle NEXTCHAR
   endif

   ! transition out of a word?
   if (inword) then
      if (thisc == endword) then
         inword = .false.
         if (thisc == ' ') endword = thisc
         wordlen = thisoff-firstoff-1
if (debug) print *, 'thisoff, firstoff, wordlen = ', thisoff, firstoff, wordlen
         if (wordlen > maxl) exit NEXTCHAR
         argword = argline(firstoff:firstoff+wordlen)
if (debug) print *, '2 arg is ', '"'//argline(firstoff:firstoff+wordlen)//'"'
      endif
      thisoff = thisoff + 1
      if (inword) then
         cycle NEXTCHAR
      else
         exit NEXTCHAR
      endif
   endif
  
enddo NEXTCHAR

endoff = thisoff ! -1?

if (wordlen > maxl) then
   write(msgstring,*) 'arg longer than max length allowed by calling code, ', maxl
   call error_handler(E_ERR,routine,msgstring,source)
endif

if (debug) print *, '3 arg is "'//trim(argword)//'"'


end subroutine get_next_arg


!--------------------------------------------------------------
! Retrieve a string column using a cached handle.
subroutine csv_get_field_char(cf, varname, varvals, context)

type(csv_file_type), intent(in)  :: cf
character(len=*),    intent(in)  :: varname
character(len=*),    intent(out) :: varvals(:)
character(len=*),    intent(in), optional :: context

character(len=*), parameter :: routine = 'csv_get_field_char'

integer                       :: fidx
integer                       :: iunit, io, iobs, nfields
character(len=MAX_FIELDS_LEN) :: line
character(len=512)            :: entries(MAX_NUM_FIELDS)

if (.not. cf%is_open) then
   string1 = 'CSV file handle has not been initialized.'
   call error_handler(E_ERR, routine, string1, context)
   return
endif

! Locate field index in cached header
fidx = csv_find_field(cf, varname)
if (fidx < 1 .or. fidx > cf%ncols) then
   string1 = 'Requested field "'//trim(varname)//'" not found in header for file "'// &
             trim(cf%filename)//'"'
   call error_handler(E_ERR, routine, string1, context)
endif

! Open file and skip header
iunit = open_file(cf%filename, action='read', form='formatted')

read(iunit, '(A)', iostat=io) line
if (io /= 0) then
   write(string1, '(A, I0)') 'Got bad read code from input file, io = ', io
   call error_handler(E_ERR, routine, string1, context)
   call close_file(iunit)
   return
endif

! Read rows
if (size(varvals) /= cf%nrows) then
   write(string1, '(A, I0, A, I0)') 'Size mismatch: varvals has ', size(varvals),    &
        ' entries but csv handle reports nrows = ', cf%nrows
   call error_handler(E_ERR, routine, string1, context)
endif

do iobs = 1, cf%nrows
   read(iunit, '(A)', iostat=io) line
   if (io /= 0) then
      write(string1, '(A, I0, A, I0)') 'Unexpected EOF or read error after ', iobs-1, &
           ' data lines in "', trim(cf%filename)//'"'
      call error_handler(E_ERR, routine, string1, context)
   endif

   call split_fields(line, cf%delim, nfields, entries)

   ! Parse the column entry. If it's _EMPTY_ then 
   ! treat it as empty string to make it MISSING 
   ! when converted to integer or real  
   if (trim(entries(fidx)) == EMPTY_ENTRY) then
      varvals(iobs) = ''        
   else
      varvals(iobs) = trim(entries(fidx))
   endif
enddo

call close_file(iunit)

end subroutine csv_get_field_char


!--------------------------------------------------------------
! Read integer data from csv file
subroutine csv_get_field_int(cf, varname, varvals, context)

type(csv_file_type), intent(in)  :: cf
character(len=*),    intent(in)  :: varname
integer,             intent(out) :: varvals(:)
character(len=*),    intent(in), optional :: context

character(len=*), parameter :: routine = 'csv_get_field_int'

integer            :: i
character(len=512) :: strvals(size(varvals))

! Read it as string first
call csv_get_field_char(cf, varname, strvals, context)

! Then convert to int
do i = 1, size(varvals)
   varvals(i) = string_to_integer(strvals(i))
enddo

end subroutine csv_get_field_int


!--------------------------------------------------------------
! Read real data from csv file
subroutine csv_get_field_real(cf, varname, varvals, context)

type(csv_file_type), intent(in)  :: cf
character(len=*),    intent(in)  :: varname
real(r8),            intent(out) :: varvals(:)
character(len=*),    intent(in), optional :: context

character(len=*), parameter :: routine = 'csv_get_field_real'

integer            :: i
character(len=512) :: strvals(size(varvals))

call csv_get_field_char(cf, varname, strvals, context)

! To real for all column data
do i = 1, size(varvals)
   varvals(i) = string_to_real(strvals(i))
enddo

end subroutine csv_get_field_real


!----------------------------------------------------
! Get number of rows in the csv file (without header line) 
integer function csv_get_obs_num(fname, context) result(nrows)

character(len=*), intent(in) :: fname
character(len=*), intent(in), optional :: context

character(len=*), parameter  :: routine = 'csv_get_obs_num'

! Number of data entries in file
call find_textfile_dims(fname, nrows)

! Omit the header
nrows = nrows - 1

! Check if no data in file
if (nrows <= 0) then 
   write(string1, '(4A, i0)') 'Input CSV file "', trim(fname), '" contains ', &
           'no data. The number of data lines is ', nrows
   call error_handler(E_ERR, routine, string1, context)
endif

end function csv_get_obs_num


!---------------------------------------------------
! Adapt get_args_from_string after adjusting delims
subroutine split_fields(line, delim, nfields, fields)

character(len=*), intent(in)  :: line
character,        intent(in)  :: delim
integer,          intent(out) :: nfields
character(len=*), intent(out) :: fields(:)

character(len=MAX_FIELDS_LEN) :: work

! Clean the line then parse it 
work = normalize_delims(line, delim)
call get_args_from_string(work, nfields, fields)

end subroutine split_fields


!----------------------------------------------------------------------
! Replace ',' and ';' with blanks to use above parsers. 
! We also need to treat empty fields so that we don't
! collapse with the spaces and cause any column drifts. 
! This serves as a wrapper for 'get_args_from_string'
! Example: 
! A;B;;;;C;; --> A B _EMPTY_ _EMPTY_ _EMPTY_ C _EMPTY_ _EMPTY_
function normalize_delims(line, delim) result(out_line)

character(len=*), intent(in)  :: line
character,        intent(in)  :: delim

character(len=MAX_FIELDS_LEN) :: out_line
integer                       :: i, j, L, k, lee
logical                       :: prev_is_delim

! Start as with a delimiter 
out_line      = ' '
prev_is_delim = .true.

j = 1
L = len_trim(line)

lee = len(EMPTY_ENTRY)

! Go over the line 1 character at a time
do i = 1, L
   if (line(i:i) == char(13)) cycle
   if (line(i:i) == delim) then
      ! Found a delim
      if (prev_is_delim) then
         ! insert placeholder + 1 space
         out_line(j:j+lee-1) = EMPTY_ENTRY
         j = j+lee
         out_line(j:j) = ' '

         j = j+1
      else
         ! normal delimiter
         out_line(j:j) = ' '
         j = j+1
      endif
      prev_is_delim = .true.
      if (j > MAX_FIELDS_LEN - 64) exit ! prevent overflow; 64 is a small cushion
   else
      out_line(j:j) = line(i:i) 

      j = j+1
      prev_is_delim = .false.
      if (j > MAX_FIELDS_LEN - 64) exit
   endif
enddo

! Trailing empty field: line ends with a delimiter (or several)
if (L > 0 .and. line(L:L) == delim) then
   out_line(j:j+lee-1) = EMPTY_ENTRY
   j = j + lee
endif

! Trim right spaces
k = j - 1
do while (k >= 1 .and. out_line(k:k) == ' ')
   k = k - 1
enddo

if (k < 1) then
   out_line = ''
else
   out_line = out_line(1:k)
endif

end function normalize_delims


!---------------------------------------------------
! Find field index using cached header in csv_file_type.
integer function csv_find_field(cf, key) result(idx)

type(csv_file_type), intent(in) :: cf
character(len=*),    intent(in) :: key

integer           :: i
character(len=64) :: field_name, field_key

field_key = adjustl(trim(key))
call to_upper(field_key)

idx = -1

do i = 1, cf%ncols
   field_name = adjustl(trim(cf%fields(i)))
   call to_upper(field_name)
   if (field_name == field_key) then
      idx = i
      return
   endif
enddo

end function csv_find_field


!-------------------------------------------------
! Detect the delimiter of a CSV file. One should expect a ','
! but certain data files can come with a ';' so
! let's take care of both cases
function detect_delim(line) result(delim)

character(len=*), intent(in) :: line

character :: delim
integer   :: i, nsemi, ncomma, L

L = len_trim(line)

nsemi  = 0
ncomma = 0

do i = 1, L
   if (line(i:i) == ';') nsemi  = nsemi  + 1
   if (line(i:i) == ',') ncomma = ncomma + 1
enddo

! Is there a better check? 
if (nsemi >= ncomma) then
   delim = ';'
else
  delim = ','
endif

end function detect_delim


!--------------------------------------------------------------
! Help function to get field/column index of requested variable
integer function csv_get_field_index(cf, varname) result(idx)

type(csv_file_type), intent(in) :: cf
character(len=*),    intent(in) :: varname

idx = csv_find_field(cf, varname)

end function csv_get_field_index


!---------------------------------------------
! Check if a field exist in the csv file
logical function csv_field_exists(cf, varname) result(found)

type(csv_file_type), intent(in) :: cf
character(len=*),    intent(in) :: varname

integer :: idx

idx = csv_find_field(cf, varname)
found = (idx > 0)

end function csv_field_exists


!----------------------------------------
! Print the fields names in the csv file
subroutine csv_print_header(cf)

type(csv_file_type), intent(in) :: cf

integer :: i

write(*, '(2X, A, I0, A)') 'CSV Header Fields (', cf%ncols, ' columns):'
do i = 1, cf%ncols
   print '(I5,2X,A)', i, trim(cf%fields(i))
end do
end subroutine csv_print_header


!-------------------------------------------
! Open a CSV handle: cache header/dims.
! By doing so, we won't need to open the file 
! every time to read header or get dimensions. 
subroutine csv_open(fname, cf, context)

character(len=*),    intent(in)  :: fname
type(csv_file_type), intent(out) :: cf
character(len=*),    intent(in), optional :: context

character(len=*), parameter :: routine = 'csv_open'

integer                       :: io, iunit
character(len=MAX_FIELDS_LEN) :: line

! Reset
cf%filename = fname
cf%nrows    = 0
cf%ncols    = 0
cf%delim    = ','
cf%fields   = ''
cf%is_open  = .false.

! Number of rows (excluding header)
cf%nrows = csv_get_obs_num(fname, context)

! Read header and delimiter
iunit = open_file(fname, action='read', form='formatted')

read(iunit, '(A)', iostat=io) line
if (io /= 0) then
   write(string1, '(A, I0)') 'Got bad read code from input file, io = ', io
   call error_handler(E_ERR, routine, string1, context)
   call close_file(iunit)
   return
endif

cf%delim = detect_delim(line)

call split_fields(line, cf%delim, cf%ncols, cf%fields)
call close_file(iunit)

cf%is_open = .true.

end subroutine csv_open


!---------------------------------------------------
! Just invalidate/reset the csv handle.
subroutine csv_close(cf)

type(csv_file_type), intent(inout) :: cf

cf%filename = ''
cf%nrows    = 0
cf%ncols    = 0
cf%delim    = ','
cf%fields   = ''
cf%is_open  = .false.

end subroutine csv_close

end module parse_args_mod


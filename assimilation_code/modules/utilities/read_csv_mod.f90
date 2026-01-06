! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! Utility routines for reading simple ASCII, CSV-like tabular data.
!
! This module provides a lightweight interface to work with
! non-NetCDF data products. It supports files with a single
! header row and delimited fields (comma, semicolon, or whitespace).
!
! Some features:
!  - Cached CSV file handle with header, delimiter, and row & column count.
!  - Column-based access by field name (character, integer, real).
!  - Automatic delimiter detection with optional override.
!
! Files are opened once using csv_open() and closed explicitly using
! csv_close(). Internally, the file is rewound as needed to reread data.
!
! Example usage:
!
!   type(csv_file_type) :: cf
!   real(r8), allocatable :: lat(:)
!
!   call csv_open('input.csv', cf)
!   allocate(lat(csv_get_nrows(cf)))
!   call csv_get_field(cf, 'lat', lat)
!   call csv_close(cf)

module read_csv_mod

use types_mod,      only : r8
use utilities_mod,  only : error_handler, E_ERR, find_textfile_dims, &
                           open_file, close_file, to_upper,          &
                           string_to_real, string_to_integer

implicit none
private

public :: csv_get_nrows,       &
          csv_get_field,       &
          csv_field_exists,    &
          csv_print_header,    & 
          csv_get_field_index, &
          csv_file_type,       & 
          get_csv_words_from_string, &
          csv_open,            &
          csv_close

interface csv_get_field
    module procedure csv_get_field_char
    module procedure csv_get_field_int
    module procedure csv_get_field_real
end interface csv_get_field

character(len=*), parameter :: source         = 'read_csv_mod.f90'
character(len=*), parameter :: EMPTY_ENTRY    = '_EMPTY_'            ! Used to fill in missing data in the raw file
integer,          parameter :: MAX_FIELDS_LEN = 15000
integer,          parameter :: MAX_NUM_FIELDS = 1000

! a csv file structure
type csv_file_type
   private 
   character(len=256) :: filename = '' 
   integer            :: nrows    = 0  
   integer            :: ncols    = 0 
   integer            :: iunit    = -1 
   character          :: delim    = ','
   character(len=512) :: fields(MAX_NUM_FIELDS)
   logical            :: is_open  = .false.
end type csv_file_type


character(len=512) :: string1

contains


!--------------------------------------------------------------
! Retrieve a string column using a cached handle.
subroutine csv_get_field_char(cf, varname, varvals, context)

type(csv_file_type), intent(inout) :: cf
character(len=*),    intent(in)    :: varname
character(len=*),    intent(out)   :: varvals(:)
character(len=*),    intent(in), optional :: context

character(len=*), parameter :: routine = 'csv_get_field_char'

integer                       :: fidx
integer                       :: io, iobs, nfields
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

! Skip header and start reading
call csv_skip_header(cf)

! Read rows
if (size(varvals) /= cf%nrows) then
   write(string1, '(A, I0, A, I0)') 'Size mismatch: varvals has ', size(varvals),    &
        ' entries but csv handle reports nrows = ', cf%nrows
   call error_handler(E_ERR, routine, string1, context)
endif

do iobs = 1, cf%nrows
   read(cf%iunit, '(A)', iostat=io) line
   if (io /= 0) then
      write(string1, '(A, I0, A, I0)') 'Unexpected EOF or read error after ', iobs-1, &
           ' data lines in "', trim(cf%filename)//'"'
      call error_handler(E_ERR, routine, string1, context)
   endif

   call get_csv_words_from_string(line, cf%delim, nfields, entries)

   ! Parse the column entry. If it's _EMPTY_ then 
   ! treat it as empty string to make it MISSING 
   ! when converted to integer or real  
   if (trim(entries(fidx)) == EMPTY_ENTRY) then
      varvals(iobs) = ''        
   else
      varvals(iobs) = trim(entries(fidx))
   endif
enddo

end subroutine csv_get_field_char


!--------------------------------------------------------------
! Read integer data from csv file
subroutine csv_get_field_int(cf, varname, varvals, context)

type(csv_file_type), intent(inout) :: cf
character(len=*),    intent(in)    :: varname
integer,             intent(out)   :: varvals(:)
character(len=*),    intent(in), optional :: context

character(len=*), parameter :: routine = 'csv_get_field_int'

integer            :: i
character(len=512) :: strvals(size(varvals))

! Read it as string first
call csv_get_field_char(cf, varname, strvals, context)

! Then convert to int.
! string_to_integer() returns MISSING_I if the number
! can't be converted.
do i = 1, size(varvals)
   varvals(i) = string_to_integer(strvals(i))
enddo

end subroutine csv_get_field_int


!--------------------------------------------------------------
! Read real data from csv file
subroutine csv_get_field_real(cf, varname, varvals, context)

type(csv_file_type), intent(inout) :: cf
character(len=*),    intent(in)    :: varname
real(r8),            intent(out)   :: varvals(:)
character(len=*),    intent(in), optional :: context

character(len=*), parameter :: routine = 'csv_get_field_real'

integer            :: i
character(len=512) :: strvals(size(varvals))

call csv_get_field_char(cf, varname, strvals, context)

! To real for all column data.
! string_to_real() returns MISSING_R8 if the number
! can't be converted.
do i = 1, size(varvals)
   varvals(i) = string_to_real(strvals(i))
enddo

end subroutine csv_get_field_real


!----------------------------------------------------
! Get number of rows in the csv file (without header line) 
integer function csv_get_nrows_from_file(fname, context) result(nrows)

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

end function csv_get_nrows_from_file


!---------------------------------------------------
! Find field index using cached header in csv_file_type.
! Returns -1 if field not found.
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
! Delimiter detection is heuristic unless forced_delim is provided.
! Default is to return ','.
function detect_delim(line, forced_delim) result(delim)

character(len=*), intent(in) :: line
character,        intent(in), optional :: forced_delim

character :: delim
integer   :: i, nsemi, ncomma, L

! Caller knows best
if (present(forced_delim)) then
   delim = forced_delim
   return
endif

L = len_trim(line)

nsemi  = 0
ncomma = 0

do i = 1, L
   if (line(i:i) == ';') nsemi  = nsemi  + 1
   if (line(i:i) == ',') ncomma = ncomma + 1
enddo

! Is there a better check? 
! Changed to return comma unless more semicolons are found.
! If equal, default to commas which is the more common format.
! (CSV == comma separated values)
if (nsemi > ncomma) then
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
! Rewind a file and skip its header
subroutine csv_skip_header(cf, context)

type(csv_file_type), intent(inout)        :: cf
character(len=*),    intent(in), optional :: context

integer                       :: io
character(len=MAX_FIELDS_LEN) :: line

character(len=*), parameter :: routine = 'csv_skip_header'

! We just need to skip the header, so make sure 
! we're back at the top
rewind(cf%iunit, iostat=io)
if (io /= 0) then
   write(string1,'(A,I0)') 'REWIND failed, io=', io
   call error_handler(E_ERR, routine, string1, context)
   return
endif

read(cf%iunit, '(A)', iostat=io) line
if (io /= 0) then
   write(string1,'(A,I0)') 'READ(header) failed, io=', io
   call error_handler(E_ERR, routine, string1, context)
   return
endif

end subroutine csv_skip_header


!-------------------------------------------
! Accessor function to get number of rows
function csv_get_nrows(cf) result(ndata)

type(csv_file_type), intent(in) :: cf
integer :: ndata

ndata = cf%nrows

end function csv_get_nrows


!-------------------------------------------
! Open a CSV handle: cache header/dims.
! By doing so, we won't need to open the file 
! every time to read header or get dimensions. 
subroutine csv_open(fname, cf, forced_delim, context)

character(len=*),    intent(in)  :: fname
type(csv_file_type), intent(out) :: cf
character(len=*),    intent(in), optional :: forced_delim
character(len=*),    intent(in), optional :: context

character(len=*), parameter :: routine = 'csv_open'

integer                       :: io
character(len=MAX_FIELDS_LEN) :: line

! Reset
cf%filename = fname
cf%nrows    = 0
cf%ncols    = 0
cf%iunit    = -1
cf%delim    = ','
cf%fields   = ''
cf%is_open  = .false.

! Number of rows (excluding header)
cf%nrows = csv_get_nrows_from_file(fname, context)

! Read header and delimiter
cf%iunit = open_file(fname, action='read', form='formatted')

read(cf%iunit, '(A)', iostat=io) line
if (io /= 0) then
   call csv_close(cf)

   write(string1, '(A, I0)') 'Got bad read code from input file, io = ', io
   call error_handler(E_ERR, routine, string1, context)
   return
endif

! Can also enforce a specific delim as a second argument
cf%delim = detect_delim(line, forced_delim)

call get_csv_words_from_string(line, cf%delim, cf%ncols, cf%fields)

cf%is_open = .true.

end subroutine csv_open


!---------------------------------------------------
! Just invalidate/reset the csv handle.
subroutine csv_close(cf)

type(csv_file_type), intent(inout) :: cf

if (cf%iunit > 0)  call close_file(cf%iunit)

cf%filename = ''
cf%nrows    = 0 
cf%ncols    = 0 
cf%iunit    = -1
cf%delim    = ',' 
cf%fields   = ''
cf%is_open  = .false.

end subroutine csv_close


!------------------------------------------------------------------------------
! parse a single string up into delimeter-separated words
!
! This code is closely related to the parse routines
! in parse_args_mod.f90.  It is customized for CSV rules
! using separator characters that are not blanks.

subroutine get_csv_words_from_string(inline, delim, wordcount, words)

 character(len=*), intent(in)  :: inline
 character,        intent(in)  :: delim
 integer,          intent(out) :: wordcount
 character(len=*), intent(out) :: words(:)

! in all these offsets, they are relative to 1, left hand char in string:
!  firstoff is offset to next delimiter character starting a word
!  thisoff  is offset to the current character
!  finaloff is offset of the last non-delimiter character in the string
! inword is a logical which toggles when inside a word or not
! maxw are the max number of words, defined by what the caller passes in
! maxl is the max length of any one word, again defined by the size of the
!  incoming array.

integer :: firstoff, finaloff, thisoff
logical :: inword
integer :: maxw, maxl
integer :: wordlen, i

character(len=len(inline)) :: wordline
character(len=512) :: msgstring, msgstring2
character :: endword, thisc
character(len=*), parameter :: routine = 'get_csv_words_from_string'


! maxw is max number of 'words' allowed
! maxl is the max length of any one 'word'

maxw = size(words)
maxl = len(words(1))

words = ''
wordcount = 0

finaloff = len_trim(inline)
if (finaloff <= 0) return

wordline = inline

firstoff = 1
thisoff  = 1
inword = .true.
wordlen = 0
endword = delim

NEXTCHAR: do
   ! end of input?
   if (thisoff > finaloff) then
      ! if currently in a word, complete it
      if (inword) then
         wordcount = wordcount + 1
         if (wordcount > maxw) exit NEXTCHAR
         wordlen = thisoff-firstoff-1
         if (wordlen > maxl) exit NEXTCHAR
         words(wordcount) = wordline(firstoff:firstoff+wordlen)
      endif
      exit NEXTCHAR
   endif

   ! next character on line
   thisc = wordline(thisoff:thisoff)

   ! this (escape by backslash) doesn't seem to be universially supported 
   ! by CSV files but i can't see that it hurts.

   ! escaped chars - backslash prevents interpretation of next char
   if (thisc == '\') then
      ! move the remainder of the string over, overwriting the \ and
      ! skipping the next char.
      do i=thisoff, finaloff-1
         wordline(i:i) = wordline(i+1:i+1)
      enddo
      wordline(finaloff:finaloff) = ' '
      finaloff = finaloff-1
      thisoff = thisoff+1
      cycle NEXTCHAR
   endif

   ! transition into a word?  this is slightly more complex than blank
   ! separated words.  in a CSV file, the delimiters separate fields, so
   ! the first one doesn't start with one, and the last field doesn't end
   ! with one.  quotes can be used immediately after a delimiter to keep
   ! field data together.  the next char after a closing quote should be
   ! the field delimieter.

   ! start of a delimiter-separated string.
   ! unlike strings of blanks, you can't skip strings of consecutive delimiters
   ! and the first and last fields aren't enclosed by delimiters.
   if (.not. inword) then 
      if (thisc == delim) then
         inword = .true.
         thisoff = thisoff+1  ! skip delimeter
         firstoff = thisoff   ! first char of field
         endword = thisc
      else
         write(msgstring, *) "error?  not in word, next char not delimiter"
         call error_handler(E_ERR,routine,msgstring,source)
      endif 
      cycle NEXTCHAR
   endif

   ! transition out of a word?
   ! also, if the first character of a word is a quote, the
   ! word continues until the closing quote.
   if (inword) then
      ! if first char of string is a quote, skip it and mark it as
      ! the new delimiter
      if ((thisoff == firstoff) .and. &
          (thisc == '"' .or. thisc == "'")) then
         endword = thisc
         thisoff = thisoff+1
         firstoff = thisoff   ! reset start of field 
         cycle NEXTCHAR
      endif
      ! if we come to a delimiter, check for quote and remove it
      ! and reset the delimiter char
      if (thisc == endword) then
         inword = .false.
         wordlen = thisoff-firstoff-1
         if (thisc == '"' .or. thisc == "'") then
            endword = delim     
            thisoff = thisoff+1  ! skip quote
         endif
         wordcount = wordcount + 1
         if (wordcount > maxw) exit NEXTCHAR
         if (wordlen > maxl) exit NEXTCHAR
         words(wordcount) = wordline(firstoff:firstoff+wordlen)
         cycle NEXTCHAR
      endif
      thisoff = thisoff + 1  ! normal case, word contents OR end of word, skip delimiter
      cycle NEXTCHAR
   endif
  
enddo NEXTCHAR

if (wordcount > maxw) then
   write(msgstring,*) 'more delimeter-separated words than max number allowed by calling code, ', maxw
   call error_handler(E_ERR,routine,msgstring,source)
endif

if (wordlen > maxl) then
   write(msgstring,*) 'one or more words longer than max length allowed by calling code, ', maxl
   call error_handler(E_ERR,routine,msgstring,source)
endif


end subroutine get_csv_words_from_string


end module read_csv_mod

! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> This module supports the obs_impact_tool and reading in the
!> table at runtime for use during the assimilation phase, to
!> alter the impact of observations on the state vector.
!>
!> these routines build a table(ntypes, nkinds) which can be
!> indexed quickly at runtime.  this happens in the assimilation
!> loop so performance matters.

module obs_impact_mod

use      types_mod, only : r8, obstypelength, missing_r8
use  utilities_mod, only : register_module, error_handler, E_ERR, E_MSG,       &
                           open_file, close_file, nc_check, get_next_filename, &
                           find_namelist_in_file, check_namelist_read,         &
                           do_nml_file, do_nml_term, nmlfileunit, to_upper
use  obs_kind_mod        ! all kinds/types, so impossible to enumerate them here
use parse_args_mod, only : get_args_from_string

implicit none
private

public :: create_impact_table, &
          allocate_impact_table, read_impact_table, free_impact_table, &
          get_impact_table_name

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
character(len=128), parameter :: id  = "$Id$"

integer :: ios, i

! make all strings the same as the max length of the parameter limit
integer, parameter :: string_length = paramname_length

integer, parameter :: MAXLINELEN = 512
character(len=MAXLINELEN) :: readbuf
character(len=6) :: readformat = '(A512)'

! global line number for error messages
integer :: linecount
character(len=MAXLINELEN+32) :: errline

integer, parameter :: MAXWORDS = 1024
integer :: wordcount
character(len=string_length+1) :: wordarray(MAXWORDS)
integer :: itemlist(MAXWORDS)

! table of contents entry types
integer, parameter :: ENTRY_UNKNOWN  = 0
integer, parameter :: ENTRY_DARTKIND = 1
integer, parameter :: ENTRY_DARTTYPE = 2
integer, parameter :: ENTRY_GROUP    = 3
integer, parameter :: ENTRY_KEYWORD  = 4

! predefined keywords, group names, and dart type/kind names
type type_entry
  character(len=string_length) :: ename
  integer :: etype
  integer :: evalue
end type type_entry

type type_toc
   integer :: toc_count
   integer :: type_count
   integer :: kind_count
   type(type_entry), allocatable :: toc_entries(:)
end type

! this is the single, global table of contents.
! FIXME: toc should maybe be dictionary (dict?)
! but for modularity, we include it as an argument
! in all the calls, generally last except for optional args.

! user-defined groups 
type type_group_entry
  integer :: member_count
  integer, allocatable :: group_members(:)
end type type_group_entry

type type_group
  integer :: group_count
  type(type_group_entry), allocatable :: groups(:)
end type type_group

! this is the single, global table of groups.
! but for modularity, we include it as an argument
! in all the calls, generally last except for optional args.

integer :: this_state
integer, parameter :: STATE_UNDEFINED = 0
integer, parameter :: STATE_DEFGROUP = 1
integer, parameter :: STATE_DEFNOTGROUP = 2
integer, parameter :: STATE_DEFIMPACT = 3

! are we processing ALL, ALLTYPES, or ALLKINDS
integer :: this_category

! define a comment character (#, %, !) which makes the rest
! of the input line a comment
character, parameter :: commentchar = '#'

character(len=string_length) :: current_group_name = ''
integer :: current_group_id = 0
logical :: otherwork

logical, save :: module_initialized = .false.

character(len=512) :: msgstring, msgstring2, msgstring3

! namelist: input/output names, values, etc
character(len=512) :: input_filename  = ''
character(len=512) :: output_filename = ''
logical :: debug = .false.  ! .true. for more output

! namelist
namelist /obs_impact_tool_nml/  &
   input_filename,  &
   output_filename, &
   debug

contains
 
!----------------------------------------------------------------------

! TOOL:
subroutine create_impact_table()

! this is the routine that reads in the config file
! and creates an output file that's suitable for reading
! at runtime.  so this is the 'create tool' part.

type(type_toc)        :: toc
type(type_group)      :: group_toc
real(r8), allocatable :: table(:,:)
integer               :: funit

! initialization and setup
call initialize_module()

! set up space for the output table
call allocate_impact_table(table)

! read in types and kinds strings, compute sizes,
! allocate space for the toc and group types

call init_strings_and_toc(toc, group_toc)


this_state = STATE_UNDEFINED

! read in ascii control file with and build table
call error_handler(E_MSG, 'obs_impact_tool', ' reading impacts from: '//trim(input_filename))


linecount = 0
!print *, 'file contents: '

! loop until read returns nothing
funit = open_file(input_filename, action='read')
readloop: do
   ! get the next line
   read(funit, readformat, iostat=ios) readbuf
   if (ios /= 0) exit readloop
   linecount = linecount + 1

!print *, 'read input line, was: '//trim(readbuf)

   ! prepare msg in case of error
   write(errline, '(A,I6,A)') 'error on line ', linecount, '; line contents: '

   ! parse it into white-space separated words
   call get_args_from_string(readbuf, wordcount, wordarray)
   if (wordcount < 1) cycle readloop

   ! count actual number of words (e.g. stop at comment char)
   ! and note which ones are already in the TOC.  do not add
   ! new ones if they do not yet exist.
   call search_toc(wordcount, wordarray, itemlist, toc)

   ! change state based on keywords, add new group names to toc
   call state_change(wordcount, wordarray, itemlist, toc, group_toc, this_state, otherwork)
 
   if (otherwork) then
      if (this_state == STATE_DEFGROUP) then
         ! add to groups if defining groups.
         call add_to_group(wordcount, itemlist, current_group_id, group_toc, toc)
      else if (this_state == STATE_DEFNOTGROUP) then
         ! remove from groups if defining exclusion groups.
         call sub_from_group(wordcount, itemlist, current_group_id, this_category, group_toc, toc)
      else if (this_state == STATE_DEFIMPACT) then
         ! build impact table based on what is going on now
         call update_impact(table, wordcount, itemlist, wordarray, toc, group_toc)
      endif
   endif

enddo readloop
call close_file(funit)

call error_handler(E_MSG, 'obs_impact_tool', ' end of input')

! anything still missing will be set to 1.0
where (table(:,:) == missing_r8) table(:,:) = 1.0_r8

call output_impact_table(output_filename, table, toc)

if (.true.) then  ! was debug
   write(msgstring, *) 'closing file ',  trim(output_filename)
   call error_handler(E_MSG, 'obs_impact_tool', msgstring)
endif

call finalize_tables(toc, group_toc, table)

end subroutine create_impact_table

!----------------------------------------------------------------------

! TOOL & RUNTIME:
subroutine allocate_impact_table(table, ntypes, nkinds)

real(r8), allocatable, intent(out) :: table(:,:)
integer, optional,     intent(out) :: ntypes
integer, optional,     intent(out) :: nkinds

integer :: kind_count, type_count

! initialization and setup
!call initialize_module()

! output table is dimensioned (numtypes, 0:numkinds)
! space for results, and initial values
! default to 'unset'.  at the end, anything unset will be
! changed to 1.0, which means impact like usual with no change.

kind_count = get_num_raw_obs_kinds() 
type_count = get_num_obs_kinds() 

allocate(table(type_count, 0:kind_count))
table(:,:) = missing_r8

if (present(ntypes)) ntypes = type_count
if (present(nkinds)) nkinds = kind_count

end subroutine allocate_impact_table

!----------------------------------------------------------------------

! this was more like the create_impact_table routine but has since
! been replaced with a routine that reads 3 words per line, errors 
! out if not 3, and then calls the set routine immediately.  input
! lines have to be: 
!  type_name kind_name real_value
! at this stage in the process.  it identifies words which aren't 
! type or kind names without a dictionary or state machine.

! RUNTIME:
subroutine read_impact_table(sourcefile, table, allow_any_values)

character(len=*), intent(in)    :: sourcefile
real(r8),         intent(inout) :: table(:,0:)
logical,          intent(in)    :: allow_any_values

integer :: i, j
integer :: funit
character(len=obstypelength) :: typename, kindname
real(r8) :: rvalue

! assumes the table is already allocated the proper size.
! could check here for that?

call error_handler(E_MSG, 'obs_impact_table', ' reading impacts from: '//trim(sourcefile))

linecount = 0
!print *, 'file contents: '

! loop until read returns nothing
funit = open_file(sourcefile, action='read')
readloop: do
   ! get the next line
   read(funit, readformat, iostat=ios) readbuf
   if (ios /= 0) exit readloop
   linecount = linecount + 1

!print *, 'read input line, was: '//trim(readbuf)

   ! prepare msg in case of error
   write(errline, '(A,I6,A)') 'error on line ', linecount, '; line contents: '

   ! parse it into 3 white-space separated words, two strings and one real
   read(readbuf, *, iostat=ios) typename, kindname, rvalue
   if (ios /= 0) then
      call error_handler(E_ERR, 'read_impact_table', &
                        'lines must contain a specific type, a generic kind, and a real value', &
                        source, revision, revdate, text2=errline, text3=readbuf)
   endif
!print *, trim(typename)//' '//trim(kindname)//' ', rvalue

   call set_impact(table, typename, kindname, rvalue, allow_any_values)

enddo readloop
call close_file(funit)

! almost ready for use - set missing entries to 1.0_r8
where (table(:,:) == missing_r8) table(:,:) = 1.0_r8

! debug only?
do j=0, ubound(table, 2)
   do i=1, ubound(table, 1)
      if (table(i, j) /= 1.0_r8) then
         write(errline, '(2A33,F12.6)') get_obs_kind_name(i), get_raw_obs_kind_name(j), table(i,j)
         call error_handler(E_MSG, 'read_impact_table', errline)
      endif
      ! if (table(i, j) /= 1.0_r8) print *, i, j, table(i,j)
   enddo
enddo

end subroutine read_impact_table

!----------------------------------------------------------------------

! RUNTIME:
subroutine set_impact(table, typename, kindname, rvalue, allow_any_values)

real(r8),         intent(inout) :: table(:,0:)
character(len=*), intent(in)    :: typename
character(len=*), intent(in)    :: kindname
real(r8),         intent(in)    :: rvalue
logical,          intent(in)    :: allow_any_values

! change this to false to allow any value between 0 and 1
logical, save :: must_be_on_or_off = .true.

integer :: index1, index2

! expect exactly 3 tokens on these lines:
!  type  kind  real(r8)

index1 = get_obs_kind_index(typename)
index2 = get_raw_obs_kind_index(kindname)
!print *, 'in set_impact, index values are: ', index1, index2

if (index1 < 0) then
   call error_handler(E_ERR, 'obs_impact', &
                     'first word on line is unrecognized specific obs TYPE', &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif
if (index2 < 0) then
   call error_handler(E_ERR, 'obs_impact', &
                     'second word on line is unrecognized generic KIND', &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif

! options for the actual impact value here could be:  
!   1. anything goes
!   2. restrict range to between 0.0 and 1.0
!   3. allow ONLY 0.0 or 1.0

if (.not. allow_any_values) then
   if (must_be_on_or_off) then
      if (rvalue /= 0.0_r8 .and. rvalue /= 1.0_r8) then 
         call error_handler(E_ERR, 'obs_impact', &
                           'impact values must be 0 or 1', &
                           source, revision, revdate, text2=readbuf)
      endif
   else 
      if (rvalue < 0.0_r8 .or. rvalue > 1.0_r8) then 
         call error_handler(E_ERR, 'obs_impact', &
                           'impact values must be between 0 and 1, inclusive', &
                           source, revision, revdate, text2=readbuf, &
                           text3='set "allow_any_impact_value=.true." in namelist to allow')
      endif
   endif
endif

! ok, we have the first, second, and value.  update the table.

if ((table(index1, index2) /= missing_r8) .and. &
    (table(index1, index2) /= rvalue)) then
   write(msgstring, *) 'error; different impact factor already set for this pair'
   write(msgstring2, *) 'previous value was', table(index1, index2), ' new value is ', rvalue
   write(msgstring3, *) trim(typename), ' and ', trim(kindname)
   call error_handler(E_ERR, 'obs_impact', msgstring, &
                     source, revision, revdate, text2=msgstring2, text3=msgstring3)
   
endif

table(index1, index2) = rvalue

end subroutine set_impact

!----------------------------------------------------------------------

! TOOL & RUNTIME:
subroutine free_impact_table(table)

real(r8), intent(inout), allocatable :: table(:,:)

deallocate(table)

end subroutine free_impact_table

!----------------------------------------------------------------------

! TOOL & RUNTIME:
function get_impact_table_name()

character(len=512) :: get_impact_table_name

get_impact_table_name = output_filename

end function get_impact_table_name

!----------------------------------------------------------------------

! TOOL & RUNTIME:
subroutine initialize_module()

integer :: funit

if (module_initialized) return

module_initialized = .true.
call register_module(id)

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_impact_tool_nml", funit)
read(funit, nml = obs_impact_tool_nml, iostat = ios)
call check_namelist_read(funit, ios, "obs_impact_tool_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_impact_tool_nml)
if (do_nml_term()) write(     *     , nml=obs_impact_tool_nml)

if (debug) then
   call error_handler(E_MSG, 'obs_impact_tool', ' debug on')
endif

end subroutine initialize_module

!----------------------------------------------------------------------

! TOOL:
subroutine init_strings_and_toc(toc, group_toc)

type(type_toc),        intent(inout) :: toc
type(type_group),      intent(inout) :: group_toc

integer :: i, type_count, kind_count, ntoc
integer :: groupsize, membersize
character(len=obstypelength), allocatable :: knownkinds(:)
character(len=obstypelength), allocatable :: knowntypes(:)

! first, get all possible generic kinds from obs_kind_mod
! BEWARE!! the first kind number is 0!!!!!
kind_count = get_num_raw_obs_kinds() 
allocate(knownkinds(0:kind_count))
do i=0, kind_count
   knownkinds(i) = get_raw_obs_kind_name(i)
   call to_upper(knownkinds(i))
enddo

! then, get all possible specific types from obs_kind_mod
type_count = get_num_obs_kinds() 
allocate(knowntypes(type_count))
do i=1, type_count
   knowntypes(i) = get_obs_kind_name(i)
   call to_upper(knowntypes(i))
enddo


! toc size  - arbitrary size, could be larger
ntoc = (type_count + kind_count) * 4

! space for groups - both sizes completely arbitrary
groupsize = kind_count * 2
membersize = kind_count * 4
group_toc%group_count = groupsize
allocate(group_toc%groups(groupsize))
do i=1, groupsize
   group_toc%groups(i)%member_count = 0
   allocate(group_toc%groups(i)%group_members(membersize))
enddo

allocate(toc%toc_entries(ntoc))

call build_toc(kind_count, knownkinds, type_count, knowntypes, toc)

deallocate(knownkinds, knowntypes)

end subroutine init_strings_and_toc

!----------------------------------------------------------------------

! TOOL:
subroutine finalize_tables(toc, group_toc, table)

type(type_toc),        intent(inout) :: toc
type(type_group),      intent(inout) :: group_toc
real(r8), allocatable, intent(inout), optional :: table(:,:)

integer :: i

do i=1, group_toc%group_count
   deallocate(group_toc%groups(i)%group_members)
enddo
deallocate(group_toc%groups)
group_toc%group_count = 0

deallocate(toc%toc_entries)
toc%toc_count = 0
toc%kind_count = 0
toc%type_count = 0

if (present(table)) call free_impact_table(table)

end subroutine finalize_tables

!----------------------------------------------------------------------

subroutine search_toc(wordcount, wordarray, itemlist, toc)

integer,          intent(inout) :: wordcount
character(len=*), intent(in)    :: wordarray(:)
integer,          intent(out)   :: itemlist(:)
type(type_toc),   intent(inout) :: toc
 
integer :: wid
character(len=string_length) :: nextword
logical :: exists

! initialize the return values so we can just return
! when we reach the end of the list
itemlist(:) = -1

! make sure there is at least one word
if (wordcount < 1) return

do wid=1, wordcount
      
   nextword = wordarray(wid)
   call to_upper(nextword)

   if (is_comment(nextword)) then
      wordcount = wid-1
      return
   endif

   ! look up id numbers for keywords which are already in the toc
   ! but do not add new ones yet
   exists = in_toc(nextword, toc, itemlist(wid))
   
enddo

end subroutine search_toc

!----------------------------------------------------------------------

subroutine build_toc(nkinds, kindnames, ntypes, typenames, toc)

integer,          intent(in)    :: nkinds
character(len=*), intent(in)    :: kindnames(0:)
integer,          intent(in)    :: ntypes
character(len=*), intent(in)    :: typenames(:)
type(type_toc),   intent(inout) :: toc

integer :: i

! add keyword entries to toc table
call add_keywords(toc)

! a bit specialized but types and kind counts seem to be
! needed everywhere, so put them where we can access them.
call add_limits(toc, ntypes, nkinds)

! at some point we should make types and kinds disjoint
! sets of integer parameter values, because then there is 
! no confusion about using one when you meant the other.
! (e.g. do not use 1->n at all; start one at 1000 and the
! other at 5000, and subtract the offset each time before
! using.)

! the first kind has a value of 0
do i=0, nkinds
   if (kindnames(i) /= 'UNKNOWN') then
      call add_toc(kindnames(i), toc, ENTRY_DARTKIND, i)
   endif
enddo


do i=1, ntypes
   if (typenames(i) /= 'UNKNOWN') then
      call add_toc(typenames(i), toc, ENTRY_DARTTYPE, i)
   endif
enddo


end subroutine build_toc

!----------------------------------------------------------------------

function in_toc(nextword, toc, t_index)

character(len=*),  intent(in)  :: nextword
type(type_toc),    intent(in)  :: toc
integer, optional, intent(out) :: t_index
logical                        :: in_toc

integer :: i

if (present(t_index)) t_index = -1
in_toc = .false.

do i=1, toc%toc_count
   if (toc%toc_entries(i)%ename == nextword) then
      if (present(t_index)) t_index = i
      in_toc = .true.
      return
   endif
enddo

! fall through intentionally and return false (and index -1) here
!print *, 'did not find '//trim(nextword)//' in toc'

end function in_toc

!----------------------------------------------------------------------

subroutine add_toc(nextword, toc, etype, evalue)

character(len=*),  intent(in)    :: nextword
type(type_toc),    intent(inout) :: toc
integer,           intent(in)    :: etype
integer, optional, intent(in)    :: evalue

integer :: i

do i=1, toc%toc_count
   if (toc%toc_entries(i)%ename == nextword) then
      call error_handler(E_ERR, 'add_toc:', 'word already in toc: '//trim(nextword), &
                       source, revision, revdate, text2='cannot have duplicate entries')
   endif
enddo

toc%toc_entries(i)%etype = etype
toc%toc_entries(i)%ename = nextword
if (present(evalue)) then
   toc%toc_entries(i)%evalue = evalue
else
   toc%toc_entries(i)%evalue = -1
endif

toc%toc_count = toc%toc_count + 1

end subroutine add_toc

!----------------------------------------------------------------------

function get_name_by_value(evalue, etype, toc)

integer,          intent(in) :: evalue
integer,          intent(in) :: etype
type(type_toc),   intent(in) :: toc
character(len=string_length) :: get_name_by_value

integer :: i

get_name_by_value = 'UNKNOWN'

ENTRIES: do i=1, toc%toc_count
   if (toc%toc_entries(i)%etype /= etype) cycle ENTRIES

   if (toc%toc_entries(i)%evalue == evalue) then
      get_name_by_value = toc%toc_entries(i)%ename
      return
   endif
enddo ENTRIES

end function get_name_by_value

!----------------------------------------------------------------------

subroutine print_toc(toc, entry)

type(type_toc),    intent(in) :: toc
integer, optional, intent(in) :: entry

integer :: i

if (present(entry)) then
   write(*,'(A34,I5)') 'table of contents: item : ', entry
   write(*, '(A32,I3,I5)') trim(toc%toc_entries(entry)%ename), &
                              toc%toc_entries(entry)%etype, toc%toc_entries(entry)%evalue
else
   write(*,'(A34,I5)') 'table of contents: item count: ', toc%toc_count
   do i=1, toc%toc_count
      write(*, '(I5,A32,I3,I5)') i, trim(toc%toc_entries(i)%ename), &
                              toc%toc_entries(i)%etype, toc%toc_entries(i)%evalue
   enddo
endif

end subroutine print_toc

!----------------------------------------------------------------------

subroutine add_keywords(toc)

type(type_toc),   intent(inout) :: toc

! if you need to add any new keywords here is
! where you do that.  you also have to add the
! code to process what they do some place else.

call add_toc("GROUP",    toc, ENTRY_KEYWORD)
call add_toc("IMPACT",   toc, ENTRY_KEYWORD)

call add_toc("END",      toc, ENTRY_KEYWORD)
call add_toc("EXCEPT",   toc, ENTRY_KEYWORD)

call add_toc("UNKNOWN",  toc, ENTRY_KEYWORD)
call add_toc("ALL",      toc, ENTRY_KEYWORD)

call add_toc("ALLKINDS", toc, ENTRY_KEYWORD)
call add_toc("ALLTYPES", toc, ENTRY_KEYWORD)


end subroutine add_keywords

!----------------------------------------------------------------------

function is_keyword(item, name, toc)
 
integer,            intent(in) :: item
character(len=*),   intent(in) :: name
type(type_toc),     intent(in) :: toc
logical                        :: is_keyword

is_keyword = .false.   

if ((get_toc_name(item, toc) == name) .and. &
    (get_toc_type(item, toc) == ENTRY_KEYWORD)) is_keyword = .true.

! returns here if not keyword

end function is_keyword

!----------------------------------------------------------------------

subroutine add_limits(toc, ntypes, nkinds)

type(type_toc),   intent(inout) :: toc
integer,          intent(in)    :: ntypes
integer,          intent(in)    :: nkinds

toc%kind_count = nkinds
toc%type_count = ntypes

end subroutine add_limits

!----------------------------------------------------------------------

subroutine state_change(wordcount, wordarray, itemlist, toc, group_toc, this_state, otherwork)

integer,          intent(inout) :: wordcount
character(len=*), intent(in)    :: wordarray(:)
integer,          intent(in)    :: itemlist(:)
type(type_toc),   intent(inout) :: toc
type(type_group), intent(inout) :: group_toc
integer,          intent(inout) :: this_state
logical,          intent(out)   :: otherwork

! wordcount has been updated to take into account comment chars
! words have been uppercased already, and itemlist is -1 if the
! word is not yet in the toc.

! if the current line is changing state, there is no other work
! to be done.  if it is not a state line, it is adding items to
! a block or defining a relationship and needs to be processed further.
otherwork = .false.

! blank lines, comment lines, etc cannot change state
if (wordcount < 1) return

if (itemlist(1) < 0) then
   call error_handler(E_ERR, 'state_change',  &
                     'first word on line is unrecognized keyword, KIND, TYPE, or group name', &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif


select case (this_state)

   case (STATE_UNDEFINED)
      if (is_keyword(itemlist(1), 'GROUP', toc)) then
         !GROUP - start definition
         if (wordcount < 2) then
            call error_handler(E_ERR, 'state_change', &
                               'GROUP keyword must be followed by a group name', &
                               source, revision, revdate, text2=errline, text3=readbuf)
         endif
         if (itemlist(2) > 0) then
            call error_handler(E_ERR, 'state_change', &
                               'group name '//trim(wordarray(2))//' already in use', &
                                source, revision, revdate, text2=errline, text3=readbuf)
         endif

         this_state = STATE_DEFGROUP
         current_group_name = wordarray(2)
         call to_upper(current_group_name)
         current_group_id = current_group_id + 1
         call add_toc(current_group_name, toc, ENTRY_GROUP, current_group_id)
!print *, 'adding group name to toc, name,id:', trim(wordarray(2)), current_group_id
            
      else if (is_keyword(itemlist(1), 'IMPACT', toc)) then
         !IMPACT - start definition
         if (wordcount > 1) then
            call error_handler(E_ERR, 'state_change', &
                               'IMPACT keyword should not be followed by additional text', & 
                               source, revision, revdate, text2=errline, text3=readbuf)
         endif

         this_state = STATE_DEFIMPACT
          
      else if (is_keyword(itemlist(1), 'END', toc)) then
         !END - error if not inside another block
         call error_handler(E_ERR, 'state_change', &
                            'END keyword found while not inside a GROUP or IMPACT block', & 
                            source, revision, revdate, text2=errline, text3=readbuf)

      else
         ! if we are already outside a block, there should not be any other text
         call error_handler(E_ERR, 'state_change', &
                            'text found while not inside a GROUP or IMPACT block', & 
                             source, revision, revdate, text2=errline, text3=readbuf)
      endif

   case (STATE_DEFGROUP, STATE_DEFNOTGROUP)
      if (is_keyword(itemlist(1), 'GROUP', toc)) then
         !GROUP while already in a group defn - error
         call error_handler(E_ERR, 'state_change', 'GROUP keyword found inside a GROUP block', &
                            source, revision, revdate, text2=errline, text3=readbuf)
      else if (is_keyword(itemlist(1), 'IMPACT', toc)) then
         !IMPACT while already in a group defn - error
         call error_handler(E_ERR, 'state_change', 'IMPACT keyword found inside a GROUP block', &
                            source, revision, revdate, text2=errline, text3=readbuf)
      else if (is_keyword(itemlist(1), 'END', toc)) then
         !END - change state if word 2 matches, else error
         if (wordcount < 2) then
            call error_handler(E_ERR, 'state_change', &
                               'END keyword found but missing second word GROUP', & 
                               source, revision, revdate, text2=errline, text3=readbuf)
         endif
         if (.not. is_keyword(itemlist(2), 'GROUP', toc)) then
            call error_handler(E_ERR, 'state_change', &
                               'END keyword found but second word not GROUP', & 
                               source, revision, revdate, text2=errline, text3=readbuf)
         endif

         this_state = STATE_UNDEFINED

      else if ((is_keyword(itemlist(1), 'ALL',      toc)) .or. &
               (is_keyword(itemlist(1), 'ALLKINDS', toc)) .or. &
               (is_keyword(itemlist(1), 'ALLTYPES', toc))) then

         !ALL and friends - add all, and then if word 2 is EXCEPT, take away
         if (this_state == STATE_DEFNOTGROUP) then
            call error_handler(E_ERR, 'state_change', &
                               'ALL keyword found but already excluding items', &
                               source, revision, revdate, text2=errline, text3=readbuf)
         endif
       
         ! FIXME
         this_category = itemlist(1)

         ! require group be empty to start.
         if (get_group_count(current_group_id, group_toc) /= 0) then
            call error_handler(E_ERR, 'state_change', &
                               'ALL EXCEPT must be first lines inside a group', & 
                               source, revision, revdate, text2=errline, text3=readbuf)
         endif
   
         ! add all kinds, types, or types+kinds to group.
         call add_all_to_group(toc, current_group_id, itemlist(1), group_toc)

         if (wordcount >= 2) then
            if (.not. is_keyword(itemlist(2), 'EXCEPT', toc)) then
               call error_handler(E_ERR, 'state_change', &
                                  'ALL keyword found but second word not EXCEPT', & 
                                  source, revision, revdate, text2=errline, text3=readbuf)
            endif
   
!print *, 'changing state to defnotgroup'
            ! finish any kinds/groups on rest of this line.
            if (wordcount > 2) &
               call sub_from_group(wordcount-2, itemlist(3:wordcount), current_group_id, itemlist(1), group_toc, toc)

         endif

         ! subsequent lines remove kinds and groups from group.
         this_state = STATE_DEFNOTGROUP

      else
         ! anything else is fine; no state change
         otherwork = .true.

      endif

   case (STATE_DEFIMPACT)
      if (is_keyword(itemlist(1), 'GROUP', toc)) then
            !GROUP while already in a impact defn - error
            call error_handler(E_ERR, 'state_change', 'GROUP keyword found inside an IMPACT block', &
                               source, revision, revdate, text2=errline, text3=readbuf)
      else if (is_keyword(itemlist(1), 'IMPACT', toc)) then
            !IMPACT while already in a impact defn - error
            call error_handler(E_ERR, 'state_change', 'IMPACT keyword found inside an IMPACT block', &
                                source, revision, revdate, text2=errline, text3=readbuf)
      else if (is_keyword(itemlist(1), 'END', toc)) then
            !END - change state if word 2 matches, else error
            if (wordcount < 2) then
               call error_handler(E_ERR, 'state_change', &
                                  'END keyword found but missing second word IMPACT', & 
                                  source, revision, revdate, text2=errline, text3=readbuf)
            endif
            if (.not. is_keyword(itemlist(2), 'IMPACT', toc)) then
               call error_handler(E_ERR, 'state_change', &
                                  'END keyword found but second word not IMPACT', & 
                                  source, revision, revdate, text2=errline, text3=readbuf)
            endif

            this_state = STATE_UNDEFINED

      else
            ! anything else is fine; no state change
            otherwork = .true.

      endif

  case default
    call error_handler(E_ERR, 'state_change:', 'internal error: bad state value', &
                       source, revision, revdate)
end select

end subroutine state_change

!----------------------------------------------------------------------

function has_nextword(wstring, nextword)

character(len=*), intent(in)    :: wstring
character(len=*), intent(out)   :: nextword
logical                         :: has_nextword

! assume there is a next word on the line
has_nextword = .true.

! upper case the string
nextword = wstring
call to_upper(nextword)

! return false if there is no next word because of comment char
if (is_comment(nextword)) has_nextword = .false.

end function has_nextword

!----------------------------------------------------------------------

function is_comment(nexttoken)

character(len=*), intent(in) :: nexttoken
logical                      :: is_comment

if (nexttoken(1:1) == commentchar) then
   is_comment = .true.
else
   is_comment = .false.
endif

end function is_comment

!----------------------------------------------------------------------

subroutine add_to_group(wordcount, itemlist, groupid, group_toc, toc)

integer,          intent(in)    :: wordcount
integer,          intent(in)    :: itemlist(:)
integer,          intent(in)    :: groupid
type(type_group), intent(inout) :: group_toc
type(type_toc),   intent(inout) :: toc

integer :: i, j, num_mem

call check_groupindex(groupid, 'add_to_group', group_toc)

do i = 1, wordcount
   if (itemlist(i) < 0) then
      call error_handler(E_ERR, 'add_to_group', &
                         'unrecognized item cannot be added to group', &
                         source, revision, revdate, text2=errline, text3=readbuf)
   endif

   ! check for dups
   do j=1, group_toc%groups(groupid)%member_count
      if (group_toc%groups(groupid)%group_members(j) == itemlist(i)) then
         call error_handler(E_ERR, 'add_to_group', &
                            'duplicate member found; '//trim(get_toc_name(itemlist(i),toc))//' already a member of group', &
                            source, revision, revdate, text2=errline, text3=readbuf)
      endif
   enddo

   ! make sure they are kinds or types, do not allow nested groups
   if (get_toc_type(itemlist(i), toc) == ENTRY_GROUP) then
      call error_handler(E_ERR, 'add_to_group', &
                         'group members must be DART KINDS or TYPES (nested groups not allowed)', &
                         source, revision, revdate, text2=errline, text3=readbuf)
   endif


   num_mem = group_toc%groups(groupid)%member_count + 1

   group_toc%groups(groupid)%member_count = num_mem
   group_toc%groups(groupid)%group_members(num_mem) = itemlist(i)

enddo


end subroutine add_to_group

!----------------------------------------------------------------------

subroutine add_all_to_group(toc, groupid, category, group_toc)

type(type_toc),   intent(in)    :: toc
integer,          intent(in)    :: groupid
integer,          intent(in)    :: category
type(type_group), intent(inout) :: group_toc

integer :: i, itemtype, num_mem
logical :: dokinds, dotypes

dokinds = .false.
dotypes = .false.

if (is_keyword(category, 'ALL', toc)) then
   dokinds = .true.
   dotypes = .true.
else if (is_keyword(category, 'ALLKINDS', toc)) then
   dokinds = .true.
else if (is_keyword(category, 'ALLTYPES', toc)) then
   dotypes = .true.
else
   call error_handler(E_ERR, 'add_all_to_group', &
                      'internal error: ALLxxx keyword unrecognized', &
                      source, revision, revdate)
endif

call check_groupindex(groupid, 'add_all_to_group', group_toc)

! for all kinds and/or types in toc, add them as requested
toc_loop: do i = 1, toc%toc_count

   itemtype = get_toc_type(i, toc)
   if ((itemtype == ENTRY_DARTKIND .and. dokinds) .or. &
       (itemtype == ENTRY_DARTTYPE .and. dotypes)) then

      num_mem = group_toc%groups(groupid)%member_count + 1

      !print *, 'adding index i, with name,val ', trim(get_toc_name(i, toc)), get_toc_value(i, toc)
      group_toc%groups(groupid)%member_count = num_mem
      group_toc%groups(groupid)%group_members(num_mem) = i
   endif

enddo toc_loop

end subroutine add_all_to_group

!----------------------------------------------------------------------

subroutine sub_from_group(wordcount, itemlist, groupid, category, group_toc, toc)

integer,          intent(in)    :: wordcount
integer,          intent(in)    :: itemlist(:)
integer,          intent(in)    :: groupid
integer,          intent(in)    :: category
type(type_group), intent(inout) :: group_toc
type(type_toc),   intent(inout) :: toc

integer :: i, j, itemtype
logical :: isgroup
integer :: itemval, itemcount, subgroupid
logical :: dokinds, dotypes

dokinds = .false.
dotypes = .false.

if (is_keyword(category, 'ALL', toc)) then
   dokinds = .true.
   dotypes = .true.
else if (is_keyword(category, 'ALLKINDS', toc)) then
   dokinds = .true.
else if (is_keyword(category, 'ALLTYPES', toc)) then
   dotypes = .true.
else
   call error_handler(E_ERR, 'sub_from_group', &
                      'internal error: ALLxxx keyword unrecognized', &
                      source, revision, revdate)
endif

call check_groupindex(groupid, 'sub_from_group', group_toc)

do i = 1, wordcount
   if (itemlist(i) < 0) then
      call error_handler(E_ERR, 'sub_from_group', &
                         'EXCEPT items must be known KINDS, TYPES, or groups', &
                         source, revision, revdate, text2=errline, text3=readbuf)

   endif

   ! see if itemlist(i) is an item or a group.  if item, remove it.
   ! if group, iterate through the group and remove the member items.
   ! for now, do NOT allow nested groups, although at some point 
   ! i guess we could start recursing.  but it makes the
   ! error checking a nightmare, and hard to detect logical loops.
   ! for now, just say no.

   call iteminfo(itemlist(i), itemval, isgroup, itemcount, toc, group_toc, kindonly = .false.)

   if (isgroup) subgroupid = itemval

   do j=1, itemcount
      
      if (isgroup) then
  
         itemval = get_group_member_by_index(subgroupid, j, group_toc)

         itemtype = get_toc_type(itemval, toc)
         if (.not. (itemtype == ENTRY_DARTKIND .and. dokinds) .and. &
             .not. (itemtype == ENTRY_DARTTYPE .and. dotypes)) then
            call error_handler(E_ERR, 'sub_from_group', &
                               'error using EXCEPT; item '//trim(get_toc_name(itemval, toc))//' not a group member', &
                               source, revision, revdate, text2=errline, text3=readbuf)
         endif

      else
         itemval = itemlist(i)
         if (.not. is_group_member(groupid, itemval, group_toc)) then
            call error_handler(E_ERR, 'sub_from_group', &
                               'error using EXCEPT; item '//trim(get_toc_name(itemval, toc))//' not a member of group', &
                               source, revision, revdate, text2=errline, text3=readbuf)
         endif
      endif
   
      call del_group_member(groupid, itemval, group_toc)

   enddo

enddo

end subroutine sub_from_group

!----------------------------------------------------------------------

function is_group_member(groupid, itemnum, group_toc)

integer,          intent(in) :: groupid
integer,          intent(in) :: itemnum
type(type_group), intent(in) :: group_toc
logical                      :: is_group_member

is_group_member = .false.

call check_groupindex(groupid, 'is_group_member', group_toc)

if (group_toc%groups(groupid)%member_count < 0) then
   call error_handler(E_ERR, 'is_group_member', &
                     'no members in requested group', &
                     source, revision, revdate)
endif

do i=1, group_toc%groups(groupid)%member_count
   if (group_toc%groups(groupid)%group_members(i) == itemnum) then
      is_group_member = .true.
      return
   endif
enddo

end function is_group_member

!----------------------------------------------------------------------

function get_group_member_by_index(groupid, membernum, group_toc)

integer,          intent(in) :: groupid
integer,          intent(in) :: membernum
type(type_group), intent(in) :: group_toc
integer                      :: get_group_member_by_index

call check_groupindex(groupid, 'get_group_member_by_index', group_toc)

if (group_toc%groups(groupid)%member_count <= 0) then
   call error_handler(E_ERR, 'get_group_member_by_index', &
                     'no members in requested group', &
                     source, revision, revdate)
endif

if (membernum > group_toc%groups(groupid)%member_count) then
   call error_handler(E_ERR, 'get_group_member_by_index', &
                     'membernum larger than number of members in this group', &
                     source, revision, revdate)
endif

get_group_member_by_index = group_toc%groups(groupid)%group_members(membernum)

end function get_group_member_by_index

!----------------------------------------------------------------------

subroutine del_group_member(groupid, itemnum, group_toc)

integer,          intent(in)    :: groupid
integer,          intent(in)    :: itemnum
type(type_group), intent(inout) :: group_toc

integer :: i, membernum
integer :: num_mem

call check_groupindex(groupid, 'del_group_member', group_toc)

membernum = -1
FINDME: do i=1, group_toc%groups(groupid)%member_count
   if (group_toc%groups(groupid)%group_members(i) == itemnum) then
      membernum = i
      exit FINDME
   endif
enddo FINDME

if (membernum < 0 .or. membernum > group_toc%groups(groupid)%member_count) then
   call error_handler(E_ERR, 'del_group_member', &
                      'membernum not a valid member number', &
                      source, revision, revdate)
endif

num_mem = group_toc%groups(groupid)%member_count - 1
   
do i=membernum, num_mem
   group_toc%groups(groupid)%group_members(i) = group_toc%groups(groupid)%group_members(i+1) 
enddo
   
group_toc%groups(groupid)%member_count = num_mem

end subroutine del_group_member

!----------------------------------------------------------------------

function get_group_count(groupid, group_toc)

integer,          intent(in) :: groupid
type(type_group), intent(in) :: group_toc
integer                      :: get_group_count

call check_groupindex(groupid, 'get_group_count', group_toc)

if (group_toc%groups(groupid)%member_count < 0) then
   call error_handler(E_ERR, 'get_group_count', &
                     'no members in requested group', &
                     source, revision, revdate)
endif

get_group_count = group_toc%groups(groupid)%member_count

end function get_group_count

!----------------------------------------------------------------------

subroutine print_group_contents(groupid, group_toc)

integer,          intent(in) :: groupid
type(type_group), intent(in) :: group_toc

integer :: i, num

call check_groupindex(groupid, 'print_group_contents', group_toc)

num = group_toc%groups(groupid)%member_count

write(*,*) 'group id number ', groupid
write(*,*) 'member count ', num

do i=1, num
   write(*, '(A14,I5,I5)') 'member, id ', i, group_toc%groups(groupid)%group_members(i)
enddo

end subroutine print_group_contents

!----------------------------------------------------------------------

subroutine check_groupindex(groupid, caller, group_toc)

integer,          intent(in) :: groupid
character(len=*), intent(in) :: caller
type(type_group), intent(in) :: group_toc

if (groupid <= 0) then
   call error_handler(E_ERR, caller, &
                     'internal error: group index must be 1 or larger', &
                     source, revision, revdate)
endif

if (groupid > group_toc%group_count) then
   call error_handler(E_ERR, caller, &
                     'internal error: group index larger than number of groups', &
                     source, revision, revdate)
endif

end subroutine check_groupindex

!----------------------------------------------------------------------

subroutine find_all_types_for_kind(givenkind, typecount, typesfound, typeidlist)
 
integer, intent(in)  :: givenkind
integer, intent(in)  :: typecount
integer, intent(out) :: typesfound
integer, intent(out) :: typeidlist(:)

integer :: i, kindindex

! given a generic kind, find all the specific types which
! map to this kind.  i believe this has to be done by iterating
! through every type to see what the kind is.

! this routine is dealing with the raw indices into the kinds
! and types array - it no longer is using the table of contents values.

typesfound = 0
typeidlist(:) = -1

do i=1, typecount
   kindindex = get_obs_kind_var_type(i)
   if (kindindex == givenkind) then
      typesfound = typesfound + 1
      typeidlist(typesfound) = i
   endif
enddo

end subroutine find_all_types_for_kind

!----------------------------------------------------------------------

subroutine update_impact(table, wordcount, itemlist, wordarray, toc, group_toc)

real(r8),         intent(inout) :: table(:,0:)
integer,          intent(in)    :: wordcount
integer,          intent(in)    :: itemlist(:)
character(len=*), intent(in)    :: wordarray(:)
type(type_toc),   intent(in)    :: toc
type(type_group), intent(in)    :: group_toc

real(r8) :: rvalue
integer :: itemval1, itemval2, groupid1, groupid2
integer :: groupitem1, groupitem2
integer :: itemcount1, itemcount2
logical :: isgroup1, isgroup2
integer :: i, j

! expect exactly 3 tokens on these lines:
!  kind_or_type_or_group   kind_or_group   real(r8)
! (where 2nd group cannot contain types)

call check_impact_line(wordcount, itemlist, wordarray, rvalue)

! ok, we have the first, second, and value.  update the table.
! any combination of kind names or group names is allowed.

call iteminfo(itemlist(1), itemval1, isgroup1, itemcount1, toc, group_toc, kindonly=.false.)
call iteminfo(itemlist(2), itemval2, isgroup2, itemcount2, toc, group_toc, kindonly=.true.)

! both single items
if (.not. isgroup1 .and. .not. isgroup2) then

   call fill_impact_table_values(table, itemval1, itemval2, rvalue, toc, itemlist(1), 'k/k')

! single item to group
else if (.not. isgroup1 .and. isgroup2) then
   groupid2 = itemval2
   do i=1, itemcount2
      groupitem2 = get_group_member_by_index(groupid2, i, group_toc)
      call require_only_kinds_in_group(groupitem2, groupid2, toc)
      itemval2 = get_toc_value(groupitem2, toc)   
      call fill_impact_table_values(table, itemval1, itemval2, rvalue, toc, itemlist(1), 'k/g')
   enddo

! group to single item
else if (isgroup1 .and. .not. isgroup2) then
   groupid1 = itemval1
   itemval2 = get_toc_value(itemlist(2), toc)
   do i=1, itemcount1
      groupitem1 = get_group_member_by_index(groupid1, i, group_toc)
      itemval1 = get_toc_value(groupitem1, toc)
      call fill_impact_table_values(table, itemval1, itemval2, rvalue, toc, groupitem1, 'g/k')
   enddo

! both groups
else if (isgroup1 .and. isgroup2) then
   groupid1 = itemval1
   groupid2 = itemval2
   do i=1, itemcount1
      groupitem1 = get_group_member_by_index(groupid1, i, group_toc)
      itemval1 = get_toc_value(groupitem1, toc)
      do j=1, itemcount2
         groupitem2 = get_group_member_by_index(groupid2, j, group_toc)
         call require_only_kinds_in_group(groupitem2, groupid2, toc)
         itemval2 = get_toc_value(groupitem2, toc)
         call fill_impact_table_values(table, itemval1, itemval2, rvalue, toc, groupitem1, 'g/g')
      enddo
   enddo

else
   call error_handler(E_ERR, 'update_impact', &
                     'first and second entries on line must be either a KIND or a group name', &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif

end subroutine update_impact

!----------------------------------------------------------------------
 
subroutine check_impact_line(wordcount, itemlist, wordarray, rvalue)

integer,          intent(in)  :: wordcount
integer,          intent(in)  :: itemlist(:)
character(len=*), intent(in)  :: wordarray(:)
real(r8),         intent(out) :: rvalue

if (wordcount /= 3) then
   call error_handler(E_ERR, 'check_impact_line', &
                     'impact lines require exactly 3 words per line', &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif
if (itemlist(1) < 0) then
   call error_handler(E_ERR, 'check_impact_line', &
                     'first word on line is unrecognized KIND, TYPE, or group name', &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif
if (itemlist(2) < 0) then
   call error_handler(E_ERR, 'check_impact_line', &
                     'second word on line is unrecognized KIND or group name', &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif

call extract_value(wordarray(3), rvalue)

! FIXME: do this at runtime now.  if you want to check values
! at table build time, uncomment this code again.
!if (.not. allow_any_factor_value) then
!   if (rvalue < 0.0_r8 .or. rvalue > 1.0_r8) then 
!      call error_handler(E_ERR, 'check_impact_line', &
!                        'impact values must be between 0 and 1, inclusive', &
!                        source, revision, revdate, text2=readbuf, &
!                        text3='set "allow_any_factor_value=.true." in namelist to allow')
!   endif
!endif

end subroutine check_impact_line

!----------------------------------------------------------------------

! TOOL:
subroutine fill_impact_table_values(table, tableindex1, tableindex2, rvalue, toc, tocindex1, label)

real(r8),         intent(inout)        :: table(:,0:)
integer,          intent(in)           :: tableindex1
integer,          intent(in)           :: tableindex2
real(r8),         intent(in)           :: rvalue
type(type_toc),   intent(in)           :: toc
integer,          intent(in)           :: tocindex1
character(len=*), intent(in), optional :: label

integer :: i, typeid(toc%type_count), typesfound
integer :: kindid

if (present(label)) then
!   print *, trim(label)
endif

! expand tableindex1 if type here and loop over list
if (get_toc_type(tocindex1, toc) == ENTRY_DARTKIND) then

   ! this routine is now dealing with the absolute indices into the types
   ! and kinds array; not the indices in the table of contents.  the values
   ! returned can be used directly to set the table entries.

   kindid = tableindex1 - 1   ! kinds must be offset down by one, 0:num_kinds
   call find_all_types_for_kind(kindid, toc%type_count, typesfound, typeid)
 
   do i=1, typesfound
      call set_impact_table_values(table, typeid(i), tableindex2, rvalue, toc, label)
   enddo

else

   call set_impact_table_values(table, tableindex1, tableindex2, rvalue, toc, label)

endif

end subroutine fill_impact_table_values

!----------------------------------------------------------------------

! TOOL:
subroutine set_impact_table_values(table, tableindex1, tableindex2, rvalue, toc, label)

real(r8),         intent(inout)        :: table(:,0:)
integer,          intent(in)           :: tableindex1
integer,          intent(in)           :: tableindex2
real(r8),         intent(in)           :: rvalue
type(type_toc),   intent(in)           :: toc
character(len=*), intent(in), optional :: label


if (present(label)) then
!   print *, trim(label)
endif

!print *, 'set_impact_table_values: ', tableindex1, tableindex2
if ((table(tableindex1, tableindex2) /= missing_r8) .and. &
    (table(tableindex1, tableindex2) /= rvalue)) then
   write(msgstring, *) 'error; different impact factor already set for this pair'
   write(msgstring2, *) 'previous value was', table(tableindex1, tableindex2), ' new value is ', rvalue
   write(msgstring3, *) trim(get_name_by_value(tableindex1, ENTRY_DARTTYPE, toc)), ' and ',  &
                        trim(get_name_by_value(tableindex2, ENTRY_DARTKIND, toc))
   call error_handler(E_ERR, 'update_impact', msgstring, &
                     source, revision, revdate, text2=msgstring2, text3=msgstring3)
   
endif

table(tableindex1, tableindex2) = rvalue

end subroutine set_impact_table_values

!----------------------------------------------------------------------

subroutine iteminfo(item, itemid, isgroup, itemcount, toc, group_toc, kindonly)

integer,          intent(in)  :: item
integer,          intent(out) :: itemid
logical,          intent(out) :: isgroup
integer,          intent(out) :: itemcount
type(type_toc),   intent(in)  :: toc
type(type_group), intent(in)  :: group_toc
logical,          intent(in)  :: kindonly

integer :: entryinfo

entryinfo = get_toc_type(item, toc) 

! our usage of this is asymmetric, but it is not clear how
! much the tool should know about this.
if(entryinfo == ENTRY_GROUP) then
   itemid = get_toc_value(item, toc)
   isgroup = .true.
   itemcount = get_group_count(itemid, group_toc)
   ! FIXME: same for group members for item 2 - not symmetric? check here?

else if (kindonly .and. entryinfo == ENTRY_DARTTYPE) then
   call error_handler(E_ERR, 'iteminfo', &
                     'entry on right, '//trim(get_toc_name(item,toc))//' can only be a DART KIND', &
                     source, revision, revdate, text2=errline, text3=readbuf)

else if (kindonly .and. entryinfo == ENTRY_DARTKIND) then
   itemid = get_toc_value(item, toc)
   isgroup = .false.
   itemcount = 1

else if (.not. kindonly .and. (entryinfo == ENTRY_DARTKIND .or. entryinfo == ENTRY_DARTTYPE)) then
   itemid = get_toc_value(item, toc)
   isgroup = .false.
   itemcount = 1

else

   call error_handler(E_ERR, 'iteminfo', &
                     'unrecognized input; expecting a DART KIND, TYPE, or GROUP name', &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif


end subroutine iteminfo

!----------------------------------------------------------------------

subroutine require_only_kinds_in_group(item, groupid, toc)

integer,        intent(in) :: item
integer,        intent(in) :: groupid
type(type_toc), intent(in) :: toc

character(len=string_length) :: gname

if (get_toc_type(item, toc) /= ENTRY_DARTKIND) then

   gname = get_name_by_value(groupid, ENTRY_GROUP, toc)
   call error_handler(E_ERR, 'require_only_kinds_in_group', &
                     'if in righthand column, group '//trim(gname)//' must include only DART KINDS', &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif

end subroutine require_only_kinds_in_group

!----------------------------------------------------------------------

subroutine must_be_specific_type(item, toc)

integer,        intent(in) :: item
type(type_toc), intent(in) :: toc

character(len=string_length) :: ename

if (get_toc_type(item, toc) /= ENTRY_DARTTYPE) then

   ename = get_toc_name(item, toc)
   call error_handler(E_ERR, 'must_be_specific_type', &
                     'lefthand column must be a DART TYPE, not '//trim(ename), &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif

end subroutine must_be_specific_type

!----------------------------------------------------------------------

subroutine must_be_generic_kind(item, toc)

integer,        intent(in) :: item
type(type_toc), intent(in) :: toc

character(len=string_length) :: ename

if (get_toc_type(item, toc) /= ENTRY_DARTKIND) then

   ename = get_toc_name(item, toc)
   call error_handler(E_ERR, 'must_be_generic_kind', &
                     'righthand column must be a DART KIND, not '//trim(ename), &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif

end subroutine must_be_generic_kind

!----------------------------------------------------------------------

function get_toc_type(item, toc)
 
integer,        intent(in) :: item
type(type_toc), intent(in) :: toc
integer                    :: get_toc_type

if (item > toc%toc_count) then
   call error_handler(E_ERR, 'get_toc_type', &
                     'called with an index larger than the number of items in toc', &
                     source, revision, revdate)
endif

!print *, 'get_toc_type: item = ', item
get_toc_type = toc%toc_entries(item)%etype

end function get_toc_type

!----------------------------------------------------------------------

function get_toc_value(item, toc)
 
integer,        intent(in) :: item
type(type_toc), intent(in) :: toc
integer                    :: get_toc_value

if (item > toc%toc_count) then
   call error_handler(E_ERR, 'get_toc_value', &
                     'called with an index larger than the number of items in toc', &
                     source, revision, revdate)
endif

get_toc_value = toc%toc_entries(item)%evalue

end function get_toc_value

!----------------------------------------------------------------------

function get_toc_name(item, toc)
 
integer,        intent(in)   :: item
type(type_toc), intent(in)   :: toc
character(len=string_length) :: get_toc_name

if (item > toc%toc_count) then
   call error_handler(E_ERR, 'get_toc_name', &
                     'called with an index larger than the number of items in toc', &
                     source, revision, revdate)
endif

get_toc_name = toc%toc_entries(item)%ename

end function get_toc_name

!----------------------------------------------------------------------

subroutine extract_value(wordarray, rvalue)

character(len=*), intent(in)  :: wordarray
real(r8),         intent(out) :: rvalue


! how to validate this?  look for chars other than 0-9 . - e d ?
read(wordarray, *, iostat=ios) rvalue
if (ios /= 0) then
   call error_handler(E_ERR, 'extract_value', &
                     'bad numeric value, expecting third item on line to be a number', &
                     source, revision, revdate, text2=errline, text3=readbuf)
endif

end subroutine extract_value

!----------------------------------------------------------------------

subroutine output_impact_table(oname, table, toc)
 
character(len=*), intent(in) :: oname
real(r8),         intent(in) :: table(:,0:)
type(type_toc),   intent(in) :: toc

integer :: i, j, funit
integer :: ntypes, nkinds

funit = open_file(oname, 'formatted', action='write')

ntypes = ubound(table, 1)
nkinds = ubound(table, 2)

do j=0, nkinds
   do i=1, ntypes
      if (table(i, j) /= 1.0_r8) then
print *, i, j, table(i, j)
         write(funit, '(A34,A34,F18.6)') &
                        trim(get_name_by_value(i, ENTRY_DARTTYPE, toc)), &
                        trim(get_name_by_value(j, ENTRY_DARTKIND, toc)), &
                        table(i, j)
      endif
   enddo
enddo

call close_file(funit)

end subroutine output_impact_table

!----------------------------------------------------------------------

end module obs_impact_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

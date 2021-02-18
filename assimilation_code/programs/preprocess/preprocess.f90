! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

!> Define the observation types and obs/state quantities understood by DART.
!>
!> Takes a list of observation type module path names and physical quantity
!> module path names and creates two output files that are compiled into
!> the rest of the DART executables.  They define the quantities and observation
!> types that are recognized and processed in any way. 
!>
!> The output module for observation types can contain multiple fragments of 
!> standard F90 that may be required to implement forward observation operators 
!> for DART. The sections are retrieved from the files by this program and 
!> inserted into the appropriate blanks in the DEFAULT_obs_def_mod.F90 and 
!> DEFAULT_obs_kind_mod.F90 templates. 
!>
!> The final obs_def_mod.f90 and obs_kind_mod.f90 that are created contain
!> the default code plus all the code required from the selected observation
!> type modules. Preprocess also inserts the required identifiers and strings
!> for the corresponding observation quantities.
!>
!> @todo FIXME: the module name is set by whatever is in the DEFAULT_xxx file
!> which is still obs_kind_mod instead of _qty_.  also i have changed the default
!> filenames here to be _qty_(HK?) but in reality all the pathnames files for the
!> build system expect _kind_ and the names are always specified in our input.nml
!> namelist files, so we aren't depending on the defaults. (maybe we should?)
!> and finally, this code generates module 'use' statements, and those, 
!> to keep things working, are still outputing "use obs_kind_mod" instead 
!> of "use obs_qty_mod".   this should change but i'm not sure when.
!>

program preprocess

use types_mod,      only : r8, MISSING_R8     ! @todo FIXME r8 is needed; is MISSING_R8?
use utilities_mod,  only : error_handler, E_ERR, E_MSG,   &
                           file_exist, open_file, close_file,              &
                           initialize_utilities, do_nml_file, do_nml_term, &
                           find_namelist_in_file, check_namelist_read,     &
                           finalize_utilities, log_it
use parse_args_mod, only : get_args_from_string, get_name_val_pairs_from_string, get_next_arg

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source = 'preprocess.f90'

! Pick something ridiculously large and forget about it (lazy)
integer, parameter   :: MAX_TYPES = 5000, MAX_QTYS = 5000

! max valid number of tokens per line.  allocate more
! to handle error conditions.
integer, parameter   :: MAX_TOKENS = 20
integer              :: ntokens
character(len=256)   :: token(MAX_TOKENS)
character(len=256)   :: valtokens(MAX_TOKENS)


character(len=256)   :: line, test, t_string
integer              :: iunit, io, i, j, k, l
integer              :: linenum1, linenum2, linenum3, linenum4
integer              :: num_types_found, num_qtys_found
logical              :: duplicate, qty_found, temp_user, is_more, match
character(len=512)   :: err_string, err_string2, err_string3
character(len=6)     :: full_line_in  = '(A256)'
character(len=3)     :: full_line_out = '(A)'

logical :: DEBUG = .false.

integer, parameter :: MAX_NAME_LEN = 32
character(len=MAX_NAME_LEN) :: temp_type, temp_qty

type qty_info_type
   character(len=MAX_NAME_LEN) :: name                 = 'null'
   integer                     :: num_nameval_pairs    = 0
   character(len=256)          :: namepair(MAX_TOKENS) = ''
   character(len=256)          :: valpair(MAX_TOKENS)  = ''
end type qty_info_type

type(qty_info_type) :: qty_info(0:MAX_QTYS)


type obs_info_type
   character(len=MAX_NAME_LEN) :: name         = 'null'
   character(len=MAX_NAME_LEN) :: qty          = 'null'
   logical                     :: has_usercode = .false.

end type obs_info_type

type(obs_info_type) :: obs_info(1:MAX_TYPES)

! specific marker strings
!                                                                1         2         3         4         5         6
!                                                       123456789012345678901234567890123456789012345678901234567890
character(len=*),parameter :: TYP_STARTDEFN_STRING  = '! BEGIN DART PREPROCESS TYPE DEFINITIONS'
character(len=*),parameter :: KND_STARTLIST_STRING  = '! BEGIN DART PREPROCESS KIND LIST'
character(len=*),parameter :: TYP_ENDDEFN_STRING    = '! END DART PREPROCESS TYPE DEFINITIONS'
character(len=*),parameter :: KND_ENDLIST_STRING    = '! END DART PREPROCESS KIND LIST'
character(len=*),parameter :: QTY_STARTDEFN_STRING  = '! BEGIN DART PREPROCESS QUANTITY DEFINITIONS'
character(len=*),parameter :: KND_STARTDEFN_STRING  = '! BEGIN DART PREPROCESS KIND DEFINITIONS'
character(len=*),parameter :: QTY_ENDDEFN_STRING    = '! END DART PREPROCESS QUANTITY DEFINITIONS'
character(len=*),parameter :: KND_ENDDEFN_STRING    = '! END DART PREPROCESS KIND DEFINITIONS'
character(len=*),parameter :: INSERT_INTS_STRING    = '! DART PREPROCESS INTEGER DECLARATIONS INSERTED HERE'
character(len=*),parameter :: INSERT_INIT_STRING    = '! DART PREPROCESS DERIVED TYPE INITIALIZATIONS INSERTED HERE'

! output format decorations
character(len=78) :: separator_line = &
'!---------------------------------------------------------------------------'
character(len=78) :: blank_line = &
'                                                                            '
!! currently unused, but available if wanted:
!character(len=12) :: start_line = '!  Start of '
!character(len=12) :: end_line =   '!  End of   '
!character(len=78) :: blank_comment_line = &
!'!                                                                           '

integer, parameter :: NUM_SECTIONS = 7
! List of the DART PREPROCESS strings for obs_def type files.
character(len=29) :: preprocess_string(NUM_SECTIONS) = (/ &
      'MODULE CODE                  ', &
      'USE FOR OBS_QTY_MOD          ', &
      'USE OF SPECIAL OBS_DEF MODULE', &
      'GET_EXPECTED_OBS_FROM_DEF    ', &
      'READ_OBS_DEF                 ', &
      'WRITE_OBS_DEF                ', &
      'INTERACTIVE_OBS_DEF          '/)

!! Must match the list above.  Is there a default code section that we can
!! autogenerate for this section?
!logical :: default_code(NUM_SECTIONS) = &
!   (/ .false., &
!      .false., &
!      .false., &
!      .true.,  &
!      .true.,  &
!      .true.,  &
!      .true.  /)

integer, parameter :: MODULE_ITEM = 1
integer, parameter :: QTY_ITEM = 2
integer, parameter :: USE_ITEM = 3
integer, parameter :: GET_EXPECTED_ITEM = 4
integer, parameter :: READ_ITEM = 5
integer, parameter :: WRITE_ITEM = 6
integer, parameter :: INTERACTIVE_ITEM = 7

integer :: num_obs_type_files = 0
integer :: num_quantity_files = 0
integer :: obs_def_in_unit, obs_def_out_unit
integer :: obs_qty_in_unit, obs_qty_out_unit, in_unit

integer, parameter   :: MAX_OBS_TYPE_FILES = 1000
integer, parameter   :: MAX_QUANTITY_FILES = 1000
logical :: file_has_usercode(MAX_OBS_TYPE_FILES) = .false.


integer, parameter :: NML_STRLEN = 256

! The namelist reads in a sequence of path_names that are absolute or
! relative to the working directory in which preprocess is being executed
! and these files are used to fill in observation type/qty details in
! DEFAULT_obs_def_mod.f90 and DEFAULT_obs_qty_mod.f90.

! original values:
character(len=NML_STRLEN) :: input_obs_def_mod_file = &
                        '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90'
character(len=NML_STRLEN) :: output_obs_def_mod_file = &
                        '../../../observations/forward_operators/obs_def_mod.f90'
character(len=NML_STRLEN) :: input_obs_kind_mod_file = &
                        '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90'
character(len=NML_STRLEN) :: output_obs_kind_mod_file = &
                        '../../../assimilation_code/modules/observations/obs_kind_mod.f90'
character(len=NML_STRLEN) :: input_files(MAX_OBS_TYPE_FILES) = 'null'
logical                   :: overwrite_output = .true.

! jump through hoops to maintain backwards compatibility in namelist.
! {input,output}_obs_def_mod_file, overwrite_output same as before. 
! but if these have been changed from the defaults, they override the older values.
character(len=NML_STRLEN) :: input_obs_qty_mod_file = 'was input_obs_kind_mod_file'
character(len=NML_STRLEN) :: output_obs_qty_mod_file = 'was output_obs_kind_mod_file'
character(len=NML_STRLEN) :: obs_type_files(MAX_OBS_TYPE_FILES) = 'was input_files'
character(len=NML_STRLEN) :: quantity_files(MAX_QUANTITY_FILES) = '_new_nml_item_'

namelist /preprocess_nml/ input_obs_def_mod_file, output_obs_def_mod_file,   &
                          input_obs_qty_mod_file, output_obs_qty_mod_file,   &
                          input_obs_kind_mod_file, output_obs_kind_mod_file, &  ! deprecate
                          obs_type_files, quantity_files, overwrite_output,  &
                          input_files                                           ! deprecate

!---------------------------------------------------------------------------
! start of program code

!Begin by reading the namelist
call initialize_utilities('preprocess')

! Read the namelist entry
call find_namelist_in_file("input.nml", "preprocess_nml", iunit)
read(iunit, nml = preprocess_nml, iostat = io)
call check_namelist_read(iunit, io, "preprocess_nml")

! mess with the namelists so we can rename items to be more accurate
! without breaking the world.  this should be removed when it is deemed
! ok to change namelist entries.
call ensure_backwards_compatibility()

! Output the namelist file information
call log_it('Path names of default obs_def and obs_qty modules')
call log_it(input_obs_def_mod_file)
call log_it(input_obs_qty_mod_file)

call log_it('Path names of output obs_def and obs_qty modules')
call log_it(output_obs_def_mod_file)
call log_it(output_obs_qty_mod_file)

! The name for the default files is required. Error if these are null.
call cannot_be_null(input_obs_def_mod_file, 'input_obs_def_mod_file')
call cannot_be_null(input_obs_qty_mod_file, 'input_obs_qty_mod_file')

call cannot_be_null(output_obs_def_mod_file, 'output_obs_def_mod_file')
call cannot_be_null(output_obs_qty_mod_file, 'output_obs_qty_mod_file')

call log_it('INPUT obs_def files:')

do i = 1, MAX_OBS_TYPE_FILES
   if(obs_type_files(i) == 'null') exit
   call log_it(obs_type_files(i))
   num_obs_type_files = i
enddo

call log_it('INPUT quantity files:')

do i = 1, MAX_QUANTITY_FILES
   if(quantity_files(i) == 'null') exit
   call log_it(quantity_files(i))
   num_quantity_files = i
enddo

! Open/Create the DEFAULT and OUTPUT files and set the unit numbers

call open_file_for_read(input_obs_def_mod_file, 'input file ', obs_def_in_unit)
call open_file_for_read(input_obs_qty_mod_file, 'input file ', obs_qty_in_unit)

call open_file_for_write(output_obs_def_mod_file, 'output file ', overwrite_output, obs_def_out_unit)
call open_file_for_write(output_obs_qty_mod_file, 'output file ', overwrite_output, obs_qty_out_unit)

!______________________________________________________________________________
! Preprocessing for the obs_qty module
! Read in the quantity table and build index numbers for them. 

! Copy over lines from obs_qty_template file up to the next insertion point
linenum1 = 0
call copy_until(obs_qty_in_unit,   input_obs_qty_mod_file, INSERT_INTS_STRING, linenum1, &
                obs_qty_out_unit, output_obs_qty_mod_file, .false.)

num_qtys_found = 0

! QTY_STATE_VARIABLE is always defined (and always offset 0).
qty_info(num_qtys_found)%name = 'QTY_STATE_VARIABLE'
qty_info(num_qtys_found)%namepair(1) = 'desc'
qty_info(num_qtys_found)%valpair(1)  = 'basic item in a state'
qty_info(num_qtys_found)%num_nameval_pairs = 1
num_qtys_found = num_qtys_found + 1

!> for each quantity input file, read in quantities to be used.
!> duplicates are ignored (so you can list overlapping subsets 
!> of quantities).  units and (if specified) bounds must match
!> if duplicate.

SEARCH_QUANTITY_FILES: do j = 1, num_quantity_files

   call open_file_for_read(quantity_files(j), 'quantity_files', in_unit)

   ! Read until the ! BEGIN QUANTITY DEFINITION marker string is found
   linenum2 = 0
   call read_until(in_unit, quantity_files(j), QTY_STARTDEFN_STRING, linenum2, KND_STARTDEFN_STRING)

   ! Subsequent lines can contain QTY_xxx lines or comments or
   ! the end string.
   DEFINE_QTYS: do
      call get_next_line(in_unit, full_line_in, QTY_ENDDEFN_STRING, &
                         quantity_files(j), linenum2, KND_ENDDEFN_STRING)

      ! Look for the ! END QTY LIST in the current line
      test = adjustl(line)
      if(test == QTY_ENDDEFN_STRING) exit DEFINE_QTYS
      if(test == KND_ENDDEFN_STRING) then
         ! @todo FIXME: after some adjustment time, here is 
         ! where you would print out a 'deprecated' warning
         ! if the alternate delimiter form is encountered.
         exit DEFINE_QTYS
      endif
   
      ! All lines between start/end must be comments or QTY_xxx strings
      !
      ! Formats accepted:  
      !
      ! QTY_string
      ! QTY_string name=value ...
      !
      ! QTY_string ! comments
      !       
      ! ! comment
      !
   
      ! ntokens < 0 means error
      ! ntokens == 0 means cycle without error
      ! ntokens > 0 means some work to do

      call parse_line(line, ntokens, token, err_string, pairs_expected=.true., valtokens=valtokens)
      if (ntokens < 0) &
         call quantity_error(err_string, line, quantity_files(j), linenum2)
   
      if (ntokens == 0) cycle DEFINE_QTYS

      if (token(1)(1:4) /= 'QTY_') then
         err_string = 'QTY_xxx not found as first word on line'
         call quantity_error(err_string, line, quantity_files(j), linenum2)
      endif

      if (len_trim(token(1)) > MAX_NAME_LEN) then
         write(err_string, *) 'Quantity names are limited to ', MAX_NAME_LEN, ' characters'
         call quantity_error(err_string, line, quantity_files(j), linenum2)
      endif

      ! this loop adds a new quantity to the string array if it's new.
      duplicate = .false.
      qtys: do i=0, num_qtys_found-1
         if (qty_info(i)%name == token(1)) then
            duplicate = .true.
            call resolve_duplicates(i, ntokens, token, valtokens, quantity_files(j), linenum2) 
            exit qtys
         endif
      enddo qtys
   
      if (.not. duplicate) then
         qty_info(num_qtys_found)%name = token(1)
         do i=1, ntokens-1
            qty_info(num_qtys_found)%namepair(i) = token(i+1)
            qty_info(num_qtys_found)%valpair(i)  = valtokens(i+1)
         enddo
         qty_info(num_qtys_found)%num_nameval_pairs = ntokens-1
         num_qtys_found = num_qtys_found + 1
      endif
   
   enddo DEFINE_QTYS

   ! Close this quantity file
   call close_file(in_unit)

enddo SEARCH_QUANTITY_FILES

!______________________________________________________________________________
! Preprocessing of obs_def files for info to add to the obs_qty module
! Get all the type/qty strings from all of the obs_def files 
! up front and then insert stuff.  Easier to error check and combine
! duplicate quantities.

! Initial number of types is 0, qtys is based on previous section
! (unlike quantities, types are 1 based)
num_types_found = 0

SEARCH_OBS_DEF_FILES: do j = 1, num_obs_type_files

   call open_file_for_read(obs_type_files(j), 'obs_type_files', in_unit)

   ! Read until the ! BEGIN TYPE DEFINITIONS is found
   linenum2 = 0
   call read_until(in_unit, obs_type_files(j), TYP_STARTDEFN_STRING, linenum2, KND_STARTLIST_STRING)

   ! Subsequent lines contain the type_identifier (same as type_string), and
   ! qty_string separated by commas, and optional usercode flag
   EXTRACT_TYPES: do
      call get_next_line(in_unit, full_line_in, TYP_ENDDEFN_STRING, &
                         obs_type_files(j), linenum2, KND_ENDLIST_STRING)

      ! Look for the ! END TYPE DEFINITIONS in the current line
      test = adjustl(line)
      if(test == TYP_ENDDEFN_STRING) exit EXTRACT_TYPES
      if(test == KND_ENDLIST_STRING) then
         ! FIXME: here is where you could print out a 'deprecated' warning
         ! if the alternate delimiter form is encountered.
         exit EXTRACT_TYPES
      endif

      ! All lines between start/end must be type/qty lines.
      ! Format:  ! type_string, qty_string [, COMMON_CODE]

      call parse_line(line, ntokens, token, err_string, pairs_expected=.false.)
      if (ntokens < 0) &
         call typeqty_error(err_string, line, obs_type_files(j), linenum2)
   
      if (ntokens == 0) cycle EXTRACT_TYPES

      if (ntokens < 2 .or. ntokens > 3) then
         err_string = 'expected OBS_TYPE  QTY_xxx  [COMMON_CODE]. unable to process.'
         call typeqty_error(err_string, line, obs_type_files(j), linenum2)
      endif

      if (token(2)(1:4) /= 'QTY_') then
         err_string = 'QTY_xxx not found as second word on line'
         call typeqty_error(err_string, line, obs_type_files(j), linenum2)
      endif

      ! save results in temp vars for now, so we can check for
      ! duplicates (not allowed in types) or duplicates (which are
      ! expected in quantities)
      if (len_trim(token(1)) > len(temp_type)) then
         write(err_string, *) 'Observation type names are limited to ', MAX_NAME_LEN, ' characters'
         call typeqty_error(err_string, line, obs_type_files(j), linenum2)
      endif
      temp_type = trim(token(1))

      if (len_trim(token(2)) > len(temp_qty)) then
         write(err_string, *) 'Quantity names are limited to ', MAX_NAME_LEN, ' characters'
         call typeqty_error(err_string, line, obs_type_files(j), linenum2)
      endif
      temp_qty = trim(token(2))

      if (ntokens == 3) then
         if (token(3) /= 'COMMON_CODE') then
            err_string = 'if third word present on line, it must be COMMON_CODE'
            call typeqty_error(err_string, line, obs_type_files(j), linenum2)
         endif
         temp_user = .false.
      else
         temp_user = .true.
         file_has_usercode(j) = .true.
      endif 

      ! Another type/qty line; increment the type count.  
      ! Check the qtys for unrecognized quantities.

      ! only allow duplicate definitions if the qtys match and there is
      ! no usercode in either the original definition or this one
      do i=1, num_types_found
         if (obs_info(i)%name == temp_type) then
            ! allow dups if 1) same qty 2) no usercode
            if ((obs_info(i)%qty /= temp_qty) .or. &
                (obs_info(i)%has_usercode) .or. temp_user) then
               err_string = &
                  'Duplicate! This observation type has already been processed with different options'
               call typeqty_error(err_string, line, obs_type_files(j), linenum2)
            else ! dup, can safely ignore?
               cycle EXTRACT_TYPES
            endif
         endif
      enddo

      num_types_found = num_types_found + 1
      obs_info(num_types_found)%name = trim(temp_type)
      obs_info(num_types_found)%qty = trim(temp_qty) 
      obs_info(num_types_found)%has_usercode = temp_user

      ! multiple obs types mapping to the same qty is ok.
      qty_found = .false.
      do i=0, num_qtys_found - 1
         if (qty_info(i)%name == temp_qty) then
            qty_found = .true.
            exit
         endif
      enddo
      if (.not. qty_found) then
         ! an error to specify a quantity not already defined
         err_string = 'Unrecognized QTY_xxx quantity'
         call typeqty_error(err_string, line, obs_type_files(j), linenum2)
      endif

   enddo EXTRACT_TYPES

   ! Close this obs_def file
   call close_file(in_unit)
enddo SEARCH_OBS_DEF_FILES


! Write out the integer declaration lines for qtys

call write_separator_line(obs_qty_out_unit)
call write_blank_line(obs_qty_out_unit)
write(obs_qty_out_unit, '(A)') '! Integer definitions for DART QUANTITIES'
call write_blank_line(obs_qty_out_unit)
do i = 0, num_qtys_found-1
   write(obs_qty_out_unit, '(A30,A32,A3,I5)') &
      'integer, parameter, public :: ', trim(qty_info(i)%name), ' = ', i
enddo
call write_blank_line(obs_qty_out_unit)

!--

call write_blank_line(obs_qty_out_unit)
write(line, '(A,I5)') 'integer, parameter, public :: MAX_DEFINED_QUANTITIES = ', &
   num_qtys_found-1
write(obs_qty_out_unit, '(A)') trim(line)
call write_blank_line(obs_qty_out_unit)
call write_separator_line(obs_qty_out_unit)

!--

! Write out the integer declaration lines for types
call write_separator_line(obs_qty_out_unit)
call write_blank_line(obs_qty_out_unit)
write(obs_qty_out_unit, '(A)') '! Integer definitions for DART OBS TYPES'
call write_blank_line(obs_qty_out_unit)
do i = 1, num_types_found 
   write(obs_qty_out_unit, '(A30,A32,A3,I5)') &
      'integer, parameter, public :: ', trim(obs_info(i)%name), ' = ', i
enddo
call write_blank_line(obs_qty_out_unit)

!--

! Write out the max num of obs_types, too
call write_blank_line(obs_qty_out_unit)
write(line, '(A,I5)') 'integer, parameter, public :: MAX_DEFINED_TYPES_OF_OBS = ', &
   num_types_found
write(obs_qty_out_unit, '(A)') trim(line)
call write_blank_line(obs_qty_out_unit)
call write_separator_line(obs_qty_out_unit)

! Copy over lines up to the next insertion point
call copy_until(obs_qty_in_unit,   input_obs_qty_mod_file, INSERT_INIT_STRING, linenum1, &
                obs_qty_out_unit, output_obs_qty_mod_file, .false.)

!--

! Write out the definitions of each entry of obs_qty_info
do i = 0, num_qtys_found-1
   write(line, *) 'obs_qty_info(', i, ')%index = ', i
   write(obs_qty_out_unit, '(A)') trim(line)

   write(line, *) 'obs_qty_info(', i, ')%name = ', '"'//trim(qty_info(i)%name)//'"'
   write(obs_qty_out_unit, '(A)') trim(line)

   write(line, *) 'obs_qty_info(', i, ')%nitems = ', qty_info(i)%num_nameval_pairs
   write(obs_qty_out_unit, '(A)') trim(line)

   do j=1, qty_info(i)%num_nameval_pairs
      write(line, *) 'obs_qty_info(', i, ')%itemname(', j, ') = ', '"'//trim(qty_info(i)%namepair(j))//'"'
      write(obs_qty_out_unit, '(A)') trim(line)
      write(line, *) 'obs_qty_info(', i, ')%itemvalue(', j, ') = ', '"'//trim(qty_info(i)%valpair(j))//'"'
      write(obs_qty_out_unit, '(A)') trim(line)
   enddo
   call write_blank_line(obs_qty_out_unit)  ! too many blank lines?
enddo
call write_blank_line(obs_qty_out_unit)

! Write out the definitions of each entry of obs_type_info
do i = 1, num_types_found
   write(line, '(A,I5,3A)') 'obs_type_info(', i, ') = obs_type_type(', &
      trim(obs_info(i)%name), ", &"
   write(obs_qty_out_unit, '(A)') trim(line)
   write(line, *) '   ', "'", trim(obs_info(i)%name), "', ", &
      trim(obs_info(i)%qty), ', .false., .false., .false.)'
   write(obs_qty_out_unit, '(A)') trim(line)
enddo
call write_blank_line(obs_qty_out_unit)


! Copy over rest of lines
call copy_until_end(obs_qty_in_unit, obs_qty_out_unit)

call close_file(obs_qty_in_unit)
call close_file(obs_qty_out_unit)

!______________________________________________________________________________

!______________________________________________________________________________
! Now do the obs_def files
! Read DEFAULT file line by line and copy into output file until
! Each insertion point is found. At the insertion points, copy the
! appropriate code from each requested obs_def file into the output obs_def
! file and then proceed.

! There are seven special code sections (ITEMS) in the obs_def file at present.
! Each copies code in from the special type specific obs_qty modules
linenum4 = 0
ITEMS: do i = 1, NUM_SECTIONS

   t_string = '! DART PREPROCESS ' // trim(preprocess_string(i)) // ' INSERTED HERE'
   call copy_until(obs_def_in_unit, input_obs_def_mod_file, t_string, linenum4, &
                   obs_def_out_unit, output_obs_def_mod_file, .false.)

   ! The 'USE FOR OBS_QTY_MOD' section is handled differently; lines are not
   ! copied, they are generated based on the list of types and qtys.
   if(i == QTY_ITEM) then
      ! Create use statements for both the QTY_ quantities and the individual
      ! observation type strings.
      call write_separator_line(obs_def_out_unit)
      call write_blank_line(obs_def_out_unit)
      do k = 1, num_types_found
         write(obs_def_out_unit, '(A)') &
            'use obs_kind_mod, only : ' // trim(obs_info(k)%name)     !!!!! FIXME: kind -> qty
      enddo

      call write_blank_line(obs_def_out_unit)
      do k = 0, num_qtys_found-1
         write(obs_def_out_unit, '(A)') &
            'use obs_kind_mod, only : ' // trim(qty_info(k)%name)      !!!!! FIXME: kind -> qty
      enddo
      call write_blank_line(obs_def_out_unit)

      call write_separator_line(obs_def_out_unit)
      call write_blank_line(obs_def_out_unit)
      cycle ITEMS
   endif


   ! Insert the code for this ITEM from each requested obs_def 'modules'
   do j = 1, num_obs_type_files
      if (.not. file_has_usercode(j)) then
         if (i == MODULE_ITEM) then
            call write_separator_line(obs_def_out_unit)
            write(obs_def_out_unit, '(2A)') &
               '!No module code needed for ', trim(obs_type_files(j))
            call write_separator_line(obs_def_out_unit)
         endif
         !if (i == USE_ITEM) then
         !   call write_separator_line(obs_def_out_unit)
         !   write(obs_def_out_unit, '(2A)') &
         !      '!No use statements needed for ', trim(obs_type_files(j))
         !   call write_separator_line(obs_def_out_unit)
         !endif
         cycle
      endif

      ! Since there might someday be a lot of these, 
      ! open and close them each time needed
      call open_file_for_read(obs_type_files(j), 'obs_type_files', in_unit)
      linenum3 = 0

      ! Read until the appropriate ITEM # label is found in the input 
      ! for this obs_type
      t_string = '! BEGIN DART PREPROCESS ' // trim(preprocess_string(i)) 
      call read_until(in_unit, obs_type_files(j), t_string, linenum3)
      
      ! decoration or visual separation, depending on your viewpoint
      if (i == MODULE_ITEM) then
         call write_separator_line(obs_def_out_unit)
         write(obs_def_out_unit, '(2A)') '! Start of code inserted from ', &
            trim(obs_type_files(j))
         call write_separator_line(obs_def_out_unit)
         call write_blank_line(obs_def_out_unit)
      endif

      ! Copy all code until the end of item into the output obs_def file
 
      t_string = '! END DART PREPROCESS ' // trim(preprocess_string(i))
      call copy_until(in_unit, obs_type_files(j), t_string, linenum3, &
                      obs_def_out_unit, output_obs_def_mod_file, (i /= MODULE_ITEM))
         

      ! decoration or visual separation, depending on your viewpoint
      if (i == MODULE_ITEM) then
         call write_blank_line(obs_def_out_unit)
         call write_separator_line(obs_def_out_unit)
         write(obs_def_out_unit, '(2A)') '! End of code inserted from ', &
            trim(obs_type_files(j))
         call write_separator_line(obs_def_out_unit)
      endif

      ! Got everything from this file, move along
      call close_file(in_unit)
  enddo
   
  ! Now check to see if this item has any types which are expecting us
  ! to automatically generate the code for them.
  do j = 1, num_types_found
     if (obs_info(j)%has_usercode) cycle

     select case (i)
     case (GET_EXPECTED_ITEM)
        write(obs_def_out_unit, '(A)')  '      case(' // trim(obs_info(j)%name) // ')'
        write(obs_def_out_unit, '(3A)') '         call interpolate(state_handle, ens_size, location, ', &
           trim(obs_info(j)%qty), ', expected_obs, istatus)'
     case (READ_ITEM, WRITE_ITEM, INTERACTIVE_ITEM)
        write(obs_def_out_unit, '(A)')  '   case(' // trim(obs_info(j)%name) // ')'
        write(obs_def_out_unit, '(A)')  '      continue'
     case default
       ! nothing to do for others 
     end select
  enddo

enddo ITEMS

call copy_until_end(obs_def_in_unit, obs_def_out_unit)

call close_file(obs_def_in_unit)
call close_file(obs_def_out_unit)

call error_handler(E_MSG,'preprocess','Finished successfully.')
call finalize_utilities('preprocess')

!------------------------------------------------------------------------------

contains

!> read from the given unit number until you run out of input
!> (which is a fatal error) or until the requested line contents
!> is found.  update the linenum so reasonable error messages
!> can be provided.

subroutine read_until(iunit, iname, stop_string, linenum, alt_stop_string)

integer,          intent(in) :: iunit
character(len=*), intent(in) :: iname
character(len=*), intent(in) :: stop_string
integer,          intent(inout) :: linenum
character(len=*), intent(in), optional :: alt_stop_string

call do_until(iunit, iname, stop_string, linenum, .false., 0, '', .false., alt_stop_string)

end subroutine read_until

!------------------------------------------------------------------------------

subroutine copy_until(iunit, iname, stop_string, linenum, ounit, oname, trimfirst, alt_stop_string)

integer,          intent(in) :: iunit
character(len=*), intent(in) :: iname
character(len=*), intent(in) :: stop_string
integer,          intent(inout) :: linenum
integer,          intent(in) :: ounit
character(len=*), intent(in) :: oname
logical,          intent(in) :: trimfirst
character(len=*), intent(in), optional :: alt_stop_string

call do_until(iunit, iname, stop_string, linenum, .true., ounit, oname, trimfirst, alt_stop_string)

end subroutine copy_until

!------------------------------------------------------------------------------

subroutine do_until(iunit, iname, stop_string, linenum, docopy, ounit, oname, trimfirst, alt_stop_string)

integer,          intent(in) :: iunit
character(len=*), intent(in) :: iname
character(len=*), intent(in) :: stop_string
integer,          intent(inout) :: linenum
logical,          intent(in) :: docopy
integer,          intent(in) :: ounit
character(len=*), intent(in) :: oname
logical,          intent(in) :: trimfirst
character(len=*), intent(in), optional :: alt_stop_string

integer :: ierr
character(len=256) :: line   ! this should match the format length

! Read until the given string is found
FIND_NEXT: do

   read(iunit, full_line_in, IOSTAT = ierr) line

   ! If end of file, input file is incomplete or weird stuff happened
   if(ierr /= 0) then
      write(err_string,  *) 'Did not find required line containing ', trim(stop_string)
      write(err_string2, *) 'reading file ', trim(iname)
      call error_handler(E_ERR, 'do_until', err_string, source, text2=err_string2)
   endif
   linenum = linenum + 1

   ! Look for the required string in the current line
   if(index(line, trim(stop_string)) > 0) exit FIND_NEXT
   if(present(alt_stop_string)) then
      if(index(line, trim(alt_stop_string)) > 0) exit FIND_NEXT
      ! FIXME: here is where you could print out a 'deprecated' warning
      ! if the alternate delimiter form is encountered.
   endif

   ! if doing a copy, write it verbatim to the output file
   if (docopy) then
       if (trimfirst) then
          write(ounit, full_line_out, IOSTAT = ierr) trim(line(2:))
       else
          write(ounit, full_line_out, IOSTAT = ierr) trim(line)
       endif
       if(ierr /= 0) then
          write(err_string,  *) 'Write error, returned code = ', ierr
          write(err_string2, *) 'writing file ', trim(oname)
          call error_handler(E_ERR, 'do_until', err_string, source, text2=err_string2)
       endif
   endif


enddo FIND_NEXT

end subroutine do_until

!------------------------------------------------------------------------------

subroutine copy_until_end(iunit, ounit)

integer, intent(in) :: iunit
integer, intent(in) :: ounit

integer :: ierr

DO_NEXT: do

   read(iunit, full_line_in, IOSTAT = ierr) line
   if (ierr /= 0) exit DO_NEXT

   write(ounit, full_line_out) trim(line)

enddo DO_NEXT

end subroutine copy_until_end

!------------------------------------------------------------------------------

subroutine typeqty_error(errtext, line, filename, linenum)

character(len=*), intent(in) :: errtext, line, filename
integer,          intent(in) :: linenum

write(err_string2,'(A,I5,A)') '"'//trim(filename)//'", line ', linenum, &
                              ', contents: "'//trim(line)//'"'
write(err_string3,   '(A,A)') 'expected:  ! UniqueObsType, Quantity  ',&
                              ' or  ! UniqueObsType, Quantity, COMMON_CODE'
call error_handler(E_ERR, 'typeqty_error', errtext, source, &
                   text2=err_string2, text3=err_string3)

end subroutine typeqty_error

!------------------------------------------------------------------------------

subroutine quantity_error(errtext, line, filename, linenum)

character(len=*), intent(in) :: errtext, line, filename
integer,          intent(in) :: linenum

write(err_string2, '(A,I5,A)')'"'//trim(filename)//'", has bad line (', linenum,')'
write(err_string3,      '(A)')'bad line contents: ['//trim(line)//']'
call error_handler(E_ERR, 'quantity_error', errtext, source, &
                   text2=err_string2, text3=err_string3)

end subroutine quantity_error

!-------------------------------------------------------------------------------

! get the next line from the input file, and bump up the line count
! it is an error if you get to the end of the file.

subroutine get_next_line(iunit, format, end_string, filename, linenum, alt_end_string)

integer,          intent(in) :: iunit
character(len=*), intent(in) :: format, end_string, filename
integer,          intent(inout) :: linenum
character(len=*), intent(in), optional :: alt_end_string

integer :: ierr

read(iunit, format, IOSTAT = ierr) line
if (ierr == 0) then
   linenum = linenum + 1
   return
endif

! If end of file, input file is incomplete or weird stuff happened
write(err_string, *) 'file "', trim(filename), &
                  '" does NOT contain ', trim(end_string), ' or ', trim(alt_end_string)
call error_handler(E_ERR, 'get_next_line', err_string, source)

end subroutine get_next_line

!------------------------------------------------------------------------------

! if the file exists, open for reading.
! otherwise it is an error.

subroutine open_file_for_read(filename, label, fileunit)
 character(len=*), intent(in) :: filename, label
 integer, intent(out) :: fileunit

if(file_exist(filename)) then
   fileunit = open_file(filename, action='read')
   return
endif

! If file does not exist it is an error
write(err_string, *) trim(label) // ' "' // trim(filename), & 
                     '" does NOT exist (and must).'
call error_handler(E_ERR, 'open_file_for_read', err_string, source)

end subroutine open_file_for_read

!------------------------------------------------------------------------------

! if the file does NOT exist, create.  if override flag is set,
! allow an existing file to be completely overwritten.
! otherwise it is an error.

subroutine open_file_for_write(filename, label, overwrite, fileunit)

character(len=*), intent(in)  :: filename, label
logical,          intent(in)  :: overwrite
integer,          intent(out) :: fileunit

if((.not. file_exist(filename)) .or. overwrite) then
   fileunit = open_file(filename, action='write')
   return
endif

! If file *does* exist and we haven't said ok to overwrite, error
write(err_string, *) trim(label) // ' "' // trim(filename), &
                     '" exists and will not be overwritten: Please remove or rename'
call error_handler(E_ERR, 'open_file_for_write', err_string, source)

end subroutine open_file_for_write

!------------------------------------------------------------------------------

subroutine cannot_be_null(varvalue, varname)

character(len=*), intent(in) :: varvalue, varname

if(varvalue /= 'null') return

call error_handler(E_ERR, 'cannot_be_null', &
                  'Namelist must provide ' // trim(varname),source)

end subroutine cannot_be_null

!------------------------------------------------------------------------------

subroutine write_separator_line(unitnum)

integer, intent(in) :: unitnum

write(unitnum, '(A)') separator_line

end subroutine write_separator_line

!------------------------------------------------------------------------------

subroutine write_blank_line(unitnum)

integer, intent(in) :: unitnum

write(unitnum, '(A)') blank_line

end subroutine write_blank_line

!------------------------------------------------------------------------------

! can this be common for qtys & types?  maybe...
! ntokens < 0 means error
! ntokens == 0 means cycle without error
! ntokens > 0 means some work to do

! apparently intent(out) strings must have a length.

subroutine parse_line(line, ntokens, tokens, estring, pairs_expected, valtokens)

character(len=*),   intent(in)  :: line
integer,            intent(out) :: ntokens
character(len=256), intent(out) :: tokens(:)
character(len=512), intent(out) :: estring
logical,            intent(in)  :: pairs_expected
character(len=256), intent(out), optional :: valtokens(:)

character(len=256) :: test
integer :: i
integer :: trimmed_len, start_of_comment
integer :: npairs, endoff
character(len=256) :: namepairs(MAX_TOKENS), valpairs(MAX_TOKENS)

ntokens = 0
estring = ''

! plan:
! 1) strip the required leading ! in col 1
! 2) truncate string at next !
! 3) parse into tokens by spaces, respecting " or ' strings (not yet)
! 4) strip trailing commas off tokens for now.

if (line(1:1) /= '!') then 
   ntokens = -1
   estring = 'line must begin with !'
   return
endif

test = adjustl(line(2:))
trimmed_len = len_trim(test)

start_of_comment = index(test, '!')
if (start_of_comment > 0) &
  test(start_of_comment:trimmed_len) = ' '

trimmed_len = len_trim(test)
if (trimmed_len == 0) return

! hack for now - replace all commas from line with spaces.
! this could mess things up if there is a unit type with a 
! comma in it, or a description which includes a comma.
! FIXME: reconsider this decision

do 
   i = index(test, ',')
   if (i == 0) exit
   test(i:i) = ' '
enddo

! get rid of tabs.
do 
   i = index(test, achar(9))
   if (i == 0) exit
   test(i:i) = ' '
enddo

! parse here
if (pairs_expected .and. present(valtokens)) then
   call get_next_arg(test, 1, tokens(1), endoff)
   valtokens(1) = ''
   call get_name_val_pairs_from_string(test(endoff:), npairs, namepairs, valpairs, is_more)
   ntokens = 1 + npairs
   do i=2, ntokens
      tokens(i) = namepairs(i-1)
      valtokens(i) = valpairs(i-1)
   enddo
else
   call get_args_from_string(test, ntokens, tokens)
endif

! get anything?
if (ntokens <= 0) return

if (DEBUG) then
   print *, "line: ", trim(test)
   if (pairs_expected) then
      print *, "ntokens, name/val tokens: ", ntokens
      print *, 1, trim(tokens(1))
      do i=2, ntokens
         print *, i, " ", trim(tokens(i)), ' = ', trim(valtokens(i))
      enddo
      print *, "is more?", is_more
   else
      print *, "ntokens, tokens: ", ntokens
      do i=1, ntokens
         print *, i, " ", trim(tokens(i))
      enddo
   endif
endif

end subroutine parse_line

!------------------------------------------------------------------------------

! this routine must resolve what happens if another entry for the same
! quantity is encountered.  its choices are: to fatally error out, to add
! a metadata name=val pair to the existing entry, print a warning,
! allow a duplicate if the new entry is a subset of the old -- see the
! comments in the code for what each section is doing in case someone
! wants to change the behavior, add printed warnings, etc.

subroutine resolve_duplicates(qty_indx, ntokens, tname, tval, infile, linenum)

integer,          intent(in) :: qty_indx
integer,          intent(in) :: ntokens
character(len=*), intent(in) :: tname(MAX_TOKENS)
character(len=*), intent(in) :: tval(MAX_TOKENS)
character(len=*), intent(in) :: infile
integer,          intent(in) :: linenum
 
integer :: k, l, first_t, last_t
logical :: match(MAX_TOKENS)

! two checks needed here.
!
! one is to find all existing entries in the new list and make
! sure the values match.  also track which ones we've found, so we
! can flag missing entries.
!
! two is looking for new entries which aren't in the existing list.
! make it easy to see where to change the behavior if someone
! decides that should be ok.  for now, it is an error.
! 
! in both cases, metadata items could be in a different order
! which is not an error.

if (DEBUG) print *, 'working on quantity ', trim(qty_info(qty_indx)%name)

! this loop starts at 2 because 1 is the QTY name.
! the token pairs are tname(2), tval(2), tname(3), tval(3), etc.
tokens: do l=2, ntokens
   match(l) = .false.
   ! loop over any existing name,val pairs already associated with this entry
   do k=1, qty_info(qty_indx)%num_nameval_pairs 

      if (qty_info(qty_indx)%namepair(k) /= tname(l)) cycle 

      ! l=new token index, k = existing token index
      if (qty_info(qty_indx)%valpair(k) /= tval(l)) then
         if (DEBUG) print *, 'found inconsistent value for same name: ', &
                    trim(qty_info(qty_indx)%valpair(k))," /= ",trim(tval(l))
         call incompatible_duplicates(qty_indx, ntokens, tname, tval, infile, linenum)
      endif
      if (DEBUG) print *, 'found match to existing name/value pair: ', trim(tname(l))," == ",trim(tval(l))
      match(l) = .true.
   enddo
   if (match(l)) cycle tokens

   ! we have found a new token.  for now, we aren't allowing this.
   ! if you want to add it, comment this line out and comment in the code below
   if (DEBUG) print *, 'found new name not in previous entries, "'//trim(tval(l))//'"'
   call incompatible_duplicates(qty_indx, ntokens, tname, tval, infile, linenum)

   ! here is the code if you wanted to combine new metadata items from
   ! different occurrances of the same QTY line.

   !k = qty_info(qty_indx)%num_nameval_pairs + 1
   !qty_info(qty_indx)%num_nameval_pairs = k
   !qty_info(qty_indx)%namepair(k) = tname(l)
   !qty_info(qty_indx)%valpair(k) = tval(l)
enddo tokens

! if we get here and all tokens haven't been matched, error out.
! comment this section out if you want to allow this instead.

! the valid match array values start at 2 since token 1 is the
! quantity, and go to number of pairs + 1 (since skipping first
! token).  this line is making sure all existing metadata pairs
! have also been specified here in a duplicate entry.
first_t = 2
l = qty_info(qty_indx)%num_nameval_pairs
last_t = (2+l) - 1
if ((l > 0) .and. (.not. all(match(first_t:last_t)))) then
   if (DEBUG) print *, 'not all metadata matched existing entry:'
   if (DEBUG) print *, l, match(first_t:last_t)

   call incompatible_duplicates(qty_indx, ntokens, tname, tval, infile, linenum)
endif

end subroutine resolve_duplicates

!------------------------------------------------------------------------------

! we have decided these don't match.  give the user help figuring out why.
! this routine currently errors out with a fatal error, so does not return.

subroutine incompatible_duplicates(qty_indx, ntokens, tname, tval, infile, linenum)

integer,          intent(in) :: qty_indx
integer,          intent(in) :: ntokens
character(len=*), intent(in) :: tname(MAX_TOKENS)
character(len=*), intent(in) :: tval(MAX_TOKENS)
character(len=*), intent(in) :: infile
integer,          intent(in) :: linenum

integer :: i, j
character(len=*), parameter :: routine = 'incompatible_duplicates'
character(len=1) :: empty = ''

! get existing entry from qty_info(qty_indx)
! get new info from ntokens, tname=tval
! print infile, linenum for newer duplicate - we don't have
! the original entry info anymore

call error_handler(E_MSG, empty, empty)
call error_handler(E_MSG, empty, 'Incompatible duplicate entry detected')

if (qty_info(qty_indx)%num_nameval_pairs == 0) then
   write(err_string, "(A)") 'Existing entry for '//trim(qty_info(qty_indx)%name)//' has no metadata'
   call error_handler(E_MSG, empty, err_string, source)
else
   write(err_string, "(A)") 'Existing entry for '//trim(qty_info(qty_indx)%name)//' has this metadata:'
   call error_handler(E_MSG, empty, err_string)
   do i=1, qty_info(qty_indx)%num_nameval_pairs 
      write(err_string, "(A,I3,A)") 'item ', i, ":  '"//trim(qty_info(qty_indx)%namepair(i))//"="//trim(qty_info(qty_indx)%valpair(i))//"'"
      call error_handler(E_MSG, empty, err_string, source)
   enddo
endif

call error_handler(E_MSG, empty, empty)

if (ntokens <= 1) then
   write(err_string, "(A)") 'Duplicate entry for '//trim(token(1))// ' has no metadata'
   call error_handler(E_MSG, empty, err_string, source)
else
   write(err_string, "(A)") 'Duplicate entry for '//trim(token(1))// ' has this metadata:'
   call error_handler(E_MSG, empty, err_string, source)
   do i=2, ntokens
      ! token 1 is the qty name.  
      write(err_string, "(A,I3,A)") 'item ', i-1, ":  '"//trim(tname(i))//"="//trim(tval(i))//"'"
      call error_handler(E_MSG, empty, err_string, source)
   enddo
endif
write(err_string, "(A,I5)") 'File '//trim(infile)//' at line number ', linenum
call error_handler(E_MSG, empty, err_string, source)


call error_handler(E_MSG, empty, empty)
call error_handler(E_ERR, routine, 'Resolve entries to proceed', source)

end subroutine incompatible_duplicates

!------------------------------------------------------------------------------

! mess with the namelist entries to remain backwards compatible.
! after some adjustment time, start warning about deprecated names
! and eventually error out and remove this subroutine.

subroutine ensure_backwards_compatibility()

integer :: i

! copy over the older named entries to the new names if the new name
! was not found/set by reading the namelist.

if (input_obs_qty_mod_file == 'was input_obs_kind_mod_file') &
   input_obs_qty_mod_file = input_obs_kind_mod_file

if (output_obs_qty_mod_file == 'was output_obs_kind_mod_file') &
   output_obs_qty_mod_file = output_obs_kind_mod_file

! this is trickier because it is an array of strings.  if you
! initialize it, all the entries are set to the initial string.
! the code expects 'null' to end the list.  simple solution is
! to test and copy over all input_file entries which should set
! unused entries to null.
if (obs_type_files(1) == 'was input_files') then
   do i=1, MAX_OBS_TYPE_FILES
     obs_type_files(i) = input_files(i)
   enddo
endif

! since we changed the defaults, fix up if caller did specify 
! one or more qty files.
do i=1, MAX_OBS_TYPE_FILES
   if (obs_type_files(i) == 'was input_files') obs_type_files(i) = 'null'
enddo

! for backwards compatibility, if no quantity files are given use this default file
! which contains all previously defined quantities.  again, since it is a string array
! set all entries to null then set the first one to the single default file.
if (quantity_files(1) == '_new_nml_item_') then
   quantity_files(:) = 'null'
   quantity_files(1) = default_quantity_file()
endif

! since we changed the defaults, fix up if caller did specify 
! one or more qty files.
do i=1, MAX_QUANTITY_FILES
   if (quantity_files(i) == '_new_nml_item_') quantity_files(i) = 'null'
enddo

end subroutine ensure_backwards_compatibility

!------------------------------------------------------------------------------

!> Determine the location of the file providing backwards-compatible behavior
!> if people do not supply one or more specific xxx_quantities_mod.f90 files
!> in the preprocess_nml:quantity_files namelist.
!>
!> Remove this routine once backwards compatibility issues have been fully deprecated

function default_quantity_file()

character(len=256) :: default_quantity_file

integer :: i

default_quantity_file = &
   'assimilation_code/modules/observations/default_quantities_mod.f90'

ITERATE : do i = 1,10

   ! RETURN successfully if the filename exists
   if ( file_exist(default_quantity_file) ) return

   ! Check to see if the candidate name is too long
   if (len_trim(default_quantity_file) + 4 >= len(default_quantity_file)) exit ITERATE

   write(default_quantity_file,'(A)') '../'//trim(default_quantity_file)

enddo ITERATE

write(err_string ,*)'Unable to find relative location of default_quantities_mod.f90'
write(err_string2,*)'Normally located in assimilation_code/modules/observations in DART source tree'
write(err_string3,*)'Set namelist preprocess_nml:quantity_files to appropriate quantity file(s)'
call error_handler(E_ERR, 'default_quantity_file', err_string, source, text2=err_string2, text3=err_string3)

end function default_quantity_file

!------------------------------------------------------------------------------

end program preprocess


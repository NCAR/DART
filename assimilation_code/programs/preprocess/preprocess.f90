! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Takes a list of observation type module path names and physical quantity
!> module path names and creates two output files that are compiled into
!> the rest of the DART executables.  They define the quantities and observation
!> types that are recognized and processed in any way. 
!>
!> The output module for observation types can contain multiple fragments of 
!> standard F90 that may be required to implement forward observation operators 
!> for DART. The sections are retrieved from the files by this program and 
!> inserted into the appropriate blanks in the DEFAULT_obs_def_mod.F90 and 
!> DEFAULT_obs_qty_mod.F90 templates. 
!>
!> The final obs_def_mod.f90 and obs_qty_mod.f90 that are created contain
!> the default code plus all the code required from the selected observation
!> type modules. Preprocess also inserts the required identifiers and strings
!> for the corresponding observation quantities.
!>
!> @todo FIXME: the module name is set by whatever is in the DEFAULT_xxx file
!> which is still obs_kind_mod instead of _qty_.  also i have changed the default
!> filenames here to be _qty_ but in reality all the pathnames files for the
!> build system expect _kind_ and the names are always specified in our input.nml
!> namelist files, so we aren't depending on the defaults. (maybe we should?)
!> and finally, this code generates module 'use' statements, and those, 
!> to keep things working, are still outputing "use obs_kind_mod" instead 
!> of "use obs_qty_mod".   this should change but i'm not sure when.
!>
!> @todo FIXME: this version is being changed to require at least 2 strings
!> per line - the qty name and units (which could be "none" if there really
!> are none).  
!> it can also have 2 additional optional tokens on the line which are any 
!> physical limits for the min and max values of this quantity.  for example, 
!> temperature in K can't be negative. relative humidity must be [0,1].  
!>
!> don't specify them if there aren't any limits, and use MISSING_R8 (or NA?)
!> for any bounds where there is only either a minimum or maximum 
!> but not the other. (if you specify either, two numbers will be required.)
!>
!> note to me:
!> since this is a text processor, could NA, null, or MISSING_R8 all be used?
!> the token MISSING_R8 should be output but may not ever be interpreted by 
!> this program. is r8 really needed as well?  again, text based processing.  
!> i guess for a numeric value it should be vetted to be valid here?
!>
!> @todo FIXME: should we hardcode QTY_STATE_VARIABLE as 0 here?
!>

program preprocess

use types_mod,      only : r8, MISSING_R8     ! @todo FIXME r8 is needed; is MISSING_R8?
use utilities_mod,  only : register_module, error_handler, E_ERR, E_MSG,   &
                           file_exist, open_file, close_file,              &
                           initialize_utilities, do_nml_file, do_nml_term, &
                           find_namelist_in_file, check_namelist_read,     &
                           finalize_utilities, log_it
use parse_args_mod, only : get_args_from_string

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! Pick something ridiculously large and forget about it (lazy)
integer, parameter   :: max_types = 5000, max_qtys = 5000


character(len = 256) :: line, test, test2, t_string
integer              :: iunit, ierr, io, i, j, k
integer              :: l_string, l2_string, total_len
integer              :: linenum1, linenum2, linenum3, linenum4
integer              :: num_types_found, num_qtys_found
logical              :: duplicate, qty_found, temp_user
character(len = 512) :: err_string
character(len = 6)   :: full_line_in  = '(A256)'
character(len = 3)   :: full_line_out = '(A)'

logical :: DEBUG = .false.

! max valid number of tokens per line is 4.  allocate more
! to handle error conditions.
integer, parameter   :: MAX_TOKENS = 20
integer              :: ntokens
character(len=128)   :: token(MAX_TOKENS)

integer, parameter :: MAX_NAME_LEN = 32
character(len=MAX_NAME_LEN) :: temp_type, temp_qty

type qty_info_type
   character(len=MAX_NAME_LEN) :: name          = 'null'
   character(len=128)          :: ud_units_name = 'none'
   real(r8)                    :: minbound      = MISSING_R8
   real(r8)                    :: maxbound      = MISSING_R8
end type qty_info_type

type(qty_info_type) :: qty_info(0:max_qtys)


type obs_info_type
   character(len=MAX_NAME_LEN) :: name         = 'null'
   character(len=MAX_NAME_LEN) :: qty          = 'null'
   logical                     :: has_usercode = .false.

end type obs_info_type

type(obs_info_type) :: obs_info(1:max_types)

! specific marker strings
!                                                             1         2         3         4         5         6
!                                                    123456789012345678901234567890123456789012345678901234567890
character(len=*),parameter :: qty_start_string   = '! BEGIN DART PREPROCESS QUANTITY LIST'
character(len=*),parameter :: qty_end_string     = '! END DART PREPROCESS QUANTITY LIST'
character(len=*),parameter :: qty2_start_string  = '! BEGIN DART PREPROCESS QUANTITY DEFINITIONS'
character(len=*),parameter :: qty2_end_string    = '! END DART PREPROCESS QUANTITY DEFINITIONS'
character(len=*),parameter :: insert_ints_string = '! DART PREPROCESS INTEGER DECLARATIONS INSERTED HERE'
character(len=*),parameter :: insert_init_string = '! DART PREPROCESS DERIVED TYPE INITIALIZATIONS INSERTED HERE'

! output format decorations
character(len = 78) :: separator_line = &
'!---------------------------------------------------------------------------'
character(len = 78) :: blank_line = &
'                                                                            '
!! currently unused, but available if wanted:
!character(len = 12) :: start_line = '!  Start of '
!character(len = 12) :: end_line =   '!  End of   '
!character(len = 78) :: blank_comment_line = &
!'!                                                                           '

integer, parameter :: NUM_SECTIONS = 7
! List of the DART PREPROCESS strings for obs_def type files.
character(len = 29) :: preprocess_string(NUM_SECTIONS) = (/ &
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

integer, parameter :: module_item = 1
integer, parameter :: qty_item = 2
integer, parameter :: use_item = 3
integer, parameter :: get_expected_item = 4
integer, parameter :: read_item = 5
integer, parameter :: write_item = 6
integer, parameter :: interactive_item = 7

integer :: num_obs_type_files = 0
integer :: num_quantity_files = 0
integer :: obs_def_in_unit, obs_def_out_unit
integer :: obs_qty_in_unit, obs_qty_out_unit, in_unit

integer, parameter   :: max_obs_type_files = 1000
integer, parameter   :: max_quantity_files = 1000
logical :: file_has_usercode(max_obs_type_files) = .false.

! The namelist reads in a sequence of path_names that are absolute or
! relative to the working directory in which preprocess is being executed
! and these files are used to fill in observation type/qty details in
! DEFAULT_obs_def_mod.f90 and DEFAULT_obs_qty_mod.f90.
character(len = 256) :: input_obs_def_mod_file = &
                        '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90'
character(len = 256) :: output_obs_def_mod_file = &
                        '../../../observations/forward_operators/obs_def_mod.f90'
character(len = 256) :: input_obs_qty_mod_file = &
                        '../../../assimilation_code/modules/observations/DEFAULT_obs_qty_mod.F90'
!HK this file does not exist
character(len = 256) :: output_obs_qty_mod_file = &
                        '../../../assimilation_code/modules/observations/obs_qty_mod.f90'
character(len = 256) :: obs_type_files(max_obs_type_files) = 'null'
character(len = 256) :: quantity_files(max_quantity_files) = 'null'
logical              :: overwrite_output = .true.

namelist /preprocess_nml/ input_obs_def_mod_file, input_obs_qty_mod_file,   &
                          output_obs_def_mod_file, output_obs_qty_mod_file, &
                          obs_type_files, quantity_files, overwrite_output

!---------------------------------------------------------------------------
! start of program code

!Begin by reading the namelist
call initialize_utilities('preprocess')
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "preprocess_nml", iunit)
read(iunit, nml = preprocess_nml, iostat = io)
call check_namelist_read(iunit, io, "preprocess_nml")

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

do i = 1, max_obs_type_files
   if(obs_type_files(i) == 'null') exit
   call log_it(obs_type_files(i))
   num_obs_type_files= i
enddo

call log_it('INPUT quantity files:')

do i = 1, max_quantity_files
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
call copy_until(obs_qty_in_unit,   input_obs_qty_mod_file, insert_ints_string, linenum1, &
                obs_qty_out_unit, output_obs_qty_mod_file, .false.)

num_qtys_found = 0

! QTY_STATE_VARIABLE is always defined (and always offset 0).
qty_info(num_qtys_found)%name = 'QTY_STATE_VARIABLE'
num_qtys_found = num_qtys_found + 1

!> for each quantity input file, read in quantities to be used.
!> duplicates are ignored (so you can list overlapping subsets 
!> of quantities).  units and (if specified) bounds must match
!> if duplicate.

SEARCH_QUANTITY_FILES: do j = 1, num_quantity_files

   call open_file_for_read(quantity_files(j), 'quantity_files', in_unit)

   ! Read until the ! BEGIN QUANTITY LIST marker string is found
   linenum2 = 0
   call read_until(in_unit, quantity_files(j), qty2_start_string, linenum2)

   ! Subsequent lines can contain QTY_xxx lines or comments or
   ! the end string.
   DEFINE_QTYS: do
      call get_next_line(in_unit, full_line_in, qty2_end_string, &
                         quantity_files(j), linenum2)

      ! Look for the ! END QTY LIST in the current line
      test = adjustl(line)
      if(test == qty2_end_string) exit DEFINE_QTYS
   
      ! All lines between start/end must be comments or QTY_xxx strings
      !
      ! Formats accepted:  
      !
      ! QTY_string, "units"
      ! QTY_string, "units" ! optional comment
      !
      ! QTY_string, "units", minbound, maxbound
      ! QTY_string, "units", minbound, maxbound  ! both required if either specified
      !       
      ! ! comment
      ! ! if setting only one bound, use MISSING_R8 as the value?  or NA? 
      ! !  it is a string at this point so can be anything we choose.
      ! !
   
      ! ntokens < 0 means error
      ! ntokens == 0 means cycle without error
      ! ntokens > 0 means some work to do

      call parse_line(line, ntokens, token, err_string)
      if (ntokens < 0) &
         call quantity_error(err_string, line, quantity_files(j), linenum2)
   
      if (ntokens == 0) cycle DEFINE_QTYS

      if (ntokens > 4 .or. ntokens == 3) then
         err_string = 'expected QTY_xxx units minbound maxbound. unable to process.'
         call quantity_error(err_string, line, quantity_files(j), linenum2)
      endif

      if (ntokens < 2) then
         err_string = 'expected QTY_xxx units. unable to process.'
         call quantity_error(err_string, line, quantity_files(j), linenum2)
      endif

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
      do i=0, num_qtys_found-1
         if (qty_info(i)%name == token(1)) then
            duplicate = .true.
            ! @todo FIXME test dups here for identity in all other fields
            exit
         endif
      enddo
   
      if (.not. duplicate) then
         qty_info(num_qtys_found)%name = trim(token(1))
         qty_info(num_qtys_found)%ud_units_name = trim(token(2))
         if (ntokens > 2) then
            call string_to_real(token(3), qty_info(num_qtys_found)%minbound, err_string)
            if (err_string /= '') then
               err_string = trim(err_string) // ' min bounds'
               call quantity_error(err_string, line, quantity_files(j), linenum2)
            endif
            call string_to_real(token(4), qty_info(num_qtys_found)%maxbound, err_string)
            if (err_string /= '') then
               err_string = trim(err_string) // ' max bounds'
               call quantity_error(err_string, line, quantity_files(j), linenum2)
            endif
         endif
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

   ! Read until the ! BEGIN QUANTITY LIST is found
   linenum2 = 0
   call read_until(in_unit, obs_type_files(j), qty_start_string, linenum2)

   ! Subsequent lines contain the type_identifier (same as type_string), and
   ! qty_string separated by commas, and optional usercode flag
   EXTRACT_TYPES: do
      call get_next_line(in_unit, full_line_in, qty_end_string, &
                         obs_type_files(j), linenum2)

      ! Look for the ! END QUANTITY LIST in the current line
      test = adjustl(line)
      if(test == qty_end_string) exit EXTRACT_TYPES

      ! All lines between start/end must be type/qty lines.
      ! Format:  ! type_string, qty_string [, COMMON_CODE]

      call parse_line(line, ntokens, token, err_string)
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
write(line, '(A,I5)') 'integer, parameter, public :: max_defined_quantities = ', &
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
write(line, '(A,I5)') 'integer, parameter, public :: max_defined_types_of_obs = ', &
   num_types_found
write(obs_qty_out_unit, '(A)') trim(line)
call write_blank_line(obs_qty_out_unit)
call write_separator_line(obs_qty_out_unit)

! Copy over lines up to the next insertion point
call copy_until(obs_qty_in_unit,   input_obs_qty_mod_file, insert_init_string, linenum1, &
                obs_qty_out_unit, output_obs_qty_mod_file, .false.)

!--

! Write out the definitions of each entry of obs_qty_names
call write_blank_line(obs_qty_out_unit)
write(line, '(A)') 'do i = 0, max_defined_quantities'
write(obs_qty_out_unit, '(A)') trim(line)
write(line, '(A)') '   obs_qty_names(i) = obs_qty_type(i, "UNKNOWN", "none", MISSING_R8, MISSING_R8)'
write(obs_qty_out_unit, '(A)') trim(line)
write(line, '(A)') 'enddo'
write(obs_qty_out_unit, '(A)') trim(line)
call write_blank_line(obs_qty_out_unit)


! Write out the definitions of each entry of obs_qty_names
do i = 0, num_qtys_found-1
   write(line, '(A,I5,5A)') 'obs_qty_names(', &
      i, ') = obs_qty_type(', trim(qty_info(i)%name),  ", '", &
      trim(qty_info(i)%name), "', &"
   write(obs_qty_out_unit, '(A)') trim(line)
   write(line, '(3A,2(F20.12,A))') "     '", &
      trim(qty_info(i)%ud_units_name), &
      "' ,", &
      qty_info(i)%minbound, &
      ", ", &
      qty_info(i)%maxbound, &
      ")"
   write(obs_qty_out_unit, '(A)') trim(line)
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
   if(i == qty_item) then
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
         if (i == module_item) then
            call write_separator_line(obs_def_out_unit)
            write(obs_def_out_unit, '(2A)') &
               '!No module code needed for ', trim(obs_type_files(j))
            call write_separator_line(obs_def_out_unit)
         endif
         !if (i == use_item) then
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
      if (i == module_item) then
         call write_separator_line(obs_def_out_unit)
         write(obs_def_out_unit, '(2A)') '! Start of code inserted from ', &
            trim(obs_type_files(j))
         call write_separator_line(obs_def_out_unit)
         call write_blank_line(obs_def_out_unit)
      endif

      ! Copy all code until the end of item into the output obs_def file
 
      t_string = '! END DART PREPROCESS ' // trim(preprocess_string(i))
      call copy_until(in_unit, obs_type_files(j), t_string, linenum3, &
                      obs_def_out_unit, output_obs_def_mod_file, (i /= module_item))
         

      ! decoration or visual separation, depending on your viewpoint
      if (i == module_item) then
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
     case (get_expected_item)
        write(obs_def_out_unit, '(A)')  '      case(' // trim(obs_info(j)%name) // ')'
        write(obs_def_out_unit, '(3A)') '         call interpolate(state_handle, ens_size, location, ', &
           trim(obs_info(j)%qty), ', expected_obs, istatus)'
     case (read_item, write_item, interactive_item)
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

call error_handler(E_MSG,'preprocess','Finished successfully.',source,revision,revdate)
call finalize_utilities('preprocess')

!------------------------------------------------------------------------------

contains

!> read from the given unit number until you run out of input
!> (which is a fatal error) or until the requested line contents
!> is found.  update the linenum so reasonable error messages
!> can be provided.
subroutine read_until(iunit, iname, stop_string, linenum)

integer, intent(in) :: iunit
character(len=*), intent(in) :: iname
character(len=*), intent(in) :: stop_string
integer, intent(inout) :: linenum

call do_until(iunit, iname, stop_string, linenum, .false., 0, '', .false.)

end subroutine read_until

!------------------------------------------------------------------------------

subroutine copy_until(iunit, iname, stop_string, linenum, ounit, oname, trimfirst)

integer, intent(in) :: iunit
character(len=*), intent(in) :: iname
character(len=*), intent(in) :: stop_string
integer, intent(inout) :: linenum
integer, intent(in) :: ounit
character(len=*), intent(in) :: oname
logical, intent(in) :: trimfirst

call do_until(iunit, iname, stop_string, linenum, .true., ounit, oname, trimfirst)

end subroutine copy_until

!------------------------------------------------------------------------------

subroutine do_until(iunit, iname, stop_string, linenum, docopy, ounit, oname, trimfirst)

integer, intent(in) :: iunit
character(len=*), intent(in) :: iname
character(len=*), intent(in) :: stop_string
integer, intent(inout) :: linenum
logical, intent(in) :: docopy
integer, intent(in) :: ounit
character(len=*), intent(in) :: oname
logical, intent(in) :: trimfirst

integer :: ierr
character(len=256) :: line   ! this should match the format length
character(len=512) :: err_string2

! Read until the given string is found
FIND_NEXT: do

   read(iunit, full_line_in, IOSTAT = ierr) line

   ! If end of file, input file is incomplete or weird stuff happened
   if(ierr /= 0) then
      write(err_string,  *) 'Did not find required line containing ', trim(stop_string)
      write(err_string2, *) 'reading file ', trim(iname)
      call error_handler(E_ERR, 'preprocess', err_string, &
            source, revision, revdate, text2=err_string2)
   endif
   linenum = linenum + 1

   ! Look for the required string in the current line
   if(index(line, trim(stop_string)) > 0) exit FIND_NEXT

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
          call error_handler(E_ERR, 'preprocess', err_string, &
                source, revision, revdate, text2=err_string2)
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

subroutine typeqty_error(errtext, line, file, linenum)
 character(len=*), intent(in) :: errtext, line, file
 integer, intent(in) :: linenum

call error_handler(E_MSG, 'preprocess error:', &
   'obs_def file has bad Obs Type/Qty line')
call error_handler(E_MSG, 'preprocess error:', errtext)
call error_handler(E_MSG, 'expected input:', &
   '! UniqueObsType, Quantity   or  ! UniqueObsType, Quantity, COMMON_CODE')
write(err_string, '(2A,I5)') trim(file), ", line number", linenum
call error_handler(E_MSG, 'bad file:', err_string)
call error_handler(E_MSG, 'bad line contents:', line)
write(err_string, *) 'See msg lines above for error details'
call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)

end subroutine typeqty_error

!------------------------------------------------------------------------------

subroutine quantity_error(errtext, line, file, linenum)
 character(len=*), intent(in) :: errtext, line, file
 integer, intent(in) :: linenum

call error_handler(E_MSG, 'preprocess error:', &
   'obs_qty file has bad Quantity line')
call error_handler(E_MSG, 'preprocess error:', errtext)
call error_handler(E_MSG, 'expected input:', &
   '! QTY_xxx "units"   or  ! comment line ')
call error_handler(E_MSG, 'or:', &
   '! QTY_xxx "units" min_bounds max_bounds')
call error_handler(E_MSG, 'where:', 'if unbounded on one side, use NA or MISSING_R8')
write(err_string, '(2A,I5)') trim(file), ", line number", linenum
call error_handler(E_MSG, 'bad file:', err_string)
call error_handler(E_MSG, 'bad line contents:', line)
write(err_string, *) 'See msg lines above for error details'
call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)

end subroutine quantity_error

!------------------------------------------------------------------------------

! get the next line from the input file, and bump up the line count
! it is an error if you get to the end of the file.

subroutine get_next_line(iunit, format, end_string, filename, linenum)
 integer, intent(in) :: iunit
 character(len=*), intent(in) :: format, end_string, filename
 integer, intent(inout) :: linenum

integer :: ierr

read(iunit, format, IOSTAT = ierr) line
if (ierr == 0) then
   linenum = linenum + 1
   return
endif

! If end of file, input file is incomplete or weird stuff happened
write(err_string, *) 'file ', trim(filename), &
                  ' does NOT contain ', trim(end_string)
call error_handler(E_ERR, 'preprocess', err_string, &
                   source, revision, revdate)

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
write(err_string, *) trim(label) // ' ' // trim(filename), & 
                     ' does NOT exist (and must).'
call error_handler(E_ERR, 'preprocess', err_string, &
      source, revision, revdate)

end subroutine open_file_for_read

!------------------------------------------------------------------------------

! if the file does NOT exist, create.  if override flag is set,
! allow an existing file to be completely overwritten.
! otherwise it is an error.

subroutine open_file_for_write(filename, label, overwrite, fileunit)
 character(len=*), intent(in) :: filename, label
 logical, intent(in) :: overwrite
 integer, intent(out) :: fileunit

if((.not. file_exist(filename)) .or. overwrite) then
   fileunit = open_file(filename, action='write')
   return
endif

! If file *does* exist and we haven't said ok to overwrite, error
write(err_string, *) trim(label) // ' ' // trim(filename), &
                     ' exists and will not be overwritten: Please remove or rename'
call error_handler(E_ERR, 'preprocess', err_string, &
      source, revision, revdate)

end subroutine open_file_for_write

!------------------------------------------------------------------------------

subroutine cannot_be_null(varvalue, varname)
 character(len=*), intent(in) :: varvalue, varname

if(varvalue /= 'null') return

call error_handler(E_ERR, 'preprocess', &
                  'Namelist must provide ' // trim(varname), &
                   source, revision, revdate)

end subroutine cannot_be_null

!------------------------------------------------------------------------------

! can this be common for qtys & types?  maybe...
! ntokens < 0 means error
! ntokens == 0 means cycle without error
! ntokens > 0 means some work to do

subroutine parse_line(line, ntokens, tokens, estring)
 character(len=*), intent(in)  :: line
 integer,          intent(out) :: ntokens
 character(len=*), intent(out) :: tokens(MAX_TOKENS)
 character(len=*), intent(out) :: estring

character(len=256) :: test
integer :: i
integer :: trimmed_len, start_of_comment

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

! hack for now - remove all commas from line.
! this will only mess things up if there is a
! unit type with a comma in it.

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
call get_args_from_string(test, ntokens, tokens)
if (DEBUG) then
   print *, "line: ", trim(test)
   print *, "ntokens, tokens: ", ntokens
   do i=1, ntokens
      print *, i, trim(tokens(i))
   enddo
endif

end subroutine parse_line

!------------------------------------------------------------------------------

subroutine string_to_real(s, r, estring)
 character(len=*), intent(in) :: s
 real(r8), intent(out) :: r
 character(len=*), intent(out) :: estring

estring = ''
r = MISSING_R8

if (s == 'MISSING_R8' .or. s == 'NA') return

read(s, *, iostat=ierr) r
if(ierr /= 0) then
   estring = 'unable to convert string to real number'   
   return
endif

end subroutine string_to_real

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

end program preprocess


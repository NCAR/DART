! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program preprocess

! Takes a list of observation type module path names. These modules contain
! multiple fragments of standard F90 that may be required to implement forward
! observation operators for DART. The sections are retrieved from the files
! by this program and inserted into the appropriate blanks in the
! DEFAULT_obs_def_mod.F90 and DEFAULT_obs_qty_mod.F90 templates. 
! The final obs_def_mod.f90 and obs_qty_mod.f90 that are created contain
! the default code plus all the code required from the selected observation
! type modules. Preprocess also inserts the required identifier and string
! for the corresponding observation quantities (and only those quantities).

! NEED TO ADD IN ALL THE ERROR STUFF

use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG,   &
                          file_exist, open_file,                          &
                          initialize_utilities, do_nml_file, do_nml_term, &
                          find_namelist_in_file, check_namelist_read,     &
                          finalize_utilities, log_it

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! Pick something ridiculously large and forget about it (lazy)
integer, parameter   :: max_types = 5000, max_qtys = 5000
character(len = 256) :: line, test, test2, type_string(max_types), &
                        qty_string(0:max_qtys), t_string, temp_type, temp_qty
integer              :: iunit, ierr, io, i, j, k
integer              :: l_string, l2_string, total_len, linenum1, linenum2
integer              :: num_types_found, num_qtys_found, qty_index(0:max_types)
logical              :: duplicate, qty_found, usercode(max_types), temp_user
character(len = 512) :: err_string
character(len = 6)   :: full_line_in  = '(A256)'
character(len = 3)   :: full_line_out = '(A)'

! specific marker strings
!                                                             1         2         3         4         5
!                                                    12345678901234567890123456789012345678901234567890
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

! List of the DART PREPROCESS strings for obs_def type files.
character(len = 29) :: preprocess_string(8) = (/ &
      'MODULE CODE                  ', &
      'USE FOR OBS_QTY_MOD          ', &
      'USE OF SPECIAL OBS_DEF MODULE', &
      'GET_EXPECTED_OBS_FROM_DEF    ', &
      'READ_OBS_DEF                 ', &
      'WRITE_OBS_DEF                ', &
      'INTERACTIVE_OBS_DEF          ', &
      'THE EIGHTH ONE IS UNDEFINED  '/)

!! Must match the list above.  Is there a default code section that we can
!! autogenerate for this section?
!logical :: default_code(8) = &
!   (/ .false., &
!      .false., &
!      .false., &
!      .true.,  &
!      .true.,  &
!      .true.,  &
!      .true.,  &
!      .false. /)

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

! A path for the default files is required. Have an error if these are null.
if(input_obs_def_mod_file == 'null') &
   call error_handler(E_ERR, 'preprocess', &
      'Namelist must provide input_obs_def_mod_file', &
      source, revision, revdate)
if(input_obs_qty_mod_file == 'null') &
   call error_handler(E_ERR, 'preprocess', &
      'Namelist must provide input_obs_qty_mod_file', &
      source, revision, revdate)

! A path for the output files is required. Have an error if these are null.
if(output_obs_def_mod_file == 'null') &
   call error_handler(E_ERR, 'preprocess', &
      'Namelist must provide output_obs_def_mod_file', &
      source, revision, revdate)
if(output_obs_qty_mod_file == 'null') &
   call error_handler(E_ERR, 'preprocess', &
      'Namelist must provide output_obs_qty_mod_file', &
      source, revision, revdate)

call log_it('INPUT obs_def files follow:')

do i = 1, max_obs_type_files
   if(obs_type_files(i) == 'null') exit
   call log_it(obs_type_files(i))
   num_obs_type_files= i
end do

call log_it('INPUT quantity files follow:')

do i = 1, max_quantity_files
   if(quantity_files(i) == 'null') exit
   call log_it(quantity_files(i))
   num_quantity_files = i
end do

! Try to open the DEFAULT and OUTPUT files
! DEFAULT files must exist or else an error
if(file_exist(input_obs_def_mod_file)) then
   ! Open the file for reading
   obs_def_in_unit = open_file(input_obs_def_mod_file, action='read')
else
   ! If file does not exist it is an error
   write(err_string, *) 'file ', trim(input_obs_def_mod_file), &
      ' must exist (and does not)'
   call error_handler(E_ERR, 'preprocess', err_string, &
      source, revision, revdate)
endif

if(file_exist(input_obs_qty_mod_file)) then
   ! Open the file for reading
   obs_qty_in_unit = open_file(input_obs_qty_mod_file, action='read')
else
   ! If file does not exist it is an error
   write(err_string, *) 'file ', trim(input_obs_qty_mod_file), &
      ' must exist (and does not)'
   call error_handler(E_ERR, 'preprocess', err_string, &
      source, revision, revdate)
endif

! Error if Output files EXIST, unless 'overwrite_output' is TRUE
if(.not. file_exist(output_obs_def_mod_file) .or. overwrite_output) then
   ! Open (create) the file for writing
   obs_def_out_unit = open_file(output_obs_def_mod_file, action='write')
else
   ! If file *does* exist and we haven't said ok to overwrite, error
   write(err_string, *) 'file ', trim(output_obs_def_mod_file), &
      ' exists and will not be overwritten: Please remove or rename'
   call error_handler(E_ERR, 'preprocess', err_string, &
      source, revision, revdate)
endif

if(.not. file_exist(output_obs_qty_mod_file) .or. overwrite_output) then
   ! Open (create) the file for writing
   obs_qty_out_unit = open_file(output_obs_qty_mod_file, action='write')
else
   ! If file *does* exist and we haven't said ok to overwrite, error
   write(err_string, *) 'file ', trim(output_obs_qty_mod_file), &
      ' exists and will not be overwritten: Please remove or rename'
   call error_handler(E_ERR, 'preprocess', err_string, &
      source, revision, revdate)
endif

!______________________________________________________________________________
! Preprocessing for the obs_qty module
! Read in the quantity table and build index numbers for them. 

! Copy over lines from obs_qty_template file up to the next insertion point
linenum1 = 0
call copy_until(obs_qty_in_unit,   input_obs_qty_mod_file, insert_ints_string, linenum1, &
                obs_qty_out_unit, output_obs_qty_mod_file)

num_qtys_found = 0
qty_string(:) = 'unknown'
qty_index(:) = -1

!> for each quantity input file, read in quantities to be used.
!> duplicates are ignored (so you can list overlapping subsets 
!> of quantities).

SEARCH_QUANTITY_FILES: do j = 1, num_quantity_files
   if(file_exist(quantity_files(j))) then
      ! Open the file for reading
         in_unit = open_file(quantity_files(j), action='read')
   else
      ! If file does not exist it is an error
      write(err_string, *) 'quantity_files ', trim(quantity_files(j)), &
         ' does NOT exist (and must).'
      call error_handler(E_ERR, 'preprocess', err_string, &
         source, revision, revdate)
   endif

   ! Read until the ! BEGIN QUANTITY LIST is found
   linenum2 = 0
   call read_until(in_unit, quantity_files(j), qty2_start_string, linenum2)

   ! Subsequent lines can contain QTY_xxx lines or comments or
   ! the end string.
   DEFINE_QTYS: do
      read(in_unit, 222, IOSTAT = ierr) line
      ! If end of file, input file is incomplete or weird stuff happened
      if(ierr /= 0) then
         write(err_string, *) 'file ', trim(quantity_files(j)), &
            ' does NOT contain ', qty2_end_string
         call error_handler(E_ERR, 'preprocess', err_string, &
            source, revision, revdate)
      endif
      linenum2 = linenum2 + 1

      ! Look for the ! END QTY LIST in the current line
      test = adjustl(line)
      if(test == qty2_end_string) exit DEFINE_QTYS
   
      ! All lines between start/end must be comments or QTY_xxx strings
      ! Format:  ! QTY_string  ! optional comment
      !    or :  ! comment
   
      ! Get rid of the leading comment and any subsequent whitespace
      if (line(1:1) /= '!') then 
         err_string = 'line must begin with !'
         call quantity_error(err_string, line, quantity_files(j), linenum2)
      endif
   
      test = adjustl(line(2:))
      total_len = len(test)
   
      ! find start of word.  if comment found first, punt
      do k = 1,  total_len
         l_string = k
         if(test(k:k) == '!') cycle DEFINE_QTYS
         if(test(k:k) /= ' ') exit
      end do

      ! If there is no QTY on the line, cycle
      l_string = index(test(l_string:), 'QTY_')
      if (l_string == 0) cycle DEFINE_QTYS
   
      ! find end of word by looking one space ahead
      do k = l_string+1, total_len
         l2_string = k - 1
         if(test(k:k) == ' ' .or. test(k:k) == '!') exit
         if (k == total_len) l2_string = k 
      end do
   
      temp_qty = test(l_string:l2_string)
   
      ! this loop adds a new quantity to the string array if it's new.
      duplicate = .false.
      do i=1, num_qtys_found
         if (qty_string(i) == temp_qty) then
            duplicate = .true.
            exit
         endif
      end do
   
      if (.not. duplicate) then
         num_qtys_found = num_qtys_found + 1
         qty_string(num_qtys_found) = temp_qty
         qty_index(num_qtys_found) = num_qtys_found
      endif
   
   end do DEFINE_QTYS

   ! Close this quantity file
   close(in_unit)
end do SEARCH_QUANTITY_FILES

!______________________________________________________________________________
! Preprocessing for the obs_qty module
! Get all the type/qty strings from all of the obs_def files 
! up front and then insert stuff.  Easier to error check and combine
! duplicate quantities.

! Initial number of types is 0, qtys is based on previous section
num_types_found = 0

SEARCH_OBS_DEF_FILES: do j = 1, num_obs_type_files
   if(file_exist(obs_type_files(j))) then
      ! Open the file for reading
         in_unit = open_file(obs_type_files(j), action='read')
   else
      ! If file does not exist it is an error
      write(err_string, *) 'obs_type_files ', trim(obs_type_files(j)), &
         ' does NOT exist (and must).'
      call error_handler(E_ERR, 'preprocess', err_string, &
         source, revision, revdate)
   endif

   ! Read until the ! BEGIN QUANTITY LIST is found
   linenum2 = 0
   call read_until(in_unit, obs_type_files(j), qty_start_string, linenum2)

   ! Subsequent lines contain the type_identifier (same as type_string), and
   ! qty_string separated by commas, and optional usercode flag
   EXTRACT_TYPES: do
      read(in_unit, 222, IOSTAT = ierr) line
      ! If end of file, input file is incomplete or weird stuff happened
      if(ierr /= 0) then
         write(err_string, *) 'file ', trim(obs_type_files(j)), &
            ' does NOT contain ', qty_end_string
         call error_handler(E_ERR, 'preprocess', err_string, &
            source, revision, revdate)
      endif
      linenum2 = linenum2 + 1

      ! Look for the ! END QUANTITY LIST in the current line
      test = adjustl(line)
      if(test == qty_end_string) exit EXTRACT_TYPES

      ! All lines between start/end must be type/qty lines.
      ! Format:  ! type_string, qty_string [, COMMON_CODE]

      ! Get rid of the leading comment and any subsequent whitespace
      if (line(1:1) /= '!') then 
         err_string = 'line must begin with !'
         call typeqty_error(err_string, line, obs_type_files(j), linenum2)
      endif

      test = adjustl(line(2:))
      total_len = len(test)

      ! Compute the length of the type_string by seeking comma
      do k = 1, total_len
         l_string = k - 1
         if(test(k:k) == ',') exit
      end do

      ! comma not found? (first one is required)
      if (l_string == total_len - 1) then
         err_string = 'strings must be separated by commas'
         call typeqty_error(err_string, line, obs_type_files(j), linenum2)
      endif

      ! save results in temp vars for now, so we can check for
      ! duplicates (not allowed in types) or duplicates (which are
      ! expected in quantities)
      temp_type = adjustl(test(1:l_string))

      ! check for another comma before end of line (not mandatory)
      do k = l_string+2, total_len
         l2_string = k - 1
         if(test(k:k) == ',') exit
      end do

      ! not found?  ok, then qty is remaining part of string
      if (l2_string == total_len - 1) then
         temp_qty = adjustl(test(l_string + 2:))
         temp_user = .true.
      else
         ! another comma found, need to have COMMON_CODE on rest of line
         test2 = adjustl(test(l2_string+2:))
         if (test2(1:11) /= 'COMMON_CODE') then
            err_string = 'if third word present on line, it must be COMMON_CODE'
            call typeqty_error(err_string, line, obs_type_files(j), linenum2)
         endif

         temp_qty = adjustl(test(l_string + 2:l2_string))
         temp_user = .false.
         
      endif

      if (temp_user) file_has_usercode(j) = .true.

      !FIXME: does not correctly flag: !type  qty, COMMON_CODE 
      ! as an error (note missing comma between type and qty)
      ! how to catch?  not allow spaces in type?
  
      ! Another type/qty line; increment the type count.  Check the qtys
      ! list for repeated occurances first before deciding this is
      ! an unrecognized quantity.

      ! only allow duplicate definitions if the qtys match and there is
      ! no usercode in either the original definition or this one
      do i=1, num_types_found
         if (type_string(i) == temp_type) then
            ! allow dups if 1) same qty 2) no usercode
            if ((qty_string(qty_index(i)) /= temp_qty) .or. &
                (temp_user) .or. (usercode(i))) then
               err_string = &
                  'Duplicate! This observation type has already been processed'
               call typeqty_error(err_string, line, obs_type_files(j), linenum2)
            else ! dup, can safely ignore?
               cycle EXTRACT_TYPES
            endif
         endif
      end do

      num_types_found = num_types_found + 1
      type_string(num_types_found) = temp_type
            usercode(num_types_found) = temp_user

      ! multiple obs types mapping to the same qty is ok.
      qty_found = .false.
      do i=1, num_qtys_found
         if (qty_string(i) == temp_qty) then
            qty_found = .true.
            exit
         endif
      end do
      if (.not. qty_found) then
         ! an error to specify a quantity not already defined
         err_string = 'Unrecognized QTY_xxx quantity'
         call typeqty_error(err_string, line, obs_type_files(j), linenum2)
      endif

   end do EXTRACT_TYPES

   ! Close this obs_def file
   close(in_unit)
end do SEARCH_OBS_DEF_FILES


! Write out the integer declaration lines for qtys

call write_separator_line(obs_qty_out_unit)
call write_blank_line(obs_qty_out_unit)
write(obs_qty_out_unit, '(A)') '! Integer definitions for DART QUANTITIES'
call write_blank_line(obs_qty_out_unit)
do i = 1, num_qtys_found 
   write(obs_qty_out_unit, '(A30,A32,A3,I5)') &
      'integer, parameter, public :: ', trim(qty_string(i)), ' = ', i-1
end do
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
      'integer, parameter, public :: ', trim(type_string(i)), ' = ', i
end do
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
                obs_qty_out_unit, output_obs_qty_mod_file)

!--

! Write out the definitions of each entry of obs_qty_names
call write_blank_line(obs_qty_out_unit)
write(line, '(A)') 'do i = 0, max_defined_quantities'
write(obs_qty_out_unit, '(A)') trim(line)
write(line, '(A)') '   obs_qty_names(i) = obs_qty_type(i, "UNKNOWN")'
write(obs_qty_out_unit, '(A)') trim(line)
write(line, '(A)') 'enddo'
write(obs_qty_out_unit, '(A)') trim(line)
call write_blank_line(obs_qty_out_unit)


! Write out the definitions of each entry of obs_qty_names
do i = 1, num_qtys_found
   write(line, '(A,I5,5A)') 'obs_qty_names(', i-1, ') = obs_qty_type(', &
      trim(qty_string(i)), ", '", trim(qty_string(i)), "')"
   write(obs_qty_out_unit, '(A)') trim(line)
end do
write(obs_qty_out_unit, '(A)') blank_line

! Write out the definitions of each entry of obs_type_info
do i = 1, num_types_found
   write(line, '(A,I5,3A)') 'obs_type_info(', i, ') = obs_type_type(', &
      trim(type_string(i)), ", &"
   write(obs_qty_out_unit, '(A)') trim(line)
   write(line, *) '   ', "'", trim(type_string(i)), "', ", &
      trim(qty_string(qty_index(i))), ', .false., .false., .false.)'
   write(obs_qty_out_unit, '(A)') trim(line)
end do
write(obs_qty_out_unit, '(A)') blank_line


! Copy over rest of lines
do
   read(obs_qty_in_unit, 222, IOSTAT = ierr) line
   ! Check for end of file
   if(ierr /=0) exit

   ! Write the line to the output file
   write(obs_qty_out_unit, '(A)') trim(line)
end do

close(obs_qty_out_unit)
!______________________________________________________________________________

!______________________________________________________________________________
! Now do the obs_def files
! Read DEFAULT file line by line and copy into output file until
! Each insertion point is found. At the insertion points, copy the
! appropriate code from each requested obs_def file into the output obs_def
! file and then proceed.

! There are seven special code sections (ITEMS) in the obs_def file at present.
! Each copies code in from the special type specific obs_qty modules
! Loop goes to N+1 so that lines after the last item are copied to the output.
ITEMS: do i = 1, 8
   READ_LINE: do
      read(obs_def_in_unit, 222, IOSTAT = ierr) line
      222 format(A256)

      ! Check for end of file (it's an error if this is before the 
      ! 7 DART ITEMS have been passed)
      if(ierr /= 0) then
         if(i < 8) then
            call error_handler(E_ERR, 'preprocess', &
               'Input DEFAULT obs_def file ended unexpectedly', &
               source, revision, revdate)
         else
            exit ITEMS
         endif
      endif

      ! Check to see if this line indicates the start of an insertion section
      test = adjustl(line)
      t_string = '! DART PREPROCESS ' // trim(preprocess_string(i)) // ' INSERTED HERE'
      if(test == t_string) exit READ_LINE

      ! Write this line into the output file
      write(obs_def_out_unit, '(A)') trim(line)

   end do READ_LINE

   ! The 'USE FOR OBS_QTY_MOD' section is handled differently; lines are not
   ! copied, they are generated based on the list of types and qtys.
   if(i == qty_item) then
      ! Create use statements for both the QTY_ quantities and the individual
      ! observation type strings.
      write(obs_def_out_unit, '(A)') separator_line
      write(obs_def_out_unit, '(A)') blank_line
      do k = 1, num_types_found
         write(obs_def_out_unit, '(A)') &
            'use obs_kind_mod, only : ' // trim(type_string(k))     !!!!! FIXME: kind -> qty
      end do
      write(obs_def_out_unit, '(A)') blank_line
      do k = 1, num_qtys_found
         write(obs_def_out_unit, '(A)') &
            'use obs_kind_mod, only : ' // trim(qty_string(k))      !!!!! FIXME: kind -> qty
      end do
      write(obs_def_out_unit, '(A)') blank_line
      write(obs_def_out_unit, '(A)') separator_line
      write(obs_def_out_unit, '(A)') blank_line
      cycle
   endif


   ! Insert the code for this ITEM from each requested obs_def 'modules'
   do j = 1, num_obs_type_files
      if (.not. file_has_usercode(j)) then
         if (i == module_item) then
            write(obs_def_out_unit, '(A)') separator_line
            write(obs_def_out_unit, 31) &
               '!No module code needed for ', &
               trim(obs_type_files(j))
            write(obs_def_out_unit, '(A)') separator_line
         endif
         !if (i == use_item) then
         !   write(obs_def_out_unit, '(A)') separator_line
         !   write(obs_def_out_unit, 31) &
         !      '!No use statements needed for ', &
         !      trim(obs_type_files(j))
         !   write(obs_def_out_unit, '(A)') separator_line
         !endif
         cycle
      endif

      ! Since there might someday be a lot of these, 
      ! open and close them each time needed
      if(file_exist(obs_type_files(j))) then
         ! Open the file for reading
         in_unit = open_file(obs_type_files(j), action='read')
      else
         ! If file does not exist it is an error
         write(err_string, *) 'obs_type_files ', trim(obs_type_files(j)), &
            ' does NOT exist.'
         call error_handler(E_ERR, 'preprocess', err_string, &
            source, revision, revdate)
      endif

      ! Read until the appropriate ITEM # label is found in the input 
      ! for this obs_type
      FIND_ITEM: do

         read(in_unit, 222, IOSTAT = ierr) line
         ! If end of file, input file is incomplete or weird stuff happened
         if(ierr /=0) then
            write(err_string, *) 'file ', trim(obs_type_files(j)), &
               ' does NOT contain ! BEGIN DART PREPROCESS ', &
               trim(preprocess_string(i))
            call error_handler(E_ERR, 'preprocess', err_string, &
               source, revision, revdate)
         endif

         ! Look for the ITEM flag
         test = adjustl(line)
         t_string = '! BEGIN DART PREPROCESS ' // trim(preprocess_string(i))
         if(test == t_string) exit FIND_ITEM

      end do FIND_ITEM
      
      ! decoration or visual separation, depending on your viewpoint
      if (i == module_item) then
         write(obs_def_out_unit, '(A)') separator_line
         write(obs_def_out_unit, 31) '! Start of code inserted from ', &
            trim(obs_type_files(j))
         write(obs_def_out_unit, '(A)') separator_line
         write(obs_def_out_unit, '(A)') blank_line
         31 format(2A)
      endif

      ! Copy all code until the end of item into the output obs_def file
      COPY_ITEM: do
         read(in_unit, 222, IOSTAT = ierr) line
         ! If end of file, input file is incomplete or weird stuff happened
         if(ierr /=0) then
            write(err_string, *) 'file ', trim(obs_type_files(j)), &
               ' does NOT contain ! END DART PREPROCESS ', &
               trim(preprocess_string(i))
            call error_handler(E_ERR, 'preprocess', err_string, &
               source, revision, revdate)
         endif

         ! Look for the ITEM flag
         test = adjustl(line)
         t_string = '! END DART PREPROCESS ' // trim(preprocess_string(i))
         if(test == t_string) exit COPY_ITEM

         ! Write the line to the output obs_def_mod.f90 file
         ! Module code, if present, is copied verbatim.  
         ! All other code sections are preceeded by a ! in col 1
         ! so it must be stripped off.
         if (i == module_item) then
            write(obs_def_out_unit, '(A)') trim(line)
         else
            write(obs_def_out_unit, '(A)') trim(line(2:))
         endif
      end do COPY_ITEM

      ! decoration or visual separation, depending on your viewpoint
      if (i == module_item) then
         write(obs_def_out_unit, '(A)') blank_line
         write(obs_def_out_unit, '(A)') separator_line
         write(obs_def_out_unit, 31) '! End of code inserted from ', &
            trim(obs_type_files(j))
         write(obs_def_out_unit, '(A)') separator_line
      endif

      ! Got everything from this file, move along
      close(in_unit)
  end do
   
  ! Now check to see if this item has any types which are expecting us
  ! to automatically generate the code for them.
  do j = 1, num_types_found
     if (usercode(j)) cycle

     select case (i)
     case (get_expected_item)
        write(obs_def_out_unit, '(A)')  '      case(' // trim(type_string(j)) // ')'
        write(obs_def_out_unit, '(3A)') '         call interpolate(state_handle, ens_size, location, ', &
           trim(qty_string(qty_index(j))), ', expected_obs, istatus)'
     case (read_item, write_item, interactive_item)
        write(obs_def_out_unit, '(A)')  '   case(' // trim(type_string(j)) // ')'
        write(obs_def_out_unit, '(A)')  '      continue'
     case default
       ! nothing to do for others 
     end select
  end do

end do ITEMS

close(obs_def_out_unit)

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

call do_until(iunit, iname, stop_string, linenum, .false., 0, '')

end subroutine read_until

!------------------------------------------------------------------------------

subroutine copy_until(iunit, iname, stop_string, linenum, ounit, oname)

integer, intent(in) :: iunit
character(len=*), intent(in) :: iname
character(len=*), intent(in) :: stop_string
integer, intent(inout) :: linenum
integer, intent(in) :: ounit
character(len=*), intent(in) :: oname

call do_until(iunit, iname, stop_string, linenum, .true., ounit, oname)

end subroutine copy_until

!------------------------------------------------------------------------------

subroutine do_until(iunit, iname, stop_string, linenum, docopy, ounit, oname)

integer, intent(in) :: iunit
character(len=*), intent(in) :: iname
character(len=*), intent(in) :: stop_string
integer, intent(inout) :: linenum
logical, intent(in) :: docopy
integer, intent(in) :: ounit
character(len=*), intent(in) :: oname

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
   if(line == stop_string) exit FIND_NEXT

   ! if doing a copy, write it verbatim to the output file
   if (docopy) then
       write(ounit, full_line_out, IOSTAT = ierr) trim(line)
       if(ierr /= 0) then
          write(err_string,  *) 'Write error, returned code = ', ierr
          write(err_string2, *) 'writing file ', trim(oname)
          call error_handler(E_ERR, 'preprocess', err_string, &
                source, revision, revdate, text2=err_string2)
       endif
   endif


end do FIND_NEXT

end subroutine do_until

!------------------------------------------------------------------------------

subroutine typeqty_error(errtext, line, file, linenum)
 character(len=*), intent(in) :: errtext, line, file
 integer, intent(in) :: linenum

call error_handler(E_MSG, 'preprocess error:', &
   'obs_def file has bad Type/Qty line')
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
   'obs_def file has bad Quantity line')
call error_handler(E_MSG, 'preprocess error:', errtext)
call error_handler(E_MSG, 'expected input:', &
   '! QTY_xxx    or  ! comment only ')
write(err_string, '(2A,I5)') trim(file), ", line number", linenum
call error_handler(E_MSG, 'bad file:', err_string)
call error_handler(E_MSG, 'bad line contents:', line)
write(err_string, *) 'See msg lines above for error details'
call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)

end subroutine quantity_error

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

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

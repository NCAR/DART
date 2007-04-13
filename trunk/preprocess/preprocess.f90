! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program preprocess

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Takes a list of observation kind module path names. These modules contain
! six fragements of standard F90 that may be required to implement forward 
! observation operators for DART. The six sections are retrieved from the file
! by the preprocessor and then inserted into the appropriate blank in the
! DEFAULT_obs_def_mod.F90. The final obs_def_mod.f90 that is created contains
! the default code plus all the code required from the selected observation
! kind modules. Preprocess also inserts the required identifier and string
! for this observation kind into the DEFAULT_obs_kind_mod.f90 which is written
! out to obs_kind_mod.f90.

! NEED TO ADD IN ALL THE ERROR STUFF

use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, file_exist, &
                             open_file, logfileunit, initialize_utilities, timestamp, &
                             find_namelist_in_file, check_namelist_read

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! Pick something ridiculously large and forget about it (lazy)
integer, parameter   :: max_kinds = 10000
character(len = 256) :: line, test, kind_string(max_kinds), &
                        raw_kind_item(max_kinds), t_string
integer              :: iunit, ierr, io, i, j, k, l_kind_string
integer              :: num_kinds_found
character(len = 169) :: err_string

! List of the DART PREPROCESS strings
character(len = 29) :: preprocess_string(7) = (/'USE FOR OBS_KIND_MOD         ', &
                                                'USE OF SPECIAL OBS_DEF MODULE', &
                                                'GET_EXPECTED_OBS_FROM_DEF    ', &
                                                'READ_OBS_DEF                 ', &
                                                'WRITE_OBS_DEF                ', &
                                                'INTERACTIVE_OBS_DEF          ', &
                                                'THE SEVENTH ONE IS UNDEFINED '/)

integer :: num_input_files = 0
integer :: obs_def_in_unit, obs_def_out_unit, obs_kind_in_unit, obs_kind_out_unit, in_unit

! The namelist reads in a sequence of path_names that are absolute or
! relative to the working directory in which preprocess is being executed
! and these files are used to fill in observation kind details in
! DEFAULT_obs_def_mod.f90 and DEFAULT_obs_kind_mod.f90.
integer, parameter   :: max_input_files = 1000
character(len = 129) ::   input_obs_def_mod_file = '../../../obs_def/DEFAULT_obs_def_mod.F90'
character(len = 129) ::  output_obs_def_mod_file = '../../../obs_def/obs_def_mod.f90'
character(len = 129) ::  input_obs_kind_mod_file = '../../../obs_kind/DEFAULT_obs_kind_mod.F90'
character(len = 129) :: output_obs_kind_mod_file = '../../../obs_kind/obs_kind_mod.f90'
character(len = 129) :: input_files(max_input_files) = 'null'

namelist /preprocess_nml/ input_obs_def_mod_file, input_obs_kind_mod_file, &
   output_obs_def_mod_file, output_obs_kind_mod_file, input_files

!Begin by reading the namelist
call initialize_utilities('preprocess')
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "preprocess_nml", iunit)
read(iunit, nml = preprocess_nml, iostat = io)
call check_namelist_read(iunit, io, "preprocess_nml")

! Output the namelist file information
write(logfileunit, *) 'Path names of default obs_def and obs_kind modules'
write(*, *) 'Path names of default obs_def and obs_kind modules'
write(logfileunit, *) trim(input_obs_def_mod_file)
write(logfileunit, *) trim(input_obs_kind_mod_file)
write(*, *) trim(input_obs_def_mod_file)
write(*, *) trim(input_obs_kind_mod_file)

write(logfileunit, *) 'Path names of output obs_def and obs_kind modules'
write(*, *) 'Path names of output obs_def and obs_kind modules'
write(logfileunit, *) trim(output_obs_def_mod_file)
write(logfileunit, *) trim(output_obs_kind_mod_file)
write(*, *) trim(output_obs_def_mod_file)
write(*, *) trim(output_obs_kind_mod_file)

! A path for the default files is required. Have an error is these are still null.
if(input_obs_def_mod_file == 'null') &
   call error_handler(E_ERR, 'preprocess', 'Namelist must provide input_obs_def_mod_file', &
      source, revision, revdate)
if(input_obs_kind_mod_file == 'null') &
   call error_handler(E_ERR, 'preprocess', 'Namelist must provide input_obs_kind_mod_file', &
      source, revision, revdate)

! A path for the output files is required. Have an error is these are still null.
if(output_obs_def_mod_file == 'null') &
   call error_handler(E_ERR, 'preprocess', 'Namelist must provide output_obs_def_mod_file', &
      source, revision, revdate)
if(output_obs_kind_mod_file == 'null') &
   call error_handler(E_ERR, 'preprocess', 'Namelist must provide output_obs_kind_mod_file', &
      source, revision, revdate)

write(logfileunit, *) 'INPUT obs_kind files follow:'
write(*, *) 'INPUT obs_kind files follow:'

do i = 1, max_input_files
   if(input_files(i) == 'null') exit
   write(logfileunit, *) trim(input_files(i))
   write(*, *) trim(input_files(i))
   num_input_files= i
end do

! Try to open the DEFAULT and OUTPUT files
! DEFAULT files must exist or else an error
if(file_exist(trim(input_obs_def_mod_file))) then
   ! Open the file for reading
   obs_def_in_unit = open_file(input_obs_def_mod_file)
else
   ! If file does not exist it is an error
   write(err_string, *) 'file ', trim(input_obs_def_mod_file), ' does not exist'                   
   call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)            
endif

if(file_exist(trim(input_obs_kind_mod_file))) then
   ! Open the file for reading
   obs_kind_in_unit = open_file(input_obs_kind_mod_file)
else
   ! If file does not exist it is an error
   write(err_string, *) 'file ', trim(input_obs_kind_mod_file), ' does not exist'                   
   call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)            
endif

! Output files must NOT EXIST or else an error
if(.not. file_exist(trim(output_obs_def_mod_file))) then
   ! Open the file for reading
   obs_def_out_unit = open_file(output_obs_def_mod_file)
else
   ! If file does not exist it is an error
   write(err_string, *) 'file ', trim(output_obs_def_mod_file), ' exists: Please Rename'                   
   call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)
endif

if(.not. file_exist(trim(output_obs_kind_mod_file))) then
   ! Open the file for reading
   obs_kind_out_unit = open_file(output_obs_kind_mod_file)
else
   ! If file does not exist it is an error
   write(err_string, *) 'file ', trim(output_obs_kind_mod_file), ' exists: Please Rename'                   
   call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)            
endif

!_______________________________________________________________________________________
! Now do preprocessing for the obs_kind module
! Easiest to get the three strings required from all of the obs_kind files up front and
! then insert stuff

! Initial number of kinds is 0
num_kinds_found = 0

SEARCH_INPUT_FILES: do j = 1, num_input_files
   if(file_exist(trim(input_files(j)))) then
      ! Open the file for reading
         in_unit = open_file(input_files(j))
      else
      ! If file does not exist it is an error
      write(err_string, *) 'input_files ', trim(input_files(j)), ' does NOT exist.'
      call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)            
   endif

   ! Read until the ! BEGIN KIND LIST is found
   FIND_KIND_LIST: do

      read(in_unit, 222, IOSTAT = ierr) line
      ! If end of file, then the input file is incomplete or weird stuff has happened
      if(ierr /=0) then
         write(err_string, *) 'file ', trim(input_files(j)), &
            ' does NOT contain ! BEGIN DART PREPROCESS KIND LIST'
         call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)
      endif

      ! Look for the ! BEGIN KIND LIST
      test = adjustl(line)
      if(test(1:33) == '! BEGIN DART PREPROCESS KIND LIST') exit FIND_KIND_LIST
   end do FIND_KIND_LIST

   ! Subsequent lines contain the kind_identifier (same as kind_string), and raw_kind_ident separated by commas
   EXTRACT_KINDS: do
      read(in_unit, 222, IOSTAT = ierr) line
      ! If end of file, then the input file is incomplete or weird stuff has happened
      if(ierr /=0) then
         write(err_string, *) 'file ', trim(input_files(j)), &
            ' does NOT contain ! END DART PREPROCESS KIND LIST'
         call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)
      endif

      ! Look for the ! END KIND LIST
      test = adjustl(line)
      if(test(1:31) == '! END DART PREPROCESS KIND LIST') exit EXTRACT_KINDS

      ! Found a kind; increment the count
      num_kinds_found = num_kinds_found + 1
      ! Otherwise this line should contain kind_identifier (same as kind_string),  raw_kind_item with leading comment
      ! Get rid of the leading comment and subsequent space
      test = adjustl(line(2:))
      ! Compute the length of the kind_string by seeking comma
      do k = 1, 256
         l_kind_string = k - 1
         if(test(k:k) == ',') exit 
      end do
      kind_string(num_kinds_found) = adjustl(test(1:l_kind_string))
  
      ! Finally get the raw_kind_item
      raw_kind_item(num_kinds_found) = adjustl(test(l_kind_string + 2:))
         
   end do EXTRACT_KINDS

   ! Close this obs_kind file 
   close(in_unit)
end do SEARCH_INPUT_FILES

! A list of num_kinds_found kinds has been found, now need to put code into the obs_kind_mod
! Begin by copying over lines until first insertion point is found
do
   read(obs_kind_in_unit, 222, IOSTAT = ierr) line
   ! Check for end of file 
   if(ierr /=0) then
      call error_handler(E_ERR, 'preprocess', 'Input DEFAULT obs_kind file ended unexpectedly', &
         source, revision, revdate)
   endif

   ! Is this the place to start writing preprocessor stuff
   test = adjustl(line)
   if(test(1:38) == '! DART PREPROCESS PUBLIC INSERTED HERE') exit
  
   ! Write the line to the output file
   write(obs_kind_out_unit, 21) trim(line)
end do

! Loop to write out all the public declarations
do i = 1, num_kinds_found
   write(line, *) 'public :: ', trim(kind_string(i))
   write(obs_kind_out_unit, 51) trim(adjustl(line))
   51 format(A)
end do

! Copy over lines up to the next insertion point
do
   read(obs_kind_in_unit, 222, IOSTAT = ierr) line
   ! Check for end of file 
   if(ierr /=0) then
      call error_handler(E_ERR, 'preprocess', 'Input DEFAULT obs_kind file ended unexpectedly', &
         source, revision, revdate)
   endif

   ! Is this the place to start writing preprocessor stuff
   test = adjustl(line)
   if(test(1:51) == '! DART PREPROCESS INTEGER DECLARATION INSERTED HERE') exit
  
   ! Write the line to the output file
   write(obs_kind_out_unit, 21) trim(line)
end do

! Write out the integer declaration lines
do i = 1, num_kinds_found
   write(line, *) 'integer, parameter :: ', trim(adjustl(kind_string(i))), ' = ', i
   write(obs_kind_out_unit, 51) trim(adjustl(line))
end do

! Write out the max_obs_kinds, too
write(line, *) 'integer, parameter :: max_obs_kinds = ', num_kinds_found
write(obs_kind_out_unit, 51) trim(adjustl(line))

! Copy over lines up to the next insertion point
do
   read(obs_kind_in_unit, 222, IOSTAT = ierr) line
   ! Check for end of file 
   if(ierr /=0) then
      call error_handler(E_ERR, 'preprocess', 'Input DEFAULT obs_kind file ended unexpectedly', &
         source, revision, revdate)
   endif

   ! Is this the place to start writing preprocessor stuff
   test = adjustl(line)
   if(test(1:51) == '! DART PREPROCESS OBS_KIND_INFO INSERTED HERE') exit
  
   ! Write the line to the output file
   write(obs_kind_out_unit, 21) trim(line)
end do

! Write out the definitions of each entry of obs_kind_info
do i = 1, num_kinds_found
   write(line, *) 'obs_kind_info(', i, ') = obs_kind_type(', trim(kind_string(i)), ", '", trim(kind_string(i)), "', &"
   write(obs_kind_out_unit, 21) trim(line)
   write(line, *) '   ', trim(raw_kind_item(i)), ', .false., .false.)'
   write(obs_kind_out_unit, 21) trim(line)
end do

   
! Copy over rest of lines
do
   read(obs_kind_in_unit, 222, IOSTAT = ierr) line
   ! Check for end of file 
   if(ierr /=0) exit
  
   ! Write the line to the output file
   write(obs_kind_out_unit, 21) trim(line)
end do
!_______________________________________________________________________________________

!_______________________________________________________________________________________
! Now do the obs_def files
! Read DEFAULT file line by line and copy into output file until
! Each insertion point is found. At the insertion points, copy the
! appropriate code from each requested obs_kind file into the output obs_def
! file and then proceed.

! There are five special code sections (ITEMS) in the obs_def file at present
! That copy code in from the special type specific obs_kind modules
! Loop goes to 6 so that stuff after the last item is also copied to the final obs_def_mod.f90
ITEMS: do i = 1, 7
   READ_LINE: do
      read(obs_def_in_unit, 222, IOSTAT = ierr) line
      222 format(A256)

      ! Check for end of file (it's an error if this before all 5 DART ITEMS have been passed)
      if(ierr /=0) then
         if(i < 7) then
            call error_handler(E_ERR, 'preprocess', 'Input DEFAULT obs_def file ended unexpectedly', &
               source, revision, revdate)
         else
            exit ITEMS  
         endif
      endif

      ! Write this line into the output file
      write(obs_def_out_unit, 21) trim(line)
      21 format(A)

      ! Check to see if this line indicates the start of an insertion section
      test = adjustl(line)
      t_string = '! DART PREPROCESS ' // trim(preprocess_string(i)) // ' INSERTED HERE'
      if(trim(test) == trim(t_string)) exit READ_LINE
   end do READ_LINE

   ! Insert the code for this ITEM from each of the requested obs_kind 'modules'
   ! The first Entry does NOT copy code but creates code from the kind list
   if(i /= 1) then 
      do j = 1, num_input_files
         ! Since there might someday be a lot of these, open and close them each time needed
         if(file_exist(trim(input_files(j)))) then
            ! Open the file for reading
            in_unit = open_file(input_files(j))
         else
            ! If file does not exist it is an error
            write(err_string, *) 'input_files ', trim(input_files(j)), ' does NOT exist.'
            call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)            
         endif

         ! Read until the appropriate ITEM # label is found in the input for this obs_kind
         FIND_ITEM: do

            read(in_unit, 222, IOSTAT = ierr) line
            ! If end of file, then the input file is incomplete or weird stuff has happened
            if(ierr /=0) then
               write(err_string, *) 'file ', trim(input_files(j)), &
                  ' does NOT contain ! BEGIN DART PREPROCESS ', trim(preprocess_string(i))
               call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)
            endif

            ! Look for the ITEM flag
            test = adjustl(line)
            t_string = '! BEGIN DART PREPROCESS ' // trim(preprocess_string(i))
            if(trim(test) == trim(t_string)) exit FIND_ITEM

         end do FIND_ITEM

         ! Now read in all the code until the end of this item and copy it into the output obs_def file
         COPY_ITEM: do
            read(in_unit, 222, IOSTAT = ierr) line
            ! If end of file, then the input file is incomplete or weird stuff has happened
            if(ierr /=0) then
               write(err_string, *) 'file ', trim(input_files(j)), &
                  ' does NOT contain ! END DART PREPROCESS ', trim(preprocess_string(i))
               call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)
            endif

            ! Look for the ITEM flag
            test = adjustl(line)
            t_string = '! END DART PREPROCESS ' // trim(preprocess_string(i))
            if(trim(test) == trim(t_string)) exit COPY_ITEM
     
            ! Write the line to the output obs_def_mod.f90 file
            write(obs_def_out_unit, 21) trim(line(2:))
         end do COPY_ITEM
      
         ! Got everything from this file, move along
         close(in_unit)
      end do

   ! If this is the first substitution spot in the DEFAULT file, create new code from kind list
   ! in the individual special obs_def file
   else
      do k = 1, num_kinds_found
         write(obs_def_out_unit, 21) 'use obs_kind_mod, only : ' // trim(kind_string(k))
      end do
   endif

end do ITEMS

!_______________________________________________________________________________________


call timestamp(source,revision,revdate,'end') ! closes the log file.

end program preprocess

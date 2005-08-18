! Data Assimilation Research Testbed -- DART
! Copyright 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program preprocess

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

! NEED TO ADD IN ALL THE ERROR STUFF

use        types_mod, only : r8, missing_i, missing_r8, RAD2DEG
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, file_exist, &
                             open_file, check_nml_error, logfileunit, close_file, &
                             initialize_utilities, timestamp

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

character(len = 256) :: line, test, ifdef_kind
integer :: iunit, ierr, io, i, j, status
integer :: num_kind_assimilate = 0, num_kind_evaluate = 0
character(len = 169) :: err_string

integer, parameter :: max_obs_kinds = 1000
integer :: input_unit, output_unit

! Namelist array to identify input and output file
character(len = 129) :: input_file = "no_default_input_file"
character(len = 129) :: output_file = "no_default_input_file"

namelist /preprocess_nml/ input_file, output_file

! Namelist array to turn on any requested observation types
character(len = 129) :: assimilate_these_obs_types(max_obs_kinds) = 'null'
character(len = 129) :: evaluate_these_obs_types(max_obs_kinds) = 'null'

namelist /obs_def_nml/ assimilate_these_obs_types, evaluate_these_obs_types



!Begin by reading the namelist
call initialize_utilities('preprocess')
call register_module(source, revision, revdate)

! Read the preprocess_nml to get the preprocessor observation types to keep
if(file_exist('input.nml')) then
   iunit = open_file(fname = 'input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = preprocess_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'preprocess_nml')
   enddo
 11 continue
   call close_file(iunit)
endif

! Open the input file on input_unit
if(file_exist(trim(adjustl(input_file)))) then
   input_unit = open_file(input_file)
else
   ! If file doesn't exist we have an error
   write(err_string, *) 'input_file ', trim(input_file), ' does not exist'
   call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)
endif

! Open the output file on output_unit
if(file_exist(trim(adjustl(output_file)))) then
   ! If file DOES exist we have an error; avoid over-writing without notifying
   write(err_string, *) 'output_file ', trim(output_file), ' exists: MOVE TO AVOID OVERWRITE'
   call error_handler(E_ERR, 'preprocess', err_string, source, revision, revdate)
else
   output_unit = open_file(output_file)
endif


! Read the obs_def_nml to get the preprocessor observation types to keep
if(file_exist('input.nml')) then
   iunit = open_file(fname = 'input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = obs_def_nml, iostat = io, end = 44)
      ierr = check_nml_error(io, 'obs_def_nml')
   enddo
 44 continue
   call close_file(iunit)
endif

call error_handler(E_MSG,'initialize_module','obs_def_nml values are',' ',' ',' ')
write(logfileunit, *) 'assimilate_these_obs_types'
write(*, *) 'ASSIMILATE_THESE_OBS_TYPES'
do i = 1, max_obs_kinds
   if(assimilate_these_obs_types(i) == 'null') goto 22
   write(logfileunit, *) trim(assimilate_these_obs_types(i))
   write(*, *) trim(assimilate_these_obs_types(i))
   num_kind_assimilate = i
end do
22 write(logfileunit, *) 'evaluate_these_obs_types'
write(*, *) 'EVALUATE_THESE_OBS_TYPES'
do i = 1, max_obs_kinds
   if(evaluate_these_obs_types(i) == 'null') goto 33
   write(logfileunit, *) trim(evaluate_these_obs_types(i))
   write(*, *) trim(evaluate_these_obs_types(i))
   num_kind_evaluate = i
end do
33 continue




! Meaning of the status variable
! status = 1: Reading normally
! status = 2: In the middle of an unimplemented ifdef, don't write

status = 1


! Loop to read each line of a standard input fortran 90 program
READ_LINE: do 
   read(input_unit, 222, IOSTAT = ierr) line
   222 format(A256)
 
   ! Check for end of file
   if(ierr /= 0) exit

   ! Remove any line that starts with #ifdef or #IFDEF
   test = adjustl(line)

   ! If status is 1 (normally reading)
   if(status == 1) then
      if(test(1:6) == '#IFDEF' .or. test(1:6) == '#ifdef') then
         ! Determine if the ifdef is one that is turned on for the obs_kinds in use
         ifdef_kind = trim(adjustl(test(7:)))
         do j = 1, num_kind_assimilate
            ! Is it being assimilated?
            if(ifdef_kind == trim(adjustl(assimilate_these_obs_types(j)))) then
               cycle READ_LINE    
            endif
         end do
         do j = 1, num_kind_evaluate
            ! Is it being evaluated?
            if(ifdef_kind == evaluate_these_obs_types(j)) then
               cycle READ_LINE    
            endif
         end do
          
         ! An unimplemented ifdef has been encountered, don't output ifdef line
         status = 2
      else if(test(1:6) == '#ENDIF' .or. test(1:6) == '#endif') then
         ! Change status to 1, don't output line
         status = 1
      else
         write(output_unit, 21) trim(line)
         21 format(A)
      endif
   else
      ! Otherwise status is 2; change to 1 if endif found
      if(test(1:6) == '#ENDIF' .or. test(1:6) == '#endif') status = 1
   endif
end do READ_LINE

call timestamp(source,revision,revdate,'end') ! closes the log file.

end program preprocess

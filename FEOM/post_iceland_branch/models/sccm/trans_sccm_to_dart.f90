! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program trans_sccm_to_dart

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

! Program to integrate assimilation model forward for asynchronous filter
! execution.

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time, GREGORIAN, set_calendar_type
use    utilities_mod, only : open_file, close_file, &
                             file_exist, initialize_utilities, register_module, &
                             error_handler, logfileunit, E_ERR, E_MSG, timestamp
use  assim_model_mod, only : static_init_assim_model, &
   get_model_size, open_restart_write, close_restart, awrite_state_restart

use model_mod, only : write_sccm_state, read_sccm_state

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(time_type)         :: target_time, model_time
real(r8), allocatable   :: x(:)
integer :: iunit, io, model_size
character(len=129) :: err_string, nml_string

!----------------------------------------------------------------
! Namelist input with default values
!
character(len = 129) :: sccm_file_name = "dart_data.dat", &
                        dart_file_name = 'dart_file_out'

namelist /trans_sccm_to_dart_nml/ sccm_file_name, dart_file_name
!----------------------------------------------------------------

call initialize_utilities('trans_sccm_to_dart')
call register_module(source,revision,revdate)

! Begin by reading the namelist input
if(file_exist('trans_sccm_to_dart.nml')) then
   iunit = open_file('trans_sccm_to_dart.nml', action = 'read')
   read(iunit, nml = trans_sccm_to_dart_nml, iostat = io)
   if(io /= 0) then
      ! A non-zero return means a bad entry was found for this namelist
      ! Reread the line into a string and print out a fatal error message.
      BACKSPACE iunit
      read(iunit, '(A)') nml_string
      write(err_string, *) 'INVALID NAMELIST ENTRY: ', trim(adjustl(nml_string))
      call error_handler(E_ERR, 'trans_sccm_to_dart:&trans_sccm_to_dart_nml problem', &
                         err_string, source, revision, revdate)
   endif
   call close_file(iunit)
endif

! Record the namelist values used for the run ...
call error_handler(E_MSG,'trans_sccm_to_dart','trans_sccm_to_dart_nml values are',' ',' ',' ')
write(logfileunit, nml=trans_sccm_to_dart_nml)
write(     *     , nml=trans_sccm_to_dart_nml)

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()
! Read in the target time
allocate(x(model_size))

!------------------- Read sccm file ----------------------

call set_calendar_type(GREGORIAN)
iunit = open_file(sccm_file_name)
call read_sccm_state(x, model_time, iunit)
call close_file(iunit)

! Write out a dart file
iunit = open_restart_write(dart_file_name)
call awrite_state_restart(model_time, x, iunit)
call close_restart(iunit)


end program trans_sccm_to_dart

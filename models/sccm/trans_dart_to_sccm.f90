! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program trans_dart_to_sccm

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

! Program to integrate assimilation model forward for asynchronous filter
! execution.

use        types_mod, only : r8
use time_manager_mod, only : time_type, GREGORIAN, set_calendar_type, get_time, &
                             operator(-)
use    utilities_mod, only : open_file, close_file, &
                             file_exist, initialize_utilities, register_module, &
                             error_handler, logfileunit, E_ERR, E_MSG, timestamp, get_unit
use  assim_model_mod, only : static_init_assim_model, get_model_size, &
                             open_restart_read, close_restart, aread_state_restart

use model_mod, only : write_sccm_state, read_sccm_state

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(time_type)         :: target_time, model_time, delta_time
real(r8), allocatable   :: x(:)
integer :: iunit, io, model_size, seconds, days, delta_seconds
character(len=129) :: err_string, nml_string

!----------------------------------------------------------------
! Namelist input with default values
!
character(len = 129) :: sccm_file_name = "dart_data.dat", &
                        dart_file_name = 'dart_file_in'

namelist /trans_dart_to_sccm_nml/ sccm_file_name, dart_file_name
!----------------------------------------------------------------

call initialize_utilities('trans_dart_to_sccm')
call register_module(source,revision,revdate)

! Begin by reading the namelist input
if(file_exist('trans_dart_to_sccm.nml')) then
   iunit = open_file('trans_dart_to_sccm.nml', action = 'read')
   read(iunit, nml = trans_dart_to_sccm_nml, iostat = io)
   if(io /= 0) then
      ! A non-zero return means a bad entry was found for this namelist
      ! Reread the line into a string and print out a fatal error message.
      BACKSPACE iunit
      read(iunit, '(A)') nml_string
      write(err_string, *) 'INVALID NAMELIST ENTRY: ', trim(adjustl(nml_string))
      call error_handler(E_ERR, 'static_init_model:&trans_dart_to_sccm_nml problem', &
                         err_string, source, revision, revdate)
   endif
   call close_file(iunit)
endif

! Record the namelist values used for the run ...
call error_handler(E_MSG,'trans_dart_to_sccm','trans_dart_to_sccm_nml values are',' ',' ',' ')
write(logfileunit, nml=trans_dart_to_sccm_nml)
write(     *     , nml=trans_dart_to_sccm_nml)


! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()
! Read in the target time
allocate(x(model_size))

!------------------- Read dart state file ----------------------

call set_calendar_type(GREGORIAN)

! Read a dart file
iunit = open_restart_read(dart_file_name)
call aread_state_restart(model_time, x, iunit, target_time)
call close_restart(iunit)

! Figure out the difference in seconds between model_time and target_time
delta_time = target_time - model_time
call get_time(delta_time, seconds, days)
delta_seconds = seconds + 3600 * days

!------------------- Write sccm file, include delta_t ----------
iunit = get_unit()
!!!iunit = open_file(sccm_file_name)
open(unit = iunit, file = sccm_file_name, delim = "none")
call write_sccm_state(x, model_time, delta_seconds, iunit)
call close_file(iunit)

end program trans_dart_to_sccm

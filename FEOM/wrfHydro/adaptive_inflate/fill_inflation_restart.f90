! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$ 

program fill_inflation_restart

! Write an inflation restart file with the right number of entries,
! based on a single inflate and standard deviation value read from the console.
! (alternatively we could read them from a namelist, but it seems like overkill.)

use types_mod,            only : r8

use utilities_mod,        only : error_handler, E_MSG,  &
                                 initialize_utilities, finalize_utilities

use ensemble_manager_mod, only : ensemble_type, write_ensemble_restart,       &
                                 init_ensemble_manager, end_ensemble_manager, &
                                 prepare_to_write_to_vars

use assim_model_mod,      only : static_init_assim_model

use            model_mod, only : get_model_size

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Read 2 floating point values from the console - the initial inflation value
! and the inflation standard deviation.  Will write an output file named
! 'inflation_ics' -- rename it to match what filter will read in.
! (the name of the inflation restart file is a namelist option.)

! Module storage for writing messages
character(len = 129) :: msgstring
character(len = 32)  :: out_file_name = 'inflate_ics'

type(ensemble_type)  :: ens_handle

real(r8) :: inf_initial, sd_initial
integer  :: ss_inflate_index    = 1 
integer  :: ss_inflate_sd_index = 2
integer  :: model_size

!============================================================================

call initialize_utilities('fill_inflation_restart')

msgstring = 'Both mean and sd will be read from console'
call error_handler(E_MSG, 'fill_inflate_restart', trim(msgstring), &
                   source, revision, revdate)
msgstring = 'Full mean and sd fields will be written to this restart file: ' // trim(out_file_name)
call error_handler(E_MSG, 'fill_inflate_restart', trim(msgstring), &
                   source, revision, revdate)

print *, 'Enter initial values for inflation mean and standard deviation:'
read *, inf_initial, sd_initial

write(msgstring, '(A, F12.6, 1X, F12.6)') 'mean and sd read from console as ', &
         inf_initial, sd_initial
call error_handler(E_MSG, 'fill_inflate_restart', trim(msgstring), &
                   source, revision, revdate)

! initialize the assim_model and model code
call static_init_assim_model()
model_size = get_model_size()

write(msgstring, *) 'Model size/restart data length = ', model_size
call error_handler(E_MSG,'',msgstring, source, revision, revdate)

call init_ensemble_manager(ens_handle, 2, model_size)
call prepare_to_write_to_vars(ens_handle)

ens_handle%vars(:, ss_inflate_index   ) = inf_initial
ens_handle%vars(:, ss_inflate_sd_index) =  sd_initial

call write_ensemble_restart(ens_handle, out_file_name, ss_inflate_index, &
                            ss_inflate_sd_index, force_single_file = .true.)

call end_ensemble_manager(ens_handle)

call finalize_utilities()

!========================================================================

end program fill_inflation_restart

! <next few lines under version control, do not edit>
! $URL$ 
! $Id$ 
! $Revision$ 
! $Date$ 

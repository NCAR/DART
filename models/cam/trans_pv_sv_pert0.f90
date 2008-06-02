! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program trans_pv_sv_pert0

! <next few lines under version control, do not edit>
! $URL: http://subversion.ucar.edu/DAReS/DART/trunk/models/cam/trans_pv_sv_time0.f90 $
! $Id: trans_pv_sv_time0.f90 2713 2007-03-26 04:09:04Z thoar $
! $Revision: 2713 $
! $Date: 2007-03-25 22:09:04 -0600 (Sun, 25 Mar 2007) $

!----------------------------------------------------------------------
! purpose: interface between CAM and DART, 
!          Modify field (esp. parameter field) values to have spread
!
! method: Read CAM 'initial' file (netCDF format).
!         Use ensemble member (in file 'element') to set a unique seed for init_ran1.
!         Change values by call to model_mod:pert_model_state.
!         Reform fields into a state vector.
!         Write out state vector in "proprietary" format for DART
!
! author: Kevin Raeder 3/21/07
!         based on trans_pv_sv_pert0 
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : nmlfileunit, error_handler, E_ERR, E_MSG, find_namelist_in_file, &
                             check_namelist_read, initialize_utilities, finalize_utilities
use        model_mod, only : model_type, init_model_instance, read_cam_init, &
                             prog_var_to_vector, pert_model_state

use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
   init_assim_model, get_model_size , set_model_state_vector, write_state_restart, &
   set_model_time, open_restart_read, open_restart_write, close_restart, &
   aread_state_restart
! Guam; move time stripping from advance_model to here
use time_manager_mod, only : time_type, read_time, set_time

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL: http://subversion.ucar.edu/DAReS/DART/trunk/models/cam/trans_pv_sv_time0.f90 $", &
   revision = "$Revision: 2713 $", &
   revdate  = "$Date: 2007-03-25 22:09:04 -0600 (Sun, 25 Mar 2007) $"

character (len = 128) :: file_name = 'caminput.nc', file_out = 'temp_ud'
! Hawaii;                          file_time = 'temp_ic'
! trans_pv_sv_pert0 should get its time from the namelist, 
! not from temp_ic, which came from filter_ics, which will not exist for a new
! set of fields comprising the state vector.

! Temporary allocatable storage to read in a native format for cam state
type(assim_model_type) :: x
type(model_type)       :: var
type(time_type)        :: model_time
real(r8), allocatable  :: x_state(:), x_temp(:)
integer                :: file_unit, x_size, iunit, io

! define exactly the same stuff as in the perfect_model_obs namelist.
logical :: start_from_restart = .false., output_restart = .false., &
           interf_provided
integer :: async = 0
integer :: init_time_days = 0, init_time_seconds = 0, output_interval = 1, &
           first_obs_days, first_obs_seconds, &
           last_obs_days,  last_obs_seconds
character(len = 129) :: restart_in_file_name  = 'perfect_ics',     &
                        restart_out_file_name = 'perfect_restart', &
                        obs_seq_in_file_name  = 'obs_seq.in',      &
                        obs_seq_out_file_name = 'obs_seq.out',     &
                        adv_ens_command       = './advance_model.csh'

namelist /perfect_model_obs_nml/                                            &
          start_from_restart, output_restart, async,                        &
          init_time_days, first_obs_days, first_obs_seconds, last_obs_days, &
          last_obs_seconds,init_time_seconds, output_interval,              &
          restart_in_file_name, restart_out_file_name,                      &
          obs_seq_in_file_name, obs_seq_out_file_name,                      &
          adv_ens_command

!-------------------------

call initialize_utilities('trans_sv_pv_pert0')

! Static init assim model calls static_init_model
PRINT*,'static_init_assim_model in trans_pv_sv'
call static_init_assim_model()

! Initialize the assim_model instance
call init_assim_model(x)

! Allocate the local state vector
x_size = get_model_size()
allocate(x_state(x_size), x_temp(x_size))

! Allocate the instance of the cam model type for storage
call init_model_instance(var)

! Read the file cam state fragments into var
call read_cam_init(file_name, var)

! transform fields into state vector for DART
call prog_var_to_vector(var, x_state)

! Replace mean/base value and/or add spread to (parameter) fields 
! according to the namelist entries pert_names, pert_sd, and pert_base_vals.
! Random value (< pert_sd) will be added to pert_base_vals.
! The seed for the random # will be derived from the ensemble member #,
! so it should be unique for each ensemble member (and variable).
! (This isn't always true for some random number generators).

call pert_model_state(x_state, x_temp, interf_provided)

! Put this in the structure
call set_model_state_vector(x, x_temp)

! Read the namelist entry
call find_namelist_in_file("input.nml", "perfect_model_obs_nml", iunit)
read(iunit, nml = perfect_model_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "perfect_model_obs_nml")

! Record the namelist values used for the run ...
write(nmlfileunit, nml=perfect_model_obs_nml)
write(     *     , nml=perfect_model_obs_nml)


call filter_set_initial_time

call set_model_time (x, model_time)
call close_restart(file_unit)

! Get channel for output 
file_unit = open_restart_write(file_out)
PRINT*,'In trans_pv_sv file_out unit = ',file_unit
PRINT*,' '
! write out state vector in "proprietary" format
call write_state_restart(x, file_unit)
call close_restart(file_unit)

call finalize_utilities()

!-------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------

subroutine filter_set_initial_time

if(init_time_days >= 0) then
   model_time = set_time(init_time_seconds, init_time_days)
else
   model_time = set_time(0, 0)
endif

end subroutine filter_set_initial_time

end program trans_pv_sv_pert0


! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module default_model_mod

!> bypass routines for all required entry points.
!> if a model has no need of a routine, use these instead.

use        types_mod,      only : r8, i8, i4, MISSING_R8

use time_manager_mod,      only : time_type, set_time

use     location_mod,      only : location_type, set_location, set_location_missing, &
                                  get_close_type, get_close_obs, get_close_state, &
                                  convert_vertical_obs, convert_vertical_state

use    utilities_mod,      only : register_module, error_handler, E_ERR, E_MSG, nmlfileunit, &
                                  do_output, find_namelist_in_file, check_namelist_read,     &
                                  do_nml_file, do_nml_term

use netcdf_utilities_mod,  only : nc_check

use ensemble_manager_mod,  only : ensemble_type

use dart_time_io_mod,      only : read_model_time, write_model_time

implicit none
private

public :: get_model_size, &
          adv_1step, &
          get_state_meta_data, &
          model_interpolate, &
          shortest_time_between_assimilations, &
          end_model, &
          static_init_model, &
          init_time, &
          fail_init_time, &
          init_conditions, &
          fail_init_conditions, &
          nc_write_model_atts, &
          nc_write_model_vars, &
          pert_model_copies, &
          get_close_obs, &
          get_close_state, &
          convert_vertical_obs, &
          convert_vertical_state, &
          read_model_time, &
          write_model_time


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

contains

!==================================================================

subroutine static_init_model()

end subroutine static_init_model

!------------------------------------------------------------------

subroutine init_conditions(x)
real(r8), intent(out) :: x(:)

! default
x = 0.0_r8

end subroutine init_conditions

!------------------------------------------------------------------

subroutine fail_init_conditions(x)
real(r8), intent(out) :: x(:)

call error_handler(E_ERR, 'init_conditions', 'this model cannot provide initial conditions', &
                   source, revision, revdate)

! default
x = 0.0_r8

end subroutine fail_init_conditions

!------------------------------------------------------------------

subroutine adv_1step(x, time)
real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

call error_handler(E_ERR, 'adv_1step', 'unable to advance model', &
                   source, revision, revdate)

end  subroutine adv_1step

!------------------------------------------------------------------

function get_model_size()
integer :: get_model_size

get_model_size = 1

end function get_model_size

!------------------------------------------------------------------

subroutine init_time(time)
type(time_type), intent(out) :: time

time = set_time(0, 0)

end subroutine init_time

!------------------------------------------------------------------

subroutine fail_init_time(time)
type(time_type), intent(out) :: time

call error_handler(E_ERR, 'init_time', 'this model cannot provide an initial time', &
                   source, revision, revdate)

time = set_time(0, 0)

end subroutine fail_init_time

!------------------------------------------------------------------

subroutine model_interpolate(state_handle, ens_size, location, obs_quantity, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: obs_quantity
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

! all fail
expected_obs(:) = MISSING_R8
istatus(:) = 1

end subroutine model_interpolate

!------------------------------------------------------------------

function shortest_time_between_assimilations()
type(time_type) :: shortest_time_between_assimilations

! default to 1 day
shortest_time_between_assimilations = set_time(0, 1)

end function shortest_time_between_assimilations

!------------------------------------------------------------------

subroutine get_state_meta_data(state_handle, index_in, location, var_type)

type(ensemble_type), intent(in)  :: state_handle !< some large models need this
integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

location = set_location_missing()
if (present(var_type)) var_type = 0    ! default variable type

end subroutine get_state_meta_data

!------------------------------------------------------------------

subroutine end_model()

end subroutine end_model

!------------------------------------------------------------------

subroutine nc_write_model_atts(ncid, domain_id) 

integer, intent(in) :: ncid
integer, intent(in) :: domain_id

end subroutine nc_write_model_atts

!------------------------------------------------------------------

subroutine nc_write_model_vars(ncid, domain_id, state_ens_handle, memberindex, timeindex)

integer,             intent(in) :: ncid      
integer,             intent(in) :: domain_id
type(ensemble_type), intent(in) :: state_ens_handle
integer, optional,   intent(in) :: memberindex
integer, optional,   intent(in) :: timeindex

end subroutine nc_write_model_vars

!--------------------------------------------------------------------

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,   intent(in) :: ens_size
real(r8),  intent(in) :: pert_amp
logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_copies

!===================================================================
! End of model_mod
!===================================================================
end module default_model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

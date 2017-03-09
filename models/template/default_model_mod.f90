! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module default_model_mod

! bypass routines for all required entry points.
! if a model has no need of a routine, use these instead.

use        types_mod,      only : r8, i8, i4, MISSING_R8

use time_manager_mod,      only : time_type, set_time

use     location_mod,      only : location_type, set_location, set_location_missing, &
                                  get_close_maxdist_init, get_close_obs_init, &
                                  loc_get_close_obs => get_close_obs, get_close_type

use    utilities_mod,      only : register_module, error_handler, E_ERR, E_MSG, nmlfileunit, &
                                  do_output, find_namelist_in_file, check_namelist_read,     &
                                  do_nml_file, do_nml_term, nc_check

use         obs_kind_mod,  only : RAW_STATE_VARIABLE

use ensemble_manager_mod,  only : ensemble_type

use dart_time_io_mod,      only : read_model_time, write_model_time

implicit none
private

public :: get_model_size, &
          adv_1step, &
          get_state_meta_data, &
          model_interpolate, &
          get_model_time_step, &
          end_model, &
          static_init_model, &
          init_time, &
          init_conditions, &
          nc_write_model_atts, &
          nc_write_model_vars, &
          pert_model_copies, &
          get_close_maxdist_init, &
          get_close_obs_init, &
          get_close_obs, &
          get_close_type, &
          vert_convert, &
          query_vert_localization_coord, &
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

subroutine model_interpolate(state_handle, ens_size, location, itype, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: itype
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

! all fail
expected_obs(:) = MISSING_R8
istatus(:) = 1

end subroutine model_interpolate

!------------------------------------------------------------------

function get_model_time_step()
type(time_type) :: get_model_time_step

! default to 1 day
get_model_time_step = set_time(0, 1)

end function get_model_time_step

!------------------------------------------------------------------

subroutine get_state_meta_data(state_handle, index_in, location, var_type)

type(ensemble_type), intent(in)  :: state_handle !< some large models need this
integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

location = set_location_missing()
if (present(var_type)) var_type = RAW_STATE_VARIABLE    ! default variable type

end subroutine get_state_meta_data

!------------------------------------------------------------------

subroutine end_model()

end subroutine end_model

!------------------------------------------------------------------

function nc_write_model_atts( ncFileID, model_mod_writes_state_variables ) result (ierr)

integer, intent(in)  :: ncFileID      ! netCDF file identifier
logical, intent(out) :: model_mod_writes_state_variables
integer              :: ierr          ! return value of function

ierr = 0                             ! assume normal termination
model_mod_writes_state_variables = .false.

end function nc_write_model_atts

!------------------------------------------------------------------

function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

ierr = 0                      ! assume normal termination

end function nc_write_model_vars

!--------------------------------------------------------------------

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,   intent(in) :: ens_size
real(r8),  intent(in) :: pert_amp
logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_copies

!--------------------------------------------------------------------

subroutine vert_convert(state_handle, location, obs_kind, istatus)

type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
integer,             intent(out) :: istatus

istatus = 0

end subroutine vert_convert

!--------------------------------------------------------------------

function query_vert_localization_coord()

integer :: query_vert_localization_coord

!> @TODO should define some parameters including something
!> like HAS_NO_VERT for this use.

query_vert_localization_coord = -1

end function query_vert_localization_coord

!--------------------------------------------------------------------

!> Pass through to the code in the locations module

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, &
                         obs_kind, num_close, close_ind, dist, state_handle)

type(ensemble_type),         intent(in)     :: state_handle
type(get_close_type),        intent(in)     :: gc
type(location_type),         intent(inout)  :: base_obs_loc, obs_loc(:)
integer,                     intent(in)     :: base_obs_kind, obs_kind(:)
integer,                     intent(out)    :: num_close, close_ind(:)
real(r8),                    intent(out)    :: dist(:)


call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                          num_close, close_ind, dist)

end subroutine get_close_obs

!--------------------------------------------------------------------

!===================================================================
! End of model_mod
!===================================================================
end module default_model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

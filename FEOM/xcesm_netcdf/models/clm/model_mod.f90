! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the wrapper for the CLM model mod to allow it to be
! used stand-alone without the other CESM components.
! See clm_model_mod.f90 for the actual model code.


! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, SECPERDAY, MISSING_R8, rad2deg, PI
use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date,                           &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)
use     location_mod, only : location_type, get_dist,                          &
                             get_close_maxdist_init, get_close_obs_init,       &
                             set_location, get_close_type,                     &
                             VERTISHEIGHT, get_location, vert_is_height,       &
                             vert_is_level, vert_is_surface, get_close_obs
use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text, do_output

use clm_model_mod, only :    &
   get_model_size => clm_get_model_size,       &
   adv_1step => clm_adv_1step,            &
   get_state_meta_data => clm_get_state_meta_data,  &
   model_interpolate => clm_model_interpolate,    &
   get_model_time_step => clm_get_model_time_step,  &
   static_init_model => clm_static_init_model,    &
   end_model => clm_end_model,            &
   init_time => clm_init_time,            &
   init_conditions => clm_init_conditions,      &
   nc_write_model_atts => clm_nc_write_model_atts,  &
   nc_write_model_vars => clm_nc_write_model_vars,  &
   pert_model_state => clm_pert_model_state,     &
   ens_mean_for_model => clm_ens_mean_for_model,   &
   to_dart_state_vector => clm_to_dart_state_vector, &
   sv_to_restart_file => clm_sv_to_restart_file,   &
   get_gridsize,             &
   get_clm_restart_filename, &
   get_state_time,           &
   get_grid_vertval,         &
   compute_gridcell_value,   &
   DART_get_var,             &
   get_model_time

use typesizes
use netcdf 

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          static_init_model,      &
          end_model,              &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model


! these routines are public because they are used by
! other programs (e.g. the converters).  the interfaces
! for these routines can be changed if needed.
public :: get_gridsize,              &
          sv_to_restart_file,        &
          get_clm_restart_filename,  &
          get_state_time,            &
          get_grid_vertval,          &
          compute_gridcell_value,    &
          DART_get_var,              &
          get_model_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

contains

!%! 
!%! !------------------------------------------------------------------
!%! !------------------------------------------------------------------
!%! 
!%! subroutine static_init_model()
!%! 
!%! call clm_static_init_model()
!%! 
!%! end subroutine static_init_model
!%! 
!%! !------------------------------------------------------------
!%! 
!%! subroutine init_conditions(x)
!%!  real(r8), intent(out) :: x(:)
!%! 
!%! call clm_init_conditions(x)
!%! 
!%! end subroutine init_conditions
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! subroutine adv_1step(x, time)
!%!  real(r8),        intent(inout) :: x(:)
!%!  type(time_type), intent(in)    :: time
!%! 
!%! call clm_adv_1step(x, time)
!%! 
!%! end subroutine adv_1step
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! function get_model_size()
!%!  integer :: get_model_size
!%! 
!%! get_model_size = clm_get_model_size()
!%! 
!%! end function get_model_size
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! subroutine init_time(time)
!%!  type(time_type), intent(out) :: time
!%! 
!%! call clm_init_time(time)
!%! 
!%! end subroutine init_time
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! subroutine model_interpolate(x, location, obs_type, interp_val, istatus)
!%!  real(r8),            intent(in) :: x(:)
!%!  type(location_type), intent(in) :: location
!%!  integer,             intent(in) :: obs_type
!%!  real(r8),           intent(out) :: interp_val
!%!  integer,            intent(out) :: istatus
!%! 
!%! call clm_model_interpolate(x, location, obs_type, interp_val, istatus)
!%! 
!%! end subroutine model_interpolate
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! function get_model_time_step()
!%!  type(time_type) :: get_model_time_step
!%! 
!%! get_model_time_step = clm_get_model_time_step()
!%! 
!%! end function get_model_time_step
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! subroutine get_state_meta_data(index_in, location, var_type)
!%!  integer,             intent(in)  :: index_in
!%!  type(location_type), intent(out) :: location
!%!  integer,             intent(out), optional :: var_type
!%! 
!%! call clm_get_state_meta_data(index_in, location, var_type)
!%! 
!%! end subroutine get_state_meta_data
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! subroutine end_model()
!%! 
!%! call clm_end_model()
!%! 
!%! end subroutine end_model
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! function nc_write_model_atts(ncFileID)
!%!  integer, intent(in)  :: ncFileID            ! netCDF file identifier
!%!  integer              :: nc_write_model_atts ! function return value
!%! 
!%! nc_write_model_atts = clm_nc_write_model_atts(ncFileID)
!%! 
!%! end function nc_write_model_atts
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! function nc_write_model_vars(ncFileID, statevec, copyindex, timeindex) 
!%!  integer,                intent(in) :: ncFileID            ! netCDF file identifier
!%!  real(r8), dimension(:), intent(in) :: statevec
!%!  integer,                intent(in) :: copyindex
!%!  integer,                intent(in) :: timeindex
!%!  integer                            :: nc_write_model_vars ! function return value
!%! 
!%! nc_write_model_vars = clm_nc_write_model_vars(ncFileID, statevec, copyindex, timeindex) 
!%! 
!%! end function nc_write_model_vars
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! subroutine pert_model_state(state, pert_state, interf_provided)
!%!  real(r8), intent(in)  :: state(:)
!%!  real(r8), intent(out) :: pert_state(:)
!%!  logical,  intent(out) :: interf_provided
!%! 
!%! call clm_pert_model_state(state, pert_state, interf_provided)
!%! 
!%! end subroutine pert_model_state
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! subroutine ens_mean_for_model(ens_mean)
!%!  real(r8), intent(in) :: ens_mean(:)
!%! 
!%! call clm_ens_mean_for_model(ens_mean)
!%! 
!%! end subroutine ens_mean_for_model
!%! 
!%! !------------------------------------------------------------------
!%! 
!%! subroutine sv_to_restart_file(state_vector, filename, statedate)
!%!  real(r8),         intent(in) :: state_vector(:)
!%!  character(len=*), intent(in) :: filename 
!%!  type(time_type),  intent(in) :: statedate
!%! 
!%! call clm_sv_to_restart_file(state_vector, filename, statedate)
!%! 
!%! end subroutine sv_to_restart_file
!%! 
!%! !------------------------------------------------------------------
!%! ! End of model_mod
!%! !------------------------------------------------------------------
!%! 

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

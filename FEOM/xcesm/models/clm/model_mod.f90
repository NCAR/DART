! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the wrapper for the CLM model mod to allow it to be
! used stand-alone without the other CESM components.
! See clm_model_mod.f90 for the actual model code.


use clm_model_mod

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


end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

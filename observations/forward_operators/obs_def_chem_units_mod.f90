! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! An example of a simple forward operator that involves more than
! just interpolating directly from a state vector in a model.
!
! This section defines a specific type in the left column and
! can be any string you want to use for an observation.  The
! right column must be a generic kind that already exists in
! the obs_kind/DEFAULT_obs_kind_mod.F90 file.

! BEGIN DART PREPROCESS KIND LIST
! IASI_CO      QTY_CO
! END DART PREPROCESS KIND LIST

! This section will be added to the main obs_def_mod.f90 that
! is going to be generated, to allow it to call the code we
! are defining here.

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_chem_units_mod, only : get_expected_chem_units
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! This section will be dropped into a large case statement in the
! main obs_def_mod.f90 code to control what happens with each
! observation type that is processed.

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!   case(IASI_CO)
!        call get_expected_chem_units(state_handle, ens_size, location, expected_obs, istatus) 
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! The next few sections do nothing because there is no additional
! data to read, write, or prompt for.  But there still needs to be a
! case statement in the large select, so they must be here.

! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(IASI_CO)
!     continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(IASI_CO)
!     continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(IASI_CO)
!     continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! This is the code that implements the forward operator.
! Define a module, and make public anything that will be called
! from the main obs_def_mod.f90 file.  Here it is just the
! get_expected routine.  There isn't any initialization needed
! but the stub is there; it could read a namelist if there are
! any run-time options to be set.

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_chem_units_mod

use        types_mod, only : r8, missing_r8
use    utilities_mod, only : register_module
use     location_mod, only : location_type
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : QTY_CO
use  ensemble_manager_mod, only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: get_expected_windspeed

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

real(r8) :: convert_table(max_defined_quantities) = 1.0_r8

contains

! ---------------------------------------------------

subroutine initialize_module

! Handle any module initialization tasks

if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

! ---------------------------------------------------

subroutine get_expected_chem_units(state_handle, ens_size, location, chem, istatus)  

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: chem(ens_size)
integer,             intent(out) :: istatus(ens_size)

! Forward operator for converting chemical species units.
! must match the call in the GET_EXPECTED_OBS_FROM_DEF section above.

real(r8) :: uwind(ens_size)   ! zonal wind component
real(r8) :: vwind(ens_size)   ! meridional wind component
integer  :: this_istatus(ens_size)
logical  :: return_now

if ( .not. module_initialized ) call initialize_module

istatus = 0   ! to use track_status, it must start out 0

if (first_time) then
   first_time = .false.
   ! read in table
   convert_table(:) = 1.0_r8
endif

! Zonal wind at this location - this calls the model_mod code.
call interpolate(state_handle, ens_size, location, QTY_CO, chem, this_istatus)
call track_status(ens_size, this_istatus, wspd, istatus, return_now)
if (return_now) return

! The actual forward operator computation.  This is the value that
! will be returned.  istatus (the return code) of 0 is good,
! any value > 0 indicates an error.  (values < 0 reserved for
! system use.)

where (istatus == 0) chem = convert_table(QTY_CO) * chem
    
end subroutine get_expected_chem_units

! ---------------------------------------------------

end module obs_def_chem_units_mod
! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!-----------------------------------------------------------------------------
!   DART Code:  Johnny Hendricks , hendric at ucar.edu
!   Original DART/Radar work: Nancy Collins  
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS TYPE DEFINITIONS
! DENSITY_ION_E,               QTY_DENSITY_ION_E
! MOLEC_OXYGEN_MIXING_RATIO,   QTY_MOLEC_OXYGEN_MIXING_RATIO
! ATOMIC_OXYGEN_MIXING_RATIO,  QTY_ATOMIC_OXYGEN_MIXING_RATIO
! DENSITY_ION_OP,              QTY_DENSITY_ION_OP
! ION_O_MIXING_RATIO,          QTY_ION_O_MIXING_RATIO
! ATOMIC_H_MIXING_RATIO,       QTY_ATOMIC_H_MIXING_RATIO
! END DART PREPROCESS TYPE DEFINITIONS
!-----------------------------------------------------------------------------
                             
!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!    use obs_def_ion_density_mod, only : get_expected_oxygen_ion_val
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(DENSITY_ION_OP)
!     call get_expected_oxygen_ion_val(state_handle, ens_size, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(DENSITY_ION_E, &
!        MOLEC_OXYGEN_MIXING_RATIO, &
!        ATOMIC_OXYGEN_MIXING_RATIO, &
!        DENSITY_ION_OP, &
!        ION_O_MIXING_RATIO, &
!        ATOMIC_H_MIXING_RATIO)
!      continue
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(DENSITY_ION_E, &
!        MOLEC_OXYGEN_MIXING_RATIO, &
!        ATOMIC_OXYGEN_MIXING_RATIO, &
!        DENSITY_ION_OP, &
!        ION_O_MIXING_RATIO, &
!        ATOMIC_H_MIXING_RATIO)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(DENSITY_ION_E, &
!        MOLEC_OXYGEN_MIXING_RATIO, &
!        ATOMIC_OXYGEN_MIXING_RATIO, &
!        DENSITY_ION_OP, &
!        ION_O_MIXING_RATIO, &
!        ATOMIC_H_MIXING_RATIO)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_ion_density_mod

use        types_mod, only : r8, missing_r8, PI, deg2rad
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             check_namelist_read, find_namelist_in_file,   &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, write_location, read_location, &
                             interactive_location, get_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : QTY_DENSITY_ION_E,              & ! Right QTY to use?
                             QTY_MOLEC_OXYGEN_MIXING_RATIO,  &
                             QTY_ATOMIC_OXYGEN_MIXING_RATIO, &
                             QTY_DENSITY_ION_OP,             &                  
                             QTY_ION_O_MIXING_RATIO,         & ! newly defined
                             QTY_ATOMIC_H_MIXING_RATIO,      & ! newly defined
                             QTY_TEMPERATURE,                &
                             QTY_PRESSURE

use ensemble_manager_mod,  only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

!>@todo compare get_expected_oxygen_ion_val to obs_def_upper_atm_mod.f90:get_expected_oxygen_ion_density ... identical

public :: get_expected_oxygen_ion_val, oxygen_ion_density

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_ion_density_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical :: module_initialized = .false.

!--------------------------------------------------------------
! WACCM-X; put into common/types_mod.f90?
real(r8), PARAMETER :: kboltz = 1.380648E-23_r8    ! [N*m/K]
real(r8), PARAMETER :: universal_gas_constant = 8314.0_r8 ! [J/K/kmol]

!--------------------------------------------------------------         
! WACCM-X; flag to check for the definition of QTYs needed by oxygen_ion_density
logical                 :: first_oxygen_ion_call = .true.

!--------------------------------------------------------------
! Namelist with default values
! 
! use_variable_mean_mass:
 
logical  :: use_variable_mean_mass = .false.

namelist /obs_def_ion_density_mod_nml/ use_variable_mean_mass


contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of executable routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine initialize_module

! Called once to set values and allocate space

integer :: iunit, io, rc

! Prevent multiple calls from executing this code more than once.
if (module_initialized) return

module_initialized = .true.

! Log the version of this source file.
call register_module(source, revision, revdate)

! Read the namelist entry.
call find_namelist_in_file("input.nml", "obs_def_ion_density_mod_nml", iunit)
read(iunit, nml = obs_def_ion_density_mod_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_ion_density_mod_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=obs_def_ion_density_mod_nml)
if (do_nml_term()) write(     *     , nml=obs_def_ion_density_mod_nml)

end subroutine initialize_module

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ION density section
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine get_expected_oxygen_ion_val(state_handle, ens_size, location, obs_val, istatus)

!-----------------------------------------------------------------------------
! This function was implemented for WACCM-X. 
! Check the units for use with other models.
! Given DART state vector and a location, it computes O+ density [1/cm^3].
! The istatus variable should be returned as 0 unless there is a problem.
!

type(ensemble_type),    intent(in) :: state_handle
integer,                intent(in) :: ens_size
type(location_type),    intent(in) :: location
real(r8),              intent(out) :: obs_val(ens_size)
integer,               intent(out) :: istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

istatus = 0 ! Need to initialize this to zero for track_status.
obs_val = MISSING_R8

call oxygen_ion_density(state_handle, ens_size, location, obs_val, istatus)

end subroutine get_expected_oxygen_ion_val

!----------------------------------------------------------------------

subroutine oxygen_ion_density(state_handle, ens_size, location, ion_val, istatus)

type(ensemble_type),    intent(in) :: state_handle
integer,                intent(in) :: ens_size
type(location_type),    intent(in) :: location
real(r8),              intent(out) :: ion_val(ens_size)
integer,               intent(out) :: istatus(ens_size)

logical  :: debug = .false.   ! set to .true. to enable debug printout
integer,  dimension(ens_size) :: mmr_o1_status, mmr_o2_status, mmr_n2_status, ion_op_status
integer,  dimension(ens_size) :: mmr_h1_status, mmr_op_status, p_status, t_status
real(r8), dimension(ens_size) :: mmr_o1, mmr_o2, mmr_n2, mmr_h1, mmr_op   ! mass mixing ratio 
real(r8), dimension(ens_size) :: ion_op, mbar, pressure, temperature 
real(r8), dimension(3)        :: debug_location
real(r8) :: N2_molar_mass, O_molar_mass, O2_molar_mass, H_molar_mass

logical  :: return_now

O_molar_mass  = 28.0_r8 
O2_molar_mass = 16.0_r8 
H_molar_mass  = 32.0_r8 
N2_molar_mass = 1.0_r8 

! Some models have density as part of the state.
! If it is available, just return it. If density is not state,
! then we need to create it from its constituents.
call interpolate(state_handle, ens_size, location, QTY_DENSITY_ION_OP, ion_op, ion_op_status)
call track_status(ens_size, ion_op_status, ion_val, istatus, return_now)
if (all(istatus(:) == 0 )) return

call interpolate(state_handle, ens_size, location, QTY_ATOMIC_OXYGEN_MIXING_RATIO, mmr_o1, mmr_o1_status)
call track_status(ens_size, mmr_o1_status, ion_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_MOLEC_OXYGEN_MIXING_RATIO, mmr_o2, mmr_o2_status)
call track_status(ens_size, mmr_o2_status, ion_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_ATOMIC_H_MIXING_RATIO, mmr_h1, mmr_h1_status)
call track_status(ens_size, mmr_h1_status, ion_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_ION_O_MIXING_RATIO, mmr_op, mmr_op_status)
call track_status(ens_size, mmr_op_status, ion_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_PRESSURE, pressure, p_status)
call track_status(ens_size, p_status, ion_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, QTY_TEMPERATURE, temperature, t_status)
call track_status(ens_size, t_status, ion_val, istatus, return_now)
if (return_now) return

!---------------------------------------------------------------------------------------------------
!  Need to get number density (cgs units) from mass mixing ratio (kg/kg).  
!  mbar is g/mole, same as rMass units
!       kg/kg * (g/mole)/(g/mole) * (Pa = N/m^2)/((Joules/K = N*m/K) * (K)) = m-3 * 1E-06 = cm-3
!---------------------------------------------------------------------------------------------------
! WACCM-X .i file pressure unit is Pa 

where(istatus == 0) 
mmr_n2 = 1.0_r8 - (mmr_o1 + mmr_o2 + mmr_h1)

mbar = 1.0_r8/( mmr_o1/O_molar_mass   &
              + mmr_o2/O2_molar_mass  &
              + mmr_h1/H_molar_mass   &
              + mmr_n2/N2_molar_mass)

ion_val = mmr_op * mbar/O_molar_mass * pressure/(kboltz * temperature) * 1.E-06_r8
endwhere

return

if (debug) then
   debug_location = get_location(location)
   print *, 'final ion_val: ', ion_val
   print *, 'istatus: ', istatus
endif

end subroutine oxygen_ion_density

!----------------------------------------------------------------------

end module obs_def_ion_density_mod
! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------


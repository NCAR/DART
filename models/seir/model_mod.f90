! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

use types_mod,             only : r8, i8, i4, MISSING_R8

use time_manager_mod,      only : time_type, set_time

use location_mod,          only : location_type, set_location, get_location,  &
                                  get_close_obs, get_close_state,             &
                                  convert_vertical_obs, convert_vertical_state

use utilities_mod,         only : register_module, do_nml_file, do_nml_term,    &
                                  nmlfileunit, find_namelist_in_file,           &
                                  check_namelist_read

use location_io_mod,      only :  nc_write_location_atts, nc_write_location

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_end_define_mode

use         obs_kind_mod,  only : QTY_STATE_VARIABLE

use  mpi_utilities_mod,    only : my_task_id

use random_seq_mod,        only : random_seq_type, init_random_seq, random_gaussian

use ensemble_manager_mod,  only : ensemble_type, get_my_num_vars, get_my_vars

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain

use default_model_mod,     only : end_model, nc_write_model_vars

use dart_time_io_mod,      only : read_model_time, write_model_time

implicit none
private

public :: get_model_size,       &
          get_state_meta_data,  &
          model_interpolate,    &
          shortest_time_between_assimilations, &
          static_init_model,    &
          init_conditions,      &
          init_time,            &
          adv_1step,            &
          nc_write_model_atts

public :: pert_model_copies,      &
          nc_write_model_vars,    &
          get_close_obs,          &
          get_close_state,        &
          end_model,              &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time, &
          write_model_time

character(len=*), parameter :: source   = "seir/model_mod.f90"

type(location_type) :: state_loc ! state location, compute once and store for speed
type(random_seq_type)            :: random_seq
type(time_type)                  :: time_step

! Input parameters
integer(i8) :: model_size        = 7
integer     :: time_step_days    = 0
integer     :: time_step_seconds = 3600
real(r8)    :: t_incub           = 5.6_r8      ! Incubation period (days)
real(r8)    :: t_infec           = 3.8_r8      ! Infection time (days)
real(r8)    :: t_recov           = 14.0_r8     ! Recovery period (days)
real(r8)    :: t_death           = 7.0_r8      ! Time until death (days)
real(r8)    :: alpha             = 0.007_r8    ! Vaccination rate (per day)
integer(i8) :: theta             = 12467       ! New birth and new residents (persons per day)       
real(r8)    :: mu                = 0.000025_r8 ! Natural death rate (persons per day)
real(r8)    :: sigma             = 0.05_r8     ! Vaccination inefficacy (e.g., 95% efficient)
real(r8)    :: beta              = 1.36e-9_r8  ! Transmission rate (per day)
real(r8)    :: kappa             = 0.00308     ! Mortality rate  
real(r8)    :: delta_t           = 0.04167_r8  ! Model time step; 1/24 := 1 hour
integer(i8) :: num_pop           = 331996199   ! Population (US)    
real(r8)    :: pert_size         = 1.0         ! Size of perturbation (lognormal pdf param)

real(r8)    :: gama, delta, lambda, rho 

! Other related model parameters
real(r8), parameter :: E0 = 1.0_r8   ! Exposed (not yet infected)
real(r8), parameter :: I0 = 1.0_r8   ! Infected (not yet quarantined)
real(r8), parameter :: Q0 = 1.0_r8   ! Quarantined (confirmed and infected)
real(r8), parameter :: R0 = 1.0_r8   ! Recovered
real(r8), parameter :: D0 = 1.0_r8   ! Dead
real(r8), parameter :: V0 = 1.0_r8   ! Vaccinated
real(r8)            :: S0            ! Susceptible

namelist /model_nml/ model_size, time_step_days, time_step_seconds, &
                     delta_t, num_pop, pert_size, &
                     t_incub, t_infec, t_recov, t_death, &
                     alpha, theta, beta, sigma, kappa, mu   

contains


!------------------------------------------------------------------
!

subroutine static_init_model()

real(r8) :: x_loc
integer  :: i, dom_id

! Do any initial setup needed, including reading the namelist values
call initialize()

! Define the locations of the model state variables
! The SEIR variables have no physical location, 
! and so I'm placing all 7 variables at the same 
! virtual point in space.
x_loc = 0.5_r8  
state_loc =  set_location(x_loc)

! This time is both the minimum time you can ask the model to advance
! and it sets the assimilation window. 
! All observations within +/- 1/2 this interval from the current
! model time will be assimilated.
time_step = set_time(time_step_seconds, time_step_days)

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
dom_id = add_domain(model_size)

end subroutine static_init_model

!------------------------------------------------------------------
! Advance the SEIR model using a four-step RK scheme

subroutine adv_1step(x, time)

real(r8), intent(inout)      :: x(:)
type(time_type), intent(in)  :: time

real(r8), dimension(size(x)) :: xi, x1, x2, x3, x4, dx 

! 1st step
call seir_eqns(x, dx)   
x1 = delta_t * dx

! 2nd step
xi = x + 0.5_r8 * delta_t * dx
call seir_eqns(xi, dx)  
x2 = delta_t * dx

! 3rd step
xi = x + 0.5_r8 * delta_t * dx
call seir_eqns(xi, dx)    
x3 = delta_t * dx

! 4th step
xi = x + delta_t * dx
call seir_eqns(xi, dx)  
x4 = delta_t * dx

! Compute new value for x
x = x + x1/6.0_r8 + x2/3.0_r8 + x3/3.0_r8 + x4/6.0_r8

end subroutine adv_1step

!------------------------------------------------------------------
! SEIR Model Equations
! The following extended SEIR Model with Vaccination 
! is adaopted from Ghostine et al. (2021):
! Ghostine, R., Gharamti, M., Hassrouny, S and Hoteit, I 
! "An Extended SEIR Model with Vaccination for Forecasting
! the COVID-19 Pandemic in Saudi Arabia Using an Ensemble
! Kalman Filter" Mathematics 2021, 9, 636.
! https://dx.doi.org/10.3390/math9060636 

subroutine seir_eqns(x, fx)

! State: x = [S, E, I, Q, R, D, V]

real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: fx(:)

fx(1) = theta - beta * x(1) * x(3) - alpha * x(1) - mu * x(1)
fx(2) = beta * x(1) * x(3) - gama * x(2) + sigma * beta * x(7) * x(3) - mu * x(2)
fx(3) = gama * x(2) - delta * x(3) - mu * x(3)
fx(4) = delta * x(3) - (1.0_r8 - kappa) * lambda * x(4) - kappa * rho * x(4) - mu * x(4)
fx(5) = (1.0_r8 - kappa) * lambda * x(4) - mu * x(5)
fx(6) = kappa * rho * x(4)
fx(7) = alpha * x(1) - sigma * beta * x(7) * x(3) - mu * x(7)

end subroutine seir_eqns

!------------------------------------------------------------------
! Perturbs a model state for generating initial ensembles.
! Returning interf_provided .true. means this code has
! added uniform small independent perturbations to a
! single ensemble member to generate the full ensemble.
subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,   intent(in) :: ens_size
real(r8),  intent(in) :: pert_amp
logical,  intent(out) :: interf_provided

integer                  :: i,j, num_my_grid_points
real(r8)                 :: rng

interf_provided = .true.

call init_random_seq(random_seq, my_task_id()+1)
! if we are running with more than 1 task, then
! we have all the ensemble members for a subset of
! the model state.  which variables we have are determined
! by looking at the global index number into the state vector.

num_my_grid_points = get_my_num_vars(state_ens_handle)

do i=1,num_my_grid_points
    
    ! Lognormal Distribution
    do j= 1, ens_size
         rng = pert_size * random_gaussian(random_seq, 0.0_r8, 1.0_r8)
         state_ens_handle%copies(j, i) = &
            state_ens_handle%copies(j, i) * exp(rng)
    end do
end do


end subroutine pert_model_copies

!------------------------------------------------------------------
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

subroutine init_conditions(x)

real(r8), dimension (model_size) :: x0
real(r8), intent(out)            :: x(:)

x0 = (/S0, E0, I0, Q0, R0, D0, V0/)
x  = x0

end subroutine init_conditions

!------------------------------------------------------------------
! Companion interface to init_conditions. Returns a time that is somehow 
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

subroutine init_time(time)

type(time_type), intent(out) :: time

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time


!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 
! This interface is required for all applications.

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size



!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.
! This interface is required for all applications.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations



!------------------------------------------------------------------
! Given a state handle, a location, and a model state variable quantity,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is an integer that specifies the quantity of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

subroutine model_interpolate(state_handle, ens_size, location, iqty, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: iqty
real(r8),            intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,             intent(out) :: istatus(ens_size)

! Given the nature of the SEIR model, no interpolation 
! is needed. Only identity observations are utilized.

! This should be the result of the interpolation of a
! given quantity (iqty) of variable at the given location.
expected_obs(:) = MISSING_R8

! The return code for successful return should be 0. 
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
istatus(:) = 1

end subroutine model_interpolate



!------------------------------------------------------------------
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument quantity
! (qty) can be returned if the model has more than one type of field
! (for instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

subroutine get_state_meta_data(index_in, location, qty)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: qty

! these should be set to the actual location and state quantity
location = state_loc !(index_in)
if (present(qty)) qty = QTY_STATE_VARIABLE 

end subroutine get_state_meta_data



!------------------------------------------------------------------
! Do any initialization/setup, including reading the
! namelist values.

subroutine initialize()

integer :: iunit, io

! Print module information
call register_module(source)

! Read the namelist 
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Output the namelist values if requested
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Compute other model parameters
gama   = 1.0_r8 / t_incub 
delta  = 1.0_r8 / t_infec
lambda = 1.0_r8 / t_recov
rho    = 1.0_r8 / t_death

! Compute initial value for S:
S0 = num_pop - E0 - I0 - Q0 - R0 - D0

end subroutine initialize


!------------------------------------------------------------------
! Writes model-specific attributes to a netCDF file

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid
integer, intent(in) :: domain_id
integer(i4)         :: offset

! put file into define mode.

integer :: msize

msize = int(model_size, i4)

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source)

call nc_add_global_attribute(ncid, "model", "seir")

call nc_write_location_atts(ncid, msize)
call nc_end_define_mode(ncid)

! Note that location for this model isn't used 
do offset = 1, msize
   call nc_write_location(ncid, state_loc, offset)
enddo

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!===================================================================
! End of model_mod
!===================================================================
end module model_mod


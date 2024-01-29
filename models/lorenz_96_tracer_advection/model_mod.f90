! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

use types_mod,             only : r8, i8, i4

use time_manager_mod,      only : time_type, set_time

use location_mod,          only : location_type, set_location, get_location,  &
                                  get_close_obs, get_close_state,             &
                                  convert_vertical_obs, convert_vertical_state

use utilities_mod,         only : register_module, do_nml_file, do_nml_term,    &
                                  nmlfileunit, find_namelist_in_file,           &
                                  check_namelist_read, E_ERR, error_handler

use location_io_mod,       only :  nc_write_location_atts, nc_write_location

use netcdf_utilities_mod,  only : nc_add_global_attribute, nc_synchronize_file,      &
                                  nc_add_global_creation_time, nc_begin_define_mode, &
                                  nc_end_define_mode

use obs_kind_mod,          only : QTY_STATE_VARIABLE, QTY_TRACER_SOURCE, &
                                  QTY_TRACER_CONCENTRATION, get_name_for_quantity

use  mpi_utilities_mod,    only : my_task_id

use random_seq_mod,        only : random_seq_type, init_random_seq, random_gaussian

use ensemble_manager_mod,  only : ensemble_type, get_my_num_vars, get_my_vars

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain, add_dimension_to_variable, &
                                  finished_adding_domain, state_structure_info

use default_model_mod,     only : end_model, nc_write_model_vars, &
                                  init_time

use dart_time_io_mod,      only : read_model_time, write_model_time

implicit none
private

! these routines must be public and you cannot change the
! arguments because they will be called *from* other DART code.

!> required routines with code in this module
public :: get_model_size,                      &
          get_state_meta_data,                 &
          model_interpolate,                   &
          shortest_time_between_assimilations, &
          static_init_model,                   &
          init_conditions,                     &
          adv_1step,                           &
          nc_write_model_atts

!> required routines where code is in other modules
public :: pert_model_copies,      &
          nc_write_model_vars,    &
          init_time,              &
          get_close_obs,          &
          get_close_state,        &
          end_model,              &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          write_model_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = "model_mod.f90"
character(len=32 ), parameter :: revision = ""
character(len=128), parameter :: revdate  = ""


! Namelist with default values

integer(i8) :: model_size               = 120
real(r8)    :: forcing                  = 8.00_r8
real(r8)    :: delta_t                  = 0.05_r8

! Tracer model parameters with default values
! mean velocity
real(r8)    :: mean_velocity            = 0.00_r8
! velocity normalization
real(r8)    :: pert_velocity_multiplier = 5.00_r8
! diffusion everywhere
real(r8)    :: diffusion_coef           = 0.00_r8
! include an exponential sink rate
real(r8)    :: e_folding                = 0.25_r8
! Also include a fixed sink so tracer can get to 0
real(r8)    :: sink_rate                = 0.1_r8
! Tracer source model parameters
! Amount injected per unit time; This is not currently implemented
real(r8)    :: source_rate              = 100.00_r8
real(r8)    :: point_tracer_source_rate = 5.0_r8

! Allows having negative tracer values to test bounded above filter algorithms
logical     :: positive_tracer          = .true.

! Allows testing non-zero bounds above 
logical     :: bound_above_is_one       = .false.

integer     :: time_step_days           = 0
integer     :: time_step_seconds        = 3600

namelist /model_nml/ model_size, forcing, delta_t, mean_velocity,         &
                     pert_velocity_multiplier, diffusion_coef, e_folding, &
                     sink_rate, source_rate, point_tracer_source_rate,    &
                     positive_tracer, bound_above_is_one,                 &
                     time_step_days, time_step_seconds

! number state variable quantities
integer, parameter  :: NVARS = 3 ! QTY_STATE_VARIABLE, QTY_TRACER_CONCENTRATION, QTY_TRACER_SOURCE

! module global variables
integer :: grid_size  
integer :: var_offset  
integer :: conc_offset 
integer :: source_offset

! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step

! Module storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

contains


!------------------------------------------------------------------
!> Do single time step advance for lorenz 96 model
!> using four-step rk time step.
!>
!> the 'time' argument is unused here but used in larger models.

subroutine adv_1step(x, time)

real(r8), intent(inout)        :: x(:) ! for a grid_size = 40, model_size = 120: 
                            ! positions (1-40) tracer (41-80) and source (81-120)
type(time_type), intent(inout) :: time

real(r8)    :: velocity, target_loc, frac, ratio
integer(r8) :: low, hi, up, down, i
real(r8), dimension(grid_size) :: x1, x2, x3, x4, x_new, dx, inter, q_diff, q_new, q

! Test for tracer with upper bound of 1; Subtract 1 when entering here, then add it back on
if(bound_above_is_one) x(grid_size + 1:2*grid_size) = x(grid_size + 1:2*grid_size) - 1.0_r8

q = x(grid_size + 1 :2*grid_size)  ! QTY_TRACER_CONCENTRATION
! Doing an upstream semi-lagrangian advection for q for each grid point
do i = 1, grid_size
    ! Get the target point
    velocity = (mean_velocity + x(i))*pert_velocity_multiplier

    ! Bail out if the velocity number of grid points per dt exceeds the whole domain size
    if(abs(velocity * delta_t) > grid_size) then
      call error_handler(E_ERR, 'adv_1step', 'Lagrangian Velocity ridiculously large')
    endif

    target_loc = i - velocity*delta_t
    ! Get the bounding grid point
    low = floor(target_loc)
    hi = low + 1
    frac = target_loc - low

    ! Assume for now that we are not looking upstream for multiple revolutions
    ! consistent with error failure above

    if (low < 1) then
      low = low + grid_size
    else if (low > grid_size) then
      low = low - grid_size
    end if

    if (hi < 1) then
      hi = hi + grid_size
    else if (hi > grid_size) then
      hi = hi - grid_size
    end if

    q_new(i) = (1 - frac)*q(low) + frac*q(hi)
end do

! Diffusion for smoothing and avoiding shocky behavior
do i = 1, grid_size
    down = i - 1;
    if (down < 1) then
      down = down + grid_size
    end if
    up = i + 1;
    if (up > grid_size) then
      up = up - grid_size
    end if
    ! Should be sure this is the right way to time normalize
    q_diff(i) = diffusion_coef * delta_t * (q_new(down) + q_new(up) - 2*q_new(i))
end do

q_new = q_new + q_diff*delta_t

! Add source 
q_new = x((2*grid_size)+1 : model_size)*delta_t + q_new
! Add exponential sinks at every grid point
ratio = exp((-1)*e_folding*delta_t)
q_new = ratio*q_new

! Add in an additional uniform sink so that stuff can get to zero for tests
if(positive_tracer) then
   q_new = max(0.0_r8, q_new - sink_rate * delta_t)
else
   q_new = min(0.0_r8, q_new + sink_rate * delta_t)
endif

x(grid_size+1:2*(grid_size)) = q_new

! RK4 solver for the lorenz-96 equations

! Compute first intermediate step
call comp_dt(x(1: grid_size), dx)
x1 = delta_t * dx
inter = x(1: grid_size)+ x1/2

! Compute second intermediate step
call comp_dt(inter, dx)
x2 = delta_t * dx
inter = x(1: grid_size) + x2/2

! Compute third intermediate step
call comp_dt(inter, dx)
x3 = delta_t * dx
inter = x(1: grid_size) + x3

! Compute fourth intermediate step
call comp_dt(inter, dx)
x4 = delta_t * dx

! Compute new value for x
x_new = x(1: grid_size) + x1/6 + x2/3 + x3/3 + x4/6

x(1: grid_size) = x_new

! Test for tracer with upper bound of 1; Subtract 1 when entering, then add it back on
if(bound_above_is_one) x(grid_size + 1:2*grid_size) = x(grid_size + 1:2*grid_size) + 1.0_r8

end subroutine adv_1step

!------------------------------------------------------------------
!> Computes the time tendency of the lorenz 1996 model given current state

subroutine comp_dt(x, dt)

real(r8), intent(in) :: x(grid_size)
real(r8), intent(out) :: dt(grid_size)

integer :: j, jp1, jm1, jm2

do j = 1, grid_size
      jp1 = j + 1
      if (jp1 > grid_size) jp1 = 1
      jm2 = j - 2
      if (jm2 < 1) jm2 = grid_size + jm2
      jm1 = j - 1
      if (jm1 < 1) jm1 = grid_size

      dt(j) = (x(jp1) - x(jm2)) * x(jm1) - x(j) + forcing
end do

end subroutine comp_dt

!------------------------------------------------------------------
!> This routine is called once and initializes any information needed
!> by other routines in this file.

subroutine static_init_model()

real(r8) :: delta_loc
integer  :: i, dom_id
character(20) :: string1

! Do any initial setup needed, including reading the namelist values
call initialize()

! check that model size is divisible by the number of variables
if ( mod(model_size,NVARS) /= 0 ) then 
   write(string1,'(I5)') NVARS 
   call error_handler(E_ERR, 'static_init_model', 'model_size must be a multiple of '// trim(string1))
endif
     
grid_size = model_size/NVARS
var_offset = 0 
conc_offset = grid_size
source_offset = grid_size*2

! Create storage for locations
allocate(state_loc(model_size))

! Define the locations of the model state variables
delta_loc = 1.0_r8/real(grid_size, r8)

do i = 1, model_size
  state_loc(i) =  set_location(delta_loc*mod(i-1,grid_size))
enddo

! The time_step in terms of a time type must also be initialized.
! (Need to determine appropriate non-dimensionalization conversion
! for L96 from Shree Khare.)
time_step = set_time(time_step_seconds, time_step_days)

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
dom_id = add_domain('template.nc', NVARS, &
                              (/ 'state_variable      ', &
                                 'tracer_concentration', &
                                 'source              ' /),&
                              (/QTY_STATE_VARIABLE, QTY_TRACER_CONCENTRATION, QTY_TRACER_SOURCE/))


end subroutine static_init_model

!------------------------------------------------------------------
!> Supply initial conditions for lorenz 96 if not reading restart
!> data from a file.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)


! Set all variables, winds, tracer concentration, and source to 0
x(:) = 0.0_r8
! Add a single perturbation to L96 state (winds) to generate evolution
x(1) = 0.1_r8
! For these tests, single tracer source at the first grid point
if(positive_tracer) then
   x(grid_size*2 + 1) = point_tracer_source_rate
else
   ! Make it negative if testing negative tracers
   x(grid_size*2 + 1) = -point_tracer_source_rate
endif

end subroutine init_conditions

!------------------------------------------------------------------
! Returns number of items in state vector

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size

!------------------------------------------------------------------
!> Return the minimum amount of time the model can be advanced by.
!> This is unrelated to the internal RK timesteps.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations

!------------------------------------------------------------------
!> Given a location, return the interpolated state value.
!> expected_obs() and istatus() are arrays.

subroutine model_interpolate(state_handle, ens_size, location, itype, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: itype
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer(i8) :: lower_index, upper_index
integer :: offset
real(r8) :: lctn, lctnfrac
real(r8) :: x_lower(ens_size) !< the lower piece of state vector
real(r8) :: x_upper(ens_size) !< the upper piece of state vector

character(128) :: string1

! All forward operators supported
istatus(:) = 0

! Convert location to real
lctn = get_location(location)
! Multiply by grid size assuming domain is [0, 1] cyclic
lctn = grid_size * lctn

lower_index = int(lctn) + 1
upper_index = lower_index + 1
if(lower_index > grid_size) lower_index = lower_index - grid_size
if(upper_index > grid_size) upper_index = upper_index - grid_size

lctnfrac = lctn - int(lctn)

! Now figure out which type of quantity the location indicates
if (itype == QTY_STATE_VARIABLE) then
   offset = var_offset
else if (itype == QTY_TRACER_CONCENTRATION) then
   offset = conc_offset
else if (itype == QTY_TRACER_SOURCE) then
   offset = source_offset
else
   write(string1, *) 'quantity ', itype, ' ('//trim(get_name_for_quantity(itype))// &
                  ') is not supported in model_interpolate'
   call error_handler(E_ERR,'model_interpolate',string1, source, revision, revdate)
end if

! Add the offset to the lower and the upper indices
lower_index = lower_index + offset
upper_index = upper_index + offset

! Grab the correct pieces of state vector

! Lower value
x_lower(:) = get_state(lower_index, state_handle)

! Upper value
x_upper(:) = get_state(upper_index, state_handle)

! calculate the obs value

expected_obs(:) = (1.0_r8 - lctnfrac) * x_lower(:) + lctnfrac * x_upper(:)

end subroutine model_interpolate
!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location and variable type.

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

if (present(var_type)) then
   if (index_in <= grid_size) then
      var_type = QTY_STATE_VARIABLE
      location = state_loc(index_in)
   end if

   if (grid_size < index_in .and. index_in <= (2*grid_size)) then
      var_type = QTY_TRACER_CONCENTRATION
      location = state_loc(index_in)
   end if

   if ( (2*grid_size) < index_in .and. index_in <= model_size) then
      var_type = QTY_TRACER_SOURCE
      location = state_loc(index_in)
   end if
end if

end subroutine get_state_meta_data

!------------------------------------------------------------------
!> Do any initialization/setup, including reading the
!> namelist values.

subroutine initialize()

integer :: iunit, io

! Print module information
call register_module(source, revision, revdate)

! Read the namelist
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Output the namelist values if requested
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=            model_nml)

end subroutine initialize

!------------------------------------------------------------------
!> Perturbs a model state for generating initial ensembles.
!> Returning interf_provided .true. means this code has
!> added uniform small independent perturbations to a
!> single ensemble member to generate the full ensemble.
subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,   intent(in) :: ens_size
real(r8),  intent(in) :: pert_amp
logical,  intent(out) :: interf_provided

integer :: i,j, num_my_grid_points
integer(i8), allocatable :: my_grid_points(:)
type(location_type) :: location
integer :: var_type

interf_provided = .true.

call init_random_seq(random_seq, my_task_id()+1)
! if we are running with more than 1 task, then
! we have all the ensemble members for a subset of
! the model state.  which variables we have are determined
! by looking at the global index number into the state vector.

! how many grid points does my task have to work on?
! and what are their indices into the full state vector?
num_my_grid_points = get_my_num_vars(state_ens_handle)
allocate(my_grid_points(num_my_grid_points))
call get_my_vars(state_ens_handle, my_grid_points)

do i=1,num_my_grid_points
    call get_state_meta_data(my_grid_points(i), location, var_type)

    if (var_type == QTY_STATE_VARIABLE) then
        do j=1,ens_size
            state_ens_handle%copies(j, i) = random_gaussian(random_seq, &
               state_ens_handle%copies(j, i), pert_amp)
        end do
    elseif(var_type == QTY_TRACER_CONCENTRATION) then
        ! For now, can just let all ensemble members be identical
        ! Spread will be generated by the chaotic flow field
    !Perturbing all source grid points
    else if (var_type == QTY_TRACER_SOURCE) then
        ! For now, can just keep source constant so it will not evolve
        ! Need to perturb to do source estimation
    end if
end do

deallocate(my_grid_points)
end subroutine pert_model_copies

!------------------------------------------------------------------
! Writes model-specific attributes to a netCDF file

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid
integer, intent(in) :: domain_id

integer :: msize

msize = int(model_size, i4)

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate", revdate )

call nc_add_global_attribute(ncid, "model", "Lorenz_96_Tracer_Advection")
call nc_add_global_attribute(ncid, "model_forcing", forcing )
call nc_add_global_attribute(ncid, "model_delta_t", delta_t )
call nc_add_global_attribute(ncid, "source_rate", source_rate)
call nc_add_global_attribute(ncid, "sink_rate", sink_rate)
call nc_add_global_attribute(ncid, "exponential_sink_folding", e_folding)
call nc_add_global_attribute(ncid, "mean_velocity", mean_velocity)
call nc_add_global_attribute(ncid, "pert_velocity_multiplier", pert_velocity_multiplier)
call nc_add_global_attribute(ncid, "diffusion_coef", diffusion_coef)


call nc_write_location_atts(ncid, grid_size)
call nc_end_define_mode(ncid)
call nc_write_location(ncid, state_loc(1:grid_size), grid_size)

call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!===================================================================
! End of model_mod
!===================================================================
end module model_mod


! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

use types_mod,             only : r8, i8, i4

use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use time_manager_mod,      only : time_type, set_time

use     location_mod,      only : location_type, set_location, get_location, &
                                  get_close_obs, get_close_state, &
                                  convert_vertical_obs, convert_vertical_state

use utilities_mod,         only : register_module, do_nml_file, do_nml_term,    &
                                  nmlfileunit, find_namelist_in_file, E_ERR,    &
                                  check_namelist_read, error_handler

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_end_define_mode

use location_io_mod,      only :  nc_write_location_atts, nc_get_location_varids, &
                                  nc_write_location

use default_model_mod,     only : end_model, pert_model_copies, nc_write_model_vars, &
                                  init_time

use         obs_kind_mod,  only : QTY_STATE_VARIABLE, QTY_1D_PARAMETER

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain

use dart_time_io_mod,      only : read_model_time, write_model_time

implicit none
private

! these routines must be public and you cannot change the
! arguments because they will be called *from* other DART code.

!> required routines with code in this module
public :: get_model_size, &
          get_state_meta_data,  &
          model_interpolate, &
          shortest_time_between_assimilations, &
          static_init_model, &
          init_conditions,    &
          init_time, &
          adv_1step, &
          nc_write_model_atts

!> required routines where code is in other modules
public :: pert_model_copies, &
          nc_write_model_vars, &
          get_close_obs, &
          get_close_state, &
          end_model, &
          convert_vertical_obs, &
          convert_vertical_state, &
          read_model_time, &
          write_model_time


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Namelist with default values

integer  :: num_state_vars = 40
real(r8) :: forcing    = 8.00_r8
real(r8) :: delta_t    = 0.05_r8
integer  :: time_step_days = 0
integer  :: time_step_seconds = 3600
logical  :: reset_forcing = .false.
real(r8) :: random_forcing_amplitude = 0.0_r8

namelist /model_nml/ num_state_vars, forcing, delta_t, time_step_days, &
   time_step_seconds, reset_forcing, random_forcing_amplitude

! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step
integer(i8) :: model_size

! Adding random noise
type(random_seq_type) :: random


contains


!------------------------------------------------------------------
!> initialize model.  this routine is called once.

subroutine static_init_model()

real(r8) :: x_loc
integer  :: i, iunit, io, dom_id

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Model size is twice the number of state_vars
model_size = 2 * num_state_vars

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Create storage for locations
allocate(state_loc(model_size))

! Define the locations of the model state variables
do i = 1, num_state_vars
   x_loc = (i - 1.0_r8) / num_state_vars
   state_loc(i) =  set_location(x_loc)
   ! Forcing is at same location as corresponding state variable
   state_loc(i + num_state_vars) = set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.
time_step = set_time(time_step_seconds, time_step_days)

! Initialize the random sequence
call init_random_seq(random)

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
dom_id = add_domain(model_size)

end subroutine static_init_model


!------------------------------------------------------------------
!> Computes the time tendency of the lorenz 1996 model given current state
 
subroutine comp_dt(x, dt)

real(r8), intent( in) ::  x(:)
real(r8), intent(out) :: dt(:)

integer :: j, jp1, jm1, jm2

do j = 1, num_state_vars
   jp1 = j + 1
   if(jp1 > num_state_vars) jp1 = 1
   jm2 = j - 2
   if(jm2 < 1) jm2 = num_state_vars + jm2
   jm1 = j - 1
   if(jm1 < 1) jm1 = num_state_vars
   
   dt(j) = (x(jp1) - x(jm2)) * x(jm1) - x(j) + x(num_state_vars + j)
end do

! Time tendency for the forcing variables
dt(num_state_vars + 1 : model_size) = 0.0_r8


! Try adding in some random spread; fixed across the instances (basically a global var)
if(.not. reset_forcing) &
   dt(num_state_vars + 1 : model_size) = &
      random_gaussian(random, 0.0_r8, random_forcing_amplitude)

! Try adding in some random noise to each forcing variable, completely local
!if(.not. reset_forcing) then
!   do j = num_state_vars + 1, model_size
!      dt(j) = random_gaussian(random, 0.0_r8, random_forcing_amplitude)
!   end do
!endif 
   
end subroutine comp_dt


!------------------------------------------------------------------
!> Does single time step advance for forced lorenz 96 model
!> using four-step rk time step

subroutine adv_1step(x, time)

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

real(r8), dimension(size(x)) :: x1, x2, x3, x4, dx, inter

! If reset forcing is set, grab the value from the namelist and hold it fixed
if(reset_forcing) x(num_state_vars + 1: model_size) = forcing

call comp_dt(x, dx)        !  Compute the first intermediate step
x1    = delta_t * dx
inter = x + x1 / 2.0_r8

call comp_dt(inter, dx)    !  Compute the second intermediate step
x2    = delta_t * dx
inter = x + x2 / 2.0_r8

call comp_dt(inter, dx)    !  Compute the third intermediate step
x3    = delta_t * dx
inter = x + x3

call comp_dt(inter, dx)    !  Compute fourth intermediate step
x4 = delta_t * dx

!  Compute new value for x

x = x + x1/6.0_r8 + x2/3.0_r8 + x3/3.0_r8 + x4/6.0_r8

end subroutine adv_1step


!------------------------------------------------------------------
!> Returns number of items in the state vector

function get_model_size()
integer(i8) :: get_model_size

! Number of state variables is twice the model size if forcing is assumed
! this was computed in the static_init_model() routine
get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------
!> Initial conditions for forced lorenz 96

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

! Set state variable to the value of forcing with 1st element slightly perturbed
! Set forcing parameters to 8.0 if being assimilated
x    = forcing
x(1) = 1.001_r8 * forcing

end subroutine init_conditions


!------------------------------------------------------------------
!> Interpolates from state vector x to the location.

subroutine model_interpolate(state_handle, ens_size, location, itype, expected_val, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: expected_val(ens_size)
integer,            intent(out) :: istatus(ens_size)

integer(i8) :: lower_index, upper_index
integer :: i
real(r8) :: lctn, lctnfrac

! All forward operators supported
istatus(:) = 0

! Convert location to real
lctn = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
lctn = num_state_vars * lctn

lower_index = int(lctn) + 1
upper_index = lower_index + 1
if(lower_index > num_state_vars) lower_index = lower_index - num_state_vars
if(upper_index > num_state_vars) upper_index = upper_index - num_state_vars

lctnfrac = lctn - int(lctn)
expected_val(:) = (1.0_r8 - lctnfrac) * get_state(lower_index, state_handle) + &
                            lctnfrac  * get_state(upper_index, state_handle)

if(1 == 1) return


! All the stuff below is for strange forward operator tests; not currently used
!-----------------------------------------------
!!!expected_val = expected_val ** 2
!!!if(1 == 1) return

! Temporarily add on an observation from the other side of the domain, too
lower_index = lower_index + model_size / 2
if(lower_index > model_size) lower_index = lower_index - model_size
upper_index = upper_index + model_size / 2
if(upper_index > model_size) upper_index = upper_index - model_size
expected_val(:) = expected_val(:) + &
             lctnfrac  * get_state(lower_index, state_handle) + &
   (1.0_r8 - lctnfrac) * get_state(upper_index, state_handle)
if(1 == 1) return


! Next one does an average over a range of points
expected_val(:) = 0.0_r8
lower_index = lower_index - 7
upper_index = upper_index - 7
if(lower_index < 1) lower_index = lower_index + model_size
if(upper_index < 1) upper_index = upper_index + model_size

do i = 1, 15
   if(lower_index > model_size) lower_index = lower_index - model_size
   if(upper_index > model_size) upper_index = upper_index - model_size
   expected_val(:) = expected_val(:) + &
      (1.0_r8 - lctnfrac) * get_state(lower_index, state_handle) + &
                lctnfrac  * get_state(upper_index, state_handle)
   lower_index = lower_index + 1
   upper_index = upper_index + 1
end do

end subroutine model_interpolate


!------------------------------------------------------------------
!> Return the minimum advance time of the model

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations


!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location and optionally type.

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

location = state_loc(index_in)
if (present(var_type)) then
   if(index_in <= num_state_vars) then
      var_type = QTY_STATE_VARIABLE    ! default variable type
   else
      ! If forcing parameter is being assimilated, it has this var_type
      var_type = QTY_1D_PARAMETER
   endif 
endif

end subroutine get_state_meta_data


!------------------------------------------------------------------
!> Writes the model-specific attributes to a netCDF file

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in)  :: ncid
integer, intent(in) :: domain_id

integer :: msize

! other parts of the dart system will write the state into the file
! so this routine just needs to write any model-specific
! attributes it wants to record.

msize = int(model_size, i4)

! Write Global Attributes

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate", revdate )

call nc_add_global_attribute(ncid, "model", "Forced_Lorenz_96")
call nc_add_global_attribute(ncid, "model_forcing", forcing )
call nc_add_global_attribute(ncid, "model_delta_t", delta_t )
call nc_add_global_attribute(ncid, "model_num_state_vars", num_state_vars)
call nc_add_global_attribute(ncid, "model_time_step_days", time_step_days)
call nc_add_global_attribute(ncid, "model_time_step_seconds", time_step_seconds)
call nc_add_global_attribute(ncid, "model_random_forcing_amplitude", random_forcing_amplitude)
if (reset_forcing) then
   call nc_add_global_attribute(ncid, "model_reset_forcing", "TRUE")
else
   call nc_add_global_attribute(ncid, "model_reset_forcing", "FALSE")
endif

call nc_write_location_atts(ncid, msize)
call nc_end_define_mode(ncid)
call nc_write_location(ncid, state_loc, msize)

call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

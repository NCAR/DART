! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

use types_mod,             only : r8, i8, i4

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

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain

use default_model_mod,     only : end_model, pert_model_copies, nc_write_model_vars, &
                                  init_time

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
          adv_1step, &
          nc_write_model_atts

!> required routines where code is in other modules
public :: pert_model_copies, &
          nc_write_model_vars, &
          init_time, &
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

integer(i8) :: model_size = 40
real(r8)    :: forcing    = 8.00_r8
real(r8)    :: delta_t    = 0.05_r8
integer     :: time_step_days = 0
integer     :: time_step_seconds = 3600

namelist /model_nml/ model_size, forcing, delta_t, time_step_days, time_step_seconds


! module global variables

! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step


contains


!------------------------------------------------------------------
!> Do single time step advance for lorenz 96 model
!> using four-step rk time step.
!>
!> the 'time' argument is unused here but used in larger models.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

real(r8), dimension(size(x)) :: x1, x2, x3, x4, dx, inter

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
!> Computes the time tendency of the lorenz 1996 model given current state

subroutine comp_dt(x, dt)

real(r8), intent(in)  ::  x(:)
real(r8), intent(out) :: dt(:)

integer :: j, jp1, jm1, jm2, ms

! avoid compiler bugs with long integers
! being used as loop indices.
ms = model_size
do j = 1, ms
   jp1 = j + 1
   if(jp1 > ms) jp1 = 1
   jm2 = j - 2
   if(jm2 < 1) jm2 = ms + jm2
   jm1 = j - 1
   if(jm1 < 1) jm1 = ms
   
   dt(j) = (x(jp1) - x(jm2)) * x(jm1) - x(j) + forcing
end do

end subroutine comp_dt

!------------------------------------------------------------------
!> This routine is called once and initializes any information needed
!> by other routines in this file.

subroutine static_init_model()

real(r8) :: x_loc
integer  :: i, dom_id

! Do any initial setup needed, including reading the namelist values
call initialize()

! Create storage for locations
allocate(state_loc(model_size))

! Define the locations of the model state variables
do i = 1, model_size
   x_loc = (i - 1.0_r8) / model_size
   state_loc(i) =  set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized. 
! (Need to determine appropriate non-dimensionalization conversion 
! for L96 from Shree Khare.)
time_step = set_time(time_step_seconds, time_step_days)

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
dom_id = add_domain(model_size)

end subroutine static_init_model

!------------------------------------------------------------------
!> Supply initial conditions for lorenz 96 if not reading restart
!> data from a file.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

x    = forcing
x(1) = 1.001_r8 * forcing

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
real(r8) :: lctn, lctnfrac
real(r8) :: x_lower(ens_size) !< the lower piece of state vector
real(r8) :: x_upper(ens_size) !< the upper piece of state vector

! All forward operators supported
istatus(:) = 0

! Convert location to real
lctn = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
lctn = model_size * lctn

lower_index = int(lctn) + 1
upper_index = lower_index + 1
if(lower_index > model_size) lower_index = lower_index - model_size
if(upper_index > model_size) upper_index = upper_index - model_size

lctnfrac = lctn - int(lctn)

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

location = state_loc(index_in)
if (present(var_type)) var_type = QTY_STATE_VARIABLE 

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
if (do_nml_term()) write(     *     , nml=model_nml)

end subroutine initialize

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

call nc_add_global_attribute(ncid, "model", "Lorenz_96")
call nc_add_global_attribute(ncid, "model_forcing", forcing )
call nc_add_global_attribute(ncid, "model_delta_t", delta_t )

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

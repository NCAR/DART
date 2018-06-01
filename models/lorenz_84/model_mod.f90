! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! Implements Lorenz 84 3 variable model with intermediate level
! attractor.

use        types_mod,      only : r8, i8, i4

use time_manager_mod,      only : time_type, set_time

use    utilities_mod,      only : register_module, error_handler, E_ERR, E_MSG, nmlfileunit, &
                                  do_output, find_namelist_in_file, check_namelist_read,     &
                                  do_nml_file, do_nml_term

use     location_mod,      only : location_type, set_location, get_location, &
                                  get_close_obs, get_close_state, &
                                  convert_vertical_obs, convert_vertical_state

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode

use location_io_mod,      only :  nc_write_location_atts, nc_get_location_varids, &
                                  nc_write_location

use default_model_mod,     only : end_model, pert_model_copies, nc_write_model_vars, &
                                  init_time

use         obs_kind_mod,  only : QTY_STATE_VARIABLE

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain, state_structure_info

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

!  define model parameters

! Model size is fixed
integer(i8), parameter :: model_size = 3

! Namelist with default values

real(r8) ::      a = 0.25_r8
real(r8) ::      b = 4.00_r8
real(r8) ::      f = 8.00_r8 
real(r8) ::      g = 1.25_r8
real(r8) :: deltat = 0.01_r8
integer  :: time_step_days = 0
integer  :: time_step_seconds = 3600

namelist /model_nml/ a, b, f, g, deltat, time_step_days, time_step_seconds

! Define the location of the state variables in module storage
type(location_type) :: state_loc(model_size)
type(time_type)     :: time_step


contains


!------------------------------------------------------------------
!> initializes class data - called once

subroutine static_init_model()

real(r8) :: x_loc
integer(i8) :: i
integer  :: iunit, io, dom_id

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Define the locations of the model state variables
do i = 1, model_size
   x_loc = (i - 1.0_r8) / model_size
   state_loc(i) =  set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized.
time_step = set_time(time_step_seconds, time_step_days)

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
dom_id = add_domain(model_size)

end subroutine static_init_model


!------------------------------------------------------------------
!> compute time tendency of lorenz-84 3-variable model given current state

subroutine comp_dt(x, dt)

real(r8), intent( in)  ::  x(:)
real(r8), intent(out) :: dt(:)

!  compute lorenz 84 model time tendency

dt(1) = -x(2)**2 -x(3)**2 - a*x(1) + a*f
dt(2) = x(1)*x(2) - b*x(1)*x(3) - x(2) + g
dt(3) = b*x(1)*x(2) + x(1)*x(3) - x(3)

end subroutine comp_dt


!------------------------------------------------------------------
!>  defines off-attractor initial conditions for lorenz-84

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

x(1) = 0.0_r8
x(2) = 1.0_r8
x(3) = 0.0_r8

end subroutine init_conditions


!------------------------------------------------------------------
!> does single time step advance 
subroutine adv_1step(x, time)

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

real(r8) :: fract

fract = 1.0_r8
call adv_single(x, fract)

end  subroutine adv_1step


!------------------------------------------------------------------
!>  does single time step advance for lorenz 3 variable model using two step rk

subroutine adv_single(x, fract)

real(r8), intent(inout) :: x(:)
real(r8), intent(in)    :: fract

real(r8) :: x1(3), x2(3), dx(3)

call comp_dt(x, dx)            !  compute the first intermediate step
x1 = x + fract * deltat * dx

call comp_dt(x1, dx)           !  compute the second intermediate step
x2 = x1 + fract * deltat * dx

!  new value for x is average of original value and second intermediate

x = (x + x2) / 2.0_r8

end subroutine adv_single


!------------------------------------------------------------------
!> Returns number of items in the state vector

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------
!> Interpolates an ensemble of expected values at the given location.

subroutine model_interpolate(state_handle, ens_size, location, itype, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: itype
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer(i8) :: lower_index, upper_index
real(r8) :: lctn, lctnfrac

! All interpolations okay for now
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
expected_obs(:) = (1.0_r8 - lctnfrac) * get_state(lower_index, state_handle) + &
                            lctnfrac  * get_state(upper_index, state_handle)

end subroutine model_interpolate


!------------------------------------------------------------------
!> Returns the mininum time step of the model.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations


!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location and optionally the quantity.

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

location = state_loc(index_in)
if (present(var_type)) var_type = QTY_STATE_VARIABLE    ! default variable quantity

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

call nc_add_global_attribute(ncid, "model_source", source)
call nc_add_global_attribute(ncid, "model_revision", revision)
call nc_add_global_attribute(ncid, "model_revdate", revdate)

call nc_add_global_attribute(ncid, "model", "Lorenz_84")
call nc_add_global_attribute(ncid, "model_a", a)
call nc_add_global_attribute(ncid, "model_b", b)
call nc_add_global_attribute(ncid, "model_f", f)
call nc_add_global_attribute(ncid, "model_g", g)
call nc_add_global_attribute(ncid, "model_delta_t", deltat)

call nc_write_location_atts(ncid, msize)
call nc_end_define_mode(ncid)
call nc_write_location(ncid, state_loc, msize)

call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts


!-------------------------------------------------------------------
!>  compute linear operator around state nl
subroutine linearize(nl, l)

real(r8) :: nl(3), l(3, 3)

l(1, 1) = -1.0_r8 *     a     * deltat + 1.0_r8
l(1, 2) = -1.0_r8 * nl(2)     * deltat
l(1, 3) = -1.0_r8 * nl(3)     * deltat
l(2, 1) = (nl(2) - b * nl(3)) * deltat
l(2, 2) = (-1.0_r8 + nl(1))   * deltat + 1.0_r8
l(2, 3) = -1.0_r8 * b * nl(1) * deltat
l(3, 1) = (b * nl(2) + nl(3)) * deltat
l(3, 2) =  b * nl(1)          * deltat
l(3, 3) = (nl(1) - 1.0_r8)    * deltat + 1.0_r8

end subroutine linearize

!--------------------------------------------------------------------

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

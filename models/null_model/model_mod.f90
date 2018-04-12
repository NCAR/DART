! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This module provides very simple models for evaluating filtering algorithms.
! It can provide simple linear growth around a fixed point, a random draw from
! a Gaussian, or combinations of the two. Select the model advance method by
! the namelist item 'advance_method'.

use        types_mod,      only : r8, i8, i4

use    utilities_mod,      only : register_module, error_handler, E_ERR, E_MSG, &
                                  nmlfileunit, do_output, find_namelist_in_file, &
                                  check_namelist_read, do_nml_file, do_nml_term

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_sync, &
                                 nc_add_global_creation_time, nc_redef, nc_enddef

use time_manager_mod, only : time_type, set_time

use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain

use     location_mod,      only : location_type, set_location, get_location, &
                                  get_close_obs, get_close_state, &
                                  convert_vertical_obs, convert_vertical_state

use location_io_mod,      only :  nc_write_location_atts, nc_get_location_varids, &
                                  nc_write_location

use         obs_kind_mod,  only : QTY_STATE_VARIABLE

use ensemble_manager_mod,  only : ensemble_type

use dart_time_io_mod,      only : read_model_time, write_model_time

use default_model_mod,     only : pert_model_copies, nc_write_model_vars, &
                                  init_time, init_conditions

implicit none
private

! these routines must be public and you cannot change the
! arguments because they will be called *from* other DART code.

!> required routines with code in this module
public :: get_model_size, &
          get_state_meta_data, &
          model_interpolate, &
          shortest_time_between_assimilations, &
          static_init_model, &
          adv_1step, &
          end_model, &
          nc_write_model_atts

!> required routines where code is in other modules
public :: pert_model_copies, &
          nc_write_model_vars, &
          init_time, &
          init_conditions, &
          get_close_obs, &
          get_close_state, &
          convert_vertical_obs, &
          convert_vertical_state, &
          read_model_time, &
          write_model_time


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Basic model parameters controlled by nameslist; have defaults

! Namelist with default values
! Model size can be as small as 1 here.
integer  :: model_size        = 2
real(r8) :: delta_t           = 0.05_r8
integer  :: time_step_days    = 0
integer  :: time_step_seconds = 3600
real(r8) :: noise             = 1.0_r8
real(r8) :: noise_mult_factor = 1.0_r8
character(len=64) :: advance_method = 'rk'

! valid options for advance_method are: rk, simple, random
namelist /model_nml/ model_size, delta_t, time_step_days, &
                     time_step_seconds, advance_method, &
                     noise, noise_mult_factor


! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step

! Following is for repeatable random numbers
logical :: first_ens_seq = .true.
type (random_seq_type) :: ens_seq

contains

!------------------------------------------------------------------
!> Called once to initialize information for null model

subroutine static_init_model()

real(r8) :: x_loc
integer  :: i, iunit, io, dom_id

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

allocate(state_loc(model_size))

! Define the locations of the model state variables
do i = 1, model_size
   x_loc = (i - 1.0_r8) / model_size
   state_loc(i) =  set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized. 
time_step = set_time(time_step_seconds, time_step_days)

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
dom_id = add_domain(int(model_size, i8))

end subroutine static_init_model


!------------------------------------------------------------------
!> Computes time tendency of the null model given the current state

subroutine comp_dt(x, dt)

real(r8), intent( in) ::  x(:)
real(r8), intent(out) :: dt(:)

integer :: j

if(first_ens_seq) then
   call init_random_seq(ens_seq)
   first_ens_seq = .false.
end if

do j = 1, model_size
   if (noise <= 0.0_r8) then
      dt(j) =  x(j)
   else
      ! add noise
      dt(j) = noise_mult_factor * random_gaussian(ens_seq, 0.0_r8, noise)
   endif
end do

end subroutine comp_dt



!------------------------------------------------------------------
!> does single time step advance for null model
!> using four step rk time step.  set method by namelist.
!> valid values for 'advance_method' are:  rk, simple, random

subroutine adv_1step(x, time)

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

real(r8), dimension(size(x)) :: x1, x2, x3, x4, dx, inter

if (advance_method == 'rk') then

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
   
else if (advance_method == 'simple') then

   ! simple timestepping
   call comp_dt(x, dx)
   x = x + delta_t * dx

else if (advance_method == 'random') then

   ! IDEALIZED DISTRIBUTION TEST: Just draw from a random distribution
   call comp_dt(x, dx)
   x = dx

else

   call error_handler(E_ERR, 'adv_1step', 'unrecognized advance_method: "'//trim(advance_method)//'"', &
                      source, revision, revdate, text2='valid strings are: rk, simple, random')
endif


end subroutine adv_1step


!------------------------------------------------------------------
!> Interpolates an ensemble of expected values at the given location.
!>
!> Argument itype is not used here because there is only one type of variable.

subroutine model_interpolate(state_handle, ens_size, location, itype, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer(i8)  :: lower_index, upper_index
real(r8) :: lctn, lctnfrac

! All obs okay for now
istatus = 0

! Convert location to real
lctn = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
lctn = model_size * lctn

lower_index = int(lctn) + 1
upper_index = lower_index + 1
if(lower_index > model_size) lower_index = lower_index - model_size
if(upper_index > model_size) upper_index = upper_index - model_size

lctnfrac = lctn - int(lctn)
expected_obs = (1.0_r8 - lctnfrac) * get_state(lower_index, state_handle) + lctnfrac * get_state(upper_index, state_handle)

end subroutine model_interpolate


!------------------------------------------------------------------
!> Returns number of items in the state vector

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size

!------------------------------------------------------------------
!> Returns the mininum time step of the model.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations


!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location.

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

location = state_loc(index_in)
if (present(var_type)) var_type = QTY_STATE_VARIABLE    ! default variable quantity

end subroutine get_state_meta_data


!------------------------------------------------------------------
!> Called once at the end of execution; free allocated storage

subroutine end_model()

deallocate(state_loc)

end subroutine end_model


!------------------------------------------------------------------
!> Writes the model-specific attributes to a netCDF file

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in)  :: ncid
integer, intent(in) :: domain_id

! other parts of the dart system will write the state into the file
! so this routine just needs to write any model-specific
! attributes it wants to record.


! Write Global Attributes 

call nc_redef(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate", revdate )

call nc_add_global_attribute(ncid, "model", "null")
call nc_add_global_attribute(ncid, "model_delta_t", delta_t )

!>@todo FIXME: if we keep noise, noise_mult_factor, and advance_method,
!> add them to the attributes on the diagnostic file.

call nc_write_location_atts(ncid, model_size)
call nc_enddef(ncid)
call nc_write_location(ncid, state_loc, model_size)

call nc_sync(ncid)

end subroutine nc_write_model_atts

!===================================================================

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

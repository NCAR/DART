! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the model module for the Lorenz 96 2-scale model, documented in
! Lorenz (1995).  It also has the option of the variant on the model 
! from Smith (2001), and it is invoked by setting local_y = .true. in the 
! namelist.  The time step, coupling, forcing, number of X variables, and the
! number of Ys per X are all specified in the namelist.  Defaults are chosen
! depending on whether the Lorenz or Smith option is specified in the namelist.
! Lorenz is the default model.
!
! May 06 2004, modification of T. Hoar's L96 1-scale model by Josh Hacker

use types_mod,             only : r8, i8, i4
use time_manager_mod,      only : time_type, set_time
use location_mod,          only : location_type, set_location, get_location,  &
                                  get_close_obs, get_close_state,             &
                                  convert_vertical_obs, convert_vertical_state
use utilities_mod,         only : register_module, do_nml_file, do_nml_term,    &
                                  nmlfileunit, find_namelist_in_file,           &
                                  check_namelist_read, E_ERR, error_handler
use location_io_mod,       only : nc_add_location_atts
use netcdf_utilities_mod,  only : nc_add_global_attribute, nc_synchronize_file, nc_put_variable, &
                                  nc_add_global_creation_time, nc_begin_define_mode, nc_end_define_mode, &
                                  nc_add_attribute_to_variable, nc_define_dimension, &
                                  nc_define_real_variable, nc_define_integer_variable
use         obs_kind_mod,  only : QTY_LARGE_SCALE_STATE, QTY_SMALL_SCALE_STATE
use ensemble_manager_mod,  only : ensemble_type
use distributed_state_mod, only : get_state
use state_structure_mod,   only : add_domain, add_dimension_to_variable, finished_adding_domain
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

! Basic model parameters controlled by namelist

integer  :: model_size_x = 36
integer  :: y_per_x    = 10
real(r8) :: forcing    = 15.0_r8
real(r8) :: delta_t    = 0.005_r8
real(r8) :: coupling_b = 10.0_r8
real(r8) :: coupling_c = 10.0_r8
real(r8) :: coupling_h = 1.0_r8
logical  :: local_y = .false.  ! default Lorenz' approach
integer  :: time_step_days = 0
integer  :: time_step_seconds = 3600
character(len=256) :: template_file = 'null'

namelist /model_nml/ model_size_x, y_per_x, forcing, delta_t, &
                     coupling_b, coupling_c, coupling_h, &
                     local_y, time_step_days, time_step_seconds, &
                     template_file


! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step

type static_data
  integer               :: x_size, y_size
  integer(i8)           :: model_size
  real(r8)              :: b, c, h, f
  real(r8), allocatable :: x_loc(:), y_loc(:)
end type static_data

type(static_data) :: l96

contains

!------------------------------------------------------------------
!> Initializes class data for this model. Sets the location of the
!> state variables, initializes the time type and state vector structure.

subroutine static_init_model()

integer  :: i, iunit, io, kcount, dom_id

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! constants
l96%h = coupling_h
l96%b = coupling_b
l96%c = coupling_c
l96%f = forcing

! size of domain and state
l96%x_size = model_size_x
l96%y_size = model_size_x*y_per_x
l96%model_size = l96%x_size + l96%y_size

! Create storage for locations
allocate(l96%x_loc(l96%x_size))
allocate(l96%y_loc(l96%y_size))
allocate(state_loc(l96%model_size))

! Define the locations of the model state variables
! Carrying all three for now - may not be necessary
kcount = 1
do i = 1, l96%x_size
   l96%x_loc(i) = (i - 1.0_r8) / l96%x_size
   state_loc(kcount) =  set_location(l96%x_loc(i))
   kcount = kcount+1
end do
do i = 1, l96%y_size
   l96%y_loc(i) = (i - 1.0_r8) / l96%y_size
   state_loc(kcount) =  set_location(l96%y_loc(i))
   kcount = kcount+1
end do

! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.
time_step = set_time(time_step_seconds, time_step_days)

! Define the domain so the dart library code can read/write it.

if (template_file /= "null") then
   ! If we are have specified a template file then DART will use it to discover
   ! the sizes and dimensions of the given variables.  Tell DART which NetCDF variable
   ! names to read into the model state, and what quantities they correspond to.
   ! Note that the state data will be read from the initial conditions file; 
   ! the template file is only used to determine sizes and attributes.
   dom_id = add_domain(template_file, 2, (/ "X", "Y" /), (/ QTY_LARGE_SCALE_STATE, QTY_SMALL_SCALE_STATE /))
else
   ! If we are using the built-in initial conditions instead of reading a
   ! restart file (for example, running perfect_model_obs with different model
   ! sizes and so not reading from an existing restart file), then define the
   ! variable names and sizes to DART explicitly.  The values in the state vector
   ! will come from the call to init_conditions().
   dom_id = add_domain(2, (/ "X", "Y" /), (/ QTY_LARGE_SCALE_STATE, QTY_SMALL_SCALE_STATE /))
   call add_dimension_to_variable(dom_id, 1, "Xlocation", l96%x_size)
   call add_dimension_to_variable(dom_id, 2, "Ylocation", l96%y_size)
   call finished_adding_domain(dom_id)
endif

end subroutine static_init_model


!------------------------------------------------------------------
!> subroutine comp_dt(x, dt) (note used for both x and y together)
!> Computes the time tendency of the lorenz 1996 model given current state

subroutine comp_dt(x, dt)

real(r8), intent( in) ::  x(:)
real(r8), intent(out) :: dt(:)

integer :: j, jp1, jp2, jm1, jm2 
integer :: k, kp1,      km1, km2 
integer :: jk
integer :: xs, xe, ys, ye, js, je
real(r8) :: fast_sum, c1, c2, c3
real(r8), dimension(y_per_x) :: tmp_Y, tmp_dt

c1 = l96%c * l96%b
c2 = l96%c
c3 = l96%h * l96%c / l96%b

ys = l96%x_size + 1
ye = l96%model_size
! first small-scale, then large

if ( local_y ) then

! Smith's version treats each group of Y variables independently (could be
! a mistake in the paper)

   do jk = ys, ye, y_per_x
      tmp_Y = x(jk:jk+y_per_x-1)
      do j = 1, y_per_x
         jp1 = j + 1
         if(jp1 > y_per_x) jp1 = 1
         jp2 = j + 2
         if(jp2 > y_per_x) jp2 = jp2 - y_per_x 
         jm2 = j - 2
         if(jm2 < 1) jm2 = y_per_x + jm2
         jm1 = j - 1
         if(jm1 < 1) jm1 = y_per_x
         tmp_dt(j) = c1 * tmp_Y(jp1) * (tmp_Y(jm1)-tmp_Y(jp2)) - c2 * tmp_Y(j) &
             + c3 * x(int(real(jk+j-1-ys,r8)/real(y_per_x,r8))+1)
      enddo
      dt(jk:jk+y_per_x-1) = tmp_dt
   enddo
else
! Lorenz' version treats all the Y variables together so that they
! are aware of each other across the Xs

   do j = ys, ye
      jp1 = j + 1
      if(jp1 > ye) jp1 = ys
      jp2 = j + 2
      if(jp2 > ye) jp2 = ys + jp2 - ye - 1
      jm2 = j - 2
      if(jm2 < ys) jm2 = ye + jm2 - ys + 1
      jm1 = j - 1
      if(jm1 < ys) jm1 = ye
      dt(j) = c1 * x(jp1) * (x(jm1)-x(jp2)) - c2 * x(j) &
           + c3 * x(int(real(j-ys,r8)/real(y_per_x,r8))+1)
   enddo
endif

! Now large scale
xs = 1
xe = l96%x_size
do k = xs, xe
   js = (k-1)*y_per_x + l96%x_size + 1
   je = k*y_per_x + l96%x_size 
   fast_sum = sum(x(js:je))
   kp1 = k + 1
   if(kp1 > xe) kp1 = 1
   km2 = k - 2
   if(km2 < 1) km2 = xe + km2
   km1 = k - 1
   if(km1 < 1) km1 = xe
   dt(k) = (x(kp1) - x(km2)) * x(km1) - x(k) + l96%f &
           - c3 * fast_sum
end do

end subroutine comp_dt



!------------------------------------------------------------------
!> Initial conditions for model if not starting from a restart file.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

! Need noise in both the large and small scales.  The small-scale noise
! should be present in each group of Ys in case the local_y option is selected

x(1:l96%x_size) = l96%f
x(l96%x_size+1:l96%model_size) = 0.0_r8
x(1) = 1.001_r8 * l96%f

x(l96%x_size+1:) = 0.01_r8 * x(2)
x(l96%x_size+1:l96%model_size:y_per_x) = 0.011_r8 * x(2) 

end subroutine init_conditions


!------------------------------------------------------------------
!> Does single time step advance for lorenz 96 2-scale model
!> using four-step rk time step

subroutine adv_1step(x, time)

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

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
!> Returns size of model

function get_model_size()

integer(i8) :: get_model_size

get_model_size = l96%model_size

end function get_model_size


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
integer(i8) :: base_index, top_index
real(r8) :: lctn, lctnfrac
real(r8) :: x_lower(ens_size) !< the lower piece of state vector
real(r8) :: x_upper(ens_size) !< the upper piece of state vector

! All forward operators supported
istatus(:) = 0

! Convert location to real
lctn = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
if ( itype == QTY_LARGE_SCALE_STATE ) then
   lctn = l96%x_size * lctn
   base_index = 1
   top_index = l96%x_size
elseif ( itype == QTY_SMALL_SCALE_STATE ) then
   lctn = l96%y_size * lctn
   base_index = l96%x_size + 1
   top_index = l96%model_size
else
   call error_handler(E_ERR,'model_interpolate', 'cannot handle this type', &
                     source, revision, revdate)
endif

lower_index = int(lctn) + base_index
upper_index = lower_index + 1
if(lower_index > top_index) lower_index = lower_index - top_index
if(upper_index > top_index) upper_index = upper_index - top_index

lctnfrac = lctn - int(lctn)

! Lower value
x_lower(:) = get_state(lower_index, state_handle)

! Upper value
x_upper(:) = get_state(upper_index, state_handle)

! calculate the obs value

expected_obs(:) = (1.0_r8 - lctnfrac) * x_lower(:) + lctnfrac * x_upper(:)

end subroutine model_interpolate


!------------------------------------------------------------------
!> Return the minimum amount of time the model can be advanced by.
!> This is unrelated to the internal RK timesteps.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations


!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location and optionally quantity.  THESE ARE USING
!> NON-STANDARD QUANTITY TYPES

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type                                      
location = state_loc(index_in)
if (present(var_type)) then
  var_type = QTY_LARGE_SCALE_STATE
  if ( index_in > model_size_x ) var_type = QTY_SMALL_SCALE_STATE
endif

end subroutine get_state_meta_data


!------------------------------------------------------------------
! Writes model-specific attributes to a netCDF file

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid
integer, intent(in) :: domain_id

integer :: i, xs, xe, ys, ye
integer(i8) :: indx
type(location_type) :: lctn

real(r8), allocatable :: xloc(:),yloc(:)

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate", revdate )

call nc_add_global_attribute(ncid, "model", "Lorenz_96_2scale")
call nc_add_global_attribute(ncid, "model_delta_t", delta_t )
call nc_add_global_attribute(ncid, "model_forcing", l96%f )
call nc_add_global_attribute(ncid, "model_coupling_b", l96%b )
call nc_add_global_attribute(ncid, "model_coupling_c", l96%c )
call nc_add_global_attribute(ncid, "model_coupling_h", l96%h )
call nc_add_global_attribute(ncid, "slow variables (large scale)", "X")
call nc_add_global_attribute(ncid, "fast variables (small scale)", "Y")


! Define the dimensions IDs for X and Y location dimensions
call nc_define_dimension(ncid, "Xlocation", l96%x_size)
call nc_define_dimension(ncid, "Ylocation", l96%y_size)

call nc_define_real_variable(ncid, "Xlocation", "Xlocation")
call nc_add_location_atts   (ncid, "Xlocation")
call nc_define_real_variable(ncid, "Ylocation", "Ylocation")
call nc_add_location_atts   (ncid, "Ylocation")

! Leave define mode so we can fill
call nc_end_define_mode(ncid)

! The starting and ending indices of the X vars and Y vars
! in the state vector
xs = 1
xe = l96%x_size
ys = xe + 1
ye = xe + l96%y_size

allocate(xloc(l96%x_size), yloc(l96%y_size))

! Fill the location variables.

do i = xs, xe
   indx = i ! convert to an i8
   call get_state_meta_data(indx,lctn)
   xloc(i) = get_location(lctn)
enddo

do i = ys, ye
   indx = i ! convert to an i8
   call get_state_meta_data(indx,lctn)
   yloc(i-ys+1) = get_location(lctn)
enddo

call nc_put_variable(ncid, "Xlocation", xloc, "nc_write_model_atts")
call nc_put_variable(ncid, "Ylocation", yloc, "nc_write_model_atts")
call nc_synchronize_file(ncid)

deallocate(xloc, yloc)

end subroutine nc_write_model_atts

!-------------------------------------------------------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Interface to Lorenz-04 3-variable model

module model_mod

use             types_mod, only : r8, i8, i4

use      time_manager_mod, only : time_type, set_time

use          location_mod, only : location_type, set_location, get_location, &
                                  get_close_obs, get_close_state, &
                                  convert_vertical_obs, convert_vertical_state

use         utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, nmlfileunit, &
                                  do_output, find_namelist_in_file, check_namelist_read,     &
                                  do_nml_file, do_nml_term

use  netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                  nc_add_global_creation_time, &
                                  nc_begin_define_mode, nc_end_define_mode

use          obs_kind_mod, only : QTY_STATE_VARIABLE

use  ensemble_manager_mod, only : ensemble_type

use distributed_state_mod, only : get_state

use   state_structure_mod, only : add_domain

use      dart_time_io_mod, only : read_model_time, write_model_time

use       location_io_mod, only : nc_write_location_atts, nc_get_location_varids, &
                                  nc_write_location

use     default_model_mod, only : end_model, pert_model_copies, nc_write_model_vars, &
                                  init_time

implicit none
private

! these routines must be public and you cannot change the
! arguments because they will be called *from* other DART code.

!> required routines with code in this module
public :: static_init_model, &
          init_conditions, &
          init_time, &
          adv_1step, &
          get_model_size, &
          get_state_meta_data, &
          model_interpolate, &
          shortest_time_between_assimilations, &
          nc_write_model_atts

!> required routines where code is in other modules
public :: nc_write_model_vars, &
          pert_model_copies, &
          read_model_time, &
          write_model_time, &
          get_close_obs, &
          get_close_state, &
          convert_vertical_obs, &
          convert_vertical_state, &
          end_model

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!---------------------------------------------------------------
! Basic model parameters controlled by nameslist; have defaults
!
integer(i8) :: model_size        = 960_i8
real(r8)    :: forcing           = 15.00_r8
real(r8)    :: delta_t           = 0.001_r8
real(r8)    :: space_time_scale  = 10.00_r8
real(r8)    :: coupling          = 3.00_r8
integer     :: K                 = 32
integer     :: smooth_steps      = 12
integer     :: time_step_days    = 0
integer     :: time_step_seconds = 3600
integer     :: model_number      = 3 ! (2 for single scale, 3 for 2-scale, Lorenz 05)

namelist /model_nml/ model_size, forcing, delta_t, space_time_scale, coupling, K, &
       smooth_steps, time_step_days, time_step_seconds, model_number

!---------------------------------------------------------------- 
! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step

!---------------------------------------------------------------- 
! Define the averaging function for the production of x from z (L2k4)
real(r8), allocatable :: a(:)

!---------------------------------------------------------------- 
! Define some parameters for computational efficiency
integer  :: H
integer  :: K2
integer  :: K4
integer  :: ss2
real(r8) :: sts2

contains

!==================================================================


!------------------------------------------------------------------
!> Initializes class data for this model. For now, simply outputs the
!> identity info, sets the location of the state variables, and initializes
!> the time type for the time stepping (is this general enough for time???)

subroutine static_init_model()

real(r8) :: x_loc
real(r8) :: ri
real(r8) :: alpha, beta
integer  :: i, iunit, io
integer  :: j, dom_id

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! if K is even, H = K/2 and the first/last summation 
! terms are divided by 2.  if K is odd, H = (K-1)/2 and 
! the first and last terms are taken as-is.  this code
! only implements the algorithm for even K so test for it.
if (int((K+1)/2) /= int(K/2)) then
   call error_handler(E_ERR,'static_init_model',&
         'Model only handles even values of K', &
          source, revision, revdate)
endif

! Create storage for locations
allocate(state_loc(model_size))

! Define the locations of the model state variables
do i = 1, model_size
   x_loc = (i - 1.0_r8) / model_size
   state_loc(i) =  set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.
time_step = set_time(time_step_seconds, time_step_days)

! Generate the alpha and beta parameters for the calculation of "a"
alpha = (3.0_r8*(smooth_steps**2) + 3.0_r8) &
      / (2.0_r8*(smooth_steps**3) + 4.0_r8*smooth_steps)
beta  = (2.0_r8*(smooth_steps**2) + 1.0_r8) &
      / (1.0_r8*(smooth_steps**4) + 2.0_r8*(smooth_steps**2))

! The "a" vector is a smoothing filter for the production of x and y from z
! in L2k4. Apologies for the "ri" and "j" construct
allocate(a(2*smooth_steps + 1))
ri = - smooth_steps - 1.0_r8
j = 0
do i = - smooth_steps, smooth_steps
   j = j + 1
   ri = ri + 1.0_r8
   a(j) = alpha - beta*abs(ri)
end do

! defining parameters to help reduce the number of operations in the calculation
! of dz/dt
H    = K/2
K2   = 2*K
K4   = 4*K
ss2  = 2*smooth_steps
sts2 = space_time_scale**2   

dom_id = add_domain(model_size)

end subroutine static_init_model


!------------------------------------------------------------------
!> subroutine init_conditions(x)
!>
!> Initial conditions for lorenz 04
!> It is assumed that this is called before any other routines in this
!> module. Should probably make that more formal and perhaps enforce for
!> more comprehensive models.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

x    = forcing
x(1) = 1.001_r8 * forcing

end subroutine init_conditions



!------------------------------------------------------------------
!> Does single time step advance for lorenz 04 model
!> using four-step rk time step

subroutine adv_1step(x, time)

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

real(r8), dimension(size(x)) :: x1, x2, x3, x4, dx, inter
real(r8), dimension(size(x)) :: dxt

call comp_dt(x, dx)    !  Compute the first intermediate step
x1    = delta_t * dx
inter = x + x1 / 2.0_r8

call comp_dt(inter, dx)!  Compute the second intermediate step
x2    = delta_t * dx
inter = x + x2 / 2.0_r8

call comp_dt(inter, dx)!  Compute the third intermediate step
x3    = delta_t * dx
inter = x + x3

call comp_dt(inter, dx)!  Compute fourth intermediate step
x4 = delta_t * dx

!  Compute new value for x

dxt = x1/6.0_r8 + x2/3.0_r8 + x3/3.0_r8 + x4/6.0_r8
x = x + dxt

end subroutine adv_1step


!------------------------------------------------------------------
!> Computes the time tendency of the lorenz 2004 model given current state.
!>
!> The model equations are given by
!> 
!> Model 2 (II)
!>      dX_i
!>      ---- = [X,X]_{K,i} -  X_i + F 
!>       dt                
!>                         
!>
!> Model 3 (III)
!>      dZ_i
!>      ---- = [X,X]_{K,i} + b^2 (-Y_{i-2}Y_{i-1} + Y_{i-1}Y_{i+1})
!>       dt                +  c  (-Y_{i-2}X_{i-1} + Y_{i-1}X_{i+1})
!>                         -  X_i - b Y_i + F,
!>
!> where
!>
!>     [X,X]_{K,i} = -W_{i-2K}W_{i-K} 
!>                 +  sumprime_{j=-(K/2)}^{K/2} W_{i-K+j}X_{i+K+j}/K,
!>
!>      W_i =  sumprime_{j=-(K/2)}^{K/2} X_{i-j}/K,
!>
!> and sumprime denotes a special kind of summation where the first
!> and last terms are divided by 2.
!>
!> NOTE: The equations above are only valid for K even.  If K is odd,
!> then sumprime is replaced by the traditional sum, and the K/2 limits
!> of summation are replaced by (K-1)/2. THIS CODE ONLY IMPLEMENTS THE
!> K EVEN SOLUTION!!!
!>
!> The variable that is integrated is X (model II) or Z (model III), 
!> but the integration of Z requires
!> the variables X and Y.  For model III they are obtained by
!>
!>      X_i = sumprime_{j= -J}^{J} a_j Z_{i+j}
!>      Y_i = Z_i - X_i.
!>
!> The "a" coefficients are given by
!>
!>      a_j = alpha - beta |j|,
!> 
!> where
!>
!>      alpha = (3J^2 + 3)/(2J^3 + 4J)
!>      beta  = (2J^2 + 1)/(1J^4 + 2J^2).
!>
!> This choice of alpha and beta ensures that X_i will equal Z_i
!> when Z_i varies quadratically over the interval 2J.   This choice
!> of alpha and beta means that sumprime a_j = 1 and 
!> sumprime (j^2) a_j = 0.
!>
!> Note that the impact of this filtering is to put large-scale
!> variations into the X variable, and small-scale variations into
!> the Y variable.
!> 
!> The parameter names above are based on those that appear in
!> Lorenz 04.  To map to the code below, set:
!>
!>       F = forcing 
!>       b = space_time_scale
!>       c = coupling
!>       K = K
!>       J = smooth_steps
subroutine comp_dt(z, dt) 

real(r8), intent( in)        ::  z(:)
real(r8), intent(out)        :: dt(:)
real(r8), dimension(size(z)) :: x, y
real(r8)                     :: xwrap(- K4:model_size + K4)
real(r8)                     :: ywrap(- K4:model_size + K4)
real(r8)                     ::    wx(- K4:model_size + K4)
real(r8)                     :: xx
integer                      :: i, j

! could branch this differently for more effecient model II

if ( model_number == 3 ) then
   ! Decompose z into x and y
   call z2xy(z,x,y)
elseif ( model_number == 2 ) then
   x = z
   y = 0.0_r8   ! just a dummy
else
   call error_handler(E_ERR,'comp_dt',&
         'Do not know that model number', source, revision, revdate)
endif

! Deal with cyclic boundary conditions using buffers
do i = 1, model_size
   xwrap(i) = x(i)
   ywrap(i) = y(i)
end do

! Fill the xwrap and ywrap buffers
do i = 1, K4
   xwrap(- K4 + i)       = xwrap(model_size - K4 +i)
   xwrap(model_size + i) = xwrap(i)
   ywrap(- K4 + i)       = ywrap(model_size - K4 +i)
   ywrap(model_size + i) = ywrap(i)
end do

! Calculate the W's
do i = 1, model_size
   wx(i) = xwrap(i - (-H))/2.00_r8
   do j = - H + 1, H - 1
      wx(i) = wx(i) + xwrap(i - j)
   end do
   wx(i) = wx(i) + xwrap(i - H)/2.00_r8
   wx(i) = wx(i)/K
end do

! Fill the W buffers
do i = 1, K4
   wx(- K4 + i)       = wx(model_size - K4 + i)
   wx(model_size + i) = wx(i)
end do

! Generate dz/dt
do i = 1, model_size
   xx = wx(i - K + (-H))*xwrap(i + K + (-H))/2.00_r8
   do j = - H + 1, H - 1
      xx = xx + wx(i - K + j)*xwrap(i + K + j)
   end do
   xx = xx + wx(i - K + H)*xwrap(i + K + H)/2.00_r8
   xx = - wx(i - K2)*wx(i - K) + xx/K
      
   if ( model_number == 3 ) then
     dt(i) = xx + (sts2)*( - ywrap(i - 2)*ywrap(i - 1) &
         + ywrap(i - 1)*ywrap(i + 1)) + coupling*( - ywrap(i - 2)*xwrap(i - 1) &
         + ywrap(i - 1)*xwrap(i + 1)) - xwrap(i) - space_time_scale*ywrap(i) &
         + forcing
   else ! must be model II
     dt(i) = xx - xwrap(i) + forcing
   endif

end do

end subroutine comp_dt


!------------------------------------------------------------------
!> Decomposes z into x and y for L2k4

subroutine z2xy(z,x,y)

integer :: i, j, ia
real(r8), intent( in) :: z(:)
real(r8), intent(out) :: x(:)
real(r8), intent(out) :: y(:)
real(r8)              :: zwrap(- ss2:model_size + ss2)

! Fill zwrap
do i = 1, model_size
   zwrap(i) = z(i)
end do
zwrap( - ss2) = zwrap(model_size - ss2)
do i = 1, ss2
   zwrap( - ss2 + i) = zwrap(model_size - ss2 + i)
   zwrap(model_size + i) = zwrap(i)
end do

! Generate the x variables
do i = 1, model_size
   ia = 1
   x(i) = a(ia)*zwrap(i - ( - smooth_steps))/2.00_r8
   do j = - smooth_steps + 1, smooth_steps - 1
      ia = ia + 1
      x(i) = x(i) + a(ia)*zwrap(i - j)
   end do
   ia = ia + 1
   x(i) = x(i) + a(ia)*zwrap(i - smooth_steps)/2.00_r8
end do

! Generate the y variables
do i = 1, model_size
   y(i) = z(i) - x(i)
end do

end subroutine z2xy


!------------------------------------------------------------------
!> Returns size of model

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------
!> Interpolates from state vector to the location. 

!> Argument itype is not used here because there is only one type of variable.
!> Type is needed to allow swap consistency with more complex models.

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
lctn = model_size * lctn

lower_index = int(lctn) + 1
upper_index = lower_index + 1
if(lower_index > model_size) lower_index = lower_index - model_size
if(upper_index > model_size) upper_index = upper_index - model_size

lctnfrac = lctn - int(lctn)
expected_val(:) = (1.0_r8 - lctnfrac) * get_state(lower_index, state_handle) + &
                            lctnfrac  * get_state(upper_index, state_handle)

if(1 == 1) return
! code returns here - the code below in this subroutine is not executed.

!!!expected_val(:) = expected_val(:) ** 2
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
expected_val = 0.0_r8
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
!> Returns the the time step of the model. In the long run should be repalced
!> by a more general routine that returns details of a general time-stepping
!> capability.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations


!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location and optionally quantity.


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

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate", revdate )
call nc_add_global_attribute(ncid,  "model", "Lorenz_04")

if ( model_number == 2 ) then
   call nc_add_global_attribute(ncid, "model_scale", "single")
else if ( model_number == 3 ) then
   call nc_add_global_attribute(ncid, "model_scale", "2-scale")
endif
call nc_add_global_attribute(ncid, "model_forcing", forcing )
call nc_add_global_attribute(ncid, "model_delta_t", delta_t )
call nc_add_global_attribute(ncid, "space_time_scale", space_time_scale )
call nc_add_global_attribute(ncid, "coupling", coupling )
call nc_add_global_attribute(ncid, "K", K )
call nc_add_global_attribute(ncid, "smooth_steps", smooth_steps )
call nc_add_global_attribute(ncid, "time_step_days", time_step_days )
call nc_add_global_attribute(ncid, "time_step_seconds", time_step_seconds )

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

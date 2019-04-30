! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

use types_mod,             only : r8, i8, i4

use location_mod,          only : location_type, set_location, get_location, &
                                  get_close_obs, get_close_state, &
                                  convert_vertical_obs, convert_vertical_state

use utilities_mod,         only : register_module, error_handler, E_ERR, E_MSG, do_output, &
                                  nmlfileunit, find_namelist_in_file, check_namelist_read, &
                                  do_nml_file, do_nml_term

use netcdf_utilities_mod,  only : nc_add_global_attribute, &
                                  nc_add_global_creation_time, &
                                  nc_begin_define_mode, nc_end_define_mode

use location_io_mod,       only : nc_write_location_atts, nc_get_location_varids, &
                                  nc_write_location

use obs_kind_mod,          only : QTY_STATE_VARIABLE

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain

use dart_time_io_mod,      only : read_model_time, write_model_time

use default_model_mod,     only : end_model, pert_model_copies, nc_write_model_vars, &
                                  init_time

use random_seq_mod,        only : random_seq_type, random_gaussian, &
                                  init_random_seq, several_random_gaussians

use time_manager_mod,      only : time_type, set_time

implicit none
private

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
          get_close_obs, &
          get_close_state, &
          init_time, &
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

integer(i8), parameter :: model_size = 9

!  define model parameters
! c is sqrt(0.75)

real(r8), private, parameter :: a(3) = (/  1.0_r8,  1.0_r8, 3.0_r8 /), &
                                b(3) = (/ -1.5_r8, -1.5_r8, 0.5_r8 /), &
                                f(3) = (/ 0.10_r8,  0.0_r8, 0.0_r8 /), &
                                h(3) = (/ -1.0_r8,  0.0_r8, 0.0_r8 /), &
                                nu = 1.0_r8 / 48.0_r8, kappa = nu, c = 0.8660254_r8

! Define the location of the state variables in module storage
! This is used for distance dependence stuff in more general models but is 
! currently just set to give 0 separation distance for all 9 vars here.

type(location_type) :: state_loc(model_size)
type(time_type)     :: time_step

! Need reproducible sequence of noise added so that different runs
! can be cleanly compared

logical :: first_ran_call = .true.
type(random_seq_type) :: ran_seq


! Namelist with default values

real(r8) :: g = 8.0_r8  ! lorenz default; 9.90_r8 is a higher dimension attractor

real(r8) :: deltat = 1.0_r8 / 12.0_r8     ! model time step
integer :: time_step_days = 0
integer :: time_step_seconds = 3600

logical :: add_noise = .false.

namelist /model_nml/ g, deltat, time_step_days, time_step_seconds, add_noise



contains

!----------------------------------------------------------------------
!> Initializes class data for this model. For now, simply outputs the
!> identity info, sets the location of the state variables, and initializes
!> the time type for the time stepping (is this general enough for time???)

subroutine static_init_model()

real(r8) :: x_loc
integer :: i, iunit, io, dom_id

! Register the module into the logfile
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

! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.
time_step = set_time(time_step_seconds, time_step_days)

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
dom_id = add_domain(model_size)

end subroutine static_init_model


!----------------------------------------------------------------------
!> computes time tendency for the 9 variable lorenz model given the
!> current values of the 9 variables (x, y, z each with three elements)
!> and the values of a number of parameters.

subroutine comp_dt(xxx, dxxx)

real(r8), intent(in)  :: xxx(:)
real(r8), intent(out) :: dxxx(:)

real(r8) :: x(3), y(3), z(3), dx(3), dy(3), dz(3)
real(r8) :: rnum(9)
integer  :: i, j, k

!  unpack the 9-vectors into the x, y and z 3-vectors

call unpack9var(xxx, x, y, z)

!  equations 33-35 from lorenz, 1980, jas, p.1688
!  equations are defined with cyclic indices

do i = 1, 3
   J = MOD(I, 3) + 1
   K = MOD(I + 1, 3) + 1

   !  BEGIN WITH THE FIRST EQUATION (33)

   dx(i) = (a(i)*b(i)*x(j)*x(k) - c*(a(i) - a(k))*x(j)*y(k) + &
           c*(a(i) - a(j))*y(j)*x(k) - 2*c**2*y(j)*y(k) - &
           nu*a(i)**2*x(i) + a(i)*y(i) - a(i)*z(i)) / a(i)

   !  equation (34)

   dy(i) = (-1.*a(k)*b(k)*x(j)*y(k) - a(j)*b(j)*y(j)*x(k) + &
           c*(a(k) - a(j))*y(j)*y(k) - a(i)*x(i) - nu*a(i)**2*y(i)) / a(i)

   !  equation (35)

   dz(i) = -1.*b(k)*x(j)*(z(k) - h(k)) - &
           b(j)*(z(j) - h(j))*x(k) + c*y(j)*(z(k) -h(k)) - &
           c*(z(j) - h(j))*y(k) + g*a(i)*x(i) - kappa*a(i)*z(i) + f(i)
end do   

call pack9var(dx, dy, dz, dxxx)     !  pack the results into 9 vector

! We need to initialize the random gen for experiments to be repeatable
if(first_ran_call) then
   first_ran_call = .false.
   call init_random_seq(ran_seq)
end if

! ADDITION OF SOME NOISE AT 1/10 the amplitude of DT
if (add_noise) then
   call several_random_gaussians(ran_seq, 0.0_r8, dxxx(i) / 10.0_r8, 9, rnum)
   do i = 1, 9
      dxxx(i) = dxxx(i) + rnum(i)
   end do
endif

end subroutine comp_dt


!---------------------------------------------------------------------------
!> set of routines used to switch between 3 3-variable sets and 9-variable
!> set for lorenz 9 variable pe model

subroutine pack9var(x, y, z, pert)

real(r8), intent(in)  :: x(3), y(3), z(3)
real(r8), intent(out) :: pert(9)

pert(1:3) = x
pert(4:6) = y
pert(7:9) = z

end subroutine pack9var


!---------------------------------------------------------------------------
!> inverse of pack9var above

subroutine unpack9var(pert, x, y, z)

real(r8), intent(in)  :: pert(9)
real(r8), intent(out) :: x(3), y(3), z(3)

x = pert(1:3)
y = pert(4:6)
z = pert(7:9)

end subroutine unpack9var


!-------------------------------------------------------------------------
!> does one time step advance for 9 variable model using two-step rk.

subroutine adv_1step(x, time)

implicit none

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

real(r8) :: fract = 1.0_r8

call adv_single(x, fract)

end subroutine adv_1step


!-------------------------------------------------------------------------
!> does one time step advance for 9 variable model using two-step rk

subroutine adv_single(x, fract)

real(r8), intent(inout) :: x(:)
real(r8), intent(in)    :: fract

real(r8) :: x1(9), x2(9), dx(9)

!  compute the first intermediate step

call comp_dt(x, dx)
x1 = x + fract * deltat * dx

!  compute the second intermediate step

call comp_dt(x1, dx)
x2 = x1 + fract * deltat * dx

!  new value for x is average of original value and second intermediate

x = (x + x2) / 2.0_r8

end subroutine adv_single


!----------------------------------------------------------------------
!> generates set of random off-attractor initial conditions for 9 variable

subroutine init_conditions(x)

real(r8), intent(out) ::  x(:)

x = 0.10_r8

end subroutine init_conditions


!---------------------------------------------------------------------------
!> unused?
! application of balance eqtn as initialization method
! to the 9 variable lorenz model given the
! current values of the 9 variables (x, y, z each with three elements)
! and the values of a number of parameters.

subroutine balance_init(xxx, init_xxx)

real(r8), intent(in)  :: xxx(9)
real(r8), intent(out) :: init_xxx(9)

real(r8) :: x(3), y(3), z(3)
!!!real(r8) :: wkspce(3)
real(r8) :: lhs(3,3),rhs(3)
integer  :: i, j, k
!!!integer  :: ifail

!  unpack the 9-vectors into the x, y and z 3-vectors

call unpack9var(xxx, x, y, z)

! application of balance eqtn as initialization method, curry ET.AL, tellus,
!    1995, p. 153-154, section 3.3
! equation 11 from gent & mcwilliams, 1982, jas, p.4

DO i = 1, 3
   j    = mod(i, 3) + 1
   k    = mod(i + 1, 3) + 1
   z(i) = (a(i)*y(i)-2*c**2*y(j)*y(k))/a(i)
end do

!  equation 29/30 from gent & mcwilliams, 1982, jas, p.6

do i = 1, 3
   j = mod(i, 3) + 1
   k = mod(i + 1, 3) + 1
   lhs(i,i) = ( a(i)*a(j)*a(k)*(1+g*a(i))-2*c**2*(a(j)**2* &
      b(j)*y(j)**2+a(k)**2*b(k)*y(k)**2) )
   lhs(j,i) = -( a(j)*a(k)*( y(k)*(2*c**2-a(k)*b(k))+a(i)* &
      b(k)*(z(k)-h(k)) )+2*c**2*a(i)*a(j)*b(i)*y(i) *y(j)  )
   lhs(k,i) = -( a(j)*a(k)*( y(j)*(2*c**2-a(j)*b(j))+a(i)* &
      b(j)*(z(j)-h(j)) )+2*c**2*a(i)*a(k)*b(i)*y(i) *y(k)  )
   rhs(i) = a(j)*a(k)*(  c*(a(k)-a(j))*y(j)*y(k)+c*a(i)*( &
      (z(j)-h(j))*y(k)-y(j)*(z(k)-h(k)) )+a(i)*( nu* &
      a(i)*(z(i)-y(i))-f(i) )  )-2*c**2*(  c*a(j)* &
      (a(j)-a(i))*y(i)*y(j)**2+c*a(k)*(a(i)-a(k))*y(i) &
      *y(k)**2-nu*a(j)*a(k)*(a(j)+a(k))*y(j)*y(k)  )
end do

!!!ifail=0
!  n.a.g. ROUTINE TO CALCULATE APPROXIMATE SOLUTION X TO aX=B (I.E. lhs*X=rhs)
!!!call f04arf_wrap(lhs,3,rhs,3,x,wkspce,ifail)
!        print *,'XINIT=',X

!  pack the results into 9 vector

call pack9var(x, y, z, init_xxx)

end subroutine balance_init


!-------------------------------------------------------------------------
! Returns number of items in state vector

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size



!------------------------------------------------------------------------

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations


!---------------------------------------------------------------------

subroutine model_interpolate(state_handle, ens_size, location, itype, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: itype
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer(i8)  :: lower_index, upper_index
real(r8) :: lctn, lctnfrac

! All obs okay for now
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
expected_obs = (1.0_r8 - lctnfrac) * get_state(lower_index, state_handle) + &
                         lctnfrac  * get_state(upper_index, state_handle)

end subroutine model_interpolate


!---------------------------------------------------------------------
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
!> Writes the model-specific attributes to a netCDF file

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
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

call nc_add_global_attribute(ncid, "model", "9var")

call nc_add_global_attribute(ncid, "model_g", g )
call nc_add_global_attribute(ncid, "model_deltat", deltat )
call nc_add_global_attribute(ncid, "model_a", a )
call nc_add_global_attribute(ncid, "model_b", b )
call nc_add_global_attribute(ncid, "model_f", f )
call nc_add_global_attribute(ncid, "model_h", h )
call nc_add_global_attribute(ncid, "model_nu", nu )
call nc_add_global_attribute(ncid, "model_kappa", kappa )
call nc_add_global_attribute(ncid, "model_c", c )

call nc_write_location_atts(ncid, msize)
call nc_end_define_mode(ncid)
call nc_write_location(ncid, state_loc, msize)

end subroutine nc_write_model_atts

!--------------------------------------------------------------------

!===================================================================
! End of 9var model_mod 
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

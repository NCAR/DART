module model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

use types_mod
use location_mod, only : location_type, get_dist, set_location, get_location
use time_manager_mod

use random_seq_mod,   only : random_gaussian, random_seq_type, &
                            init_random_seq, several_random_gaussians

private

public init_model, get_model_size, init_conditions, adv_1step, advance, &
   model_output, balance_init , static_init_model, init_time, &
   get_state_meta_data, get_model_time_step, end_model, model_interpolate

integer, parameter :: model_size = 9

!  define model parameters
! c is sqrt(0.75)

real(r8), private, parameter :: a(3) = (/  1.0_r8,  1.0_r8, 3.0_r8 /), &
                                b(3) = (/ -1.5_r8, -1.5_r8, 0.5_r8 /), &
                                f(3) = (/ 0.10_r8,  0.0_r8, 0.0_r8 /), &
                                h(3) = (/ -1.0_r8,  0.0_r8, 0.0_r8 /), &
                                nu = 1.0_r8 / 48.0_r8, kappa = nu, c = 0.8660254_r8


 real(r8), parameter :: g = 8.0_r8     ! lorenz default
!real(r8), parameter :: g = 9.90_r8    ! higher dimension attractor

real(r8), parameter :: deltat = 1.0_r8 / 12.0_r8     ! model time step

! Define the location of the state variables in module storage
! This is used for distance dependence stuff in more general models but is 
! currently just set to give 0 separation distance for all 9 vars here.

type(location_type) :: state_loc(model_size)
type(time_type)     :: time_step

! Need reproducible sequence of noise added so that different runs
! can be cleanly compared

logical :: first_ran_call = .true.
type(random_seq_type) :: ran_seq


contains

!======================================================================


subroutine static_init_model()
!----------------------------------------------------------------------
! subroutine static_init_model()
!
! Initializes class data for this model. For now, simply outputs the
! identity info, sets the location of the state variables, and initializes
! the time type for the time stepping (is this general enough for time???)

implicit none
real(r8) :: x_loc
integer :: i
character(len=128) :: source,revision,revdate

! let CVS fill strings ... DO NOT EDIT ...

source   = "$Source$"
revision = "$Revision$"
revdate  = "$Date$"

! Ultimately,  change output to diagnostic output block ...

write(*,*)'model attributes:'
write(*,*)'   ',source
write(*,*)'   ',revision
write(*,*)'   ',revdate

! Define the locations of the model state variables
do i = 1, model_size
   x_loc = (i - 1.0) / model_size
   state_loc(i) =  set_location(x_loc)
end do


! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.
time_step = set_time(3600, 0)


end subroutine static_init_model




  subroutine comp_dt(xxx, dxxx)
!----------------------------------------------------------------------
! subroutine comp_dt(xxx, dxxx)
!
! computes time tendency for the 9 variable lorenz model given the
! current values of the 9 variables (x, y, z each with three elements)
! and the values of a number of parameters.

implicit none

real(r8), intent(in)  :: xxx(:)
real(r8), intent(out) :: dxxx(:)

real(r8) :: x(3), y(3), z(3), dx(3), dy(3), dz(3), rnum(9)
integer  :: i, j, k

!  unpack the 9-vectors into the x, y and z 3-vectors

call unpack(xxx, x, y, z)

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

call pack(dx, dy, dz, dxxx)     !  pack the results into 9 vector

! OPTIONAL ADDItion OF NOISE
! ADDITION OF SOME NOISE AT 1/10 the amplitude of DT
! Do we need to initialize the repeatable random gen
!!!if(first_ran_call) then
!!!   first_ran_call = .false.
!!!   call init_random_seq(ran_seq)
!!!end if

!!!call several_random_gaussians(ran_seq, dble(0.0), dxxx(i) / 10.0, 9, rnum)
!!!do i = 1, 9
!!!   dxxx(i) = dxxx(i) + rnum(i)
!!!end do

end subroutine comp_dt




  subroutine pack(x, y, z, pert)
!---------------------------------------------------------------------------
! subroutine pack(x, y, z, pert)
!
! set of routines used to switch between 3 3-variable sets and 9-variable
! set for lorenz 9 variable pe model

implicit none

!  pack and unpack convert from x, y, z to full 9 vector format

real(r8), intent(in)  :: x(3), y(3), z(3)
real(r8), intent(out) :: pert(9)

pert(1:3) = x
pert(4:6) = y
pert(7:9) = z

end subroutine pack



  subroutine unpack(pert, x, y, z)
!---------------------------------------------------------------------------
! subroutine unpack(pert, x, y, z)
!
! inverse of pack above

implicit none

real(r8), intent(in)  :: pert(9)
real(r8), intent(out) :: x(3), y(3), z(3)

x = pert(1:3)
y = pert(4:6)
z = pert(7:9)

end subroutine unpack



  subroutine advance(x, num, xnew, time)
!-----------------------------------------------------------------------
! subroutine advance(x, num, xnew, time)
!
! advance advances the 9 variable model by a given number of steps

implicit none

real(r8), intent(in)  :: x(9)
integer,  intent(in)  :: num
real(r8), intent(out) :: xnew(9)
type(time_type), intent(in) :: time

integer :: i

xnew = x                 !  copy initial conditions to avoid overwrite

do i = 1, num            !  advance the appropriate number of steps
   call adv_1step(xnew, time)
end do

end subroutine advance





  subroutine adv_1step(x, time)
!-------------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! does one time step advance for 9 variable model using two-step rk.
! The Time argument is needed for compatibility with more complex models
! that need to know the time to compute their time tendency and is not
! used in L96. Is there a better way to do this in F90 than to just hang
! this argument out everywhere?

implicit none

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

real(r8) :: fract = 1.0_r8

call adv_single(x, fract)

end subroutine adv_1step




  subroutine adv_single(x, fract)
!-------------------------------------------------------------------------
! subroutine adv_single(x, fract)
!
! does one time step advance for 9 variable model using two-step rk

implicit none

real(r8), intent(inout) :: x(:)
real(r8), intent(in)    :: fract

real(r8) :: x1(9), x2(9), dx(9)
integer  :: i

!  compute the first intermediate step

call comp_dt(x, dx)
x1 = x + fract * deltat * dx

!  compute the second intermediate step

call comp_dt(x1, dx)
x2 = x1 + fract * deltat * dx

!  new value for x is average of original value and second intermediate

x = (x + x2) / 2.0_r8

end subroutine adv_single



  subroutine init_conditions(x)
!----------------------------------------------------------------------
! subroutine init_conditions(x)
!
! generates set of random off-attractor initial conditions for 9 variable

implicit none

real(r8), intent(out) ::  x(:)     ! TJH ... guessed at intent ...

integer  :: i
real(r8) :: x_loc

x = 0.10_r8

end subroutine init_conditions



  subroutine linearize(nl, l)
!----------------------------------------------------------------------
! subroutine linearize(nl, l)
!

implicit none

real(r8) :: nl(:), l(:, :)  ! TJH ... no intent ?

!  no-op subroutine header for linking in standard packages

end subroutine linearize



  subroutine balance_init(xxx, init_xxx)
!---------------------------------------------------------------------------
! subroutine balance_init(xxx, init_xxx)
!
! application of balance eqtn as initialization method
! to the 9 variable lorenz model given the
! current values of the 9 variables (x, y, z each with three elements)
! and the values of a number of parameters.

implicit none

real(r8), intent(in)  :: xxx(9)
real(r8), intent(out) :: init_xxx(9)

real(r8) :: x(3), y(3), z(3), dx(3), dy(3), dz(3), wkspce(3)
real(r8) :: lhs(3,3),rhs(3)
integer  :: i, j, k, ifail

!  unpack the 9-vectors into the x, y and z 3-vectors

call unpack(xxx, x, y, z)

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

ifail=0
!  n.a.g. ROUTINE TO CALCULATE APPROXIMATE SOLUTION X TO aX=B (I.E. lhs*X=rhs)
!!!call f04arf_wrap(lhs,3,rhs,3,x,wkspce,ifail)
!        print *,'XINIT=',X

!  pack the results into 9 vector

call pack(x, y, z, init_xxx)

end subroutine balance_init



  subroutine get_close_pts(list, num)
!-------------------------------------------------------------------------
! subroutine get_close_pts(list, num)
!

implicit none

integer, intent(in)    :: num
integer, intent(inout) :: list(model_size, num)

integer :: i, offset, index, temp

do i = 1, model_size

   do offset = -num/2, -num/2 + num - 1
      index = i + offset
      if(index > model_size) index = index - model_size
      if(index < 1         ) index = model_size + index
      list(i, offset + num/2 + 1) = index
   end do

   ! Always need the actual point first in list

   temp = list(i, 1)
   list(i, 1) =  list(i, num / 2 + 1)
   list(i, num / 2 + 1) = temp

end do

end subroutine get_close_pts



  function get_model_size()
!-------------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

integer :: get_model_size

get_model_size = model_size

end function get_model_size




function get_model_time_step()
!------------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



function model_interpolate(x, location, type)
!---------------------------------------------------------------------
!
! Interpolates from state vector x to the location. It's not particularly
! happy dumping all of this straight into the model. Eventually some
! concept of a grid underlying models but above locations is going to
! be more general. May want to wait on external infrastructure projects
! for this?

! Argument type is not used here because there is only one type of variable.
! Type is needed to allow swap consistency with more complex models.

implicit none

real(r8) :: model_interpolate
real(r8), intent(in) :: x(:)
type(location_type), intent(in) :: location
integer, intent(in) :: type

integer :: lower_index, upper_index
real(r8) :: loc, fraction

! Convert location to real
loc = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
loc = model_size * loc

lower_index = int(loc) + 1
upper_index = lower_index + 1
if(lower_index > model_size) lower_index = lower_index - model_size
if(upper_index > model_size) upper_index = upper_index - model_size

fraction = loc - int(loc)
model_interpolate = (1.0_r8 - fraction) * x(lower_index) + fraction * x(upper_index)

end function model_interpolate




subroutine get_state_meta_data(index, location)
!---------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?

implicit none

integer, intent(in) :: index
type(location_type), intent(out) :: location

location = state_loc(index)

end subroutine get_state_meta_data





  subroutine init_model()
!-------------------------------------------------------------------------
! subroutine init_model()
!
! Stub for model initialization, not needed for 9var

end subroutine init_model




subroutine init_time(time)
!----------------------------------------------------------------------
!
! Gets the initial time for a state from the model. Where should this info
! come from in the most general case?

implicit none

type(time_type), intent(out) :: time

! For now, just set to 0
time = set_time(0, 0)

end subroutine init_time



subroutine end_model()
!------------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Nothing for L96 for now.


end subroutine end_model






!===================================================================
! End of 9var model_mod 
!===================================================================

end module model_mod

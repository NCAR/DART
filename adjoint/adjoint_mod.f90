! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module adjoint_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! Does fake variational 4D fit; used to simulate optimal adjoint
! performance for comparison to other methods
! Makes use of NAG optimization routines, would have to be modified
! to run on platforms without NAG support.

use types_mod
use obs_mod,      only : take_obs
use model_mod,    only : advance, adv_1step
use nag_wrap_mod, only : e04dgf_wrap

! implicit none?
private

public min_cost, adj_forecast_eval, adj_forecast_out

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

real(r8), parameter :: small_dist = 0.001_r8, delta = 1e-10
integer, parameter :: max_1d = 100, max_iters = 100, num_refine = 20

! Parameters and storage for forecast evaluation
integer, parameter :: num_leads = 20
integer :: count = 0
real(r8) :: err(0:num_leads) = 0.0_r8

integer :: g_num_obs, g_n_times, g_obs_freq
real(r8), allocatable :: g_obs(:, :)

contains



  subroutine min_cost(x, obs, num_obs, n_times, obs_freq, xnew)
!========================================================================
! subroutine min_cost(x, obs, num_obs, n_times, obs_freq, xnew)
!
! xnew is value of best fit at LAST STEP (assimilated value)

implicit none

integer,  intent(in)  :: num_obs, n_times, obs_freq
real(r8), intent(in)  :: x(:), obs(num_obs, n_times)
real(r8), intent(out) :: xnew(:)

real(r8) :: cost, grad(size(x)), xp(size(x)), xmin(size(x)), objf
real(r8) :: objgrd(size(x)), work(13*size(x)), user(1)
integer  :: iwork(size(x) + 1), iuser(1), ifail, i, j, iter

!if(1 == 1) goto 20

ifail = -1
g_num_obs = num_obs
g_n_times = n_times
g_obs_freq = obs_freq
allocate(g_obs(num_obs, n_times))
g_obs = obs
xmin = x
write(*, *) 'xmin before eo4 ', xmin
call e04dgf_wrap(size(x), funct2, iter, objf, objgrd, &
                 xmin, iwork, work, iuser, user, ifail)
write(*, *) 'xmin after eo4 ', xmin
xnew = xmin
deallocate(g_obs)

if(1 == 0) then    ! TJH essentially comment out the following block

 20 cost = get_cost(x, obs, num_obs, n_times, obs_freq)
    write(*, *) 'initial cost is ', cost
    xnew = x
    do i = 1, max_iters
       call get_grad(xnew, obs, num_obs, n_times, obs_freq, grad) 
       call one_d_min(xnew, grad, obs, num_obs, n_times, obs_freq, xp)
       xnew = xp
       cost = get_cost(xnew, obs, num_obs, n_times, obs_freq)
       write(*, *) 'cost after iter ', i, ' is ', cost
    end do

endif

10 do i = 1, n_times - 1       ! Want to return value at last step

   write(*, *) 'advancing ', i
   do j = 1, size(xnew)
      write(*, *) 'xnew ', j, xnew(j)
   end do
   call advance(xnew, obs_freq, xnew)
   write(*, *) 'done advancing ', i

end do 
   
end subroutine min_cost


  subroutine adj_forecast_eval(x_in, adj_in, obs_freq)
!========================================================================
! subroutine adj_forecast_eval(x_in, adj_in, obs_freq)
!

implicit none

real(r8), intent(in) :: x_in(:), adj_in(:)
integer,  intent(in) :: obs_freq

real(r8) :: x(size(x_in)), adj(size(adj_in))
integer  :: i, j

count = count + 1
x     = x_in        ! Copy input
adj   = adj_in

err(0) = err(0) + sqrt(sum((x - adj)**2))     ! initial error mean 

do i = 1, num_leads          ! Compute the forecast errors
   do j = 1, obs_freq        ! WARNING: DOESN'T WORK IN IMPERFECT MODELS
      call adv_1step(adj)
      call adv_1step(x)
   end do
   err(i) = err(i) + sqrt(sum((x - adj)**2))
end do

end subroutine adj_forecast_eval



  subroutine adj_forecast_out()
!========================================================================
! subroutine adj_forecast_out()

implicit none

integer :: i

if(count == 0) return

do i = 0, num_leads
   write(*, *) 'adjoint lead ', i, ' error ', real(err(i) / count)
end do

end subroutine adj_forecast_out



  function get_cost(ic, obs, num_obs, n_times, obs_freq)
!========================================================================
! function get_cost(ic, obs, num_obs, n_times, obs_freq)
!
! Computes cost functional for fit starting with ic at first obs time
! Have a total number of n_times that obs were taken and they were taken
! every obs_freq model steps
!

implicit none

integer,  intent(in) :: num_obs, n_times, obs_freq
real(r8), intent(in) :: ic(:), obs(num_obs, n_times)
real(r8)             :: get_cost

real(r8) :: x(size(ic))
integer  :: i

x = ic
get_cost = sqrt(sum((take_obs(x) - obs(:, 1))**2))

do i = 1, n_times - 1
   call advance(x, obs_freq, x)
   get_cost = get_cost + sqrt(sum((take_obs(x) - obs(:, i + 1))**2))
end do

end function get_cost


  subroutine one_d_min(x, grad, obs, num_obs, n_times, obs_freq, xnew)
!==========================================================================
! subroutine one_d_min(x, grad, obs, num_obs, n_times, obs_freq, xnew)
!
! Starting with ic x, searches in direction grad in a stupid way for a minimum

implicit none

integer,  intent(in)  :: num_obs, n_times, obs_freq
real(r8), intent(in)  :: x(:), grad(:), obs(num_obs, n_times)
real(r8), intent(out) :: xnew(:)

real(r8) :: cost, new_cost
integer  :: i

cost = get_cost(x, obs, num_obs, n_times, obs_freq)     ! initial cost
write(*, *) 'iniital cost in one_d_min ', cost

do i = 1, max_1d
   xnew = x - 2**(i - 1) * small_dist * grad
   new_cost = get_cost(xnew, obs, num_obs, n_times, obs_freq)
   write(*, *) 'cost ', i, ' in one_d_min ', new_cost
   if(new_cost > cost) then
      if(i < 3) then
         write(*, *) 'ERROR(one_d_min): no min found'
         stop
      endif
      xnew = x - 2**(i - 2) * small_dist * grad
      call refine_one_d(x - 2**(i - 3) * small_dist * grad, &
         x - 2**(i - 1) * small_dist * grad, obs, num_obs, n_times, obs_freq, xnew)
      return
   end if
   cost = new_cost
end do   

write(*, *) 'ERROR(one_d_min): fell off the end'
stop

end subroutine one_d_min



  subroutine refine_one_d(x1_in, x2_in, obs, num_obs, n_times, obs_freq, xnew)
!==========================================================================
! subroutine refine_one_d(x1_in, x2_in, obs, num_obs, n_times, obs_freq, xnew)

implicit none

integer,  intent(in)  :: num_obs, n_times, obs_freq
real(r8), intent(in)  :: x1_in(:), x2_in(:), obs(num_obs, n_times)
real(r8), intent(out) :: xnew(:)

real(r8) :: cost1, cost2, costm
real(r8) :: xm(size(x1_in)), x1(size(x1_in)), x2(size(x1_in))
integer  :: i

x1    = x1_in
cost1 = get_cost(x1, obs, num_obs, n_times, obs_freq)
x2    = x2_in
cost2 = get_cost(x2, obs, num_obs, n_times, obs_freq)

! Know that min is between x1 and x2

do i = 1, num_refine
   xm = (x1 + x2) / 2.0
   costm = get_cost(xm, obs, num_obs, n_times, obs_freq)
   write(*, *) 'in refine ', i, ' cost is ', costm
   if(cost1 < cost2) then
      cost2 = costm
      x2 = xm
   else
      cost1 = costm
      x1 = xm
   end if
end do

end subroutine refine_one_d



  subroutine get_grad(x, obs, num_obs, n_times, obs_freq, grad)
!==========================================================================
! subroutine get_grad(x, obs, num_obs, n_times, obs_freq, grad)

implicit none

integer,  intent(in)  :: num_obs, n_times, obs_freq
real(r8), intent(in)  :: x(:), obs(num_obs, n_times)
real(r8), intent(out) :: grad(:)

real(r8) :: base_cost, xp(size(x)), xm(size(x)), cp, cm
integer  :: i

! For now do silly delta perturbations in each state direction to get gradient
! Could use adjoint for efficiency

!base_cost = get_cost (x, obs, num_obs, n_times, obs_freq)

do i = 1, size(x)
   xp = x
   xp(i) = xp(i) + delta
   cp = get_cost(xp, obs, num_obs, n_times, obs_freq)
   xm = x
   xm(i) = xm(i) - delta
   cm = get_cost(xm, obs, num_obs, n_times, obs_freq)
   grad(i) = (cp - cm) / (2.0_d0 * delta)
end do

end subroutine get_grad   


  subroutine funct2(mode, n, x_in, objf_out, objgrd_out, nstate, iuser, user)
!==========================================================================
! subroutine funct2(mode, n, x_in, objf_out, objgrd_out, nstate, iuser, user)
!
! TJH ... I had to guess at some of the "intents" ...
! TJH ... "user" not used ... as far as I can tell ...
!

implicit none

integer, intent(in)    :: n, mode, nstate, iuser(:)
real,    intent(inout) :: x_in(n)
real,    intent(out)   :: objf_out, objgrd_out(n)
real                   :: user(:)                 ! not used????

real(r8) :: x(n), objf, objgrd(n)

x = x_in

objf = get_cost(x, g_obs, g_num_obs, g_n_times, g_obs_freq)
if(mode /= 0) call get_grad(x, g_obs, g_num_obs, g_n_times, g_obs_freq, objgrd)

x_in = real(x)
objf_out = real(objf)
objgrd_out = real(objgrd)

end subroutine funct2

!==========================================================================

end module adjoint_mod

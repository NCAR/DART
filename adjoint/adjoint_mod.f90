module adjoint_mod
!
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

use obs_mod, only : take_obs

use model_mod, only : advance, adv_1step

use nag_wrap_mod, only : e04dgf_wrap

private

public min_cost, adj_forecast_eval, adj_forecast_out

double precision, parameter :: small_dist = 1e-3 , delta = 1e-10
integer, parameter :: max_1d = 100, max_iters = 100, num_refine = 20

! Parameters and storage for forecast evaluation
integer, parameter :: num_leads = 20
integer :: count = 0
double precision :: err(0:num_leads) = 0.0

integer :: g_num_obs, g_n_times, g_obs_freq
double precision, allocatable :: g_obs(:, :)


contains

!========================================================================

subroutine min_cost(x, obs, num_obs, n_times, obs_freq, xnew)

implicit none

integer, intent(in) :: num_obs, n_times, obs_freq
double precision, intent(in) :: x(:), obs(num_obs, n_times)
! xnew is value of best fit at LAST STEP (assimilated value)
double precision, intent(out) :: xnew(:)

double precision :: cost, grad(size(x)), xp(size(x)), xmin(size(x)), objf
double precision :: objgrd(size(x)), work(13*size(x)), user(1)
integer :: iwork(size(x) + 1), iuser(1), ifail, i, j, iter

!if(1 == 1) goto 20

ifail = -1
g_num_obs = num_obs
g_n_times = n_times
g_obs_freq = obs_freq
allocate(g_obs(num_obs, n_times))
g_obs = obs
xmin = x
write(*, *) 'xmin before eo4 ', xmin
call e04dgf_wrap(size(x), funct2, iter, objf, objgrd, xmin, iwork, work, iuser, &
   user, ifail)
write(*, *) 'xmin after eo4 ', xmin
xnew = xmin
deallocate(g_obs)

if(1 == 1) goto 10



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

! Want to return value at last step
10 do i = 1, n_times - 1
   write(*, *) 'advancing ', i
   do j = 1, size(xnew)
      write(*, *) 'xnew ', j, xnew(j)
   end do
   call advance(xnew, obs_freq, xnew)
   write(*, *) 'done advancing ', i
end do 
   
end subroutine min_cost

!========================================================================

subroutine adj_forecast_eval(x_in, adj_in, obs_freq)

implicit none

double precision, intent(in) :: x_in(:), adj_in(:)
integer, intent(in) :: obs_freq

double precision :: x(size(x_in)), adj(size(adj_in))
integer :: i, j

count = count + 1

! Copy input
x = x_in
adj = adj_in

! Compute initial error mean 
err(0) = err(0) + sqrt(sum((x - adj)**2))

! Compute the forecast errors
! WARNING: DOESN'T WORK IN IMPERFECT MODELS
do i = 1, num_leads
   do j = 1, obs_freq
      call adv_1step(adj)
      call adv_1step(x)
   end do
   err(i) = err(i) + sqrt(sum((x - adj)**2))
end do

end subroutine adj_forecast_eval

!========================================================================

subroutine adj_forecast_out()

implicit none

integer :: i

if(count == 0) return

do i = 0, num_leads
   write(*, *) 'adjoint lead ', i, ' error ', real(err(i) / count)
end do

end subroutine adj_forecast_out

!========================================================================

function get_cost(ic, obs, num_obs, n_times, obs_freq)

implicit none

double precision :: get_cost
integer, intent(in) :: num_obs, n_times, obs_freq
double precision, intent(in) :: ic(:), obs(num_obs, n_times)

double precision :: x(size(ic))
integer :: i

! Computes cost functional for fit starting with ic at first obs time
! Have a total number of n_times that obs were taken and they were taken
! every obs_freq model steps

x = ic
get_cost = sqrt(sum((take_obs(x) - obs(:, 1))**2))

do i = 1, n_times - 1
   call advance(x, obs_freq, x)
   get_cost = get_cost + sqrt(sum((take_obs(x) - obs(:, i + 1))**2))
end do

end function get_cost

!==========================================================================

subroutine one_d_min(x, grad, obs, num_obs, n_times, obs_freq, xnew)

implicit none

! Starting with ic x, searches in direction grad in a stupid way for a minimum
integer, intent(in) :: num_obs, n_times, obs_freq
double precision, intent(in) :: x(:), grad(:), obs(num_obs, n_times)
double precision, intent(out) :: xnew(:)

double precision :: cost, new_cost
integer :: i

! Compute initial cost
cost = get_cost(x, obs, num_obs, n_times, obs_freq)
write(*, *) 'iniital cost in one_d_min ', cost

do i = 1, max_1d
   xnew = x - 2**(i - 1) * small_dist * grad
   new_cost = get_cost(xnew, obs, num_obs, n_times, obs_freq)
   write(*, *) 'cost ', i, ' in one_d_min ', new_cost
   if(new_cost > cost) then
      if(i < 3) then
         write(*, *) 'ERROR: NO MIN FOUND in one_d_min'
         stop
      endif
      xnew = x - 2**(i - 2) * small_dist * grad
      call refine_one_d(x - 2**(i - 3) * small_dist * grad, &
         x - 2**(i - 1) * small_dist * grad, obs, num_obs, n_times, obs_freq, xnew)
      return
   end if
   cost = new_cost
end do   

write(*, *) 'error: fell off the end in one_d_min'
stop

end subroutine one_d_min

!==========================================================================

subroutine refine_one_d(x1_in, x2_in, obs, num_obs, n_times, obs_freq, xnew)

implicit none

integer, intent(in) :: num_obs, n_times, obs_freq
double precision, intent(in) :: x1_in(:), x2_in(:), obs(num_obs, n_times)
double precision, intent(out) :: xnew(:)

double precision :: cost1, cost2, costm, xm(size(x1_in)), x1(size(x1_in)), x2(size(x1_in))
integer :: i

x1 = x1_in
cost1 = get_cost(x1, obs, num_obs, n_times, obs_freq)
x2 = x2_in
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

!==========================================================================

subroutine get_grad(x, obs, num_obs, n_times, obs_freq, grad)

implicit none

integer, intent(in) :: num_obs, n_times, obs_freq
double precision, intent(in) :: x(:), obs(num_obs, n_times)
double precision, intent(out) :: grad(:)

double precision :: base_cost, xp(size(x)), xm(size(x)), cp, cm
integer i

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
   grad(i) = (cp - cm) / (2.0 * delta)
end do

end subroutine get_grad   

!==========================================================================

subroutine funct2(mode, n, x_in, objf_out, objgrd_out, nstate, iuser, user)

implicit none

integer :: n, mode, nstate, iuser(:)
real :: x_in(n), objf_out, objgrd_out(n), user(:)
double precision x(n), objf, objgrd(n)

x = dble(x_in)

objf = get_cost(x, g_obs, g_num_obs, g_n_times, g_obs_freq)
if(mode /= 0) call get_grad(x, g_obs, g_num_obs, g_n_times, g_obs_freq, objgrd)

x_in = real(x)
objf_out = real(objf)
objgrd_out = real(objgrd)

end subroutine funct2

!==========================================================================

end module adjoint_mod

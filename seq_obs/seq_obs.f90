program seq_obs


! CREATED 9 April, 2001: Uses update of observation only followed by 
! correlative update of each associated state variable. This allows
! arbitrarily fancy updates of the obs state without having to mess with
! two dimensions or anything else. Also should be faster since
! the update computations need not be done.

! This comment to test CVS checkin, 26 Feb. 2002

!----------------------------------------------------------------------

use assim_tools_mod, only : add_noise, read_restart, write_restart, &
   assim_tools_init, obs_increment, update_from_obs_inc, &
   obs_increment2, obs_increment3, obs_increment4, local_update_from_obs_inc, &
   robust_update_from_obs_inc
use model_mod, only : init_model, get_model_size, init_conditions, adv_1step, &
   adv_true_state, output, diag_output_index, model_output, state_loc 
!, balance_init
!   dp_to_grid
use obs_mod, only : num_obs, take_obs, take_single_obs, ens_ics, &
   obs_var, init_obs, get_close_state, obs_loc
use assim_diag_mod, only : diag_type, bin_num, output_diagnostics, &
   save_diagnostics, assim_diag_init
use utilities_mod, only : file_exist, open_file, check_nml_error, &
   print_version_number, close_file
!use forecast_eval_mod, only : forecast_advance, forecast_out
use sort_mod, only : sort, index_sort
use loc_and_dist_mod, only : loc_type, get_dist
use cov_cutoff_mod, only : comp_cov_factor
use correct_corr_mod, only : corr_update_from_obs_inc, corr_obs_increment

!-----------------------------------------------------------------------

implicit none


!=======================================================================

!---- namelist with default values

integer :: ens_size = 0
integer :: spin_up_steps = 0
integer :: num_steps = 0
integer :: output_start = 0
integer :: ens_spin_up_steps = 0
integer :: obs_freq = 0
double precision :: cov_inflate = 1.0
double precision :: mean_inflate = 1.0
double precision :: noise_amp = 0.0

logical :: restart = .FALSE.
logical :: add_on= .FALSE.

namelist /lim_obs_seq_nml/ ens_size, spin_up_steps, num_steps, output_start, &
   ens_spin_up_steps, obs_freq, cov_inflate, mean_inflate, restart, add_on, noise_amp

!--- module name and version number ----
character(len = 11), parameter :: module_name = 'seq_obs'
character(len = 5), parameter :: vers_num = 'x11.0'

!=======================================================================

type (diag_type) :: diag, pre_diag

character*4 :: file_name
character*40 :: diag_file

integer :: i, j, k, num, ifail, start_step, last_step, end_step, index, jjj
integer :: ind_unit(9), obs_unit(9), diag_unit, lji, model_size
integer :: unit, ierr, io

!---------------------------------------------------------------------------

double precision, allocatable :: x(:), ens(:, :), ens_mean(:), ens_obs(:)
double precision, allocatable :: obs_inc(:), wt(:), ens_inc(:), spread_inc(:)
double precision, allocatable :: obs(:), obs_variance(:)
double precision :: sx, s_x2, err_var, mean_inc

! Storage for getting indices of nearest state variables for obs
integer, allocatable :: close_state_list(:), obs_index(:)
integer :: num_close

!--------------------------------------------------------------------------

! Storage for growth rate inflation test
double precision :: dif, rate

! Storage for ramped covariance tests
double precision :: cov_factor, dist, cutoff

!--------------------------------------------------------------------------


! Read namelist for run time control
if(file_exist('input.nml')) then
   unit = open_file(file = 'input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(unit, nml = lim_obs_seq_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'lim_obs_seq_nml')
   enddo
 11 continue
   call close_file(unit)
endif

! Write the namelist to a log file
unit = open_file(file = 'logfile.out', action = 'write')
call print_version_number(unit, module_name, vers_num)
write(unit, nml = lim_obs_seq_nml)
call close_file(unit)

write(*, *) 'namelist inputs are'
write(*, *) 'ens_size ', ens_size
write(*, *) 'spin_up_steps ', spin_up_steps
write(*, *) 'num_steps ', num_steps
write(*, *) 'output_start ', output_start
write(*, *) 'ens_spin_up_steps ', ens_spin_up_steps
write(*, *) 'obs_freq ', obs_freq
write(*, *) 'cov_inflate ', cov_inflate
write(*, *) 'mean_inflate ', mean_inflate
write(*, *) 'restart ', restart
write(*, *) 'add_on ', add_on
write(*, *) 'noise_amp ', noise_amp

!--------------------------------------------------------------------------

! Initialize the model 
call init_model()

! Initialize assim tools module
call assim_tools_init()

! Get model size and allocate model size dependent storage
model_size = get_model_size()

allocate(x(model_size), ens(model_size, ens_size), ens_mean(model_size), &
   close_state_list(model_size), ens_inc(ens_size), &
   ens_obs(ens_size), obs_inc(ens_size), wt(ens_size), obs_index(ens_size), &
   spread_inc(ens_size))

! SNZ: You should read in your conditions here, not by adding extra write.
call init_conditions(x)

! Begin by initializing the observations
call init_obs()
! Allocate storage that depends on the number of observations
allocate(obs(num_obs), obs_variance(num_obs))

write(*, *) 'obs_freq is ', obs_freq
write(*, *) 'model size is ', model_size, 'ensemble size is ', ens_size
write(*, *) 'mean_inflate is ', mean_inflate, 'cov_inflate is ', cov_inflate

! Initialize output diagnostics data structures
call assim_diag_init(diag, model_size, ens_size, diag_output_index, &
   size(diag_output_index), add_on)
call assim_diag_init(pre_diag, model_size, ens_size, diag_output_index, &
   size(diag_output_index), add_on)

!  OPEN UNITS FOR OUTPUT OF INDIVIDUAL MODEL COMPONENTS
do i = 1, min(model_size, 9)
   write(file_name, 61) 'out', i
 61 format(a3, i1)
   write(*, *) 'file_name is ', file_name, ' unit is ', ind_unit(i)
! Get an available unit for this file
   if(add_on) then
      ind_unit(i) = open_file(file = file_name, action = 'append')
   else
      ind_unit(i) = open_file(file = file_name, action = 'write')
   end if
   write(file_name, 61) 'obs', i
   write(*, *) 'file_name is ', file_name, ' unit is ', obs_unit(i)
   if(add_on) then
      obs_unit(i) = open_file(file = file_name, action = 'append')
   else
      obs_unit(i) = open_file(file = file_name, action = 'write')
   end if
end do

! Get a file for diagnostic output
write(diag_file, 31) ens_size, obs_freq, cov_inflate
31 format('ens', i3.3, '_obsint', i3.3, '_inf', f5.3)

if(restart) then
   call read_restart(x, ens, last_step)
   start_step = last_step + 1
   end_step = start_step + num_steps - 1
else

! Let model go for a while to get out onto attractor
   write(*, *) 'doing advance of ', spin_up_steps, ' steps'
   do i = 1, spin_up_steps
      call adv_true_state(x)
   end do

   start_step = 1
   end_step = num_steps

! Set up the first guess for assimilated ensemble
   call ens_ics(x, ens)

! Advance control and ensemble a little ways
   do i = 1, ens_spin_up_steps
      write(*, *) 'ensemble spinup steps is ', i
      call adv_true_state(x)
      do j = 1, ens_size
         call adv_1step(ens(:, j))
      end do
   end do
endif

!  Main loop to advance the perfect model and assimilation in time
do i = start_step, end_step
!   write(*, *) 'advancing 1 step ', i
   call adv_true_state(x)
   do j = 1, ens_size
      call adv_1step(ens(:, j))
   end do

!  Form a new sample if obs
   if(i / obs_freq * obs_freq == i) then
      write(*, *) 'doing new obs i = ', i

      obs = add_noise(take_obs(x))

! Compute the ensemble mean of assimilation
      ens_mean = sum(ens, dim=2) / ens_size

! Temporary barot output and animation blocks could go here (see below) 

! Do pre-assimilation diagnostics
      if(i >= output_start) then
         call save_diagnostics(ens, ens_mean, x, pre_diag)
! Write out the observations for individual components
         do j = 1, min(model_size, 9)
            write(obs_unit(j), 41) i, obs(diag_output_index(j))
         end do
      endif
 41   format(i8, 1x, e12.6)

! Compute and print the distance between the ensemble mean and truth
      write(*, *) 'scaled distance to truth is ', real(sqrt(sum((x - ens_mean) * (x - ens_mean))))

! get observational variance (assumed diagnonal for now)
         obs_variance = obs_var()

! Covariance inflate put into prior ensemble estimate here
      do j = 1, ens_size
         ens(:, j) = sqrt(cov_inflate) * (ens(:, j) - ens_mean) + ens_mean
      end do

! Bounding covariance inflate tests
!      do k = 1, model_size
!         dif = sum(dabs(ens(k, :) - ens_mean(k))) / ens_size
!!!         if(dif < 5e5) then
!!!            rate = 1.0 + (sqrt(cov_inflate) - 1.0) * ((5e5 - dif)/ 5e5)
!         if(dif < 8) then
!            rate = 1.0 + (sqrt(cov_inflate) - 1.0) * ((8 - dif)/ 8)
!         else
!            rate = 1.0
!         end if
!         write(*, *) k, dif, rate
!         do j = 1, ens_size
!            ens(k, j) = rate * (ens(k, j) - ens_mean(k)) + ens_mean(k)
!         end do
!      end do



!------------------------------------------------------------------------

! Loop through each observation available at this time
      do lji = 1, num_obs
! Negative obs variance means to not use this ob
         if(obs_variance(lji) > 0.0) then
            call get_close_state(lji, close_state_list, model_size, num_close)

! Get ensemble for this observation
            do j = 1, ens_size
               ens_obs(j) = take_single_obs(ens(:, j), lji)
            end do


! Do an index sort of the obs value for use in local regression
! Take out when not needed, this is expensive
!!!            call index_sort(ens_obs, obs_index, ens_size)



! Get the update increments for this observed variable
            call obs_increment(ens_obs, ens_size, obs(lji), obs_variance(lji), &
               obs_inc)
!!!            call corr_obs_increment(ens_obs, ens_size, obs(lji), obs_variance(lji), &
!!!               mean_inc, spread_inc)

! Section to do adjustment point by point
            do k = 1, num_close
               j = close_state_list(k)

! ADDITION OF DISTANCE WEIGHTING COVARIANCE

               dist = get_dist(state_loc(j), obs_loc(lji))
!              write(*, *) 'state ', j, 'obs ', lji, real(dist)
               cutoff = 0.20
               cov_factor = comp_cov_factor(dist, cutoff)



! Temporary kluge to avoid the distance dependent covariance envelope
!              if(cov_factor > 0.0) cov_factor = 1.0
!            if(dist < 0.025) then

               cov_factor = 1.0

!            else if(dist > 0.1) then 
!               cov_factor = 0.0
!            else 
!               cov_factor = 1.0 - ((dist - 0.025) / 0.075)
!            endif
!             if(cov_factor > 0.0) write(*, *) 'cov factor ', k, dist, cov_factor




               if(1 == 1) then
!!!               if(cov_factor > 0.0) then
!                  write(*, *) 'ob state ', lji, j
! The new update method
                  call update_from_obs_inc(ens_obs, obs_inc, ens(j, :), &
                     ens_size, ens_inc, cov_factor)

!!!                  call corr_update_from_obs_inc(ens_obs, mean_inc, spread_inc, ens(j, :), &
!!!                     ens_size, ens_inc, cov_factor)
!                  do jjj = 1, ens_size
!                     write(*, *) jjj, ens_inc(jjj)
!                  end do
! Test of local version;REMEMBER TO COMMENT OUT SORTING FOR THIS ABOVE
!!!                 call local_update_from_obs_inc(ens_obs, obs_inc, ens(j, :), &
!!!                     ens_size, ens_inc, cov_factor, 9, obs_index)
! Test of robust regression methods; probably not applicable except for speed
!!!                 call robust_update_from_obs_inc(ens_obs, obs_inc, ens(j, :), &
!!!                     ens_size, ens_inc, cov_factor, obs_index)
!                  do jjj = 1, ens_size
!                     write(*, *) jjj, ens_inc(jjj)
!                  end do

! NEW TWO PART OBS UPDATE!!!
                  ens(j, :) = ens(j, :) + ens_inc

               endif
            end do
         endif

      end do

! TEMPORARY INCLUSION OF BALANCE CALC
!      do j = 1, ens_size
!         call balance_init(ens(:, j), ens(:, j))
!      end do

! Compute ensemble mean of current assimilated state
   ens_mean = sum(ens, dim=2) / ens_size

! Do diagnostics
      if(i >= output_start) then
         call save_diagnostics(ens, ens_mean, x, diag)
! First stab at evaluating some forecasts
!         call forecast_advance(x, ens, i, diag_output_index, &
!            size(diag_output_index), add_on)
      endif

      write(*, *) 'scaled distance to truth after is ', real(sqrt(sum((x - ens_mean) * (x - ens_mean))))   

   endif           ! End of computations done at observation times only

! Compute ensemble mean of current assimilated state; NEED TO MIN ENS_MEAN CALLS
   ens_mean = sum(ens, dim=2) / ens_size

   ! Output the first 9 variables

   if(i >= output_start) then

      do j = 1, min(model_size, 9)

         index = diag_output_index(j)

         ! Temporary computation of 'error variance' to compare to Evensen

         sx = sum(ens(index, :))
	 s_x2 = sum(ens(index, :) * ens(index, :))
	 err_var = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

         write(ind_unit(j), 21) i, x(index), ens_mean(index), err_var, ens(index, 1:min(ens_size,10))

      end do

   end if

 21 format(i8, 1x, 13(e12.6, 1x))

end do

! Final diagnostic output follows
diag_unit = open_file(file = diag_file, action = 'write')
call output_diagnostics(diag, diag_unit)

write(*, *) 'STARTING PRE_ASSIM DIAGNOSTIC OUTPUT'
call output_diagnostics(pre_diag, diag_unit)

!write(*, *) 'STARTING FORECAST DIAGNOSTICS'
!call forecast_out(diag_unit, diag_output_index)

! Writing a final restart
write(*, *) 'writing a restart'
call write_restart(x, ens, end_step)

end program seq_obs

!-----------------------------------------------------

subroutine correl(ens, model_size, ens_size, cor)

implicit none

integer, intent(in) :: model_size, ens_size
double precision, intent(in) :: ens(model_size, ens_size)
double precision, intent(out) :: cor(model_size, model_size)

integer :: i, j
double precision :: x(ens_size), y(ens_size), sx, sy, sx2, sy2, sxy

! Temporary to look at correlations in ensemble
do i = 1, model_size
   do j = 1, model_size
      x = ens(i, :)
      y = ens(j, :)
      sx = sum(x)
      sy = sum(y)
      sx2 = sum(x**2)
      sy2 = sum(y**2)
      sxy = sum(x * y)
      cor(i, j) = (ens_size * sxy - sx * sy) / sqrt((ens_size * sx2 - sx**2) * &
                                                    (ens_size * sy2 - sy**2))
   end do
end do      

end subroutine correl

!===========================================================================================
 
function get_correl(x, y, n)
 
implicit none
 
integer, intent(in) :: n
double precision, intent(in) :: x(n), y(n)
double precision :: get_correl
double precision :: sx, sy, sx2, sy2, sxy
 
sx = sum(x)
sy = sum(y)
sx2 = sum(x**2)
sy2 = sum(y**2)
sxy = sum(x * y)
 
get_correl = (n * sxy - sx * sy) / sqrt((n * sx2 - sx**2) * (n * sy2 - sy**2) )
 
if(dabs(get_correl)  > 1.01) then
 
   write(*, *) 'get_correl is ', real(get_correl)
   write(*, *) 'n is ', n
   write(*, *) 'x is ', x
   write(*, *) 'y is ', y
endif
 
end function get_correl
 
!===========================================================================================


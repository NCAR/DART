! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module assim_tools_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! A variety of operations required by assimilation.

use types_mod, only :: r8
!!!use obs_mod,        only : num_obs, obs_var
use utilities_mod,  only : file_exist, open_file, check_nml_error, &
                           close_file
!                           print_version_number, close_file
use sort_mod,       only : index_sort

! Added 22 January, 2001 to duplicate observations no matter what else is
! done with random number generator. Allows clear enkf_2d comparisons.

use random_seq_mod, only : random_seq_type, random_gaussian, &
                           init_random_seq, random_uniform

logical :: first_ran_call = .true., first_inc_ran_call = .true.
logical :: first_obs18_call = .true.

type (random_seq_type) :: ran_seq, inc_ran_seq

implicit none
private

public read_restart, write_restart, assim_tools_init, &
   obs_increment, obs_increment2, obs_increment3, obs_increment4, &
   update_from_obs_inc, local_update_from_obs_inc, robust_update_from_obs_inc, &
   obs_inc_index, obs_inc4_index, obs_increment5, obs_increment6, &
   obs_increment7, obs_increment8, obs_increment9, obs_increment10, &
   obs_increment11, obs_increment12, obs_increment13, obs_increment14, &
   obs_increment15, obs_increment16, obs_increment17, obs_increment18, &
   obs_increment_group, obs_increment_particle, &
   linear_obs_increment, linear_update_from_obs_inc, look_for_bias

!============================================================================

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!---- namelist with default values

real(r8) :: cor_cutoff = 0.0_r8

namelist / assim_tools_nml / cor_cutoff

character(len=128), parameter :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

!============================================================================

contains



subroutine assim_tools_init()
!============================================================================
! subroutine assim_tools_init()
!

implicit none

integer :: unit, ierr, io

! Read namelist for run time control

if(file_exist('input.nml')) then
   unit = open_file('input.nml', action = 'read')
   ierr = 1

! TJH Coding Standard does not allow use of the "end=" construct.
!
!  do while(ierr /= 0)
!     read(unit, nml = assim_tools_nml, iostat = io, end = 11)
!     ierr = check_nml_error(io, 'assim_tools_nml')
!  enddo
!11 continue

   READBLOCK: do while(ierr /= 0)
      read(unit, nml = assim_tools_nml, iostat = io)
      if ( io < 0 ) exit READBLOCK          ! end-of-file
      ierr = check_nml_error(io, 'assim_tools_nml')
   enddo READBLOCK
 
   call close_file(unit)
endif

! TJH 14.03.2002 What do we do if the namelist does not exist?

! Write the namelist and module info to a log file

unit = open_file('logfile.out', action = 'append')
write(unit,*)'assim_tools attributes:'
write(unit,*)'   ',trim(adjustl(source))
write(unit,*)'   ',trim(adjustl(revision))
write(unit,*)'   ',trim(adjustl(revdate))
write(unit, nml = assim_tools_nml)
call close_file(unit)

write(*, *) 'Correlation cutoff in fast_seq_non_identity_prod is ', cor_cutoff

end subroutine assim_tools_init




subroutine read_restart(x, ens, last_step)
!============================================================================
! subroutine read_restart(x, ens, last_step)
!

implicit none

real(r8), intent(inout) :: x(:), ens(:, :)
integer,  intent(inout) :: last_step

integer :: chan

chan = open_file('restart.in', action = 'read')
read(chan, *) last_step
read(chan, *) ens
read(chan, *) x
write(*, *) 'successfully read x'
call close_file(chan)

end subroutine read_restart



subroutine write_restart(x, ens, last_step)
!========================================================================
! subroutine write_restart(x, ens, last_step)
!

implicit none

real(r8), intent(in) :: x(:), ens(:, :)
integer,  intent(in) :: last_step

integer :: chan

chan = open_file('restart.out', action = 'write')
write(chan, *) last_step
write(chan, *) ens
write(chan, *) x
call close_file(chan)

end subroutine write_restart



subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc, a)
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
real(r8), intent(out) :: a

real(r8) :: prior_mean, new_mean, prior_var, var_ratio, sum_x

! Compute prior variance and mean from sample
sum_x = sum(ens)
prior_mean = sum_x / ens_size
prior_var = (sum(ens * ens) - sum_x**2 / ens_size) / (ens_size - 1)

var_ratio = obs_var / (prior_var + obs_var)
new_mean = var_ratio * (prior_mean  + prior_var*obs / obs_var)

a = sqrt(var_ratio)

obs_inc = a * (ens - prior_mean) + new_mean - ens

end subroutine obs_increment


subroutine obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc)
!------------------------------------------------------------------------
!
! A observation space only particle filter implementation for a
! two step sequential update filter. Second version, 2 October, 2003.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: weight(ens_size), rel_weight(ens_size), cum_weight(0:ens_size)
real(r8) :: base, frac, new_val(ens_size), weight_sum
integer :: i, j, index(ens_size), ens_index(ens_size), new_index(ens_size)

! Begin by computing a weight for each of the prior ensemble members
do i = 1, ens_size
   weight(i) = exp(-1.0 * (ens(i) - obs)**2 / (2.0 * obs_var))
end do

! Compute relative weight for each ensemble member
weight_sum = sum(weight)
do i = 1, ens_size
   rel_weight(i) = weight(i) / weight_sum
end do

! Compute cumulative weights at boundaries
cum_weight(0) = 0
do i = 1, ens_size
   cum_weight(i) = cum_weight(i - 1) + rel_weight(i)
!   write(*, 131) i, weight(i), rel_weight(i), cum_weight(i)
131 format(1x, i3, 3(e10.4, 1x))
end do
! Fix up for round-off error if any
cum_weight(ens_size) = 1.0

! Do a deterministic implementation: just divide interval into ens_size parts and see
! which interval this is in (careful to offset; not start at 0)
base = 1.0 / (ens_size * 2.0)
do i = 1, ens_size
   frac = base + (i - 1.0) / ens_size
! Now search in the cumulative range to see where this frac falls
! Can make this search more efficient by limiting base
   do j = 1, ens_size
      if(cum_weight(j - 1) < frac .and. frac < cum_weight(j)) then
         index(i) = j
!         write(*, *) i, frac, 'gets index ', j
         goto 111
      end if
   end do
111 continue
end do

! Set the new values for the ensemble members
do i = 1, ens_size
   new_val(i) = ens(index(i))
!   write(*, *) 'new_val ', i, new_val(i)
end do

! Try sorting to make increments as small as possible
call index_sort(ens, ens_index, ens_size)
call index_sort(new_val, new_index, ens_size)
do i = 1, ens_size
   obs_inc(ens_index(i)) = new_val(new_index(i)) - ens(ens_index(i))
end do

end subroutine obs_increment_particle




subroutine obs_increment18(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment

! Working with sampling error correction in obs space from 23 Sept. 2003

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8) :: prior_cov, sx, s_x2
real(r8) :: base_fract, new_fract, final_mean
real(r8) :: ratio(189), factor(189), var_ratio, this_factor

integer :: i, index

! Need to read in files for factor if this is first time through
if(first_obs18_call) then
   first_obs18_call = .false.
   open(unit = 51, file = 'inc18_factor')
   do i = 1, 189
      read(51, *) ratio(i), factor(i)
!      write(*, *) 'reading ', i, ratio(i), factor(i)
   end do
endif

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

! Find the factor for this var ratio
var_ratio = prior_cov / obs_var
if(var_ratio > ratio(189)) then
   this_factor = 1.0
else
   do i = 1, 189
      if(var_ratio < ratio(i)) then
         this_factor = factor(i)
         goto 111
      endif
   end do
endif
   
111 continue
!write(*, *) 'var_ratio, this_factor ', var_ratio, this_factor
!!!if(this_factor < 1.01) this_factor = 1.015

prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

! From work on 22 Sept. 2003 (see sys_sim203) to handle sampling bias
a = sqrt(this_factor * new_cov * prior_cov_inv)

! Also want to increase fraction moved towards ob by same factor
base_fract = (new_mean - prior_mean) / (obs - prior_mean)
new_fract = this_factor * base_fract
final_mean = new_fract * (obs - prior_mean) + prior_mean

obs_inc = a * (ens - prior_mean) + final_mean - ens

end subroutine obs_increment18




subroutine obs_increment_group(ens, ens_size, obs, obs_var, obs_inc, &
   num_groups)
!========================================================================
!
! EAKF version of obs increment for groups of ensembles

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
integer, intent(in) :: num_groups

real(r8) :: a, a1(num_groups), obs_var_inv
real(r8) :: new_mean(num_groups)
real(r8) :: prior_mean(num_groups), prior_cov(num_groups)
real(r8) :: sx, s_x2
real(r8) :: mean_mean, mean_var, inf_prior_cov_inv, inf_prior_cov
real(r8) :: inf_new_cov(num_groups)
real(r8) :: new_ens(ens_size) 

integer :: grp_size, group, grp_bot, grp_top

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute the prior_mean and variance for each group
! Divide ensemble into num_groups groups
grp_size = ens_size / num_groups
do group = 1, num_groups
   grp_bot = (group - 1) * grp_size + 1
   grp_top = grp_bot + grp_size - 1

! Compute prior covariance and mean from sample
   sx = sum(ens(grp_bot:grp_top))
   s_x2 = sum(ens(grp_bot:grp_top) * ens(grp_bot:grp_top))
   prior_mean(group) = sx / grp_size
   prior_cov(group) = (s_x2 - sx**2 / grp_size) / (grp_size - 1)
end do

! Compute the variance of the prior mean estimates between groups
sx = sum(prior_mean)
s_x2 = sum(prior_mean * prior_mean)
mean_mean = sx / num_groups
mean_var = (s_x2 - sx**2 / num_groups) / (num_groups - 1)

!!!write(*, *) 'mean var is ', mean_var

! Can use incremented variance for computing new means only

do group = 1, num_groups
   grp_bot = (group - 1) * grp_size + 1
   grp_top = grp_bot + grp_size - 1

! Look out for next; probably not exactly right
   inf_prior_cov = prior_cov(group) + mean_var
   a1(group) = sqrt(inf_prior_cov / prior_cov(group))
!   write(*, *) 'a1 ', a1(group), 1.0 / a1(group)
! Also need to inflate the ensemble beforehand
   new_ens(grp_bot:grp_top) = a1(group) * (ens(grp_bot:grp_top) - prior_mean(group)) + &
      prior_mean(group)

   inf_prior_cov_inv = 1.0 / inf_prior_cov
   inf_new_cov(group) = 1.0 / (inf_prior_cov_inv + obs_var_inv)
   a = sqrt(inf_new_cov(group) * inf_prior_cov_inv)
!   write(*, *) 'first a ', a

   new_mean(group) = inf_new_cov(group) * (inf_prior_cov_inv * prior_mean(group) + obs / obs_var)
   
   new_ens(grp_bot:grp_top) = a * (new_ens(grp_bot:grp_top) - prior_mean(group)) &
      + new_mean(group)
end do


! Quick, more stable cheat??? Just restore the factor added in
!!!do group = 1, num_groups
!!!   grp_bot = (group - 1) * grp_size + 1
!!!   grp_top = grp_bot + grp_size - 1
!!!   new_ens(grp_bot:grp_top) = (1.0 / a1(group)) * &
!!!      (new_ens(grp_bot:grp_top) - new_mean(group)) + new_mean(group)
!!!   obs_inc(grp_bot:grp_top) = new_ens(grp_bot:grp_top) - ens(grp_bot:grp_top)
!!!end do

!!!if(1 == 1) return


! Now, have to remove that portion of the variance that is held in the varying means
! Compute the variance of the posterior mean estimates between groups
sx = sum(new_mean)
s_x2 = sum(new_mean * new_mean)
mean_mean = sx / num_groups
mean_var = (s_x2 - sx**2 / num_groups) / (num_groups - 1)

! Now, remove this part of variance from the individual groups of ensembles
do group = 1, num_groups
   grp_bot = (group - 1) * grp_size + 1
   grp_top = grp_bot + grp_size - 1
   if(mean_var > 1.00 * inf_new_cov(group)) then
      a = 0.0
   else
      a = sqrt((inf_new_cov(group) - mean_var) / inf_new_cov(group))
   endif
!   write(*, *) 'inf second a ', a
   new_ens(grp_bot:grp_top) =  a * (new_ens(grp_bot:grp_top) - new_mean(group)) &
      + new_mean(group)
   obs_inc(grp_bot:grp_top) = new_ens(grp_bot:grp_top) - ens(grp_bot:grp_top)
end do

end subroutine obs_increment_group




subroutine third_increment_group(ens, ens_size, obs, obs_var, obs_inc, &
   num_groups)
!========================================================================
!
! EAKF version of obs increment for groups of ensembles

! Approach here failed. Be sure to understand why.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
integer, intent(in) :: num_groups

real(r8) :: a(num_groups), prior_mean(num_groups), new_mean(num_groups)
real(r8) :: ens2(1, ens_size), obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_cov_inv, new_cov
real(r8) :: prior_cov, sx, s_x2, new_global_mean, mean_sum
real(r8) :: mean_mean, a_min

integer :: i, grp_size, group, grp_bot, grp_top

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var


! Now do the variance from the individual groups
! Divide ensemble into num_groups groups
grp_size = ens_size / num_groups
do group = 1, num_groups
   grp_bot = (group - 1) * grp_size + 1
   grp_top = grp_bot + grp_size - 1

! Compute prior covariance and mean from sample
   sx = sum(ens(grp_bot:grp_top))
   s_x2 = sum(ens(grp_bot:grp_top) * ens(grp_bot:grp_top))
   prior_mean(group) = sx / grp_size
   prior_cov = (s_x2 - sx**2 / grp_size) / (grp_size - 1)

   prior_cov_inv = 1.0 / prior_cov
   new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

   new_mean(group) = new_cov * (prior_cov_inv * prior_mean(group) + obs / obs_var)

   a(group) = sqrt(new_cov * prior_cov_inv)

end do

write(*, *) 'a are ', a
a_min = maxval(a)
!a_min = sum(a) / num_groups


do group = 1, num_groups
   grp_bot = (group - 1) * grp_size + 1
   grp_top = grp_bot + grp_size - 1

   obs_inc(grp_bot:grp_top) = a_min * (ens(grp_bot:grp_top) - prior_mean(group)) &
      + new_mean(group) - ens(grp_bot:grp_top)

end do

end subroutine third_increment_group




subroutine bad_increment_group(ens, ens_size, obs, obs_var, obs_inc, &
   num_groups)
!========================================================================
!
! EAKF version of obs increment for groups of ensembles

! Approach here failed. Be sure to understand why.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
integer, intent(in) :: num_groups

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8) :: prior_cov, sx, s_x2, new_global_mean, mean_sum
real(r8) :: mean_mean

integer :: i, grp_size, group, grp_bot, grp_top

! Accumulate individual ensemble means for grand mean
mean_sum = 0.0
! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! First compute the new mean from larger group statistics
! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

new_global_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)


! Now do the variance from the individual groups
! Divide ensemble into num_groups groups
grp_size = ens_size / num_groups
do group = 1, num_groups
   grp_bot = (group - 1) * grp_size + 1
   grp_top = grp_bot + grp_size - 1

! Compute prior covariance and mean from sample
   sx = sum(ens(grp_bot:grp_top))
   s_x2 = sum(ens(grp_bot:grp_top) * ens(grp_bot:grp_top))
   prior_mean = sx / grp_size
   prior_cov = (s_x2 - sx**2 / grp_size) / (grp_size - 1)

   prior_cov_inv = 1.0 / prior_cov
   new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

   new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)
   mean_sum = mean_sum + new_mean

   a = sqrt(new_cov * prior_cov_inv)

   obs_inc(grp_bot:grp_top) = a * (ens(grp_bot:grp_top) - prior_mean) &
      + new_mean - ens(grp_bot:grp_top)

end do

mean_mean = mean_sum / num_groups
!!!write(*, *) 'global and mean mean ', new_global_mean, mean_mean

! Add on the difference to all ensemble members
obs_inc = obs_inc + (new_global_mean - mean_mean)

end subroutine bad_increment_group




subroutine obs_increment12(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment
! Look at adding additional variance from systematic error. Add same
! amount to each observation. 11 April, 2003.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2, inf_prior_cov, inf_prior_cov_inv
! WARNING, CLEARLY ONLY WORKS FOR ONLY SURFACE PRESSURE OBS RIGHT NOW
real(r8), parameter :: add_cov = 100.0

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

! Do the mean with the uninflated values???
inf_prior_cov = prior_cov + add_cov
inf_prior_cov_inv = 1.0 / inf_prior_cov
new_cov = 1.0 / (inf_prior_cov_inv + obs_var_inv)
if(new_cov > prior_cov) new_cov = prior_cov
write(*, *) 'prior, new cov ', prior_cov, new_cov

a = sqrt(new_cov * prior_cov_inv)

obs_inc = a * (ens - prior_mean) + new_mean - ens

end subroutine obs_increment12



subroutine obs_increment11(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment

! Ideas from 9 April, 2003. Single global tuning factor that reduces
! amount by which covariance is decreased by constant global factor.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8), parameter :: red_factor = 0.1

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

! TEMPORARY LOOK AT INFLATING HERE; see notes from 12 Sept. 2001
!!!cov = cov * 1.01 OR prior_cov = prior_cov * 1.01

prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

a = sqrt(new_cov * prior_cov_inv)
write(*, *) 'red_factor,  original a', red_factor, a

! Adjust a to reduce it
a = 1.0 - red_factor * (1.0 - a)
write(*, *) '1 - a ', 1.0 - a
write(*, *) 'dif term ', red_factor * (1.0 - a)
write(*, *) 'reduced a', a

obs_inc = a * (ens - prior_mean) + new_mean - ens

end subroutine obs_increment11





subroutine obs_increment14(ens, ens_size, obs, obs_var, obs_inc)
! in its present form.
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment
! ent14 explores changing the shape of the relation between the ratio
! of the error and the expected error_sd and the factor used to increase
! the spread of the prior and obs variance.
! 
! Extending initial success with obs_increment7. In this case, begin 
! by doing the standard update for the mean for the observation. Then
! go back and evaluate the ratio of the prior difference between the
! ensemble mean and the obs and the expected value given the prior
! sample variance and the specified obs variance. If the ratio
! exceeds a threshold, then adjust the prior variance and obs variance
! by the same factor to reduce the ratio to the threshold (could also
! consider damping it toward the threshold by some factor to make
! the shocks less?). Then, recompute the updated variance using this
! and adjust the ensemble spread only. This means that a bad observation
! will increase the spread some but will not immediately pull the mean
! off toward the bad observation any further.

! ORIGINALLY identical to obs_increment7. Looking at changing inflate
! ratio for state versus obs.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8) :: error, diff_sd, ratio, inf_obs_var, inf_ens(ens_size)
real(r8) :: factor
real(r8) :: inf_obs_var_inv, inf_prior_cov, inf_prior_cov_inv
real(r8), parameter :: ratio_threshold = 2.0

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

! Compute the updated mean and covariance in the standard fashion
prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)
new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

! AT THIS POINT CAN EVALUATE INCONSISTENCY
! CHECK ALL THIS
error = prior_mean - obs
diff_sd = sqrt(obs_var + prior_cov)
ratio = abs(error / diff_sd)

! Only modify if the ratio exceeds the threshold value
!!!if(ratio < 1.0) factor = 1.0
!!!if(ratio > 1.0) factor = 1.0 + (ratio - 1.0) / 2
factor = 1.0 + (ratio - 1.0) / 4.0

! Can now inflate by this ratio and then do adjustment
inf_obs_var = factor**2 * obs_var
inf_obs_var_inv = 1.0 / inf_obs_var
! Form inflated ensemble
inf_ens = prior_mean + factor * (ens - prior_mean)
inf_prior_cov = factor**2 * prior_cov


inf_prior_cov_inv = 1.0 / inf_prior_cov
new_cov = 1.0 / (inf_prior_cov_inv + inf_obs_var_inv)

a = sqrt(new_cov * inf_prior_cov_inv)

obs_inc = a * (inf_ens - prior_mean) + new_mean - ens

end subroutine obs_increment14





subroutine obs_increment15(ens, ens_size, obs, obs_var, obs_inc)
! in its present form.
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment
! ent15 explores changing the shape of the relation between the ratio
! of the error and the expected error_sd and the factor used to increase
! the spread of the prior and obs variance. As per exp14, it uses
! an inflation factor of 1.0 + 0.5 * (ratio - 1) BUT, it only applies
! this if the ratio is > 1. Otherwise, factor is left at 1.0.
! 
! Extending initial success with obs_increment7. In this case, begin 
! by doing the standard update for the mean for the observation. Then
! go back and evaluate the ratio of the prior difference between the
! ensemble mean and the obs and the expected value given the prior
! sample variance and the specified obs variance. If the ratio
! exceeds a threshold, then adjust the prior variance and obs variance
! by the same factor to reduce the ratio to the threshold (could also
! consider damping it toward the threshold by some factor to make
! the shocks less?). Then, recompute the updated variance using this
! and adjust the ensemble spread only. This means that a bad observation
! will increase the spread some but will not immediately pull the mean
! off toward the bad observation any further.

! ORIGINALLY identical to obs_increment7. Looking at changing inflate
! ratio for state versus obs.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8) :: error, diff_sd, ratio, inf_obs_var, inf_ens(ens_size)
real(r8) :: factor
real(r8) :: inf_obs_var_inv, inf_prior_cov, inf_prior_cov_inv
real(r8), parameter :: ratio_threshold = 2.0

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

! Compute the updated mean and covariance in the standard fashion
prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)
new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

! AT THIS POINT CAN EVALUATE INCONSISTENCY
! CHECK ALL THIS
error = prior_mean - obs
diff_sd = sqrt(obs_var + prior_cov)
ratio = abs(error / diff_sd)

! Only modify if the ratio exceeds the threshold value
if(ratio > 1.0) then
   factor = 1.0 + (ratio - 1.0) / 2
else
   factor = 1.0
endif

! Can now inflate by this ratio and then do adjustment
inf_obs_var = factor**2 * obs_var
inf_obs_var_inv = 1.0 / inf_obs_var
! Form inflated ensemble
inf_ens = prior_mean + factor * (ens - prior_mean)
inf_prior_cov = factor**2 * prior_cov


inf_prior_cov_inv = 1.0 / inf_prior_cov
new_cov = 1.0 / (inf_prior_cov_inv + inf_obs_var_inv)

a = sqrt(new_cov * inf_prior_cov_inv)

obs_inc = a * (inf_ens - prior_mean) + new_mean - ens

end subroutine obs_increment15





subroutine obs_increment16(ens, ens_size, obs, obs_var, obs_inc)
! in its present form.
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment
! ent15 explores changing the shape of the relation between the ratio
! of the error and the expected error_sd and the factor used to increase
! the spread of the prior and obs variance. As per exp14, it uses
! an inflation factor of 1.0 + 0.5 * (ratio - 1) BUT, it only applies
! this if the ratio is > 1. Otherwise, factor is left at 1.0.
! 
! Extending initial success with obs_increment7. In this case, begin 
! by doing the standard update for the mean for the observation. Then
! go back and evaluate the ratio of the prior difference between the
! ensemble mean and the obs and the expected value given the prior
! sample variance and the specified obs variance. If the ratio
! exceeds a threshold, then adjust the prior variance and obs variance
! by the same factor to reduce the ratio to the threshold (could also
! consider damping it toward the threshold by some factor to make
! the shocks less?). Then, recompute the updated variance using this
! and adjust the ensemble spread only. This means that a bad observation
! will increase the spread some but will not immediately pull the mean
! off toward the bad observation any further.

! ORIGINALLY identical to obs_increment7. Looking at changing inflate
! ratio for state versus obs.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8) :: error, diff_sd, ratio, inf_obs_var, inf_ens(ens_size)
real(r8) :: factor
real(r8) :: inf_obs_var_inv, inf_prior_cov, inf_prior_cov_inv
real(r8), parameter :: ratio_threshold = 2.0

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

! Compute the updated mean and covariance in the standard fashion
prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)
new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

! AT THIS POINT CAN EVALUATE INCONSISTENCY
! CHECK ALL THIS
error = prior_mean - obs
diff_sd = sqrt(obs_var + prior_cov)
ratio = abs(error / diff_sd)

! Only modify if the ratio exceeds the threshold value
if(ratio > 1.0) then
   factor = 1.0 + (ratio - 1.0) / 2.0 
else if (ratio < 1.0) then
   factor = 1.0 + (ratio - 1.0) / 6
endif

! Can now inflate by this ratio and then do adjustment
inf_obs_var = factor**2 * obs_var
inf_obs_var_inv = 1.0 / inf_obs_var
! Form inflated ensemble
inf_ens = prior_mean + factor * (ens - prior_mean)
inf_prior_cov = factor**2 * prior_cov


inf_prior_cov_inv = 1.0 / inf_prior_cov
new_cov = 1.0 / (inf_prior_cov_inv + inf_obs_var_inv)

a = sqrt(new_cov * inf_prior_cov_inv)

obs_inc = a * (inf_ens - prior_mean) + new_mean - ens

end subroutine obs_increment16






subroutine obs_increment17(ens, ens_size, obs, obs_var, obs_inc, slope, a)
! in its present form.
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment
! ent15 explores changing the shape of the relation between the ratio
! of the error and the expected error_sd and the factor used to increase
! the spread of the prior and obs variance. As per exp14, it uses
! an inflation factor of 1.0 + 0.5 * (ratio - 1) BUT, it only applies
! this if the ratio is > 1. Otherwise, factor is left at 1.0.
! 
! Extending initial success with obs_increment7. In this case, begin 
! by doing the standard update for the mean for the observation. Then
! go back and evaluate the ratio of the prior difference between the
! ensemble mean and the obs and the expected value given the prior
! sample variance and the specified obs variance. If the ratio
! exceeds a threshold, then adjust the prior variance and obs variance
! by the same factor to reduce the ratio to the threshold (could also
! consider damping it toward the threshold by some factor to make
! the shocks less?). Then, recompute the updated variance using this
! and adjust the ensemble spread only. This means that a bad observation
! will increase the spread some but will not immediately pull the mean
! off toward the bad observation any further.

! ORIGINALLY identical to obs_increment7. Looking at changing inflate
! ratio for state versus obs.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
real(r8), intent(in) :: slope
real(r8), intent(out) :: a

real(r8) :: obs_var_inv
real(r8) :: prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8) :: error, diff_sd, ratio, inf_obs_var, inf_ens(ens_size)
real(r8) :: factor, b, c
real(r8) :: inf_obs_var_inv, inf_prior_cov, inf_prior_cov_inv

real(r8), parameter :: threshold = 1.0
integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

! Compute the updated mean and covariance in the standard fashion
prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)
new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

! THIS POINT CAN EVALUATE INCONSISTENCY
error = prior_mean - obs
diff_sd = sqrt(obs_var + prior_cov)
ratio = abs(error / diff_sd)

! FIRST PASS GROSS OBS ERROR CHECK;SIMPLEST POSSIBLE
if(ratio > 2.5) then
   write(*, *) 'DISCARDING OBS ____________________________'
   obs_inc = 0.0 
   a = 0.0
   return
endif

! Only modify if the ratio exceeds the threshold value
if(ratio > threshold .and. slope > 0.0) then
   b = (1.0 / slope) ** (1.0 / (slope - 1.0))
   c = -1.0 * b**slope
   factor = ratio / ((ratio - threshold + b)**slope + c + threshold)
!!!   factor = ratio / (0.5 + (ratio - 0.75)**0.5)


!!!   factor = 1.0 + (ratio - 1.0)**0.2 * slope
else 
   factor = 1.0
endif

!!!if(factor > 1.0) write(*, 11) ratio, factor, ratio/factor
!!!11 format(1x, 'ratio, factor, adj ', 3(e10.4, 1x))

! Can now inflate by this ratio and then do adjustment
inf_obs_var = factor**2 * obs_var
inf_obs_var_inv = 1.0 / inf_obs_var
! Form inflated ensemble
inf_ens = prior_mean + factor * (ens - prior_mean)
inf_prior_cov = factor**2 * prior_cov

inf_prior_cov_inv = 1.0 / inf_prior_cov
new_cov = 1.0 / (inf_prior_cov_inv + inf_obs_var_inv)

a = sqrt(new_cov * inf_prior_cov_inv)
!write(*, *) 'a in 17 is ', a

obs_inc = a * (inf_ens - prior_mean) + new_mean - ens

end subroutine obs_increment17




subroutine obs_increment13(ens, ens_size, obs, obs_var, obs_inc)
! in its present form.
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment
! 
! Extending initial success with obs_increment7. In this case, begin 
! by doing the standard update for the mean for the observation. Then
! go back and evaluate the ratio of the prior difference between the
! ensemble mean and the obs and the expected value given the prior
! sample variance and the specified obs variance. If the ratio
! exceeds a threshold, then adjust the prior variance and obs variance
! by the same factor to reduce the ratio to the threshold (could also
! consider damping it toward the threshold by some factor to make
! the shocks less?). Then, recompute the updated variance using this
! and adjust the ensemble spread only. This means that a bad observation
! will increase the spread some but will not immediately pull the mean
! off toward the bad observation any further.

! ORIGINALLY identical to obs_increment7. Looking at changing inflate
! ratio for state versus obs.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8) :: error, diff_sd, ratio, inf_obs_var, inf_ens(ens_size)
real(r8) :: inf_obs_var_inv, inf_prior_cov, inf_prior_cov_inv
real(r8), parameter :: ratio_threshold = 2.0

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

! Compute the updated mean and covariance in the standard fashion
prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)
new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

! AT THIS POINT CAN EVALUATE INCONSISTENCY
! CHECK ALL THIS
error = prior_mean - obs
diff_sd = sqrt(obs_var + prior_cov)
ratio = abs(error / (ratio_threshold * diff_sd))

! Only modify if the ratio exceeds the threshold value
if(ratio < 1.0) ratio = 1.0

! Can now inflate by this ratio and then do adjustment
inf_obs_var = ratio**2 * obs_var
inf_obs_var_inv = 1.0 / inf_obs_var
! Form inflated ensemble
inf_ens = prior_mean + 1.2 * ratio * (ens - prior_mean)
inf_prior_cov = 1.2**2 * ratio**2 * prior_cov


inf_prior_cov_inv = 1.0 / inf_prior_cov
new_cov = 1.0 / (inf_prior_cov_inv + inf_obs_var_inv)

a = sqrt(new_cov * inf_prior_cov_inv)

obs_inc = a * (inf_ens - prior_mean) + new_mean - ens

end subroutine obs_increment13





subroutine obs_increment7(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment
! Third try messing around on 24 March, inflate obs_var and prior_var
! by same factor to make them consistent with difference between mean
! and obs and then proceed.

! In cases with no systematic error, this worked extremely well, slightly
! reducing error and vastly increasing spread. The test is exp_test2 was
! done with the threshold for correction beginning when the error is 
! more than 2 sd's of the expected.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8) :: error, diff_sd, ratio, inf_obs_var, inf_ens(ens_size)
real(r8) :: inf_obs_var_inv, inf_prior_cov, inf_prior_cov_inv

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)


! AT THIS POINT CAN EVALUATE INCONSISTENCY
error = prior_mean - obs
diff_sd = sqrt(obs_var + prior_cov)
ratio = abs(error / (2.0 * diff_sd))

! Only work with ratio's larger than 1?
if(ratio < 1.0) ratio = 1.0

! Can now inflate by this ratio and then do adjustment
inf_obs_var = ratio**2 * obs_var
inf_obs_var_inv = 1.0 / inf_obs_var
! Form inflated ensemble
inf_ens = prior_mean + ratio * (ens - prior_mean)
inf_prior_cov = ratio**2 * prior_cov

! CHECK ALL THIS
inf_prior_cov_inv = 1.0 / inf_prior_cov
new_cov = 1.0 / (inf_prior_cov_inv + inf_obs_var_inv)

new_mean = new_cov * (inf_prior_cov_inv * prior_mean + obs / inf_obs_var)

a = sqrt(new_cov * inf_prior_cov_inv)

obs_inc = a * (inf_ens - prior_mean) + new_mean - ens

end subroutine obs_increment7





subroutine obs_increment8(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! Uses a semi-formal model of systematic error. Here, assumes this
! is only in model. See notes from 1 April, 2003. Idea is that a
! distribution is hypothesized for an inflation factor, 1+alpha,
! and then the product of the probablility of a given value of
! alpha and the probability of the observation given this alpha
! is maximized by linear search.

! In cases with no systematic error, this worked extremely well, slightly
! reducing error and vastly increasing spread. The test is exp_test2 was
! done with the threshold for correction beginning when the error is 
! more than 2 sd's of the expected.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8) :: inf_ens(ens_size), alpha
! The tuning knob on the systematic error model, see 1 April, 01
real(r8), parameter :: sigma_alpha = 0.25

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

! Use maximum likelihood bias model to get inflation for model only
! Eventually want to split between model and obs
call bias_max_likelihood(prior_cov, sigma_alpha, abs(obs - prior_mean), &
   obs_var, alpha)

! Can now inflate by factor (1 + alpha)
! Form inflated ensemble
inf_ens = prior_mean + (1.0 + alpha) * (ens - prior_mean)
prior_cov = (1.0 + alpha)**2 * prior_cov


! CHECK ALL THIS

prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

a = sqrt(new_cov * prior_cov_inv)

obs_inc = a * (inf_ens - prior_mean) + new_mean - ens

end subroutine obs_increment8





subroutine obs_increment9(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! Uses a semi-formal model of systematic error. Here, assumes this
! is only in model. See notes from 1 April, 2003. Idea is that a
! distribution is hypothesized for an inflation factor, 1+alpha,
! and then the product of the probablility of a given value of
! alpha and the probability of the observation given this alpha
! is maximized by linear search.

! In cases with no systematic error, this worked extremely well, slightly
! reducing error and vastly increasing spread. The test is exp_test2 was
! done with the threshold for correction beginning when the error is 
! more than 2 sd's of the expected.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8) :: inf_ens(ens_size), alpha, inf_obs_var
! The tuning knob on the systematic error model, see 1 April, 01
real(r8), parameter :: sigma_alpha = 0.350

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

! Use maximum likelihood bias model to get inflation for model only
! Eventually want to split between model and obs
call bias_max_likelihood2(prior_cov, sigma_alpha, abs(obs - prior_mean), &
   obs_var, alpha)

! Can now inflate by factor (1 + alpha)
inf_obs_var = (1.0 + alpha)**2 * obs_var
obs_var_inv = 1.0 / inf_obs_var
! Form inflated ensemble
inf_ens = prior_mean + (1.0 + alpha) * (ens - prior_mean)
prior_cov = (1.0 + alpha)**2 * prior_cov


! CHECK ALL THIS

prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

new_mean = new_cov * (prior_cov_inv * prior_mean + obs / inf_obs_var)

a = sqrt(new_cov * prior_cov_inv)

obs_inc = a * (inf_ens - prior_mean) + new_mean - ens

end subroutine obs_increment9





subroutine obs_increment10(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! Uses a semi-formal model of systematic error. Here, assumes this
! is only in model. See notes from 1 April, 2003. Idea is that a
! distribution is hypothesized for an inflation factor, 1+alpha,
! and then the product of the probablility of a given value of
! alpha and the probability of the observation given this alpha
! is maximized by linear search.
! This version only inflates the variance. It adjusts the mean with
! the unmodified algorithm. Hope is to avoid being knocked around
! by sampling error in obs error while still dealing with persistent
! systematic error by increasing spread.

! In cases with no systematic error, this worked extremely well, slightly
! reducing error and vastly increasing spread. The test is exp_test2 was
! done with the threshold for correction beginning when the error is 
! more than 2 sd's of the expected.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8) :: inf_ens(ens_size), alpha, inf_obs_var, new_ens(ens_size)
real(r8) :: inf_obs_var_inv, inf_prior_cov, inf_prior_cov_inv, inf_new_cov
! The tuning knob on the systematic error model, see 1 April, 01
real(r8), parameter :: sigma_alpha = 0.50


integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)
prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)
new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

! Use maximum likelihood bias model to get inflation for model only
! Eventually want to split between model and obs
call bias_max_likelihood2(prior_cov, sigma_alpha, abs(obs - prior_mean), &
   obs_var, alpha)

!write(*, *) 'orig mean and cov ', prior_mean, prior_cov

! Can now inflate by factor (1 + alpha)
inf_obs_var = (1.0 + alpha)**2 * obs_var
inf_obs_var_inv = 1.0 / inf_obs_var
! Form inflated ensemble
inf_prior_cov = (1.0 + alpha)**2 * prior_cov

! CHECK ALL THIS
inf_prior_cov_inv = 1.0 / inf_prior_cov
inf_new_cov = 1.0 / (inf_prior_cov_inv + inf_obs_var_inv)

! Don't let obs increase uncertainty
if(inf_new_cov > 1.0 * prior_cov) inf_new_cov = 1.0 * prior_cov

a = sqrt(inf_new_cov * prior_cov_inv)
obs_inc = a * (ens - prior_mean) + new_mean - ens
write(*, *) 'sd factor a is ', a

end subroutine obs_increment10




subroutine bias_max_likelihood(prior_var, sigma_alpha, delta, obs_var, max_alpha)

implicit none

real(r8), intent(in) :: prior_var, sigma_alpha, delta, obs_var
real(r8), intent(out) :: max_alpha

real(r8) :: prior_sd, obs_sd
real(r8) :: alpha, sigma_t, prob_obs, prob_alpha, prob
real(r8) :: max_prob, max_prob_obs, max_prob_alpha
real(r8), parameter :: pi = 3.14159
integer :: i, max_i

write(*, *) '-----------------------------------------------'
write(*, *) sqrt(prior_var), delta, sqrt(obs_var)

! Zero out the max_prob initial
max_prob = 0.0

! Convert variance to SD
prior_sd = sqrt(prior_var)
obs_sd = sqrt(obs_var)

! Do an outward search on alpha for a maximum likelihood state
do i = 0, 100
   alpha = i * 0.02

! Compute the expected variance of the difference
   sigma_t = sqrt(((1.0 + alpha) * prior_sd)**2 + obs_sd**2)

! Probability that obs would be taken is
   prob_obs = (1.0 / (sigma_t * sqrt(2.0) * pi)) * exp((delta / sigma_t)**2 / (-2.0))
!  NOTE: MULTIPLYING BY 2.0 BECAUSE DISTRIBUTION IS ONE-SIDED; CHECK THIS
! TURNS OUT THAT THIS DOESN"T MATTER FOR MINIMIZING???
   prob_alpha = 2.0 * (1.0 / (sigma_alpha * sqrt(2.0) * pi)) * exp((alpha / sigma_alpha)**2 / (-2.0))
   prob = prob_obs * prob_alpha

!   write(*, 11) i, alpha, prob_obs, prob_alpha, prob
! 11 format(1x, i3, 1x, 4(e10.4, 1x))

! Keep the maximum likelihood state
   if(max_prob < prob) then
      max_prob = prob
      max_prob_obs = prob_obs
      max_prob_alpha = prob_alpha
      max_alpha = alpha
      max_i = i
   end if
end do

!write(*, *) 'max       i     alpha        prob_obs       prob_alpha      prob'
write(*, 21) max_alpha, max_prob_obs, max_prob_alpha, max_prob
21 format(1x, 4(e10.4, 1x))

end subroutine bias_max_likelihood





subroutine bias_max_likelihood2(prior_var, sigma_alpha, delta, obs_var, max_alpha)

implicit none

real(r8), intent(in) :: prior_var, sigma_alpha, delta, obs_var
real(r8), intent(out) :: max_alpha

real(r8) :: prior_sd, obs_sd
real(r8) :: alpha, sigma_t, prob_obs, prob_alpha, prob
real(r8) :: max_prob, max_prob_obs, max_prob_alpha
real(r8), parameter :: pi = 3.14159
integer :: i, max_i

write(*, *) '-----------------------------------------------'
write(*, *) sqrt(prior_var), delta, sqrt(obs_var)

! Zero out the max_prob initial
max_prob = 0.0

! Convert variance to SD
prior_sd = sqrt(prior_var)
obs_sd = sqrt(obs_var)

! Do an outward search on alpha for a maximum likelihood state
do i = 0, 100
   alpha = i * 0.02

! Compute the expected variance of the difference
   sigma_t = sqrt(((1.0 + alpha) * prior_sd)**2 + &
       ((1.0 + alpha) * obs_sd)**2)

! Probability that obs would be taken is
   prob_obs = (1.0 / (sigma_t * sqrt(2.0) * pi)) * exp((delta / sigma_t)**2 / (-2.0))
!  NOTE: MULTIPLYING BY 2.0 BECAUSE DISTRIBUTION IS ONE-SIDED; CHECK THIS
! TURNS OUT THAT THIS DOESN"T MATTER FOR MINIMIZING???
   prob_alpha = 2.0 * (1.0 / (sigma_alpha * sqrt(2.0) * pi)) * exp((alpha / sigma_alpha)**2 / (-2.0))
   prob = prob_obs * prob_alpha

!   write(*, 11) i, alpha, prob_obs, prob_alpha, prob
! 11 format(1x, i3, 1x, 4(e10.4, 1x))

! Keep the maximum likelihood state
   if(max_prob < prob) then
      max_prob = prob
      max_prob_obs = prob_obs
      max_prob_alpha = prob_alpha
      max_alpha = alpha
      max_i = i
   end if
end do

!write(*, *) 'max       i     alpha        prob_obs       prob_alpha      prob'
write(*, 21) max_alpha, max_prob_obs, max_prob_alpha, max_prob
21 format(1x, 4(e10.4, 1x))

end subroutine bias_max_likelihood2





subroutine obs_increment6(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment
! Second try on 24 March to deal with bias. this time don't reduce
! spread when there seems to be a large bias.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8) :: error, diff_sd, ratio, inv_ratio, delta_cov, rev_new_cov

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

! TEMPORARY LOOK AT INFLATING HERE; see notes from 12 Sept. 2001
!!!cov = cov * 1.01 OR prior_cov = prior_cov * 1.01

prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

! Look for signs of bias
error = prior_mean - obs
diff_sd = sqrt(obs_var + prior_cov)
ratio = abs(error / diff_sd)
if(ratio > 1.0) then
   inv_ratio = 1.0 / ratio
   delta_cov = new_cov - prior_cov
   delta_cov = inv_ratio * delta_cov
   rev_new_cov = prior_cov + delta_cov
else
   rev_new_cov = new_cov
endif


!write(*, *) 'error, diff_sd ', error, diff_sd
!write(*, *) 'inv_ratio ', inv_ratio
!write(*, *) 'new, rev_new_cov ', new_cov, rev_new_cov

 



new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

a = sqrt(rev_new_cov * prior_cov_inv)

obs_inc = a * (ens - prior_mean) + new_mean - ens

end subroutine obs_increment6





subroutine obs_increment5(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc)
!
! EAKF version of obs increment
! Test version from 24 March, 2003 looking at local explicit handling
! of prior bias.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8):: prior_cov, sx, s_x2
real(r8) :: error, diff_sd, bias

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior covariance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

! TEMPORARY LOOK AT INFLATING HERE; see notes from 12 Sept. 2001
!!!cov = cov * 1.01 OR prior_cov = prior_cov * 1.01

prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

a = sqrt(new_cov * prior_cov_inv)

obs_inc = a * (ens - prior_mean) + new_mean - ens


! Look for bias based inconsistency in the original distribution and
! do explicit correction
error = prior_mean - obs
diff_sd = sqrt(obs_var + prior_cov)
if(abs(error) > 3.0 * diff_sd) then
   bias = (abs(error) - 3.0 * diff_sd) * (error / abs(error))
else
   bias = 0.0
endif

obs_inc = obs_inc - bias

!write(*, *) 'prior_mean, obs, error', prior_mean, obs, error
!write(*, *) 'obs_sd, prior_sd, diff_sd ', sqrt(obs_var), sqrt(prior_cov), diff_sd
!write(*, *) 'computed bias ', bias

end subroutine obs_increment5



subroutine obs_inc_index(ens, ens_size, obs, obs_var, obs_inc, index)
!========================================================================
! subroutine obs_inc_index(ens, ens_size, obs, obs_var, obs_inc, index)
!

! EAKF version of obs increment with closest obs neighbor indexing
! See notes from 14 November, 2001

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
integer, intent(out) :: index(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8) :: prior_cov, sx, s_x2
real(r8) :: new_obs(ens_size), min_dist

integer :: i, j

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute the regression directly
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

new_mean = new_cov * (prior_cov_inv * prior_mean + obs / obs_var)

a = sqrt(new_cov * prior_cov_inv)

new_obs = a * (ens - prior_mean) + new_mean
! For each update ensemble member, find closest original member and an increment
do i = 1, ens_size
   min_dist = dabs(new_obs(i) - ens(1))
   index(i) = 1
   do j = 2, ens_size
      if(dabs(new_obs(i) - ens(j)) < min_dist) then
         min_dist = dabs(new_obs(i) - ens(j))
         index(i) = j
      end if
   end do

! Need signed increment
   obs_inc(i) = new_obs(i) - ens(index(i))
!   write(*, *) i, index(i), real(obs_inc(i)), real(new_obs(i) - ens(i))
end do


end subroutine obs_inc_index



subroutine obs_increment4(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment4(ens, ens_size, obs, obs_var, obs_inc)
!

! ENKF version of obs increment

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean(ens_size)
real(r8) :: sx, s_x2, prior_cov
real(r8) :: temp_mean, temp_obs(ens_size)
integer :: ens_index(ens_size), new_index(ens_size)

integer :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior mean and covariance
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

! If this is first time through, need to initialize the random sequence
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq)
   first_inc_ran_call = .false.
endif

! Generate perturbed obs
do i = 1, ens_size
    temp_obs(i) = random_gaussian(inc_ran_seq, obs, sqrt(obs_var))
end do

! Move this so that it has original obs mean
temp_mean = sum(temp_obs) / ens_size
temp_obs(:) = temp_obs(:) - temp_mean + obs

! Loop through pairs of priors and obs and compute new mean
do i = 1, ens_size
   new_mean(i) = new_cov * (prior_cov_inv * ens(i) + temp_obs(i) / obs_var)
   obs_inc(i) = new_mean(i) - ens(i)
end do

! Try sorting to make increments as small as possible
!call index_sort(ens, ens_index, ens_size)
!call index_sort(new_mean, new_index, ens_size)
!do i = 1, ens_size
!   obs_inc(ens_index(i)) = new_mean(new_index(i)) - ens(ens_index(i))
!end do

end subroutine obs_increment4



subroutine obs_inc4_index(ens, ens_size, obs, obs_var, obs_inc, index)
!========================================================================
! subroutine obs_inc4_index(ens, ens_size, obs, obs_var, obs_inc, index)

! ENKF version of obs increment
! This is the version that uses local linear updates by finding nearest
! neighbor and increment from that neighbor.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
integer, intent(out) :: index(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean(ens_size)
real(r8) :: sx, s_x2, prior_cov
real(r8) :: temp_mean, temp_obs(ens_size), min_dist
integer :: ens_index(ens_size), new_index(ens_size)

integer :: i, j

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Test of computing directly
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

prior_cov_inv = 1.0 / prior_cov
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

! If this is first time through, need to initialize the random sequence
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq)
   first_inc_ran_call = .false.
endif

! Generate perturbed obs
do i = 1, ens_size
    temp_obs(i) = random_gaussian(inc_ran_seq, obs, sqrt(obs_var))
end do
! Move this so that it has original obs mean
temp_mean = sum(temp_obs) / ens_size
temp_obs(:) = temp_obs(:) - temp_mean + obs

! Loop through pairs of priors and obs and compute new mean
do i = 1, ens_size
   new_mean(i) = new_cov * (prior_cov_inv * ens(i) + temp_obs(i) / obs_var)
end do

! For each update ensemble member, find closest original member and an increment
do i = 1, ens_size
   min_dist = dabs(new_mean(i) - ens(1))
   index(i) = 1
   do j = 2, ens_size
      if(dabs(new_mean(i) - ens(j)) < min_dist) then
         min_dist = dabs(new_mean(i) - ens(j))
         index(i) = j
      end if
   end do

! Need signed increment
   obs_inc(i) = new_mean(i) - ens(index(i))
!   write(*, *) i, index(i), real(obs_inc(i)), real(new_mean(i) - ens(i))
end do

end subroutine obs_inc4_index



subroutine obs_increment3(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment3(ens, ens_size, obs, obs_var, obs_inc)
!

! Kernel version of obs increment

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, prior_cov
real(r8) :: sx, s_x2
real(r8) :: weight(ens_size), new_mean(ens_size)
real(r8) :: cum_weight, total_weight, cum_frac(ens_size)
real(r8) :: unif, norm, new_member(ens_size)

integer :: i, j, kernel, ens_index(ens_size), new_index(ens_size)

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Compute prior mean and covariance
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)


! For kernels, scale the prior covariance
prior_cov = prior_cov / 10.0

prior_cov_inv = 1.0 / prior_cov

! Compute new covariance once for these kernels
new_cov = 1.0 / (prior_cov_inv + obs_var_inv)

! New mean is computed ens_size times as is weight
do i = 1, ens_size
   new_mean(i) = new_cov*(prior_cov_inv * ens(i) + obs / obs_var)
   weight(i) =  2.71828 ** (-0.5 * (ens(i)**2 * prior_cov_inv + &
      obs**2 * obs_var_inv - new_mean(i)**2 / new_cov))
end do

! Compute total weight
total_weight = sum(weight)
cum_weight = 0.0
do i = 1, ens_size
   cum_weight = cum_weight + weight(i)
   cum_frac(i) = cum_weight / total_weight
end do

! If this is first time through, need to initialize the random sequence
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq)
   first_inc_ran_call = .false.
endif

! Generate a uniform random number and a Gaussian for each new member
do i = 1, ens_size
   unif = random_uniform(inc_ran_seq)
! Figure out which kernel it's in
   do j = 1, ens_size
      if(unif < cum_frac(j)) then
         kernel = j
         goto 10
      end if
   end do
10 continue

! Next calculate a unit normal in this kernel
   norm = random_gaussian(inc_ran_seq, dble(0.0), sqrt(new_cov))
! Now generate the new ensemble member
   new_member(i) = new_mean(kernel) + norm
end do

! Try sorting to make increments as small as possible
call index_sort(ens, ens_index, ens_size)
call index_sort(new_member, new_index, ens_size)

do i = 1, ens_size
   obs_inc(ens_index(i)) = new_member(new_index(i)) - ens(ens_index(i))
end do

end subroutine obs_increment3



subroutine obs_increment2(ens, ens_size, obs, obs_var_in, obs_inc)
!========================================================================
! subroutine obs_increment2(ens, ens_size, obs, obs_var_in, obs_inc)
!

! This is a research version of obs increment using generalized 
! distributions. Should not be used in present from.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var_in
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: ens2(1, ens_size), a, obs_var_inv, cov(1, 1)
real(r8) :: mean(1), prior_mean, prior_cov_inv, new_cov, new_mean
real(r8) :: prior_cov, sx, s_x2

integer, parameter :: num = 10000
integer :: i, indx, cur_index, j , ind_sort(ens_size)
real(r8) :: tot_dense, lctn(ens_size), target
real(r8) :: obs_dense(num), state_dense(num), dense(num)
real(r8) :: cum_dense(0:num) 
real(r8) :: kernel_var, mxmum, mnmum, rnge, x(num)
real(r8) :: min_ens, min_obs, max_ens, max_obs


real(r8) :: obs_var


!!! TEST, GET RID OF THIS
obs_var = 1e10


! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0 / obs_var

! Need to compute the prior covariance and mean; copy to 2d for now but
! get rid of this step at some point; just fits interface
!ens2(1, :) = ens
!call sample_cov(ens2, 1, ens_size, cov, xmean = mean)
!prior_mean = mean(1)
!prior_cov = cov(1, 1)


! Test of computing directly
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_cov = (s_x2 - sx**2 / ens_size) / (ens_size - 1)


prior_cov_inv = 1.0 / prior_cov

! Need to fix a better kernel width later!!!!!
kernel_var = prior_cov / 5.0

! WARNING: NEED TO SORT AND GO

! Dumbest possible way, just partition and search
! For now, do 1 sd past max and min
mxmum = maxval(ens) + sqrt(prior_cov)
mnmum = minval(ens) - sqrt(prior_cov)

!!!max_ens = maxval(ens) + sqrt(prior_cov)
!!!max_obs = obs + 2.0 * sqrt(obs_var)
!!!if(max_ens > max_obs) then
!!!   mxmum = max_ens
!!!else
!!!   mxmum = max_obs
!!!endif

!!!min_ens = minval(ens) - sqrt(prior_cov)
!!!min_obs = obs - 2.0 * sqrt(obs_var)
!!!if(min_ens < min_obs) then
!!!   mnmum = min_ens
!!!else
!!!   mnmum = min_obs
!!!endif

rnge = mxmum - mnmum

!write(*, *) 'min, max, range ', real(mnmum), real(mxmum), real(rnge)
!write(*, *) 'obs . var ', obs, obs_var

cum_dense(0) = 0

do i = 1, num
   x(i) = mnmum + (i - 1.0) * rnge / (num - 1.0)
   obs_dense(i) = exp((x(i) - obs)**2 / (-2. * obs_var))
   state_dense(i) = 0.0
   do j = 1, ens_size
      state_dense(i) = state_dense(i) + exp((x(i) - ens(j))**2 / (-2. * kernel_var))
   end do
   dense(i) = obs_dense(i) * state_dense(i)
   cum_dense(i) = cum_dense(i - 1) + dense(i)
!   write(*, 11) i, real(x(i)), real(obs_dense(i)), real(state_dense(i)), real(dense(i)), cum_dense(i)
11 format(1x, i3, 1x, 5(e10.4, 1x))
end do
tot_dense = cum_dense(num)

! Figure out where the new ensemble members should be located
! Problem, outside of boxes? Error for now
target = (tot_dense / (ens_size + 1))
if(target < cum_dense(1)) then
   write(*, *) 'ERROR:first target outside of box in increment_obs'
end if
target = ens_size * (tot_dense / (ens_size + 1))
if(target > cum_dense(num)) then
   write(*, *) 'ERROR: last target outside of box in increment_obs'
endif

indx = 1
do i = 2, num
! At what cumulative density should next point go
   cur_index = indx
   do j = cur_index, ens_size
      target = indx * (tot_dense / (ens_size + 1))
! Search to see if some should go between next indices
      if(target > cum_dense(i)) goto 10
! Otherwise, in this bin, linearly interpolate
      lctn(indx) = x(i - 1) + (target - cum_dense(i - 1)) / &
         (cum_dense(i) - cum_dense(i - 1)) * (x(i) - x(i - 1))
!      write(*, *) 'lctn ', indx, target, lctn(indx)
! Increment index for position being sought, if all found exit outer loop
      indx = indx + 1
      if(indx > ens_size) goto 20
   end do
10 continue
end do

20 continue

! Sort the ensemble
call index_sort(ens, ind_sort, ens_size)
do i = 1, ens_size
!   write(*, *) i, ens(ind_sort(i)), lctn(i)
end do

do i = 1, ens_size
   obs_inc(ind_sort(i)) = lctn(i) - ens(ind_sort(i))
end do

do i = 1, ens_size
!   write(*, *) i, ens(i), obs_inc(i)
end do

if(1 == 1) stop

end subroutine obs_increment2



subroutine update_from_obs_inc(obs, obs_inc, state, ens_size, &
               a, state_inc, cov_factor)
!========================================================================
! subroutine update_from_obs_inc(obs, obs_inc, state, ens_size, &
!                state_inc, cov_factor)
!

! Does linear regression of a state variable onto an observation and
! computes state variable increments from observation increments

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: obs(ens_size), obs_inc(ens_size)
real(r8), intent(in) :: state(ens_size), cov_factor
real(r8), intent(inout) :: a
real(r8), intent(out) :: state_inc(ens_size)

real(r8) :: sum_x, t(ens_size), sum_t2, sum_ty, reg_coef
real(r8) :: mean_inc, inf_state(ens_size)

real(r8) :: sum_y, sum_y2, sum_x2, correl
real(r8) :: factor, state_var_norm


! For efficiency, just compute regression coefficient here
sum_x  = sum(obs)
t      = obs - sum_x/ens_size
sum_t2 = sum(t * t)
sum_ty = sum(t * state)

if (sum_t2 /= 0.0_r8) then
   reg_coef = sum_ty/sum_t2
else
   reg_coef = 0.0_r8
endif

! Compute the sample correlation
sum_y = sum(state)
sum_y2 = sum(state**2)
state_var_norm = sum_y2 - sum_y**2 / ens_size
correl = reg_coef * sqrt(sum_t2 / state_var_norm)

! The original form follows
state_inc = cov_factor * reg_coef * obs_inc

if(1 == 1) return

! What was the 'a'; reduction in standard deviation for the ob
! Remove mean obs_inc first?
! WARNING: WOULD BE MORE ROBUST TO PASS A IN
!a = 1.0 - abs((obs_inc(1) - sum(obs_inc) / ens_size) / (obs(1) - sum_x / ens_size))

!write(*, *) 'a in update_from is ', a

! Add in a slope factor for continuity
! Following flat triad works fairly well for base L96 cases
! FOLLOWING VALUES WORK FOR 20 and 40 ENSEMBLE MEMBERS BUT NOT FOR 10
if(abs(correl) < 0.3) then
   factor = 1.0 / (1.0 - 1.00 * (1.0 - a) / 30.0) 
else if(abs(correl) < 0.6) then
   factor = 1.0 / (1.0 - 1.00 * (1.0 - a) / 35.0) 
else
   factor = 1.0
endif

!FOLLOWING ONE WORKS PRETTY WELL FOR 10 MEMBER ENSEMBLES IN L96!
!if(abs(correl) < 0.40) then
!   factor = 1.0 / (1.0 - 1.00 * (1.0 - a) / 9.0) 
!else if(abs(correl) < 0.7) then
!   factor = 1.0 / (1.0 - 1.00 * (1.0 - a) / 35.0) 
!else
!   factor = 1.0 
!endif

if(abs(correl) < 1.00) then
   inf_state = factor * (state - sum_y / ens_size) + (sum_y / ens_size)
   state_inc = state_inc + inf_state - state
endif


end subroutine update_from_obs_inc



subroutine robust_update_from_obs_inc(obs, obs_inc, state, ens_size, &
                                      state_inc, cov_factor, index_in)
!========================================================================
! subroutine robust_update_from_obs_inc(obs, obs_inc, state, ens_size, &
!                                       state_inc, cov_factor, index_in)
!

! Does a robust linear regression of a state variable onto an observation
! variable. Not rigorously tested at this point.

implicit none

integer, intent(in) :: ens_size
integer, intent(in), optional :: index_in(ens_size)
real(r8), intent(in) :: obs(ens_size), obs_inc(ens_size)
real(r8), intent(in) :: state(ens_size), cov_factor
real(r8), intent(out) :: state_inc(ens_size)

real(r8) :: ens(2, ens_size), cov(2, 2), reg
real(r8) :: y_lo, y_hi, x_lo, x_hi
integer :: index_lo(ens_size), i
integer :: lo_start, lo_end, hi_start, hi_end, lo_size, hi_size

! Need to sort the observations to be able to find local structure
! For efficiency, would like to sort this only once in main program
if(present(index_in)) then
   index_lo = index_in
else
   call index_sort(obs, index_lo, ens_size)
endif

! Create combined matrix for covariance
do i = 1, ens_size
   ens(1, i) = state(index_lo(i))
   ens(2, i) =   obs(index_lo(i))
end do

! Get mean values of lower and upper half of data sorted by index
!lo_start = 1
lo_start = ens_size / 8
lo_end = ens_size / 2
!lo_end = 3 * ens_size / 8
lo_size = lo_end - lo_start + 1
hi_start = ens_size / 2 + 1
!hi_start = 5 * ens_size / 8
!hi_end = ens_size
hi_end = 7 * ens_size / 8
hi_size = hi_end - hi_start + 1

! Method 1: use means of data halves sorted by obs
!y_lo = sum(ens(2, lo_start:lo_end)) / lo_size
!y_hi = sum(ens(2, hi_start:hi_end)) / hi_size
x_lo = sum(ens(1, lo_start:lo_end)) / lo_size
x_hi = sum(ens(1, hi_start:hi_end)) / hi_size



! Method 2: median values of data halves???;  unstable when used for x
y_lo = ens(2, ens_size / 4)
y_hi = ens(2, 3 * ens_size / 4)
!x_lo = ens(1, ens_size / 4)
!x_hi = ens(1, 3 * ens_size / 4)

reg = (x_hi - x_lo) / (y_hi - y_lo)

state_inc = cov_factor * reg * obs_inc

end subroutine robust_update_from_obs_inc



subroutine local_update_from_obs_inc(obs, obs_inc, state, ens_size, &
                     state_inc, cov_factor, half_num_neighbors, index_in)
!========================================================================
! subroutine local_update_from_obs_inc(obs, obs_inc, state, ens_size, &
!                      state_inc, cov_factor, half_num_neighbors, index_in)
!

! First stab at doing local linear regression on a set of num_neighbors
! points around a particular obs.

implicit none

integer, intent(in) :: ens_size, half_num_neighbors
integer, intent(in), optional :: index_in(ens_size)
real(r8), intent(in) :: obs(ens_size), obs_inc(ens_size)
real(r8), intent(in) :: state(ens_size), cov_factor
real(r8), intent(out) :: state_inc(ens_size)

real(r8) :: ens(2, ens_size), cov(2, 2), y2(ens_size), xy(ens_size)
real(r8) :: sx, sy, s_y2, sxy, reg
integer :: index_lo(ens_size), i, lower, upper, num_neighbors

! Make sure num_neighbors isn't too big to make sense
num_neighbors = 2 * half_num_neighbors
if(num_neighbors > ens_size - 1) then
   write(*, *) 'num_neighbors too big in local_update_from_obs_inc'
   stop
endif

! Need to sort the observations to be able to find local structure
! For efficiency, would like to sort this only once in main program
if(present(index_in)) then
   index_lo = index_in
else
   call index_sort(obs, index_lo, ens_size)
endif


! Temporary look at sorting on data, not obs
! Seems to work much better for original L96 cases
call index_sort(state, index_lo, ens_size)




! Load up the ensemble array to have state and obs sorted by obs
do i = 1, ens_size
   ens(1, i) = state(index_lo(i))
   ens(2, i) =   obs(index_lo(i))
!   write(*, *) 'state obs ', i, real(ens(1, i)), real(ens(2, i))
end do

! Do initial preparation for repeated regressions
xy = ens(1, :) * ens(2, :)
y2 = ens(2, :) ** 2

! Loop through to use nearest num_neighbors points for regression
do i = 1, ens_size
   if(i == 1 .or. (i > half_num_neighbors .and. &
      i <= ens_size - half_num_neighbors)) then
      lower = i - half_num_neighbors
      if(lower < 1) lower = 1
      if(lower > ens_size - num_neighbors) lower = ens_size - num_neighbors
      upper = lower + num_neighbors
!      write(*, *) 'upper lower ', i,lower, upper

      sx = sum(ens(1, lower:upper))
      sy = sum(ens(2, lower:upper))
      s_y2 = sum(y2(lower:upper))
      sxy = sum(xy(lower:upper))
      reg = (sxy - sx * sy / (num_neighbors + 1)) / &
             (s_y2 - sy**2 / (num_neighbors + 1))
!      write(*, *) 'reg ', i, real(reg)
   endif
   state_inc(index_lo(i)) = (cov_factor * reg) * obs_inc(index_lo(i))

end do

end subroutine local_update_from_obs_inc


!========================================================================

subroutine linear_obs_increment(ens, ens_size, obs, obs_var, var_inc, mean_inc, sd_ratio)

! EAKF version of obs increment, obs_update_inflate inflates the variance after
! computation of the mean. (See notes from 17-19 Sept. 2002)
! Uses additional linear corrections from notes in early Dec. 2002 plus modified
! mean versus variance delta's from that analysis.

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: var_inc(ens_size), mean_inc, sd_ratio

real(r8) :: a, prior_mean, prior_var, new_mean, new_var, sx, s_x2

! Compute prior variance and mean from sample
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / ens_size
prior_var = (s_x2 - sx**2 / ens_size) / (ens_size - 1)

!!!write(*, *) 'in obs_inc prior, obs_var ', real(prior_var), real(obs_var)

! Compute updated variance and mean
new_var = 1.0 / (1.0 / prior_var + 1.0 / obs_var)
new_mean = new_var * (prior_mean / prior_var + obs / obs_var)

a = sqrt(new_var / prior_var)

!!!write(*, *) 'a in obs_increment is ', a

! Compute the increments for each ensemble member
var_inc = a * (ens - prior_mean) + prior_mean - ens
mean_inc = new_mean - prior_mean

! Return a as sd_ratio
sd_ratio = a

end subroutine linear_obs_increment

!========================================================================

subroutine linear_update_from_obs_inc(obs, var_inc, mean_inc, state, &
   ens_size, state_inc, cov_factor, sd_ratio)

! Does linear regression of a state variable onto an observation and
! computes state variable increments from observation increments

implicit none

integer, intent(in) :: ens_size
real(r8), intent(in) :: obs(ens_size), var_inc(ens_size), mean_inc
real(r8), intent(in) :: state(ens_size), cov_factor, sd_ratio
real(r8), intent(out) :: state_inc(ens_size)

real(r8) :: sum_x, sum_y, sum_xy, sum_x2, reg_coef
real(r8) :: sum_y2, correl, linear_factor

real(r8) :: new_reg, ybar, yvar, xbar, xyvar



! For efficiency, just compute regression coefficient here
sum_x = sum(obs)
sum_y = sum(state)
sum_xy = sum(obs * state)
sum_x2 = sum(obs * obs)


!------------------------------------------------
! Temporary computation of correlation
!sum_y2 = sum(state * state)

!correl = (ens_size * sum_xy - sum_x * sum_y) / &
!   sqrt((ens_size * sum_x2 - sum_x**2) * (ens_size * sum_y2 - sum_y**2))
!write(*, *) 'correl/factor in update_from is ', real(correl)
!----------------------------------------------------

reg_coef = (ens_size * sum_xy - sum_x * sum_y) / (ens_size * sum_x2 - sum_x**2)

!!!write(*, *) 'update_from_obs_inc reg_coef is ', reg_coef

! Use the full ratio expression from 3 Dec. '02 Notes
if(cov_factor > 0.99999) then
   linear_factor = cov_factor
else
   linear_factor = (sqrt(cov_factor * sd_ratio**2 - cov_factor + 1.0) - 1.0) / &
      (sd_ratio - 1.0)
endif

!!!write(*, *) 'sd ratio ', sd_ratio
!!!write(*, *) 'cov_factor, linear_factor ', cov_factor, linear_factor
!!!write(*, *) 'ratio ', cov_factor / linear_factor

! Following line moves mean more (may not be right)
!!!state_inc = linear_factor * reg_coef * var_inc + cov_factor * reg_coef * mean_inc

! Test of moving mean the same as the rest of things
! Appears to be correct solution
!!!state_inc = linear_factor * reg_coef * var_inc + linear_factor * reg_coef * mean_inc

! TEST OF TAKING SQUARE OF REGRESSION COEFFICIENT ???
! WORKS GREAT IN INITIAL TESTS!!!
!!!reg_coef = reg_coef * sqrt(abs(reg_coef))
! Try multiplying by correl to reduce impact
!reg_coef = reg_coef * abs(correl)
!reg_coef = reg_coef * 0.01
state_inc = linear_factor * reg_coef * var_inc + linear_factor * reg_coef * mean_inc



end subroutine linear_update_from_obs_inc



!-----------------------------------------------------------------------

subroutine look_for_bias(ens, n, obs, obs_var, var_ratio)

implicit none

integer, intent(in) :: n
real(r8), intent(in) :: ens(n), obs, obs_var
real(r8), intent(out) :: var_ratio

real(r8) :: sx, s_x2, prior_mean, prior_var, sq_err, tot_var

! Compute variance of the ensemble prior for this obs
sx = sum(ens)
s_x2 = sum(ens * ens)
prior_mean = sx / n
prior_var = (s_x2 - sx**2 / n) / (n - 1)

! Variance of difference between obs and mean should be sum of variances
!!!sq_err = (obs - prior_mean)**2
sq_err = sum((obs - ens)**2) / n
tot_var = obs_var + prior_var
var_ratio = sq_err / tot_var

end subroutine look_for_bias

!========================================================================



!========================================================================
! end module assim_tools_mod
!========================================================================

end module assim_tools_mod

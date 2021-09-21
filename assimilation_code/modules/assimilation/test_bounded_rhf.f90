program test_bounded_rhf

! Test specific cases for obs_increment_bounded_norm_rhf
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities
use types_mod,         only : r8
use assim_tools_mod,   only : obs_increment_bounded_norm_rhf, &
                              obs_increment_rank_histogram
use random_seq_mod,    only : random_seq_type, init_random_seq, random_gaussian

implicit none

real(r8), parameter :: truth = 0.99_r8
integer,  parameter :: ens_size = 80
integer,  parameter :: num_steps = 10000
integer,  parameter :: bounds_case = 4         ! Case 4 is doubly bounded
integer,  parameter :: init_conditions_case = 2
real(r8), parameter :: bound(2) = (/0.0_r8, 1.0_r8/)

real(r8)              :: ens(ens_size), obs_inc(ens_size) 
real(r8)              :: obs, obs_var, obs_std, prior_mean, prior_var
logical               :: is_bounded(4, 2)
type(random_seq_type) :: r_seq
integer               :: i, j

! Initialize for error handler
call initialize_mpi_utilities('test_bounded_rhf')

! Initialize a repeating random sequence, seed fixed at 1
call init_random_seq(r_seq, 1)

! The 4 possible bounding cases: none, left, right, both, selected by bounds_case
is_bounded(1:4, 1) = (/.false., .true., .false., .true./)
is_bounded(1:4, 2) = (/.false., .false., .true., .true./)

! The three possible initial conditions cases, selected by init_conditions_case
! 1 -> initial ensemble is uniform with no members at bounds
! 2 -> initial ensemble is uniform with a member at each of the bounds
! 3 -> initial ensemble is two delta functions, half of members are lower bound, half are upper

! Set the observation error variance
if(truth > 0.1_r8) then 
   obs_std = 0.15 * truth
else
   ! NOTE: Go back to 0.15 for small truth to compare 0 and 1 truth cases
   obs_std = 0.15;
endif
obs_var = obs_std**2

! Initial ensemble is uninformative over [0 1]
do i = 1, ens_size
   if(init_conditions_case == 1) then
      ! Uniform ensemble that avoids the bounds
      ens(i) = (i * 1.0_r8) / (ens_size + 1.0_r8)
   elseif(init_conditions_case == 2) then
      ! Uniform ensemble that includes the bounds
      ens(i) = (i - 1.0_r8) / (ens_size - 1.0_r8)
   elseif(init_conditions_case == 3) then
      ! Gnarly prior that has either 0 or 1
      ens = 0.0_r8;
      ens(ens_size / 2 + 1:ens_size) = 1.0_r8
   else
   ! Die an ignominius death
      write(*, *) 'init_conditions_case must be 1, 2, or 3'
      stop
   endif
end do

! Loop through a sequence of observations
do i = 1, num_steps

   ! Generate truncated normal observation
   obs = 99.0_r8
   do while(obs >= 1.0_r8 .or. obs <= 0.0_r8)
      obs = random_gaussian(r_seq, truth, sqrt(obs_var))
   end do

   ! Compute prior variance and mean from sample
   prior_mean = sum(ens) / ens_size
   prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)

   ! Naive output of ensemble mean and variance
   write(41, *) prior_mean, prior_var

   ! Output of entire ensemble
   do j = 1, ens_size
      write(42, *) j, ens(j)
   end do

   ! Use the bounded rhf
   call obs_increment_bounded_norm_rhf(ens, ens_size, prior_var, obs, obs_var, &
      obs_inc, is_bounded(bounds_case, :), bound)

   ! Baseline check for original rhf
   !call obs_increment_rank_histogram(ens, ens_size, prior_var, obs, obs_var, &
      !obs_inc)

   ! Update the ensemble
   ens = ens + obs_inc

end do

! Metadata for entire ensemble
write(42, *) ens_size, truth

call finalize_mpi_utilities()

end program test_bounded_rhf

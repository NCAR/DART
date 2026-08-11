! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! Test the array level routines in perturb_mod.
!
! which take the copies array and the ensemble distribution as data (not ens_handle),
! so a serial program can stand in for any number of tasks.  No MPI required.

program perturb_test

use types_mod,     only : r8, i8
use utilities_mod, only : initialize_utilities, finalize_utilities
use perturb_mod,   only : perturb_uniform, perturb_bitwise

use test,          only : ok, plan, isabs, pass

implicit none

! Deliberately more copies than ens_size.  filter carries extras (mean, sd,
! inflation) beyond the members, so a routine that assumed ens_size was
! size(copies,1) would trample them.
integer,  parameter :: ENS_SIZE   = 20
integer,  parameter :: NUM_COPIES = ENS_SIZE + 3
integer,  parameter :: NUM_VARS   = 40
real(r8), parameter :: AMPLITUDE  = 0.1_r8
real(r8), parameter :: TOL        = 1.0e-14_r8

! output_flag = .false. keeps DART's start/finish banners off stdout.  
! to stop 'ill formed plan' when using tapview
call initialize_utilities('perturb_test', output_flag=.false.)

call plan(14)

call test_uniform()
call test_bitwise_task_invariance()
call test_no_vars_on_this_task()

call finalize_utilities()

contains

!------------------------------------------------------------------
!> A single instance copied to every member,
!> need to not have the sd the same everywhere, so that perturb_uniform can be checked.

subroutine fill_single_instance(copies)

real(r8), intent(out) :: copies(:,:)

integer :: i

do i = 1, size(copies, 2)
   copies(:, i) = 3.0_r8 * i - 17.0_r8
enddo

end subroutine fill_single_instance

!------------------------------------------------------------------
!> Sample standard deviation over the ensemble members at one element.

function member_sd(copies, ivar)

real(r8), intent(in) :: copies(:,:)
integer,  intent(in) :: ivar
real(r8)             :: member_sd

real(r8) :: mean

mean      = sum(copies(1:ENS_SIZE, ivar)) / real(ENS_SIZE, r8)
member_sd = sqrt(sum((copies(1:ENS_SIZE, ivar) - mean)**2) / real(ENS_SIZE - 1, r8))

end function member_sd

!------------------------------------------------------------------
!> The point of the uniform perturbation: the spread is the same at every
!> state element, so the covariance of an observation with any state
!> variable is identical and localization is the only thing left varying.
!>
!> Also checks that ens_size is honoured as an argument rather than taken
!> from size(copies,1), which would overwrite filter's extra copies.

subroutine test_uniform()

real(r8) :: copies(NUM_COPIES, NUM_VARS)
real(r8) :: sds(NUM_VARS)
real(r8) :: expected_sd, dev1, worst_dev
integer  :: i, j

call fill_single_instance(copies)
call perturb_uniform(copies, ENS_SIZE, AMPLITUDE)

! spread is amplitude*sqrt(N(N+1)/12), and the same everywhere
expected_sd = AMPLITUDE * sqrt(real(ENS_SIZE, r8) * real(ENS_SIZE + 1, r8) / 12.0_r8)
do i = 1, NUM_VARS
   sds(i) = member_sd(copies, i)
enddo
call isabs(minval(sds), expected_sd, TOL, 'uniform: smallest sd is amplitude*sqrt(N(N+1)/12)')
call isabs(maxval(sds), expected_sd, TOL, 'uniform: largest sd is amplitude*sqrt(N(N+1)/12)')

! the member offsets must not depend on the state index
worst_dev = 0.0_r8
do j = 1, ENS_SIZE
   dev1 = copies(j, 1) - copies(1, 1)
   do i = 1, NUM_VARS
      worst_dev = max(worst_dev, abs((copies(j, i) - copies(1, i)) - dev1))
   enddo
enddo
call isabs(worst_dev, 0.0_r8, TOL, 'uniform: member offsets do not depend on the state index')

! member 1 is the instance that was read in, it must not move
call ok(all([(abs(copies(1, i) - (3.0_r8 * i - 17.0_r8)) < TOL, i = 1, NUM_VARS)]), &
        'uniform: member 1 is unchanged')

! the copies above ens_size belong to filter, not to us
call ok(all([((abs(copies(j, i) - (3.0_r8 * i - 17.0_r8)) < TOL, &
               j = ENS_SIZE + 1, NUM_COPIES), i = 1, NUM_VARS)]), &
        'uniform: copies beyond ens_size are untouched')

end subroutine test_uniform

!------------------------------------------------------------------
!> perturb_bitwise must give the same answer however the state is split
!> across tasks.  Run it single task (whole state), then split the state the
!> way ensemble_manager does (round robin) for several task counts and
!> check every owned element against the one task answer.

subroutine test_bitwise_task_invariance()

real(r8)    :: initial(NUM_COPIES, NUM_VARS)
real(r8)    :: reference(NUM_COPIES, NUM_VARS)
real(r8)    :: mine(NUM_COPIES, NUM_VARS)
integer(i8) :: all_vars(NUM_VARS), my_vars(NUM_VARS)
integer     :: ntasks, itask, i, j, my_num_vars
logical     :: matches
character(len=64) :: label

! perturb_bitwise overwrites every member including member 1, so keep the
! state as it was read in to seed each pretend task from
call fill_single_instance(initial)

! one task owns the whole state - this is the reference answer
do i = 1, NUM_VARS
   all_vars(i) = int(i, i8)
enddo
reference = initial
!                   (copies,    ens_size, my_vars,  num_vars,          amplitude)
call perturb_bitwise(reference, ENS_SIZE, all_vars, int(NUM_VARS, i8), AMPLITUDE)

! a routine that quietly did nothing would pass every comparison below
call ok(abs(reference(2,1) - reference(1,1)) > TOL, &
        'bitwise: the reference run actually perturbed the state')

do ntasks = 2, 7

   matches = .true.

   do itask = 0, ntasks - 1

      ! round robin following what ensemble_manager does
      my_num_vars = 0
      do i = itask + 1, NUM_VARS, ntasks
         my_num_vars = my_num_vars + 1
         my_vars(my_num_vars) = int(i, i8)
      enddo

      ! this task holds only its own elements, as they were read in
      do i = 1, my_num_vars
         mine(:, i) = initial(:, my_vars(i))
      enddo

      call perturb_bitwise(mine(:, 1:my_num_vars), ENS_SIZE, &
                           my_vars(1:my_num_vars), int(NUM_VARS, i8), AMPLITUDE)

      do i = 1, my_num_vars
         do j = 1, ENS_SIZE
            if (mine(j, i) /= reference(j, my_vars(i))) matches = .false.
         enddo
      enddo

   enddo

   write(label, '(A,I2,A)') 'bitwise: ', ntasks, ' tasks give the 1 task answer, bit for bit'
   call ok(matches, trim(label))

enddo

end subroutine test_bitwise_task_invariance

!------------------------------------------------------------------
!> With more tasks than state elements a task can own nothing.  Indexing
!> copies with the local index would run off the end of a zero size array.

subroutine test_no_vars_on_this_task()

real(r8)    :: copies(NUM_COPIES, 0)
integer(i8) :: my_vars(0)

call perturb_bitwise(copies, ENS_SIZE, my_vars, int(NUM_VARS, i8), AMPLITUDE)
call pass('bitwise: a task owning no state returns without faulting')

call perturb_uniform(copies, ENS_SIZE, AMPLITUDE)
call pass('uniform: a task owning no state returns without faulting')

end subroutine test_no_vars_on_this_task

!------------------------------------------------------------------

end program perturb_test

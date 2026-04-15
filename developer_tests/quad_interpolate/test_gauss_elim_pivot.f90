! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
!> Unit tests for the gauss_elim_pivot subroutine.
!>

program test_gauss_elim_pivot

use types_mod,     only : r8, MISSING_R8
use utilities_mod, only : initialize_utilities, finalize_utilities, &
                          error_handler, E_MSG, E_ERR
use quad_utils_mod, only : gauss_elim_pivot

implicit none

! -----------------------------------------------------------------------
! Counters
! -----------------------------------------------------------------------
integer :: tests_run    = 0
integer :: tests_passed = 0

! -----------------------------------------------------------------------
! Working variables shared across tests
! -----------------------------------------------------------------------
real(r8) :: A(3,3), b(3), x(3), residual(3)
integer  :: stat
real(r8), parameter :: TOL = 1.0e-10_r8   ! acceptable residual threshold

! -----------------------------------------------------------------------
call initialize_utilities('test_gauss_elim_pivot')


! ===================================================================== !
! Test 1: diagonal system – no pivoting required                        !
! [ 2  0  0 ] [x]   [ 4]      solution: x = [2, 3, 4]                  !
! [ 0  3  0 ] [y] = [ 9]                                                !
! [ 0  0  4 ] [z]   [16]                                                !
! ===================================================================== !
A = reshape([2.0_r8, 0.0_r8, 0.0_r8, &
             0.0_r8, 3.0_r8, 0.0_r8, &
             0.0_r8, 0.0_r8, 4.0_r8], shape=[3,3], order=[2,1])
b = [4.0_r8, 9.0_r8, 16.0_r8]
call gauss_elim_pivot(A, b, x, stat)
call check_solution('Test 1 (diagonal)', A, b, x, stat, &
                    [2.0_r8, 3.0_r8, 4.0_r8], expect_fail=.false.)


! ===================================================================== !
! Test 2: full system, no pivoting needed                               !
! [ 2  1 -1 ] [x]   [ 8]      solution: x = [2, 3, -1]                 !
! [-3 -1  2 ] [y] = [-11]                                               !
! [-2  1  2 ] [z]   [-3]                                                !
! ===================================================================== !
A = reshape([ 2.0_r8,  1.0_r8, -1.0_r8, &
             -3.0_r8, -1.0_r8,  2.0_r8, &
             -2.0_r8,  1.0_r8,  2.0_r8], shape=[3,3], order=[2,1])
b = [8.0_r8, -11.0_r8, -3.0_r8]
call gauss_elim_pivot(A, b, x, stat)
call check_solution('Test 2 (full, no pivot)', A, b, x, stat, &
                    [2.0_r8, 3.0_r8, -1.0_r8], expect_fail=.false.)


! ===================================================================== !
! Test 3: system that requires row swaps (partial pivoting)             !
! [ 0  2  1 ] [x]   [ 5]      solution: x = [0, 2, 1]                  !
! [ 1  1 -1 ] [y] = [ 1]                                                !
! [ 1  3  2 ] [z]   [ 8]                                                !
! ===================================================================== !
A = reshape([0.0_r8, 2.0_r8,  1.0_r8, &
             1.0_r8, 1.0_r8, -1.0_r8, &
             1.0_r8, 3.0_r8,  2.0_r8], shape=[3,3], order=[2,1])
b = [5.0_r8, 1.0_r8, 8.0_r8]
call gauss_elim_pivot(A, b, x, stat)
call check_solution('Test 3 (row swap required)', A, b, x, stat, &
                    [0.0_r8, 2.0_r8, 1.0_r8], expect_fail=.false.)


! ===================================================================== !
! Test 4: residual check on a random-looking non-trivial system         !
! [ 3  1  2 ] [x]   [11]      solution: x = [1, 2, 3]                  !
! [ 1  4  3 ] [y] = [18]                                                !
! [ 2  2  5 ] [z]   [21]                                                !
! ===================================================================== !
A = reshape([3.0_r8, 1.0_r8, 2.0_r8, &
             1.0_r8, 4.0_r8, 3.0_r8, &
             2.0_r8, 2.0_r8, 5.0_r8], shape=[3,3], order=[2,1])
b = [11.0_r8, 18.0_r8, 21.0_r8]
call gauss_elim_pivot(A, b, x, stat)
call check_solution('Test 4 (residual check)', A, b, x, stat, &
                    [1.0_r8, 2.0_r8, 3.0_r8], expect_fail=.false.)


! ===================================================================== !
! Test 5: singular matrix (two identical rows) → stat must be 1         !
! [ 1  2  3 ]                                                           !
! [ 1  2  3 ]   (rows 1 and 2 identical)                                !
! [ 0  1  1 ]                                                           !
! ===================================================================== !
A = reshape([1.0_r8, 2.0_r8, 3.0_r8, &
             1.0_r8, 2.0_r8, 3.0_r8, &
             0.0_r8, 1.0_r8, 1.0_r8], shape=[3,3], order=[2,1])
b = [6.0_r8, 6.0_r8, 2.0_r8]
call gauss_elim_pivot(A, b, x, stat)
call check_singular('Test 5 (singular – identical rows)', stat)


! ===================================================================== !
! Test 6: all-zero matrix → stat must be 1                              !
! ===================================================================== !
A = 0.0_r8
b = [1.0_r8, 2.0_r8, 3.0_r8]
call gauss_elim_pivot(A, b, x, stat)
call check_singular('Test 6 (all-zero matrix)', stat)


! ===================================================================== !
! Test 7: rank-deficient (third row is sum of first two) → stat = 1    !
! [ 1  0  0 ]                                                           !
! [ 0  1  0 ]                                                           !
! [ 1  1  0 ]  (row3 = row1 + row2, last column all-zero)              !
! ===================================================================== !
A = reshape([1.0_r8, 0.0_r8, 0.0_r8, &
             0.0_r8, 1.0_r8, 0.0_r8, &
             1.0_r8, 1.0_r8, 0.0_r8], shape=[3,3], order=[2,1])
b = [1.0_r8, 2.0_r8, 3.0_r8]
call gauss_elim_pivot(A, b, x, stat)
call check_singular('Test 7 (rank-deficient)', stat)


! ===================================================================== !
! Test 8: negative and mixed-sign values – checks correctness and sign  !
! [-1  2  3 ] [x]   [-1]      solution: x = [1, 1, -1]                 !
! [ 4 -2  1 ] [y] = [  1]                                               !
! [ 2  3 -5 ] [z]   [ -2]                                               !
! ===================================================================== !
A = reshape([-1.0_r8,  2.0_r8,  3.0_r8, &
              4.0_r8, -2.0_r8,  1.0_r8, &
              2.0_r8,  3.0_r8, -5.0_r8], shape=[3,3], order=[2,1])
b = [-1.0_r8 + 2.0_r8 - 3.0_r8, &   ! = -2  ... compute explicitly below
      4.0_r8 - 2.0_r8 - 1.0_r8, &   ! =  1
      2.0_r8 + 3.0_r8 + 5.0_r8]     ! = 10  ... wait, correct by hand

! b = A * [1, 1, -1]:
!  row1: -1*1 + 2*1 + 3*(-1) = -1+2-3 = -2
!  row2:  4*1 - 2*1 + 1*(-1) =  4-2-1 =  1
!  row3:  2*1 + 3*1 - 5*(-1) =  2+3+5 = 10
b = [-2.0_r8, 1.0_r8, 10.0_r8]
call gauss_elim_pivot(A, b, x, stat)
call check_solution('Test 8 (mixed signs)', A, b, x, stat, &
                    [1.0_r8, 1.0_r8, -1.0_r8], expect_fail=.false.)


! ===================================================================== !
! Test 9: well-scaled system where pivot choice matters for accuracy    !
! Large off-diagonal entry forces a pivot swap                          !
! [ 0.001  1.0   0.0 ] [x]   [1.001]     solution: x = [1, 1, 1]       !
! [ 1000.0  0.0  1.0 ] [y] = [1001.0]                                   !
! [ 0.0    1.0   1.0 ] [z]   [2.0]                                      !
! ===================================================================== !
A = reshape([  0.001_r8, 1.0_r8  , 0.0_r8, &
            1000.0_r8  , 0.0_r8  , 1.0_r8, &
               0.0_r8  , 1.0_r8  , 1.0_r8], shape=[3,3], order=[2,1])
! b = A * [1, 1, 1]:
!  row1: 0.001+1+0  = 1.001
!  row2: 1000+0+1   = 1001
!  row3: 0+1+1      = 2
b = [1.001_r8, 1001.0_r8, 2.0_r8]
call gauss_elim_pivot(A, b, x, stat)
call check_solution('Test 9 (large off-diagonal – pivoting accuracy)', A, b, x, stat, &
                    [1.0_r8, 1.0_r8, 1.0_r8], expect_fail=.false.)


! -----------------------------------------------------------------------
! Summary
! -----------------------------------------------------------------------
write(*, '(A)')         ''
write(*, '(A,I0,A,I0)') 'Tests passed: ', tests_passed, ' / ', tests_run
if (tests_passed == tests_run) then
   write(*, '(A)') 'All tests PASSED.'
else
   write(*, *) tests_run - tests_passed, ' test(s) FAILED.'
   write(*, '(A)') 'One or more tests failed.'
endif

call finalize_utilities('test_gauss_elim_pivot')


! ======================================================================
contains
! ======================================================================

!----------------------------------------------------------------------
subroutine check_solution(label, A_orig, b_orig, x, stat, x_expected, expect_fail)
!> Verify that:
!>  (a) stat == 0  (routine considers system non-singular)
!>  (b) the returned x equals x_expected within TOL
!>  (c) the residual A*x - b is within TOL element-wise

character(len=*), intent(in) :: label
real(r8),         intent(in) :: A_orig(3,3), b_orig(3), x(3), x_expected(3)
integer,          intent(in) :: stat
logical,          intent(in) :: expect_fail

real(r8) :: residual(3)
integer  :: i
logical  :: ok

tests_run = tests_run + 1
ok = .true.

if (stat /= 0) then
   write(*, '(2A)') 'FAIL ', label
   write(*, '(A,I0)') '  stat = ', stat, ' (expected 0)'
   ok = .false.
endif

if (ok) then
   ! solution accuracy
   do i = 1, 3
      if (abs(x(i) - x_expected(i)) > TOL) then
         if (ok) write(*, '(2A)') 'FAIL ', label
         write(*, '(A,I0,3(A,ES12.5))') '  x(', i, ') = ', x(i), &
            '  expected ', x_expected(i), '  diff ', abs(x(i) - x_expected(i))
         ok = .false.
      endif
   enddo
endif

if (ok) then
   ! residual check: r = A_orig * x - b_orig
   do i = 1, 3
      residual(i) = A_orig(i,1)*x(1) + A_orig(i,2)*x(2) + A_orig(i,3)*x(3) - b_orig(i)
      if (abs(residual(i)) > TOL) then
         if (ok) write(*, '(2A)') 'FAIL ', label
         write(*, '(A,I0,A,ES12.5)') '  residual(', i, ') = ', residual(i)
         ok = .false.
      endif
   enddo
endif

if (ok) then
   tests_passed = tests_passed + 1
   write(*, '(2A)') 'PASS ', label
endif

end subroutine check_solution


!----------------------------------------------------------------------
subroutine check_singular(label, stat)
!> Verify that stat /= 0 for a singular or rank-deficient system.

character(len=*), intent(in) :: label
integer,          intent(in) :: stat

tests_run = tests_run + 1

if (stat /= 0) then
   tests_passed = tests_passed + 1
   write(*, '(2A)') 'PASS ', label
else
   write(*, '(2A)') 'FAIL ', label
   write(*, '(A)')  '  expected stat /= 0 for singular system, got stat = 0'
endif

end subroutine check_singular

end program test_gauss_elim_pivot

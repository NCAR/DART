! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_corr

! test of twod_gaussian distribution, and two different numerical
! methods of computing the variance and correlation.  one is more
! stable in the face of 8 byte reals being redefined as 4 byte values.

use      types_mod, only : r8
use  utilities_mod, only : check_namelist_read, find_namelist_in_file, &
                           register_module, error_handler, E_ERR,      &
                           nmlfileunit, do_nml_file, do_nml_term,      &
                           initialize_utilities, finalize_utilities
use random_seq_mod, only : random_seq_type, init_random_seq, twod_gaussians

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (random_seq_type) :: r
real(r8), allocatable :: rnum(:, :)
real(r8) :: c(2, 2), mean(2)
integer :: i, j, iunit, io

real(r8) :: mean_correlA, varA, meanA, meanc_2A, sample_correlA
real(r8) :: mean_correlB, varB, meanB, meanc_2B, sample_correlB

! namelist variables
integer  :: n_pairs = 10
integer  :: n_cycles = 100
real(r8) :: desired_correl = 0.75_r8

namelist /test_corr_nml/ n_pairs, desired_correl, n_cycles


! start of program

call initialize_utilities('test_corr')
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "test_corr_nml", iunit)
read(iunit, nml = test_corr_nml, iostat = io)
call check_namelist_read(iunit, io, "test_corr_nml")

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=test_corr_nml)
if (do_nml_term()) write(     *     , nml=test_corr_nml)


if (n_cycles < 2 .or. n_pairs < 2) then
   call error_handler(E_ERR, 'test_corr', 'both n_pairs and n_cycles must be >= 2', &
      source, revision, revdate)
endif

mean_correlA = 0.0
mean_correlB = 0.0
varA = 0.0
varB = 0.0

! Allocate space for the pairs of data
allocate(rnum(2, n_pairs))

call init_random_seq(r)

! Mean for generated
mean = 0.0
meanc_2A = 0.0
meanc_2B = 0.0

! Generate the covariance matrix
c(1, 1) = 1.0
c(2, 2) = 1.0
c(1, 2) = desired_correl
c(2, 1) = desired_correl

do j = 1, n_cycles
   
   do i = 1, n_pairs
      call twod_gaussians(r, mean, c, rnum(:, i))
   end do

   ! Compute the correlation of the random numbers two different ways
   call comp_correlA(rnum, n_pairs, sample_correlA)
   call comp_correlB(rnum, n_pairs, sample_correlB)

   mean_correlA = mean_correlA + sample_correlA
   meanc_2A = meanc_2A + sample_correlA**2
   varA = varA + (sample_correlA - desired_correl) ** 2

   mean_correlB = mean_correlB + sample_correlB
   meanc_2B = meanc_2B + sample_correlB**2
   varB = varB + (sample_correlB - desired_correl) ** 2
end do

mean_correlA = mean_correlA / n_cycles
mean_correlB = mean_correlB / n_cycles

write(*, *) ''
write(*, *) 'target mean correl is                ', desired_correl
write(*, *) 'computed mean correl, version A, is  ', mean_correlA
write(*, *) 'computed mean correl, version B, is  ', mean_correlB

write(*, *) ''
write(*, *) 'predicted variance is                ', (1.0_r8 - desired_correl**2)**2 / (n_pairs - 1._r8)
write(*, *) 'computed variance, version A, is     ', (meanc_2A - n_cycles * mean_correlA**2) / (n_cycles - 1._r8)
write(*, *) 'computed variance, version B, is     ', (meanc_2B - n_cycles * mean_correlB**2) / (n_cycles - 1._r8)

call finalize_utilities('test_corr')

contains

!-----------------------------------------------------
 
subroutine comp_correlA(ens, n, correl)
 
integer,  intent(in)  :: n
real(r8), intent(in)  :: ens(2, n)
real(r8), intent(out) :: correl

real(r8) :: sum_x, sum_y, sum_xy, sum_x2, sum_y2

sum_x = sum(ens(2, :))
sum_y = sum(ens(1, :))
sum_xy = sum(ens(2, :) * ens(1, :))
sum_x2 = sum(ens(2, :) * ens(2, :))
 
! Computation of correlation 
sum_y2 = sum(ens(1, :) * ens(1, :))

! debug:
!print *, 'x,y,xy,x2,y2: ', sum_x, sum_y, sum_xy, sum_x2, sum_y2

correl = (n * sum_xy - sum_x * sum_y) / &
   sqrt((n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2))
 
end subroutine comp_correlA

!-----------------------------------------------------
 
subroutine comp_correlB(ens, n, correl)
 
! alternative way to comput covariance and variance - more numerically
! stable in the face of 4 byte floats:

integer,  intent(in)  :: n
real(r8), intent(in)  :: ens(2, n)
real(r8), intent(out) :: correl

real(r8) :: sum_x, sum_y
real(r8) :: x_mean, y_mean, x_var, y_var, cov

sum_x = sum(ens(2, :))
sum_y = sum(ens(1, :))
 
x_mean = sum_x / n
y_mean = sum_y / n
x_var = sum((ens(2, :) - x_mean)**2) / (n - 1)
y_var = sum((ens(1, :) - y_mean)**2) / (n - 1)
cov = sum((ens(2, :) - x_mean) * (ens(1, :) - y_mean)) / (n - 1)
correl = cov / sqrt(x_var * y_var)

! debug:
!print *, x_mean, y_mean, x_var, y_var, cov, abs(x_var * y_var), sqrt(x_var*y_var), correl

end subroutine comp_correlB

!-----------------------------------------------------
 
end program test_corr

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

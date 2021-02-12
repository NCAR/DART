! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Correct covariances for fixed ensemble sizes.
!> See Anderson, J. L., 2011: Localization and Sampling Error Correction
!>   in Ensemble Kalman Filter Data Assimilation. 
!> Submitted for publication, Jan 2011.  Contact author.

program full_error

! This version of the program reads the ensemble size and base filename
! for the output from a namelist.

use      types_mod, only : r8
use  utilities_mod, only : open_file, close_file, error_handler,       &
                           find_namelist_in_file, check_namelist_read, &
                           do_nml_file, do_nml_term, nmlfileunit,      &
                           initialize_utilities, finalize_utilities, E_ERR
use random_seq_mod, only : random_seq_type, init_random_seq, &          
                           twod_gaussians

implicit none

character(len=*), parameter :: source = 'full_error.f90'

integer, parameter :: num_times   = 1
integer, parameter :: num_samples = 100000000

! ---------------
! namelist items

integer            :: ens_size = 80
character(len=256) :: output_filename = 'final_full'

namelist /full_error_nml/ ens_size, output_filename


type (random_seq_type) :: ran_id
real(r8), allocatable  :: pairs(:,:), temp(:)

real(r8) :: zero_2(2) = 0.0, cov(2, 2)
real(r8) :: t_correl, sample_correl, alpha, beta
real(r8) :: s_mean(2), s_var(2), reg_mean, reg_sd, t_sd1, t_sd2, true_correl_mean
real(r8) :: tcorrel_sum(201), reg_sum(201), reg_2_sum(201)

integer  :: i, j, k, bin_num, bin_count(201)
integer  :: iunit, io

character(len=16)  :: formstring
character(len=256) :: outname
character(len=512) :: errstring

!
! start of executable code
!

call initialize_utilities('full_error')

! Read the namelist entry
call find_namelist_in_file("input.nml", "full_error_nml", iunit)
read(iunit, nml = full_error_nml, iostat = io)
call check_namelist_read(iunit, io, "full_error_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=full_error_nml)
if (do_nml_term()) write(     *     , nml=full_error_nml)

call close_file(iunit)

call init_random_seq(ran_id)

write(*, *) 'stats for ensemble size ', ens_size

allocate(pairs(2, ens_size), temp(ens_size))

bin_count   = 0
tcorrel_sum = 0.0_r8
reg_sum     = 0.0_r8
reg_2_sum   = 0.0_r8

! Uniformly distributed sequence of true correlations  
do j = 1, num_samples
   ! Assume that true correlation is uniform [-1,1]
   t_correl = -1.0_r8 + 2.0_r8 * ((j - 1) / (num_samples - 1.0_r8))
   if(j > num_samples) t_correl = 0.0_r8

   do k = 1, num_times

      t_sd1 = 1.0_r8
      t_sd2 = 1.0_r8

      ! Generate the covariance matrix for this correlation
      cov(1, 1) = t_sd1**2
      cov(2, 2) = t_sd2**2
      cov(1, 2) = t_correl * t_sd1 * t_sd2
      cov(2, 1) = cov(1, 2)
   
      ! Loop to generate an ensemble size sample from this correl
      ! Generate a random sample of size ens_size from something with this correlation
      do i = 1, ens_size
         call twod_gaussians(ran_id, zero_2, cov, pairs(:, i))
      enddo
   
      ! Compute the sample correlation
      call comp_correl(pairs, ens_size, sample_correl)

      ! Bin statistics for each 0.01 sample correlation bin; 
      !  start with just mean true correlation
      bin_num = (sample_correl + 1.0_r8) / 0.01_r8   + 1
      ! Keep a sum of the sample_correl and squared to compute standard deviation
      tcorrel_sum(bin_num) = tcorrel_sum(bin_num) + t_correl
      !-----------------
      ! Also interested in finding out what the spurious variance reduction factor is
      ! First need to compute sample s.d. for the obs and unobs variable
      do i = 1, 2
         temp = pairs(i, :)
         call sample_mean_var(temp, ens_size, s_mean(i), s_var(i))
      enddo

      !-----------------
      reg_sum(bin_num) = reg_sum(bin_num) + t_correl * sqrt(s_var(2) / s_var(1))
      reg_2_sum(bin_num) = reg_2_sum(bin_num) + (t_correl * sqrt(s_var(2) / s_var(1)))**2

      bin_count(bin_num) = bin_count(bin_num) + 1

   enddo
enddo

! make the size of the integer output .X, .XX, .XXX, etc
! depending on the number of decimal digits in the value.
if (ens_size < 10) then
   formstring = '(2A,I1)'
else if (ens_size < 100) then
   formstring = '(2A,I2)'
else if (ens_size < 1000) then
   formstring = '(2A,I3)'
else if (ens_size < 10000) then
   formstring = '(2A,I4)'
else if (ens_size < 100000) then
   formstring = '(2A,I5)'
else 
   formstring = '(2A,I8)'
endif

! filename
write(outname, formstring) trim(output_filename), '.', ens_size

! text file, overwrite existing file if present
iunit = open_file(outname, 'formatted', 'write')

! always generate exactly 200 entries in the output file
do i = 1, 200
   ! must have at least 2 counts to compute a std dev
   if(bin_count(i) <= 1) then
      write(errstring, *) 'Bin ', i, ' has ', bin_count(i), ' counts'
      call error_handler(E_ERR,'full_error', errstring, &
         source, text2="All bins must have at least 2 counts")
   endif
   
   ! Compute the standard deviation of the true correlations
   true_correl_mean = tcorrel_sum(i) / bin_count(i)
   reg_mean = reg_sum(i) / bin_count(i)
   reg_sd = sqrt((reg_2_sum(i) - bin_count(i) * reg_mean**2) / (bin_count(i) - 1))

   if(reg_sd <= 0.0_r8) then
      alpha = 1.0_r8
   else
      !!!beta = reg_mean**2 / reg_sd**2
   ! Correct for bias in the standard deviation for very small ensembles, too 
      beta = reg_mean**2 / (reg_sd**2 * (1.0_r8 + 1.0_r8 / ens_size))
      alpha = beta / (1.0_r8 + beta)
   endif

   write(*, *) 'bin, count, mean ', i, bin_count(i), true_correl_mean, alpha
   write(iunit, 10)   i, bin_count(i), true_correl_mean, alpha
enddo

10 format (I4,I9,2G25.14)

call close_file(iunit)

! print out just to stdout how many counts, if any, fell beyond bin 200
write(*, *) 'bin, count, mean ', 201, bin_count(201), 0, 0

deallocate(pairs, temp)

call finalize_utilities()

! end of main program

contains

!-----------------------------------------------------

subroutine comp_correl(ens, n, correl)

integer,  intent(in)  :: n
real(r8), intent(in)  :: ens(2, n)
real(r8), intent(out) :: correl

real(r8) :: sum_x, sum_y, sum_xy, sum_x2, sum_y2


sum_x  = sum(ens(2, :))
sum_y  = sum(ens(1, :))
sum_xy = sum(ens(2, :) * ens(1, :))
sum_x2 = sum(ens(2, :) * ens(2, :))

! Computation of correlation
sum_y2 = sum(ens(1, :) * ens(1, :))

correl = (n * sum_xy - sum_x * sum_y) / &
   sqrt((n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2))

end subroutine comp_correl

!----------------------------------------------------------------

subroutine sample_mean_var(x, n, mean, var)

integer,  intent(in)  :: n
real(r8), intent(in)  :: x(n)
real(r8), intent(out) :: mean, var

real(r8) :: sx, s_x2

sx   = sum(x)
s_x2 = sum(x * x)
mean = sx / n
var  = (s_x2 - sx**2 / n) / (n - 1)

end subroutine sample_mean_var


end program full_error


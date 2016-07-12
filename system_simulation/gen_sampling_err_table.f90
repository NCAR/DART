! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Correct covariances for fixed ensemble sizes.
!> See Anderson, J. L., 2011: Localization and Sampling Error Correction
!>   in Ensemble Kalman Filter Data Assimilation. 
!> Submitted for publication, Jan 2011.  Contact author.

!> this version of the program creates entries for any ensemble size
!> from 2 to 1000 and outputs the values in netcdf file format.
!> it has no namelist - the output file is named 'sampling_error_correction_table.nc'

program gen_sampling_err_table

use types_mod,      only : r8
use utilities_mod,  only : error_handler, E_ERR, nc_check,      &
                           initialize_utilities, finalize_utilities
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian, &          
                           twod_gaussians, random_uniform

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


integer, parameter :: num_times   = 1
integer, parameter :: num_samples = 100000000

integer, parameter :: min_ens_size = 4
!!!!!integer, parameter :: max_ens_size = 2000
integer, parameter :: max_ens_size = 10

integer, parameter :: nentries = 200

character(len=128) :: output_filename = 'sampling_error_correction_table.nc'

type (random_seq_type) :: ran_id
real(r8), allocatable  :: pairs(:,:), temp(:)

real(r8) :: zero_2(2) = 0.0, cov(2, 2)
real(r8) :: t_correl, correl_mean, sample_correl, alpha(nentries), beta
real(r8) :: s_mean(2), s_var(2), reg_mean, reg_sd, t_sd1, t_sd2, true_correl_mean(nentries)
real(r8) :: tcorrel_sum(nentries), reg_sum(nentries), reg_2_sum(nentries)

integer  :: i, j, k, bin_num, bin_count(0:nentries+1)
integer  :: iunit, io, ncid, ens_size
integer  :: count_id, corrmean_id, alpha_id

character(len=512) :: errstring, msgstring

!
! start of executable code
!

call initialize_utilities('gen_sampling_err_table')

call init_random_seq(ran_id)

ncid = open_output_file()

call setup_output_file(ncid, count_id, corrmean_id, alpha_id)

do ens_size=min_ens_size, max_ens_size

   write(*,*) 'computing for ensemble size of ', ens_size
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
         end do
      
         ! Compute the sample correlation
         call comp_correl(pairs, ens_size, sample_correl)
   
         ! Bin statistics for each 0.01 sample correlation bin; 
         !  start with just mean true correlation
         bin_num = (sample_correl + 1.0_r8) / 0.01_r8   + 1
         if (bin_num < 1) then
            bin_count(0) = bin_count(0) + 1
            cycle
         endif
         if (bin_num > nentries) then
            bin_count(nentries+1) = bin_count(nentries+1) + 1
            cycle
         endif

         ! Keep a sum of the sample_correl and squared to compute standard deviation
         tcorrel_sum(bin_num) = tcorrel_sum(bin_num) + t_correl
         !-----------------
         ! Also interested in finding out what the spurious variance reduction factor is
         ! First need to compute sample s.d. for the obs and unobs variable
         do i = 1, 2
            temp = pairs(i, :)
            call sample_mean_var(temp, ens_size, s_mean(i), s_var(i))
         end do
   
         !-----------------
         reg_sum(bin_num) = reg_sum(bin_num) + t_correl * sqrt(s_var(2) / s_var(1))
         reg_2_sum(bin_num) = reg_2_sum(bin_num) + (t_correl * sqrt(s_var(2) / s_var(1)))**2
   
         bin_count(bin_num) = bin_count(bin_num) + 1
   
      end do
   end do
   
   ! print out just to stdout how many counts, if any, fell below bin 0
   write(*, *) 'bin, count, mean ', 0, bin_count(0), 0, 0

   ! always generate exactly 'nentries' entries in the output file
   do i = 1, nentries
      ! must have at least 2 counts to compute a std dev
      if(bin_count(i) <= 1) then
         write(*, *) 'Bin ', i, ' has ', bin_count(i), ' counts'
         cycle
         write(errstring, *) 'Bin ', i, ' has ', bin_count(i), ' counts'
         call error_handler(E_ERR,'gen_sampling_err_table', errstring, &
            source, revision, revdate, text2='All bins must have at least 2 counts')
      endif
      
      ! Compute the standard deviation of the true correlations
      true_correl_mean(i) = tcorrel_sum(i) / bin_count(i)
      reg_mean = reg_sum(i) / bin_count(i)
      reg_sd = sqrt((reg_2_sum(i) - bin_count(i) * reg_mean**2) / (bin_count(i) - 1))
   
      if(reg_sd <= 0.0_r8) then
         alpha(i) = 1.0_r8
      else
         !!!beta = reg_mean**2 / reg_sd**2
      ! Correct for bias in the standard deviation for very small ensembles, too 
         beta = reg_mean**2 / (reg_sd**2 * (1.0_r8 + 1.0_r8 / ens_size))
         alpha(i) = beta / (1.0_r8 + beta)
      endif
   
      write(*, *) 'bin, count, mean ', i, bin_count(i), true_correl_mean(i), alpha(i)
!      write(iunit, 10)   i, bin_count(i), true_correl_mean(i), alpha(i)
10 format (I4,I9,2G25.14)
   end do
   
   ! print out just to stdout how many counts, if any, fell beyond bin 'nentries'
   write(*, *) 'bin, count, mean ', nentries+1, bin_count(nentries+1), 0, 0
   
   call write_output_file(ncid, ens_size, count_id, bin_count(1:nentries), &
                                          corrmean_id, true_correl_mean, &
                                          alpha_id, alpha)
   deallocate(pairs, temp)
   
enddo  ! ens_size

call close_output_file(ncid)

call finalize_utilities()

! end of main program

contains

!----------------------------------------------------------------

!----------------------------------------------------------------
! stat routines
!----------------------------------------------------------------

!----------------------------------------------------------------

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

!----------------------------------------------------------------

!----------------------------------------------------------------
! main netcdf i/o routines
!----------------------------------------------------------------

!----------------------------------------------------------------

!> create a netcdf file and add some global attributes

function open_output_file()

integer :: open_output_file

integer :: rc, ncid

rc = nf90_create(output_filename, NF90_CLOBBER, ncid)
call nc_check(rc, 'open_output_file', 'creating '//trim(output_filename))

write(msgstring, *) 'supported ensemble sizes from ', min_ens_size, ' to ', max_ens_size
rc = nf90_put_att(ncid, NF90_GLOBAL, 'ensemble_sizes', trim(msgstring))
call nc_check(rc, 'open_output_file', 'adding global attributes')

write(msgstring, *) 'each ensemble size entry has ', nentries, ' bins'
rc = nf90_put_att(ncid, NF90_GLOBAL, 'bin_info', trim(msgstring))
call nc_check(rc, 'open_output_file', 'adding global attributes')

rc = nf90_put_att(ncid, NF90_GLOBAL, 'min_ens_size', min_ens_size)
call nc_check(rc, 'open_output_file', 'adding global attributes')

! you can get the max from the dimension declaration

open_output_file = ncid

end function open_output_file

!----------------------------------------------------------------

! define 2 dims and 3 arrays in output file
! counts on knowing global info
! returns 3 handles to use for the 3 arrays

subroutine setup_output_file(ncid, id1, id2, id3)

integer, intent(in)  :: ncid
integer, intent(out) :: id1, id2, id3

integer :: rc
integer :: nbinsDimID, nensDimID

call setup_sec_dim(ncid, 'bins', nentries,     nbinsDimID)
call setup_sec_dim(ncid, 'ens',  max_ens_size, nensDimID)

call setup_sec_data_int (ncid, 'count',          nbinsDimID, nensDimID, id1)
call setup_sec_data_real(ncid, 'true_corr_mean', nbinsDimID, nensDimID, id2)
call setup_sec_data_real(ncid, 'alpha',          nbinsDimID, nensDimID, id3)

rc = nf90_enddef(ncid)
call nc_check(rc, 'setup_output_file', 'enddef')

end subroutine setup_output_file

!----------------------------------------------------------------

! write 3 arrays to file

subroutine write_output_file(ncid, col, id1, a1, id2, a2, id3, a3)

integer,  intent(in)  :: ncid
integer,  intent(in)  :: col
integer,  intent(in)  :: id1
integer,  intent(in)  :: a1(:)
integer,  intent(in)  :: id2
real(r8), intent(in)  :: a2(:)
integer,  intent(in)  :: id3
real(r8), intent(in)  :: a3(:)

call write_sec_data_int (ncid, col, 'count',          id1, a1)
call write_sec_data_real(ncid, col, 'true_corr_mean', id2, a2)
call write_sec_data_real(ncid, col, 'alpha',          id3, a3)

end subroutine write_output_file

!----------------------------------------------------------------

subroutine close_output_file(ncid)

integer, intent(in) :: ncid

integer :: rc

rc = nf90_sync(ncid)
call nc_check(rc, 'close_output_file', 'syncing '//trim(output_filename))

rc = nf90_close(ncid)
call nc_check(rc, 'close_output_file', 'closing '//trim(output_filename))

end subroutine close_output_file

!----------------------------------------------------------------

!----------------------------------------------------------------
! helper routines for above code
!----------------------------------------------------------------

!----------------------------------------------------------------

! define a dimension

subroutine setup_sec_dim(ncid, c1, n1, id1)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: c1
integer,          intent(in)  :: n1
integer,          intent(out) :: id1

integer :: rc

rc = nf90_def_dim(ncid, name=c1, len=n1, dimid=id1)
call nc_check(rc, 'setup_sec_dim', 'adding dimension '//trim(c1))

end subroutine setup_sec_dim

!----------------------------------------------------------------

! given a name, return id for a 2d integer array

subroutine setup_sec_data_int(ncid, c1, d1, d2, id1)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: c1
integer,          intent(in)  :: d1, d2
integer,          intent(out) :: id1

integer :: rc

rc = nf90_def_var(ncid, name=c1, xtype=nf90_int, dimids=(/ d1, d2 /), varid=id1)
call nc_check(rc, 'setup_sec_data', 'defining variable '//trim(c1))

end subroutine setup_sec_data_int

!----------------------------------------------------------------

! given a name, return id for a 2d real array

subroutine setup_sec_data_real(ncid, c1, d1, d2, id1)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: c1
integer,          intent(in)  :: d1, d2
integer,          intent(out) :: id1

integer :: rc

rc = nf90_def_var(ncid, name=c1, xtype=nf90_real, dimids=(/ d1, d2 /), varid=id1)
call nc_check(rc, 'setup_sec_data', 'defining variable '//trim(c1))

end subroutine setup_sec_data_real

!----------------------------------------------------------------

subroutine write_sec_data_int(ncid, col, c1, id1, a1)

integer,          intent(in) :: ncid
integer,          intent(in) :: col
character(len=*), intent(in) :: c1
integer,          intent(in) :: id1
integer,          intent(in) :: a1(:)

integer :: rc

rc = nf90_put_var(ncid, id1, a1, start=(/ 1, col /), count=(/ size(a1), 1 /) )
call nc_check(rc, 'write_sec_data', 'writing variable "'//trim(c1)//'"')

end subroutine write_sec_data_int

!----------------------------------------------------------------

subroutine write_sec_data_real(ncid, col, c1, id1, a1) 
integer,          intent(in) :: ncid
integer,          intent(in) :: col
character(len=*), intent(in) :: c1
integer,          intent(in) :: id1
real(r8),         intent(in) :: a1(:)

integer :: rc

rc = nf90_put_var(ncid, id1, a1, start=(/ 1, col /), count=(/ size(a1), 1 /) )
call nc_check(rc, 'write_sec_data', 'writing variable "'//trim(c1)//'"')

end subroutine write_sec_data_real

!----------------------------------------------------------------


end program gen_sampling_err_table

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

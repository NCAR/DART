module assim_diag_mod

! Does a variety of diagnostics to evaluate the performance of the assimilation
! algorithms. A defined type is used to keep track of the RMS errors, bias,
! and order statistic bins.

! Temporary version created on 13 October, 2000. Many shortcuts have been
! taken to allow an initial high resolution model version to begin running.
! Need to generalize this code. Note that multiple calls to output_diag
! must be in same order as calls to assim_diag_init.

use types_mod
use chisq_mod, only: chsone
use utilities_mod, only : open_file

implicit none

private

public bin_num, save_diagnostics, output_diagnostics, assim_diag_init

type, public :: diag_type
   private
   integer  :: count, model_size, ens_size, num_var
   real(r8) :: total_rms_mean, mean_total_rms, mean_bias
   real(r8) :: sum_rms_mean, sum_mean_rms
   integer,  pointer :: bin(:, :), var_list(:)
   real(r8), pointer :: rms_mean(:), mean_rms(:), bias(:), sq_bias(:)
end type diag_type

logical :: first = .true.
integer :: bin_unit, obs_bin_unit, read_chan = -1, write_chan

contains



  subroutine assim_diag_init(d, model_size, ens_size, var_list, num_var, add_on)
!========================================================================
! subroutine assim_diag_init(d, model_size, ens_size, var_list, num_var, add_on)
!

implicit none

integer,          intent(in) :: model_size, ens_size, num_var, var_list(num_var)
logical,          intent(in) :: add_on
type (diag_type), intent(inout) :: d

! First thing that should be done is to output version number to logfile

! Initialize the diag_type structure

d%count      = 0
d%model_size = model_size
d%ens_size   = ens_size
d%num_var    = num_var

allocate( d%bin(num_var, 0:ens_size), &
          d%var_list(num_var), d%rms_mean(num_var), &
          d%mean_rms(num_var), d%bias(num_var), d%sq_bias(num_var))

! Initialize the indices of the variables to have diagnostics

d%var_list = var_list

! Initialize all accumulating diagnostic quantities to zero

d%bin            = 0
d%total_rms_mean = 0.0_r8
d%mean_total_rms = 0.0_r8
d%mean_bias      = 0.0_r8
d%sum_rms_mean   = 0.0_r8
d%sum_mean_rms   = 0.0_r8
d%rms_mean       = 0.0_r8
d%mean_rms       = 0.0_r8
d%bias           = 0.0_r8
d%sq_bias        = 0.0_r8

! See if a diagnostic restart should be read in and hasn't already been opened

if(add_on) then

   if(read_chan < 0) &
      read_chan = open_file(file = 'diag_restart.in', action = 'read')

   read(read_chan, *) d%count, d%model_size, d%ens_size, d%num_var
   read(read_chan, *) d%total_rms_mean, d%mean_total_rms, d%mean_bias
   read(read_chan, *) d%sum_rms_mean, d%sum_mean_rms, d%bin, d%var_list
   read(read_chan, *) d%rms_mean, d%mean_rms, d%bias, d%sq_bias
end if
   
end subroutine assim_diag_init



  subroutine save_diagnostics(x, ens_mean, truth, d)
!========================================================================
! subroutine save_diagnostics(x, ens_mean, truth, d)
!
! Saves information about the assimilation given the ensemble assimilated
! estimate, x, the corresponding observations, obs, and the truth. A 
! diag_type, d, is used to store all this.
!
! x(model_size, ens_size):	Ensemble of assimilated states
! ens_mean(model_size):		Ensemble mean of x
! truth(model_size):		True state of system
! d:				Diag_type to accumulate diagnostics
!------------------------------------------------------------------------

implicit none

real(r8),         intent(in)    :: x(:, :), truth(:), ens_mean(:)
type (diag_type), intent(inout) :: d

real(r8) :: squared_dif
integer  :: i, j, ind, num


d%count = d%count + 1         ! Increment the number of obs times

! Compute the total rms distances (not components)

d%total_rms_mean = d%total_rms_mean + sqrt(sum((ens_mean - truth)**2))

do i = 1, d%ens_size
   d%mean_total_rms = d%mean_total_rms + sqrt(sum((x(:, i) - truth)**2))
end do

! Compute the rms distance of individual components

do i = 1, d%num_var

   ind = d%var_list(i)
   d%rms_mean(i) = d%rms_mean(i) + abs(ens_mean(ind) - truth(ind))
   d%bias(i)     = d%bias(i)     +     ens_mean(ind) - truth(ind)

   ! Interested in relative square departure after correction of bias (see
   ! Notes from 31 August, 2001) relative to prediction from ensemble.

   squared_dif = 0.0_r8
   do j = 1, d%ens_size
      squared_dif = squared_dif + (x(ind, j) - ens_mean(ind))**2
   end do
   squared_dif = squared_dif / d%ens_size
   d%sq_bias(i) = d%sq_bias(i) + (ens_mean(ind) - truth(ind))**2 / squared_dif

   do j = 1, d%ens_size
      d%mean_rms(i) = d%mean_rms(i) + dabs(x(ind, j) - truth(ind))
   end do
end do
 
! Compute mean of bias for all variables (meaningful?)

d%mean_bias = d%mean_bias + sum(ens_mean - truth)

! 4 April, 2001: AN ERROR IN EARLIER VERSIONS LED TO AN INCORRECT
! COMPUTATION OF THE RMS ERROR FOR SUMMARY PURPOSES. THIS DOES NOT
! APPEAR TO HAVE BEEN USED ANYWHERE.

d%sum_rms_mean = d%sum_rms_mean + &
                 sqrt(sum((ens_mean - truth)**2) / d%model_size)

do j = 1, d%ens_size
   d%sum_mean_rms = d%sum_mean_rms + &
                    sqrt(sum((x(:, j) - truth)**2) / d%model_size)
end do

do i = 1, d%num_var
   ind = d%var_list(i)

   ! Do perfect bins first; this is checking the order statistics partitioning

   num = bin_num(x(ind, :), truth(ind))
   d%bin(i, num) = d%bin(i, num) + 1
end do

end subroutine save_diagnostics



  subroutine output_diagnostics(d, chan)
!===========================================================================
! subroutine output_diagnostics(d, chan)
!
! Output diagnostics output from structure d to channel chan.  
! 
! d:		diagnostics type containing accumulated diagnostics
! chan:		integer output channel number
!
! TJH 28 Mar 2002 This routine contained a lot of real/int operations. 
!                 rather than rely on the compiler to promote the integers,
!                 I "ensured" they would get promoted by multiplying by a r8.
!                 Suns compilers _used_to_be_ (?) bad...
!----------------------------------------------------------------------------

implicit none

type (diag_type), intent(in) :: d
integer, intent(in) :: chan

real(r8), dimension(0:d%ens_size) :: expected, num_in_bin
real(r8) :: df, chsq, chsq_prob, desired, trms_mean
real(r8) :: tmean_rms, tbias, tsq_bias
integer  :: i, j
real     :: dummy    ! intentionally single precision, used for printing only

! If no assims were done, just return for now

if(d%count == 0) return

! Correct computation of ensemble mean and mean ensemble RMS

dummy = d%total_rms_mean / d%count*1.0_r8
write(   *, *)'new rms of ensemble mean ', dummy
write(chan, *)'new rms of ensemble mean ', dummy

dummy = d%mean_total_rms / (d%count * d%ens_size)*1.0_r8
write(   *, *)'new ensemble mean of rms ', dummy
write(chan, *)'new ensemble mean of rms ', dummy

! Compute the ratio of the rms of mean / mean rms to the desired ratio

desired = sqrt((d%ens_size + 1) / (2.0_r8 * d%ens_size))
dummy   = (d%total_rms_mean / (d%mean_total_rms / d%ens_size*1.0_r8)) / desired
write(   *, *)'normalized ratio is ', dummy
write(chan, *)'normalized ratio is ', dummy

! Compute global means of these two quantities

write(*, *)
write(   *, *)'rms of ensemble mean is ', real(d%sum_rms_mean / (d%count*1.0_r8))
write(chan, *)'rms of ensemble mean is ', real(d%sum_rms_mean / (d%count*1.0_r8))
write(   *, *)'mean of rms is ', real(d%sum_mean_rms / (d%ens_size * d%count*1.0_r8))
write(chan, *)'mean of rms is ', real(d%sum_mean_rms / (d%ens_size * d%count*1.0_r8))

! Compute the ratio of the rms of mean / mean rms to the desired ratio

desired = sqrt((d%ens_size + 1) / (2.0_r8 * d%ens_size))
dummy   = (d%sum_rms_mean / (d%sum_mean_rms / d%ens_size*1.0_r8)) / desired
write(   *, *)'normalized ratio is ', dummy
write(chan, *)'normalized ratio is ', dummy

write(   *, *) 'mean bias is ', real(d%mean_bias / (d%model_size * d%count*1.0_r8))
write(chan, *) 'mean bias is ', real(d%mean_bias / (d%model_size * d%count*1.0_r8))

! Compute chi-square for truth bins

expected = (d%count * 1.0_r8) / (d%ens_size + 1.0_r8)
do j = 1, d%num_var
   num_in_bin = d%bin(j, :)*1.0_r8
   call chsone(num_in_bin, expected, 1, df, chsq, chsq_prob)
   write(   *, *)'variable ', d%var_list(j), ' chsq ', chsq, chsq_prob
   write(chan, *)'variable ', d%var_list(j), ' chsq ', chsq, chsq_prob
end do

! Open the bin files if this is first call and restart file

if(first) then
   bin_unit   = open_file(file = 'bin',              action = 'write')
   write_chan = open_file(file = 'diag_restart.out', action = 'write')
   first = .false.
end if

! Write out bins 

do i = 0, d%ens_size
   write(       *, 31) i, (d%bin(j, i), j = 1, d%num_var)
   write(    chan, 31) i, (d%bin(j, i), j = 1, d%num_var)
   write(bin_unit, 31) i, (d%bin(j, i), j = 1, d%num_var)
 31 format(10(i5, 1x))
end do

! Write out rms errors and ratios for selected vars

do i = 1, d%num_var
   trms_mean = d%rms_mean(i) / (d%count*1.0_r8)
   tmean_rms = d%mean_rms(i) / (d%count * d%ens_size*1.0_r8)
   tbias     = d%bias(i)     / (d%count*1.0_r8)
   tsq_bias  = d%sq_bias(i)  / (d%count*1.0_r8)
   write(   *, 41) i, trms_mean, tmean_rms, (trms_mean / tmean_rms) / desired, &
      tbias, tsq_bias
   write(chan, 41) i, trms_mean, tmean_rms, (trms_mean / tmean_rms) / desired, &
      tbias, tsq_bias
end do
41 format(1x, 'var ', i2, ' rms ', 3(e10.4, 1x), 'bias ', 2(e10.4, 1x))

! Also, output the structure to a file for possible continuation

write(write_chan, *) d%count, d%model_size, d%ens_size, d%num_var
write(write_chan, *) d%total_rms_mean, d%mean_total_rms, d%mean_bias
write(write_chan, *) d%sum_rms_mean, d%sum_mean_rms, d%bin, d%var_list
write(write_chan, *) d%rms_mean, d%mean_rms, d%bias, d%sq_bias

end subroutine output_diagnostics



  function bin_num(x, a)
!==========================================================================
! function bin_num(x, a)
!
! Determine rank order of a in array x (what bin is a in?)
! 
! x:	Array of arbitrary size
! a:	Value to be rank ordered
!--------------------------------------------------------------------------

implicit none

real(r8), intent(in) :: x(:), a
integer              :: bin_num

integer :: i

bin_num = 0
do i = 1, size(x)
   if(a > x(i)) bin_num = bin_num + 1
end do

end function bin_num



!=========================================================================
! End of assim_diag/assim_diag_mod.f90
!=========================================================================

end module assim_diag_mod

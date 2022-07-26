! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! A variety of PDFs, CDFs, quantile functions and other tools for working with distributions
! to implement quantile conserving filters in observation space and regression in quantile space.

module quantile_distributions_mod

use types_mod, only : r8, digits12, PI
implicit none
private

public :: norm_cdf, norm_inv, weighted_norm_inv, convert_to_probit, convert_from_probit, dist_param_type, &
   convert_all_to_probit, convert_all_from_probit


type dist_param_type
   real(r8), allocatable :: params(:)
end type


contains

!------------------------------------------------------------------------

subroutine convert_all_to_probit(ens_size, num_vars, state_ens, var_kind, p, probit_ens)

integer, intent(in)                  :: ens_size
integer, intent(in)                  :: num_vars
real(r8), intent(in)                 :: state_ens(:, :)
type(dist_param_type), intent(inout) :: p(num_vars)
integer, intent(in)                  :: var_kind(num_vars)
real(r8), intent(out)                :: probit_ens(:, :)

! Note that the input and output arrays may have extra copies (first subscript). Passing sections of a
! leading index could be inefficient for time and storage, so avoiding that for now.

integer  :: i

do i = 1, num_vars
   call convert_to_probit(ens_size, num_vars, state_ens(1:ens_size, :), var_kind(i), p(i), probit_ens(1:ens_size, :))
end do

end subroutine convert_all_to_probit

!------------------------------------------------------------------------

subroutine convert_to_probit(ens_size, num_vars, state_ens, var_kind, p, probit_ens)

integer, intent(in)                  :: ens_size
integer, intent(in)                  :: num_vars
real(r8), intent(in)                 :: state_ens(ens_size)
type(dist_param_type), intent(inout) :: p
integer, intent(in)                  :: var_kind
real(r8), intent(out)                :: probit_ens(ens_size)

! Note that the input and output arrays may have extra copies (first subscript). Passing sections of a
! leading index could be inefficient for time and storage, so avoiding that for now.

real(r8) :: mean, sd

! Initial test is just a bogus thing for normals which require two parameters, mean and sd
mean = sum(state_ens) / ens_size
sd  = sqrt(sum((state_ens - mean)**2) / (ens_size - 1))
! Do the probit transform for the normal
probit_ens = (state_ens - mean) / sd

! Store these for the inversion
allocate(p%params(2))
p%params(1) = mean
p%params(2) = sd

end subroutine convert_to_probit

!------------------------------------------------------------------------

subroutine convert_all_from_probit(ens_size, num_vars, probit_ens, var_kind, p, state_ens)

integer, intent(in)                  :: ens_size
integer, intent(in)                  :: num_vars
real(r8), intent(in)                 :: probit_ens(:, :)
type(dist_param_type), intent(inout) :: p(num_vars)
integer, intent(in)                  :: var_kind(num_vars)
real(r8), intent(out)                :: state_ens(:, :)

! Convert back to the orig
integer  :: i

do i = 1, num_vars
   call convert_from_probit(ens_size, 1, probit_ens(1:ens_size, i), var_kind(i), p(i), state_ens(1:ens_size, i))
end do

end subroutine convert_all_from_probit

!------------------------------------------------------------------------

subroutine convert_from_probit(ens_size, num_vars, probit_ens, var_kind, p, state_ens)

integer, intent(in)                  :: ens_size
integer, intent(in)                  :: num_vars
real(r8), intent(in)                 :: probit_ens(ens_size)
type(dist_param_type), intent(inout) :: p
integer, intent(in)                  :: var_kind
real(r8), intent(out)                :: state_ens(ens_size)

! Convert back to the orig
real(r8) :: mean, sd

mean = p%params(1)
sd   = p%params(2)
state_ens = probit_ens * sd + mean

! Free the storage
deallocate(p%params)

end subroutine convert_from_probit

!------------------------------------------------------------------------

function norm_cdf(x_in, mean, sd)

! Approximate cumulative distribution function for normal
! with mean and sd evaluated at point x_in
! Only works for x>= 0.

real(r8)             :: norm_cdf
real(r8), intent(in) :: x_in, mean, sd

real(digits12) :: x, p, b1, b2, b3, b4, b5, t, density, nx

! Convert to a standard normal
nx = (x_in - mean) / sd

x = abs(nx)


! Use formula from Abramowitz and Stegun to approximate
p = 0.2316419_digits12
b1 = 0.319381530_digits12
b2 = -0.356563782_digits12
b3 = 1.781477937_digits12
b4 = -1.821255978_digits12
b5 = 1.330274429_digits12

t = 1.0_digits12 / (1.0_digits12 + p * x)

density = (1.0_digits12 / sqrt(2.0_digits12 * PI)) * exp(-x*x / 2.0_digits12)

norm_cdf = 1.0_digits12 - density * &
   ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t

if(nx < 0.0_digits12) norm_cdf = 1.0_digits12 - norm_cdf

!write(*, *) 'cdf is ', norm_cdf

end function norm_cdf

!------------------------------------------------------------------------

subroutine weighted_norm_inv(alpha, mean, sd, p, x)

! Find the value of x for which the cdf of a N(mean, sd) multiplied times
! alpha has value p.

real(r8), intent(in)  :: alpha, mean, sd, p
real(r8), intent(out) :: x

real(r8) :: np

! Can search in a standard normal, then multiply by sd at end and add mean
! Divide p by alpha to get the right place for weighted normal
np = p / alpha

! Find spot in standard normal
call norm_inv(np, x)

! Add in the mean and normalize by sd
x = mean + x * sd

end subroutine weighted_norm_inv


!------------------------------------------------------------------------

subroutine norm_inv(p, x)

real(r8), intent(in)  :: p
real(r8), intent(out) :: x

! normal inverse
! translate from http://home.online.no/~pjacklam/notes/invnorm
! a routine written by john herrero

real(r8) :: p_low,p_high
real(r8) :: a1,a2,a3,a4,a5,a6
real(r8) :: b1,b2,b3,b4,b5
real(r8) :: c1,c2,c3,c4,c5,c6
real(r8) :: d1,d2,d3,d4
real(r8) :: q,r
a1 = -39.69683028665376_digits12
a2 =  220.9460984245205_digits12
a3 = -275.9285104469687_digits12
a4 =  138.357751867269_digits12
a5 = -30.66479806614716_digits12
a6 =  2.506628277459239_digits12
b1 = -54.4760987982241_digits12
b2 =  161.5858368580409_digits12
b3 = -155.6989798598866_digits12
b4 =  66.80131188771972_digits12
b5 = -13.28068155288572_digits12
c1 = -0.007784894002430293_digits12
c2 = -0.3223964580411365_digits12
c3 = -2.400758277161838_digits12
c4 = -2.549732539343734_digits12
c5 =  4.374664141464968_digits12
c6 =  2.938163982698783_digits12
d1 =  0.007784695709041462_digits12
d2 =  0.3224671290700398_digits12
d3 =  2.445134137142996_digits12
d4 =  3.754408661907416_digits12
p_low  = 0.02425_digits12
p_high = 1_digits12 - p_low
! Split into an inner and two outer regions which have separate fits
if(p < p_low) then
   q = sqrt(-2.0_digits12 * log(p))
   x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
      ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_digits12)
else if(p > p_high) then
   q = sqrt(-2.0_digits12 * log(1.0_digits12 - p))
   x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
      ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_digits12)
else
   q = p - 0.5_digits12
   r = q*q
   x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q / &
      (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0_digits12)
endif

end subroutine norm_inv

!------------------------------------------------------------------------


end module quantile_distributions_mod

module cov_cutoff_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

use types_mod

contains

!======================================================================



function comp_cov_factor(z_in, c)
!----------------------------------------------------------------------
! function comp_cov_factor(z_in, c)
!
! Computes a covariance cutoff function from Gaspari and Cohn
! QJRMS, 125, 723-757.  (their eqn. 4.10)
!
! z_in is the distance while c is the cutoff distance. 
! For distances greater than 2c, the cov_factor returned goes to 0.

implicit none

real(r8), intent(in) :: z_in, c
real(r8)             :: comp_cov_factor

real(r8) :: z, r

z = abs(z_in)
r = z / c

if( z >= c*2.0_r8 ) then

   comp_cov_factor = 0.0_r8

else if( z >= c .and. z < c*2.0_r8 ) then

   comp_cov_factor = r**5 / 12.0_r8  -  &
                     r**4 / 2.0_r8   +  &
                     r**3 * 5.0_r8 / 8.0_r8 + &
                     r**2 * 5.0_r8 / 3.0_r8 - &
                     5.0_r8*r + 4.0_r8 - (c * 2.0_r8) / (3.0_r8 * z) 
else
   comp_cov_factor = r**5 * (-1.0_r8 / 4.0_r8 ) + &
                     r**4 / 2.0_r8 +              &
                     r**3 * 5.0_r8/8.0_r8 -       &
                     r**2 * 5.0_r8/3.0_r8 + 1.0_r8
endif

end function comp_cov_factor

end module cov_cutoff_mod

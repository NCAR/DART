module cov_cutoff_mod

! Computes a covariance cutoff function from Gaspari and Cohn (their eqn. 4.10)
! QJRMS, 125, 723-757.

! z_in is the distance while c is the cutoff distance. For distances greater
! than 2c, the cov_factor returned goes to 0.

contains

function comp_cov_factor(z_in, c)

implicit none

double precision :: comp_cov_factor
double precision, intent(in) :: z_in, c
double precision :: z, r

z = dabs(z_in)
r = z / c

if(z >= 2*c) then
   comp_cov_factor = 0.0
else if(z >= c .and. z < 2*c) then
   comp_cov_factor = r**5 / 12. - r**4 / 2. + r**3 * 5./8. + r**2 * 5./3. - &
      5*r + 4. - (2 * c) / (3 * z) 
else
   comp_cov_factor = r**5 * (-1./4.) + r**4 / 2. + r**3 * 5./8. - r**2 * 5./3. + 1.
endif

end function comp_cov_factor

end module cov_cutoff_mod

program l96_obs_random

! Needs to be modified to use new random number packages.

use nag_wrap_mod, only : g05caf_wrap

! Places a given number of observations randomly uniformly on the sphere
! Used to test barotropic model ability to deal with increasingly sparse
! observations. Also tacks on the observational variance.

implicit none

integer :: num_obs
real, parameter :: variance = 4.0
real :: x, y, z, lon, length, lat
integer :: i

write(*, *) 'input the number of observations'
read(*, *) num_obs

write(*, *) num_obs
do i = 1, num_obs
! Compute a random point in a volume and then compute direction to surface
   x = real(g05caf_wrap(dble(0.0)))
   write(*, *) x, variance
end do

end program l96_obs_random

program barot_obs_random

use nag_wrap_mod, only : g05ddf_wrap

! Currently uses NAG, must be modified to use available random sequence
! generators.

! Places a given number of observations randomly uniformly on the sphere
! Used to test barotropic model ability to deal with increasingly sparse
! observations. Also tacks on the observational variance.

implicit none

integer :: num_obs
real, parameter :: variance = (1e6)**2
real :: x, y, z, lon, length, lat
integer :: i

write(*, *) 'input the number of observations'
read(*, *) num_obs

write(*, *) num_obs
do i = 1, num_obs
! Compute a random point in a volume and then compute direction to surface
 11   x = real(g05ddf_wrap(dble(0.0), dble(1.0)))
   y = real(g05ddf_wrap(dble(0.0), dble(1.0)))
   z = real(g05ddf_wrap(dble(0.0), dble(1.0)))
! Begin by computing longitude in degrees
   lon = atan2(y, x) * 360.0 / (2.0 * 3.14159) + 180.0
   if(lon < 0.0 .or. lon > 360.0) then
      write(*, *) 'longitude out of bound'
      stop
   endif
   length = sqrt(x**2 + y**2)
   lat = atan2(z, length) * 360.0 / (2.0 * 3.14159)
   if(lat < -90.0 .or. lat > 90.0) then
      write(*, *) 'latitude out of bound'
      stop
   endif
! For barotropic model, don't want to go too close to poles for now
   if(lat < -83 .or. lat > 83) goto 11

!----------------------------------------------------------
! Option for doing NH only
!   if(lat < 0 .or. lat > 85) goto 11
!----------------------------------------------------------
! Option for doing extended WH only
!    if(lon > 60 .and. lon < 160) goto 11
!----------------------------------------------------------
! Option for doing NA only
     if(lon < 180 .or. lon > 300 .or. lat > 80 .or. lat < 20) goto 11
!----------------------------------------------------------
   write(*, *) lon, lat, variance
end do

end program barot_obs_random

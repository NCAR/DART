! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
 
program barot_obs_random

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use    types_mod, only : r8
use nag_wrap_mod, only : g05ddf_wrap

! Currently uses NAG, must be modified to use available random sequence
! generators.

! Places a given number of observations randomly uniformly on the sphere
! Used to test barotropic model ability to deal with increasingly sparse
! observations. Also tacks on the observational variance.

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer :: num_obs
real(r8), parameter :: variance = (1e6)**2_r8
real(r8) :: x, y, z, lon, length, lat
integer :: i

write(*, *) 'input the number of observations'
read(*, *) num_obs

write(*, *) num_obs
do i = 1, num_obs
! Compute a random point in a volume and then compute direction to surface
 11   x = g05ddf_wrap(0.0_r8, 1.0_r8)
   y = g05ddf_wrap(0.0_r8, 1.0_r8)
   z = g05ddf_wrap(0.0_r8, 1.0_r8)
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

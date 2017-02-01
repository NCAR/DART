! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
program barot_obs_random

use    types_mod, only : r8
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

! this routine used to call a special NAG subroutine but we have a
! replacement one in the system now.

! Places a given number of observations randomly uniformly on the sphere
! Used to test barotropic model ability to deal with increasingly sparse
! observations. Also tacks on the observational variance.

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer :: num_obs
real(r8), parameter :: variance = (1e6)**2_r8
real(r8) :: x, y, z, lon, length, lat
integer :: i
type(random_seq_type) :: s

call init_random_sequence(s)

write(*, *) 'input the number of observations'
read(*, *) num_obs

write(*, *) num_obs
do i = 1, num_obs
! Compute a random point in a volume and then compute direction to surface
11 x = random_gaussian(0.0_r8, 1.0_r8)
   y = random_gaussian(0.0_r8, 1.0_r8)
   z = random_gaussian(0.0_r8, 1.0_r8)
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

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

module loc_and_dist_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! Two d spherical location and distance module. Inputs are latitude and 
! longitude in degrees.

private

public loc_type, get_dist, set_loc, get_loc

type loc_type
   private
   double precision :: lon, lat
end type loc_type

contains

!----------------------------------------------------------------------------

subroutine set_loc(loc, lon, lat)

! Given a location type and a double precision value between 0 and 1
! puts this value into the location.

implicit none

type (loc_type), intent(out) :: loc
double precision, intent(in) :: lon, lat

if(lon < 0.0 .or. lon > 360.0 .or. lat < -90.0 .or. lat > 90.0) then
   write(*, *) 'Error in set_loc: value of lon or lat is out range'
   write(*, *) lon, lat
   stop
endif
loc%lon = lon
loc%lat = lat

end subroutine set_loc

!---------------------------------------------------------------------------

subroutine get_loc(loc, lon, lat)

! Given a location type, return lon and lat

implicit none

type(loc_type), intent(in) :: loc
double precision, intent(out) :: lon, lat

lon = loc%lon
lat = loc%lat

end subroutine get_loc

!----------------------------------------------------------------------------

function get_dist(a, b)

! computes pseudo-distance between two locations; just relative distance
! matters so can compute Euclidean distance in 3d rather than on surface
! of sphere

implicit none

type (loc_type), intent(in) :: a, b
double precision :: get_dist

double precision :: xa, xb, ya, yb, za, zb, conv_to_rad


!write(*, *) 'a lon lat ', real(a%lon), real(a%lat)
!write(*, *) 'b lon lat ', real(b%lon), real(b%lat)
conv_to_rad = 2 * 3.14159 / 360.0

xa = cos(a%lon * conv_to_rad) * cos(a%lat * conv_to_rad)
xb = cos(b%lon * conv_to_rad) * cos(b%lat * conv_to_rad)
ya = sin(a%lon * conv_to_rad) * cos(a%lat * conv_to_rad)
yb = sin(b%lon * conv_to_rad) * cos(b%lat * conv_to_rad)
za = sin(a%lat * conv_to_rad)
zb = sin(b%lat * conv_to_rad)

get_dist = (xa - xb)**2 + (ya - yb)**2 + (za-zb)**2
end function get_dist

!----------------------------------------------------------------------------

end module loc_and_dist_mod

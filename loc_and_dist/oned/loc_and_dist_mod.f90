module loc_and_dist_mod

! CVS ID $ID$ 
! CVS Name $Name$ 
! CVS Revision $Revision$ 
! A simple prototype for a state variable location module. For one-d
! model cases this is trivial. Assumes that domain runs from 0.0 to 1.0
! to avoid the need for a domain setting initialization call. This is a
! separate module that can be used by both the model and obs modules since
! distance definitions may be shared by different models (barotropic and
! pe on sphere for instance).

private

public loc_type, get_dist, set_loc, get_loc

type loc_type
   private
   double precision :: x
end type loc_type

contains

!----------------------------------------------------------------------------

subroutine set_loc(loc, x)

! Given a location type and a double precision value between 0 and 1
! puts this value into the location.

implicit none

type (loc_type), intent(out) :: loc
double precision, intent(in) :: x

if(x < 0.0 .or. x > 1.0) then
   write(*, *) 'Error in set_loc: value of x is out of 0->1 range'
   stop
endif
loc%x = x

end subroutine set_loc

!---------------------------------------------------------------------------

subroutine get_loc(loc, x)

! Given a location type, return the x coordinate

implicit none

type(loc_type), intent(in) :: loc
double precision, intent(out) :: x

x = loc%x

end subroutine get_loc

!----------------------------------------------------------------------------

function get_dist(loc1, loc2)

implicit none

type(loc_type), intent(in) :: loc1, loc2
double precision :: get_dist

! Reentrant domain, if distance is greater than half wraparound the other way.
get_dist = dabs(loc1%x - loc2%x)
if(get_dist > 0.5) get_dist = 1.0 - get_dist

end function get_dist

!----------------------------------------------------------------------------

end module loc_and_dist_mod

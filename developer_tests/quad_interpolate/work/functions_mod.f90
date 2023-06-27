module functions_mod

use types_mod, only : r8, DEG2RAD

implicit none

private
public :: sine

contains

pure function sine(x,y)

real(r8) :: sine
real(r8), intent(in) :: x,y


sine = sin(sqrt((x*DEG2RAD)**2 + (y*DEG2RAD)**2))

end function

end module functions_mod

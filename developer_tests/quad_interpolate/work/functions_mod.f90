! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module functions_mod

use types_mod, only : r8, DEG2RAD

implicit none

private
public :: f_sine, f_row, f_col, f_sum, f_rand

contains

!--------------------------------------------------------------------
pure function f_sine(x,y)

real(r8) :: f_sine
real(r8), intent(in) :: x,y

f_sine = sin(sqrt((x*DEG2RAD)**2 + (y*DEG2RAD)**2))

end function

!--------------------------------------------------------------------
pure function f_row(x,y)

real(r8) :: f_row
real(r8), intent(in) :: x,y

f_row = x

end function

!--------------------------------------------------------------------
pure function f_col(x,y)

real(r8) :: f_col
real(r8), intent(in) :: x,y

f_col = y

end function

!--------------------------------------------------------------------
pure function f_sum(x,y)

real(r8) :: f_sum
real(r8), intent(in) :: x,y

f_sum = x + y

end function

!--------------------------------------------------------------------
pure function f_rand(x,y)

real(r8) :: f_rand
real(r8), intent(in) :: x,y


f_rand = mod(x*y, 100.0)

end function

!--------------------------------------------------------------------

end module functions_mod

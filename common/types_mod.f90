
MODULE types_mod

IMPLICIT NONE
SAVE

!----------------------------------------------------------------------------
! Attributes for variable kinds -- no need to rely on -r8 switch in compiler
! all real variables are 64bit, ditto for all integer variables.  
!----------------------------------------------------------------------------
! These are TJH's favorites
! integer, parameter :: i4 = SELECTED_INT_KIND(8)
! integer, parameter :: i8 = SELECTED_INT_KIND(17)
! integer, parameter :: r4 = SELECTED_REAL_KIND(6,30)
! integer, parameter :: c4 = SELECTED_REAL_KIND(6,30)
! integer, parameter :: r8 = SELECTED_REAL_KIND(12,100)
! integer, parameter :: c8 = SELECTED_REAL_KIND(12,100)

! These comply with the CCM4 standard, as far as I can tell.

  integer, parameter :: i4 = SELECTED_INT_KIND(8)
  integer, parameter :: i8 = SELECTED_INT_KIND(13)
  integer, parameter :: r4 = SELECTED_REAL_KIND(6,30)
  integer, parameter :: c4 = SELECTED_REAL_KIND(6,30)
  integer, parameter :: r8 = SELECTED_REAL_KIND(12)
  integer, parameter :: c8 = SELECTED_REAL_KIND(12)

!----------------------------------------------------------------------------
! Constants ... 
!----------------------------------------------------------------------------

real(kind=r8), parameter :: conv_to_rad = 2.0_r8 * 3.14159_r8 / 360.0_r8


END MODULE types_mod

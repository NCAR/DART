! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

MODULE types_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 

implicit none
private 

public :: r8, PI, DEG2RAD, RAD2DEG, MISSING_R4, MISSING_R8
public :: MISSING_I, MISSING_DATA

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

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

real(kind=r8), parameter :: PI = 3.1415926535897932346_r8
real(kind=r8), parameter :: DEG2RAD = PI / 180.0_r8
real(kind=r8), parameter :: RAD2DEG = 180.0_r8 / PI

integer,       PARAMETER ::  MISSING_I    = -888888
integer,       PARAMETER ::  MISSING_DATA = -88
real(kind=r4), PARAMETER ::  MISSING_R4   = -888888.0_r4
real(kind=r8), PARAMETER ::  MISSING_R8   = -888888.0_r8

END MODULE types_mod

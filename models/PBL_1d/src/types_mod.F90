! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

MODULE types_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

implicit none
private 

public :: i4, i8, r4, c4, r8, c8, digits12
public :: PI, DEG2RAD, RAD2DEG, MISSING_R4, MISSING_R8
public :: MISSING_I, MISSING_DATA
public :: t_kelvin, es_alpha, es_beta, es_gamma, gas_constant_v, gas_constant
public :: L_over_Rv, ps0, earth_radius, gravity

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
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
! comment in only one of the following lines, not both.   the first is the default.
  integer, parameter :: r8 = SELECTED_REAL_KIND(12)   ! real r8
!!integer, parameter :: r8 = r4                       ! alias r8 to r4
  integer, parameter :: c8 = SELECTED_REAL_KIND(12)

! 'digits12' is reserved for real variables that MUST retain 64 bits of
! precision. DO NOT CHANGE '12' to a smaller number. BAD BAD BAD things happen.
! This is a small subset of the variables. Changing this will ruin the ability
! to distinguish timesteps that are a few seconds apart, for instance.

  integer, parameter :: digits12 = SELECTED_REAL_KIND(12)

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

real(r8), PARAMETER :: t_kelvin       = 273.15_r8
real(r8), PARAMETER :: es_alpha       = 611.2_r8
real(r8), PARAMETER :: es_beta        = 17.67_r8
real(r8), PARAMETER :: es_gamma       = 243.5_r8
real(r8), PARAMETER :: gas_constant_v = 461.6_r8
real(r8), PARAMETER :: gas_constant   = 287.0_r8
real(r8), PARAMETER :: L_over_Rv      = 5418.12_r8
real(r8), PARAMETER :: ps0            = 100000.0_r8    ! Base sea level pressure
real(r8), PARAMETER :: earth_radius   = 6370.0_r8      ! km, consistant with WRF
real(r8), PARAMETER :: gravity        = 9.81_r8

END MODULE types_mod

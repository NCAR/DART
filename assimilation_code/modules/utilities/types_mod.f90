! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Should only have fortran "kinds" - number of bytes/precision - for
!> integer, real, etc.  the constants should go into their own module, 
!> or this should be renamed 'constants' or something more descriptive 
!> than 'types' (which overlaps with dart's use of kinds and types).

module types_mod

implicit none
private 

public :: i2, i4, i8, r4, c4, r8, c8, digits12
public :: PI, DEG2RAD, RAD2DEG, MISSING_R4, MISSING_R8
public :: MISSING_I, MISSING_I8, MISSING_DATA
public :: SECPERDAY
public :: t_kelvin, es_alpha, es_beta, es_gamma, gas_constant_v, gas_constant
public :: L_over_Rv, ps0, earth_radius, gravity
public :: metadatalength, obstypelength, varnamelength, vtablenamelength
public :: MAX_NUM_DOMS, MAX_FILES

!----------------------------------------------------------------------------
! constants that need to be shared - really has nothing to do with types ...
!----------------------------------------------------------------------------

integer, parameter :: metadatalength   = 64   ! 75 is max w/out wrapping in ncdump
integer, parameter :: varnamelength    = 31   ! max name of any fortran variable
                                              ! in F90/95. increased to 63 in F2003
integer, parameter :: obstypelength    = 31   ! because of variable name limits.
                                              ! use this for obs types/state quantities
integer, parameter :: vtablenamelength = 64   ! use this for the model_mod variable table items
                                           
integer, parameter :: MAX_NUM_DOMS     = 10   ! max num domains. this is arbitrarily 
                                              ! based on WRF's maximum number of domains

integer, parameter :: MAX_FILES        = 1000 ! maximum number of files

!----------------------------------------------------------------------------
! Attributes for variable kinds -- declaring sizes explicitly means we know
! exactly what precision we are using and are not relying on compiler flags
! to set the defaults for real and int.
!
! the call is: selected_int_kind(p, r)
!  where p = decimal_precision,  and r = decimal_exponent_range
!
! below, "external32_size" = number of bytes 
! 
!
! from the MPI documentation, they say:
! (this block of text is originally from the MPI documentation but it seems to apply
! to the compiler intrinsic functions SELECTED_REAL_KIND and SELECTED_INT_KIND)
!
! for reals, SELECTED_REAL_KIND(p, r) values result in:
!
!   if      (p > 33) or (r > 4931) then  external32 representation 
!                                        is undefined   
!   else if (p > 15) or (r >  307) then  external32_size = 16 
!   else if (p >  6) or (r >   37) then  external32_size =  8 
!   else                                 external32_size =  4 
!
!
! for integers, SELECTED_INT_KIND(r) results in:
!
!   if      (r > 38) then  external32 representation is undefined 
!   else if (r > 18) then  external32_size =  16  
!   else if (r >  9) then  external32_size =  8  
!   else if (r >  4) then  external32_size =  4 
!   else if (r >  2) then  external32_size =  2  
!   else                   external32_size =  1  
!
!----------------------------------------------------------------------------

! integer precision:
integer, parameter :: i2 = SELECTED_INT_KIND(3)
integer, parameter :: i4 = SELECTED_INT_KIND(8)
integer, parameter :: i8 = SELECTED_INT_KIND(13)

! real precision:
! TO RUN WITH REDUCED PRECISION REALS (and use correspondingly less memory)
! comment OUT the r8 definition below and use the second one:
integer, parameter :: r4 = SELECTED_REAL_KIND(6,30)
integer, parameter :: r8 = SELECTED_REAL_KIND(12)   ! 8 byte reals
!integer, parameter :: r8 = r4                      ! alias r8 to r4

! complex precision:
integer, parameter :: c4 = SELECTED_REAL_KIND(6,30)
integer, parameter :: c8 = SELECTED_REAL_KIND(12)



! 'digits12' is reserved for real variables that MUST retain 64 bits of
! precision. DO NOT CHANGE '12' to a smaller number. BAD BAD BAD things happen.
! This is a small subset of the variables. Changing this will ruin the ability
! to distinguish timesteps that are a few seconds apart, for instance.

integer, parameter :: digits12 = SELECTED_REAL_KIND(12)

!----------------------------------------------------------------------------
! Constants, some of which are model dependent and shouldn't be in this
! module.
!----------------------------------------------------------------------------

integer,       parameter ::  SECPERDAY = 24*60*60

real(kind=r8), parameter :: PI = 3.14159265358979323846_r8
real(kind=r8), parameter :: DEG2RAD = PI / 180.0_r8
real(kind=r8), parameter :: RAD2DEG = 180.0_r8 / PI

integer,          parameter ::  MISSING_I    = -888888
integer(kind=i8), parameter ::  MISSING_I8   = -888888_I8
integer,          parameter ::  MISSING_DATA = -88
real(kind=r4),    parameter ::  MISSING_R4   = -888888.0_r4
real(kind=r8),    parameter ::  MISSING_R8   = -888888.0_r8

real(r8), parameter :: t_kelvin       = 273.15_r8
real(r8), parameter :: es_alpha       = 611.2_r8
real(r8), parameter :: es_beta        = 17.67_r8
real(r8), parameter :: es_gamma       = 243.5_r8
real(r8), parameter :: gas_constant_v = 461.6_r8
real(r8), parameter :: gas_constant   = 287.0_r8 ! wikipedia has 287.06, WRF has 287.05 ...
real(r8), parameter :: L_over_Rv      = 5418.12_r8
real(r8), parameter :: ps0            = 100000.0_r8    ! Base sea level pressure
real(r8), parameter :: earth_radius   = 6370.0_r8      ! km, consistant with WRF
real(r8), parameter :: gravity        = 9.81_r8        ! wikipedia has 9.80665

end module types_mod


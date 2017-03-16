module constants
!$$$   module documentation block
!                .      .    .                                       .
! module:    constants
!   prgmmr: treadon          org: np23                date: 2003-09-25
!
! abstract:  This module contains the definition of various constants
!            used in the gsi code
!
! program history log:
!   2003-09-25  treadon - original code
!   2004-03-02  treadon - allow global and regional constants to differ
!   2004-06-16  treadon - update documentation
!   2004-10-28  treadon - replace parameter tiny=1.e-12 with tiny_r_kind
!                         and tiny_single
!   2004-11-16  treadon - add huge_single, huge_r_kind parameters
!   2005-08-24  derber   - move cg_term to constants from qcmod
!   2006-03-07  treadon  - add rd_over_cp_mass
!   2006-05-18  treadon  - add huge_i_kind
!   2006-06-06       su  - add var-qc wgtlim, change value to 0.25 (ECMWF)
!   2006-07-28  derber   - add r1000
!   2007-03-20  rancic   - add r3600
!   2009-02-05  cucurull - modify refractive indexes for gpsro data
!   2010-08-25  cucurull - add constants to compute compressibility factor
!                        - add option to use Rueger/Bevis refractive index coeffs
!   2010-12-20 pagowski  - add max_varname_length=12
!   2010-04-01 li        - add maximum diurnal thermocline thickness
!   2011-10-27 Huang     - add i_missing and r_missing to detect missing values
!   2011-11-01 eliu      - add minimum value for cloud water mixing ratio 
!   2012-03-07 todling   - define lower bound for trace-gases (arbitrary unit as long as small)
!
! Subroutines Included:
!   sub init_constants_derived - compute derived constants
!   sub init_constants         - set regional/global constants
!   sub gps_constants          - set Rueger/Bevis refractive index coefficients
!
! Variable Definitions:
!   see below
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$ end documentation block

  use kinds, only: r_single,r_kind,i_kind,r_quad,i_long
  implicit none

! set default as private
  private
! set subroutines as public
  public :: init_constants_derived
  public :: init_constants
  public :: gps_constants
! set passed variables to public
  public :: one,two,half,zero,deg2rad,pi,three,quarter,one_tenth
  public :: rad2deg,zero_quad,r3600,r1000,r60inv,five,four,rd_over_cp,grav
  public :: rd,rv,rozcon,rearth_equator,zero_single,tiny_r_kind,tiny_single,ten
  public :: omega,rcp,rearth,fv,h300,cp,cg_term,tpwcon,xb,ttp,psatk,xa,tmix
  public :: xai,xbi,psat,eps,omeps,wgtlim,one_quad,two_quad,epsq,climit,epsm1,hvap
  public :: hsub,cclimit,el2orc,elocp,h1000,cpr,pcpeff0,pcpeff2,delta,pcpeff1
  public :: factor1,c0,pcpeff3,factor2,dx_inv,dx_min,rhcbot,rhctop,hfus,ke2
  public :: rrow,cmr,cws,r60,huge_i_kind,huge_r_kind,t0c,rd_over_cp_mass
  public :: somigliana,grav_equator,grav_ratio,flattening,semi_major_axis
  public :: n_b,n_a,eccentricity,huge_single,constoz,g_over_rd,amsua_clw_d2
  public :: amsua_clw_d1,n_c,rd_over_g,zero_ilong
  public :: r10,r100,sqrt_tiny_r_kind,r2000,r4000
  public :: r0_01,r0_02,r0_03,r0_04,r0_05,r400,r2400
  public :: cpf_a0, cpf_a1, cpf_a2, cpf_b0, cpf_b1, cpf_c0, cpf_c1, cpf_d, cpf_e
  public :: psv_a, psv_b, psv_c, psv_d
  public :: ef_alpha, ef_beta, ef_gamma
  public :: max_varname_length
  public :: z_w_max,tfrozen
  public :: qmin,qcmin,tgmin
  public :: i_missing, r_missing

! Declare derived constants
  integer(i_kind):: huge_i_kind
  integer(i_kind), parameter :: max_varname_length=12
  real(r_single):: tiny_single, huge_single
  real(r_kind):: xai, xa, xbi, xb, dldt, rozcon,ozcon,fv, tpwcon,eps, rd_over_g
  real(r_kind):: el2orc, g_over_rd, rd_over_cp, cpr, omeps, epsm1, factor2
  real(r_kind):: factor1, huge_r_kind, tiny_r_kind, deg2rad, pi, rad2deg, cg_term
  real(r_kind):: eccentricity_linear, cv, rv, rd_over_cp_mass, cliq, rd, cp_mass
  real(r_kind):: eccentricity, grav, rearth, r60inv
  real(r_kind):: sqrt_tiny_r_kind
  real(r_kind):: n_a, n_b, n_c

! Define constants common to global and regional applications
  real(r_kind),parameter::  rearth_equator= 6.37813662e6_r_kind  ! equatorial earth radius          (m)
  real(r_kind),parameter::  omega  = 7.2921e-5_r_kind            !  angular velocity of earth       (1/s)
  real(r_kind),parameter::  cp     = 1.0046e+3_r_kind            !  specific heat of air @pressure  (J/kg/K)
  real(r_kind),parameter::  cvap   = 1.8460e+3_r_kind            !  specific heat of h2o vapor      (J/kg/K)
  real(r_kind),parameter::  csol   = 2.1060e+3_r_kind            !  specific heat of solid h2o (ice)(J/kg/K)
  real(r_kind),parameter::  hvap   = 2.5000e+6_r_kind            !  latent heat of h2o condensation (J/kg)
  real(r_kind),parameter::  hfus   = 3.3358e+5_r_kind            !  latent heat of h2o fusion       (J/kg)
  real(r_kind),parameter::  psat   = 6.1078e+2_r_kind            !  pressure at h2o triple point    (Pa)
  real(r_kind),parameter::  t0c    = 2.7315e+2_r_kind            !  temperature at zero celsius     (K)
  real(r_kind),parameter::  ttp    = 2.7316e+2_r_kind            !  temperature at h2o triple point (K)
  real(r_kind),parameter::  jcal   = 4.1855e+0_r_kind            !  joules per calorie              ()
  real(r_kind),parameter::  stndrd_atmos_ps = 1013.25e2_r_kind   ! 1976 US standard atmosphere ps   (Pa)

! Numeric constants

  integer(i_long),parameter::  zero_ilong = 0_i_long

  real(r_single),parameter::  zero_single= 0.0_r_single

  real(r_kind),parameter::  zero      = 0.0_r_kind
  real(r_kind),parameter::  r0_01     = 0.01_r_kind
  real(r_kind),parameter::  r0_02     = 0.02_r_kind
  real(r_kind),parameter::  r0_03     = 0.03_r_kind
  real(r_kind),parameter::  r0_04     = 0.04_r_kind
  real(r_kind),parameter::  r0_05     = 0.05_r_kind
  real(r_kind),parameter::  one_tenth = 0.10_r_kind
  real(r_kind),parameter::  quarter   = 0.25_r_kind
  real(r_kind),parameter::  one       = 1.0_r_kind
  real(r_kind),parameter::  two       = 2.0_r_kind
  real(r_kind),parameter::  three     = 3.0_r_kind
  real(r_kind),parameter::  four      = 4.0_r_kind
  real(r_kind),parameter::  five      = 5.0_r_kind
  real(r_kind),parameter::  ten       = 10.0_r_kind
  real(r_kind),parameter::  r10       = 10.0_r_kind
  real(r_kind),parameter::  r60       = 60._r_kind
  real(r_kind),parameter::  r100      = 100.0_r_kind
  real(r_kind),parameter::  r400      = 400.0_r_kind
  real(r_kind),parameter::  r1000     = 1000.0_r_kind
  real(r_kind),parameter::  r2000     = 2000.0_r_kind
  real(r_kind),parameter::  r2400     = 2400.0_r_kind
  real(r_kind),parameter::  r4000     = 4000.0_r_kind
  real(r_kind),parameter::  r3600     = 3600.0_r_kind

  real(r_kind),parameter:: z_w_max    = 30.0_r_kind     ! maximum diurnal thermocline thickness
  real(r_kind),parameter:: tfrozen    = 271.2_r_kind    ! sea water frozen point temperature

  real(r_quad),parameter::  zero_quad = 0.0_r_quad
  real(r_quad),parameter::  one_quad  = 1.0_r_quad
  real(r_quad),parameter::  two_quad  = 2.0_r_quad

! Constants for compressibility factor (Davis et al 1992)
  real(r_kind),parameter::  cpf_a0 =  1.58123e-6_r_kind ! K/Pa
  real(r_kind),parameter::  cpf_a1 = -2.9331e-8_r_kind  ! 1/Pa
  real(r_kind),parameter::  cpf_a2 =  1.1043e-10_r_kind ! 1/K 1/Pa
  real(r_kind),parameter::  cpf_b0 =  5.707e-6_r_kind   ! K/Pa
  real(r_kind),parameter::  cpf_b1 = -2.051e-8_r_kind   ! 1/Pa
  real(r_kind),parameter::  cpf_c0 =  1.9898e-4_r_kind  ! K/Pa
  real(r_kind),parameter::  cpf_c1 = -2.376e-6_r_kind   ! 1/Pa
  real(r_kind),parameter::  cpf_d  =  1.83e-11_r_kind   ! K2/Pa2
  real(r_kind),parameter::  cpf_e  = -0.765e-8_r_kind   ! K2/Pa2

! Constants for vapor pressure at saturation
  real(r_kind),parameter::  psv_a =  1.2378847e-5_r_kind       !  (1/K2)
  real(r_kind),parameter::  psv_b = -1.9121316e-2_r_kind       !  (1/K)
  real(r_kind),parameter::  psv_c = 33.93711047_r_kind         !
  real(r_kind),parameter::  psv_d = -6.3431645e+3_r_kind       !  (K)

! Constants for enhancement factor to calculating the mole fraction of water vapor
  real(r_kind),parameter::  ef_alpha = 1.00062_r_kind           !
  real(r_kind),parameter::  ef_beta  = 3.14e-8_r_kind           !  (1/Pa)
  real(r_kind),parameter::  ef_gamma = 5.6e-7_r_kind            !  (1/K2)

! Parameters below from WGS-84 model software inside GPS receivers.
  real(r_kind),parameter::  semi_major_axis = 6378.1370e3_r_kind     !                     (m)
  real(r_kind),parameter::  semi_minor_axis = 6356.7523142e3_r_kind  !                     (m)
  real(r_kind),parameter::  grav_polar      = 9.8321849378_r_kind    !                     (m/s2)
  real(r_kind),parameter::  grav_equator    = 9.7803253359_r_kind    !                     (m/s2) 
  real(r_kind),parameter::  earth_omega     = 7.292115e-5_r_kind     !                     (rad/s)
  real(r_kind),parameter::  grav_constant   = 3.986004418e14_r_kind  !                     (m3/s2)

! Derived geophysical constants
  real(r_kind),parameter::  flattening = (semi_major_axis-semi_minor_axis)/semi_major_axis
  real(r_kind),parameter::  somigliana = &
       (semi_minor_axis/semi_major_axis) * (grav_polar/grav_equator) - one
  real(r_kind),parameter::  grav_ratio = (earth_omega*earth_omega * &
       semi_major_axis*semi_major_axis * semi_minor_axis) / grav_constant 

! Derived thermodynamic constants
  real(r_kind),parameter::  dldti = cvap-csol
  real(r_kind),parameter::  hsub = hvap+hfus
  real(r_kind),parameter::  psatk = psat*0.001_r_kind
  real(r_kind),parameter::  tmix = ttp-20._r_kind
  real(r_kind),parameter::  elocp = hvap/cp
  real(r_kind),parameter::  rcp  = one/cp

! Constants used in GFS moist physics
  real(r_kind),parameter::  h300 = 300._r_kind
  real(r_kind),parameter::  half = 0.5_r_kind
  real(r_kind),parameter::  cclimit = 0.001_r_kind
  real(r_kind),parameter::  climit = 1.e-20_r_kind
  real(r_kind),parameter::  epsq = 2.e-12_r_kind
  real(r_kind),parameter::  h1000 = r1000
  real(r_kind),parameter::  rhcbot=0.85_r_kind
  real(r_kind),parameter::  rhctop=0.85_r_kind
  real(r_kind),parameter::  dx_max=-8.8818363_r_kind
  real(r_kind),parameter::  dx_min=-5.2574954_r_kind
  real(r_kind),parameter::  dx_inv=one/(dx_max-dx_min)
  real(r_kind),parameter::  c0=0.002_r_kind
  real(r_kind),parameter::  delta=0.6077338_r_kind
  real(r_kind),parameter::  pcpeff0=1.591_r_kind
  real(r_kind),parameter::  pcpeff1=-0.639_r_kind
  real(r_kind),parameter::  pcpeff2=0.0953_r_kind
  real(r_kind),parameter::  pcpeff3=-0.00496_r_kind
  real(r_kind),parameter::  cmr = one/0.0003_r_kind
  real(r_kind),parameter::  cws = 0.025_r_kind
  real(r_kind),parameter::  ke2 = 0.00002_r_kind
  real(r_kind),parameter::  row = r1000
  real(r_kind),parameter::  rrow = one/row
! real(r_kind),parameter::  qmin = 1.e-7_r_kind  !lower bound on ges_q

! Constant used to process ozone
  real(r_kind),parameter::  constoz = 604229.0_r_kind

! Constants used in cloud liquid water correction for AMSU-A
! brightness temperatures
  real(r_kind),parameter::  amsua_clw_d1 = 0.754_r_kind
  real(r_kind),parameter::  amsua_clw_d2 = -2.265_r_kind

! Constants used for variational qc
  real(r_kind),parameter::  wgtlim = quarter  ! Cutoff weight for concluding that obs has been
                                     ! rejected by nonlinear qc. This limit is arbitrary
                                     ! and DOES NOT affect nonlinear qc. It only affects
                                     ! the printout which "counts" the number of obs that
                                     ! "fail" nonlinear qc.  Observations counted as failing
                                     ! nonlinear qc are still assimilated.  Their weight
                                     ! relative to other observations is reduced. Changing
                                     ! wgtlim does not alter the analysis, only
                                     ! the nonlinear qc data "count"

! Minimum values for water vapor, cloud water mixing ratio, and trace gases
  real(r_kind),parameter:: qmin   = 1.e-07_r_kind   ! lower bound on ges_q
  real(r_kind),parameter:: qcmin  = 0.0_r_kind      ! lower bound on ges_cw
  real(r_kind),parameter:: tgmin  = 1.e-15_r_kind   ! lower bound on trace gases

! Constant used to detect missing input value
  integer(i_kind),parameter:: i_missing=-9999
  integer(r_kind),parameter:: r_missing=-9999._r_kind

contains

  subroutine init_constants_derived
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_constants_derived          set derived constants
!     prgmmr:    treadon          org: np23           date: 2004-12-02
! ! abstract:  This routine sets derived constants
!
! program history log:
!   2004-12-02  treadon
!   2005-03-03  treadon - add implicit none
!   2008-06-04  safford - rm unused vars
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!$$$
    implicit none

!   Trigonometric constants
    pi      = acos(-one)
    deg2rad = pi/180.0_r_kind
    rad2deg = one/deg2rad
    cg_term = (sqrt(two*pi))/two                  ! constant for variational qc
    tiny_r_kind = tiny(zero)
    sqrt_tiny_r_kind = r10*sqrt(tiny_r_kind)
    huge_r_kind = huge(zero)
    tiny_single = tiny(zero_single)
    huge_single = huge(zero_single)
    huge_i_kind = huge(0)
    r60inv=one/r60

!   Geophysical parameters used in conversion of geopotential to
!   geometric height
    eccentricity_linear = sqrt(semi_major_axis**2 - semi_minor_axis**2)
    eccentricity = eccentricity_linear / semi_major_axis

    return
  end subroutine init_constants_derived

  subroutine init_constants(regional)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_constants       set regional or global constants
!     prgmmr:    treadon          org: np23           date: 2004-03-02
!
! abstract:  This routine sets constants specific to regional or global
!            applications of the gsi
!
! program history log:
!   2004-03-02  treadon
!   2004-06-16  treadon, documentation
!   2004-10-28  treadon - use intrinsic TINY function to set value
!                         for smallest machine representable positive
!                         number
!   2004-12-03  treadon - move derived constants to init_constants_derived
!   2005-03-03  treadon - add implicit none
!
!   input argument list:
!     regional - if .true., set regional gsi constants;
!                otherwise (.false.), use global constants
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!$$$
    implicit none

    logical,intent(in   ) :: regional

    real(r_kind) reradius,g,r_d,r_v,cliq_wrf

!   Define regional constants here
    if (regional) then

!      Name given to WRF constants
       reradius = one/6370.e03_r_kind
       g        = 9.81_r_kind
       r_d      = 287.04_r_kind
       r_v      = 461.6_r_kind
       cliq_wrf = 4190.0_r_kind
       cp_mass  = 1004.67_r_kind

!      Transfer WRF constants into unified GSI constants
       rearth = one/reradius
       grav   = g
       rd     = r_d
       rv     = r_v
       cv     = cp-r_d
       cliq   = cliq_wrf
       rd_over_cp_mass = rd / cp_mass

!   Define global constants here
    else
       rearth = 6.3712e+6_r_kind
       grav   = 9.80665e+0_r_kind
       rd     = 2.8705e+2_r_kind
       rv     = 4.6150e+2_r_kind
       cv     = 7.1760e+2_r_kind
       cliq   = 4.1855e+3_r_kind
       cp_mass= zero
       rd_over_cp_mass = zero
    endif


!   Now define derived constants which depend on constants
!   which differ between global and regional applications.

!   Constants related to ozone assimilation
    ozcon = grav*21.4e-9_r_kind
    rozcon= one/ozcon

!   Constant used in vertical integral for precipitable water
    tpwcon = 100.0_r_kind/grav

!   Derived atmospheric constants
    fv         = rv/rd-one    ! used in virtual temperature equation 
    dldt       = cvap-cliq
    xa         = -(dldt/rv)
    xai        = -(dldti/rv)
    xb         = xa+hvap/(rv*ttp)
    xbi        = xai+hsub/(rv*ttp)
    eps        = rd/rv
    epsm1      = rd/rv-one
    omeps      = one-eps
    factor1    = (cvap-cliq)/rv
    factor2    = hvap/rv-factor1*t0c
    cpr        = cp*rd
    el2orc     = hvap*hvap/(rv*cp)
    rd_over_g  = rd/grav
    rd_over_cp = rd/cp
    g_over_rd  = grav/rd

    return
  end subroutine init_constants

  subroutine gps_constants(use_compress)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gps_constants     set Bevis or Rueger refractive index coeff
!     prgmmr:    cucurull          org: np23           date: 2010-08-25
!
! abstract:  This routine sets constants for the refractivity equation. GSI uses Bevis 
!            coefficients when the compressibility factors option is turned off 
!            and uses Rueger coefficients otherwise.
!
! program history log:
!   2010-08-25  cucurull
!   2010-08-25  cucurull, documentation
!
!   input argument list:
!     compress - if .true., set Rueger coefficients;
!                otherwise (.false.), use Bevis coefficients
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!$$$
    implicit none

    logical,intent(in   ) :: use_compress

!   Define refractive index coefficients here
    if (use_compress) then

       ! Constants for gpsro data (Rueger 2002)
       n_a = 77.6890_r_kind   ! K/mb
       n_b = 3.75463e+5_r_kind  ! K^2/mb
       n_c = 71.2952_r_kind   ! K/mb
    else
       ! Constants for gpsro data (Bevis et al 1994)
       n_a = 77.60_r_kind     ! K/mb
       n_b = 3.739e+5_r_kind  ! K^2/mb
       n_c = 70.4_r_kind      ! K/mb
    endif

    return
  end subroutine gps_constants

end module constants

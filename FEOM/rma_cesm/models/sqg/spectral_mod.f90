! This code is not protected by the DART copyright agreement.
! DART $Id$

!============================================================
! a module to be used with spectral surface-based models.
!============================================================

MODULE spectral_mod

    implicit none

    real,    parameter :: pi = 4.0*atan(1.0)    ! pi = 3.14159265358979323846264338327950288419 

    integer, parameter :: model = 1             ! 0:2D; 1:2sQG; 2:sQG-trop; 3:sQG-sfc; 4:HsQG
    real,    parameter :: dt    = 0.0108        ! model time step

    integer, parameter :: ntims = 1000001       ! number of model time steps
    integer, parameter :: iplot = 1000          ! plots every iplot time steps
    integer, parameter :: imean = int(10.e3/dt) ! compute area mean every imean

    logical, parameter :: linear = .FALSE.      ! flag for linear advection
    logical, parameter :: trop   = .FALSE.      ! flag for tropopause geometry
    logical, parameter :: grow   = .FALSE.      ! grow a normal mode (.TRUE./.FALSE.)
    logical, parameter :: inorm  = .FALSE.      ! flag for computing norm

    logical, parameter :: iterr = .FALSE.       ! flag for terrain
    real,    parameter :: hamp  = 0.6           ! height of terrain (0.6)
    real,    parameter :: asx   = 0.6           ! gaussian terrian major axis (1.0)
    real,    parameter :: asy   = 0.6           ! gaussian terrian minor axis (1.0)
    real,    parameter :: asig  = 1.0           ! eccentricity of gaussian terrain (1.0)

    integer, parameter :: kmax = 128/2          ! number of x waves
    integer, parameter :: lmax = 64/2           ! number of y waves

    logical, parameter :: hw    = .TRUE.        ! Hoskins-West jet on/off (.TRUE./.FALSE.)
    real,    parameter :: hwp   = 5.539118      ! Hoskins-West jet width parameter
    real,    parameter :: amu   = 1.0           ! Hoskins-West jet concentration parameter (0->1)
    real,    parameter :: shear = 1.0           ! shear parameter (1 for HW jet)
    real,    parameter :: trl   = 15.0          ! jet relaxation parameter [e-folding time] (20.0, 40.0, 1000.0)

    real,    parameter :: Ross  = 0.0           ! Rossby number (0.0, 0.1)
    real,    parameter :: gamma = 0.0           ! Ekman parameter (0.0, 0.075, 0.075*1.5)

    integer, parameter :: n   = 8               ! diffusion parameter
    real,    parameter :: tau = 20.0*dt         ! diffusion time scale 

    real,    parameter :: XL = 20.0             ! x domain length (2*pi, 20.0, 40.0)
    real,    parameter :: YL = 2*hwp            ! y domain length (2*pi, 2*5.539118)
    real,    parameter :: H  = 1.0              ! z domain length (0.01, 0.1, 1.0, 10.0)

    integer, parameter :: verbose = 2           ! 0:no mesgs; 1:imp only; 2:all

    real,    parameter :: ZH    = H             ! z domain length ( needed for pinv )
    integer, parameter :: pmax  = 11            ! number of vertical levels ( == " == )
    integer, parameter :: order = 2             ! 0:leading order only; 1-first order only; 2-full version ( == " == )

    !================================
    ! do not modify below this point:
    !================================

    integer, parameter :: kmid=kmax/2,     lmid=lmax/2,     &
                          mmax=3.125*kmax, nmax=3.125*lmax, &
                          k2=mmax-kmax,    l2=nmax-lmax,    &
                          kmaxp1=kmax+1,   lmaxp1=lmax+1
    real,    parameter :: facx=2.0*pi/XL,  facy=2.0*pi/YL
    real,    parameter :: ryl=2.0*pi/hwp

    type derivative_operators
        complex, dimension(2*kmax,2*lmax) :: dx,dy,dz,dzo,iz,izo,Id
    end type derivative_operators
    type(derivative_operators) :: d_oper

END MODULE spectral_mod

!============================================================

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

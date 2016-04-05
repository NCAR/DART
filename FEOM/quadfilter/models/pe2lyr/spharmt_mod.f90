! This code is not protected by the DART copyright agreement.
! DART $Id$

 module spharmt_mod

! fortran95 spherical harmonic module.
! For gaussian grids and triangular truncation only.
!
! version 1.1  9/30/2003 (second version - now all fortran using fft99)
! Jeff Whitaker <Jeffrey.S.Whitaker@noaa.gov>

!-----------------------------------------------------------------------------
! to use this module, put "use spharmt" at top of main program
! (just after program statement).
! to initialize arrays needed to compute spherical harmonic transforms
! at T(ntrunc) resolution on a nlon x nlat gaussian grid, 
! do something like this:
! 
! type (sphere) :: sphere_dat
! parameter (nlon=192,nlat=nlon/2)
! parameter (ntrunc=62)
! rsphere = 6.3712e6
! call spharmt_init(sphere_dat,nlon,nlat,ntrunc,re)
 
! derived data type "sphere" contains stuff needed for legendre
! transforms and fft.  By creating multiple instances of sphere
! data types, spherical harmonic transforms on multiple grids can
! be done easily in the same code.
!
! Components of sphere derived data type are:
!
! NLON (integer) - number of longitudes
! NLAT (integer) - number of latitudes
! NTRUNC (integer) - triangular truncation limit
! RE (real) - radius of sphere (m).
! PNM (real pointer array dimension ((ntrunc+1)*(ntrunc+2)/2, nlat) - 
! Associated legendre polynomials.
! HNM (real pointer array, same size as pnm) = -(1 - x*x) * d(pnm)/dx  
! at x = sin(latitude).
! GAULATS (real pointer array dimension nlat) - sin(gaussian latitudes).
! WEIGHTS (real pointer array dimension nlat) - gaussian weights.
! GWRC (real pointer array dimension nlat) - gaussian weights divided by 
! rsphere*(1-x**2).
! INDXM (integer pointer array dimension (ntrunc+1)*(ntrunc+2)/2) - value of
! zonal wavenumber m corresponding to spectral array index 
! nm=1,(ntrunc+1)*(ntrunc+2)/2)
! INDXN (integer pointer array same size as indxm) - value of
! spherical harmonic degree n corresponding to spectral array index 
! nm=1,(ntrunc+1)*(ntrunc+2)/2)
! LAP (real pointer array dimension (ntrunc+1)*(ntrunc+2)/2) - 
! lapacian operator in spectral space = 
! -n*(n+1)/rsphere**2, where n is degree of spherical harmonic.
! ILAP (real pointer array same size as lap) - inverse laplacian operator in 
! spectral space.
! TRIGS, IFAX: arrays needed by Temperton FFT.
! ISINITIALIZED (logical) - true if derived data type has been
! initialized by call to spharm_init, false if not initialized.
! 

! public routines:
!
! SPHARMT_INIT(sphere_dat,nlon,nlat,ntrunc,re):  initialize a sphere object
! (sphere_dat - derived data type "sphere"). inputs are nlon (number of unique
! longitudes), nlat (number of gaussian latitudes), and re 
! (radius of sphere in m). Must be called before anything else.
!
! SPHARMT_DESTROY(sphere_dat):  cleans up pointer arrays allocated by 
! spharmt_init.
!
! SPHARM(sphere_dat, ugrid, anm, idir):  spherical harmonic transform 
! (forward, i.e. grid to spectral, for idir=1 and backward for idir=-1). 
! Arguments are sphere derived data type (sphere_dat), gridded data (ugrid),
! complex spectral coefficients (anm), and flag specifying direction of 
! transform (idir).  See "Import Details" below for information 
! about indexing of grid and spectral arrays.

! GAULW(gaulats,weights): compute sin(gaussian latitudes) and gaussian
! weights.  Number of latitudes determined by size of input arrays.

! LEGEND(x,pnm,hnm): compute associated legendre functions (pnm) and
! their derivates hnm = (X**2-1)*D(PNM)/DX at x = sin(latitude). The
! input arrays pnm and hnm should have dimension (ntrunc+1)*(ntrunc+2)/2, 
! where ntrunc is triangular truncation limit (see Import Details below
! for a description of the spectral indexing).

! GETUV(sphere_dat,vrtnm,divnm,ug,vg): get U,V (u*cos(lat),v*cos(lat) on grid) 
! from spectral coefficients of vorticity and divergence. 
! Input:  sphere_dat, spectral coeffs. of vort and div (vrtnm,divnm).
! Output: gridded U,V (ug,vg).
!
! GETVRTDIV(sphere_dat,vrtnm,divnm,ug,vg): get spectral coefficients of 
! vorticity and divergence (vrtnm,divnm) from gridded U,V (ug,vg).
!
! COSGRAD(sphere_dat,chinm,uchig,vchig): compute cos(lat)*vector gradient.
! Inputs: sphere_dat, spectral coefficient array (chinm).
! Outputs: longitudinal and latitudinal components of gradient on the gaussian
! grid (uchig,vchig).

! SPECSMOOTH(sphere_dat,datagrid,smooth): isotropic spectral smoothing.
! Inputs: sphere_dat, smooth(ntrunc+1) (smoothing factor as a function
! of spherical harmonic degree), datagrid (gridded data to be smoothed).
! Outputs: datagrid (smoothed gridded data).
!
! SUMNM(sphere_dat,am,bm,anm,isign1,isign2):  
! given the arrays of fourier coeffs, am and bm,
! compute the complex spherical harmonic coeffs (anm) of:  
! isign1*( (1./rsphere*(1.-x**2))*d(ag)/d(lon) + (isign2/rsphere)*d(bg)/dx )
! where ag and bg are the grid pt counterparts of am,bm,
! isign1,isign2 are +1 or -1, rsphere is radius of sphere, x=sin(lat))
! am,bm can be computed from gridded data (ag,bg) using RFFT.
!
! RFFT(sphere_dat, data, coeff, idir): computes fourier harmonics in zonal
! direction of a gridded array. idir=+1 for forward (grid
! to fourier) and -1 for backward (fourier to grid) transform.
! data (nlon,nlat) contains gridded data, coeff (ntrunc+1,nlat) contains
! complex zonal fourier harmonics.  Uses FFTW library.
! 
! Important Details:
!
! The gridded data is assumed to be oriented such that i=1 is the Greenwich
! meridian and j=1 is the northernmost point. Grid indices increase eastward 
! and southward. The grid increment in longitude is 2*pi/nlon radians.
! For example nlon = 72 for a five degree grid. nlon must be greater than or
! equal to 4. The efficiency of the computation is improved when nlon is a
! product of small prime numbers.

! The spectral data is assumed to be in a complex array of dimension
! (NTRUNC+1)*(NTRUNC+2)/2. NTRUNC is the triangular truncation limit (NTRUNC=42
! for T42). NTRUNC must be <= nlon/2. Coefficients are ordered so that first
! (nm=1) is m=0,n=0, second is m=0,n=1, nm=mtrunc is m=0,n=mtrunc, nm=mtrunc+1 
! is m=1,n=1, etc. In Fortran90 syntax, values of m (degree) and n (order)
! corresponding as a function of the index nm are:

! indxm = (/((m,n=m,mtrunc),m=0,mtrunc)/)
! indxn = (/((n,n=m,mtrunc),m=0,mtrunc)/)

! Conversely, the index nm as a function of m and n is:

! nm = sum((/(i,i=mtrunc+1,mtrunc-m+2,-1)/))+n-m+1

! The associated legendre polynomials are normalized so that the integral  
! (pbar(n,m,theta)**2)*sin(theta) on the interval theta=0 to theta=pi is 1, 
! where: pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
! *sin(theta)**m/(2**n*factorial(n)) times the (n+m)th derivative of 
! (x**2-1)**n with respect to x=cos(theta).

! note: theta = 0.5*pi - phi, where phi is latitude and theta is colatitude,
! Therefore, cos(theta) = sin(phi) and sin(theta) = cos(phi).

! Note that pbar(0,0,theta)=sqrt(2)/2, 
! and pbar(1,0,theta)=0.5*sqrt(6.)*sin(lat).
!----------------------------------------------------------------------------

! explicit typing
! everything private to module, unless otherwise specified.

 implicit none
 private

 public :: sphere,spharmt_init,spharmt_destroy,spharm,rfft,cosgrad,&
           specsmooth,getuv,getvrtdiv,sumnm,gaulw,legend
 
 type sphere

! rsphere is radius of sphere in m.
   real    :: rsphere 

! nlons is number of longitudes (not including cyclic point).
   integer :: nlons

! nlats is number of gaussian latitudes.
   integer :: nlats

! for fft.
   real, dimension(:), pointer :: trigs
   integer, dimension(:), pointer :: ifax

! ntrunc is triangular truncation limit (e.g. 42 for T42 truncation)
   integer :: ntrunc

! gaulats is sin(gaussian lats), weights are gaussian weights,
! gwrc is gaussian weights divided by re*cos(lat)**2, lap is 
! the laplacian operatore -n*(n+1)/re**2 (where n is the degree 
! of the spherical harmonic),
! ilap is the inverse laplacian (1/lap, set to zero for n=0).
   real, dimension(:), pointer :: gaulats, weights, gwrc, lap, ilap

! pnm is associated legendre polynomials, hnm is -(1-x**2)d(pnm)/dx,
! where x = gaulats.
   real, dimension(:,:), pointer :: pnm, hnm

! indxm is zonal wavenumber index, indxn is spherical harmonic degree index.
   integer, dimension(:), pointer :: indxm, indxn
   logical :: isinitialized

 end type sphere

 contains

 subroutine spharmt_init(sphere_dat,nlon,nlat,ntrunc,re)

! initialize a sphere object.

   integer, intent(in) :: nlon,nlat,ntrunc
   real, intent(in) :: re
   type (sphere), intent(inout) :: sphere_dat
   real, dimension(:), allocatable :: pnm_tmp,hnm_tmp
   double precision, dimension(:), allocatable :: gaulats_tmp,weights_tmp
   integer :: nmdim,m,n,j

   sphere_dat%nlons = nlon
   sphere_dat%nlats = nlat
   sphere_dat%ntrunc = ntrunc
   sphere_dat%rsphere = re

   nmdim = (ntrunc+1)*(ntrunc+2)/2

   allocate(gaulats_tmp(nlat))
   allocate(weights_tmp(nlat))
   call gaulw(gaulats_tmp,weights_tmp)
   allocate(sphere_dat%gwrc(nlat))
   sphere_dat%gwrc = weights_tmp/(dble(re)*(1.d0-gaulats_tmp**2))
   allocate(sphere_dat%gaulats(nlat))
   allocate(sphere_dat%weights(nlat))
   sphere_dat%gaulats = gaulats_tmp
   sphere_dat%weights = weights_tmp
   deallocate(weights_tmp)

   allocate(sphere_dat%indxm(nmdim))
   allocate(sphere_dat%indxn(nmdim))
   allocate(sphere_dat%lap(nmdim))
   allocate(sphere_dat%ilap(nmdim))

   sphere_dat%indxm = (/((m,n=m,ntrunc),m=0,ntrunc)/)
   sphere_dat%indxn = (/((n,n=m,ntrunc),m=0,ntrunc)/)
   sphere_dat%lap(:)=-real(sphere_dat%indxn(:))*real(sphere_dat%indxn(:)+1)/re**2
   sphere_dat%ilap(1) = 0.
   sphere_dat%ilap(2:nmdim) = 1./sphere_dat%lap(2:nmdim)

   allocate(sphere_dat%pnm(nmdim,nlat))
   allocate(sphere_dat%hnm(nmdim,nlat))
   allocate(pnm_tmp(nmdim))
   allocate(hnm_tmp(nmdim))
   do j=1,nlat
      call legend(gaulats_tmp(j),pnm_tmp,hnm_tmp)
      sphere_dat%pnm(:,j) = pnm_tmp(:)
      sphere_dat%hnm(:,j) = hnm_tmp(:)
   enddo
   deallocate(gaulats_tmp)
   deallocate(pnm_tmp)
   deallocate(hnm_tmp)

   allocate(sphere_dat%trigs((3*nlon/2)+1))
   allocate(sphere_dat%ifax(13))
   call set99(sphere_dat%trigs,sphere_dat%ifax,nlon)

   sphere_dat%isinitialized = .TRUE.

 end subroutine spharmt_init

 subroutine spharmt_destroy(sphere_dat)

! deallocate pointers in sphere object.

   type (sphere), intent(inout) :: sphere_dat

   if (.not. sphere_dat%isinitialized) return

   deallocate(sphere_dat%gaulats)
   deallocate(sphere_dat%weights)
   deallocate(sphere_dat%gwrc)
   deallocate(sphere_dat%pnm)
   deallocate(sphere_dat%hnm)
   deallocate(sphere_dat%indxm)
   deallocate(sphere_dat%indxn)
   deallocate(sphere_dat%lap)
   deallocate(sphere_dat%ilap)
   deallocate(sphere_dat%trigs)
   deallocate(sphere_dat%ifax)

 end subroutine spharmt_destroy

 subroutine gaulw(sinlats, wts)

! compute sin of gaussian latitudes and gaussian weights.
! uses the iterative method presented in Numerical Recipes.

 double precision, intent (inout), dimension(:) :: sinlats, wts

 integer :: itermax
 integer :: i, iter, j, nlat, nprec
 double precision :: pi, pp, p1, p2, p3, z, z1, converg, ten


 ten = 10.d0
 converg = ten * EPSILON(ten)

 nprec = precision(converg)
 converg = .1**nprec

 itermax = 10

 pi = 4.0*ATAN(1.0)

 nlat = size(sinlats)
 if (size(sinlats) .ne. size(wts)) then
   print *, 'sinlats and wts must be same size in gaulw!'
   stop
 end if

 do i=1,nlat
   z = cos(pi*(i - 0.25)/(nlat + 0.5))
   do iter=1,itermax
     p1 = 1.0
     p2 = 0.0

     do j=1,nlat
        p3 = p2
        p2 = p1
        p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j
     end do

     pp = nlat*(z*p1 - p2)/(z*z - 1.0E+00)
     z1 = z
     z  = z1 - p1/pp
     if(ABS(z - z1) .LT. converg) go to 10
   end do
   print *, 'abscissas failed to converge in itermax iterations'
   print *, 'stopping in gaulw!'
   stop

10 continue

   sinlats (i)     = z
   wts (i)     = 2.0/((1.0 - z*z)*pp*pp)

 end do

 return
 end subroutine gaulw

 SUBROUTINE LEGEND(X,PMN,HMN)
!
!     THIS SUBROUTINE COMPUTES ASSOCIATED LEGENDRE  
!     FUNCTIONS, PMN, AND THE DERIVATIVE QUANTITY 
!     HMN = -(1 - X*X) * D(PMN)/DX  AT X = COS( COLATITUDE )  
!     GIVEN SOME COLATITUDE IN THE DOMAIN.
!
!     ACCURACY: 
!
!     THE RECURSION RELATION EMPLOYED IS LESS ACCURATE
!     NEAR THE POLES FOR M + N .GE. ROUGHLY 75 BECAUSE THE RECURSION
!     INVOLVES THE DIFFERENCE OF NEARLY EQUAL NUMBERS, SO
!     THE ASSOCIATED LEGENDRE FUNCTIONS ARE COMPUTED USING DOUBLE
!     PRECISION ACCUMULATORS. 
!     SEE BELOUSOV, TABLES OF NORMALIZED ASSOCIATED LEGENDRE  
!     POLYNOMIALS, C. 1962, FOR INFORMATION ON THE ACCURATE 
!     COMPUTATION OF ASSOCIATED LEGENDRE FUNCTIONS.
 
      real, dimension(:), intent(inout) ::  PMN,  HMN
      double precision, intent(in) :: x
      integer :: m,n,nm,i,nmax,np1,nmstrt,j,nmdim,ntrunc
      double precision :: A, B, PROD, SINSQ, &
      eps,epsnmp1,EPSPMN, PMNJ, PMNJM1, PMNJM2
 
 
!**** SET PARAMETERS FOR ENTRY INTO THE RECURSIVE FORMULAE. 
 
 
      SINSQ = 1.D0 - X * X
 
      A     = 1.D0
      B     = 0.D0
      PROD  = 1.D0
 
      nmdim = size(pmn)
      ntrunc = nint((-1.+sqrt(1+8*float(nmdim)))/2.)-1
      if (size(pmn) .ne. size(hmn)) then
         print *, 'pnm and hnm must be same size in subroutine legend!'
         stop
      end if
 
!**** LOOP FOR THE 'M' INDEX. 
 
      NMSTRT = 0
      DO I = 1, NTRUNC+1
 
         M    = (I - 1)
         NMAX = NTRUNC+1-M
 
!        WHEN M=0 (I=1), STANDARD LEGENDRE POLYNOMIALS ARE
!        GENERATED. 
 
         IF (M .NE. 0) THEN
               A    = A + 2.D0
               B    = B + 2.D0
               PROD = PROD * SINSQ * A / B 
         END IF
 
!****    GENERATE PMN AND HMN FOR J = 1 AND 2.
 
         PMNJM2   = SQRT(0.5D0 * PROD)
         NM = NMSTRT + 1
         PMN(NM) = PMNJM2
 
         PMNJM1   = SQRT( DBLE(2 * M + 3) ) * X * PMNJM2
         IF (NM .NE. NMDIM) PMN(NM+1) = PMNJM1
 
         NP1 = M + 1  
         EPSNMP1 = SQRT( DBLE(NP1*NP1 - M*M) / DBLE(4*NP1*NP1 - 1) )
         EPSPMN   = X * PMNJM1 - EPSNMP1 * PMNJM2 
 
         HMN(NM) = DBLE(M) * EPSNMP1 * PMNJM1 
         IF (NM .NE. NMDIM) &
         HMN(NM+1) = DBLE(M+1) * EPSPMN  -  &
         DBLE(M+2) * EPSNMP1 * PMNJM2
 
!****    LOOP FOR THE 'N' INDEX.
!        NOW APPLY THE RECURSION FORMULAE FOR J .GE. 3.
 
         DO J    = 3, NMAX  
            N = M + J - 1
            NM = NMSTRT + J
            EPS = SQRT( DBLE(N*N - M*M) / DBLE(4*N*N - 1) )
            PMNJ     = EPSPMN / EPS 
            PMN(NM) = PMNJ
 
!        COMPUTE EPS * PMN FOR J+1.
 
            EPSPMN   = X * PMNJ - EPS * PMNJM1  
 
!        COMPUTE THE DERIVATIVE.
 
            HMN(NM) = DBLE(N) * EPSPMN -  &
            DBLE(N+1) * EPS * PMNJM1 
 
            PMNJM2   = PMNJM1 
            PMNJM1   = PMNJ
         ENDDO  
         NMSTRT = NMSTRT + NMAX
     ENDDO

 END SUBROUTINE LEGEND

 subroutine rfft(sphere_dat, data, coeff, idir)

! real multiple fft (uses temperton fft991)

 type (sphere), intent(in) :: sphere_dat
 
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(inout) :: data
 real, dimension(sphere_dat%nlats*(sphere_dat%nlons+2)) :: wrk1
 real, dimension(sphere_dat%nlats*(sphere_dat%nlons+1)) :: wrk2
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats), intent(inout) :: coeff
 integer ::  nlons,nlats,ntrunc,mwaves,i,j,m,n
 integer, intent(in) :: idir

 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in rfft!'
    stop
 end if


 nlons = sphere_dat%nlons
 nlats = sphere_dat%nlats
 ntrunc = sphere_dat%ntrunc
 if (ntrunc .gt. nlons/2) then
    print *, 'ntrunc must be less than or equal to nlons in rfft'
    stop
 end if

 mwaves = ntrunc+1

!==> forward transform.

 if (idir .eq. +1) then
 
!==> copy the data into the work array.
!    transforms are computed pairwise using a complex fft.
     
   n = 0
   wrk1 = 0.
   do j=1,nlats
   do i=1,nlons+2  
      n = n + 1
      wrk1(n) = 0.0
      if (i .le. nlons) then
         wrk1(n) = data(i,j)
      end if
   enddo
   enddo
 
   call fft991(wrk1,wrk2,sphere_dat%trigs,sphere_dat%ifax,1,nlons+2,nlons,nlats,-1)
 
   n = -1
   do j=1,nlats
   do m=1,(nlons/2)+1
      n = n + 2
      if (m .le. mwaves) then
        coeff(m,j) = cmplx(wrk1(n),wrk1(n+1)) 
      end if
   enddo
   enddo

!==> inverse transform.

 else if (idir .eq. -1) then

   wrk1 = 0.
   n = -1
   do j=1,nlats
   do m=1,(nlons/2)+1
      n = n + 2
      if (m .le. mwaves) then
         wrk1(n) = real(coeff(m,j))
         wrk1(n+1) = aimag(coeff(m,j))
      end if
   enddo
   enddo
 
   call fft991(wrk1,wrk2,sphere_dat%trigs,sphere_dat%ifax,1,nlons+2,nlons,nlats,1)
 
   n = 0
   do j=1,nlats
   do i=1,nlons+2  
      n = n + 1
      if (i .le. nlons) then
         data(i,j) = wrk1(n)
      end if
   enddo
   enddo

 else
    write(6,*) ' idir must be +1 or -1 in RFFT!'
    write(6,*) ' execution terminated.'
    stop
 end if

 end subroutine rfft

 subroutine spharm(sphere_dat, ugrid, anm, idir)

! spherical harmonic transform

 type (sphere), intent(in) :: sphere_dat
 
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(inout) :: ugrid
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2), intent(inout) :: anm
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats) :: am
 integer ::  nlats,ntrunc,mwaves,nmstrt,nm,m,n,j
 integer, intent(in) :: idir

 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in spharm!'
    stop
 end if

 nlats = sphere_dat%nlats
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1

 IF (IDIR .eq. +1) then
 
!==>  GRID SPACE TO SPECTRAL SPACE TRANSFORMATION 
!     FIRST, INITIALIZE ARRAY.  
 
   anm = 0.
 
!==> perform ffts on each latitude.

   call rfft(sphere_dat, ugrid, am, 1)
 
!==>  SUM OVER ALL GAUSSIAN LATITUDES FOR EACH MODE AND EACH WAVE TO
!     OBTAIN THE TRANSFORMED VARIABLE IN SPECTRAL SPACE.
 
   do j=1,nlats
     NMSTRT = 0
     DO m = 1, mwaves
       DO n = 1, mwaves-m+1
          NM = NMSTRT + n
          anm(NM)=anm(NM)+sphere_dat%pnm(nm,j)*sphere_dat%weights(j)*am(m,j)
       enddo
       nmstrt = nmstrt + mwaves-m+1
     enddo
   enddo
 
!==>  SPECTRAL SPACE TO GRID SPACE TRANSFORMATION.
 
 else if (idir .eq. -1) then
 
 DO J = 1, NLATS 
 
!==>  INVERSE LEGENDRE TRANSFORM TO GET VALUES OF THE ZONAL FOURIER
!     TRANSFORM AT LATITUDE j.
 
!==>  SUM THE VARIOUS MERIDIONAL MODES TO BUILD THE FOURIER SERIES
!     COEFFICIENT FOR ZONAL WAVENUMBER M=I-1 AT the GIVEN LATITUDE. 
 
    NMSTRT = 0
    do m = 1, MWAVES
      am(m,j) = cmplx(0., 0.)
      DO n = 1, mwaves-m+1
         NM = NMSTRT + n
         am(m,j) = am(m,j)  +  anm(NM) * sphere_dat%pnm(NM,j) 
      enddo
      NMSTRT = NMSTRT + mwaves-m+1
    enddo
 
 ENDDO

!==>  FOURIER TRANSFORM TO COMPUTE THE VALUES OF THE VARIABLE IN GRID 
!     SPACE at THE J-TH LATITUDE. 

   call rfft(sphere_dat, ugrid, am, -1)
 
 else
   print *, 'error in spharm: idir must be -1 or +1!'
   print *, 'execution terminated in subroutine spharm'
 end if

 end subroutine spharm

 subroutine getuv(sphere_dat,vrtnm,divnm,ug,vg)

! compute U,V (winds times cos(lat)) from vrtnm,divnm
! (spectral coeffs of vorticity and divergence).

 type (sphere), intent(in) :: sphere_dat
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(out) ::  ug,vg
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2), intent(in) :: vrtnm,divnm
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats)  :: um,vm
 integer :: nlats,ntrunc,mwaves,m,j,n,nm,nmstrt
 real ::  rm

 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in getuv!'
    stop
 end if


 nlats = sphere_dat%nlats
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1

 do j=1,nlats
    nmstrt = 0
    do m=1,mwaves
       rm = m-1
       um(m,j) = cmplx(0.,0.)
       vm(m,j) = cmplx(0.,0.)
       do n=1,mwaves-m+1
          nm = nmstrt + n
          um(m,j) = um(m,j) + (sphere_dat%ilap(nm)/sphere_dat%rsphere)*( &
                    cmplx(0.,rm)*divnm(nm)*sphere_dat%pnm(nm,j) + &
                    vrtnm(nm)*sphere_dat%hnm(nm,j) )
          vm(m,j) = vm(m,j) + (sphere_dat%ilap(nm)/sphere_dat%rsphere)*( &
                    cmplx(0.,rm)*vrtnm(nm)*sphere_dat%pnm(nm,j) - &
                    divnm(nm)*sphere_dat%hnm(nm,j) )
       enddo
       nmstrt = nmstrt + mwaves-m+1
    enddo
 enddo

 call rfft(sphere_dat, ug, um, -1)
 call rfft(sphere_dat, vg, vm, -1)
 
 end subroutine getuv

 subroutine getvrtdiv(sphere_dat,vrtnm,divnm,ug,vg)

! compute vrtnm,divnm (spectral coeffs of vorticity and
! divergence) from U,V (winds time cos(lat)).

 type (sphere), intent(in) :: sphere_dat
 
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(inout) ::  ug,vg
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2), intent(in) :: vrtnm,divnm
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats) :: um,vm
 integer :: nlats,ntrunc,mwaves

 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in getvrtdiv!'
    stop
 end if

 ntrunc = sphere_dat%ntrunc
 nlats = sphere_dat%nlats
 mwaves = ntrunc+1

 call rfft(sphere_dat, ug, um, 1)
 call rfft(sphere_dat, vg, vm, 1)

 call sumnm(sphere_dat,um,vm,divnm,1,1)
 call sumnm(sphere_dat,vm,um,vrtnm,1,-1)
 
 end subroutine getvrtdiv

 subroutine sumnm(sphere_dat,am,bm,anm,isign1,isign2)
!
!  given the arrays of fourier coeffs, am and bm,
!  compute the complex spherical harmonic coeffs of:  
!
!  isign1*( (1./re*(1.-x**2))*d(ag)/d(lon) + (isign2/re)*d(bg)/dx )
!
!  where x = sin(lat), isign1 and isign2 are either +1 or -1,
!  ag and bg are the physical space values of am,bm and re
!  is the radius of the sphere.
!
!  the result is returned in anm.
!
!  for example on how to use this routine, see subroutine getvrtdiv.
!
!
 type (sphere), intent(in) :: sphere_dat
 
 integer, intent(in) :: isign1,isign2
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2) :: anm
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats), intent(in) :: am,bm
 integer :: nlats,ntrunc,mwaves,j,m,n,nm,nmstrt
 real ::  sign1,sign2,rm

 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in sumnm!'
    stop
 end if


 sign1 = float(isign1)
 sign2 = float(isign2)
 if (isign2 .ne. 1 .and. isign2 .ne. -1) then
    print *, ' isign2 must be +1 or -1 in sumnm!'
    print *, ' execution terminated in sumnm'
 end if
 if (isign1 .ne. 1 .and. isign1 .ne. -1) then
    print *, ' isign1 must be +1 or -1 in sumnm!'
    print *, ' execution terminated in sumnm'
    stop
 end if
 nlats = sphere_dat%nlats
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 
 anm = 0.
 DO J=1,NLATS
    NMSTRT = 0 
    DO m = 1, MWAVES 
       RM = M-1 
       DO n   = 1, mwaves-m+1
          NM = NMSTRT + n
          aNM(NM) = aNM(NM) + sign1*sphere_dat%GWrC(j)*(CMPLX(0.,RM) &
                    * sphere_dat%PNM(NM,J) * am(m,j) &
                    + sign2 * sphere_dat%HNM(NM,J) * bm(m,j))
       enddo
       nmstrt = nmstrt + mwaves - m +1
    enddo
 enddo

 end subroutine sumnm

 subroutine cosgrad(sphere_dat,divnm,ug,vg)

! compute coslat * gradient of spectral coefficients (divnm)
! vector gradient returned on grid as (ug,vg)

 type (sphere), intent(in) :: sphere_dat
 
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(out) ::  ug,vg
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2), intent(in) :: divnm
 complex, dimension(sphere_dat%ntrunc+1,sphere_dat%nlats) :: um,vm
 integer :: nlats,ntrunc,mwaves,j,m,n,nm,nmstrt
 real :: rm

 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in cosgrad!'
    stop
 end if


 nlats = sphere_dat%nlats
 ntrunc = sphere_dat%ntrunc
 mwaves = ntrunc+1
 
 do j=1,nlats
    nmstrt = 0
    do m=1,mwaves
       rm = (m-1)
       um(m,j) = cmplx(0.,0.)
       vm(m,j) = cmplx(0.,0.)
       do n=1,mwaves-m+1
          nm = nmstrt + n
          um(m,j) = um(m,j) + (1./sphere_dat%rsphere)* &
                    cmplx(0.,rm)*divnm(nm)*sphere_dat%pnm(nm,j) 
          vm(m,j) = vm(m,j) - (1./sphere_dat%rsphere)* &
                    divnm(nm)*sphere_dat%hnm(nm,j) 
       enddo
       nmstrt = nmstrt + mwaves - m +1
    enddo
 enddo

 call rfft(sphere_dat, ug, um, -1)
 call rfft(sphere_dat, vg, vm, -1)

 end subroutine cosgrad

 subroutine specsmooth(sphere_dat,datagrid,smooth)
 
! isotropic spectral smoothing of datagrid.
! input: smooth(sphere_dat%ntrunc+1) - smoothing factor as a
! function of degree (sphere_dat%indxn).

 type (sphere), intent(in) :: sphere_dat
 real, dimension(sphere_dat%nlons,sphere_dat%nlats), intent(inout) :: datagrid
 real, dimension(sphere_dat%ntrunc+1), intent(in) :: smooth
 complex, dimension((sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2) :: dataspec
 integer :: n,nm,nmdim

 if (.not. sphere_dat%isinitialized) then
    print *, 'uninitialized sphere object in specsmooth!'
    stop
 end if


 nmdim = (sphere_dat%ntrunc+1)*(sphere_dat%ntrunc+2)/2

 call spharm(sphere_dat, datagrid, dataspec, 1)

 do nm=1,nmdim
    n = sphere_dat%indxn(nm)
    dataspec(nm) = dataspec(nm)*smooth(n+1)
 enddo

 call spharm(sphere_dat, datagrid, dataspec, -1)

 end subroutine specsmooth

    subroutine fft99 (a,work,trigs,ifax,inc,jump,n,lot,isign) 

! purpose      performs multiple fast fourier transforms.  this package
!              will perform a number of simultaneous real/half-complex
!              periodic fourier transforms or corresponding inverse
!              transforms, i.e.  given a set of real data vectors, the
!              package returns a set of 'half-complex' fourier
!              coefficient vectors, or vice versa.  the length of the
!              transforms must be an even number greater than 4 that has
!              no other factors except possibly powers of 2, 3, and 5.
!              this is an all-fortran version of a optimized routine
!              fft99 written for xmp/ymps by dr. clive temperton of
!              ecmwf.
!
!              the package fft99f contains several user-level routines:
!
!            subroutine set99
!                an initialization routine that must be called once
!                before a sequence of calls to the fft routines
!                (provided that n is not changed).
!
!            subroutines fft99 and fft991
!                two fft routines that return slightly different
!                arrangements of the data in gridpoint space.
!
! usage        let n be of the form 2**p * 3**q * 5**r, where p .ge. 1,
!              q .ge. 0, and r .ge. 0.  then a typical sequence of
!              calls to transform a given set of real vectors of length
!              n to a set of 'half-complex' fourier coefficient vectors
!              of length n is
!
!                   dimension ifax(13),trigs(3*n/2+1),a(m*(n+2)),
!                  +          work(m*(n+1))
!
!                   call set99 (trigs, ifax, n)
!                   call fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
!
!              see the individual write-ups for set99, fft99, and
!              fft991 below, for a detailed description of the
!              arguments.
!
! history      the package was written by clive temperton at ecmwf in
!              november, 1978.  it was modified, documented, and tested
!              for ncar by russ rew in september, 1980.
!
!-----------------------------------------------------------------------
!
! subroutine set99 (trigs, ifax, n)
!
! purpose      a set-up routine for fft99 and fft991.  it need only be
!              called once before a sequence of calls to the fft
!              routines (provided that n is not changed).
!
! argument     ifax(13),trigs(3*n/2+1)
! dimensions
!
! arguments
!
! on input     trigs
!               a floating point array of dimension 3*n/2 if n/2 is
!               even, or 3*n/2+1 if n/2 is odd.
!
!              ifax
!               an integer array.  the number of elements actually used
!               will depend on the factorization of n.  dimensioning
!               ifax for 13 suffices for all n less than a million.
!
!              n
!               an even number greater than 4 that has no prime factor
!               greater than 5.  n is the length of the transforms (see
!               the documentation for fft99 and fft991 for the
!               definitions of the transforms).
!
! on output    ifax
!               contains the factorization of n/2.  ifax(1) is the
!               number of factors, and the factors themselves are stored
!               in ifax(2),ifax(3),...  if set99 is called with n odd,
!               or if n has any prime factors greater than 5, ifax(1)
!               is set to -99.
!
!              trigs
!               an array of trigonometric function values subsequently
!               used by the fft routines.
!
!-----------------------------------------------------------------------
!
! subroutine fft991 (a,work,trigs,ifax,inc,jump,n,m,isign)
!                       and
! subroutine fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
!
! purpose      perform a number of simultaneous real/half-complex
!              periodic fourier transforms or corresponding inverse
!              transforms, using ordinary spatial order of gridpoint
!              values (fft991) or explicit cyclic continuity in the
!              gridpoint values (fft99).  given a set
!              of real data vectors, the package returns a set of
!              'half-complex' fourier coefficient vectors, or vice
!              versa.  the length of the transforms must be an even
!              number that has no other factors except possibly powers
!              of 2, 3, and 5.  this is an all-fortran version of 
!              optimized routine fft991 written for xmp/ymps by
!              dr. clive temperton of ecmwf.
!
! argument     a(m*(n+2)), work(m*(n+1)), trigs(3*n/2+1), ifax(13)
! dimensions
!
! arguments
!
! on input     a
!               an array of length m*(n+2) containing the input data
!               or coefficient vectors.  this array is overwritten by
!               the results.
!
!              work
!               a work array of dimension m*(n+1)
!
!              trigs
!               an array set up by set99, which must be called first.
!
!              ifax
!               an array set up by set99, which must be called first.
!
!              inc
!               the increment (in words) between successive elements of
!               each data or coefficient vector (e.g.  inc=1 for
!               consecutively stored data).
!
!              jump
!               the increment (in words) between the first elements of
!               successive data or coefficient vectors.  on crays, 
!               try to arrange data so that jump is not a multiple of 8
!               (to avoid memory bank conflicts).  for clarification of
!               inc and jump, see the examples below.
!
!              n
!               the length of each transform (see definition of
!               transforms, below).
!
!              m
!               the number of transforms to be done simultaneously.
!
!              isign
!               = +1 for a transform from fourier coefficients to
!                    gridpoint values.
!               = -1 for a transform from gridpoint values to fourier
!                    coefficients.
!
! on output    a
!               if isign = +1, and m coefficient vectors are supplied
!               each containing the sequence:
!
!               a(0),b(0),a(1),b(1),...,a(n/2),b(n/2)  (n+2 values)
!
!               then the result consists of m data vectors each
!               containing the corresponding n+2 gridpoint values:
!
!               for fft991, x(0), x(1), x(2),...,x(n-1),0,0.
!               for fft99, x(n-1),x(0),x(1),x(2),...,x(n-1),x(0).
!                   (explicit cyclic continuity)
!
!               when isign = +1, the transform is defined by:
!                 x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
!                 where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
!                 and i=sqrt (-1)
!
!               if isign = -1, and m data vectors are supplied each
!               containing a sequence of gridpoint values x(j) as
!               defined above, then the result consists of m vectors
!               each containing the corresponding fourier cofficients
!               a(k), b(k), 0 .le. k .le n/2.
!
!               when isign = -1, the inverse transform is defined by:
!                 c(k)=(1/n)*sum(j=0,...,n-1)(x(j)*exp(-2*i*j*k*pi/n))
!                 where c(k)=a(k)+i*b(k) and i=sqrt(-1)
!
!               a call with isign=+1 followed by a call with isign=-1
!               (or vice versa) returns the original data.
!
!               note: the fact that the gridpoint values x(j) are real
!               implies that b(0)=b(n/2)=0.  for a call with isign=+1,
!               it is not actually necessary to supply these zeros.
!
! examples      given 19 data vectors each of length 64 (+2 for explicit
!               cyclic continuity), compute the corresponding vectors of
!               fourier coefficients.  the data may, for example, be
!               arranged like this:
!
! first data   a(1)=    . . .                a(66)=             a(70)
! vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
!
! second data  a(71)=   . . .                                  a(140)
! vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
!
!               and so on.  here inc=1, jump=70, n=64, m=19, isign=-1,
!               and fft99 should be used (because of the explicit cyclic
!               continuity).
!
!               alternatively the data may be arranged like this:
!
!                first         second                          last
!                data          data                            data
!                vector        vector                          vector
!
!                 a(1)=         a(2)=                           a(19)=
!
!                 x(63)         x(63)       . . .               x(63)
!        a(20)=   x(0)          x(0)        . . .               x(0)
!        a(39)=   x(1)          x(1)        . . .               x(1)
!                  .             .                               .
!                  .             .                               .
!                  .             .                               .
!
!               in which case we have inc=19, jump=1, and the remaining
!               parameters are the same as before.  in either case, each
!               coefficient vector overwrites the corresponding input
!               data vector.
!
!-----------------------------------------------------------------------
    integer, intent(in)    :: inc,jump,n,lot,isign
    integer, intent(in) :: ifax(:)
    real,    intent(in)    :: trigs(:)
    real,    intent(inout) :: a(*),work(*)

!     dimension a(n),work(n),trigs(n),ifax(1)
!
!     subroutine "fft99" - multiple fast real periodic transform
!     corresponding to old scalar routine fft9
!     procedure used to convert to half-length complex transform
!     is given by cooley, lewis and welch (j. sound vib., vol. 12
!     (1970), 315-337)
!
!     a is the array containing input and output data
!     work is an area of size (n+1)*lot
!     trigs is a previously prepared list of trig function values
!     ifax is a previously prepared list of factors of n/2
!     inc is the increment within each data 'vector'
!         (e.g. inc=1 for consecutively stored data)
!     jump is the increment between the start of each data vector
!     n is the length of the data vectors
!     lot is the number of data vectors
!     isign = +1 for transform from spectral to gridpoint
!           = -1 for transform from gridpoint to spectral
!
!     ordering of coefficients:
!         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
!         where b(0)=b(n/2)=0; (n+2) locations required
!
!     ordering of data:
!         x(n-1),x(0),x(1),x(2),...,x(n),x(0)
!         i.e. explicit cyclic continuity; (n+2) locations required
!
!     vectorization is achieved on cray by doing the transforms in
!     parallel
!
!     *** n.b. n is assumed to be an even number
!
!     definition of transforms:
!     -------------------------
!
!     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
!         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
!
!     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
!               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
!

    integer :: nfax, nx, nh, ink, igo, ibase, jbase
    integer :: i, j, k, L, m, ia, la, ib

      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30

!   if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=inc+1
      jbase=1
    do L=1,lot
      i=ibase
      j=jbase
      do m=1,n
        work(j)=a(i)
        i=i+inc
        j=j+1
      enddo
      ibase=ibase+jump
      jbase=jbase+nx
    enddo

      igo=60
      go to 40

!   preprocessing (isign=+1)
!   ------------------------

   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60

!   complex transform
!   -----------------

   40 continue
      ia=inc+1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm (a(ia),a(ia+inc),work(1),work(2),trigs, &
                   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm (work(1),work(2),a(ia),a(ia+inc),trigs, &
                   2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue

    if (isign.eq.-1) go to 130

! if necessary, transfer data from work area

    if (mod(nfax,2).ne.1) then
      ibase=1
      jbase=ia
      do L=1,lot
        i=ibase
        j=jbase
        do m=1,n
          a(j)=work(i)
          i=i+1
          j=j+inc
        enddo
        ibase=ibase+nx
        jbase=jbase+jump
      enddo
    endif

!   fill in cyclic boundary points
      ia=1
      ib=n*inc+1
      do L=1,lot
        a(ia)=a(ib)
        a(ib+inc)=a(ia+inc)
        ia=ia+jump
        ib=ib+jump
      enddo
      go to 140

!   postprocessing (isign=-1):
!   --------------------------

  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)

  140 continue

    end subroutine fft99

!##########################################################################

    subroutine fft99a (a,work,trigs,inc,jump,n,lot)
    integer, intent(in)    :: inc,jump,n,lot
    real,    intent(in)    :: trigs(:)
    real,    intent(inout) :: a(*),work(*)

!     dimension a(n),work(n),trigs(n)
!
!     subroutine fft99a - preprocessing step for fft99, isign=+1
!     (spectral to gridpoint transform)

    integer :: nh, nx, ink, k, L
    integer :: ia, ib, ja, jb, iabase, ibbase, jabase, jbbase
    real    :: c, s

      nh=n/2
      nx=n+1
      ink=inc+inc

!   a(0) and a(n/2)
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
    do L=1,lot
      work(ja)=a(ia)+a(ib)
      work(jb)=a(ia)-a(ib)
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
    enddo
 
!   remaining wavenumbers
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1

    do k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
      do L=1,lot
        work(ja)=(a(ia)+a(ib))- &
            (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
        work(jb)=(a(ia)+a(ib))+ &
            (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
        work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+ &
            (a(ia+inc)-a(ib+inc))
        work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))- &
            (a(ia+inc)-a(ib+inc))
        ia=ia+jump
        ib=ib+jump
        ja=ja+nx
        jb=jb+nx
      enddo
      iabase=iabase+ink
      ibbase=ibbase-ink
      jabase=jabase+2
      jbbase=jbbase-2
    enddo

!   wavenumber n/4 (if it exists)
    if (iabase.eq.ibbase) then
      ia=iabase
      ja=jabase
      do L=1,lot
        work(ja)=2.0*a(ia)
        work(ja+1)=-2.0*a(ia+inc)
        ia=ia+jump
        ja=ja+nx
      enddo
    endif

    end subroutine fft99a

!##########################################################################

    subroutine fft99b (work,a,trigs,inc,jump,n,lot)
    integer, intent(in)    :: inc,jump,n,lot
    real,    intent(in)    :: trigs(:)
    real,    intent(inout) :: a(*),work(*)

!     dimension work(n),a(n),trigs(n)
!
!     subroutine fft99b - postprocessing step for fft99, isign=-1
!     (gridpoint to spectral transform)

    integer :: nh, nx, ink, k, L
    integer :: ia, ib, ja, jb, iabase, ibbase, jabase, jbbase
    real    :: scale, c, s

      nh=n/2
      nx=n+1
      ink=inc+inc

!   a(0) and a(n/2)
      scale=1.0/real(n)
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
    do L=1,lot
      a(ja)=scale*(work(ia)+work(ib))
      a(jb)=scale*(work(ia)-work(ib))
      a(ja+inc)=0.0
      a(jb+inc)=0.0
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
    enddo

!   remaining wavenumbers
      scale=0.5*scale
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1

    do k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
      do L=1,lot
        a(ja)=scale*((work(ia)+work(ib)) &
           +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
        a(jb)=scale*((work(ia)+work(ib)) &
           -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
        a(ja+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
            +(work(ib+1)-work(ia+1)))
        a(jb+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
            -(work(ib+1)-work(ia+1)))
        ia=ia+nx
        ib=ib+nx
        ja=ja+jump
        jb=jb+jump
      enddo
      iabase=iabase+2
      ibbase=ibbase-2
      jabase=jabase+ink
      jbbase=jbbase-ink
    enddo

!   wavenumber n/4 (if it exists)
    if (iabase.eq.ibbase) then
      ia=iabase
      ja=jabase
      scale=2.0*scale
      do L=1,lot
        a(ja)=scale*work(ia)
        a(ja+inc)=-scale*work(ia+1)
        ia=ia+nx
        ja=ja+jump
      enddo
    endif

    end subroutine fft99b

!##########################################################################

    subroutine fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
    integer, intent(in)    :: inc,jump,n,lot,isign
    integer, intent(in) :: ifax(:)
    real,    intent(in)    :: trigs(:)
    real,    intent(inout) :: a(*),work(*)

!     dimension a(n),work(n),trigs(n),ifax(1)
!
!     subroutine "fft991" - multiple real/half-complex periodic
!     fast fourier transform
!
!     same as fft99 except that ordering of data corresponds to
!     that in mrfft2
!
!     procedure used to convert to half-length complex transform
!     is given by cooley, lewis and welch (j. sound vib., vol. 12
!     (1970), 315-337)
!
!     a is the array containing input and output data
!     work is an area of size (n+1)*lot
!     trigs is a previously prepared list of trig function values
!     ifax is a previously prepared list of factors of n/2
!     inc is the increment within each data 'vector'
!         (e.g. inc=1 for consecutively stored data)
!     jump is the increment between the start of each data vector
!     n is the length of the data vectors
!     lot is the number of data vectors
!     isign = +1 for transform from spectral to gridpoint
!           = -1 for transform from gridpoint to spectral
!
!     ordering of coefficients:
!         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
!         where b(0)=b(n/2)=0; (n+2) locations required
!
!     ordering of data:
!         x(0),x(1),x(2),...,x(n-1)
!
!     vectorization is achieved on cray by doing the transforms in
!     parallel
!
!     *** n.b. n is assumed to be an even number
!
!     definition of transforms:
!     -------------------------
!
!     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
!         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
!
!     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
!               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
!

    integer :: nfax, nx, nh, ink, igo, ibase, jbase
    integer :: i, j, k, L, m, ia, la, ib


      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30

!     if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=1
      jbase=1
    do L=1,lot
      i=ibase
      j=jbase
      do m=1,n
        work(j)=a(i)
        i=i+inc
        j=j+1
      enddo
      ibase=ibase+jump
      jbase=jbase+nx
    enddo
!
      igo=60
      go to 40
!
!   preprocessing (isign=+1)
!   ------------------------
!
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
!
!   complex transform
!   -----------------
!
   40 continue
      ia=1
      la=1
    do k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
        call vpassm (a(ia),a(ia+inc),work(1),work(2),trigs, &
                     ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
        call vpassm (work(1),work(2),a(ia),a(ia+inc),trigs, &
                     2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
    enddo

    if (isign.eq.-1) go to 130

! if necessary, transfer data from work area
    if (mod(nfax,2).ne.1) then
      ibase=1
      jbase=1
      do L=1,lot
        i=ibase
        j=jbase
        do m=1,n
          a(j)=work(i)
          i=i+1
          j=j+inc
        enddo
        ibase=ibase+nx
        jbase=jbase+jump
      enddo
    endif

!   fill in zeros at end
      ib=n*inc+1
      do L=1,lot
        a(ib)=0.0
        a(ib+inc)=0.0
        ib=ib+jump
      enddo
      go to 140

!     postprocessing (isign=-1):
!     --------------------------

  130 continue
      call fft99b (work,a,trigs,inc,jump,n,lot)

  140 continue

    end subroutine fft991

!##########################################################################

    subroutine set99 (trigs, ifax, n)
    integer, intent(in)  :: n
    integer, intent(out) :: ifax(:)
    real,    intent(out) :: trigs(:)

!     dimension ifax(13),trigs(1)
!
! mode 3 is used for real/half-complex transforms.  it is possible
! to do complex/complex transforms with other values of mode, but
! documentation of the details were not available when this routine
! was written.
!
    integer :: mode = 3
    integer :: i

      call fax (ifax, n, mode)
      i = ifax(1)
      if (ifax(i+1) .gt. 5 .or. n .le. 4) ifax(1) = -99
      if (ifax(1) .le. 0 ) then 
        write(6,*) ' set99 -- invalid n'
        stop
      endif
      call fftrig (trigs, n, mode)

    end subroutine set99

!##########################################################################

    subroutine fax (ifax,n,mode)
    integer, intent(out) :: ifax(:)
    integer, intent(in)  :: n, mode

    integer :: nn, k, L, inc, nfax, ii, istop, i, item

      nn=n
      if (iabs(mode).eq.1) go to 10
      if (iabs(mode).eq.8) go to 10
      nn=n/2
      if ((nn+nn).eq.n) go to 10
      ifax(1)=-99
      return
   10 k=1
!     test for factors of 4
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
!     test for extra factor of 2
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
!     test for factors of 3
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
!     now find remaining factors
   50 L=5
      inc=2
!     inc alternately takes on values 2 and 4
   60 if (mod(nn,L).ne.0) go to 70
      k=k+1
      ifax(k)=L
      nn=nn/L
      if (nn.eq.1) go to 80
      go to 60
   70 L=L+inc
      inc=6-inc
      go to 60
   80 ifax(1)=k-1
!     ifax(1) contains number of factors
      nfax=ifax(1)
!     sort factors into ascending order
      if (nfax.eq.1) go to 110
      do 100 ii=2,nfax
      istop=nfax+2-ii
      do 90 i=2,istop
      if (ifax(i+1).ge.ifax(i)) go to 90
      item=ifax(i)
      ifax(i)=ifax(i+1)
      ifax(i+1)=item
   90 continue
  100 continue
  110 continue

    end subroutine fax

!##########################################################################

    subroutine fftrig (trigs,n,mode)
    real,    intent(out) :: trigs(:)
    integer, intent(in)  :: n, mode

    real    :: del, angle,pi
    integer :: imode, nn, nh, i, L, la

      pi = 4.*atan(1.0)
      imode=iabs(mode)
      nn=n
      if (imode.gt.1.and.imode.lt.6) nn=n/2
      del=(pi+pi)/real(nn)
      L=nn+nn
      do i=1,L,2
        angle=0.5*real(i-1)*del
        trigs(i)=cos(angle)
        trigs(i+1)=sin(angle)
      enddo
      if (imode.eq.1) return
      if (imode.eq.8) return

      del=0.5*del
      nh=(nn+1)/2
      L=nh+nh
      la=nn+nn
      do i=1,L,2
        angle=0.5*real(i-1)*del
        trigs(la+i)=cos(angle)
        trigs(la+i+1)=sin(angle)
      enddo
      if (imode.le.3) return

      del=0.5*del
      la=la+nn
    if (mode.ne.5) then
      do i=2,nn
        angle=real(i-1)*del
        trigs(la+i)=2.0*sin(angle)
      enddo
      return
    endif

      del=0.5*del
      do i=2,n
        angle=real(i-1)*del
        trigs(la+i)=sin(angle)
      enddo

    end subroutine fftrig

!##########################################################################

    subroutine vpassm (a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
    integer, intent(in)  :: inc1, inc2, inc3, inc4, lot, n, ifac, la
    real,    intent(in)  :: a(*),b(*),trigs(*)
    real,    intent(out) :: c(*),d(*)
!
! the following two lines are the original declarations of the
! first 5 arguments.  this routine is called in at least one place
! with arrays which are larger in size than n, and in this subroutine
! indices larger than n are generated.  this causes bounds warnings
! on some compilers and core dumps on others.  when i changed the arrays
! to be assumed shape (*) the problem went away, and the results look ok.
! n. collins  15jun06
!
!   real,    intent(in)  :: a(n),b(n),trigs(n)
!   real,    intent(out) :: c(n),d(n)

!
!     subroutine "vpassm" - multiple version of "vpassa"
!     performs one pass through data
!     as part of multiple complex fft routine
!     a is first real input vector
!     b is first imaginary input vector
!     c is first real output vector
!     d is first imaginary output vector
!     trigs is precalculated table of sines " cosines
!     inc1 is addressing increment for a and b
!     inc2 is addressing increment for c and d
!     inc3 is addressing increment between a"s & b"s
!     inc4 is addressing increment between c"s & d"s
!     lot is the number of vectors
!     n is length of vectors
!     ifac is current factor of n
!     la is product of previous factors
!

    real :: sin36=0.587785252292473
    real :: cos36=0.809016994374947
    real :: sin72=0.951056516295154
    real :: cos72=0.309016994374947
    real :: sin60=0.866025403784437

    integer :: i, j, k, L, m, iink, jink, jump, ibase, jbase, igo, ijk, la1
    integer :: ia, ja, ib, jb, kb, ic, jc, kc, id, jd, kd, ie, je, ke
    real    :: c1, s1, c2, s2, c3, s3, c4, s4

      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo.gt.4) return
!del  go to (10,50,90,130),igo

  select case (igo)

!   coding for factor 2

    case (1)
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      do 20 L=1,la
      i=ibase
      j=jbase
      do 15 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
      i=i+inc3
      j=j+inc4
   15 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   20 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 40 k=la1,m,la
      kb=k+k-2
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      do 30 L=1,la
      i=ibase
      j=jbase
      do 25 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
      i=i+inc3
      j=j+inc4
   25 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   30 continue
      jbase=jbase+jump
   40 continue
!     return

!   coding for factor 3

    case (2)
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      do 60 L=1,la
      i=ibase
      j=jbase
      do 55 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
      i=i+inc3
      j=j+inc4
   55 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   60 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 80 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      do 70 L=1,la
      i=ibase
      j=jbase
      do 65 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=                                                           &
          c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
         -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=                                                           &
          s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
         +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=                                                           &
          c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
         -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=                                                           &
          s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
         +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
   65 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   70 continue
      jbase=jbase+jump
   80 continue
!     return

!   coding for factor 4

    case (3)
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      do 100 L=1,la
      i=ibase
      j=jbase
      do 95 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
      i=i+inc3
      j=j+inc4
   95 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  100 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 120 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      do 110 L=1,la
      i=ibase
      j=jbase
      do 105 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)=                                     &
          c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
         -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)=                                     &
          s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
         +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)=                                     &
          c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
         -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=                                     &
          s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
         +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=                                     &
          c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
         -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=                                     &
          s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
         +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  105 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  110 continue
      jbase=jbase+jump
  120 continue
!     return

!   coding for factor 5

    case (4)
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
      do 140 L=1,la
      i=ibase
      j=jbase
      do 135 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
        -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
        +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
        +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
        -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
        -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
        +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
        +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
        -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  135 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  140 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 160 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      do 150 L=1,la
      i=ibase
      j=jbase
      do 145 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=                                                          &
          c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(jb+j)=                                                          &
          s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(je+j)=                                                          &
          c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(je+j)=                                                          &
          s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(jc+j)=                                                          &
          c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jc+j)=                                                          &
          s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      c(jd+j)=                                                          &
          c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jd+j)=                                                          &
          s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      i=i+inc3
      j=j+inc4
  145 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  150 continue
      jbase=jbase+jump
  160 continue

  end select

    end subroutine vpassm

 end module spharmt_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

      SUBROUTINE HINES_DSP1 (DRAG_U,DRAG_V,HEAT,DIFFCO,FLUX_U,FLUX_V,
     1                       VEL_U,VEL_V,BVFREQ,DENSITY,VISC_MOL,ALT,
     2                       RMS_WIND,K_ALPHA,V_ALPHA,M_ALPHA,
     3                       SIGMA_ALPHA,SIGSQH_ALPHA,AK_ALPHA,
     &                       SPECFAC,DO_ALPHA,DRAG,DRAGIL,  
     4                       MMIN_ALPHA,I_ALPHA,SIGMA_T,DENSB,BVFB,
     5                       UBOT,VBOT,IHEATCAL,ICUTOFF,IPRNT1,IPRNT2,
     6                       NSMAX,SMCO,ALT_CUTOFF,KSTAR,M_MIN,SLOPE,
     7                       F1,F2,F3,F5,F6,NAZ,IL1,IL2,
     8                       LEVBOT,LEVTOP,NLONS,NLEVS,NAZMTH)
C
C  Main routine for Hines Doppler spread gravity wave parameterization
C  scheme which calculates zonal and meridional components of gravity 
C  wave drag, heating rates and diffusion coefficient on a longitude 
C  by altitude grid.
C
C  Aug. 13/95 - C. McLandress
C
C  Modifications:
C  --------------
C  Feb. 2/96 - C. McLandress (added: minimum cutoff wavenumber M_MIN;
C                             12 and 16 azimuths; inclusion of introwaves
C                             by UBOT and VBOT; diffusion and heating
C                             calculated at half level; logical flags;
C                             V_ALPHA array only at single level;
C                             new print routine) 
C
C  Output arguements:
C  ------------------
C
C     * DRAG_U = zonal component of gravity wave drag (m/s^2).
C     * DRAG_V = meridional component of gravity wave drag (m/s^2).
C     * HEAT   = gravity wave heating (K/sec).
C     * DIFFCO = diffusion coefficient (m^2/sec)
C     * FLUX_U = zonal component of vertical momentum flux (Pascals)
C     * FLUX_V = meridional component of vertical momentum flux (Pascals)
C
C  Input arguements:
C  -----------------
C
C     * VEL_U      = background zonal wind component (m/s) (positive is
C     *              eastward).
C     * VEL_V      = background meridional wind component (m/s) (positive
C     *              northward).
C     * BVFREQ     = background Brunt Vassala frequency (radians/sec).
C     * DENSITY    = background density (kg/m^3) 
C     * VISC_MOL   = molecular viscosity (m^2/s)
C     * ALT        = altitude of momentum, density, buoyancy levels (m)
C     *              (levels must be ordered so ALT(I,LEVTOP) > ALT(I,LEVBOT))
C     * RMS_WIND   = root mean square gravity wave wind at bottom (reference)
C     *              level (m/s).
C     * K_ALPHA    = horizontal wavenumber of each azimuth (1/m).
C     * IHEATCAL   = 1 to calculate heating rates and diffusion coefficient.
C     * IPRNT1     = 1 to print out flux, drag arrays at specified longitudes.
C     * IPRNT2     = 1 to print out azimuthal arrays at specified longitudes.
C     * ICUTOFF    = 1 to exponentially damp drag, heating and diffusion 
C     *              arrays above the altitude ALT_CUTOFF.
C     * ALT_CUTOFF = altitude in meters above which exponential decay applied.
C     * SMCO       = smoothing factor used to smooth cutoff vertical 
C     *              wavenumbers and total rms winds in vertical direction
C     *              before calculating drag or heating
C     *              (SMCO >= 1 ==>  a 1:SMCO:1 three-point stencil is used).
C     * NSMAX      = number of times smoother applied ( >= 1),
C     *            = 0 means no smoothing performed.
C     * KSTAR      = typical gravity wave horizontal wavenumber (1/m).
C     * M_MIN      = minimum allowable cutoff vertical wavenumber, 
C     *              e.g., 1/(3km). This is used only for a spectral slope
C     *              of one ==> for slope of 1.5 or 2 then M_MIN = 0. 
C     * SLOPE      = slope of incident vertical wavenumber spectrum
C     *              (SLOPE must equal 1., 1.5 or 2.).
C     * F1 to F6   = Hines's fudge factors (F4 not needed since used for
C     *              vertical flux of vertical momentum).
C     * NAZ        = number of horizontal azimuths used (4, 8, 12 or 16).
C     * IL1        = first longitudinal index to use (IL1 >= 1).
C     * IL2        = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * LEVBOT     = index of bottom (reference) drag level.
C     * LEVTOP     = index of top drag level (if LEVBOT > LEVTOP then the
C     *              vertical indexes increase from top down, otherwise
C     *              indexes increase from bottom up).
C     * NLONS      = number of longitudes.
C     * NLEVS      = number of vertical levels.
C     * NAZMTH     = azimuthal array dimension (NAZMTH >= NAZ).
C 
C  Input work arrays:
C  ------------------
C
C     * M_ALPHA      = cutoff vertical wavenumber (1/m).
C     * V_ALPHA      = wind component at each azimuth (m/s).
C     * SIGMA_ALPHA  = total rms wind in each azimuth (m/s).
C     * SIGSQH_ALPHA = portion of wind variance from waves having wave
C     *                normals in the alpha azimuth (m/s).
C     * SIGMA_T      = total rms horizontal wind (m/s).
C     * AK_ALPHA     = spectral amplitude factor at each azimuth 
C     *                (i.e.,{AjKj}) in m^4/s^2.
C     * SPECFAC      = AK_ALPHA * K_ALPHA.
C     * DO_ALPHA     = logical flag indicating azimuths and longitudes
C     *                where calculations to be performed.
C     * DRAG         = logical flag indicating longitudes where calculations
C     *                to be performed.
C     * DRAGIL       = logical flag indicating longitudes and levels where 
C     *                calculations to be performed.
C     * I_ALPHA      = Hines' integral.
C     * MMIN_ALPHA   = minimum value of cutoff wavenumber.
C     * DENSB        = background density at reference level.
C     * BVFB         = buoyancy frequency at bottom reference and
C     *                work array for ICUTOFF = 1.
C     * UBOT         = background zonal wind component at reference level.
C     * VBOT         = background meridional wind component at reference level.
C
      INTEGER  NAZ, NLONS, NLEVS, NAZMTH, IL1, IL2, LEVBOT, LEVTOP
      INTEGER  ICUTOFF, NSMAX, IHEATCAL, IPRNT1, IPRNT2
      LOGICAL  DRAGIL(NLONS,NLEVS), DRAG(NLONS), DO_ALPHA(NLONS,NAZMTH)
      REAL  KSTAR, M_MIN, F1, F2, F3, F5, F6, SLOPE, ALT_CUTOFF, SMCO
      REAL  DRAG_U(NLONS,NLEVS),   DRAG_V(NLONS,NLEVS) 
      REAL  HEAT(NLONS,NLEVS),     DIFFCO(NLONS,NLEVS)
      REAL  FLUX_U(NLONS,NLEVS),   FLUX_V(NLONS,NLEVS)
      REAL  VEL_U(NLONS,NLEVS),    VEL_V(NLONS,NLEVS)
      REAL  BVFREQ(NLONS,NLEVS),   DENSITY(NLONS,NLEVS)
      REAL  VISC_MOL(NLONS,NLEVS), ALT(NLONS,NLEVS)
      REAL  RMS_WIND(NLONS),       BVFB(NLONS),   DENSB(NLONS)
      REAL  UBOT(NLONS),           VBOT(NLONS)
      REAL  SIGMA_T(NLONS,NLEVS)
      REAL  SIGMA_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  SIGSQH_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH), V_ALPHA(NLONS,NAZMTH)
      REAL  AK_ALPHA(NLONS,NAZMTH),      K_ALPHA(NLONS,NAZMTH)
      REAL  MMIN_ALPHA(NLONS,NAZMTH) ,   I_ALPHA(NLONS,NAZMTH)
      REAL  SPECFAC(NLONS,NAZMTH)
C
C  Internal variables.
C  -------------------
C
      INTEGER  IERROR, I, N, L, LEV1P, LEV2M, NUNIT
      INTEGER  ILPRT1, ILPRT2, LEV1, LEV2, IORDER, ICOUNT
      REAL  ZERO, RMS_MIN
      DATA  ZERO / 0. /, RMS_MIN / 0.001 / 
C----------------------------------------------------------------------- 
C
C  Check that things set up correctly, abort if not.
C 
      IERROR = 0
      IF (NAZ .GT. NAZMTH)                                  IERROR = 10
      IF (NAZ.NE.4 .AND. NAZ.NE.8 
     & .AND. NAZ.NE.12 .AND. NAZ.NE.16 )                    IERROR = 20
      IF (SLOPE.NE.1. .AND. SLOPE.NE.1.5 .AND. SLOPE.NE.2.) IERROR = 30
      IF (IERROR.NE.0)  THEN
        WRITE (6,*) 'aborting in HINES_DSP1: IERROR=',IERROR
        STOP
      END IF 
C
C  Ordering of levels.
C
      IF (LEVBOT.LT.LEVTOP)  THEN
        LEV1 = LEVBOT
        LEV2 = LEVTOP
        IORDER = -1
      ELSE
        LEV1 = LEVTOP
        LEV2 = LEVBOT
        IORDER = 1
      END IF
C
      LEV1P = LEV1 + 1
      LEV2M = LEV2 - 1
C
C  Initialize output and some work arrays.
C
      DO 5 L = 1,NLEVS
      DO 5 I = IL1,IL2
        DRAG_U(I,L) = ZERO
        DRAG_V(I,L) = ZERO
        HEAT(I,L)   = ZERO
        DIFFCO(I,L) = ZERO
        FLUX_U(I,L) = ZERO
        FLUX_V(I,L) = ZERO
        SIGMA_T(I,L) = ZERO
 5    CONTINUE
C
C  Initialize cutoff wavenumber array to minimum value. 
C
      DO 10 N = 1,NAZ
      DO 10 L = 1,NLEVS
      DO 10 I = IL1,IL2
        M_ALPHA(I,L,N) = M_MIN
 10   CONTINUE
C
C  Longitudes where drag to be calculated.
C
      ICOUNT = 0
      DO 15 I = IL1,IL2
        DRAG(I) = .FALSE.
        IF (RMS_WIND(I).GE.RMS_MIN)  THEN
          DRAG(I) = .TRUE.
          ICOUNT = ICOUNT + 1
        END IF
 15   CONTINUE
C
C  Return to calling program if no drag.
C
      IF (ICOUNT.EQ.0)  RETURN
C
C  Buoyancy, density and winds at bottom level.
C
      DO 20 I = IL1,IL2
        IF (DRAG(I))  THEN
          BVFB(I)  = BVFREQ(I,LEVBOT)
          DENSB(I) = DENSITY(I,LEVBOT)
          UBOT(I)  = VEL_U(I,LEVBOT)
          VBOT(I)  = VEL_V(I,LEVBOT)
        END IF
 20   CONTINUE
C
C  Calculate cutoff vertical wavenumber and velocity variances.
C
      CALL HINES_WAVNUM ( M_ALPHA, SIGMA_ALPHA, SIGSQH_ALPHA, SIGMA_T,
     ^                    AK_ALPHA, VEL_U, VEL_V, UBOT, VBOT, 
     ^                    VISC_MOL, DENSITY, DENSB, BVFREQ, BVFB, 
     ^                    RMS_WIND, V_ALPHA, I_ALPHA, MMIN_ALPHA,
     ^                    DO_ALPHA, DRAGIL, DRAG, M_MIN, KSTAR, SLOPE, 
     ^                    F1, F2, F3, NAZ, LEVBOT, LEVTOP, 
     ^                    IL1, IL2,  NLONS, NLEVS, NAZMTH)
C
C  Multiplicative spectral factor at each azimuth.
C
      DO 50 N = 1,NAZ
      DO 50 I = IL1,IL2
        IF (DRAG(I))  THEN
          SPECFAC(I,N) = AK_ALPHA(I,N) * K_ALPHA(I,N)
        END IF
 50   CONTINUE
C
C  Smooth cutoff wavenumbers and total rms velocity in the vertical 
C  direction NSMAX times, using FLUX_U as temporary work array.
C   
      IF (NSMAX.GT.0)  THEN
        DO 80 N = 1,NAZ
          CALL HINES_SMOOTH ( M_ALPHA(1,1,N), 
     ^                        FLUX_U, DRAGIL, SMCO, NSMAX,
     ^                        IL1, IL2, LEV1, LEV2, NLONS, NLEVS )
 80     CONTINUE
        CALL HINES_SMOOTH ( SIGMA_T, 
     ^                      FLUX_U, DRAGIL, SMCO, NSMAX,
     ^                      IL1, IL2, LEV1, LEV2, NLONS, NLEVS )
      END IF
C
C  Calculate zonal and meridional components of the
C  momentum flux and drag.
C
      CALL HINES_DRAG ( FLUX_U, FLUX_V, DRAG_U, DRAG_V, 
     ^                  ALT, DENSITY, DENSB, M_ALPHA, 
     ^                  SPECFAC, DRAGIL, M_MIN, SLOPE, NAZ,
     ^                  IL1, IL2, LEV1, LEV2, NLONS, NLEVS, NAZMTH )
C
C  Heating rate and diffusion coefficient at midpoint between momentum levels.
C
      IF (IHEATCAL.EQ.1)  THEN            
        CALL HINES_HEAT ( HEAT, DIFFCO, 
     ^                    ALT, M_ALPHA, SPECFAC, DRAGIL, 
     ^                    BVFREQ, DENSITY, DENSB, SIGMA_T, VISC_MOL, 
     ^                    KSTAR, SLOPE, F2, F3, F5, F6, NAZ, IL1, IL2,
     ^                    LEVBOT, LEVTOP, NLONS, NLEVS, NAZMTH)
      END IF
C
C  Apply exponential decay to drag, heating and diffusion above 
C  ALT_CUTOFF (use BVFB as temporary work array).
C
      IF (ICUTOFF.EQ.1)  THEN		
        CALL HINES_EXP ( DRAG_U, 
     ^                   BVFB, ALT, ALT_CUTOFF, IORDER,
     ^                   IL1, IL2, LEV1, LEV2, NLONS, NLEVS )
        CALL HINES_EXP ( DRAG_V, 
     ^                   BVFB, ALT, ALT_CUTOFF, IORDER,
     ^                   IL1, IL2, LEV1, LEV2, NLONS, NLEVS )
        IF (IHEATCAL.EQ.1)  THEN		
          CALL HINES_EXP ( HEAT,
     ^                     BVFB, ALT, ALT_CUTOFF, IORDER,
     ^                     IL1, IL2, LEV1, LEV2, NLONS, NLEVS )
          CALL HINES_EXP ( DIFFCO, 
     ^                     BVFB, ALT, ALT_CUTOFF, IORDER,
     ^                     IL1, IL2, LEV1, LEV2, NLONS, NLEVS )
        END IF   
      END IF   
C
C  Print out flux, drag, etc arrays for diagnostic purposes.
C
      IF (IPRNT1.EQ.1)  THEN
        ILPRT1 = 1
        ILPRT2 = 1
        NUNIT  = 11
        CALL HINES_PRNT1 ( FLUX_U, FLUX_V, DRAG_U, DRAG_V, VEL_U, VEL_V, 
     ^                     ALT, SIGMA_T, SIGMA_ALPHA, M_ALPHA,
     ^                     1, 1, NUNIT, ILPRT1, ILPRT2, LEV1, LEV2,
     ^                     NAZ, NLONS, NLEVS, NAZMTH)
      END IF
C
C  Print out azimuthal arrays for diagnostic purposes.
C
      IF (IPRNT2.EQ.1)  THEN
        ILPRT1 = 1
        ILPRT2 = 1
        NUNIT  = 11
        CALL HINES_PRNT2 (M_ALPHA, SIGMA_ALPHA, VEL_U, VEL_V,
     ^                    UBOT, VBOT, ALT, DRAGIL, V_ALPHA,
     ^                    NUNIT, ILPRT1, ILPRT2, LEV1, LEV2,
     ^                    NAZ, NLONS, NLEVS, NAZMTH)
      END IF
C
C  Finished.
C
      RETURN
C-----------------------------------------------------------------------
      END


      SUBROUTINE HINES_WAVNUM (M_ALPHA,SIGMA_ALPHA,SIGSQH_ALPHA,SIGMA_T,
     1                         AK_ALPHA,VEL_U,VEL_V,UBOT,VBOT,VISC_MOL,
     2                         DENSITY,DENSB,BVFREQ,BVFB,
     3                         RMS_WIND,V_ALPHA,I_ALPHA,MMIN_ALPHA,
     4                         DO_ALPHA,DRAGIL,DRAG,M_MIN,KSTAR,SLOPE,
     5                         F1,F2,F3,NAZ,LEVBOT,LEVTOP,IL1,IL2,
     6                         NLONS,NLEVS,NAZMTH)
C
C  This routine calculates the cutoff vertical wavenumber and velocity
C  variances on a longitude by altitude grid needed for the Hines' Doppler 
C  spread gravity wave drag parameterization scheme.
C
C  Aug. 10/95 - C. McLandress
C
C  Modifications 
C  -------------
C  Feb. 2/96 - C. McLandress (added: minimum cutoff wavenumber M_MIN;
C                             12 and 16 azimuths; inclusion of introwaves
C                             by UBOT and VBOT; diffusion and heating
C                             calculated at half level; logical flags;
C                             new print routine) 
C
C  Output arguements:
C  ------------------
C
C     * M_ALPHA      = cutoff wavenumber at each azimuth (1/m).
C     * SIGMA_ALPHA  = total rms wind in each azimuth (m/s).
C     * SIGSQH_ALPHA = portion of wind variance from waves having wave
C     *                normals in the alpha azimuth (m/s).
C     * SIGMA_T      = total rms horizontal wind (m/s).
C     * AK_ALPHA     = spectral amplitude factor at each azimuth 
C     *                (i.e.,{AjKj}) in m^4/s^2.
C     * DRAGIL       = logical flag indicating longitudes and levels where 
C     *                calculations to be performed.
C
C  Input arguements:
C  -----------------
C
C     * VEL_U    = background zonal wind component (m/s).
C     * VEL_V    = background meridional wind component (m/s).
C     * UBOT     = background zonal wind component at bottom level.
C     * VBOT     = background meridional wind component at bottom level.
C     * VISC_MOL = molecular viscosity (m^2/s)
C     * DENSITY  = background density (kg/m^3).
C     * DENSB    = background density at model bottom (kg/m^3).
C     * BVFREQ   = background Brunt Vassala frequency (radians/sec).
C     * BVFB     = background Brunt Vassala frequency at model bottom.
C     * RMS_WIND = root mean square gravity wave wind at lowest level (m/s).
C     * DRAG     = logical flag indicating longitudes where calculations
C     *            to be performed.
C     * KSTAR    = typical gravity wave horizontal wavenumber (1/m).
C     * M_MIN    = minimum allowable cutoff vertical wavenumber, 
C     *            e.g., 1/(3km). This is used only for a spectral slope
C     *            (SLOPE) of one ==> for slope of 1.5 or 2 then M_MIN = 0. 
C     * SLOPE    = slope of incident vertical wavenumber spectrum
C     *            (SLOPE = 1., 1.5 or 2.).
C     * F1,F2,F3 = Hines's fudge factors.
C     * NAZ      = number of horizontal azimuths used (4, 8, 12 or 16).
C     * LEVBOT   = index of lowest vertical level.
C     * LEVTOP   = index of highest vertical level. 
C     * IL1      = first longitudinal index to use (IL1 >= 1).
C     * IL2      = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * NLONS    = number of longitudes.
C     * NLEVS    = number of vertical levels.
C     * NAZMTH   = azimuthal array dimension (NAZMTH >= NAZ).
C
C  Input work arrays:
C  -----------------
C
C     * V_ALPHA    = wind component at each azimuth (m/s). 
C     * I_ALPHA    = Hines' integral at a single level.
C     * MMIN_ALPHA = minimum value of cutoff wavenumber.
C     * DO_ALPHA   = logical flag indicating azimuths and longitudes
C     *              where calculations to be performed.
C
      INTEGER  NAZ, LEVBOT, LEVTOP, IL1, IL2, NLONS, NLEVS, NAZMTH
      LOGICAL DRAG(NLONS), DO_ALPHA(NLONS,NAZMTH), DRAGIL(NLONS,NLEVS)
      REAL  SLOPE, M_MIN, KSTAR, F1, F2, F3
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  SIGMA_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  SIGSQH_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  SIGMA_T(NLONS,NLEVS)
      REAL  AK_ALPHA(NLONS,NAZMTH)
      REAL  VISC_MOL(NLONS,NLEVS)
      REAL  VEL_U(NLONS,NLEVS),    VEL_V(NLONS,NLEVS)
      REAL  UBOT(NLONS),           VBOT(NLONS)
      REAL  DENSITY(NLONS,NLEVS),  DENSB(NLONS)
      REAL  BVFREQ(NLONS,NLEVS),   BVFB(NLONS),  RMS_WIND(NLONS)
      REAL  I_ALPHA(NLONS,NAZMTH), MMIN_ALPHA(NLONS,NAZMTH)
      REAL  V_ALPHA(NLONS,NAZMTH)
C
C Internal variables.
C -------------------
C
      INTEGER  IBIG
      PARAMETER  ( IBIG = 1000 )
      REAL  N_OVER_M(IBIG), SIGFAC(IBIG)
      INTEGER  I, L, N, LSTART, LINCR, LBELOW, ICOUNT
      REAL  M_SUB_M_TURB, M_SUB_M_MOL, M_TRIAL, MMSQ
      REAL  VISC, VISC_MIN, AZFAC, SP1
      DATA  VISC_MIN / 1.E-10 /
C-----------------------------------------------------------------------     
C
C  Quit if work arrays not big enough.
C
      IF (IBIG.LT.NLONS)  THEN
        WRITE (6,*)
        WRITE (6,*) ' HINES_WAVNUM: increase IBIG'
        STOP
      END IF
C
      SP1 = SLOPE + 1.
      MMSQ = M_MIN**2
C
C  Indices of levels to process.
C
      IF (LEVBOT.GT.LEVTOP)  THEN
        LSTART = LEVBOT - 1     
        LINCR  = -1
      ELSE
        LSTART = LEVBOT + 1
        LINCR  = 1
      END IF
C
C  Initialize logical flags and determine number of longitudes
C  at bottom where calculations to be done.
C
      DO 5 L = 1,NLEVS
      DO 5 I = IL1,IL2
        DRAGIL(I,L) = .FALSE.
 5    CONTINUE
      ICOUNT = 0
      DO 10 I = IL1,IL2
        IF (DRAG(I))  ICOUNT = ICOUNT + 1
        DRAGIL(I,LEVBOT) = DRAG(I)
 10   CONTINUE
      DO 15  N = 1,NAZMTH
      DO 15  I = IL1,IL2
        DO_ALPHA(I,N) = DRAG(I)
 15   CONTINUE
C
C  Use horizontal isotropy to calculate azimuthal variances at bottom level.
C
      AZFAC = 1. / FLOAT(NAZ)
      DO 20 N = 1,NAZ
      DO 20 I = IL1,IL2
        IF (DRAG(I))  THEN
          SIGSQH_ALPHA(I,LEVBOT,N) = AZFAC * RMS_WIND(I)**2
        END IF
 20   CONTINUE
C
C  Velocity variances at bottom level.
C
      CALL HINES_SIGMA ( SIGMA_T, SIGMA_ALPHA, 
     ^                   SIGSQH_ALPHA, DRAG, NAZ, LEVBOT, 
     ^                   IL1, IL2, NLONS, NLEVS, NAZMTH)
C
C  Calculate cutoff wavenumber and spectral amplitude factor 
C  at bottom level where it is assumed that background winds vanish
C  and also initialize minimum value of cutoff wavnumber.
C
        IF (SLOPE.EQ.1.)  THEN
          DO 30 N = 1,NAZ
          DO 30 I = IL1,IL2
            IF (DRAG(I))  THEN
              M_ALPHA(I,LEVBOT,N) =  BVFB(I) / 
     ^                           ( F1 * SIGMA_ALPHA(I,LEVBOT,N) 
     ^                           + F2 * SIGMA_T(I,LEVBOT) )
              AK_ALPHA(I,N) = 2. * SIGSQH_ALPHA(I,LEVBOT,N) /
     ^                        ( M_ALPHA(I,LEVBOT,N)**2 - MMSQ ) 
              MMIN_ALPHA(I,N) = M_ALPHA(I,LEVBOT,N)
            END IF
 30       CONTINUE
        ELSE
          DO 40 N = 1,NAZ
          DO 40 I = IL1,IL2
            IF (DRAG(I))  THEN
              M_ALPHA(I,LEVBOT,N) =  BVFB(I) / 
     ^                           ( F1 * SIGMA_ALPHA(I,LEVBOT,N) 
     ^                           + F2 * SIGMA_T(I,LEVBOT) )
              AK_ALPHA(I,N)   = SIGSQH_ALPHA(I,LEVBOT,N) 
     ^                      / ( M_ALPHA(I,LEVBOT,N)**SP1 / SP1 )
              MMIN_ALPHA(I,N) = M_ALPHA(I,LEVBOT,N)
            END IF
 40       CONTINUE
        END IF
C
C  Calculate quantities from the bottom upwards, 
C  starting one level above bottom.
C
      DO 150 L = LSTART,LEVTOP,LINCR
C
C  Return to calling program if no more longitudes left to do.
C
        IF (ICOUNT.EQ.0)  GO TO 160
C
C  Level beneath present level.
C
        LBELOW = L - LINCR 
C
C  Compute azimuthal wind components from zonal and meridional winds.
C
        CALL HINES_WIND ( V_ALPHA, 
     ^                    VEL_U(1,L), VEL_V(1,L), UBOT, VBOT, 
     ^                    DRAGIL(1,LBELOW), NAZ, IL1, IL2, 
     ^                    NLONS, NAZMTH )
C
C  Calculate N/m_M where m_M is maximum permissible value of the vertical
C  wavenumber (i.e., m > m_M are obliterated) and N is buoyancy frequency.
C  m_M is taken as the smaller of the instability-induced 
C  wavenumber (M_SUB_M_TURB) and that imposed by molecular viscosity
C  (M_SUB_M_MOL). Since variance at this level is not yet known
C  use value at level below.
C
        DO 50 I = IL1,IL2
          IF (DRAGIL(I,LBELOW))  THEN
            VISC = AMAX1 ( VISC_MOL(I,L), VISC_MIN )
            M_SUB_M_TURB = BVFREQ(I,L) / ( F2 * SIGMA_T(I,LBELOW) )
            M_SUB_M_MOL  = (BVFREQ(I,L)*KSTAR/VISC)**0.33333333/F3
            IF (M_SUB_M_TURB .LT. M_SUB_M_MOL)  THEN
              N_OVER_M(I) = F2 * SIGMA_T(I,LBELOW)
            ELSE
              N_OVER_M(I) = BVFREQ(I,L) / M_SUB_M_MOL 
            END IF
          END IF
  50    CONTINUE
C
C  Calculate cutoff wavenumber at this level.
C
        DO 60 N = 1,NAZ
        DO 60 I = IL1,IL2
C
          IF (DO_ALPHA(I,N))  THEN 
C
C  Calculate trial value (since variance at this level is not yet known
C  use value at level below). If trial value is negative or if it exceeds 
C  minimum value found at lower levels (not permitted) then set it 
C  to this minimum value. 
C
             denom =  F1 * SIGMA_ALPHA(I,LBELOW,N)  
     ^              + N_OVER_M(I) + V_ALPHA(I,N)
             if(abs(denom).ge.1.e-5) then
                M_TRIAL = BVFB(I) / denom
                IF (M_TRIAL.LE.0. .OR. M_TRIAL.GT.MMIN_ALPHA(I,N))  THEN
                   M_TRIAL = MMIN_ALPHA(I,N)
                END IF
             else
                M_TRIAL = MMIN_ALPHA(I,N)
             end if

            M_ALPHA(I,L,N) = M_TRIAL
C
C  Do not permit cutoff wavenumber to be less than minimum allowable value.
C
            IF (M_ALPHA(I,L,N) .LT. M_MIN)  M_ALPHA(I,L,N) = M_MIN
C
C  Reset minimum value of cutoff wavenumber if necessary.
C
            IF (M_ALPHA(I,L,N) .LT. MMIN_ALPHA(I,N))  THEN
              MMIN_ALPHA(I,N) = M_ALPHA(I,L,N)
            END IF
C
          ELSE
C
            M_ALPHA(I,L,N) = M_MIN
C
          END IF          
C
 60     CONTINUE
C
C  Calculate the Hines integral at this level.
C
        CALL HINES_INTGRL ( I_ALPHA, 
     ^                      V_ALPHA, M_ALPHA, BVFB, 
     ^                      DRAGIL(1,L), DO_ALPHA, M_MIN, SLOPE, NAZ, 
     ^                      L, IL1, IL2, NLONS, NLEVS, NAZMTH )
C
C  Calculate the velocity variances at this level.
C
        DO 80 I = IL1,IL2
          IF (DRAGIL(I,L))  THEN
            SIGFAC(I) = DENSB(I) / DENSITY(I,L) 
     ^                  * BVFREQ(I,L) / BVFB(I) 
          END IF
 80     CONTINUE
C
C  Calculate the velocity variances at this level.
C
        DO 90 N = 1,NAZ
        DO 90 I = IL1,IL2
          IF (DRAGIL(I,L))  THEN
            SIGSQH_ALPHA(I,L,N) = SIGFAC(I) * AK_ALPHA(I,N) 
     ^                            * I_ALPHA(I,N)
          END IF
  90    CONTINUE
        CALL HINES_SIGMA ( SIGMA_T, SIGMA_ALPHA, 
     ^                     SIGSQH_ALPHA, DRAGIL(1,L), NAZ, L, 
     ^                     IL1, IL2, NLONS, NLEVS, NAZMTH )
C
C  If total rms wind zero (no more drag) then set DRAG to false 
C  and update longitude counter.
C
        DO 110 I = IL1,IL2 
          IF (SIGMA_T(I,L) .EQ. 0.)  THEN
            DRAGIL(I,L) = .FALSE. 
            ICOUNT = ICOUNT - 1   
          END IF
 110    CONTINUE
C
C  End of level loop.
C
 150  CONTINUE
C
 160  CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END


      SUBROUTINE HINES_WIND (V_ALPHA,VEL_U,VEL_V,UBOT,VBOT,DRAG,
     1                       NAZ,IL1,IL2,NLONS,NAZMTH)
C
C  This routine calculates the azimuthal horizontal background wind components 
C  on a longitude at a single altitude for the case of 4, 8, 12 or 16 equally
C  spaced azimuths needed for the Hines' Doppler spread GWD parameterization 
C  scheme.
C
C  Aug. 7/95 - C. McLandress
C
C  Modifications:
C  --------------
C  Feb. 2/96 - C. McLandress (added: 12 and 16 azimuths; logical flags;
C                             only single level calculation; removed UMIN)  
C
C  Output arguement:
C  ----------------
C
C     * V_ALPHA = background wind component at each azimuth (m/s). 
C     *           (note: first azimuth is in eastward direction
C     *            and rotate in counterclockwise direction.)
C
C  Input arguements:
C  ----------------
C
C     * VEL_U  = background zonal wind component (m/s).
C     * VEL_V  = background meridional wind component (m/s).
C     * UBOT   = background zonal wind component at bottom level.
C     * VBOT   = background meridional wind component at bottom level.
C     * DRAG   = logical flag indicating longitudes where calculations
C     *          to be performed.
C     * NAZ    = number of horizontal azimuths used (4, 8, 12 or 16).
C     * IL1    = first longitudinal index to use (IL1 >= 1).
C     * IL2    = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * NLONS  = number of longitudes.
C     * NAZMTH = azimuthal array dimension (NAZMTH >= NAZ).
C
C  Constants in DATA statements.
C  ----------------------------
C
C     * COS45 = cosine of 45 degrees. 		
C     * COS30 = cosine of 30 degrees. 		
C     * SIN30 = sine of 30 degrees. 		
C     * COS22 = cosine of 22.5 degrees. 		
C     * SIN22 = sine of 22.5 degrees. 		
C
C  Subroutine arguements.
C  ---------------------
C
      INTEGER  NAZ, IL1, IL2, NLONS, NAZMTH
      LOGICAL  DRAG(NLONS)
      REAL  V_ALPHA(NLONS,NAZMTH)
      REAL  VEL_U(NLONS), VEL_V(NLONS)
      REAL  UBOT(NLONS),  VBOT(NLONS)
C
C  Internal variables.
C  -------------------
C
      INTEGER  I
      REAL U, V, COS45, COS30, SIN30, COS22, SIN22
C
      DATA  COS45 / 0.7071068 /
      DATA  COS30 / 0.8660254 /, SIN30 / 0.5       /
      DATA  COS22 / 0.9238795 /, SIN22 / 0.3826834 / 
C-----------------------------------------------------------------------     
C
C  Case with 4 azimuths.
C
      IF (NAZ.EQ.4)  THEN
        DO 20 I = IL1,IL2
          IF (DRAG(I))  THEN
            U = VEL_U(I) - UBOT(I)
            V = VEL_V(I) - VBOT(I)
            V_ALPHA(I,1) = U 
            V_ALPHA(I,2) = V
            V_ALPHA(I,3) = - U
            V_ALPHA(I,4) = - V
          END IF
 20     CONTINUE
      END IF
C
C  Case with 8 azimuths.
C
      IF (NAZ.EQ.8)  THEN
        DO 30 I = IL1,IL2
          IF (DRAG(I))  THEN
            U = VEL_U(I) - UBOT(I)
            V = VEL_V(I) - VBOT(I)
            V_ALPHA(I,1) = U 
            V_ALPHA(I,2) = COS45 * ( V + U )
            V_ALPHA(I,3) = V
            V_ALPHA(I,4) = COS45 * ( V - U )
            V_ALPHA(I,5) = - U
            V_ALPHA(I,6) = - V_ALPHA(I,2)
            V_ALPHA(I,7) = - V
            V_ALPHA(I,8) = - V_ALPHA(I,4)
          END IF
 30     CONTINUE
      END IF
C
C  Case with 12 azimuths.
C
      IF (NAZ.EQ.12)  THEN
        DO 40 I = IL1,IL2
          IF (DRAG(I))  THEN
            U = VEL_U(I) - UBOT(I)
            V = VEL_V(I) - VBOT(I)
            V_ALPHA(I,1)  = U 
            V_ALPHA(I,2)  = COS30 * U + SIN30 * V
            V_ALPHA(I,3)  = SIN30 * U + COS30 * V
            V_ALPHA(I,4)  = V
            V_ALPHA(I,5)  = - SIN30 * U + COS30 * V
            V_ALPHA(I,6)  = - COS30 * U + SIN30 * V
            V_ALPHA(I,7)  = - U
            V_ALPHA(I,8)  = - V_ALPHA(I,2)
            V_ALPHA(I,9)  = - V_ALPHA(I,3)
            V_ALPHA(I,10) = - V
            V_ALPHA(I,11) = - V_ALPHA(I,5)
            V_ALPHA(I,12) = - V_ALPHA(I,6)
          END IF
 40     CONTINUE
      END IF
C
C  Case with 16 azimuths.
C
      IF (NAZ.EQ.16)  THEN
        DO 50 I = IL1,IL2
          IF (DRAG(I))  THEN
            U = VEL_U(I) - UBOT(I)
            V = VEL_V(I) - VBOT(I)
            V_ALPHA(I,1)  = U 
            V_ALPHA(I,2)  = COS22 * U + SIN22 * V
            V_ALPHA(I,3)  = COS45 * ( U + V )
            V_ALPHA(I,4)  = COS22 * V + SIN22 * U
            V_ALPHA(I,5)  = V
            V_ALPHA(I,6)  = COS22 * V - SIN22 * U
            V_ALPHA(I,7)  = COS45 * ( V - U )
            V_ALPHA(I,8)  = - COS22 * U + SIN22 * V
            V_ALPHA(I,9)  = - U
            V_ALPHA(I,10) = - V_ALPHA(I,2)
            V_ALPHA(I,11) = - V_ALPHA(I,3)
            V_ALPHA(I,12) = - V_ALPHA(I,4)
            V_ALPHA(I,13) = - V
            V_ALPHA(I,14) = - V_ALPHA(I,6)
            V_ALPHA(I,15) = - V_ALPHA(I,7)
            V_ALPHA(I,16) = - V_ALPHA(I,8)
          END IF
 50     CONTINUE
      END IF
C
      RETURN
C-----------------------------------------------------------------------
      END


      SUBROUTINE HINES_DRAG (FLUX_U,FLUX_V,DRAG_U,DRAG_V,ALT,DENSITY,
     1                       DENSB,M_ALPHA,SPECFAC,DRAGIL,M_MIN,SLOPE,
     2                       NAZ,IL1,IL2,LEV1,LEV2,NLONS,NLEVS,NAZMTH)
C
C  Calculate zonal and meridional components of the vertical flux 
C  of horizontal momentum and corresponding wave drag (force per unit mass)
C  on a longitude by altitude grid needed for the Hines Doppler spread 
C  gravity wave parameterization scheme.
C
C  Aug. 6/95 - C. McLandress
C
C  Modifications:
C  --------------
C  Feb. 2/96 - C. McLandress (added: minimum cutoff wavenumber M_MIN;
C                             12 and 16 azimuths; logical flags)
C
C  Output arguements:
C  ------------------
C
C     * FLUX_U  = zonal component of vertical momentum flux (Pascals)
C     * FLUX_V  = meridional component of vertical momentum flux (Pascals)
C     * DRAG_U  = zonal component of drag (m/s^2).
C     * DRAG_V  = meridional component of drag (m/s^2).
C
C  Input arguements:
C  -----------------
C
C     * ALT     = altitudes (m).
C     * DENSITY = background density (kg/m^3).
C     * DENSB   = background density at bottom level (kg/m^3).
C     * M_ALPHA = cutoff vertical wavenumber (1/m).
C     * SPECFAC = AK_ALPHA * K_ALPHA.
C     * DRAGIL  = longitudes and levels where drag to be calculated.
C     * DRAGIL  = logical flag indicating longitudes and levels where 
C     *           calculations to be performed.
C     * M_MIN   = minimum allowable cutoff vertical wavenumber, 
C     *           e.g., 1/(3km). This is used only for a spectral slope
C     *           (SLOPE) of one ==> for slope of 1.5 or 2 then M_MIN = 0. 
C     * SLOPE   = slope of incident vertical wavenumber spectrum.
C     * NAZ     = number of horizontal azimuths used (4, 8, 12 or 16).
C     * IL1     = first longitudinal index to use (IL1 >= 1).
C     * IL2     = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * LEV1    = first altitude level to use (LEV1 >=1). 
C     * LEV2    = last altitude level to use (LEV1 < LEV2 <= NLEVS).
C     * NLONS   = number of longitudes.
C     * NLEVS   = number of vertical levels.
C     * NAZMTH  = azimuthal array dimension (NAZMTH >= NAZ).
C
C  Constants in DATA statement.
C  ----------------------------
C
C     * COS45 = cosine of 45 degrees. 		
C     * COS30 = cosine of 30 degrees. 		
C     * SIN30 = sine of 30 degrees.
C     * COS22 = cosine of 22.5 degrees. 		
C     * SIN22 = sine of 22.5 degrees. 		
C
C  Subroutine arguements.
C  ----------------------
C
      INTEGER  NAZ, IL1, IL2, LEV1, LEV2
      INTEGER  NLONS, NLEVS, NAZMTH
      LOGICAL  DRAGIL(NLONS,NLEVS)
      REAL  M_MIN, SLOPE
      REAL  FLUX_U(NLONS,NLEVS), FLUX_V(NLONS,NLEVS)
      REAL  DRAG_U(NLONS,NLEVS), DRAG_V(NLONS,NLEVS)
      REAL  ALT(NLONS,NLEVS),    DENSITY(NLONS,NLEVS), DENSB(NLONS)
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  SPECFAC(NLONS,NAZMTH)
C
C  Internal variables.
C  -------------------
C
      INTEGER  I, L, LEV1P, LEV2M
      REAL  COS45, COS30, SIN30, COS22, SIN22, DDZ, DDZ2, ZERO
      REAL  FLUX2,  FLUX3,  FLUX4,  FLUX5,  FLUX6,  FLUX7,  FLUX8
      REAL  FLUX9,  FLUX10, FLUX11, FLUX12, FLUX14, FLUX15, FLUX16
      DATA  COS45 / 0.7071068 /
      DATA  COS30 / 0.8660254 /, SIN30 / 0.5       /
      DATA  COS22 / 0.9238795 /, SIN22 / 0.3826834 /
      DATA  ZERO / 0. /   
C-----------------------------------------------------------------------
C
      LEV1P = LEV1 + 1
      LEV2M = LEV2 - 1
C
C  Sum over azimuths for case where SLOPE = 1 (NOTE: this is the only
C  value of the slope for which a nonzero M_MIN is used).
C
      IF (SLOPE.EQ.1.)  THEN
C
C  Case with 4 azimuths.
C
        IF (NAZ.EQ.4)  THEN
          DO 10 L = LEV1,LEV2
          DO 10 I = IL1,IL2
            IF (DRAGIL(I,L))  THEN
              FLUX_U(I,L) = DENSB(I) * (
     ^                      SPECFAC(I,1) * ( M_ALPHA(I,L,1) - M_MIN )
     ^                    - SPECFAC(I,3) * ( M_ALPHA(I,L,3) - M_MIN ) )
              FLUX_V(I,L) = DENSB(I) * (
     ^                      SPECFAC(I,2) * ( M_ALPHA(I,L,2) - M_MIN )
     ^                    - SPECFAC(I,4) * ( M_ALPHA(I,L,4) - M_MIN ) )
            END IF
 10       CONTINUE
        END IF
C
C  Case with 8 azimuths.
C
        IF (NAZ.EQ.8)  THEN
          DO 20 L = LEV1,LEV2
          DO 20 I = IL1,IL2
            IF (DRAGIL(I,L))  THEN
              FLUX2 = SPECFAC(I,2) * ( M_ALPHA(I,L,2) - M_MIN )
              FLUX4 = SPECFAC(I,4) * ( M_ALPHA(I,L,4) - M_MIN )
              FLUX6 = SPECFAC(I,6) * ( M_ALPHA(I,L,6) - M_MIN )
              FLUX8 = SPECFAC(I,8) * ( M_ALPHA(I,L,8) - M_MIN )
              FLUX_U(I,L) = DENSB(I) * (
     ^                      SPECFAC(I,1) * ( M_ALPHA(I,L,1) - M_MIN )
     ^                    - SPECFAC(I,5) * ( M_ALPHA(I,L,5) - M_MIN )
     ^                    + COS45 * ( FLUX2 - FLUX4 - FLUX6 + FLUX8 ) )
              FLUX_V(I,L) = DENSB(I) * (
     ^                      SPECFAC(I,3) * ( M_ALPHA(I,L,3) - M_MIN )
     ^                    - SPECFAC(I,7) * ( M_ALPHA(I,L,7) - M_MIN )
     ^                    + COS45 * ( FLUX2 + FLUX4 - FLUX6 - FLUX8 ) )
            END IF
 20       CONTINUE
        END IF
C
C  Case with 12 azimuths.
C
        IF (NAZ.EQ.12)  THEN
          DO 30 L = LEV1,LEV2
          DO 30 I = IL1,IL2
            IF (DRAGIL(I,L))  THEN
              FLUX2  = SPECFAC(I,2)  * ( M_ALPHA(I,L,2)  - M_MIN )
              FLUX3  = SPECFAC(I,3)  * ( M_ALPHA(I,L,3)  - M_MIN )
              FLUX5  = SPECFAC(I,5)  * ( M_ALPHA(I,L,5)  - M_MIN )
              FLUX6  = SPECFAC(I,6)  * ( M_ALPHA(I,L,6)  - M_MIN )
              FLUX8  = SPECFAC(I,8)  * ( M_ALPHA(I,L,8)  - M_MIN )
              FLUX9  = SPECFAC(I,9)  * ( M_ALPHA(I,L,9)  - M_MIN )
              FLUX11 = SPECFAC(I,11) * ( M_ALPHA(I,L,11) - M_MIN )
              FLUX12 = SPECFAC(I,12) * ( M_ALPHA(I,L,12) - M_MIN )
              FLUX_U(I,L) = DENSB(I) * (
     ^                      SPECFAC(I,1) * ( M_ALPHA(I,L,1) - M_MIN )
     ^                    - SPECFAC(I,7) * ( M_ALPHA(I,L,7) - M_MIN )
     ^                    + COS30 * ( FLUX2 - FLUX6 - FLUX8 + FLUX12 ) 
     ^                    + SIN30 * ( FLUX3 - FLUX5 - FLUX9 + FLUX11 ) )
              FLUX_V(I,L) = DENSB(I) * (
     ^                      SPECFAC(I,4)  * ( M_ALPHA(I,L,4)  - M_MIN )
     ^                    - SPECFAC(I,10) * ( M_ALPHA(I,L,10) - M_MIN )
     ^                    + COS30 * ( FLUX3 + FLUX5 - FLUX9 - FLUX11 ) 
     ^                    + SIN30 * ( FLUX2 + FLUX6 - FLUX8 - FLUX12 ) )
            END IF
 30       CONTINUE
        END IF
C
C  Case with 16 azimuths.
C
        IF (NAZ.EQ.16)  THEN
          DO 40 L = LEV1,LEV2
          DO 40 I = IL1,IL2
            IF (DRAGIL(I,L))  THEN
              FLUX2  = SPECFAC(I,2)  * ( M_ALPHA(I,L,2)  - M_MIN )
              FLUX3  = SPECFAC(I,3)  * ( M_ALPHA(I,L,3)  - M_MIN )
              FLUX4  = SPECFAC(I,4)  * ( M_ALPHA(I,L,4)  - M_MIN )
              FLUX6  = SPECFAC(I,6)  * ( M_ALPHA(I,L,6)  - M_MIN )
              FLUX7  = SPECFAC(I,7)  * ( M_ALPHA(I,L,7)  - M_MIN )
              FLUX8  = SPECFAC(I,8)  * ( M_ALPHA(I,L,8)  - M_MIN )
              FLUX10 = SPECFAC(I,10) * ( M_ALPHA(I,L,10) - M_MIN )
              FLUX11 = SPECFAC(I,11) * ( M_ALPHA(I,L,11) - M_MIN )
              FLUX12 = SPECFAC(I,12) * ( M_ALPHA(I,L,12) - M_MIN )
              FLUX14 = SPECFAC(I,14) * ( M_ALPHA(I,L,14) - M_MIN )
              FLUX15 = SPECFAC(I,15) * ( M_ALPHA(I,L,15) - M_MIN )
              FLUX16 = SPECFAC(I,16) * ( M_ALPHA(I,L,16) - M_MIN )
              FLUX_U(I,L) = DENSB(I) * (
     ^                      SPECFAC(I,1) * ( M_ALPHA(I,L,1) - M_MIN )
     ^                    - SPECFAC(I,9) * ( M_ALPHA(I,L,9) - M_MIN )
     ^                   + COS22 * ( FLUX2 - FLUX8 - FLUX10 + FLUX16 ) 
     ^                   + COS45 * ( FLUX3 - FLUX7 - FLUX11 + FLUX15 ) 
     ^                   + SIN22 * ( FLUX4 - FLUX6 - FLUX12 + FLUX14 ) )
              FLUX_V(I,L) = DENSB(I) * (
     ^                      SPECFAC(I,5)  * ( M_ALPHA(I,L,5)  - M_MIN )
     ^                    - SPECFAC(I,13) * ( M_ALPHA(I,L,13) - M_MIN )
     ^                   + COS22 * ( FLUX4 + FLUX6 - FLUX12 - FLUX14 ) 
     ^                   + COS45 * ( FLUX3 + FLUX7 - FLUX11 - FLUX15 ) 
     ^                   + SIN22 * ( FLUX2 + FLUX8 - FLUX10 - FLUX16 ) )
            END IF
 40       CONTINUE
        END IF
C
      END IF
C
C  Sum over azimuths for case where SLOPE not equal to 1
C  (NOTE: minimum wavenumber cutoff is zero).
C
      IF (SLOPE.NE.1.)  THEN
C
C  Case with 4 azimuths.
C
        IF (NAZ.EQ.4)  THEN
          DO 50 L = LEV1,LEV2
          DO 50 I = IL1,IL2
            IF (DRAGIL(I,L))  THEN
              FLUX_U(I,L) = DENSB(I) / SLOPE * (
     ^                      SPECFAC(I,1) * M_ALPHA(I,L,1)**SLOPE
     ^                    - SPECFAC(I,3) * M_ALPHA(I,L,3)**SLOPE )
              FLUX_V(I,L) = DENSB(I) / SLOPE * (
     ^                      SPECFAC(I,2) * M_ALPHA(I,L,2)**SLOPE
     ^                    - SPECFAC(I,4) * M_ALPHA(I,L,4)**SLOPE )
            END IF
 50       CONTINUE
        END IF
C
C  Case with 8 azimuths.
C
        IF (NAZ.EQ.8)  THEN
          DO 60 L = LEV1,LEV2
          DO 60 I = IL1,IL2
            IF (DRAGIL(I,L))  THEN
              FLUX2 = SPECFAC(I,2) * M_ALPHA(I,L,2)**SLOPE
              FLUX4 = SPECFAC(I,4) * M_ALPHA(I,L,4)**SLOPE
              FLUX6 = SPECFAC(I,6) * M_ALPHA(I,L,6)**SLOPE
              FLUX8 = SPECFAC(I,8) * M_ALPHA(I,L,8)**SLOPE
              FLUX_U(I,L) = DENSB(I) / SLOPE * (
     ^                     SPECFAC(I,1) * M_ALPHA(I,L,1)**SLOPE
     ^                   - SPECFAC(I,5) * M_ALPHA(I,L,5)**SLOPE
     ^                   + COS45 * ( FLUX2 - FLUX4 - FLUX6 + FLUX8 ) )
              FLUX_V(I,L) =  DENSB(I) / SLOPE * (
     ^                     SPECFAC(I,3) * M_ALPHA(I,L,3)**SLOPE
     ^                   - SPECFAC(I,7) * M_ALPHA(I,L,7)**SLOPE
     ^                   + COS45 * ( FLUX2 + FLUX4 - FLUX6 - FLUX8 ) )
            END IF
 60       CONTINUE
        END IF
C
C  Case with 12 azimuths.
C
        IF (NAZ.EQ.12)  THEN
          DO 70 L = LEV1,LEV2
          DO 70 I = IL1,IL2
            IF (DRAGIL(I,L))  THEN
              FLUX2  = SPECFAC(I,2)  * M_ALPHA(I,L,2)**SLOPE
              FLUX3  = SPECFAC(I,3)  * M_ALPHA(I,L,3)**SLOPE
              FLUX5  = SPECFAC(I,5)  * M_ALPHA(I,L,5)**SLOPE
              FLUX6  = SPECFAC(I,6)  * M_ALPHA(I,L,6)**SLOPE
              FLUX8  = SPECFAC(I,8)  * M_ALPHA(I,L,8)**SLOPE
              FLUX9  = SPECFAC(I,9)  * M_ALPHA(I,L,9)**SLOPE
              FLUX11 = SPECFAC(I,11) * M_ALPHA(I,L,11)**SLOPE
              FLUX12 = SPECFAC(I,12) * M_ALPHA(I,L,12)**SLOPE
              FLUX_U(I,L) = DENSB(I) / SLOPE * (
     ^                      SPECFAC(I,1) * M_ALPHA(I,L,1)**SLOPE
     ^                    - SPECFAC(I,7) * M_ALPHA(I,L,7)**SLOPE
     ^                    + COS30 * ( FLUX2 - FLUX6 - FLUX8 + FLUX12 ) 
     ^                    + SIN30 * ( FLUX3 - FLUX5 - FLUX9 + FLUX11 ) )
              FLUX_V(I,L) =  DENSB(I) / SLOPE * (
     ^                      SPECFAC(I,4)  * M_ALPHA(I,L,4)**SLOPE
     ^                    - SPECFAC(I,10) * M_ALPHA(I,L,10)**SLOPE
     ^                    + COS30 * ( FLUX3 + FLUX5 - FLUX9 - FLUX11 ) 
     ^                    + SIN30 * ( FLUX2 + FLUX6 - FLUX8 - FLUX12 ) )
            END IF
 70       CONTINUE
        END IF
C
C  Case with 16 azimuths.
C
        IF (NAZ.EQ.16)  THEN
          DO 80 L = LEV1,LEV2
          DO 80 I = IL1,IL2
            IF (DRAGIL(I,L))  THEN
              FLUX2  = SPECFAC(I,2)  * M_ALPHA(I,L,2)**SLOPE
              FLUX3  = SPECFAC(I,3)  * M_ALPHA(I,L,3)**SLOPE
              FLUX4  = SPECFAC(I,4)  * M_ALPHA(I,L,4)**SLOPE
              FLUX6  = SPECFAC(I,6)  * M_ALPHA(I,L,6)**SLOPE
              FLUX7  = SPECFAC(I,7)  * M_ALPHA(I,L,7)**SLOPE
              FLUX8  = SPECFAC(I,8)  * M_ALPHA(I,L,8)**SLOPE
              FLUX10 = SPECFAC(I,10) * M_ALPHA(I,L,10)**SLOPE
              FLUX11 = SPECFAC(I,11) * M_ALPHA(I,L,11)**SLOPE
              FLUX12 = SPECFAC(I,12) * M_ALPHA(I,L,12)**SLOPE
              FLUX14 = SPECFAC(I,14) * M_ALPHA(I,L,14)**SLOPE
              FLUX15 = SPECFAC(I,15) * M_ALPHA(I,L,15)**SLOPE
              FLUX16 = SPECFAC(I,16) * M_ALPHA(I,L,16)**SLOPE
              FLUX_U(I,L) = DENSB(I) / SLOPE * (
     ^                      SPECFAC(I,1) * M_ALPHA(I,L,1)**SLOPE
     ^                    - SPECFAC(I,9) * M_ALPHA(I,L,9)**SLOPE
     ^                   + COS22 * ( FLUX2 - FLUX8 - FLUX10 + FLUX16 ) 
     ^                   + COS45 * ( FLUX3 - FLUX7 - FLUX11 + FLUX15 ) 
     ^                   + SIN22 * ( FLUX4 - FLUX6 - FLUX12 + FLUX14 ) )
              FLUX_V(I,L) =  DENSB(I) / SLOPE * (
     ^                      SPECFAC(I,5)  * M_ALPHA(I,L,5)**SLOPE
     ^                    - SPECFAC(I,13) * M_ALPHA(I,L,13)**SLOPE
     ^                   + COS22 * ( FLUX4 + FLUX6 - FLUX12 - FLUX14 ) 
     ^                   + COS45 * ( FLUX3 + FLUX7 - FLUX11 - FLUX15 ) 
     ^                   + SIN22 * ( FLUX2 + FLUX8 - FLUX10 - FLUX16 ) )
            END IF
 80       CONTINUE
        END IF
C
      END IF
C
C  Calculate drag at intermediate levels using centered differences.
C      
      DO 90 L = LEV1P,LEV2M
      DO 90 I = IL1,IL2
        IF (DRAGIL(I,L))  THEN
          DDZ2 = DENSITY(I,L) * ( ALT(I,L+1) - ALT(I,L-1) )
          DRAG_U(I,L) = - ( FLUX_U(I,L+1) - FLUX_U(I,L-1) ) / DDZ2
          DRAG_V(I,L) = - ( FLUX_V(I,L+1) - FLUX_V(I,L-1) ) / DDZ2
        END IF
 90   CONTINUE
C
C  Drag at first and last levels using one-side differences.
C 
      DO 100 I = IL1,IL2
        IF (DRAGIL(I,LEV1))  THEN
          DDZ = DENSITY(I,LEV1) * ( ALT(I,LEV1P) - ALT(I,LEV1) ) 
          DRAG_U(I,LEV1) = - ( FLUX_U(I,LEV1P) - FLUX_U(I,LEV1) ) / DDZ
          DRAG_V(I,LEV1) = - ( FLUX_V(I,LEV1P) - FLUX_V(I,LEV1) ) / DDZ
        END IF
 100  CONTINUE
      DO 110 I = IL1,IL2
        IF (DRAGIL(I,LEV2))  THEN
          DDZ = DENSITY(I,LEV2) * ( ALT(I,LEV2) - ALT(I,LEV2M) )
          DRAG_U(I,LEV2) = - ( FLUX_U(I,LEV2) - FLUX_U(I,LEV2M) ) / DDZ
          DRAG_V(I,LEV2) = - ( FLUX_V(I,LEV2) - FLUX_V(I,LEV2M) ) / DDZ
        END IF
 110  CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END


      SUBROUTINE HINES_HEAT (HEAT,DIFFCO,ALT,M_ALPHA,SPECFAC,
     1                       DRAGIL,BVFREQ,DENSITY,DENSB,SIGMA_T,
     2                       VISC_MOL,KSTAR,SLOPE,F2,F3,F5,F6,NAZ,
     3                       IL1,IL2,LEVBOT,LEVTOP,NLONS,NLEVS,NAZMTH)
C
C  This routine calculates the gravity wave induced heating and eddy 
C  diffusion coefficient on a longitude by altitude grid for the Hines 
C  Doppler spread parameterization scheme. The output is placed on the 
C  intermediate levels such that the highest (and lowest) levels for 
C  diffusion and heating are 1/2 grid step above (and below) the highest 
C  (and lowest) momentum levels. This routine can be used for nonzero 
C  minimum cutoff wavenumber (M_MIN) only in the case of spectral SLOPE = 1,
C  in which case M_MIN is not needed since its vertical derivative is zero.
C
C  Aug. 6/95 - C. McLandress
C
C  Modifications:
C  --------------
C  Feb. 2/96 - C. McLandress (diffusion and heating calculated at half levels;
C                             logical flags added; vertical derivative of
C                             cutoff wavenumber calculated in this routine;
C                             molecular viscosity not used in calculation
C                             of eddy diffusion coefficient) 
C  Jul. 1/96 - C. McLandress (R1**R2 for R1<0 avoided in diffusion calc.)
C
C  Output arguements:
C  ------------------
C     * HEAT   = gravity wave heating (K/sec) (must be initialized to zero).
C     * DIFFCO = diffusion coefficient (m^2/sec).
C
C  Input arguements:
C  -----------------
C     * ALT         = altitude (m) of levels.
C     * M_ALPHA     = cutoff vertical wavenumber (1/m) specified on
C     *               the momentum levels.
C     * SPECFAC     = AK_ALPHA * K_ALPHA.
C     * DRAGIL      = longitudes and levels where drag to be calculated.
C     * BVFREQ      = background Brunt Vassala frequency (rad/sec).
C     * DENSITY     = background density (kg/m^3).
C     * DENSB       = background density at bottom level (kg/m^3).
C     * SIGMA_T     = total rms horizontal wind (m/s).
C     * VISC_MOL    = molecular viscosity (m^2/s).
C     * KSTAR       = typical gravity wave horizontal wavenumber (1/m).
C     * SLOPE       = slope of incident vertical wavenumber spectrum.
C     * F2,F3,F5,F6 = Hines's fudge factors.
C     * NAZ         = number of horizontal azimuths used.
C     * IL1         = first longitudinal index to use (IL1 >= 1).
C     * IL2         = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * LEVBOT      = index of bottom drag level.
C     * LEVTOP      = index of top drag level.
C     * NLONS       = number of longitudes.
C     * NLEVS       = number of vertical levels.
C     * NAZMTH      = azimuthal array dimension (NAZMTH >= NAZ).
C
      INTEGER  NAZ, IL1, IL2, LEVBOT, LEVTOP, NLONS, NLEVS, NAZMTH
      LOGICAL  DRAGIL(NLONS,NLEVS)
      REAL  KSTAR, SLOPE, F2, F3, F5, F6
      REAL  HEAT(NLONS,NLEVS), DIFFCO(NLONS,NLEVS), ALT(NLONS,NLEVS)
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH), SPECFAC(NLONS,NAZMTH)
      REAL  BVFREQ(NLONS,NLEVS), DENSITY(NLONS,NLEVS),  DENSB(NLONS) 
      REAL  SIGMA_T(NLONS,NLEVS), VISC_MOL(NLONS,NLEVS)
C
C Internal variables.
C -------------------
C
      INTEGER  I, L, N, L1, LH, LP1
      INTEGER  LEV1, LEV2, LEV1P, LEV2M
      REAL  M_SUB_M_TURB, M_SUB_M_MOL, M_SUB_M, HEATNG
      REAL  VISC, VISC_MIN, CPGAS, SM1, DMDZ, MALP, SIGMA, BVF, DENS
C
      DATA  CPGAS / 1004. /         ! specific heat at constant pressure
      DATA  VISC_MIN / 1.E-10 /     ! minimum permissible viscosity
C-----------------------------------------------------------------------     
C
C  For MAM level indices which increase from top down then levels of
C  diffusion coefficient and heating defined so that half level L=2 is 
C  half way between fulle levels L=1 (the model top)  and L=2, etc. 
C  For other models in which the indices increase from bottom up then 
C  first intermediate half level designated L=1.
C
      L1 = 0
      LEV1 = LEVBOT
      LEV2 = LEVTOP - 1
      IF (LEVBOT.GT.LEVTOP)  THEN
        L1 = 1
        LEV1 = LEVTOP
        LEV2 = LEVBOT - 1
      END IF
C 
C  Calculate sum of azimuths for case where SLOPE = 1.
C
      IF (SLOPE.EQ.1.)  THEN
        DO 10 N = 1,NAZ
        DO 10 L = LEV1,LEV2
          DO 5 I = IL1,IL2
            IF (DRAGIL(I,L))  THEN
              HEAT(I,L+L1) = HEAT(I,L+L1) + SPECFAC(I,N) 
     ^                       * ( M_ALPHA(I,L+1,N) - M_ALPHA(I,L,N) ) 
     ^                         / ( ALT(I,L+1) - ALT(I,L) )
            END IF
   5      CONTINUE
  10    CONTINUE
      END IF
C
C  For SLOPE not 1.
C
      IF (SLOPE.NE.1.)  THEN
        SM1 = SLOPE - 1.
        DO 20 N = 1,NAZ
        DO 20 L = LEV1,LEV2
          LH  = L + L1
          LP1 = L + 1           
          DO 15 I = IL1,IL2
            IF (DRAGIL(I,L))  THEN
              DMDZ = ( M_ALPHA(I,LP1,N) - M_ALPHA(I,L,N) ) 
     ^               / ( ALT(I,LP1) - ALT(I,L) )
              MALP = ( M_ALPHA(I,LP1,N) + M_ALPHA(I,L,N) ) / 2. 
              HEAT(I,LH) = HEAT(I,LH) + SPECFAC(I,N) * MALP**SM1 * DMDZ 
            END IF
 15       CONTINUE
 20     CONTINUE
      END IF
C
C  Avoid round off problems that may result in A**B where A < 0 by ensuring
C  that array HEAT does not change sign (which must be since dm/dz <= 0).
C  Note that sign of summed term HEAT is forced to be positive.
C
      DO 22 L = LEV1,LEV2
      DO 22 I = IL1,IL2
        HEAT(I,L+L1) = AMIN1 ( HEAT(I,L+L1), 0.0 )
        HEAT(I,L+L1) = ABS (HEAT(I,L+L1))
 22   CONTINUE 
C
C  Heating and diffusion.
C
      DO 30 L = LEV1,LEV2
C
        DO 25 I = IL1,IL2
C
          IF (DRAGIL(I,L))  THEN
C
C  Interpolate quantities at half levels.
C
            BVF   = ( BVFREQ(I,L)   + BVFREQ(I,L+1)   ) / 2.
            DENS  = ( DENSITY(I,L)  + DENSITY(I,L+1)  ) / 2.
            SIGMA = ( SIGMA_T(I,L)  + SIGMA_T(I,L+1)  ) / 2.
            VISC  = ( VISC_MOL(I,L) + VISC_MOL(I,L+1) ) / 2.
C
C  For heating allow maximum permissible value of cutoff wavenumber to
C  be the smaller of the instability-induced wavenumber (M_SUB_M_TURB) 
C  and that imposed by molecular viscosity (M_SUB_M_MOL), while for
C  diffusion coefficient ignore molecular viscosity contribution
C  since want the turbulent generated diffusion.
C
            VISC = AMAX1 ( VISC, VISC_MIN )
            M_SUB_M_TURB = BVF / ( F2 * SIGMA )
            M_SUB_M_MOL  = ( BVF * KSTAR / VISC )**0.33333333 / F3
C
C  Turbulent diffusion coefficient.
C
            HEATNG = HEAT(I,L+L1) 
     ^               * F5 * BVF / M_SUB_M_TURB * DENSB(I) / DENS
            DIFFCO(I,L+L1) = F6 * HEATNG**0.33333333 
     ^                       / M_SUB_M_TURB**1.33333333
C
C  Heating rate.
C
            M_SUB_M = AMIN1 ( M_SUB_M_TURB, M_SUB_M_MOL )
            HEATNG  = HEAT(I,L+L1) 
     ^                * F5 * BVF / M_SUB_M * DENSB(I) / DENS
            HEAT(I,L+L1) = HEATNG / CPGAS
C
          END IF
C
 25     CONTINUE
C
 30   CONTINUE
C
C  Calculate diffusion and heating at top level for the case
C  where indices of vertical levels increase from top down.
C
      IF (LEVBOT.GT.LEVTOP)  THEN
        LEV1P = LEVTOP + 1
        DO 40 I = IL1,IL2
          IF (DRAGIL(I,LEV1P))  THEN
            DIFFCO(I,LEVTOP) = DIFFCO(I,LEV1P) 
            HEAT(I,LEVTOP)   = HEAT(I,LEV1P) 
          END IF
 40     CONTINUE
        RETURN
      END IF
C
C  Calculate diffusion and heating at top level for the case
C  where indices of vertical levels increase from bottom up.
C
      IF (LEVBOT.LT.LEVTOP)  THEN
        LEV2M = LEVTOP - 1
        DO 50 I = IL1,IL2
          IF (DRAGIL(I,LEV2M))  THEN
            DIFFCO(I,LEVTOP) = DIFFCO(I,LEV2M) 
            HEAT(I,LEVTOP)   = HEAT(I,LEV2M) 
          END IF
 50     CONTINUE
      END IF
C
      RETURN
C-----------------------------------------------------------------------
      END


      SUBROUTINE HINES_SIGMA (SIGMA_T,SIGMA_ALPHA,SIGSQH_ALPHA,DRAG,
     1                        NAZ,LEV,IL1,IL2,NLONS,NLEVS,NAZMTH)
C
C  This routine calculates the total rms and azimuthal rms horizontal 
C  velocities at a given level on a longitude by altitude grid for 
C  the Hines' Doppler spread GWD parameterization scheme.
C  NOTE: only 4, 8, 12 or 16 azimuths can be used.
C
C  Aug. 7/95 - C. McLandress
C
C  Modifications:
C  --------------
C  Feb. 2/96 - C. McLandress (added: 12 and 16 azimuths; logical flags) 
C
C  Output arguements:
C  ------------------
C
C     * SIGMA_T     = total rms horizontal wind (m/s).
C     * SIGMA_ALPHA = total rms wind in each azimuth (m/s).
C
C  Input arguements:
C  -----------------
C
C     * SIGSQH_ALPHA = portion of wind variance from waves having wave
C     *                normals in the alpha azimuth (m/s).
C     * DRAG      = logical flag indicating longitudes where calculations
C     *             to be performed.
C     * NAZ       = number of horizontal azimuths used (4, 8, 12 or 16).
C     * LEV       = altitude level to process.
C     * IL1       = first longitudinal index to use (IL1 >= 1).
C     * IL2       = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * NLONS     = number of longitudes.
C     * NLEVS     = number of vertical levels.
C     * NAZMTH    = azimuthal array dimension (NAZMTH >= NAZ).
C
C  Subroutine arguements.
C  ---------------------
C
      INTEGER  LEV, NAZ, IL1, IL2
      INTEGER  NLONS, NLEVS, NAZMTH
      LOGICAL  DRAG(NLONS)
      REAL  SIGMA_T(NLONS,NLEVS)
      REAL  SIGMA_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  SIGSQH_ALPHA(NLONS,NLEVS,NAZMTH)
C
C  Internal variables.
C  -------------------
C
      INTEGER  I, N
      REAL  SUM_2468,   SUM_1357 
      REAL  SUM_26812,  SUM_35911,  SUM_1379
      REAL  SUM_461012, SUM_24810,  SUM_15711  
      REAL  SUM_281016, SUM_371115, SUM_461214, SUM_13911  
      REAL  SUM_481216, SUM_571315, SUM_241012, SUM_15913
      REAL  SUM_681416, SUM_351113, SUM_261014, SUM_17915
      REAL  C22SQ, C67SQ
      DATA  C22SQ / 0.8535534 /, C67SQ / 0.1464466 /  
C-----------------------------------------------------------------------     
C
C  Calculate azimuthal rms velocity for the 4 azimuth case.
C
      IF (NAZ.EQ.4)  THEN
        DO 10 I = IL1,IL2
          IF (DRAG(I))  THEN
            SIGMA_ALPHA(I,LEV,1) = SQRT ( SIGSQH_ALPHA(I,LEV,1)
     ^                                  + SIGSQH_ALPHA(I,LEV,3) )
            SIGMA_ALPHA(I,LEV,2) = SQRT ( SIGSQH_ALPHA(I,LEV,2)
     ^                                  + SIGSQH_ALPHA(I,LEV,4) )
            SIGMA_ALPHA(I,LEV,3) = SIGMA_ALPHA(I,LEV,1)
            SIGMA_ALPHA(I,LEV,4) = SIGMA_ALPHA(I,LEV,2)
          END IF
 10     CONTINUE
      END IF
C
C  Calculate azimuthal rms velocity for the 8 azimuth case.
C
      IF (NAZ.EQ.8)  THEN
        DO 20 I = IL1,IL2
          IF (DRAG(I))  THEN
            SUM_1357 = SIGSQH_ALPHA(I,LEV,1) + SIGSQH_ALPHA(I,LEV,3) 
     ^               + SIGSQH_ALPHA(I,LEV,5) + SIGSQH_ALPHA(I,LEV,7) 
            SUM_2468 = SIGSQH_ALPHA(I,LEV,2) + SIGSQH_ALPHA(I,LEV,4)
     ^               + SIGSQH_ALPHA(I,LEV,6) + SIGSQH_ALPHA(I,LEV,8)
            SIGMA_ALPHA(I,LEV,1) = SQRT ( SIGSQH_ALPHA(I,LEV,1) 
     ^                           + SIGSQH_ALPHA(I,LEV,5) + SUM_2468/2. )
            SIGMA_ALPHA(I,LEV,2) = SQRT ( SIGSQH_ALPHA(I,LEV,2) 
     ^                           + SIGSQH_ALPHA(I,LEV,6) + SUM_1357/2. )
            SIGMA_ALPHA(I,LEV,3) = SQRT ( SIGSQH_ALPHA(I,LEV,3) 
     ^                           + SIGSQH_ALPHA(I,LEV,7) + SUM_2468/2. )
            SIGMA_ALPHA(I,LEV,4) = SQRT ( SIGSQH_ALPHA(I,LEV,4) 
     ^                           + SIGSQH_ALPHA(I,LEV,8) + SUM_1357/2. )
            SIGMA_ALPHA(I,LEV,5) = SIGMA_ALPHA(I,LEV,1)
            SIGMA_ALPHA(I,LEV,6) = SIGMA_ALPHA(I,LEV,2)
            SIGMA_ALPHA(I,LEV,7) = SIGMA_ALPHA(I,LEV,3)
            SIGMA_ALPHA(I,LEV,8) = SIGMA_ALPHA(I,LEV,4)
          END IF
 20     CONTINUE
      END IF
C
C  Calculate azimuthal rms velocity for the 12 azimuth case.
C
      IF (NAZ.EQ.12)  THEN
        DO 30 I = IL1,IL2
          IF (DRAG(I))  THEN
            SUM_26812  = SIGSQH_ALPHA(I,LEV, 2) + SIGSQH_ALPHA(I,LEV, 6) 
     ^                 + SIGSQH_ALPHA(I,LEV, 8) + SIGSQH_ALPHA(I,LEV,12) 
            SUM_35911  = SIGSQH_ALPHA(I,LEV, 3) + SIGSQH_ALPHA(I,LEV, 5) 
     ^                 + SIGSQH_ALPHA(I,LEV, 9) + SIGSQH_ALPHA(I,LEV,11) 
            SUM_1379   = SIGSQH_ALPHA(I,LEV, 1) + SIGSQH_ALPHA(I,LEV, 3) 
     ^                 + SIGSQH_ALPHA(I,LEV, 7) + SIGSQH_ALPHA(I,LEV, 9) 
            SUM_461012 = SIGSQH_ALPHA(I,LEV, 4) + SIGSQH_ALPHA(I,LEV, 6) 
     ^                 + SIGSQH_ALPHA(I,LEV,10) + SIGSQH_ALPHA(I,LEV,12) 
            SUM_24810  = SIGSQH_ALPHA(I,LEV, 2) + SIGSQH_ALPHA(I,LEV, 4) 
     ^                 + SIGSQH_ALPHA(I,LEV, 8) + SIGSQH_ALPHA(I,LEV,10) 
            SUM_15711  = SIGSQH_ALPHA(I,LEV, 1) + SIGSQH_ALPHA(I,LEV, 5) 
     ^                 + SIGSQH_ALPHA(I,LEV, 7) + SIGSQH_ALPHA(I,LEV,11)
            SIGMA_ALPHA(I,LEV,1)  = SQRT ( SIGSQH_ALPHA(I,LEV,1) 
     ^                            + SIGSQH_ALPHA(I,LEV,7) 
     ^                            + 0.75*SUM_26812 + 0.25*SUM_35911 )
            SIGMA_ALPHA(I,LEV,2)  = SQRT ( SIGSQH_ALPHA(I,LEV,2) 
     ^                            + SIGSQH_ALPHA(I,LEV,8) 
     ^                            + 0.75*SUM_1379 + 0.25*SUM_461012 )
            SIGMA_ALPHA(I,LEV,3)  = SQRT ( SIGSQH_ALPHA(I,LEV,3) 
     ^                            + SIGSQH_ALPHA(I,LEV,9) 
     ^                            + 0.75*SUM_24810 + 0.25*SUM_15711 )
            SIGMA_ALPHA(I,LEV,4)  = SQRT ( SIGSQH_ALPHA(I,LEV,4) 
     ^                            + SIGSQH_ALPHA(I,LEV,10) 
     ^                            + 0.75*SUM_35911 + 0.25*SUM_26812 )
            SIGMA_ALPHA(I,LEV,5)  = SQRT ( SIGSQH_ALPHA(I,LEV,5) 
     ^                            + SIGSQH_ALPHA(I,LEV,11) 
     ^                            + 0.75*SUM_461012 + 0.25*SUM_1379 )
            SIGMA_ALPHA(I,LEV,6)  = SQRT ( SIGSQH_ALPHA(I,LEV,6) 
     ^                            + SIGSQH_ALPHA(I,LEV,12) 
     ^                            + 0.75*SUM_15711 + 0.25*SUM_24810 )
            SIGMA_ALPHA(I,LEV,7)  = SIGMA_ALPHA(I,LEV,1)
            SIGMA_ALPHA(I,LEV,8)  = SIGMA_ALPHA(I,LEV,2)
            SIGMA_ALPHA(I,LEV,9)  = SIGMA_ALPHA(I,LEV,3)
            SIGMA_ALPHA(I,LEV,10) = SIGMA_ALPHA(I,LEV,4)
            SIGMA_ALPHA(I,LEV,11) = SIGMA_ALPHA(I,LEV,5)
            SIGMA_ALPHA(I,LEV,12) = SIGMA_ALPHA(I,LEV,6)
          END IF
 30    CONTINUE
      END IF
C
C  Calculate azimuthal rms velocity for the 16 azimuth case.
C
      IF (NAZ.EQ.16)  THEN
        DO 40 I = IL1,IL2
          IF (DRAG(I))  THEN
            SUM_281016 = SIGSQH_ALPHA(I,LEV, 2) + SIGSQH_ALPHA(I,LEV, 8) 
     ^                 + SIGSQH_ALPHA(I,LEV,10) + SIGSQH_ALPHA(I,LEV,16) 
            SUM_371115 = SIGSQH_ALPHA(I,LEV, 3) + SIGSQH_ALPHA(I,LEV, 7) 
     ^                 + SIGSQH_ALPHA(I,LEV,11) + SIGSQH_ALPHA(I,LEV,15) 
            SUM_461214 = SIGSQH_ALPHA(I,LEV, 4) + SIGSQH_ALPHA(I,LEV, 6) 
     ^                 + SIGSQH_ALPHA(I,LEV,12) + SIGSQH_ALPHA(I,LEV,14) 
            SUM_13911  = SIGSQH_ALPHA(I,LEV, 1) + SIGSQH_ALPHA(I,LEV, 3) 
     ^                 + SIGSQH_ALPHA(I,LEV, 9) + SIGSQH_ALPHA(I,LEV,11) 
            SUM_481216 = SIGSQH_ALPHA(I,LEV, 4) + SIGSQH_ALPHA(I,LEV, 8) 
     ^                 + SIGSQH_ALPHA(I,LEV,12) + SIGSQH_ALPHA(I,LEV,16) 
            SUM_571315 = SIGSQH_ALPHA(I,LEV, 5) + SIGSQH_ALPHA(I,LEV, 7) 
     ^                 + SIGSQH_ALPHA(I,LEV,13) + SIGSQH_ALPHA(I,LEV,15)
            SUM_241012 = SIGSQH_ALPHA(I,LEV, 2) + SIGSQH_ALPHA(I,LEV, 4) 
     ^                 + SIGSQH_ALPHA(I,LEV,10) + SIGSQH_ALPHA(I,LEV,12) 
            SUM_15913  = SIGSQH_ALPHA(I,LEV, 1) + SIGSQH_ALPHA(I,LEV, 5) 
     ^                 + SIGSQH_ALPHA(I,LEV, 9) + SIGSQH_ALPHA(I,LEV,13) 
            SUM_681416 = SIGSQH_ALPHA(I,LEV, 6) + SIGSQH_ALPHA(I,LEV, 8) 
     ^                 + SIGSQH_ALPHA(I,LEV,14) + SIGSQH_ALPHA(I,LEV,16) 
            SUM_351113 = SIGSQH_ALPHA(I,LEV, 3) + SIGSQH_ALPHA(I,LEV, 5) 
     ^                 + SIGSQH_ALPHA(I,LEV,11) + SIGSQH_ALPHA(I,LEV,13) 
            SUM_261014 = SIGSQH_ALPHA(I,LEV, 2) + SIGSQH_ALPHA(I,LEV, 6) 
     ^                 + SIGSQH_ALPHA(I,LEV,10) + SIGSQH_ALPHA(I,LEV,14) 
            SUM_17915  = SIGSQH_ALPHA(I,LEV, 1) + SIGSQH_ALPHA(I,LEV, 7) 
     ^                 + SIGSQH_ALPHA(I,LEV, 9) + SIGSQH_ALPHA(I,LEV,15)
            SIGMA_ALPHA(I,LEV,1)  = SQRT (  
     ^                   SIGSQH_ALPHA(I,LEV, 1) + SIGSQH_ALPHA(I,LEV, 9) 
     ^                 + C22SQ * SUM_281016     + 0.5 * SUM_371115 
     ^                 + C67SQ * SUM_461214 ) 
            SIGMA_ALPHA(I,LEV,2)  = SQRT (  
     ^                   SIGSQH_ALPHA(I,LEV, 2) + SIGSQH_ALPHA(I,LEV,10) 
     ^                 + C22SQ * SUM_13911      + 0.5 * SUM_481216 
     ^                 + C67SQ * SUM_571315 ) 
            SIGMA_ALPHA(I,LEV,3)  = SQRT (  
     ^                   SIGSQH_ALPHA(I,LEV, 3) + SIGSQH_ALPHA(I,LEV,11) 
     ^                 + C22SQ * SUM_241012     + 0.5 * SUM_15913
     ^                 + C67SQ * SUM_681416 ) 
            SIGMA_ALPHA(I,LEV,4)  = SQRT (  
     ^                   SIGSQH_ALPHA(I,LEV, 4) + SIGSQH_ALPHA(I,LEV,12) 
     ^                 + C22SQ * SUM_351113     + 0.5 * SUM_261014
     ^                 + C67SQ * SUM_17915 ) 
            SIGMA_ALPHA(I,LEV,5)  = SQRT (  
     ^                   SIGSQH_ALPHA(I,LEV, 5) + SIGSQH_ALPHA(I,LEV,13) 
     ^                 + C22SQ * SUM_461214     + 0.5 * SUM_371115
     ^                 + C67SQ * SUM_281016 ) 
            SIGMA_ALPHA(I,LEV,6)  = SQRT (  
     ^                   SIGSQH_ALPHA(I,LEV, 6) + SIGSQH_ALPHA(I,LEV,14) 
     ^                 + C22SQ * SUM_571315     + 0.5 * SUM_481216
     ^                 + C67SQ * SUM_13911 ) 
            SIGMA_ALPHA(I,LEV,7)  = SQRT (  
     ^                   SIGSQH_ALPHA(I,LEV, 7) + SIGSQH_ALPHA(I,LEV,15) 
     ^                 + C22SQ * SUM_681416     + 0.5 * SUM_15913
     ^                 + C67SQ * SUM_241012 ) 
            SIGMA_ALPHA(I,LEV,8)  = SQRT (  
     ^                   SIGSQH_ALPHA(I,LEV, 8) + SIGSQH_ALPHA(I,LEV,16) 
     ^                 + C22SQ * SUM_17915      + 0.5 * SUM_261014
     ^                 + C67SQ * SUM_351113 ) 
            SIGMA_ALPHA(I,LEV,9)  = SIGMA_ALPHA(I,LEV,1)
            SIGMA_ALPHA(I,LEV,10) = SIGMA_ALPHA(I,LEV,2)
            SIGMA_ALPHA(I,LEV,11) = SIGMA_ALPHA(I,LEV,3)
            SIGMA_ALPHA(I,LEV,12) = SIGMA_ALPHA(I,LEV,4)
            SIGMA_ALPHA(I,LEV,13) = SIGMA_ALPHA(I,LEV,5)
            SIGMA_ALPHA(I,LEV,14) = SIGMA_ALPHA(I,LEV,6)
            SIGMA_ALPHA(I,LEV,15) = SIGMA_ALPHA(I,LEV,7)
            SIGMA_ALPHA(I,LEV,16) = SIGMA_ALPHA(I,LEV,8)
          END IF
 40    CONTINUE
      END IF
C
C  Initialize rms wind.
C
      DO 50 I = IL1,IL2
        SIGMA_T(I,LEV) = 0.
 50   CONTINUE
C
C  Calculate total rms wind.
C
      DO 60 N = 1,NAZ
      DO 60 I = IL1,IL2
        IF (DRAG(I))  THEN
          SIGMA_T(I,LEV) = SIGMA_T(I,LEV) + SIGSQH_ALPHA(I,LEV,N)
        END IF   
 60   CONTINUE
      DO 70 I = IL1,IL2
        IF (DRAG(I))  THEN
          SIGMA_T(I,LEV) = SQRT ( SIGMA_T(I,LEV) )
        END IF   
 70   CONTINUE
C
      RETURN
C-----------------------------------------------------------------------     
      END


      SUBROUTINE HINES_INTGRL (I_ALPHA,V_ALPHA,M_ALPHA,BVFB,
     1                         DRAG,DO_ALPHA,M_MIN,SLOPE,NAZ,LEV,
     2                         IL1,IL2,NLONS,NLEVS,NAZMTH)
C
C  This routine calculates the vertical wavenumber integral
C  for a single vertical level at each azimuth on a longitude grid
C  for the Hines' Doppler spread GWD parameterization scheme.
C  NOTE: (1) only spectral slopes of 1, 1.5 or 2 are permitted.
C        (2) the integral is written in terms of the product QM
C            which by construction is always less than 1. Series
C            solutions are used for small |QM| and analytical solutions
C            for remaining values.
C
C  Aug. 8/95 - C. McLandress
C
C  Modifications:
C  --------------
C  Feb. 2/96 - C. McLandress (added minimum cutoff wavenumber M_MIN and
C                             logical flags; I_ALPHA < 0 disallowed)
C
C  Output arguements:
C  ------------------
C
C     * I_ALPHA  = Hines' integral.
C     * DRAG     = logical flag indicating longitudes where calculations
C     *            to be performed.
C     * DO_ALPHA = logical flag indicating azimuths and longitudes
C     *            where calculations to be performed (for SLOPE=1).
C
C  Input arguements:
C  -----------------
C
C     * V_ALPHA  = azimuthal wind component (m/s) at this level.
C     * M_ALPHA  = azimuthal cutoff vertical wavenumber (1/m).
C     * BVFB     = background Brunt Vassala frequency at model bottom.
C     * DO_ALPHA = logical flag indicating azimuths and longitudes
C     *            where calculations to be performed.
C     * M_MIN    = minimum allowable cutoff vertical wavenumber, 
C     *            e.g., 1/(3km). This is used only for a spectral slope
C     *            (SLOPE) of one ==> for slope of 1.5 or 2 then M_MIN = 0. 
C     * SLOPE    = slope of initial vertical wavenumber spectrum 
C     *            (must use SLOPE = 1, 1.5 or 2)
C     * NAZ      = number of horizontal azimuths used.
C     * LEV      = altitude level to process.
C     * IL1      = first longitudinal index to use (IL1 >= 1).
C     * IL2      = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * NLONS    = number of longitudes.
C     * NLEVS    = number of vertical levels.
C     * NAZMTH   = azimuthal array dimension (NAZMTH >= NAZ).
C
C  Constants in DATA statements:
C  ----------------------------
C
C     * QMIN = minimum value of Q_ALPHA (avoids indeterminant form of integral)
C     * QM_MIN = minimum value of Q_ALPHA * M_ALPHA (used to avoid numerical
C     *          problems).
C
      INTEGER  LEV, NAZ, IL1, IL2, NLONS, NLEVS, NAZMTH
      LOGICAL  DRAG(NLONS), DO_ALPHA(NLONS,NAZMTH)
      REAL  I_ALPHA(NLONS,NAZMTH), V_ALPHA(NLONS,NAZMTH)
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  BVFB(NLONS), M_MIN, SLOPE
C
C  Internal variables.
C  -------------------
C
      INTEGER  I, N
      REAL  Q_ALPHA, QM, QMM, SQRTQM, Q_MIN, QM_MIN, IVAL, ZERO
C
      DATA  Q_MIN / 1.0 /, QM_MIN / 0.01 /, ZERO / 0. /
C-----------------------------------------------------------------------     
C
C  Initialize DRAG array.
C
      DO 5 I = IL1,IL2
        DRAG(I) = .FALSE.
 5    CONTINUE
C
C  For integer value SLOPE = 1.
C
      IF (SLOPE .EQ. 1.)  THEN
C
        DO 10 N = 1,NAZ
        DO 10 I = IL1,IL2
C
C  Calculate integral only in regions where cutoff wavenumber
C  is greater than minimum allowable value, otherwise set integral
C  to zero and turn off flag to calculate this azimuth.
C
          IF (M_ALPHA(I,LEV,N).GT.M_MIN)  THEN
C
            DRAG(I) = .TRUE.
            Q_ALPHA = V_ALPHA(I,N) / BVFB(I)
            QM      = Q_ALPHA * M_ALPHA(I,LEV,N)
            QMM     = Q_ALPHA * M_MIN 
C
C  If |QM| is small then use first 4 terms series of Taylor series
C  expansion of integral in order to avoid indeterminate form of integral,
C  otherwise use analytical form of integral.
C
            IF (ABS(Q_ALPHA).LT.Q_MIN .OR. ABS(QM).LT.QM_MIN)  THEN  
              IF (Q_ALPHA .EQ. ZERO)  THEN
                IVAL = ( M_ALPHA(I,LEV,N)**2 - M_MIN**2 ) / 2.
              ELSE
                IVAL = ( QM**2/2.  + QM**3/3.  + QM**4/4.  + QM**5/5. 
     ^                 - QMM**2/2. - QMM**3/3. - QMM**4/4. - QMM**5/5. )
     ^                 / Q_ALPHA**2
              END IF
            ELSE
              IVAL = - ( ALOG(1.-QM) - ALOG(1.-QMM) + QM - QMM ) 
     ^                 / Q_ALPHA**2
            END IF
C
C  If I_ALPHA negative (due to round off error) then set it to zero.
C
            I_ALPHA(I,N) = AMAX1 ( IVAL, ZERO )
C
          ELSE
C
            I_ALPHA(I,N)  = ZERO
            DO_ALPHA(I,N) = .FALSE.
C
          END IF
C
 10     CONTINUE
C
      END IF
C
C  For integer value SLOPE = 2 (NOTE: only for M_MIN = 0!).
C
      IF (SLOPE .EQ. 2.)  THEN
C
        DO 20 N = 1,NAZ
        DO 20 I = IL1,IL2
C
          IF (DO_ALPHA(I,N))  THEN
C
            DRAG(I) = .TRUE.
            Q_ALPHA = V_ALPHA(I,N) / BVFB(I)
            QM = Q_ALPHA * M_ALPHA(I,LEV,N)
C
C  If |QM| is small then use first 4 terms series of Taylor series
C  expansion of integral in order to avoid indeterminate form of integral,
C  otherwise use analytical form of integral.
C
            IF (ABS(Q_ALPHA).LT.Q_MIN .OR. ABS(QM).LT.QM_MIN)  THEN  
              IF (Q_ALPHA .EQ. ZERO)  THEN
                IVAL = M_ALPHA(I,LEV,N)**3 / 3.
              ELSE
                IVAL = ( QM**3/3. + QM**4/4. + QM**5/5. + QM**6/6. ) 
     ^                 / Q_ALPHA**3
              END IF
            ELSE
              IVAL = - ( ALOG(1.-QM) + QM + QM**2/2.) / Q_ALPHA**3
            END IF
C
C  If I_ALPHA negative (due to round off error) then set it to zero.
C
            I_ALPHA(I,N) = AMAX1 ( IVAL, ZERO )
C
          ELSE
C
            I_ALPHA(I,N) = ZERO
C
          END IF
C
 20     CONTINUE
C
      END IF
C
C  For real value SLOPE = 1.5 (NOTE: only for M_MIN = 0!).
C
      IF (SLOPE .EQ. 1.5)  THEN
C
        DO 30 N = 1,NAZ
        DO 30 I = IL1,IL2
C
          IF (DO_ALPHA(I,N))  THEN
C
            DRAG(I) = .TRUE.
            Q_ALPHA = V_ALPHA(I,N) / BVFB(I)
            QM = Q_ALPHA * M_ALPHA(I,LEV,N)       
C
C  If |QM| is small then use first 4 terms series of Taylor series
C  expansion of integral in order to avoid indeterminate form of integral,
C  otherwise use analytical form of integral.
C
            IF (ABS(Q_ALPHA).LT.Q_MIN .OR. ABS(QM).LT.QM_MIN)  THEN  
              IF (Q_ALPHA .EQ. ZERO)  THEN
                IVAL = M_ALPHA(I,LEV,N)**2.5 / 2.5
              ELSE
                IVAL = ( QM/2.5    + QM**2/3.5 
     ^                 + QM**3/4.5 + QM**4/5.5 ) 
     ^                 * M_ALPHA(I,LEV,N)**1.5 / Q_ALPHA
              END IF
            ELSE
              QM     = ABS(QM)
              SQRTQM = SQRT(QM)
              IF (Q_ALPHA .GE. ZERO)  THEN
                IVAL = ( ALOG( (1.+SQRTQM)/(1.-SQRTQM) )
     ^                 -2.*SQRTQM*(1.+QM/3.) ) / Q_ALPHA**2.5
              ELSE
                IVAL = 2. * ( ATAN(SQRTQM) + SQRTQM*(QM/3.-1.) )
     ^                 / ABS(Q_ALPHA)**2.5
              END IF
            END IF
C
C  If I_ALPHA negative (due to round off error) then set it to zero.
C
            I_ALPHA(I,N) = AMAX1 ( IVAL, ZERO )
C
          ELSE
C
            I_ALPHA(I,N) = ZERO
C
          END IF
C
 30     CONTINUE
C
      END IF
C
      RETURN
C-----------------------------------------------------------------------
      END


      SUBROUTINE HINES_SMOOTH (DATA,WORK,DRAGIL,COEFF,NSMOOTH,
     1                         IL1,IL2,LEV1,LEV2,NLONS,NLEVS)
C
C  Smooth a longitude by altitude array in the vertical over a
C  specified number of levels using a three point smoother. 
C
C  NOTE: input array DATA is modified on output!
C
C  Aug. 3/95 - C. McLandress
C
C  Modifications:
C  --------------
C  Feb. 2/96 - C. McLandress (added logical flag)
C
C  Output arguement:
C  ----------------
C
C     * DATA    = smoothed array (on output).
C
C  Input arguements:
C  -----------------
C
C     * DATA    = unsmoothed array of data (on input).
C     * WORK    = work array of same dimension as DATA.
C     * DRAGIL  = logical flag indicating longitudes and levels where 
C     *           calculations to be performed.
C     * COEFF   = smoothing coefficient for a 1:COEFF:1 stencil.
C     *           (e.g., COEFF = 2 will result in a smoother which
C     *           weights the level L gridpoint by two and the two 
C     *           adjecent levels (L+1 and L-1) by one).
C     * NSMOOTH = number of times to smooth in vertical.
C     *           (e.g., NSMOOTH=1 means smoothed only once, 
C     *           NSMOOTH=2 means smoothing repeated twice, etc.)
C     * IL1     = first longitudinal index to use (IL1 >= 1).
C     * IL2     = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * LEV1    = first altitude level to use (LEV1 >=1). 
C     * LEV2    = last altitude level to use (LEV1 < LEV2 <= NLEVS).
C     * NLONS   = number of longitudes.
C     * NLEVS   = number of vertical levels.
C
C  Subroutine arguements.
C  ----------------------
C
      INTEGER  NSMOOTH, IL1, IL2, LEV1, LEV2, NLONS, NLEVS
      LOGICAL DRAGIL(NLONS,NLEVS)
      REAL  COEFF
      REAL  DATA(NLONS,NLEVS), WORK(NLONS,NLEVS)
C
C  Internal variables.
C  -------------------
C
      INTEGER  I, L, NS, LEV1P, LEV2M
      REAL  SUM_WTS
C-----------------------------------------------------------------------     
C
C  Calculate sum of weights.
C
      SUM_WTS = COEFF + 2.
C
      LEV1P = LEV1 + 1
      LEV2M = LEV2 - 1
C
C  Smooth NSMOOTH times
C
      DO 50 NS = 1,NSMOOTH
C
C  Copy data into work array.
C
        DO 20 L = LEV1,LEV2
        DO 20 I = IL1,IL2
          WORK(I,L) = DATA(I,L)
 20     CONTINUE
C
C  Smooth array WORK in vertical direction and put into DATA.
C
        DO 30 L = LEV1P,LEV2M
        DO 30 I = IL1,IL2
          IF (DRAGIL(I,L))  THEN
            DATA(I,L) = ( WORK(I,L+1) + COEFF*WORK(I,L) + WORK(I,L-1) ) 
     &                    / SUM_WTS 
          END IF
 30     CONTINUE
C
 50   CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END


      SUBROUTINE HINES_EXP (DATA,DATA_ZMAX,ALT,ALT_EXP,IORDER,
     1                      IL1,IL2,LEV1,LEV2,NLONS,NLEVS)
C
C  This routine exponentially damps a longitude by altitude array 
C  of data above a specified altitude.
C
C  Aug. 13/95 - C. McLandress
C
C  Output arguement:
C  -----------------
C
C     * DATA = modified data array.
C
C  Input arguements:
C  -----------------
C
C     * DATA    = original data array.
C     * ALT     = altitudes.
C     * ALT_EXP = altitude above which exponential decay applied.
C     * IORDER	= 1 means vertical levels are indexed from top down 
C     *           (i.e., highest level indexed 1 and lowest level NLEVS);
C     *           .NE. 1 highest level is index NLEVS.
C     * IL1     = first longitudinal index to use (IL1 >= 1).
C     * IL2     = last longitudinal index to use (IL1 <= IL2 <= NLONS).
C     * LEV1    = first altitude level to use (LEV1 >=1). 
C     * LEV2    = last altitude level to use (LEV1 < LEV2 <= NLEVS).
C     * NLONS   = number of longitudes.
C     * NLEVS   = number of vertical
C
C  Input work arrays:
C  ------------------
C
C     * DATA_ZMAX = data values just above altitude ALT_EXP.
C
      INTEGER  IORDER, IL1, IL2, LEV1, LEV2, NLONS, NLEVS
      REAL  ALT_EXP
      REAL  DATA(NLONS,NLEVS), DATA_ZMAX(NLONS), ALT(NLONS,NLEVS)
C
C  Internal variables.
C  -------------------
C
      INTEGER  LEVBOT, LEVTOP, LINCR, I, L
      REAL  HSCALE
      DATA  HSCALE / 5.E3 /
C-----------------------------------------------------------------------     
C
C  Index of lowest altitude level (bottom of drag calculation).
C
      LEVBOT = LEV2
      LEVTOP = LEV1
      LINCR  = 1
      IF (IORDER.NE.1)  THEN
        LEVBOT = LEV1
        LEVTOP = LEV2
        LINCR  = -1
      END IF
C
C  Data values at first level above ALT_EXP.
C
      DO 20 I = IL1,IL2
        DO 10 L = LEVTOP,LEVBOT,LINCR
          IF (ALT(I,L) .GE. ALT_EXP)  THEN
            DATA_ZMAX(I) = DATA(I,L) 
          END IF	   
 10     CONTINUE
 20   CONTINUE
C
C  Exponentially damp field above ALT_EXP to model top at L=1.
C
      DO 40 L = 1,LEV2 
        DO 30 I = IL1,IL2
          IF (ALT(I,L) .GE. ALT_EXP)  THEN
            DATA(I,L) = DATA_ZMAX(I) * EXP( (ALT_EXP-ALT(I,L))/HSCALE )
          END IF
 30     CONTINUE
 40   CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END


      SUBROUTINE HINES_PRNT1 (FLUX_U,FLUX_V,DRAG_U,DRAG_V,VEL_U,VEL_V,
     1                        ALT,SIGMA_T,SIGMA_ALPHA,M_ALPHA,
     2                        IU_PRINT,IV_PRINT,NMESSG,
     3                        ILPRT1,ILPRT2,LEVPRT1,LEVPRT2,
     4                        NAZ,NLONS,NLEVS,NAZMTH)
C
C  Print out altitude profiles of various quantities from
C  Hines Doppler spread gravity wave parameterization scheme.
C  (NOTE: only for NAZ = 4, 8 or 12). 
C
C  Aug. 8/95 - C. McLandress
C
C  Modifications:
C  --------------
C  Feb. 2/96 - C. McLandress (12 and 16 azimuths)
C
C  Input arguements:
C  -----------------
C
C     * IU_PRINT = 1 to print out values in east-west direction.
C     * IV_PRINT = 1 to print out values in north-south direction.
C     * NMESSG   = unit number for printed output.
C     * ILPRT1   = first longitudinal index to print.
C     * ILPRT2   = last longitudinal index to print.
C     * LEVPRT1  = first altitude level to print.
C     * LEVPRT2  = last altitude level to print.
C
      INTEGER  NAZ, ILPRT1, ILPRT2, LEVPRT1, LEVPRT2
      INTEGER  NLONS, NLEVS, NAZMTH
      INTEGER  IU_PRINT, IV_PRINT, NMESSG
      REAL  FLUX_U(NLONS,NLEVS), FLUX_V(NLONS,NLEVS)
      REAL  DRAG_U(NLONS,NLEVS), DRAG_V(NLONS,NLEVS)
      REAL  VEL_U(NLONS,NLEVS), VEL_V(NLONS,NLEVS)
      REAL  ALT(NLONS,NLEVS), SIGMA_T(NLONS,NLEVS)
      REAL  SIGMA_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH)
C
C  Internal variables.
C  -------------------
C
      INTEGER  N_EAST, N_WEST, N_NORTH, N_SOUTH
      INTEGER  I, L
C-----------------------------------------------------------------------
C
C  Azimuthal indices of cardinal directions.
C
      N_EAST = 1
      IF (NAZ.EQ.4)  THEN
        N_NORTH = 2
        N_WEST  = 3       
        N_SOUTH = 4       
      ELSE IF (NAZ.EQ.8)  THEN
        N_NORTH = 3
        N_WEST  = 5       
        N_SOUTH = 7       
      ELSE IF (NAZ.EQ.12)  THEN
        N_NORTH = 4
        N_WEST  = 7       
        N_SOUTH = 10       
      ELSE IF (NAZ.EQ.16)  THEN
        N_NORTH = 5
        N_WEST  = 9       
        N_SOUTH = 13       
      END IF
C
C  Print out values for range of longitudes.
C
      DO 100 I = ILPRT1,ILPRT2
C
C  Print east-west wind, sigmas, cutoff wavenumbers, flux and drag.
C
        IF (IU_PRINT.EQ.1)  THEN
          WRITE (NMESSG,*) 
          WRITE (NMESSG,6001) I
          WRITE (NMESSG,6005) 
 6001     FORMAT ( 'Hines GW (east-west) at longitude I =',I3)
 6005     FORMAT (15x,' U ',2x,'sig_E',2x,'sig_T',3x,'m_E',
     &            4x,'m_W',4x,'fluxU',5x,'gwdU')
          DO 10 L = LEVPRT1,LEVPRT2
            WRITE (NMESSG,6701) ALT(I,L)/1.E3, VEL_U(I,L),
     &                          SIGMA_ALPHA(I,L,N_EAST), SIGMA_T(I,L),
     &                          M_ALPHA(I,L,N_EAST)*1.E3, 
     &                          M_ALPHA(I,L,N_WEST)*1.E3,
     &                          FLUX_U(I,L)*1.E5, DRAG_U(I,L)*24.*3600.
  10      CONTINUE
 6701     FORMAT (' z=',f7.2,1x,3f7.1,2f7.3,f9.4,f9.3)
        END IF
C
C  Print north-south winds, sigmas, cutoff wavenumbers, flux and drag.
C
        IF (IV_PRINT.EQ.1)  THEN
          WRITE(NMESSG,*) 
          WRITE(NMESSG,6002) 
 6002     FORMAT ( 'Hines GW (north-south) at longitude I =',I3)
          WRITE(NMESSG,6006) 
 6006     FORMAT (15x,' V ',2x,'sig_N',2x,'sig_T',3x,'m_N',
     &            4x,'m_S',4x,'fluxV',5x,'gwdV')
          DO 20 L = LEVPRT1,LEVPRT2
            WRITE (NMESSG,6701) ALT(I,L)/1.E3, VEL_V(I,L),
     &                          SIGMA_ALPHA(I,L,N_NORTH), SIGMA_T(I,L),
     &                          M_ALPHA(I,L,N_NORTH)*1.E3, 
     &                          M_ALPHA(I,L,N_SOUTH)*1.E3,
     &                          FLUX_V(I,L)*1.E5, DRAG_V(I,L)*24.*3600.
 20       CONTINUE
        END IF
C
 100  CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END


      SUBROUTINE HINES_PRNT2 (M_ALPHA,SIGMA_ALPHA,VEL_U,VEL_V,
     1                        UBOT,VBOT,ALT,DRAGIL,V_ALPHA,
     2                        NMESSG,ILPRT1,ILPRT2,LEVPRT1,LEVPRT2,
     3                        NAZ,NLONS,NLEVS,NAZMTH)
C
C  Print out altitude profiles of cutoff wavenumbers, rms winds and
C  background winds at each horizontal azimuth for the Hines Doppler spread 
C  gravity wave parameterization scheme.
C
C  Feb. 2/96 - C. McLandress
C
C  Input arguements:
C  -----------------
C
C     * M_ALPHA      = cutoff wavenumber at each azimuth (1/m).
C     * SIGMA_ALPHA  = total rms wind in each azimuth (m/s).
C     * VEL_U    = background zonal wind component (m/s).
C     * VEL_V    = background meridional wind component (m/s).
C     * UBOT     = background zonal wind component at bottom level.
C     * VBOT     = background meridional wind component at bottom level.
C     * ALT      = altitude (m).
C     * DRAGIL   = logical flag indicating longitudes and levels where 
C     *            calculations to be performed.
C     * NMESSG   = unit number for printed output.
C     * ILPRT1   = first longitudinal index to print.
C     * ILPRT2   = last longitudinal index to print.
C     * LEVPRT1  = first altitude level to print.
C     * LEVPRT2  = last altitude level to print.
C     * NAZ      = actual number of horizontal azimuths used.
C     * NLONS    = number of longitudes.
C     * NLEVS    = number of vertical levels.
C     * NAZMTH   = azimuthal array dimension (NAZMTH >= NAZ).
C
C  Input work arrays:
C  ------------------
C
C     * V_ALPHA  = wind component at each azimuth (m/s). 
C
      INTEGER  NAZ, ILPRT1, ILPRT2, LEVPRT1, LEVPRT2, NMESSG
      INTEGER  NLONS, NLEVS, NAZMTH
      LOGICAL DRAGIL(NLONS,NLEVS)
      REAL  M_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  SIGMA_ALPHA(NLONS,NLEVS,NAZMTH)
      REAL  VEL_U(NLONS,NLEVS),    VEL_V(NLONS,NLEVS)
      REAL  ALT(NLONS,NLEVS), UBOT(NLONS),  VBOT(NLONS)
      REAL  V_ALPHA(NLONS,NAZMTH)
C
C Internal variables.
C -------------------
C
      INTEGER  IBIG, I, L, N, NAZ1
      PARAMETER  ( IBIG = 50 )
      REAL  ZKM, WORK(IBIG)
C-----------------------------------------------------------------------     
      NAZ1 = NAZ
      IF (NAZ.GT.12)  NAZ1 = 12
C
C  Print out values for range of longitudes.
C
      DO 100 I = ILPRT1,ILPRT2
C
C Print cutoff wavenumber at all azimuths.
C
        WRITE (NMESSG,*) 
        WRITE (NMESSG,6001) I
        WRITE (NMESSG,*) 
 6001   FORMAT ('Cutoff wavenumber (X 1.E3) at longitude I =',I3)
        DO 10 L = LEVPRT1,LEVPRT2
          ZKM = ALT(I,L)/1.E3
          DO 5 N = 1,NAZ1
            WORK(N) = M_ALPHA(I,L,N) * 1.E3
  5       CONTINUE
          WRITE (NMESSG,6100) ZKM, (WORK(N),N=1,NAZ1)
 10     CONTINUE
        IF (NAZ.GT.12)  THEN
          DO 11 L = LEVPRT1,LEVPRT2
            ZKM = ALT(I,L)/1.E3
            DO 6 N = 13,NAZ
              WORK(N) = M_ALPHA(I,L,N) * 1.E3
  6         CONTINUE
            WRITE (NMESSG,6100) ZKM, (WORK(N),N=13,NAZ)
 11       CONTINUE
        END IF
        WRITE (NMESSG,*) 
 6100   FORMAT (F5.1,'km',12F6.2)
C
C Print rms wind at all azimuths.
C
        WRITE (NMESSG,*) 
        WRITE (NMESSG,6002) I
        WRITE (NMESSG,*) 
 6002   FORMAT ('RMS wind (m/s) at longitude I =',I3)
        DO 20 L = LEVPRT1,LEVPRT2
          ZKM = ALT(I,L)/1.E3
          WRITE (NMESSG,6110) ZKM, (SIGMA_ALPHA(I,L,N),N=1,NAZ1)
 20     CONTINUE
        IF (NAZ.GT.12)  THEN
          DO 21 L = LEVPRT1,LEVPRT2
            ZKM = ALT(I,L)/1.E3
            WRITE (NMESSG,6110) ZKM, (SIGMA_ALPHA(I,L,N),N=13,NAZ)
 21       CONTINUE
        END IF
        WRITE (NMESSG,*) 
 6110   FORMAT (F5.1,'km',12F6.1)
C
C Print background wind at all azimuths.
C
        WRITE (NMESSG,*) 
        WRITE (NMESSG,6003) I
        WRITE (NMESSG,*) 
 6003   FORMAT ('Background wind (m/s) at longitude I =',I3)
        DO 30 L = LEVPRT1,LEVPRT2
          ZKM = ALT(I,L)/1.E3
          CALL HINES_WIND ( V_ALPHA, 
     ^                      VEL_U(1,L), VEL_V(1,L), UBOT, VBOT, 
     ^                      DRAGIL(1,L), NAZ, I, I, NLONS, NAZMTH )
          WRITE (NMESSG,6110) ZKM, (V_ALPHA(I,N),N=1,NAZ1)
 30     CONTINUE
        IF (NAZ.GT.12)  THEN
          DO 31 L = LEVPRT1,LEVPRT2
            ZKM = ALT(I,L)/1.E3
            CALL HINES_WIND ( V_ALPHA, 
     ^                        VEL_U(1,L), VEL_V(1,L), UBOT, VBOT, 
     ^                        DRAGIL(1,L), NAZ, I, I, NLONS, NAZMTH )
            WRITE (NMESSG,6110) ZKM, (V_ALPHA(I,N),N=13,NAZ)
 31       CONTINUE
        END IF
C
 100  CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
      END


      SUBROUTINE HINES_SETUP (IERROR,NAZ,SLOPE,F1,F2,F3,F5,F6,KSTAR,
     1                        M_MIN,ICUTOFF,ALT_CUTOFF,SMCO,NSMAX,
     2                        IHEATCAL,NMESSG,NAZMTH)
C
C  This routine specifies various parameters needed for the
C  the Hines Doppler spread gravity wave drag parameterization scheme.
C
C  Aug. 8/95 - C. McLandress
C
C  Modifications:
C  --------------
C  Feb. 2/96 - C. McLandress (minimum cutoff wavenumber M_MIN;
C                             12 and 16 azimuths; calculation of arrays
C                             K_ALPHA no longer done here)
C
C  Output arguements:
C  ------------------
C
C     * IERROR     = error flag.
C     *            = 0 no errors.
C     *            = 10 ==> NAZ > NAZMTH
C     *            = 20 ==> invalid number of azimuths (NAZ must be 4 or 8).
C     *            = 30 ==> invalid slope (SLOPE must be 1., 1.5 or 2.).
C     *            = 40 ==> invalid smoother (SMCO must be >= 1.)
C     * NAZ        = number of horizontal azimuths used (4, 8, 12 or 16).
C     * SLOPE      = slope of incident vertical wavenumber spectrum
C     *              (code set up presently for SLOPE 1., 1.5 or 2.).
C     * F1         = "fudge factor" used in calculation of trial value of
C     *              azimuthal cutoff wavenumber M_ALPHA (1.2 <= F1 <= 1.9).
C     * F2         = "fudge factor" used in calculation of maximum
C     *              permissible instabiliy-induced cutoff wavenumber 
C     *              M_SUB_M_TURB (0.1 <= F2 <= 1.4).
C     * F3         = "fudge factor" used in calculation of maximum 
C     *              permissible molecular viscosity-induced cutoff wavenumber 
C     *              M_SUB_M_MOL (0.1 <= F2 <= 1.4).
C     * F5         = "fudge factor" used in calculation of heating rate
C     *              (1 <= F5 <= 3).
C     * F6         = "fudge factor" used in calculation of turbulent 
C     *              diffusivity coefficient.
C     * KSTAR      = typical gravity wave horizontal wavenumber (1/m)
C     *              used in calculation of M_SUB_M_TURB.
C     * M_MIN      = minimum allowable cutoff vertical wavenumber, 
C     *              e.g., 1/(3km). This is used only for a spectral slope
C     *              (SLOPE) of one ==> for slope of 1.5 or 2 then M_MIN = 0. 
C     * ICUTOFF    = 1 to exponentially damp off GWD, heating and diffusion 
C     *              arrays above ALT_CUTOFF; otherwise arrays not modified.
C     * ALT_CUTOFF = altitude in meters above which exponential decay applied.
C     * SMCO       = smoother used to smooth cutoff vertical wavenumbers
C     *              and total rms winds before calculating drag or heating.
C     *              (==> a 1:SMCO:1 stencil used; SMCO >= 1.).
C     * NSMAX      = number of times smoother applied ( >= 1),
C     *            = 0 means no smoothing performed.
C     * IHEATCAL   = 1 to calculate heating rates and diffusion coefficient.
C     *            = 0 means only drag and flux calculated.
C
C  Input arguements:
C  -----------------
C
C     * NMESSG  = output unit number where messages to be printed.
C     * NAZMTH  = azimuthal array dimension (NAZMTH >= NAZ).
C
      INTEGER  NAZ, NAZMTH, IHEATCAL, ICUTOFF
      INTEGER  NMESSG, NSMAX, IERROR
      REAL  KSTAR, M_MIN, SLOPE, F1, F2, F3, F5, F6, ALT_CUTOFF, SMCO
C-----------------------------------------------------------------------     
C
C  Specify constants.
C
      NAZ   = 8
c      NAZ   = 12
c      SLOPE = 1.
c      SLOPE = 1.5
      SLOPE = 2.
      F1    = 1.5 
      F2    = 0.3 
      F3    = 1.0 
      F5    = 1.0
!      F6    = 0.1
      F6    = 0.04
      M_MIN = 0.
      IF (SLOPE.EQ.1.)  M_MIN = 1./3.E3
      KSTAR = 7.E-6
c      KSTAR = 14.E-6
c      ICUTOFF    = 0   
      ICUTOFF    = 1   
      ALT_CUTOFF = 105.E3
      SMCO       = 2.0 
      NSMAX      = 3
c      NSMAX      = 0
      IHEATCAL   = 1
c      IHEATCAL   = 0
C
C  Print information to output file.
C
c      WRITE (NMESSG,6000)
 6000 FORMAT (/' Subroutine HINES_SETUP:')
c      WRITE (NMESSG,*)  '  SLOPE = ', SLOPE
c      WRITE (NMESSG,*)  '  NAZ = ', NAZ
c      WRITE (NMESSG,*)  '  F1,F2,F3  = ', F1, F2, F3
c      WRITE (NMESSG,*)  '  F5,F6     = ', F5, F6
c      WRITE (NMESSG,*)  '  KSTAR     = ', KSTAR
c      WRITE (NMESSG,*)  '  M_MIN     = ', M_MIN
c      IF (ICUTOFF .EQ. 1)  THEN
c        WRITE (NMESSG,*) '  Drag exponentially damped above ',
c     &                       ALT_CUTOFF/1.E3
c      END IF
c      IF (NSMAX.LT.1 )  THEN
c        WRITE (NMESSG,*) '  No smoothing of cutoff wavenumbers, etc'
c      ELSE
c        WRITE (NMESSG,*) '  Cutoff wavenumbers and sig_t smoothed:'
c        WRITE (NMESSG,*) '    SMCO  =', SMCO
c        WRITE (NMESSG,*) '    NSMAX =', NSMAX
c      END IF
C
C  Check that things are setup correctly and log error if not
C
      IERROR = 0
      IF (NAZ .GT. NAZMTH)                                  IERROR = 10
      IF (NAZ.NE.4 .AND. NAZ.NE.8 
     & .AND. NAZ.NE.12 .AND. NAZ.NE.16 )                    IERROR = 20
      IF (SLOPE.NE.1. .AND. SLOPE.NE.1.5 .AND. SLOPE.NE.2.) IERROR = 30
      IF (SMCO .LT. 1.)                                     IERROR = 40
C
      RETURN
C-----------------------------------------------------------------------
      END

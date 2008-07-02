!WRF:MODEL_MP:PHYSICS
!
MODULE module_mp_etanew
      USE module_misc, only : wrf_dm_on_monitor
!
!-----------------------------------------------------------------------
      REAL,PRIVATE,SAVE ::  ABFR, CBFR, CIACW, CIACR, C_N0r0,           &
     &  CN0r0, CN0r_DMRmin, CN0r_DMRmax, CRACW, CRAUT, ESW0,            &
     &  RFmax, RQR_DR1, RQR_DR2, RQR_DR3, RQR_DRmin,                    &
     &  RQR_DRmax, RR_DRmin, RR_DR1, RR_DR2, RR_DR3, RR_DRmax
!
      INTEGER, PRIVATE,PARAMETER :: MY_T1=1, MY_T2=35
      REAL,PRIVATE,DIMENSION(MY_T1:MY_T2),SAVE :: MY_GROWTH
!
      REAL, PRIVATE,PARAMETER :: DMImin=.05e-3, DMImax=1.e-3,           &
     &      DelDMI=1.e-6,XMImin=1.e6*DMImin
      INTEGER, PUBLIC,PARAMETER :: XMImax=1.e6*DMImax, XMIexp=.0536,    &
     &                             MDImin=XMImin, MDImax=XMImax
      REAL, PRIVATE,DIMENSION(MDImin:MDImax) ::                         &
     &      ACCRI,SDENS,VSNOWI,VENTI1,VENTI2
!
      REAL, PRIVATE,PARAMETER :: DMRmin=.05e-3, DMRmax=.45e-3,          &
     &      DelDMR=1.e-6,XMRmin=1.e6*DMRmin, XMRmax=1.e6*DMRmax
      INTEGER, PRIVATE,PARAMETER :: MDRmin=XMRmin, MDRmax=XMRmax                   
      REAL, PRIVATE,DIMENSION(MDRmin:MDRmax)::                          &
     &      ACCRR,MASSR,RRATE,VRAIN,VENTR1,VENTR2
!
      INTEGER, PRIVATE,PARAMETER :: Nrime=40
      REAL, DIMENSION(2:9,0:Nrime),PRIVATE,SAVE :: VEL_RF
!
      INTEGER,PARAMETER :: NX=7501
      REAL, PARAMETER :: XMIN=180.0,XMAX=330.0
      REAL, DIMENSION(NX),PRIVATE,SAVE :: TBPVS,TBPVS0
      REAL, PRIVATE,SAVE :: C1XPVS0,C2XPVS0,C1XPVS,C2XPVS
!
      REAL, PRIVATE,PARAMETER ::                                        &
!--- Physical constants follow:
     &   CP=1004.6, EPSQ=1.E-12, GRAV=9.806, RHOL=1000., RD=287.04      &
     &  ,RV=461.5, T0C=273.15, XLS=2.834E6                              &
!--- Derived physical constants follow:
     &  ,EPS=RD/RV, EPS1=RV/RD-1., EPSQ1=1.001*EPSQ                     &
     &  ,RCP=1./CP, RCPRV=RCP/RV, RGRAV=1./GRAV, RRHOL=1./RHOL          &
     &  ,XLS1=XLS*RCP, XLS2=XLS*XLS*RCPRV, XLS3=XLS*XLS/RV              &
!--- Constants specific to the parameterization follow:
!--- CLIMIT/CLIMIT1 are lower limits for treating accumulated precipitation
     &  ,CLIMIT=10.*EPSQ, CLIMIT1=-CLIMIT                               &
     &  ,C1=1./3.                                                       &
     &  ,DMR1=.1E-3, DMR2=.2E-3, DMR3=.32E-3                            &
     &  ,XMR1=1.e6*DMR1, XMR2=1.e6*DMR2, XMR3=1.e6*DMR3
      INTEGER, PARAMETER :: MDR1=XMR1, MDR2=XMR2, MDR3=XMR3
!
! ======================================================================
!--- Important tunable parameters that are exported to other modules
!  * RHgrd - threshold relative humidity for onset of condensation
!  * T_ICE - temperature (C) threshold at which all remaining liquid water
!            is glaciated to ice
!  * T_ICE_init - maximum temperature (C) at which ice nucleation occurs
!  * NLImax - maximum number concentrations (m**-3) of large ice (snow/graupel/sleet) 
!  * NLImin - minimum number concentrations (m**-3) of large ice (snow/graupel/sleet) 
!  * N0r0 - assumed intercept (m**-4) of rain drops if drop diameters are between 0.2 and 0.45 mm
!  * N0rmin - minimum intercept (m**-4) for rain drops 
!  * NCW - number concentrations of cloud droplets (m**-3)
!  * FLARGE1, FLARGE2 - number fraction of large ice to total (large+snow) ice 
!          at T>0C and in presence of sublimation (FLARGE1), otherwise in
!          presence of ice saturated/supersaturated conditions
! ======================================================================
      REAL, PUBLIC,PARAMETER ::                                         &
     &  RHgrd=1.                                                        &
     & ,T_ICE=-30.                                                      &
     & ,T_ICEK=T0C+T_ICE                                                &
     & ,T_ICE_init=-5.                                                  &
     & ,NLImax=5.E3                                                     &
     & ,NLImin=1.E3                                                     &
     & ,N0r0=8.E6                                                       &
     & ,N0rmin=1.E4                                                     &
     & ,NCW=100.E6                                                      &
     & ,FLARGE1=1.                                                      &
     & ,FLARGE2=.2
!--- Other public variables passed to other routines:
      REAL,PUBLIC,SAVE ::  QAUT0
      REAL, PUBLIC,DIMENSION(MDImin:MDImax) :: MASSI
!
!
      CONTAINS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE ETAMP_NEW (itimestep,DT,DX,DY,                         &
     &                      dz8w,rho_phy,p_phy,pi_phy,th_phy,qv,qt,     &
     &                      LOWLYR,SR,                                  &
     &                      F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY,           &
     &                      QC,QR,QS,                                   &
     &                      mp_restart_state,tbpvs_state,tbpvs0_state,  &
     &                      RAINNC,RAINNCV,                             &
     &                      ids,ide, jds,jde, kds,kde,		        &
     &                      ims,ime, jms,jme, kms,kme,		        &
     &                      its,ite, jts,jte, kts,kte )
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: ITLO=-60, ITHI=40

      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE                     &
     &                     ,ITIMESTEP

      REAL, INTENT(IN) 	    :: DT,DX,DY
      REAL, INTENT(IN),     DIMENSION(ims:ime, kms:kme, jms:jme)::      &
     &                      dz8w,p_phy,pi_phy,rho_phy
      REAL, INTENT(INOUT),  DIMENSION(ims:ime, kms:kme, jms:jme)::      &
     &                      th_phy,qv,qt
      REAL, INTENT(INOUT),  DIMENSION(ims:ime, kms:kme, jms:jme ) ::    &
     &                      qc,qr,qs
      REAL, INTENT(INOUT),  DIMENSION(ims:ime, kms:kme, jms:jme ) ::    &
     &                      F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY
      REAL, INTENT(INOUT),  DIMENSION(ims:ime,jms:jme)           ::     &
     &                                                   RAINNC,RAINNCV
      REAL, INTENT(OUT),    DIMENSION(ims:ime,jms:jme):: SR
!
      REAL,DIMENSION(*),INTENT(INOUT) :: MP_RESTART_STATE
!
      REAL,DIMENSION(nx),INTENT(INOUT) :: TBPVS_STATE,TBPVS0_STATE
!
      INTEGER, DIMENSION( ims:ime, jms:jme ),INTENT(INOUT) :: LOWLYR

!-----------------------------------------------------------------------
!     LOCAL VARS
!-----------------------------------------------------------------------

!     NSTATS,QMAX,QTOT are diagnostic vars

      INTEGER,DIMENSION(ITLO:ITHI,4) :: NSTATS
      REAL,   DIMENSION(ITLO:ITHI,5) :: QMAX
      REAL,   DIMENSION(ITLO:ITHI,22):: QTOT

!     SOME VARS WILL BE USED FOR DATA ASSIMILATION (DON'T NEED THEM NOW). 
!     THEY ARE TREATED AS LOCAL VARS, BUT WILL BECOME STATE VARS IN THE 
!     FUTURE. SO, WE DECLARED THEM AS MEMORY SIZES FOR THE FUTURE USE

!     TLATGS_PHY,TRAIN_PHY,APREC,PREC,ACPREC,SR are not directly related 
!     the microphysics scheme. Instead, they will be used by Eta precip 
!     assimilation.

      REAL,  DIMENSION( ims:ime, kms:kme, jms:jme ) ::                  &
     &       TLATGS_PHY,TRAIN_PHY
      REAL,  DIMENSION(ims:ime,jms:jme):: APREC,PREC,ACPREC
      REAL,  DIMENSION(its:ite, kts:kte, jts:jte):: t_phy

      INTEGER :: I,J,K,KFLIP
      REAL :: WC
!
!-----------------------------------------------------------------------
!**********************************************************************
!-----------------------------------------------------------------------

      STOP 'ETANEW NOT IMPLEMENTED'
!
  END SUBROUTINE ETAMP_NEW
      END MODULE module_mp_etanew

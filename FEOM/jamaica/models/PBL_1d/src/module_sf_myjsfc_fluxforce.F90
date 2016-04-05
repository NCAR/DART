!----------------------------------------------------------------------
!
MODULE MODULE_SF_MYJSFC_fluxforce
!
!----------------------------------------------------------------------
!
      USE MODULE_MODEL_CONSTANTS
      use module_sf_myjsfc

! removed paraneters from module_sf_myjsfc
!
!----------------------------------------------------------------------
!
CONTAINS
!
!----------------------------------------------------------------------
  SUBROUTINE MYJSFC_fluxforce(ITIMESTEP,HT,DZ                     & 
     &                 ,PMID,PINT,TH,T,QV,QC,U,V,Q2                    &
     &                 ,TSK,QSFC,THZ0,QZ0,UZ0,VZ0                      &
     &                 ,LOWLYR,XLAND                                   &
     &                 ,USTAR,ZNT,Z0BASE,PBLH,MAVAIL,RMOL              &
     &                 ,AKHS,AKMS                                      &
     &                 ,CHS,CHS2,CQS2,HFX,QFX,FLX_LH,FLHC,FLQC         &
     &                 ,QGH,CPM,CT                                     &
     &          ,U10,V10,T02,TH02,TSHLTR,TH10,Q02,QSHLTR,Q10,PSHLTR    &
     &     ,calc_ust,z_o,z_t,z_q&
     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                        &
     &                 ,IMS,IME,JMS,JME,KMS,KME                        &
     &                 ,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: ITIMESTEP
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: LOWLYR
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: HT,MAVAIL      &
     &                                             ,XLAND,Z0BASE
!
!

      LOGICAL :: calc_ust

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(IN) :: DZ         &
     &                                                     ,PMID,PINT  &
     &                                                     ,Q2,QC,QV   &
     &                                                     ,T,TH       &
     &                                                     ,U,V   
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: FLX_LH,PSHLTR &
     &                                              ,Q10,QSHLTR    &
     &                                              ,TH10,TSHLTR,T02       &
     &                                              ,U10,V10,TH02,Q02

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(in) :: HFX,QFX

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(inout) :: TSK

!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: AKHS,AKMS       &
     &                                                ,PBLH,QSFC
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: QZ0,RMOL,THZ0   &
     &                                                ,USTAR,UZ0,VZ0   &
     &                                                ,ZNT
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: CHS,CHS2,CQS2     &
     &                                              ,CPM,CT,FLHC,FLQC  &
     &                                              ,QGH 

      REAL :: z_o,z_t,z_q

!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: I,J,K,KFLIP,LMH,LPBL,NTSD
!
      REAL :: A,APESFC,B,BTGX,CWMLOW                                   &
     &       ,DQDT,DTDIF,DTDT,DUDT,DVDT                                &
     &       ,FIS                                                      &
     &       ,P02P,P10P,PLOW,PSFC,PTOP,QLOW,QS02,QS10                  &
     &       ,RAPA,RAPA02,RAPA10,RATIOMX,RDZ,SEAMASK,SM                &
     &       ,T02P,T10P,TEM,TH02P,TH10P,THLOW,THELOW,THM               &
     &       ,TLOW,TZ0,ULOW,VLOW,ZSL
!
      REAL,DIMENSION(KTS:KTE) :: CWMK,PK,Q2K,QK,THEK,THK,TK,UK,VK
!
      REAL,DIMENSION(KTS:KTE-1) :: EL,ELM

!
      REAL,DIMENSION(KTS:KTE+1) :: ZHK
!
      REAL,DIMENSION(ITS:ITE,JTS:JTE) :: THSK
!
      REAL,DIMENSION(ITS:ITE,KTS:KTE+1,JTS:JTE) :: ZINT
!
!----------------------------------------------------------------------
!**********************************************************************
      REAL :: rlow

!----------------------------------------------------------------------
!
!***  MAKE PREPARATIONS
!
!----------------------------------------------------------------------
      DO J=JTS,JTE
      DO K=KTS,KTE+1
      DO I=ITS,ITE
        ZINT(I,K,J)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        ZINT(I,KTE+1,J)=HT(I,J)     ! Z at bottom of lowest sigma layer
        PBLH(I,J)=-1.
!
!!!!!!!!!
!!!!!! UNCOMMENT THESE LINES IF USING ETA COORDINATES
!!!!!!!!!
!!!!!!  ZINT(I,KTE+1,J)=1.E-4       ! Z of bottom of lowest eta layer
!!!!!!  ZHK(KTE+1)=1.E-4            ! Z of bottom of lowest eta layer
!
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO K=KTE,KTS,-1
        KFLIP=KTE+1-K
        DO I=ITS,ITE
          ZINT(I,K,J)=ZINT(I,K+1,J)+DZ(I,KFLIP,J)
        ENDDO
      ENDDO
      ENDDO
!
      NTSD=ITIMESTEP
!
      IF(NTSD==1)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          USTAR(I,J)=0.1
          FIS=HT(I,J)*G
          SM=XLAND(I,J)-1.
!!!       Z0(I,J)=SM*Z0SEA+(1.-SM)*(Z0(I,J)*Z0MAX+FIS*FCM+Z0LAND)
!!!       ZNT(I,J)=SM*Z0SEA+(1.-SM)*(ZNT(I,J)*Z0MAX+FIS*FCM+Z0LAND)
        ENDDO
        ENDDO
      ENDIF
!
!!!!  IF(NTSD==1)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          CT(I,J)=0.
        ENDDO
        ENDDO
!!!!  ENDIF
!
!----------------------------------------------------------------------
        setup_integration:  DO J=JTS,JTE
!----------------------------------------------------------------------
!
        DO I=ITS,ITE


!
!***  LOWEST LAYER ABOVE GROUND MUST BE FLIPPED
!
          LMH=KTE-LOWLYR(I,J)+1
!
          PTOP=PINT(I,KTE+1,J)      ! KTE+1=KME
          PSFC=PINT(I,LOWLYR(I,J),J)
! Define THSK here (for first timestep mostly)
          THSK(I,J)=TSK(I,J)/(PSFC*1.E-5)**CAPA
!
!***  CONVERT LAND MASK (1 FOR SEA; 0 FOR LAND)
!
          SEAMASK=XLAND(I,J)-1.
!
!***  FILL 1-D VERTICAL ARRAYS
!***  AND FLIP DIRECTION SINCE MYJ SCHEME
!***  COUNTS DOWNWARD FROM THE DOMAIN'S TOP
!
          DO K=KTE,KTS,-1
            KFLIP=KTE+1-K
            THK(K)=TH(I,KFLIP,J)
            TK(K)=T(I,KFLIP,J)
            RATIOMX=QV(I,KFLIP,J)
            QK(K)=RATIOMX/(1.+RATIOMX)
            PK(K)=PMID(I,KFLIP,J)
            CWMK(K)=QC(I,KFLIP,J)
            THEK(K)=(CWMK(K)*(-ELOCP/TK(K))+1.)*THK(K)
            Q2K(K)=2.*Q2(I,KFLIP,J)
!
!
!***  COMPUTE THE HEIGHTS OF THE LAYER INTERFACES
!
            ZHK(K)=ZINT(I,K,J)
!
          ENDDO
          ZHK(KTE+1)=HT(I,J)          ! Z at bottom of lowest sigma layer
!
          DO K=KTE,KTS,-1
            KFLIP=KTE+1-K
            UK(K)=U(I,KFLIP,J)
            VK(K)=V(I,KFLIP,J)
          ENDDO
!
!***  FIND THE HEIGHT OF THE PBL
          LPBL=LMH
          DO K=LMH-1,1,-1
            IF(Q2K(K)<=EPSQ2*FH) THEN
               LPBL=K
               GO TO 110
            ENDIF
!

          ENDDO

          LPBL=1
!
!--------------THE HEIGHT OF THE PBL------------------------------------
!-----------------------------------------------------------------------
!
 110      PBLH(I,J)=ZHK(LPBL)-ZHK(LMH+1)

!
!----------------------------------------------------------------------
!***
!***  FIND THE SURFACE EXCHANGE COEFFICIENTS
!***
!----------------------------------------------------------------------

          PLOW=PK(LMH)
          TLOW=TK(LMH)
          THLOW=THK(LMH)
          THELOW=THEK(LMH)
          QLOW=QK(LMH)
          CWMLOW=CWMK(LMH)
          ULOW=UK(LMH)
          VLOW=VK(LMH)
          ZSL=(ZHK(LMH)-ZHK(LMH+1))*0.5
          APESFC=(PSFC*1.E-5)**CAPA
          TZ0=THZ0(I,J)*APESFC


          CALL SFCDIF_fluxforce(NTSD,SEAMASK,THSK(I,J),QSFC(I,J),PSFC  &
     &               ,UZ0(I,J),VZ0(I,J),TZ0,THZ0(I,J),QZ0(I,J)         &
     &               ,USTAR(I,J),ZNT(I,J),Z0BASE(I,J),CT(I,J),RMOL(I,J) &
     &               ,AKMS(I,J),AKHS(I,J),PBLH(I,J),MAVAIL(I,J)        &
     &               ,CHS(I,J),CHS2(I,J),CQS2(I,J)                     &
     &               ,HFX(I,J),QFX(I,J),FLX_LH(I,J)                    &
     &               ,FLHC(I,J),FLQC(I,J),QGH(I,J),CPM(I,J)            &
     &               ,ULOW,VLOW,TLOW,THLOW,THELOW,QLOW,CWMLOW          &
     &               ,ZSL,PLOW                                         &
     &               ,U10(I,J),V10(I,J),TSHLTR(I,J),TH10(I,J)          &
     &               ,QSHLTR(I,J),Q10(I,J),PSHLTR(I,J)                 &
     &,calc_ust,z_o,z_t,z_q&
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &               ,IMS,IME,JMS,JME,KMS,KME                          &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE,i,j,ZHK(LMH+1))

!
!***  REMOVE SUPERATURATION AT 2M AND 10M
!

          tsk(i,j)=thsk(i,j)*(psfc/p1000mb)**rcp

          RAPA=APESFC
          TH02P=TSHLTR(I,J)
          TH10P=TH10(I,J)
          TH02(I,J)=TSHLTR(I,J)
!
          RAPA02=RAPA-GOCP02/TH02P
          RAPA10=RAPA-GOCP10/TH10P
!
          T02P=TH02P*RAPA02
          T10P=TH10P*RAPA10
! 1 may 06 tgs         T02(I,J) = T02P
          T02(I,J) = TH02(I,J)*APESFC
!
          P02P=(RAPA02**RCAP)*1.E5
          P10P=(RAPA10**RCAP)*1.E5
!
          QS02=PQ0/P02P*EXP(A2*(T02P-A3)/(T02P-A4))
          QS10=PQ0/P10P*EXP(A2*(T10P-A3)/(T10P-A4))
!
          IF(QSHLTR(I,J)>QS02)QSHLTR(I,J)=QS02
          IF(Q10   (I,J)>QS10)Q10   (I,J)=QS10
          Q02(I,J)=QSHLTR(I,J)/(1.-QSHLTR(I,J))
!----------------------------------------------------------------------
!
        ENDDO
!
!----------------------------------------------------------------------
      ENDDO setup_integration
!----------------------------------------------------------------------
 
    END SUBROUTINE MYJSFC_fluxforce
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
   SUBROUTINE SFCDIF_fluxforce(NTSD,SEAMASK,THS,QS,PSFC              &
     &                 ,UZ0,VZ0,TZ0,THZ0,QZ0                           &
     &                 ,USTAR,Z0,Z0BASE,CT,RLMO,AKMS,AKHS,PBLH,WETM    &
     &                 ,CHS,CHS2,CQS2,HFX,QFX,FLX_LH,FLHC,FLQC,QGH,CPM &
     &                 ,ULOW,VLOW,TLOW,THLOW,THELOW,QLOW,CWMLOW        &
     &                 ,ZSL,PLOW                                       &
     &                 ,U10,V10,TH02,TH10,Q02,Q10,PSHLTR               &
     &                 ,calc_ust,z_o,z_t,z_q&
     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                        &
     &                 ,IMS,IME,JMS,JME,KMS,KME                        &
     &                 ,ITS,ITE,JTS,JTE,KTS,KTE,i,j,ZSFC)
!     ****************************************************************
!     *                                                              *
!     *                       SURFACE LAYER                          *
!     *                                                              *
!     ****************************************************************
!----------------------------------------------------------------------
     USE module_ideal, ONLY : qsat

!
      IMPLICIT NONE
!
!----------------------------------------------------------------------


      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE,i,j
!
      INTEGER,INTENT(IN) :: NTSD
!

      LOGICAL :: calc_ust

      REAL,INTENT(IN) :: CWMLOW,PBLH,PLOW,QLOW,PSFC,SEAMASK,ZSFC       &
     &                  ,THELOW,THLOW,TLOW,TZ0,ULOW,VLOW,WETM,ZSL  &
     &                  ,Z0BASE

!
      REAL,INTENT(OUT) :: CHS,CHS2,CPM,CQS2,CT,FLHC,FLQC,FLX_LH    &
     &                   ,PSHLTR,Q02,Q10,QGH,RLMO,TH02,TH10,U10,V10
!
      REAL,INTENT(inOUT) :: THS

      REAL,INTENT(in) :: HFX,QFX
!
      REAL,INTENT(INOUT) :: AKHS,AKMS,QZ0,THZ0,USTAR,UZ0,VZ0,Z0,QS

      REAL :: z_o,z_t,z_q

!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: ITR,K
!
      REAL :: A,B,BTGH,BTGX,CXCHL,CXCHS,DTHV,DU2,ELFC,FCT              &
     &       ,HLFLX,HSFLX,HV,PSH02,PSH10,PSHZ,PSHZL,PSM10,PSMZ,PSMZL   &
     &       ,RDZ,RDZT,RIB,RLMA,RLMN,RLMP                              &
     &       ,RLOGT,RLOGU,RWGH,RZ,RZST,RZSU,SIMH,SIMM,TEM,THM          &
     &       ,UMFLX,USTARK,VMFLX,WGHT,WGHTT,WGHTQ,WSTAR2               &
     &       ,X,XLT,XLT4,XLU,XLU4,XT,XT4,XU,XU4,ZETALT,ZETALU          &
     &       ,ZETAT,ZETAU,ZQ,ZSLT,ZSLU,ZT,ZU,TOPOTERM,ZZIL
!
!***  DIAGNOSTICS
!
      REAL :: AKHS02,AKHS10,AKMS02,AKMS10,EKMS10,QSAT10,QSAT2          &
     &       ,RLNT02,RLNT10,RLNU10,SIMH02,SIMH10,SIMM10,T02,T10        &
     &       ,TERM1,RLOW,U10E,V10E,WSTAR,XLT02,XLT024,XLT10            &
     &       ,XLT104,XLU10,XLU104,XU10,XU104,ZT02,ZT10,ZTAT02,ZTAT10   &
     &       ,ZTAU,ZTAU10,ZU10,ZUUZ
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------

      REAL :: tskv,thskv,rhox,molength,zeta,zeta2,zeta10,yy,psih,xx,psim,&
           &hvflx

      REAL, PARAMETER :: alphastable= 5., alphaunstable=16.,&
           &a_stable = 1.,b_stable = 2./3.,c_stable = 5., d_stable = .35,&
           &epsilonmol=1.e-7


      RDZ=1./ZSL
      CXCHL=EXCML*RDZ
      CXCHS=EXCMS*RDZ
!
      BTGX=G/THLOW
      ELFC=VKARMAN*BTGX
!
      IF(PBLH>1000.)THEN
        BTGH=BTGX*PBLH
      ELSE
        BTGH=BTGX*1000.
      ENDIF



!
!----------------------------------------------------------------------
!
!***  SEA POINTS
!
!----------------------------------------------------------------------
!
      IF(SEAMASK>0.5)THEN 
!
!----------------------------------------------------------------------
        DO ITR=1,ITRMX
!----------------------------------------------------------------------
          Z0=MAX(USTFC*USTAR*USTAR,1.59E-5)
!
!***  VISCOUS SUBLAYER, JANJIC MWR 1994
!
!----------------------------------------------------------------------
          IF(USTAR<USTC)THEN
!----------------------------------------------------------------------
!
            IF(USTAR<USTR)THEN
!
              IF(NTSD==1)THEN
                 AKMS=CXCHS
                 AKHS=CXCHS
             ENDIF
!

              ZU=FZU1*SQRT(SQRT(Z0*USTAR*RVISC))/USTAR
              WGHT=AKMS*ZU*RVISC
              RWGH=WGHT/(WGHT+1.)
              UZ0=(ULOW*RWGH+UZ0)*0.5
              VZ0=(VLOW*RWGH+VZ0)*0.5

              z_o=zu
!
              ZT=FZT1*ZU
              ZQ=FZQ1*ZT
              z_t=zt
              z_q=zq
              WGHTT=AKHS*ZT*RTVISC
              WGHTQ=AKHS*ZQ*RQVISC

!
              IF(NTSD>1)THEN
!                THZ0=((WGHTT*THLOW+THS)/(WGHTT+1.)+THZ0)*0.5
!                QZ0=((WGHTQ*QLOW+QS)/(WGHTQ+1.)+QZ0)*0.5
              ELSE
                THZ0=(WGHTT*THLOW+THS)/(WGHTT+1.)
                QZ0=(WGHTQ*QLOW+QS)/(WGHTQ+1.)
              ENDIF
!
           ENDIF
!

           IF(USTAR>=USTR.AND.USTAR<USTC)THEN
              ZU=Z0
              UZ0=0.
              VZ0=0.
              z_o=zu
!
              ZT=FZT2*SQRT(SQRT(Z0*USTAR*RVISC))/USTAR
              ZQ=FZQ2*ZT
              z_t=zt
              z_q=zq
              WGHTT=AKHS*ZT*RTVISC
              WGHTQ=AKHS*ZQ*RQVISC
!
              IF(NTSD>1)THEN
!                THZ0=((WGHTT*THLOW+THS)/(WGHTT+1.)+THZ0)*0.5
!                QZ0=((WGHTQ*QLOW+QS)/(WGHTQ+1.)+QZ0)*0.5
              ELSE
                THZ0=(WGHTT*THLOW+THS)/(WGHTT+1.)
                QZ0=(WGHTQ*QLOW+QS)/(WGHTQ+1.)
              ENDIF
!
           ENDIF
!----------------------------------------------------------------------
        ELSE
!----------------------------------------------------------------------
            ZU=Z0
            UZ0=0.
            VZ0=0.
!
            ZT=Z0
            THZ0=THS
            QZ0=QS

!
            ZQ=Z0

            z_o=zu
            z_t=zt
            z_q=zq
            
!----------------------------------------------------------------------
         ENDIF
!----------------------------------------------------------------------

         IF (ITR==1) THEN
            RLOW=PLOW/(R_D*TLOW)
            HSFLX=-HFX/(rlow*cp)
         ENDIF
         
            z_t=zt
            z_q=z_t

         tskv=tlow*(1.+p608*qlow)
         rhox= psfc / (r_d * tskv)
!         rhox=rlow
         thskv=thlow*(1.+p608*qlow)

         hvflx=hfx/(cp*rhox)+p608*thskv*qfx/rhox
         molength=-ustar**3*thlow/&
              &(vkarman*g*hvflx)
         IF (ABS(molength) < epsilonmol) THEN
            molength=SIGN(1.,molength)*epsilonmol
         ENDIF

         rlmo=1./molength


!
!----------------------------------------------------------------------
!***  1./MONIN-OBUKHOV LENGTH
!----------------------------------------------------------------------
!
!            RLMO=ELFC*AKHS*DTHV/USTAR**3
!
            ZETALU=ZSLU*RLMO
            ZETALT=ZSLT*RLMO
            ZETAU=ZU*RLMO
            ZETAT=ZT*RLMO
!
            ZETALU=MIN(MAX(ZETALU,ZTMIN1),ZTMAX1)
            ZETALT=MIN(MAX(ZETALT,ZTMIN1),ZTMAX1)
            ZETAU=MIN(MAX(ZETAU,ZTMIN1/RZSU),ZTMAX1/RZSU)
            ZETAT=MIN(MAX(ZETAT,ZTMIN1/RZST),ZTMAX1/RZST)
!
!----------------------------------------------------------------------
!***   WATER FUNCTIONS
!----------------------------------------------------------------------
!
            RZ=(ZETAU-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZ=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
!
            RZ=(ZETALU-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZL=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
!
            SIMM=PSMZL-PSMZ+RLOGU
!
            RZ=(ZETAT-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZ=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
            RZ=(ZETALT-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZL=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
            SIMH=(PSHZL-PSHZ+RLOGT)*FH01
!----------------------------------------------------------------------
            USTARK=USTAR*VKARMAN
            AKMS=MAX(USTARK/SIMM,CXCHS)
            AKHS=MAX(USTARK/SIMH,CXCHS)



!
!----------------------------------------------------------------------
!***  BELJAARS CORRECTION FOR USTAR
!----------------------------------------------------------------------
!
            IF (hvflx > 0) THEN !- why not
               WSTAR2=WWST2*ABS(BTGH*hvflx)**(2./3.)
            ELSE
               WSTAR2=0.
            ENDIF


!----------------------------------------------------------------------
!         ENDIF  !  End of turbulent branch
!----------------------------------------------------------------------
!
        ENDDO  !  End of the iteration loop over sea points
!
!----------------------------------------------------------------------
!
!***  LAND POINTS
!
!----------------------------------------------------------------------
!
     ELSE  
!
!----------------------------------------------------------------------
!

        
!
        ZU=Z0
        UZ0=0.
        VZ0=0.
!
        ZT=ZU*ZTFC

        IF(NTSD>1)THEN
        ELSE
           THZ0=THS
           QZ0=QS
        ENDIF
!
!        THZ0=THS
!        QZ0=QS

        z_o=zu
        z_t=zt
        z_q=zq
        
        RLOW=PLOW/(R_D*TLOW)
!        HSFLX=-HFX/(rlow*cp)


        ZSLU=ZSL+ZU
        RZSU=ZSLU/ZU
        RLOGU=LOG(RZSU)
        ZSLT=ZSL+ZU


!mp   Topo modification of ZILFC term
!
          TOPOTERM=TOPOFAC*ZSFC**2.
          TOPOTERM=MAX(TOPOTERM,3.0)
!

          IF(hvflx < 0.)THEN
            ZZIL=ZILFC*TOPOTERM
          ELSE
            ZZIL=ZILFC
          ENDIF
!
!-


!----------------------------------------------------------------------
!
          land_point_iteration: DO ITR=1,ITRMX
!
!----------------------------------------------------------------------
!***  ZILITINKEVITCH FIX FOR ZT
!----------------------------------------------------------------------
!
! oldform   ZT=MAX(EXP(ZZIL*SQRT(USTAR*ZU))*ZU,EPSZT)
            ZT=MAX(EXP(ZZIL*SQRT(USTAR*Z0BASE))*Z0BASE,EPSZT)
            RZST=ZSLT/ZT
            RLOGT=LOG(RZST)

            z_t=zt
            z_q=z_t

         tskv=tlow*(1.+p608*qlow)
         rhox= psfc / (r_d * tskv)
!         rhox=rlow
         thskv=thlow*(1.+p608*qlow)

         hvflx=hfx/(cp*rhox)+p608*thskv*qfx/rhox

         molength=-ustar**3*thlow/&
              &(vkarman*g*hvflx)
         IF (ABS(molength) < epsilonmol) THEN
            molength=SIGN(1.,molength)*epsilonmol
         ENDIF

         rlmo=1./molength

!----------------------------------------------------------------------
!***  1./MONIN-OBUKHOV LENGTH-SCALE
!----------------------------------------------------------------------
!
!            RLMO=ELFC*AKHS*DTHV/USTAR**3

            ZETALU=ZSLU*RLMO
            ZETALT=ZSLT*RLMO
            ZETAU=ZU*RLMO
            ZETAT=ZT*RLMO

!
            ZETALU=MIN(MAX(ZETALU,ZTMIN2),ZTMAX2)
            ZETALT=MIN(MAX(ZETALT,ZTMIN2),ZTMAX2)
            ZETAU=MIN(MAX(ZETAU,ZTMIN2/RZSU),ZTMAX2/RZSU)
            ZETAT=MIN(MAX(ZETAT,ZTMIN2/RZST),ZTMAX2/RZST)
!
!----------------------------------------------------------------------
!***  LAND FUNCTIONS
!----------------------------------------------------------------------
!
            RZ=(ZETAU-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZ=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
!
            RZ=(ZETALU-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZL=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
!
            SIMM=PSMZL-PSMZ+RLOGU
!
            RZ=(ZETAT-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZ=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
            RZ=(ZETALT-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZL=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
            SIMH=(PSHZL-PSHZ+RLOGT)*FH02
!----------------------------------------------------------------------
            USTARK=USTAR*VKARMAN
            AKMS=MAX(USTARK/SIMM,CXCHL)
            AKHS=MAX(USTARK/SIMH,CXCHL)
            

!----------------------------------------------------------------------
!***  BELJAARS CORRECTION FOR USTAR
!----------------------------------------------------------------------
!

            IF (hvflx > 0) THEN !- why not
               WSTAR2=WWST2*ABS(BTGH*hvflx)**(2./3.)
            ELSE
               WSTAR2=0.
            ENDIF
               
            USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)
            
!
!----------------------------------------------------------------------
         ENDDO land_point_iteration  !  End of iteration for land points
!----------------------------------------------------------------------
!
!       ENDIF  !  End of turbulant branch over land
!
!----------------------------------------------------------------------
!
      ENDIF  !  End of land/sea branch

!no need to recalculate skin temp for water


!----------------------------------------------------------------------
!***  COUNTERGRADIENT FIX
!----------------------------------------------------------------------
!      HV=hvflx
!      IF(HV>0.)THEN
!        FCT=-10.*(BTGX)**(-1./3.)
!        CT=FCT*(HV/(PBLH*PBLH))**(2./3.)
!      ELSE
        CT=0.
!      ENDIF

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!***  THE FOLLOWING DIAGNOSTIC BLOCK PRODUCES 2-m and 10-m VALUES
!***  FOR TEMPERATURE, MOISTURE, AND WINDS.  IT IS DONE HERE SINCE
!***  THE VARIOUS QUANTITIES NEEDED FOR THE COMPUTATION ARE LOST
!***  UPON EXIT FROM THE ROTUINE.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      WSTAR=SQRT(WSTAR2)/WWST
!
      UMFLX=AKMS*(ULOW -UZ0 )
      VMFLX=AKMS*(VLOW -VZ0 )
      HSFLX=-HFX/(rhox*cp)
      HLFLX=-QFX/rhox
      
      AKHS=HSFLX/(THLOW-THZ0)
      qz0=qlow-HLFLX/akhs
      if (qlow==qsat(plow,tlow)) qz0=qlow
!      PRINT *,AKHS*(THLOW-THZ0),'111',thz0
!      PRINT *,AKHS*(qlow-qz0),'111',qz0

!      thz0=thlow-HSFLX/AKHS
!      ths=thz0
!      PRINT *,thz0,akhs,'111'
!      ths=THLOW-HSFLX/akhs
!      thz0=THLOW-HSFLX/akhs
!      dthv=HLFLX/akhs
!      ths=(thlow*(1.+p608*qlow)-dthv)/(1.+p608*qlow)

!      PRINT *,HSFLX,'111'
!      HSFLX=AKHS*(THLOW-THZ0)
!      PRINT *,HSFLX,'222'
!      HLFLX=AKHS*(QLOW -QZ0 )
!----------------------------------------------------------------------
!     IF(RIB>=RI)CTHEN
!----------------------------------------------------------------------
!       IF(SEAMASK>0.5)THEN
!         AKMS10=MAX( VISC/10.,CXCHS)
!         AKHS02=MAX(TVISC/02.,CXCHS)
!         AKHS10=MAX(TVISC/10.,CXCHS)
!       ELSE
!         AKMS10=MAX( VISC/10.,CXCHL)
!         AKHS02=MAX(TVISC/02.,CXCHL)
!         AKHS10=MAX(TVISC/10.,CXCHL)
!       ENDIF
!----------------------------------------------------------------------
!     ELSE
!----------------------------------------------------------------------
        ZU10=ZU+10.
        ZT02=ZT+02.
        ZT10=ZT+10.
!
        RLNU10=LOG(ZU10/ZU)
        RLNT02=LOG(ZT02/ZT)
        RLNT10=LOG(ZT10/ZT)

        ZTAU10=ZU10*RLMO
        ZTAT02=ZT02*RLMO
        ZTAT10=ZT10*RLMO
!
!----------------------------------------------------------------------
!***  SEA
!----------------------------------------------------------------------
!
        IF(SEAMASK>0.5)THEN
!
!----------------------------------------------------------------------
          ZTAU10=MIN(MAX(ZTAU10,ZTMIN1),ZTMAX1)
          ZTAT02=MIN(MAX(ZTAT02,ZTMIN1),ZTMAX1)
          ZTAT10=MIN(MAX(ZTAT10,ZTMIN1),ZTMAX1)
!----------------------------------------------------------------------
          RZ=(ZTAU10-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSM10=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
!
          SIMM10=PSM10-PSMZ+RLNU10
!
          RZ=(ZTAT02-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH02=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
          SIMH02=(PSH02-PSHZ+RLNT02)*FH01
!
          RZ=(ZTAT10-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH10=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
          SIMH10=(PSH10-PSHZ+RLNT10)*FH01
!
          AKMS10=MAX(USTARK/SIMM10,CXCHS)
          AKHS02=MAX(USTARK/SIMH02,CXCHS)
          AKHS10=MAX(USTARK/SIMH10,CXCHS)
!
!----------------------------------------------------------------------
!***  LAND
!----------------------------------------------------------------------
!
        ELSE
!
!----------------------------------------------------------------------
          ZTAU10=MIN(MAX(ZTAU10,ZTMIN2),ZTMAX2)
          ZTAT02=MIN(MAX(ZTAT02,ZTMIN2),ZTMAX2)
          ZTAT10=MIN(MAX(ZTAT10,ZTMIN2),ZTMAX2)
!----------------------------------------------------------------------
          RZ=(ZTAU10-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSM10=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
!
          SIMM10=PSM10-PSMZ+RLNU10
!
          RZ=(ZTAT02-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH02=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
          SIMH02=(PSH02-PSHZ+RLNT02)*FH02
!
          RZ=(ZTAT10-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH10=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
          SIMH10=(PSH10-PSHZ+RLNT10)*FH02
!
          AKMS10=MAX(USTARK/SIMM10,CXCHL)
          AKHS02=MAX(USTARK/SIMH02,CXCHL)
          AKHS10=MAX(USTARK/SIMH10,CXCHL)
!----------------------------------------------------------------------
        ENDIF
!----------------------------------------------------------------------
!     ENDIF
!----------------------------------------------------------------------
      U10 =UMFLX/AKMS10+UZ0
      V10 =VMFLX/AKMS10+VZ0
      TH02=HSFLX/AKHS02+THZ0
      TH10=HSFLX/AKHS10+THZ0
      Q02 =HLFLX/AKHS02+QZ0
      Q10 =HLFLX/AKHS10+QZ0
      TERM1=-0.068283/TLOW
      PSHLTR=PSFC*EXP(TERM1)
!
!----------------------------------------------------------------------
!***  COMPUTE "EQUIVALENT" Z0 TO APPROXIMATE LOCAL SHELTER READINGS.
!----------------------------------------------------------------------
!
      U10E=U10
      V10E=V10
!
      IF(SEAMASK<0.5)THEN
!
        ZUUZ=AMIN1(ZU*0.50,0.18)
        ZU=AMAX1(ZU*0.35,ZUUZ)

!
        ZU10=ZU+10.
        RZSU=ZU10/ZU
        RLNU10=LOG(RZSU)
!
        ZETAU=ZU*RLMO
        ZTAU10=ZU10*RLMO
!
        ZTAU10=MIN(MAX(ZTAU10,ZTMIN2),ZTMAX2)
        ZETAU=MIN(MAX(ZETAU,ZTMIN2/RZSU),ZTMAX2/RZSU)
!
        RZ=(ZTAU10-ZTMIN2)/DZETA2
        K=INT(RZ)
        RDZT=RZ-REAL(K)
        K=MIN(K,KZTM2)
        K=MAX(K,0)
        PSM10=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
        SIMM10=PSM10-PSMZ+RLNU10
        EKMS10=MAX(USTARK/SIMM10,CXCHL)
!
        U10E=UMFLX/EKMS10+UZ0
        V10E=VMFLX/EKMS10+VZ0
!
     ENDIF
!
      U10=U10E
      V10=V10E
!
!----------------------------------------------------------------------
!***  SET OTHER WRF DRIVER ARRAYS
!----------------------------------------------------------------------
!




      CHS=AKHS
      CHS2=AKHS02
      CQS2=AKHS02
!      HFX=-RLOW*CP*HSFLX
!      QFX=-RLOW*HLFLX*WETM

!no diagnostics for 2 and 10m for now
1000  CONTINUE

      FLX_LH=XLV*QFX
!      FLHC=RLOW*CP*AKHS
      FLHC=rhox*CP*AKHS
!      FLQC=RLOW*AKHS*WETM
      FLQC=rhox*AKHS*WETM
!!!   QGH=PQ0/PSHLTR*EXP(A2S*(TSK-A3S)/(TSK-A4S))
      QGH=((1.-SEAMASK)*PQ0+SEAMASK*PQ0SEA)                            &
     &     /PLOW*EXP(A2S*(TLOW-A3S)/(TLOW-A4S))
      QGH=QGH/(1.-QGH)    !Convert to mixing ratio
      CPM=CP*(1.+0.8*QLOW)
!
!***  DO NOT COMPUTE QS OVER LAND POINTS HERE SINCE IT IS
!***  A PROGNOSTIC VARIABLE THERE.  IT IS OKAY TO USE IT 
!***  AS A DIAGNOSTIC OVER WATER SINCE IT WILL CAUSE NO
!***  INTERFERENCE BEFORE BEING RECOMPUTED IN MYJPBL.
!
      IF(SEAMASK>0.5)THEN
        QS=QLOW+QFX/(RLOW*AKHS)
        QS=QS/(1.-QS)
      ENDIF
!----------------------------------------------------------------------
!
    END SUBROUTINE SFCDIF_FLUXFORCE
!
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
    END MODULE MODULE_SF_MYJSFC_FLUXFORCE
!
!----------------------------------------------------------------------

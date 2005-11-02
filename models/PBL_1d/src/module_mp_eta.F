!WRF:MODEL_MP:PHYSICS
!
MODULE module_mp_eta

CONTAINS
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
   SUBROUTINE ETAMP    (itimestep,DT,XLAND,pi3d,p8w,P_PHY,TH,QV,QL,  &
                        RAINNC,PREC,T0ETA,Q0ETA,P0ETA,               &
                        IDS,IDE,JDS,JDE,KDS,KDE,                     &
                        IMS,IME,JMS,JME,KMS,KME,                     &
                        ITS,ITE,JTS,JTE,KTS,KTE                      )
!---------------------------------------------------------------------
      IMPLICIT NONE
!---------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                  &
                           ,IMS,IME,JMS,JME,KMS,KME                  &
                           ,ITS,ITE,JTS,JTE,KTS,KTE                  &
                           ,itimestep

      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme)::p8w,p_phy,pi3d
      REAL, INTENT(IN), DIMENSION(ims:ime, jms:jme):: XLAND
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme):: P0ETA
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme):: &
                                                    QV,QL,Q0ETA,T0ETA,TH
      REAL, INTENT(INOUT), DIMENSION(ims:ime, jms:jme):: PREC,RAINNC

      REAL, INTENT(IN)                     :: DT
      REAL, DIMENSION(its:ite, kts:kte, jts:jte)::TFLIP,QFLIP,QLFLIP,PFLIP,T
      REAL, DIMENSION(its:ite, kts:kte+1, jts:jte)::P8WFLIP
      REAL, DIMENSION(its:ite, kts:kte, jts:jte)::T0FLIP,Q0FLIP,P0FLIP

      INTEGER :: I,J,K,KFLIP

!!! for itimestep = 1, T0ETA,Q0ETA,P0ETA are approximation only

      DO J = jts,jte
      DO K = KTS,KTE
      DO I = its,ite
         T(I,K,J) = TH(i,k,j)*pi3d(i,k,j)
      ENDDO
      ENDDO
      ENDDO

      if (itimestep .eq. 1)  then
      DO J = jts,jte
      DO K = KTS,KTE
      DO I = its,ite
         T0ETA(I,K,J) = T(I,K,J)
         Q0ETA(I,K,J) = MAX(0., QV(I,K,J))
         P0ETA(I,K,J) = P_PHY(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      endif


      DO K = KTS,KTE
         KFLIP=KTE+1-K
         DO J = jts,jte
         DO I = its,ite
            TFLIP (I,K,J) = T(I,KFLIP,J)
            QFLIP (I,K,J) = MAX(0.,QV(I,KFLIP,J)/(1.+QV(I,KFLIP,J)))
            QLFLIP(I,K,J) = QL(I,KFLIP,J)
            PFLIP (I,K,J) = P_PHY(I,KFLIP,J)
            T0FLIP(I,K,J) = T0ETA(I,KFLIP,J)
            Q0FLIP(I,K,J) = MAX(0.,Q0ETA(I,KFLIP,J)/(1.+Q0ETA(I,KFLIP,J)))
            P0FLIP(I,K,J) = P0ETA(I,KFLIP,J)
         ENDDO
         ENDDO
      ENDDO

      DO K = KMS,KME
         KFLIP=KME+1-K
         DO J = jts,jte
         DO I = its,ite
            P8WFLIP(I,K,J) = P8W(I,KFLIP,J)
         ENDDO
         ENDDO
      ENDDO

      CALL GSCOND (DT,1,XLAND,PFLIP,P0FLIP,                           &
                   TFLIP,QFLIP,QLFLIP,T0FLIP,Q0FLIP,                  &
                   ids,ide, jds,jde, kds,kde,                         &
                   ims,ime, jms,jme, kms,kme,                         &
                   its,ite, jts,jte, kts,kte                          )

      CALL PRECPD (DT,1,QFLIP,QLFLIP,Q0FLIP,T0FLIP,TFLIP,XLAND,PFLIP, &
                   P8WFLIP,PREC,RAINNC,                               &
                   ids,ide, jds,jde, kds,kde,                         &
                   ims,ime, jms,jme, kms,kme,                         &
                   its,ite, jts,jte, kts,kte                          )
   
      DO J = jts,jte
      DO K = KTS,KTE
         KFLIP=KTE+1-K
         DO I = its,ite
            T(I,KFLIP,J)     = TFLIP (I,K,J) 
            QV(I,KFLIP,J)    = MAX(0., QFLIP(I,K,J)/(1. - QFLIP(I,K,J)))
            QL(I,KFLIP,J)    = QLFLIP(I,K,J) 
            T0ETA(I,KFLIP,J) = T0FLIP(I,K,J) 
            Q0ETA(I,KFLIP,J) = MAX(0., Q0FLIP(I,K,J)/(1. - Q0FLIP(I,K,J)))
            P0ETA(I,K,J)     = P_PHY(I,K,J)
         ENDDO
      ENDDO
      ENDDO

      DO K = KTS,KTE
      DO J = jts,jte
      DO I = its,ite
         TH(I,K,J) = T(i,k,j)/pi3d(i,k,j)
      ENDDO
      ENDDO
      ENDDO

 
   END SUBROUTINE ETAMP 

!----------------------------------------------------------------------
 SUBROUTINE GSCOND(DT,NPHS,XLAND,PFLIP,P0FLIP,                        &
                   T,Q,CWM,T0,Q0,                                     &
                   ids,ide, jds,jde, kds,kde,                         &
                   ims,ime, jms,jme, kms,kme,                         &
                   its,ite, jts,jte, kts,kte                          )
!----------------------------------------------------------------------
 IMPLICIT NONE
!----------------------------------------------------------------------
       REAL, PARAMETER :: A2=17.2693882,A3=273.16,A4=35.86,           &
                          PQ0=379.90516, CP=1004.6,ELWV=2.50E6,       &
                          ELIV=2.834E6,EPSQ=2.E-12,R=287.04,CPR=CP*R
       REAL, PARAMETER :: RCP=1./CP
!---------------------------------------------------------------------- 
       INTEGER, INTENT(IN)        :: ids,ide, jds,jde, kds,kde ,      &
                                     ims,ime, jms,jme, kms,kme ,      &
                                     its,ite, jts,jte, kts,kte

       REAL, INTENT(INOUT), DIMENSION(its:ite, kts:kte, jts:jte):: &
                                                        Q,CWM,Q0,T0,T
       REAL, INTENT(IN), DIMENSION(its:ite, kts:kte, jts:jte)::PFLIP,P0FLIP
       REAL, INTENT(IN), DIMENSION(ims:ime, jms:jme)::XLAND

       INTEGER, INTENT(IN)                     :: NPHS
       REAL,    INTENT(IN)                     :: DT

       INTEGER, DIMENSION(its:ite, jts:jte):: LMH
       REAL,    DIMENSION(its:ite, jts:jte):: U00
       REAL,    DIMENSION(2*kte)   :: UL

!      INTEGER, DIMENSION(ims:ime, kms;kme, jms:jme):: HTM
!----------------------------------------------------------------------
       INTEGER:: IW(kts:kte)
       REAL,DIMENSION(kts:kte, its:ite, jts:jte):: T_T,T0_T,Q_T,Q0_T, &
                                                   CWM_T
!                                                  ,HTM_T
!      REAL,DIMENSION(its:ite,jts:jte)          :: PDSL
!      REAL,DIMENSION(its:ite,jts:jte)          :: RES
!---------------------------------------------------------------------- 
       INTEGER :: LM,J,I,L
       INTEGER :: LMHIJ,LML,IWKL
       REAL    :: MR,KE,DTPH,RDTPH,TWODT,RTWODT,C0,C1,C2,US,EPS,CCLIMIT
       REAL    :: CLIMIT,UTIM,U00IJ,P0IJ,RESIJ,QKL,E0,COND,CWMKL
       REAL    :: PDSLIJ,TKL,TMT0,TMT15,AI,BIQW,QI,QINT,U00KL,FITHH,PP
       REAL    :: PP0,AT,AQ,AP,FIW,EC,CCRKL,RQKLL,CCR,RQKL,QC,ELV,BI,QW,FI
       REAL    :: THH,US00,CCRKL1,AA,AB,AC,AD,AE,AF,AG
       REAL    :: CONDK,QTEMP,RQTMP,CONE0
!-------------------------------------------------------------
!-----------------------------------------------------------------------
!--------------PREPARATORY CALCULATIONS---------------------------------

      LM = kte
      DTPH =float(NPHS)*DT
      RDTPH=1./DTPH
      TWODT=DTPH
      RTWODT=1./TWODT                                                   
      C0=1.5E-4                                                         
      C1=300.                                                          
      C2=0.5                                                            
      MR=3.0E-4                                                         
      KE=2.0E-5                                                         
      US=1.
      EPS=0.622
      CCLIMIT=1.0E-3                                                    
      CLIMIT =1.0E-20                                                   

      do j = jts,jte
      do i = its,ite
         LMH(i,j) = kme-1
!        RES(i,j) = 1.
      enddo
      enddo

!     do j = jts,jte
!     do i = its,ite
!     do k = kms,kme
!        HTM(i,k,j) = 1.
!     enddo
!     enddo
!     enddo

      do j = jts,jte
      do i = its,ite
         U00(I,J)=(2.-XLAND(i,j))*0.75+(XLAND(i,j) -1.)*0.80
      enddo
      enddo

      DO L=1,2*LM
      IF(L.GE.LM-10.AND.L.LE.LM)THEN
        UL(L)=0.1*FLOAT(L-LM+10)
      ELSE
        UL(L)=0.
      ENDIF
      ENDDO


!-----------------------------------------------------------------------
!------------------PADDING SPECIFIC HUMIDITY & CWM IF TOO SMALL---------
!-----------------------------------------------------------------------
      DO J=jts,jte
      DO L=1,LM                              
      DO I=its,ite
        IF(Q(I,L,J).LT.EPSQ)Q(I,L,J)=EPSQ  ! *HTM(I,L,J)                       
        IF(CWM(I,L,J).LT.CLIMIT)CWM(I,L,J)=CLIMIT  ! *HTM(I,L,J)                   
      ENDDO
      ENDDO
      ENDDO
!
!     DO J=jts,jte 
!     DO I=its,ite
!       PDSL(I,J)=RES(I,J)*PD(I,J)                                              
!     ENDDO
!     ENDDO
!
      IW(1)=0
      UTIM=1.
!
!-----------------------------------------------------------------------
!*************BEGINNING OF GRID-SCALE CONDENSATION/EVAP. LOOP***********
!-----------------------------------------------------------------------
!***
!***  TRANSPOSE ARRAYS
!***
!     CALL SGETMO(T,LDA,LDA,LM,T_T,LM)
!     CALL SGETMO(Q,LDA,LDA,LM,Q_T,LM)
!     CALL SGETMO(HTM,LDA,LDA,LM,HTM_T,LM)
!     CALL SGETMO(CWM,LDA,LDA,LM,CWM_T,LM)
!     CALL SGETMO(T0,LDA,LDA,LM,T0_T,LM)
!     CALL SGETMO(Q0,LDA,LDA,LM,Q0_T,LM)

      CALL SGETMO(T,T_T,                   &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )
      CALL SGETMO(Q,Q_T,                   &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )
!     CALL SGETMO(HTM,HTM_T,               &
!                 ids,ide,jds,jde,kds,kde, &
!                 ims,ime,jms,jme,kms,kme, &
!                 its,ite,jts,jte,kts,kte  )
      CALL SGETMO(CWM,CWM_T,               &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )

      CALL SGETMO(T0,T0_T,                 &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )
      CALL SGETMO(Q0,Q0_T,                 &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )
!
!-----------------------------------------------------------------------
!------------------QW, QI AND QINT--------------------------------------
!-----------------------------------------------------------------------
      DO 100 J=jts,jte
      DO 100 I=its,ite
!-----------------------------------------------------------------------
!
      LMHIJ=LMH(I,J)
!     HBM2IJ=HBM2(I,J)
      U00IJ=U00(I,J)
!     P0IJ=P0(I,J)
!     RESIJ=RES(I,J)
!     PDSLIJ=PDSL(I,J)
!
      DO 90 L=2,LM
!
      TKL=T_T(L,I,J)
      QKL=Q_T(L,I,J)
      CWMKL=CWM_T(L,I,J)
!
      COND=0.
      E0=0.
      LML=LM-LMHIJ
!     HH=HTM_T(L,I,J)*HBM2IJ
!     HH=HBM2IJ
      TMT0=(TKL-273.16) ! *HH                                                
      TMT15=AMIN1(TMT0,-15.) ! *HH                                            
      AI=0.008855
      BI=1.
!
      IF(TMT0.LT.-20.)THEN
        AI=0.007225
        BI=0.9674
      ENDIF
!
!     QW=HH*PQ0/(PDSLIJ*AETA(L)+PT(i,j))*EXP(HH*A2*(T(I,L,J)-A3)/(T(I,L,J)-A4))                  
      QW=PQ0/PFLIP(I,L,J)*EXP(A2*(T(I,L,J)-A3)/(T(I,L,J)-A4))                  
      QI=QW*(BI+AI*AMIN1(TMT0,0.))                               
      QINT=QW*(1.-0.00032*TMT15*(TMT15+15.))                        
      IF(TMT0.LE.-40.)QINT=QI                                    
!-----------------------------------------------------------------------
!-------------------ICE-WATER ID NUMBER IW------------------------------
!-----------------------------------------------------------------------
      IF(TMT0.LT.-15.)THEN
        U00KL=U00IJ+UL(L+LML)*(0.95-U00IJ)*UTIM
        FI=Q(I,L,J)-U00KL*QI
        IF(FI.GT.0..OR.CWMKL.GT.CLIMIT)THEN                    
          IW(L)=1                                                   
        ELSE                                                           
          IW(L)=0                                                   
        ENDIF                                                         
      ENDIF
!
      IF(TMT0.GE.0.)THEN
        IW(L)=0                                                      
      ENDIF
!
      IF(TMT0.LT.0.0.AND.TMT0.GE.-15.)THEN
        IW(L)=0
        IF(IW(L-1).EQ.1.AND.CWMKL.GT.CLIMIT)IW(L)=1 
      ENDIF
!-----------------------------------------------------------------------
!--------------CONDENSATION AND EVAPORATION OF CLOUD--------------------
!------------------------AT, AQ AND DP/DT-------------------------------
!-----------------------------------------------------------------------
      THH=TWODT ! *HH                                                      
      PP=PFLIP(I,L,J)
!     PP0=P0IJ*RESIJ*AETA(L)+PT(i,j)
      PP0=P0FLIP(I,L,J)
      AT=(TKL-T0_T(L,I,J))*RTWODT                                        
      AQ=(QKL-Q0_T(L,I,J))*RTWODT                                        
      AP=(PP-PP0)*RTWODT         
      IWKL=IW(L)                                                      
      U00KL=U00IJ+UL(L+LML)*(0.95-U00IJ)*UTIM
!-----------------------------------------------------------------------
!----------------THE SATUATION SPECIFIC HUMIDITY------------------------
!-----------------------------------------------------------------------
      FIW=FLOAT(IWKL)                                                   
      ELV=(1.-FIW)*ELWV+FIW*ELIV                                     
      QC =(1.-FIW)*QINT+FIW*QI
!-----------------------------------------------------------------------
!----------------THE RELATIVE HUMIDITY----------------------------------
!-----------------------------------------------------------------------
      IF(QC.LE.0.)THEN                                                
        RQKL=0.
      ELSE                                                             
        RQKL=QKL/QC                                                 
      ENDIF                                                             
!-----------------------------------------------------------------------
!----------------CLOUD COVER RATIO CCR----------------------------------
!-----------------------------------------------------------------------
      IF(RQKL.LE.U00KL)THEN                                      
        CCR=0.
      ELSE                                                            
        RQKLL=AMIN1(US,RQKL)                                         
        CCR=1.-SQRT((US-RQKLL)/(US-U00KL))                     
      ENDIF                                                             
!-----------------------------------------------------------------------
!-----------CORRECT CCR IF IT IS TOO SMALL IN LARGE CWM REGIONS--------
!-----------------------------------------------------------------------
      IF(CCR.GE.0.01.AND.CCR.LE.0.2.AND.CWMKL.GE.0.2E-3)THEN
        CCR=AMIN1(1.,CWMKL*1.E3)
      ENDIF
!
      CCRKL=CCR
!-----------------------------------------------------------------------
!-------GIVE UP THIS POINT  IF NO CLOUD NOR CONDENSATION EXIST--------- 
!-----------------------------------------------------------------------
      IF(CCRKL.LE.CCLIMIT.AND.CWMKL.LE.CLIMIT)GO TO 90
!-----------------------------------------------------------------------
!----------------EVAPORATION OF CLOUD WATER-----------------------------
!-----------------------------------------------------------------------
        EC=0.
        IF(CCRKL.LE.CCLIMIT.AND.CWMKL.GT.CLIMIT)THEN 
          EC=QC*(U00KL-RQKL)*RTWODT                           
          E0=AMAX1(EC,0.0)                                             
          E0=AMIN1(CWMKL*RTWODT,E0) ! *HH                           
          E0=AMAX1(0.,E0)                                          
        ENDIF
!-----------------------------------------------------------------------
!----------------CONDENSATION OF CLOUD----------------------------------
!-----------------------------------------------------------------------
      IF(CCRKL.LE.0.20.OR.QC.LE.EPSQ)THEN                        
!     IF(CCRKL.LE.CCLIMIT.OR.QC.LE.EPSQ)THEN                        
        COND=0.                                                    
        GO TO 80                                                      
      ENDIF                                                            
!-----------------------------------------------------------------------
!-----------THE EQS. FOR COND. HAS BEEN REORGANIZED TO REDUCE CPU------
!-----------------------------------------------------------------------
      US00=US-U00KL 
      CCRKL1=1.-CCRKL
      AA=EPS*ELV*PP*QKL
      AB=CCRKL*CCRKL1*QC*US00
      AC=AB+0.5*CWMKL
      AD=AB*CCRKL1
      AE=CPR*TKL*TKL
      AF=AE*PP
      AG=AA*ELV
      AI=CP*AA
      COND=(AC-AD)*(AF*AQ-AI*AT+AE*QKL*AP)/(AC*(AF+AG))
!-----------------------------------------------------------------------
!-----------CHECK & CORRECT IF OVER CONDENSATION OCCURS-----------------
!-----------------------------------------------------------------------
      CONDK=(QKL-U00KL*QC*0.1)*RTWODT                             
!     CONDK=(QKL-U00KL*QC*0.6)*RTWODT                             
      IF(COND.GT.CONDK)THEN
        COND=CONDK                           
      ENDIF
!-----------------------------------------------------------------------
!----------CHECK & CORRECT IF SUPERSATUATION IS TOO HIGH----------------
!-----------------------------------------------------------------------
      QTEMP=QKL-AMAX1(0.,(COND-E0))*THH
      RQTMP=QTEMP/QC
      IF(RQTMP.GE.1.10)THEN
        COND=(QKL-1.10*QC)*RTWODT
      ENDIF
!-----------------------------------------------------------------------
      IF(COND.LT.0.)THEN                                           
        COND=0.
      ENDIF                                                             
!-----------------------------------------------------------------------
!-------------------UPDATE OF T, Q AND CWM------------------------------
!-----------------------------------------------------------------------
   80 CONTINUE                                                          
      CONE0=COND-E0
      CWM_T(L,I,J)=CONE0*THH+CWMKL                                          
!     
!-----------------------------------------------------------------------
!     ACCUMULATE LATENT HEATING DUE TO GRID-SCALE PRECIP/EVAP.
!     SCALE BY THE RECIPROCAL OF THE PERIOD AT WHICH THIS ROUTINE
!     IS CALLED.  THIS PERIOD IS THE PHYSICS TIMESTEP.
!-----------------------------------------------------------------------
!
      T_T(L,I,J)=ELV*RCP*CONE0*THH+TKL                                      
      Q_T(L,I,J)=-CONE0*THH+QKL                                             
      IF(CWM_T(L,I,J).LE.0.)CWM_T(L,I,J)=0.                            
!
   90 CONTINUE
!
  100 CONTINUE                                     
!-----------------------------------------------------------------------
!-------------------SAVE T, Q AND P FOR THIS STEP-----------------------
!-----------------------------------------------------------------------
!***
!***  TRANSPOSE BACK THE NEEDED ARRAYS
!***
!     CALL SGETMO(T_T,LM,LM,LDA,T,LDA)
!     CALL SGETMO(Q_T,LM,LM,LDA,Q,LDA)
!     CALL SGETMO(CWM_T,LM,LM,LDA,CWM,LDA)

      CALL SGETMO_T(T_T,T,                   &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )
      CALL SGETMO_T(Q_T,Q,                   &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )


      CALL SGETMO_T(CWM_T,CWM,               &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )
!
      DO J=jts,jte
      DO L=1,LM
      DO I=its,ite
        Q0(I,L,J)=Q(I,L,J)                                                   
        T0(I,L,J)=T(I,L,J)                                                   
      ENDDO
      ENDDO
      ENDDO
!
!     DO J=jts,jte
!     DO I=its,ite
!       P0(I,J)=PD(I,J)
!     ENDDO
!     ENDDO
!-----------------------------------------------------------------------
 END SUBROUTINE GSCOND                          

!=======================================================================
 SUBROUTINE SGETMO(A,ACOL_T, ids,ide, jds,jde, kds,kde,                &
                             ims,ime, jms,jme, kms,kme,                &
                             its,ite, jts,jte, kts,kte                 )
!-----------------------------------------------------------------------
 IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER,INTENT(IN) ::  ids,ide, jds,jde, kds,kde,                &
                             ims,ime, jms,jme, kms,kme,                &
                             its,ite, jts,jte, kts,kte

      REAL, INTENT(IN) ,DIMENSION(its:ite,kts:kte,jts:jte) :: A
      REAL, INTENT(OUT),DIMENSION(kts:kte,its:ite,jts:jte) :: ACOL_T
      INTEGER::i,j,k 
!------------------------------------------------------
      DO j=jts,jte
      DO k=kts,kte
      DO i=its,ite
        ACOL_T(K,I,J)=A(I,K,J)
      ENDDO
      ENDDO
      ENDDO
!
 END SUBROUTINE SGETMO


!=======================================================================
 SUBROUTINE SGETMO_T(ACOL_T, A, ids,ide, jds,jde, kds,kde,             &
                             ims,ime, jms,jme, kms,kme,                &
                             its,ite, jts,jte, kts,kte                 )
!-----------------------------------------------------------------------
 IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER,INTENT(IN) ::  ids,ide, jds,jde, kds,kde,                &
                             ims,ime, jms,jme, kms,kme,                &
                             its,ite, jts,jte, kts,kte

      REAL, INTENT(OUT), DIMENSION(its:ite,kts:kte,jts:jte) :: A
      REAL, INTENT(IN ), DIMENSION(kts:kte,its:ite,jts:jte) :: ACOL_T
      INTEGER::i,j,k
!------------------------------------------------------
      DO j=jts,jte
      DO k=kts,kte
      DO i=its,ite
        A(I,K,J)=ACOL_T(K,I,J)
      ENDDO
      ENDDO
      ENDDO
!
 END SUBROUTINE SGETMO_T

!----------------------------------------------------------------------
 SUBROUTINE PRECPD (DT,NPHS,Q,CWM,Q0,T0,T,XLAND,PFLIP,                &
                   P8WFLIP,PREC,ACPREC,                               &
                   ids,ide, jds,jde, kds,kde,                         &
                   ims,ime, jms,jme, kms,kme,                         &
                   its,ite, jts,jte, kts,kte                          )
!----------------------------------------------------------------------
 IMPLICIT NONE
!----------------------------------------------------------------------
       REAL,PARAMETER :: A2=17.2693882,A3=273.16,A4=35.86,  &
                         PQ0=379.90516,R=287.04,C0=0.15,    &
                         CP=1004.6,ELWV=2.50E6,ELIV=2.834E6,ROW=1.E3, &
                         G=9.8,EPSQ=2.E-12,ELIW=ELIV-ELWV         
!
       REAL,PARAMETER :: RCP=1./CP,RROW=1./ROW
!---------------------------------------------------------------------- 
       INTEGER, INTENT(IN)        :: ids,ide, jds,jde, kds,kde ,      &
                                     ims,ime, jms,jme, kms,kme ,      &
                                     its,ite, jts,jte, kts,kte
!----------------------------------------------------------------------
      INTEGER, INTENT(IN)                     :: NPHS
      REAL,    INTENT(IN)                     :: DT
      REAL, INTENT(IN), DIMENSION(ims:ime, jms:jme):: XLAND
      REAL, INTENT(INOUT), DIMENSION(ims:ime, jms:jme):: PREC,ACPREC

! ACPREC : Accumulation precipitation from time 0 for grid scale only
! PREC   : Accumulation precipitation within one time step for grid scale only

      REAL, INTENT(INOUT), DIMENSION(its:ite, kts:kte, jts:jte)::CWM, &
                                                             Q,Q0,T0,T
!     REAL, INTENT(IN), DIMENSION(its:ite, jts:jte):: PD

      REAL, DIMENSION((ite-its+1)*(jte-jts+1))   :: PRECSL1,PRECRL1
      INTEGER,DIMENSION((ite-its+1)*(jte-jts+1)) :: IWL1
      INTEGER,DIMENSION((ite-its+1)*(jte-jts+1)) :: IPREC,JPREC
      REAL,DIMENSION(kts:kte, its:ite, jts:jte):: T_T,T0_T,Q_T,Q0_T, &
                                                  CWM_T

      REAL, INTENT(IN), DIMENSION(its:ite, kts:kte, jts:jte)::PFLIP
      REAL, INTENT(IN), DIMENSION(its:ite, kts:kte+1, jts:jte)::P8WFLIP
      INTEGER, DIMENSION(its:ite, jts:jte):: LMH,SR
      REAL,    DIMENSION(its:ite, jts:jte):: U00
      REAL,    DIMENSION(2*kte)   :: UL

      INTEGER :: LM,I,J,L,MYJS,MYJE,MYIS,MYIE,N,IMJM_LOC
      INTEGER :: MYJS2,MYJE2,NPRE,LML,IWL,IWK
      REAL    :: DTPH,RDTPH,TWODT,RTWODT,KE,US,EPS,CCLIMIT,CLIMIT,CWS, &
                 CSM1,CRS1,CRS2,CR,MI0,AA2,AVRAIN,ARATIM,UTIM,TTEMP,   &
                 WFIX,WMIN
      REAL    :: PDSL,CONST,U00IJ,ULL,PRECRL,PRECSL,PRAUT,PSAUT, &
                 PRACW,CONDE,RCONDE,TT,QQ,WW,U00KL,WMINK,PRECRK,PRECSK
      REAL    :: ERS,ERR,PSACI,PSM,PSM1,PSM2,PPR,PPS,CPDR,PID,TK,QK,  &
                 TMT0,TMT15,AI,BI,QW,QI,QINT,FI,FIW,QC,RQ,CCR,RQKLL,CWMK
      REAL    :: EXPF,AA1,AMAXCM,CS,TMT0K,U00KLT,AMAXRQ,HHT,QTEMP,TMT0T,&
                 TMT15T,QWT,QIT,QINTT,QCT,RQT,RQTT,ERK,RPRS,ERRT,ERST, &
                 AMAXPS,PRECRS,PRECSS,APREC,TOTPPT
!-----------------------------------------------------------------------
!***********************************************************************
!--------------PREPARATORY CALCULATIONS---------------------------------

     
      MYJS=jts
      MYJE=jte
      LM = kte
      MYIS=its
      MYIE=ite
      MYJS2=jts
      MYJE2=jte
      IMJM_LOC = (ite-its+1)*(jte-jts+1)

      DTPH  = float(NPHS) * DT
      RDTPH = 1. / DTPH
      TWODT= DTPH
      RTWODT=1./TWODT                                                   
      KE=2.0E-5                                                         
      US=1.
      EPS=0.622E0                                                       
      CCLIMIT=1.0E-3                                                    
      CLIMIT=1.0E-20                                                    
      CWS=0.025                                                         
      CSM1=5.0000E-8                                                    
      CRS1=5.00000E-6                                                   
      CRS2=6.66600E-10                                                  
      CR=5.0E-4                                                         
      MI0=5.0E-4                                                        
      AA2=1.25E-3                                                       
!     AVRAIN=AVRAIN+1.
!     ARATIM=ARATIM+1

      do j = jts,jte
      do i = its,ite
         LMH(i,j) = kme-1
      enddo
      enddo

      do j = jts,jte
      do i = its,ite
         U00(I,J)=(2.-XLAND(i,j))*0.75+(XLAND(i,j) -1.)*0.80
      enddo
      enddo

      DO L=1,2*LM
      IF(L.GE.LM-10.AND.L.LE.LM)THEN
        UL(L)=0.1*FLOAT(L-LM+10)
      ELSE
        UL(L)=0.
      ENDIF
      ENDDO

!-------------------PADDING CLOUD MIXING RATIO IF TOO SMALL-------------
      DO 20 L=1,LM                                                      
      DO J=MYJS,MYJE
      DO I=MYIS,MYIE
!       CWM(I,L,J)=CWM(I,L,J)*HTM(I,L,J)*HBM2(I,J) 
!       CWM(I,L,J)=CWM(I,L,J) ! *HBM2(I,J)                              
        IF(CWM(I,L,J).LT.0.)CWM(I,L,J)=0.
!------------------PADDING SPECIFIC HUMIDITY IF TOO SMALL---------------
        IF(Q(I,L,J).LT.EPSQ)Q(I,L,J)=EPSQ ! *HTM(I,L,J)                       
      ENDDO
      ENDDO
   20 CONTINUE                                                          
!
      UTIM=1.
!-----------------------------------------------------------------------
      DO N=1,IMJM_LOC
        IWL1(N)=0
        PRECRL1(N)=0.
        PRECSL1(N)=0.
      ENDDO
!------------CHOOSE THE COLUMNS WHERE PREC CAN BE PRODUCED--------------
      NPRE=0
      DO 35 J=MYJS2,MYJE2
      DO 35 I=MYIS,MYIE
!
      DO L=2,LM
        TTEMP=0.025*(T(I,L,J)-273.16)
        WFIX=0.9814*EXP(0.01873*L)
        WMIN=0.1E-3*EXP(TTEMP)*WFIX
        IF(CWM(I,L,J).GT.WMIN)GO TO 33
      ENDDO
!
      GO TO 35
   33 NPRE=NPRE+1
      IPREC(NPRE)=I
      JPREC(NPRE)=J
   35 CONTINUE
!------------------------------------------------------------------------
!***
!***  TRANSPOSE ARRAYS
!***
      CALL SGETMO(T,T_T,                   &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )
      CALL SGETMO(Q,Q_T,                   &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )
      CALL SGETMO(CWM,CWM_T,               &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )

!     CALL SGETMO(T,LDA,LDA,LM,T_T,LM)
!     CALL SGETMO(Q,LDA,LDA,LM,Q_T,LM)
!     CALL SGETMO(CWM,LDA,LDA,LM,CWM_T,LM)
!     CALL SGETMO(HTM,LDA,LDA,LM,HTM_T,LM)
!     CALL SGETMO(TRAIN,LDA,LDA,LM,TRAIN_T,LM)
!-----------------------------------------------------------------------
!-----------------BEGINING OF PRECIPITATION CALCULATION-----------------
!-----------------------------------------------------------------------
!***
!***  LOOP OVER ALL POSSIBLE PRECIPITATION POINTS
!***
!!omp parallel do 
!!omp& private(aai,aetal,ai,amaxcm,amaxps,amaxrq,aprec,bi,ccr,conde,
!!omp&         const,cpdr,cs,cwmk,detal,erk,err,errt,ers,erst,expf,fi,
!!omp&         fiw,hbm2k,hh,htmk,i,iwl,j,lml,mi0,pdsl,pid,ppr,pps,
!!omp&         pracw,praut,precrk,precrl,precsk,precsl,precss,psaci,
!!omp&         psaut,psm,psm1,psm2,qc,qct,qi,qint,qintt,qit,qk,qq,
!!omp&         qtemp,qw,qwt,rconde,rprs,rq,rqkll,rqt,rqtt,tk,tmt0,
!!omp&         tmt0k,tmt0t,tmt15,tmt15t,totppt,tt,ttemp,u00ij,u00kl,
!!omp&         u00klt,ull,wfix,wmink,ww)
!
      DO 300 N=1,NPRE
!
      I=IPREC(N)
      J=JPREC(N)
!     HBM2K=HBM2(I,J)
!     PDSL=PD(I,J)                                              
!     CONST=PDSL/G*TWODT                                        
      LML=LM-LMH(I,J)
      U00IJ=U00(I,J)
!
      DO 180 L=2,LM  
!
!     DETAL=DETA(L)
      ULL=UL(L)
!     AETAL=AETA(L)
      WFIX=0.9814*EXP(0.01873*L)
!
      PRECRL=0.                                                      
      PRECSL=0.                                                      
      PRAUT=0.                                                      
      PSAUT=0.                                                      
      PRACW=0.                                                      
      PSACI=0.                                                      
      ERR  =0.                                                      
      ERS  =0.                                                      
      PSM  =0.                                                      
      PSM1 =0.                                                      
      PSM2 =0.                                                      
      PPR  =0.
      PPS  =0.
      CPDR =0.
!     HH   =0.
      PID  =0.
      IWL  =0.
      CONDE=0.
      RCONDE=0.
!-----------------------------------------------------------------------
      TT=T_T(L,I,J)
      QQ=Q_T(L,I,J)
      WW=CWM_T(L,I,J)
!     HTMK=HTM_T(L,I,J)
!-----------------------------------------------------------------------
      U00KL=U00IJ+UL(L+LML)*(0.95-U00IJ)*UTIM
      TTEMP=0.025*(TT-273.16)
      WMINK=0.1E-3*EXP(TTEMP)*WFIX
!-----------------------------------------------------------------------
!----------CHOOSE THE POINTS WHERE PRECIPITATION CAN BE PRODUCED--------
!-----------------------------------------------------------------------
      PRECRK=AMAX1(0.,PRECRL1(N))                                         
      PRECSK=AMAX1(0.,PRECSL1(N))                                         
!     HH=HTMK*HBM2K
      IF(WW.LT.WMINK.AND.(PRECRK+PRECSK).EQ.0.)THEN
        PID=0.
      ELSE
!       PID=HH
        PID=1
      ENDIF
!-----------------------------------------------------------------------
!-------------------QW, QI AND QINT-------------------------------------
!-----------------------------------------------------------------------
      IF(PID.EQ.1.)THEN
!       CONST=PDSL/G*TWODT                                        
!       CONDE=CONST*DETAL                                               
        CONDE=(P8WFLIP(I,L+1,J)-P8WFLIP(I,L,J))/G*TWODT
        RCONDE=1./CONDE                                                   
        TK=TT
        QK=QQ
        TMT0=(TK-273.16) ! *HH
        TMT15=AMIN1(TMT0,-15.) ! *HH
        AI=0.008855
        BI=1.
!
        IF(TMT0.LT.-20.)THEN
          AI=0.007225
          BI=0.9674
        ENDIF
!
!       QW=HH*PQ0/(PDSL*AETAL+PT)  &
!                    *EXP(HH*A2*(TK-A3)/(TK-A4))                  
        QW=PQ0/PFLIP(I,L,J)  &
                     *EXP(A2*(TK-A3)/(TK-A4))                  
        QI=QW*(BI+AI*AMIN1(TMT0,0.))                                 
        QINT=QW*(1.-0.00032*TMT15*(TMT15+15.))                          
        IF(TMT0.LE.-40.)QINT=QI
!-------------------ICE-WATER ID NUMBER IW------------------------------
        IF(TMT0.LT.-15.)THEN                                             
          FI=QK-U00KL*QI
          IF(FI.GT.0..OR.WW.GT.CLIMIT) THEN                    
            IWL=1                                                   
          ELSE                                                           
            IWL=0                                                   
          ENDIF                                                         
        ENDIF                                                            
!
        IF(TMT0.LT.0.0.AND.TMT0.GE.-15.0)THEN                            
          IWL=0                                                   
          IF(IWL1(N).EQ.1.AND.WW.GT.CLIMIT)IWL=1 
        ENDIF                                                            
!
        IF(TMT0.GE.0.)THEN                                               
          IWL=0                                                      
        ENDIF                                                            
!----------------THE SATUATION SPECIFIC HUMIDITY------------------------
        FIW=FLOAT(IWL)                                                
        QC=(1.-FIW)*QINT+FIW*QI
!----------------THE RELATIVE HUMIDITY----------------------------------
        IF(QC.LE.0.)THEN                                                
          RQ=1.E-10
        ELSE                                                             
          RQ=QK/QC
        ENDIF                                                             
!----------------CLOUD COVER RATIO CCR----------------------------------
        IF(RQ.LE.U00KL)THEN                                      
          CCR=0.
        ELSE                                                            
          RQKLL=AMIN1(US,RQ)                                         
          CCR=1.-SQRT((US-RQKLL)/(US-U00KL))                     
        ENDIF                                                             
!-----------CORRECT CCR IF IT IS TOO SMALL IN LARGE CWM REGIONS--------
        IF(CCR.GE.0.01.AND.CCR.LE.0.2.AND. WW.GE.0.2E-3)THEN
          CCR=AMIN1(1.,WW*1.0E3)
        ENDIF
      ENDIF
   60              CONTINUE
!-----------------------------------------------------------------------
!------------------PRECIPITATION PRODUCTION RATES-----------------------
!------------------AUTO-CONVERT RATES-----------------------------------
!-----------------------------------------------------------------------
      IF(PID.EQ.1.)THEN
        IWK=IWL
        CWMK=AMAX1(0.,WW-CLIMIT)                                 
        MI0=WMINK
!
        IF(IWK.EQ.1)THEN                                                 
          EXPF=EXP(0.025*TMT0)                                           
          AA1=1.E-3*EXPF                                                 
          PSAUT=AA1*AMAX1(0.,CWMK-MI0)                                
          CPDR=-PSAUT*TWODT                                           
          IF(-CPDR.GE.CWMK)THEN                                         
            CPDR=-CWMK                                                 
            PSAUT=-CPDR*RTWODT                                      
          ENDIF                                                         
        ELSE                                                              
          AMAXCM=AMAX1(0.,CWMK-MI0)                                      
          PRAUT=C0*AMAXCM*AMAXCM                                   
          CPDR=-PRAUT*TWODT                                           
          IF(-CPDR.GE.CWMK)THEN                                         
            CPDR=-CWMK                                                  
            PRAUT=-CPDR*RTWODT                                       
          ENDIF                                                         
        ENDIF                                                            
        PPR=PRAUT*CONDE
        PPS=PSAUT*CONDE
      ENDIF
!
      IF(PID.EQ.1.)THEN
!       WW=CPDR*HH+WW                                   
        WW=CPDR+WW                                   
        PRECRL=PRECRL1(N)+PPR ! *HH
        PRECSL=PRECSL1(N)+PPS ! *HH
      ENDIF
!-----------------------------------------------------------------------
!-----------------------ACCRETIONS--------------------------------------
!-----------------------------------------------------------------------
      IF(PID.EQ.1.)THEN
        IWK=IWL
        CWMK=WW
        PRECRK=AMAX1(0.,PRECRL1(N))                                         
        PRECSK=AMAX1(0.,PRECSL1(N))                                         
        IF(IWK.EQ.1)THEN                                                 
          EXPF=EXP(0.025*TMT0)                                           
          CS=AA2*EXPF                                                    
          PSACI=CS*AMAX1(0.,CWMK)*PRECSK                          
          CPDR=-PSACI*TWODT                                           
          IF(-CPDR.GE.CWMK)THEN                                         
            CPDR=-CWMK                                                 
            PSACI=-CPDR*RTWODT                                      
          ENDIF                                                         
        ELSE                                                              
          PRACW=CR*AMAX1(0.,CWMK)*(PRECRK+PRECSK)             
          CPDR=-PRACW*TWODT                                           
          IF(-CPDR.GE.CWMK)THEN                                         
            CPDR=-CWMK                                                  
            PRACW=-CPDR*RTWODT                                       
          ENDIF                                                         
        ENDIF                                                            
        PPR=PRACW*CONDE
        PPS=PSACI*CONDE
      ENDIF
!
      IF(PID.EQ.1.)THEN
!       WW=CPDR*HH+WW                                   
        WW=CPDR+WW                                   
        PRECRL=PRECRL+PPR ! *HH                                       
        PRECSL=PRECSL+PPS ! *HH                                       
      ENDIF
!-----------------------------------------------------------------------
!-----EVAPORATION/CONDENSATION OF PRECIPITATION-------------------------
!***** ERR & ERS POSITIVE--EVAPORATION                                  
!***** ERR & ERS NEGTIVE---CONDENSATION                                 
!-----------------------------------------------------------------------
      IF(PID.EQ.1.0)THEN
        QK=QQ
        TMT0K=TMT0
        IF(TMT0K.LT.-30.)TMT0K=-30.                                        
        PRECRK=AMAX1(0.,PRECRL)                                         
        PRECSK=AMAX1(0.,PRECSL)                                         
!---------------------------------------------------------------------- 
! INCREASE THE EVAPORATION/CONDENSATION FOR STRONG/LIGHT PREC           
!---------------------------------------------------------------------- 
        U00KLT=U00KL
        AMAXRQ=AMAX1(0.,U00KL-RQ)   
        ERR=KE*AMAXRQ*PRECRK**0.5                                      
!
        IF(TMT0.GE.0.)THEN                                               
          ERS=0.                                                         
        ELSE                                                              
          ERS=(CRS1+CRS2*TMT0K)*AMAXRQ*PRECSK/U00KLT                      
        ENDIF                                                            
!
        IF(ERR+ERS.LE.1.E-20) GO TO 125
!---------------CORRECT IF OVER-EVAPO./COND. OCCURS-------------------- 
!       HHT=HH*TWODT
        HHT=TWODT
        TTEMP=TT-RCP*(ELWV*ERR+ELIV*ERS)*HHT
        QTEMP=QQ+HHT*(ERR+ERS)
        TMT0T=(TTEMP-273.16)  ! *HH
        IF(TMT0T.LT.-30.)TMT0T=-30.                                        
        TMT15T=AMIN1(TMT0T,-15.)  ! *HH
        AI=0.008855
        BI=1.
!
        IF(TMT0T.LT.-20.)THEN
          AI=0.007225
          BI=0.9674
        ENDIF
!
!       QWT=HH*PQ0/(PDSL*AETAL+PT)*EXP(HH*A2*(TTEMP-A3)/(TTEMP-A4))                  
        QWT=PQ0/PFLIP(I,L,J)*EXP(A2*(TTEMP-A3)/(TTEMP-A4))                  
        QIT=QWT*(BI+AI*AMIN1(TMT0T,0.))                                 
        QINTT=QWT*(1.-0.00032*TMT15T*(TMT15T+15.))                          
        IF(TMT0T.LE.-40.)QINTT=QIT                                    
        FIW=FLOAT(IWL)
        QCT=(1.-FIW)*QINTT+FIW*QIT
!
        IF(QCT.LE.1.E-10) THEN
          RQT=1.E-10
          RQTT=1.E-10
        ELSE
          RQT=QTEMP/QCT
          RQTT=QQ/QCT
        ENDIF
!
        IF(RQT.LE.U00KL) GO TO 125
!
        ERK=(U00KL-RQTT)*QCT*RTWODT
        RPRS=ERK/(PRECRK+PRECSK)                                          
        ERRT=PRECRK*RPRS                                                
        ERST=PRECSK*RPRS                                                
        ERR=AMAX1(0.,0.5*(ERR+ERRT))
        ERS=AMAX1(0.,0.5*(ERS+ERST))
!         
 125    CONTINUE
!
        PPR=-ERR*CONDE                                                 
        PPS=-ERS*CONDE                                                 
!
        IF(-PPR.GE.PRECRK)THEN                                         
          PPR=-PRECRK                                                  
          ERR=-PPR*RCONDE
        ENDIF                                                            
!
        IF(-PPS.GE.PRECSK)THEN                                         
          PPS=-PRECSK                                                  
          ERS=-PPS*RCONDE
        ENDIF                                                            
!
      ENDIF                                                            
!
      IF(PID.EQ.1.)THEN
        PRECRL=PRECRL+PPR ! *HH
        PRECSL=PRECSL+PPS ! *HH
      ENDIF                                                            
!-----------------------------------------------------------------------
!--------------------MELTING OF THE SNOW--------------------------------
!-----------------------------------------------------------------------
      IF(PID.EQ.1.)THEN
        CWMK=WW                                                     
        AMAXPS=AMAX1(0.,PRECSL)                                      
!
        IF(TMT0.GT.0.)THEN                                               
          PSM1=CSM1*TMT0*TMT0*AMAXPS                                  
          PSM2=CWS*CR*CWMK*AMAXPS                                     
          PSM=PSM1+PSM2                                         
        ELSE                                                              
          PSM1=0.
          PSM2=0.
          PSM=0.
        ENDIF                                                            
!
        PPR=PSM*CONDE
        PPS=-PSM*CONDE
!
        IF(-PPS.GE.AMAXPS)THEN                                         
          PPS=-AMAXPS                                                  
          PPR=AMAXPS                                                   
          PSM1=-PPS*RCONDE                                            
          PSM2=0.
          PSM=PSM1
        ENDIF                                                            
!
      ENDIF                                                            
!
      IF(PID.EQ.1.)THEN
        PRECRL=PRECRL+PPR ! *HH                                       
        PRECSL=PRECSL+PPS ! *HH                                       
      ENDIF                                                            
!-----------------------------------------------------------------------
!---------------UPDATE T AND Q------------------------------------------
!-----------------------------------------------------------------------
      IF(PID.EQ.1.)THEN
!       HHT=HH*TWODT                                        
        HHT=TWODT                                        
        TT=-RCP*(ELWV*ERR+ELIV*ERS+ELIW*PSM1)*HHT+TT
        QQ=(ERR+ERS)*HHT+QQ
      ENDIF                                                            
!
!     IF(HH.EQ.1.)THEN
        IWL1(N)=IWL
        PRECRL1(N)=PRECRL
        PRECSL1(N)=PRECSL
!     ENDIF
!     
!     ACCUMULATE LATENT HEATING DUE TO GRID-SCALE PRECIP/EVAP.
!     SCALE BY THE RECIPROCAL OF THE PERIOD AT WHICH THIS ROUTINE
!     IS CALLED.  THIS PERIOD IS THE PHYSICS TIMESTEP.
!
!     TRAIN_T(L,I,J)=TRAIN_T(L,I,J)+(TT-T_T(L,I,J))*RDTPH
      T_T(L,I,J)=TT
      Q_T(L,I,J)=QQ
      CWM_T(L,I,J)=WW
  180                      CONTINUE
!-----------------------------------------------------------------------
!-------------------THE PRECIPITATION ON SFC----------------------------
!-----------------------------------------------------------------------
      PRECRS=PRECRL1(N)*RROW                                         
      PRECSS=PRECSL1(N)*RROW                                         
!
!     APREC=PRECRS+PRECSS

! from m to mm

      APREC=(PRECRS+PRECSS)*1000.
!     PREC(I,J)=PREC(I,J)+PRECRS+PRECSS                              
      PREC(I,J)=APREC
      ACPREC(I,J)=ACPREC(I,J)+APREC
!-----------------------------------------------------------------------
!---------------THE SNOW AND RAIN RATIO OF SFC PREC---------------------
!----SR IS THE RATIO OF SNOW TO THE TOTAL PRECIP------------------------
!----IF TOTAL PRECIP IS ZERO, SR IS ZERO--------------------------------
!-----------------------------------------------------------------------
      TOTPPT=PRECRS+PRECSS
      IF (TOTPPT.GT.1.E-8) THEN
       SR(I,J)=PRECSS/TOTPPT
      ELSE
       SR(I,J)=0.
      ENDIF
!-----------------------------------------------------------------------
  300 CONTINUE
!-----------------------------------------------------------------------
!***
!***  TRANSPOSE BACK
!***

      CALL SGETMO_T(T_T,T,                   &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )
      CALL SGETMO_T(Q_T,Q,                   &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )


      CALL SGETMO_T(CWM_T,CWM,               &
                  ids,ide,jds,jde,kds,kde, &
                  ims,ime,jms,jme,kms,kme, &
                  its,ite,jts,jte,kts,kte  )

!     CALL SGETMO(T_T,LM,LM,LDA,T,LDA)
!     CALL SGETMO(Q_T,LM,LM,LDA,Q,LDA)
!     CALL SGETMO(CWM_T,LM,LM,LDA,CWM,LDA)
!     CALL SGETMO(TRAIN_T,LM,LM,LDA,TRAIN,LDA)
!-----------------------------------------------------------------------

 END SUBROUTINE PRECPD 

END MODULE module_mp_eta

      SUBROUTINE GTD6(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
C        Neutral Atmosphere Empirical Model from the surface to lower
C          exosphere  MSISE90 (JGR, 96, 1159-1172, 1991)
C         A.E.Hedin 4/24/90
C           See subroutine GHP6 to specify a pressure rather than
C           altitude.
C     INPUT:
C        IYD - YEAR AND DAY AS YYDDD or DDD (day of year from 1 to 365)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C     Note:  Ut, Local Time, and Longitude are used independently in the
C            model and are not of equal importance for every situation.  
C            For the most physically realistic calculation these three
C            variables should be consistent (STL=SEC/3600+GLONG/15).
C            F107, F107A, and AP effects are not large below 80 km 
C            and these can be set to 150., 150., and 4. respectively.
C     OUTPUT:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)                       
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
C      TO GET OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.) 
C
C      O, H, and N set to zero below 72.5 km
C      Exospheric temperature set to average for altitudes below 120 km.
C        
C           The following is for test and special purposes:
C            TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW)
C               WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
C               FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C               FOR THE FOLLOWING VARIATIONS
C               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
C               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
C               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
C               7 - DIURNAL               8 - SEMIDIURNAL
C               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
C              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
C              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
C              15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
C              16 - ALL TINF VAR         17 - ALL TLB VAR
C              18 - ALL TN1 VAR           19 - ALL S VAR
C              20 - ALL TN2 VAR           21 - ALL NLB VAR
C              22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR
C
C              To get current values of SW: CALL TRETRV(SW)
C
      DIMENSION D(8),T(2),AP(1),D6(8),T6(2)
      DIMENSION ZN3(5),ZN2(4)
      COMMON/GTS3C/TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
     & ,G0,RL,DD,DB14,TR12
      COMMON/MESO6/TN1(5),TN2(4),TN3(5),TGN1(2),TGN2(2),TGN3(2)
      COMMON/LOWER6/PTM(10),PDM(10,8)
      COMMON/PARM6/PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),
     $ PMA(100,10)
      COMMON/DATIM6/ISD(3),IST(2),NAM(2)
      COMMON/DATIME/ISDATE(3),ISTIME(2),NAME(2)
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/MAVG6/PAVGM(10)
      COMMON/DMIX/DM04,DM16,DM28,DM32,DM40,DM01,DM14
      COMMON/PARMB/GSURF,RE
      COMMON/METSEL/IMR
      EXTERNAL GTD6BK
      DATA MN3/5/,ZN3/32.5,20.,15.,10.,0./
      DATA MN2/4/,ZN2/72.5,55.,45.,32.5/
      DATA ZMIX/62.5/,ALAST/99999./,MSSL/-999/

c	write(6,*) 'gtd6'
c	write(6,*) IYD
c	write(6,*) MASS
c	write(6,*) ALT
c	write(6,*) GLAT
c	write(6,*) GLONG
c	write(6,*) STL
c	write(6,*) F107A
c	write(6,*) F107
c	write(6,*) AP
c	write(6,*) SEC
	
C      Put identification data into common/datime/
      DO 1 I=1,3
        ISDATE(I)=ISD(I)
    1 CONTINUE
      DO 2 I=1,2
        ISTIME(I)=IST(I)
        NAME(I)=NAM(I)
    2 CONTINUE
C
Ce        Test for changed input
      V1=VTST(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,1)
C       Latitude variation of gravity (none for SW(2)=0)
      XLAT=GLAT
      IF(SW(2).EQ.0) XLAT=45.
      CALL GLATF(XLAT,GSURF,RE)
C
      XMM=PDM(5,3)
C
C       THERMOSPHERE/UPPER MESOSPHERE [above ZN2(1)]
      ALTT=AMAX1(ALT,ZN2(1))
      MSS=MASS
Ce       Only calculate N2 in thermosphere if alt in mixed region
      IF(ALT.LT.ZMIX.AND.MASS.GT.0) MSS=28
Ce       Only calculate thermosphere if input parameters changed
Ce         or altitude above ZN2(1) in mesosphere
      IF(V1.EQ.1..OR.ALT.GT.ZN2(1).OR.ALAST.GT.ZN2(1).OR.MSS.NE.MSSL)
     $ THEN
        CALL GTS6(IYD,SEC,ALTT,GLAT,GLONG,STL,F107A,F107,AP,MSS,D6,T6)
        DM28M=DM28
C         metric adjustment
        IF(IMR.EQ.1) DM28M=DM28*1.E6
        MSSL=MSS
      ENDIF
      T(1)=T6(1)
      T(2)=T6(2)
      IF(ALT.GE.ZN2(1)) THEN
        DO 5 J=1,8
          D(J)=D6(J)
    5   CONTINUE
        GOTO 10
      ENDIF
C
C       LOWER MESOSPHERE/UPPER STRATOSPHERE [between ZN3(1) and ZN2(1)]
C         Temperature at nodes and gradients at end nodes
C         Inverse temperature a linear function of spherical harmonics
Ce        Only calculate nodes if input changed
       IF(V1.EQ.1..OR.ALAST.GE.ZN2(1)) THEN
        TGN2(1)=TGN1(2)
        TN2(1)=TN1(5)
        TN2(2)=PMA(1,1)*PAVGM(1)/(1.-SW(20)*GLOB6S(PMA(1,1)))
        TN2(3)=PMA(1,2)*PAVGM(2)/(1.-SW(20)*GLOB6S(PMA(1,2)))
        TN2(4)=PMA(1,3)*PAVGM(3)/(1.-SW(20)*SW(22)*GLOB6S(PMA(1,3)))
        TGN2(2)=PAVGM(9)*PMA(1,10)*(1.+SW(20)*SW(22)*GLOB6S(PMA(1,10)))
     $  *TN2(4)*TN2(4)/(PMA(1,3)*PAVGM(3))**2
        TN3(1)=TN2(4)
       ENDIF
       IF(ALT.GE.ZN3(1)) GOTO 6
C
C       LOWER STRATOSPHERE AND TROPOSPHERE [below ZN3(1)]
C         Temperature at nodes and gradients at end nodes
C         Inverse temperature a linear function of spherical harmonics
Ce        Only calculate nodes if input changed
        IF(V1.EQ.1..OR.ALAST.GE.ZN3(1)) THEN
         TGN3(1)=TGN2(2)
         TN3(2)=PMA(1,4)*PAVGM(4)/(1.-SW(22)*GLOB6S(PMA(1,4)))
         TN3(3)=PMA(1,5)*PAVGM(5)/(1.-SW(22)*GLOB6S(PMA(1,5)))
         TN3(4)=PMA(1,6)*PAVGM(6)/(1.-SW(22)*GLOB6S(PMA(1,6)))
         TN3(5)=PMA(1,7)*PAVGM(7)/(1.-SW(22)*GLOB6S(PMA(1,7)))
         TGN3(2)=PMA(1,8)*PAVGM(8)*(1.+SW(22)*GLOB6S(PMA(1,8)))
     $   *TN3(5)*TN3(5)/(PMA(1,7)*PAVGM(7))**2
        ENDIF
    6   CONTINUE
        IF(MASS.EQ.0) GOTO 50
Ce          Linear transition to full mixing at ZMIX from almost
Ce            full mixing at ZN2(1) to improve efficiency
        DMC=0
        IF(ALT.GT.ZMIX) DMC=1.-(ZN2(1)-ALT)/(ZN2(1)-ZMIX)
        DZ28=D6(3)
C      ***** N2 DENSITY ****
        DMR=D6(3)/DM28M-1.
        D(3)=DENSM(ALT,DM28M,XMM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
        D(3)=D(3)*(1.+DMR*DMC)
C      ***** HE DENSITY ****
        D(1)=0
        IF(MASS.NE.4.AND.MASS.NE.48) GOTO 204
          DMR=D6(1)/(DZ28*PDM(2,1))-1.
          D(1)=D(3)*PDM(2,1)*(1.+DMR*DMC)
  204   CONTINUE
C      **** O DENSITY ****
        D(2)=0
  216   CONTINUE
C      ***** O2 DENSITY ****
        D(4)=0
        IF(MASS.NE.32.AND.MASS.NE.48) GOTO 232
          DMR=D6(4)/(DZ28*PDM(2,4))-1.
          D(4)=D(3)*PDM(2,4)*(1.+DMR*DMC)
  232   CONTINUE
C      ***** AR DENSITY ****
        D(5)=0
        IF(MASS.NE.40.AND.MASS.NE.48) GOTO 240
          DMR=D6(5)/(DZ28*PDM(2,5))-1.
          D(5)=D(3)*PDM(2,5)*(1.+DMR*DMC)
  240   CONTINUE
C      ***** HYDROGEN DENSITY ****
        D(7)=0
C      ***** ATOMIC NITROGEN DENSITY ****
        D(8)=0
C
C       TOTAL MASS DENSITY
C
        IF(MASS.EQ.48)
     $   D(6) = 1.66E-24*(4.*D(1)+16.*D(2)+28.*D(3)+32.*D(4)+40.*D(5)+
     &       D(7)+14.*D(8))  
        T(2)=TZ
   10 CONTINUE
      GOTO 90
   50 CONTINUE
      DD=DENSM(ALT,1.,0,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)                
      T(2)=TZ
   90 CONTINUE
      ALAST=ALT
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GHP6(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,
     $  D,T,PRESS)
C       FIND ALTITUDE OF PRESSURE SURFACE (PRESS) FROM GTD6
C     INPUT:
C        IYD - YEAR AND DAY AS YYDDD
C        SEC - UT(SEC)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
C                    TO CURRENT TIME
C        PRESS - PRESSURE LEVEL(MB)
C     OUTPUT:
C        ALT - ALTITUDE(KM) 
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
      COMMON/PARMB/GSURF,RE
      COMMON/METSEL/IMR
      DIMENSION D(8),T(2)
      DATA BM/1.3806E-19/,RGAS/831.4/
      DATA TEST/.00043/
      PL=ALOG10(PRESS)
C      Initial altitude estimate
      IF(PL.GE.-5.) THEN
         IF(PL.GT.2.5) ZI=18.06*(3.00-PL)
         IF(PL.GT..75.AND.PL.LE.2.5) ZI=14.98*(3.08-PL)
         IF(PL.GT.-1..AND.PL.LE..75) ZI=17.8*(2.72-PL)
         IF(PL.GT.-2..AND.PL.LE.-1.) ZI=14.28*(3.64-PL)
         IF(PL.GT.-4..AND.PL.LE.-2.) ZI=12.72*(4.32-PL)
         IF(PL.LE.-4.) ZI=25.3*(.11-PL)
         IDAY=MOD(IYD,1000)
         CL=GLAT/90.
         CL2=CL*CL
         IF(IDAY.LT.182) CD=1.-IDAY/91.25
         IF(IDAY.GE.182) CD=IDAY/91.25-3.
         CA=0
         IF(PL.GT.-1.11.AND.PL.LE.-.23) CA=1.0
         IF(PL.GT.-.23) CA=(2.79-PL)/(2.79+.23)
         IF(PL.LE.-1.11.AND.PL.GT.-3.) CA=(-2.93-PL)/(-2.93+1.11)
         Z=ZI-4.87*CL*CD*CA-1.64*CL2*CA+.31*CA*CL
      ENDIF
      IF(PL.LT.-5.) Z=22.*(PL+4.)**2+110
      L=0
C      ITERATION LOOP
   10 CONTINUE
        L=L+1
        CALL GTD6(IYD,SEC,Z,GLAT,GLONG,STL,F107A,F107,AP,48,D,T)
        XN=D(1)+D(2)+D(3)+D(4)+D(5)+D(7)+D(8)
        P=BM*XN*T(2)
        IF(IMR.EQ.1) P=P*1.E-6
        DIFF=PL-ALOG10(P)
        IF(ABS(DIFF).LT.TEST .OR. L.EQ.6) GOTO 20
        XM=D(6)/XN/1.66E-24
        G=GSURF/(1.+Z/RE)**2
        SH=RGAS*T(2)/(XM*G)
C         New altitude estimate using scale height
        Z=Z-SH*DIFF*2.302
        GOTO 10
   20 CONTINUE
      IF(L.EQ.6) WRITE(6,100) PRESS,DIFF
  100 FORMAT(1X,29HGHP6 NOT CONVERGING FOR PRESS,1P2E12.2)
      ALT=Z
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GLATF(LAT,GV,REFF)
C      CALCULATE LATITUDE VARIABLE GRAVITY (GV) AND EFFECTIVE
C      RADIUS (REFF)
      REAL LAT,LATL
      DATA DGTR/1.74533E-2/,LATL/-999./
      IF(LAT.NE.LATL) C2 = COS(2.*DGTR*LAT)
      LATL=LAT
      GV = 980.616*(1.-.0026373*C2)
      REFF = 2.*GV/(3.085462E-6 + 2.27E-9*C2)*1.E-5
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION VTST(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,IC)
C       Test if geophysical variables or switches changed and save
C       Return 0 if unchanged and 1 if changed
      DIMENSION AP(7),IYDL(2),SECL(2),GLATL(2),GLL(2),STLL(2)
      DIMENSION FAL(2),FL(2),APL(7,2),SWL(25,2),SWCL(25,2)
      COMMON/CSW/SW(25),ISW,SWC(25)
      DATA IYDL/2*-999/,SECL/2*-999./,GLATL/2*-999./,GLL/2*-999./
      DATA STLL/2*-999./,FAL/2*-999./,FL/2*-999./,APL/14*-999./
      DATA SWL/50*-999./,SWCL/50*-999./
      VTST=0
      IF(IYD.NE.IYDL(IC)) GOTO 10
      IF(SEC.NE.SECL(IC)) GOTO 10
      IF(GLAT.NE.GLATL(IC)) GOTO 10
      IF(GLONG.NE.GLL(IC)) GOTO 10
      IF(STL.NE.STLL(IC)) GOTO 10
      IF(F107A.NE.FAL(IC)) GOTO 10
      IF(F107.NE.FL(IC)) GOTO 10
      j = 7
      if (SW(9) .ne. -1.) j = 1
      DO 5 I=1,j
        IF(AP(I).NE.APL(I,IC)) GOTO 10
    5 CONTINUE
      DO 7 I=1,25
        IF(SW(I).NE.SWL(I,IC)) GOTO 10
        IF(SWC(I).NE.SWCL(I,IC)) GOTO 10
    7 CONTINUE
      GOTO 20
   10 CONTINUE
      VTST=1
      IYDL(IC)=IYD
      SECL(IC)=SEC
      GLATL(IC)=GLAT
      GLL(IC)=GLONG
      STLL(IC)=STL
      FAL(IC)=F107A
      FL(IC)=F107
      DO 15 I=1,j
        APL(I,IC)=AP(I)
   15 CONTINUE
      DO 16 I=1,25
        SWL(I,IC)=SW(I)
        SWCL(I,IC)=SWC(I)
   16 CONTINUE
   20 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GTS6(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
C        Neutral Thermosphere Model above 72.5 km for MSISE-90
C         A.E.Hedin 3/9/90
C         Coefficients not changed for 120km and above, but results may differ
C        by a few percent from MSIS-86 (GTS5) with introduction of a
C        latitude dependent accel. of gravity.
C         Lower thermosphere reformulated for better continuation into
C        lower atmosphere.
C        For efficiency:
C         Exospheric temperature left at average value for alt below 120km;
C         120 km gradient left at average value for alt below 72 km;
C     INPUT:
C        IYD - YEAR AND DAY AS YYDDD
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM) (GREATER THAN 72.5 KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C     Note:  Ut, Local Time, and Longitude are used independently in the
C            model and are not of equal importance for every situation.  
C            For the most physically realistic calculation these three
C            variables should be consistent (STL=SEC/3600+GLONG/15).
C     OUTPUT:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
C           The following is for test and special purposes:
C           (1) LOWER BOUND QUANTITIES IN COMMON/GTS3C/
C           (2) TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW)
C               WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
C               FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C               FOR THE FOLLOWING VARIATIONS
C               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
C               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
C               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
C               7 - DIURNAL               8 - SEMIDIURNAL
C               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
C              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
C              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
C              15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
C              16 - ALL TINF VAR         17 - ALL TLB VAR
C              18 - ALL TN1 VAR           19 - ALL S VAR
C              20 - ALL TN2 VAR           21 - ALL NLB VAR
C              22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR
C
C              To get current values of SW: CALL TRETRV(SW)
C
      LOGICAL METER
      DIMENSION ZN1(5)
      COMMON/GTS3C/TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
     & ,G0,RL,DD,DB14,TR12
      COMMON/MESO6/TN1(5),TN2(4),TN3(5),TGN1(2),TGN2(2),TGN3(2)
      DIMENSION D(8),T(2),MT(10),AP(1),ALTL(8)
      COMMON/LOWER6/PTM(10),PDM(10,8)
      COMMON/PARM6/PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),
     $ PMA(100,10)
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/TTEST/TINFG,GB,ROUT,TT(15)
      COMMON/DMIX/DM04,DM16,DM28,DM32,DM40,DM01,DM14
      COMMON/METSEL/IMR
      DATA MT/48,0,4,16,28,32,40,1,49,14/
      DATA ALTL/200.,400.,160.,200.,240.,450.,320.,450./
      DATA MN1/5/,ZN1/120.,110.,100.,90.,72.5/
      DATA DGTR/1.74533E-2/,DR/1.72142E-2/,ALAST/-999./
Ce        Test for changed input
      V2=VTST(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,2)
C
      YRD=IYD
      ZA=PDL(16,2)
      ZN1(1)=ZA
      DO 2 J=1,8
        D(J)=0.
    2 CONTINUE
Ce       TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1)
      IF(ALT.GT.ZN1(1)) THEN
        IF(V2.EQ.1..OR.ALAST.LE.ZN1(1)) TINF=PTM(1)*PT(1)
     $  *(1.+SW(16)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PT))
      ELSE
        TINF=PTM(1)*PT(1)
      ENDIF
      T(1)=TINF
Ce         GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
      IF(ALT.GT.ZN1(5)) THEN
        IF(V2.EQ.1.OR.ALAST.LE.ZN1(5)) G0=PTM(4)*PS(1)
     $   *(1.+SW(19)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PS))
      ELSE
        G0=PTM(4)*PS(1)
      ENDIF
Ce      Calculate these temperatures only if input changed
      IF(V2.EQ.1.)
     $  TLB=PTM(2)*(1.+SW(17)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,
     $  F107A,F107,AP,PD(1,4)))*PD(1,4)
       S=G0/(TINF-TLB)
Ce       Lower thermosphere temp variations not significant for
Ce        density above 300 km
       IF(ALT.LT.300.) THEN
        IF(V2.EQ.1..OR.ALAST.GE.300.) THEN
         TN1(2)=PTM(7)*PTL(1,1)/(1.-SW(18)*GLOB6S(PTL(1,1)))
         TN1(3)=PTM(3)*PTL(1,2)/(1.-SW(18)*GLOB6S(PTL(1,2)))
         TN1(4)=PTM(8)*PTL(1,3)/(1.-SW(18)*GLOB6S(PTL(1,3)))
         TN1(5)=PTM(5)*PTL(1,4)/(1.-SW(18)*SW(20)*GLOB6S(PTL(1,4)))
         TGN1(2)=PTM(9)*PMA(1,9)*(1.+SW(18)*SW(20)*GLOB6S(PMA(1,9)))
     $   *TN1(5)*TN1(5)/(PTM(5)*PTL(1,4))**2
        ENDIF
       ELSE
        TN1(2)=PTM(7)*PTL(1,1)
        TN1(3)=PTM(3)*PTL(1,2)
        TN1(4)=PTM(8)*PTL(1,3)
        TN1(5)=PTM(5)*PTL(1,4)
        TGN1(2)=PTM(9)*PMA(1,9)
     $  *TN1(5)*TN1(5)/(PTM(5)*PTL(1,4))**2
       ENDIF
C
      Z0=ZN1(4)
      T0=TN1(4)
      ZLB=PTM(6)
      TR12=1.
C
      IF(MASS.EQ.0) GO TO 50
C       N2 variation factor at Zlb
      G28=SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107, 
     & AP,PD(1,3))
C        Variation of turbopause height
      DAY=AMOD(YRD,1000.)
      ZHF=PDL(25,2)
     $    *(1.+SW(5)*PDL(25,1)*SIN(DGTR*GLAT)*COS(DR*(DAY-PT(14))))
C
      YRD=IYD
      T(1)=TINF
      XMM=PDM(5,3)
C
      DO 10 J = 1,10
      IF(MASS.EQ.MT(J))   GO TO 15
   10 CONTINUE
      WRITE(6,100) MASS
      GO TO 90
   15 IF(ALT.GT.ALTL(6).AND.MASS.NE.28.AND.MASS.NE.48) GO TO 17
C
C       **** N2 DENSITY ****
C
C      Diffusive density at Zlb
      DB28 = PDM(1,3)*EXP(G28)*PD(1,3)
C      Diffusive density at Alt
      D(3)=DENSU(ALT,DB28,TINF,TLB, 28.,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
      DD=D(3)
C      Turbopause
      ZH28=PDM(3,3)*ZHF
      ZHM28=PDM(4,3)*PDL(6,2) 
      XMD=28.-XMM
C      Mixed density at Zlb
      B28=DENSU(ZH28,DB28,TINF,TLB,XMD,-1.,TZ,ZLB,S,MN1,ZN1,TN1,TGN1)
      IF(ALT.GT.ALTL(3).OR.SW(15).EQ.0.) GO TO 17
C      Mixed density at Alt
      DM28=DENSU(ALT,B28,TINF,TLB,XMM,0.,TZ,ZLB,S,MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(3)=DNET(D(3),DM28,ZHM28,XMM,28.)
   17 CONTINUE
      GO TO (20,50,20,25,90,35,40,45,25,48),  J
   20 CONTINUE
C
C       **** HE DENSITY ****
C
C       Density variation factor at Zlb
      G4 = SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,1))
C      Diffusive density at Zlb
      DB04 = PDM(1,1)*EXP(G4)*PD(1,1)
C      Diffusive density at Alt
      D(1)=DENSU(ALT,DB04,TINF,TLB, 4.,-.4,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
      DD=D(1)
      IF(ALT.GT.ALTL(1).OR.SW(15).EQ.0.) GO TO 24
C      Turbopause
      ZH04=PDM(3,1)
      ZHM04=ZHM28
C      Mixed density at Zlb
      B04=DENSU(ZH04,DB04,TINF,TLB,4.-XMM,-1.4,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM04=DENSU(ALT,B04,TINF,TLB,XMM,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(1)=DNET(D(1),DM04,ZHM04,XMM,4.)
C      Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,1)/B04)
      ZC04=PDM(5,1)*PDL(1,2)
      HC04=PDM(6,1)*PDL(2,2)
C      Net density corrected at Alt
      D(1)=D(1)*CCOR(ALT,RL,HC04,ZC04)
   24 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   25 CONTINUE
C
C      **** O DENSITY ****
C
C       Density variation factor at Zlb
      G16= SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,2))
C      Diffusive density at Zlb
      DB16 =  PDM(1,2)*EXP(G16)*PD(1,2)
C       Diffusive density at Alt
      D(2)=DENSU(ALT,DB16,TINF,TLB, 16.,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
      DD=D(2)
      IF(ALT.GT.ALTL(2).OR.SW(15).EQ.0.) GO TO 34
C  Corrected from PDM(3,1) to PDM(3,2)  12/2/85
C       Turbopause
      ZH16=PDM(3,2)
      ZHM16=ZHM28
C      Mixed density at Zlb
      B16=DENSU(ZH16,DB16,TINF,TLB,16-XMM,-1.,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM16=DENSU(ALT,B16,TINF,TLB,XMM,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(2)=DNET(D(2),DM16,ZHM16,XMM,16.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,2)*ABS(PDL(17,2))/B16)
      HC16=PDM(6,2)*PDL(4,2)
      ZC16=PDM(5,2)*PDL(3,2)
      D(2)=D(2)*CCOR(ALT,RL,HC16,ZC16)
C       Chemistry correction
      HCC16=PDM(8,2)*PDL(14,2)
      ZCC16=PDM(7,2)*PDL(13,2)
      RC16=PDM(4,2)*PDL(15,2)
C      Net density corrected at Alt
      D(2)=D(2)*CCOR(ALT,RC16,HCC16,ZCC16)
   34 CONTINUE
      IF(MASS.NE.48 .AND. MASS.NE.49) GO TO 90
   35 CONTINUE
C
C       **** O2 DENSITY ****
C
C       Density variation factor at Zlb
      G32= SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,5))
C      Diffusive density at Zlb
      DB32 = PDM(1,4)*EXP(G32)*PD(1,5)
C       Diffusive density at Alt
      D(4)=DENSU(ALT,DB32,TINF,TLB, 32.,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
      IF(MASS.EQ.49) THEN
         DD=DD+2.*D(4)
      ELSE
         DD=D(4)
      ENDIF
      IF(ALT.GT.ALTL(4).OR.SW(15).EQ.0.) GO TO 39
C       Turbopause
      ZH32=PDM(3,4)
      ZHM32=ZHM28
C      Mixed density at Zlb
      B32=DENSU(ZH32,DB32,TINF,TLB,32.-XMM,-1.,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM32=DENSU(ALT,B32,TINF,TLB,XMM,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(4)=DNET(D(4),DM32,ZHM32,XMM,32.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,4)/B32)
      HC32=PDM(6,4)*PDL(8,2)
      ZC32=PDM(5,4)*PDL(7,2)
C      Net density corrected at Alt
      D(4)=D(4)*CCOR(ALT,RL,HC32,ZC32)
   39 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   40 CONTINUE
C
C       **** AR DENSITY ****
C
C       Density variation factor at Zlb
      G40= SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,6))
C      Diffusive density at Zlb
      DB40 = PDM(1,5)*EXP(G40)*PD(1,6)
C       Diffusive density at Alt
      D(5)=DENSU(ALT,DB40,TINF,TLB, 40.,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
      DD=D(5)
      IF(ALT.GT.ALTL(5).OR.SW(15).EQ.0.) GO TO 44
C       Turbopause
      ZH40=PDM(3,5)
      ZHM40=ZHM28
C      Mixed density at Zlb
      B40=DENSU(ZH40,DB40,TINF,TLB,40.-XMM,-1.,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM40=DENSU(ALT,B40,TINF,TLB,XMM,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(5)=DNET(D(5),DM40,ZHM40,XMM,40.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,5)/B40)
      HC40=PDM(6,5)*PDL(10,2)
      ZC40=PDM(5,5)*PDL(9,2)
C      Net density corrected at Alt
      D(5)=D(5)*CCOR(ALT,RL,HC40,ZC40)
   44 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   45 CONTINUE
C
C        **** HYDROGEN DENSITY ****
C
C       Density variation factor at Zlb
      G1 = SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,7))
C      Diffusive density at Zlb
      DB01 = PDM(1,6)*EXP(G1)*PD(1,7)
C       Diffusive density at Alt
      D(7)=DENSU(ALT,DB01,TINF,TLB,1.,-.4,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
      DD=D(7)
      IF(ALT.GT.ALTL(7).OR.SW(15).EQ.0.) GO TO 47
C       Turbopause
      ZH01=PDM(3,6)
      ZHM01=ZHM28
C      Mixed density at Zlb
      B01=DENSU(ZH01,DB01,TINF,TLB,1.-XMM,-1.4,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM01=DENSU(ALT,B01,TINF,TLB,XMM,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(7)=DNET(D(7),DM01,ZHM01,XMM,1.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,6)*ABS(PDL(18,2))/B01)
      HC01=PDM(6,6)*PDL(12,2)
      ZC01=PDM(5,6)*PDL(11,2)
      D(7)=D(7)*CCOR(ALT,RL,HC01,ZC01)
C       Chemistry correction
      HCC01=PDM(8,6)*PDL(20,2)
      ZCC01=PDM(7,6)*PDL(19,2)
      RC01=PDM(4,6)*PDL(21,2)
C      Net density corrected at Alt
      D(7)=D(7)*CCOR(ALT,RC01,HCC01,ZCC01)
   47 CONTINUE
   48 CONTINUE
C
C        **** ATOMIC NITROGEN DENSITY ****
C
C       Density variation factor at Zlb
      G14 = SW(21)*GLOBE6(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,8))
C      Diffusive density at Zlb
      DB14 = PDM(1,7)*EXP(G14)*PD(1,8)
C       Diffusive density at Alt
      D(8)=DENSU(ALT,DB14,TINF,TLB,14.,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
      DD=D(8)
      IF(ALT.GT.ALTL(8).OR.SW(15).EQ.0.) GO TO 49
C       Turbopause
      ZH14=PDM(3,7)
      ZHM14=ZHM28
C      Mixed density at Zlb
      B14=DENSU(ZH14,DB14,TINF,TLB,14.-XMM,-1.,
     $  T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM14=DENSU(ALT,B14,TINF,TLB,XMM,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
C      Net density at Alt
      D(8)=DNET(D(8),DM14,ZHM14,XMM,14.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,7)*ABS(PDL(3,1))/B14)
      HC14=PDM(6,7)*PDL(2,1)
      ZC14=PDM(5,7)*PDL(1,1)
      D(8)=D(8)*CCOR(ALT,RL,HC14,ZC14)
C       Chemistry correction
      HCC14=PDM(8,7)*PDL(5,1)
      ZCC14=PDM(7,7)*PDL(4,1)
      RC14=PDM(4,7)*PDL(6,1)
C      Net density corrected at Alt
      D(8)=D(8)*CCOR(ALT,RC14,HCC14,ZCC14)
   49 CONTINUE
      IF(MASS.NE.48) GO TO 90
C
C       TOTAL MASS DENSITY
C
      D(6) = 1.66E-24*(4.*D(1)+16.*D(2)+28.*D(3)+32.*D(4)+40.*D(5)+
     &       D(7)+14.*D(8))
      DB48=1.66E-24*(4.*DB04+16.*DB16+28.*DB28+32.*DB32+40.*DB40+DB01+
     &        14.*DB14)
      GO TO 90
C       TEMPERATURE AT ALTITUDE
   50 CONTINUE
      DDUM=DENSU(ALT,1.,TINF,TLB,0.,0.,T(2),ZLB,S,MN1,ZN1,TN1,TGN1)
      GO TO 90
   90 CONTINUE
C       ADJUST DENSITIES FROM CGS TO KGM
      IF(IMR.EQ.1) THEN
        DO 95 I=1,8
          D(I)=D(I)*1.E6
   95   CONTINUE
        D(6)=D(6)/1000.
      ENDIF
      ALAST=ALT
      RETURN
  100 FORMAT(1X,'MASS', I5, '  NOT VALID')
      ENTRY METER6(METER)
      IMR=0
      IF(METER) IMR=1
      END
C-----------------------------------------------------------------------
      FUNCTION GLOBE6(YRD,SEC,LAT,LONG,TLOC,F107A,F107,AP,P)
C       CALCULATE G(L) FUNCTION 
C       Upper Thermosphere Parameters
      REAL LAT, LONG,LONGL
      DIMENSION P(150),SV(25),AP(1)
      COMMON/TTEST/TINF,GB,ROUT,T(15)
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/LPOLY/PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     $ IYR,DAY,DF,DFA,APD,APDF,APT(4),XLONG,CLONG,SLONG
      DATA DGTR/1.74533E-2/,DR/1.72142E-2/, XL/1000./,TLL/1000./
      DATA SW9/1./,DAYL/-1./,P14/-1000./,P18/-1000./,P32/-1000./
      DATA HR/.2618/,SR/7.2722E-5/,SV/25*1./,NSW/14/,P39/-1000./
      DATA LONGL/-999./
C       3hr Magnetica activity functions
      G0(A)=(A-4.+(P(26)-1.)*(A-4.+(EXP(-ABS(P(25))*(A-4.))-1.)/ABS(P(25
     *))))
       SUMEX(EX)=1.+(1.-EX**19)/(1.-EX)*EX**(.5)
      SG0(EX)=(G0(AP(2))+(G0(AP(3))*EX+G0(AP(4))*EX*EX+G0(AP(5))*EX**3
     $ +(G0(AP(6))*EX**4+G0(AP(7))*EX**12)*(1.-EX**8)/(1.-EX))
     $ )/SUMEX(EX)
C
      IF(ISW.NE.64999) CALL TSELEC(SV)
      DO 10 J=1,14
       T(J)=0
   10 CONTINUE
      IF(SW(9).GT.0) SW9=1.
      IF(SW(9).LT.0) SW9=-1.
      IYR = YRD/1000.
      DAY = YRD - IYR*1000.
      XLONG=LONG
      IF(XL.EQ.LAT)   GO TO 15
C          CALCULATE LEGENDRE POLYNOMIALS
      C = SIN(LAT*DGTR)
      S = COS(LAT*DGTR)
      C2 = C*C
      C4 = C2*C2
      S2 = S*S
      PLG(2,1) = C
      PLG(3,1) = 0.5*(3.*C2 -1.)
      PLG(4,1) = 0.5*(5.*C*C2-3.*C)
      PLG(5,1) = (35.*C4 - 30.*C2 + 3.)/8.
      PLG(6,1) = (63.*C2*C2*C - 70.*C2*C + 15.*C)/8.
      PLG(7,1) = (11.*C*PLG(6,1) - 5.*PLG(5,1))/6.
C     PLG(8,1) = (13.*C*PLG(7,1) - 6.*PLG(6,1))/7.
      PLG(2,2) = S
      PLG(3,2) = 3.*C*S
      PLG(4,2) = 1.5*(5.*C2-1.)*S
      PLG(5,2) = 2.5*(7.*C2*C-3.*C)*S
      PLG(6,2) = 1.875*(21.*C4 - 14.*C2 +1.)*S
      PLG(7,2) = (11.*C*PLG(6,2)-6.*PLG(5,2))/5.
C     PLG(8,2) = (13.*C*PLG(7,2)-7.*PLG(6,2))/6.
C     PLG(9,2) = (15.*C*PLG(8,2)-8.*PLG(7,2))/7.
      PLG(3,3) = 3.*S2
      PLG(4,3) = 15.*S2*C
      PLG(5,3) = 7.5*(7.*C2 -1.)*S2
      PLG(6,3) = 3.*C*PLG(5,3)-2.*PLG(4,3)
      PLG(7,3)=(11.*C*PLG(6,3)-7.*PLG(5,3))/4.
      PLG(8,3)=(13.*C*PLG(7,3)-8.*PLG(6,3))/5.
      PLG(4,4) = 15.*S2*S
      PLG(5,4) = 105.*S2*S*C 
      PLG(6,4)=(9.*C*PLG(5,4)-7.*PLG(4,4))/2.
      PLG(7,4)=(11.*C*PLG(6,4)-8.*PLG(5,4))/3.
      XL=LAT
   15 CONTINUE
      IF(TLL.EQ.TLOC)   GO TO 16
      IF(SW(7).EQ.0.AND.SW(8).EQ.0.AND.SW(14).EQ.0) GOTO 16
      STLOC = SIN(HR*TLOC)
      CTLOC = COS(HR*TLOC)
      S2TLOC = SIN(2.*HR*TLOC)
      C2TLOC = COS(2.*HR*TLOC)
      S3TLOC = SIN(3.*HR*TLOC)
      C3TLOC = COS(3.*HR*TLOC)
      TLL = TLOC
   16 CONTINUE
      IF(LONG.NE.LONGL) THEN
        CLONG=COS(DGTR*LONG)
        SLONG=SIN(DGTR*LONG)
      ENDIF
      LONGL=LONG
      IF(DAY.NE.DAYL.OR.P(14).NE.P14) CD14=COS(DR*(DAY-P(14)))
      IF(DAY.NE.DAYL.OR.P(18).NE.P18) CD18=COS(2.*DR*(DAY-P(18)))
      IF(DAY.NE.DAYL.OR.P(32).NE.P32) CD32=COS(DR*(DAY-P(32)))
      IF(DAY.NE.DAYL.OR.P(39).NE.P39) CD39=COS(2.*DR*(DAY-P(39)))
      DAYL = DAY
      P14 = P(14)
      P18 = P(18)
      P32 = P(32)
      P39 = P(39)
C         F10.7 EFFECT
      DF = F107 - F107A
      DFA=F107A-150.
      T(1) =  P(20)*DF + P(21)*DF*DF + P(22)*DFA
     & + P(30)*DFA**2
      F1 = 1. + (P(48)*DFA +P(20)*DF+P(21)*DF*DF)*SWC(1)
      F2 = 1. + (P(50)*DFA+P(20)*DF+P(21)*DF*DF)*SWC(1)
C        TIME INDEPENDENT
      T(2) =
     1  (P(2)*PLG(3,1) + P(3)*PLG(5,1)+P(23)*PLG(7,1))
     & +(P(15)*PLG(3,1))*DFA*SWC(1)
     2 +P(27)*PLG(2,1)
C        SYMMETRICAL ANNUAL
      T(3) =
     1 (P(19) )*CD32
C        SYMMETRICAL SEMIANNUAL
      T(4) =
     1 (P(16)+P(17)*PLG(3,1))*CD18
C        ASYMMETRICAL ANNUAL
      T(5) =  F1*
     1  (P(10)*PLG(2,1)+P(11)*PLG(4,1))*CD14
C         ASYMMETRICAL SEMIANNUAL
      T(6) =    P(38)*PLG(2,1)*CD39
C        DIURNAL
      IF(SW(7).EQ.0) GOTO 200
      T71 = (P(12)*PLG(3,2))*CD14*SWC(5)
      T72 = (P(13)*PLG(3,2))*CD14*SWC(5)
      T(7) = F2*
     1 ((P(4)*PLG(2,2) + P(5)*PLG(4,2) + P(28)*PLG(6,2)
     2 + T71)*CTLOC
     4 + (P(7)*PLG(2,2) + P(8)*PLG(4,2) +P(29)*PLG(6,2)
     5 + T72)*STLOC)
  200 CONTINUE
C        SEMIDIURNAL
      IF(SW(8).EQ.0) GOTO 210
      T81 = (P(24)*PLG(4,3)+P(36)*PLG(6,3))*CD14*SWC(5) 
      T82 = (P(34)*PLG(4,3)+P(37)*PLG(6,3))*CD14*SWC(5)
      T(8) = F2*
     1 ((P(6)*PLG(3,3) + P(42)*PLG(5,3) + T81)*C2TLOC
     3 +(P(9)*PLG(3,3) + P(43)*PLG(5,3) + T82)*S2TLOC)
  210 CONTINUE
C        TERDIURNAL
      IF(SW(14).EQ.0) GOTO 220
      T(14) = F2*
     1 ((P(40)*PLG(4,4)+(P(94)*PLG(5,4)+P(47)*PLG(7,4))*CD14*SWC(5))*
     $ S3TLOC
     2 +(P(41)*PLG(4,4)+(P(95)*PLG(5,4)+P(49)*PLG(7,4))*CD14*SWC(5))*
     $ C3TLOC)
  220 CONTINUE
C          MAGNETIC ACTIVITY BASED ON DAILY AP

      IF(SW9.EQ.-1.) GO TO 30
      APD=(AP(1)-4.)
      P44=P(44)
      P45=P(45)
      IF(P44.LT.0) P44=1.E-5
      APDF = (APD+(P45-1.)*(APD+(EXP(-P44  *APD)-1.)/P44  ))
      IF(SW(9).EQ.0) GOTO 40
      T(9)=APDF*(P(33)+P(46)*PLG(3,1)+P(35)*PLG(5,1)+
     $ (P(101)*PLG(2,1)+P(102)*PLG(4,1)+P(103)*PLG(6,1))*CD14*SWC(5)+
     $ (P(122)*PLG(2,2)+P(123)*PLG(4,2)+P(124)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(125))))
      GO TO 40
   30 CONTINUE
      IF(P(52).EQ.0) GO TO 40
      EXP1 = EXP(-10800.*ABS(P(52))/(1.+P(139)*(45.-ABS(LAT))))
      IF(EXP1.GT..99999) EXP1=.99999
      EXP2 = EXP(-10800.*ABS(P(54)))
      IF(EXP2.GT..99999) EXP2=.99999
      IF(P(25).LT.1.E-4) P(25)=1.E-4
      APT(1)=SG0(EXP1)
      APT(3)=SG0(EXP2)
      IF(SW(9).EQ.0) GOTO 40
      T(9) = APT(1)*(P(51)+P(97)*PLG(3,1)+P(55)*PLG(5,1)+
     $ (P(126)*PLG(2,1)+P(127)*PLG(4,1)+P(128)*PLG(6,1))*CD14*SWC(5)+
     $ (P(129)*PLG(2,2)+P(130)*PLG(4,2)+P(131)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(132))))
  40  CONTINUE
      IF(SW(10).EQ.0.OR.LONG.LE.-1000.) GO TO 49
C        LONGITUDINAL
      IF(SW(11).EQ.0) GOTO 230
      T(11)= (1.+P(81)*DFA*SWC(1))*
     $((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $ +P(104)*PLG(2,2)+P(105)*PLG(4,2)+P(106)*PLG(6,2)
     $ +SWC(5)*(P(110)*PLG(2,2)+P(111)*PLG(4,2)+P(112)*PLG(6,2))*CD14)*
     $     CLONG
     $ +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $ +P(107)*PLG(2,2)+P(108)*PLG(4,2)+P(109)*PLG(6,2)
     $ +SWC(5)*(P(113)*PLG(2,2)+P(114)*PLG(4,2)+P(115)*PLG(6,2))*CD14)*
     $  SLONG)
  230 CONTINUE
C        UT AND MIXED UT,LONGITUDE
      IF(SW(12).EQ.0) GOTO 240
      T(12)=(1.+P(96)*PLG(2,1))*(1.+P(82)*DFA*SWC(1))*
     $(1.+P(120)*PLG(2,1)*SWC(5)*CD14)*
     $((P(69)*PLG(2,1)+P(70)*PLG(4,1)+P(71)*PLG(6,1))*
     $     COS(SR*(SEC-P(72))))
      T(12)=T(12)+SWC(11)*
     $ (P(77)*PLG(4,3)+P(78)*PLG(6,3)+P(79)*PLG(8,3))*
     $     COS(SR*(SEC-P(80))+2.*DGTR*LONG)*(1.+P(138)*DFA*SWC(1))
  240 CONTINUE
C        UT,LONGITUDE MAGNETIC ACTIVITY
      IF(SW(13).EQ.0) GOTO 48
      IF(SW9.EQ.-1.) GO TO 45
      T(13)= APDF*SWC(11)*(1.+P(121)*PLG(2,1))*
     $((P( 61)*PLG(3,2)+P( 62)*PLG(5,2)+P( 63)*PLG(7,2))*
     $     COS(DGTR*(LONG-P( 64))))
     $ +APDF*SWC(11)*SWC(5)*
     $ (P(116)*PLG(2,2)+P(117)*PLG(4,2)+P(118)*PLG(6,2))*
     $     CD14*COS(DGTR*(LONG-P(119)))
     $ + APDF*SWC(12)*
     $ (P( 84)*PLG(2,1)+P( 85)*PLG(4,1)+P( 86)*PLG(6,1))*
     $     COS(SR*(SEC-P( 76)))
      GOTO 48
   45 CONTINUE
      IF(P(52).EQ.0) GOTO 48
      T(13)=APT(1)*SWC(11)*(1.+P(133)*PLG(2,1))*
     $((P(53)*PLG(3,2)+P(99)*PLG(5,2)+P(68)*PLG(7,2))*
     $     COS(DGTR*(LONG-P(98))))
     $ +APT(1)*SWC(11)*SWC(5)*
     $ (P(134)*PLG(2,2)+P(135)*PLG(4,2)+P(136)*PLG(6,2))*
     $     CD14*COS(DGTR*(LONG-P(137)))
     $ +APT(1)*SWC(12)*
     $ (P(56)*PLG(2,1)+P(57)*PLG(4,1)+P(58)*PLG(6,1))*
     $     COS(SR*(SEC-P(59)))
   48 CONTINUE
   49 CONTINUE
      TINF=P(31)
      DO 50 I = 1,NSW
   50 TINF = TINF + ABS(SW(I))*T(I)
      GLOBE6 = TINF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE TSELEC(SV)
C        SET SWITCHES
C        SW FOR MAIN TERMS, SWC FOR CROSS TERMS
      DIMENSION SV(1),SAV(25),SVV(1)
      COMMON/CSW/SW(25),ISW,SWC(25)
      DO 100 I = 1,25
        SAV(I)=SV(I)
        SW(I)=AMOD(SV(I),2.)
        IF(ABS(SV(I)).EQ.1.OR.ABS(SV(I)).EQ.2.) THEN
          SWC(I)=1.
        ELSE
          SWC(I)=0.
        ENDIF
  100 CONTINUE
      ISW=64999
      RETURN
      ENTRY TRETRV(SVV)
      DO 200 I=1,25
        SVV(I)=SAV(I)
  200 CONTINUE
      END
C-----------------------------------------------------------------------
      FUNCTION GLOB6S(P)
C      VERSION OF GLOBE FOR LOWER ATMOSPHERE 1/17/90
      REAL LONG
      COMMON/LPOLY/PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     $ IYR,DAY,DF,DFA,APD,APDF,APT(4),LONG,CLONG,SLONG
      COMMON/CSW/SW(25),ISW,SWC(25)
      DIMENSION P(100),T(14)
      DATA DR/1.72142E-2/,DGTR/1.74533E-2/
      DATA DAYL/-1./,P32,P18,P14,P39/4*-1000./
      DO 10 J=1,14
        T(J)=0.
   10 CONTINUE
      IF(DAY.NE.DAYL.OR.P32.NE.P(32)) CD32=COS(DR*(DAY-P(32)))
      IF(DAY.NE.DAYL.OR.P18.NE.P(18)) CD18=COS(2.*DR*(DAY-P(18)))       
      IF(DAY.NE.DAYL.OR.P14.NE.P(14)) CD14=COS(DR*(DAY-P(14)))
      IF(DAY.NE.DAYL.OR.P39.NE.P(39)) CD39=COS(2.*DR*(DAY-P(39)))
      DAYL=DAY
      P32=P(32)
      P18=P(18)
      P14=P(14)
      P39=P(39)
C
C       F10.7
      T(1)=P(22)*DFA
C       TIME INDEPENDENT
      T(2)=P(2)*PLG(3,1)+P(3)*PLG(5,1)+P(23)*PLG(7,1)
     $     +P(27)*PLG(2,1)+P(28)*PLG(4,1)+P(29)*PLG(6,1)
C       SYMMETRICAL ANNUAL
      T(3)=(P(19)+P(48)*PLG(3,1)+P(30)*PLG(5,1))*CD32
C       SYMMETRICAL SEMIANNUAL
      T(4)=(P(16)+P(17)*PLG(3,1)+P(31)*PLG(5,1))*CD18
C       ASYMMETRICAL ANNUAL
      T(5)=(P(10)*PLG(2,1)+P(11)*PLG(4,1)+P(36)*PLG(6,1))*CD14
C       ASYMMETRICAL SEMIANNUAL
      T(6)=(P(38)*PLG(2,1))*CD39
C        DIURNAL
      IF(SW(7).EQ.0) GOTO 200
      T71 = P(12)*PLG(3,2)*CD14*SWC(5)
      T72 = P(13)*PLG(3,2)*CD14*SWC(5)
      T(7) = 
     1 ((P(4)*PLG(2,2) + P(5)*PLG(4,2)
     2 + T71)*CTLOC
     4 + (P(7)*PLG(2,2) + P(8)*PLG(4,2)
     5 + T72)*STLOC)
  200 CONTINUE
C        SEMIDIURNAL
      IF(SW(8).EQ.0) GOTO 210
      T81 = (P(24)*PLG(4,3)+P(47)*PLG(6,3))*CD14*SWC(5) 
      T82 = (P(34)*PLG(4,3)+P(49)*PLG(6,3))*CD14*SWC(5)
      T(8) = 
     1 ((P(6)*PLG(3,3) + P(42)*PLG(5,3) + T81)*C2TLOC
     3 +(P(9)*PLG(3,3) + P(43)*PLG(5,3) + T82)*S2TLOC)
  210 CONTINUE
C        TERDIURNAL
      IF(SW(14).EQ.0) GOTO 220
      T(14) = P(40)*PLG(4,4)*S3TLOC
     $ +P(41)*PLG(4,4)*C3TLOC
  220 CONTINUE
C       MAGNETIC ACTIVITY
      IF(SW(9).EQ.0) GOTO 40
      IF(SW(9).EQ.1)
     $ T(9)=APDF*(P(33)+P(46)*PLG(3,1)*SWC(2))
      IF(SW(9).EQ.-1)
     $ T(9)=(P(51)*APT(3)+P(97)*PLG(3,1)*APT(3)*SWC(2))
   40 CONTINUE
      IF(SW(10).EQ.0.OR.SW(11).EQ.0.OR.LONG.LE.-1000.) GO TO 49
C        LONGITUDINAL
      T(11)= (1.+PLG(2,1)*(P(81)*SWC(5)*COS(DR*(DAY-P(82)))
     $           +P(86)*SWC(6)*COS(2.*DR*(DAY-P(87))))
     $        +P(84)*SWC(3)*COS(DR*(DAY-P(85)))
     $           +P(88)*SWC(4)*COS(2.*DR*(DAY-P(89))))
     $ *((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $   +P(75)*PLG(2,2)+P(76)*PLG(4,2)+P(77)*PLG(6,2)
     $    )*CLONG
     $  +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $   +P(78)*PLG(2,2)+P(79)*PLG(4,2)+P(80)*PLG(6,2)
     $    )*SLONG)
   49 CONTINUE
      TT=0.
      DO 50 I=1,14
   50 TT=TT+ABS(SW(I))*T(I)
      GLOB6S=TT
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION DENSU(ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,
     $  MN1,ZN1,TN1,TGN1)
C       Calculate Temperature and Density Profiles for MSIS models
C       New lower thermo polynomial 10/30/89
      DIMENSION ZN1(MN1),TN1(MN1),TGN1(2),XS(5),YS(5),Y2OUT(5)
      COMMON/PARMB/GSURF,RE
      COMMON/LSQV/MP,II,JG,LT,QPB(50),IERR,IFUN,N,J,DV(60)
      DATA RGAS/831.4/
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
CCCCCCWRITE(6,*) 'DB',ALT,DLB,TINF,TLB,XM,ALPHA,ZLB,S2,MN1,ZN1,TN1
      DENSU=1.
C        Joining altitude of Bates and spline
      ZA=ZN1(1)
      Z=AMAX1(ALT,ZA)
C      Geopotential altitude difference from ZLB
      ZG2=ZETA(Z,ZLB)
C      Bates temperature
      TT=TINF-(TINF-TLB)*EXP(-S2*ZG2)
      TA=TT
      TZ=TT
      DENSU=TZ
      IF(ALT.GE.ZA) GO TO 10
C
C       CALCULATE TEMPERATURE BELOW ZA
C      Temperature gradient at ZA from Bates profile
      DTA=(TINF-TA)*S2*((RE+ZLB)/(RE+ZA))**2
      TGN1(1)=DTA 
      TN1(1)=TA
      Z=AMAX1(ALT,ZN1(MN1))
      MN=MN1
      Z1=ZN1(1)
      Z2=ZN1(MN)
      T1=TN1(1)
      T2=TN1(MN)
C      Geopotental difference from Z1
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 20 K=1,MN
        XS(K)=ZETA(ZN1(K),Z1)/ZGDIF
        YS(K)=1./TN1(K)
   20 CONTINUE
C        End node derivatives
      YD1=-TGN1(1)/(T1*T1)*ZGDIF
      YD2=-TGN1(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)
C       temperature at altitude
      TZ=1./Y
      DENSU=TZ
   10 IF(XM.EQ.0.) GO TO 50
C
C      CALCULATE DENSITY ABOVE ZA
      GLB=GSURF/(1.+ZLB/RE)**2
      GAMMA=XM*GLB/(S2*RGAS*TINF)
      EXPL=EXP(-S2*GAMMA*ZG2)
      IF(EXPL.GT.50.OR.TT.LE.0.) THEN
        EXPL=50.
      ENDIF
C       Density at altitude
      DENSA=DLB*(TLB/TT)**(1.+ALPHA+GAMMA)*EXPL
      DENSU=DENSA
      IF(ALT.GE.ZA) GO TO 50
C
C      CALCULATE DENSITY BELOW ZA
      GLB=GSURF/(1.+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C       integrate spline temperatures
      CALL SPLINI(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50..OR.TZ.LE.0.) THEN
        EXPL=50.
      ENDIF
C       Density at altitude
      DENSU=DENSU*(T1/TZ)**(1.+ALPHA)*EXP(-EXPL)
   50 CONTINUE
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION DENSM(ALT,D0,XM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
C       Calculate Temperature and Density Profiles for lower atmos.
      DIMENSION ZN3(MN3),TN3(MN3),TGN3(2),XS(10),YS(10),Y2OUT(10)
      DIMENSION ZN2(MN2),TN2(MN2),TGN2(2)
      COMMON/PARMB/GSURF,RE
      COMMON/FIT/TAF
      COMMON/LSQV/MP,II,JG,LT,QPB(50),IERR,IFUN,N,J,DV(60)
      DATA RGAS/831.4/
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
      DENSM=D0
      IF(ALT.GT.ZN2(1)) GOTO 50
C      STRATOSPHERE/MESOSPHERE TEMPERATURE
      Z=AMAX1(ALT,ZN2(MN2))
      MN=MN2
      Z1=ZN2(1)
      Z2=ZN2(MN)
      T1=TN2(1)
      T2=TN2(MN)
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 210 K=1,MN
        XS(K)=ZETA(ZN2(K),Z1)/ZGDIF
        YS(K)=1./TN2(K)
  210 CONTINUE
      YD1=-TGN2(1)/(T1*T1)*ZGDIF
      YD2=-TGN2(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)
C       Temperature at altitude
      TZ=1./Y
      IF(XM.EQ.0.) GO TO 20
C
C      CALCULATE STRATOSPHERE/MESOSPHERE DENSITY 
      GLB=GSURF/(1.+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C       Integrate temperature profile
      CALL SPLINI(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.) EXPL=50.
C       Density at altitude
      DENSM=DENSM*(T1/TZ)*EXP(-EXPL)
   20 CONTINUE
      IF(ALT.GT.ZN3(1)) GOTO 50
C
C      TROPOSPHERE/STRATOSPHERE TEMPERATURE
      Z=ALT
      MN=MN3
      Z1=ZN3(1)
      Z2=ZN3(MN)
      T1=TN3(1)
      T2=TN3(MN)
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 220 K=1,MN
        XS(K)=ZETA(ZN3(K),Z1)/ZGDIF
        YS(K)=1./TN3(K)
  220 CONTINUE
      YD1=-TGN3(1)/(T1*T1)*ZGDIF
      YD2=-TGN3(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)
C       temperature at altitude
      TZ=1./Y
      IF(XM.EQ.0.) GO TO 30
C
C      CALCULATE TROPOSPHERIC/STRATOSPHERE DENSITY 
C     
      GLB=GSURF/(1.+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C        Integrate temperature profile
      CALL SPLINI(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.) EXPL=50.
C        Density at altitude
      DENSM=DENSM*(T1/TZ)*EXP(-EXPL)
   30 CONTINUE
   50 CONTINUE
      IF(XM.EQ.0) DENSM=TZ
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
C        CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
C        ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
C        X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        N: SIZE OF ARRAYS X,Y
C        YP1,YPN: SPECIFIED DERIVATIVES AT X(1) AND X(N); VALUES
C                 >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
C        Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
      PARAMETER (NMAX=100)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF(YP1.GT..99E30) THEN
        Y2(1)=0
        U(1)=0
      ELSE
        Y2(1)=-.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     $    /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
   11 CONTINUE
      IF(YPN.GT..99E30) THEN
        QN=0
        UN=0
      ELSE
        QN=.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
   12 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
C        CALCULATE CUBIC SPLINE INTERP VALUE
C        ADAPTED FROM NUMBERICAL RECIPES BY PRESS ET AL.
C        XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        Y2A: ARRAY OF SECOND DERIVATIVES
C        N: SIZE OF ARRAYS XA,YA,Y2A
C        X: ABSCISSA FOR INTERPOLATION
C        Y: OUTPUT VALUE
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
    1 CONTINUE
      IF(KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X) THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF(H.EQ.0) WRITE(6,*) 'BAD XA INPUT TO SPLINT'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     $  ((A*A*A-A)*Y2A(KLO)+(B*B*B-B)*Y2A(KHI))*H*H/6.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINI(XA,YA,Y2A,N,X,YI)
C       INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
C        XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        Y2A: ARRAY OF SECOND DERIVATIVES
C        N: SIZE OF ARRAYS XA,YA,Y2A
C        X: ABSCISSA ENDPOINT FOR INTEGRATION
C        Y: OUTPUT VALUE
      DIMENSION XA(N),YA(N),Y2A(N)
      YI=0
      KLO=1
      KHI=2
    1 CONTINUE
      IF(X.GT.XA(KLO).AND.KHI.LE.N) THEN
        XX=X
        IF(KHI.LT.N) XX=AMIN1(X,XA(KHI))
        H=XA(KHI)-XA(KLO)
        A=(XA(KHI)-XX)/H
        B=(XX-XA(KLO))/H
        A2=A*A
        B2=B*B
        YI=YI+((1.-A2)*YA(KLO)/2.+B2*YA(KHI)/2.+
     $     ((-(1.+A2*A2)/4.+A2/2.)*Y2A(KLO)+
     $     (B2*B2/4.-B2/2.)*Y2A(KHI))*H*H/6.)*H
        KLO=KLO+1
        KHI=KHI+1
        GOTO 1
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION DNET(DD,DM,ZHM,XMM,XM)
C       TURBOPAUSE CORRECTION FOR MSIS MODELS
C         Root mean density
C       8/20/80
C          DD - diffusive density
C          DM - full mixed density
C          ZHM - transition scale length
C          XMM - full mixed molecular weight
C          XM  - species molecular weight
C          DNET - combined density
      A=ZHM/(XMM-XM)
      IF(DM.GT.0.AND.DD.GT.0) GOTO 5
        WRITE(6,*) 'DNET LOG ERROR',DM,DD,XM
        IF(DD.EQ.0.AND.DM.EQ.0) DD=1.
        IF(DM.EQ.0) GOTO 10
        IF(DD.EQ.0) GOTO 20
    5 CONTINUE
      YLOG=A*ALOG(DM/DD)
      IF(YLOG.LT.-10.) GO TO 10
      IF(YLOG.GT.10.)  GO TO 20
        DNET=DD*(1.+EXP(YLOG))**(1/A)
        GO TO 50
   10 CONTINUE
        DNET=DD
        GO TO 50
   20 CONTINUE
        DNET=DM
        GO TO 50
   50 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION  CCOR(ALT, R,H1,ZH)
C        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
C        ALT - altitude
C        R - target ratio
C        H1 - transition scale length
C        ZH - altitude of 1/2 R
      E=(ALT-ZH)/H1
      IF(E.GT.70.) GO TO 20
      IF(E.LT.-70.) GO TO 10
        EX=EXP(E)
        CCOR=R/(1.+EX)
        GO TO 50
   10   CCOR=R
        GO TO 50
   20   CCOR=0.
        GO TO 50
   50 CONTINUE
      CCOR=EXP(CCOR)
       RETURN
      END
C-----------------------------------------------------------------------
      BLOCK DATA GTD6BK
C          MSISE 90 12-MAR-90   
      COMMON/PARM6/PT1(50),PT2(50),PT3(50),PA1(50),PA2(50),PA3(50),
     $ PB1(50),PB2(50),PB3(50),PC1(50),PC2(50),PC3(50),
     $ PD1(50),PD2(50),PD3(50),PE1(50),PE2(50),PE3(50),
     $ PF1(50),PF2(50),PF3(50),PG1(50),PG2(50),PG3(50),
     $ PH1(50),PH2(50),PH3(50),PI1(50),PI2(50),PI3(50),
     $ PJ1(50),PJ2(50),PJ3(50),PK1(50),PL1(50),PL2(50),
     $ PM1(50),PM2(50),PN1(50),PN2(50),PO1(50),PO2(50),
     $ PP1(50),PP2(50),PQ1(50),PQ2(50),PR1(50),PR2(50),
     $ PS1(50),PS2(50),PU1(50),PU2(50),PV1(50),PV2(50),
     $ PW1(50),PW2(50),PX1(50),PX2(50),PY1(50),PY2(50),
     $ PZ1(50),PZ2(50)
      COMMON/LOWER6/PTM(10),PDM(10,8)
      COMMON/MAVG6/PAVGM(10)
      COMMON/DATIM6/ISDATE(3),ISTIME(2),NAME(2)
      COMMON/METSEL/IMR
      DATA IMR/0/
      DATA ISDATE/'12-M','AR-9','0   '/,ISTIME/'15:0','9:04'/
      DATA NAME/'MSIS','E 90'/
C         TEMPERATURE
      DATA PT1/
     *  9.96040E-01, 3.85528E-02, 3.03445E-03,-1.05531E-01,-6.07134E-03,
     * -5.16278E-04,-1.15622E-01, 2.02240E-03, 9.90156E-03,-1.27371E-01,
     * -3.02449E-02, 1.23512E-02,-5.26277E-03,-8.45398E+00, 0.00000E+00,
     *  1.42370E-02, 0.00000E+00, 1.25818E+02, 8.05486E-03, 1.64419E-03,
     * -6.21452E-06, 3.11701E-03, 0.00000E+00, 3.86578E-03, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,-6.41110E-06,
     *  0.00000E+00, 3.00150E+01, 5.33297E-03, 3.89146E-03, 2.04725E-03,
     *  0.00000E+00, 0.00000E+00,-1.92645E-02, 2.75905E+00, 1.47284E-03,
     *  3.41345E-04,-1.17388E-03,-3.54589E-04, 1.13139E-01, 1.69134E-01,
     *  5.08295E-03, 3.65016E-05, 4.26385E-03, 1.15102E-04, 5.11819E-03/
      DATA PT2/
     *  6.09108E-03, 4.04995E-05, 1.53049E-03, 2.41470E-05, 2.30764E-03,
     *  1.55267E-03, 1.33722E-03,-1.82318E-03,-2.63007E+02, 0.00000E+00,
     *  1.37337E-03, 9.95774E-04, 0.00000E+00,-1.08983E+02, 5.62606E-03,
     *  5.94053E-03, 1.09358E-03, 0.00000E+00,-1.33410E-02,-2.43409E-02,
     * -1.35688E-02, 3.11370E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -2.83023E+03, 8.45583E-04, 5.38706E-04, 0.00000E+00, 2.47956E+02,
     *  2.92246E-03, 0.00000E+00, 0.00000E+00, 7.47703E-05, 8.87993E-04,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -1.16540E-02,-4.49173E-03,-3.53189E-04,-1.73933E-04,-1.53218E-04,
     * -5.65411E-01, 7.77272E-03,-9.11784E+01, 6.45187E-04, 0.00000E+00/
      DATA PT3/
     * -8.37685E-04, 2.42318E-03, 4.73796E-03,-3.01801E-03,-4.23564E-03,
     * -2.48289E-03, 9.19286E-04, 2.16372E-03, 8.63968E-04, 1.89689E-03,
     *  4.15654E-03, 0.00000E+00, 1.18068E-02, 3.31190E-03, 0.00000E+00,
     *  1.20222E-03, 0.00000E+00, 0.00000E+00,-3.07246E+00, 0.00000E+00,
     *  0.00000E+00, 6.72403E-04, 1.08930E-03, 9.72278E-04, 4.68242E+00,
     * -3.15034E-04, 4.00059E-03, 5.15036E-03, 1.62989E-03, 1.08824E-03,
     *  9.95261E-04, 4.18955E+00,-3.64059E-01, 1.70182E-03, 0.00000E+00,
     *  0.00000E+00,-3.20120E+00, 0.00000E+00, 5.80206E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         HE DENSITY
      DATA PA1/
     *  1.04934E+00,-2.88362E-02,-2.07095E-01,-1.03314E-01,-7.02373E-03,
     *  1.29664E-02, 4.08853E-01,-9.19895E-03,-1.88660E-02, 1.40927E+00,
     *  1.75033E-01, 1.87351E-02, 1.10979E-01,-7.42871E+00, 0.00000E+00,
     *  2.67143E-01,-5.95979E-02, 1.05038E+02,-8.40963E-02,-6.97632E-04,
     *  2.06521E-06, 7.65306E-04, 0.00000E+00, 0.00000E+00, 1.26762E-01,
     *  1.28876E-01,-5.04479E-02,-1.30735E-02,-2.24348E-02, 0.00000E+00,
     *  0.00000E+00,-1.50832E+02,-6.29928E-03, 0.00000E+00,-4.07760E-03,
     *  0.00000E+00, 0.00000E+00, 5.25725E-02,-3.11486E+01,-3.13351E-03,
     *  2.75838E-03, 0.00000E+00, 0.00000E+00, 1.11247E-01, 1.08815E-01,
     * -4.66713E-02, 0.00000E+00,-3.29329E-03, 0.00000E+00, 1.67838E-03/
      DATA PA2/
     * -9.16691E-03, 3.45044E-05,-9.71806E-03, 0.00000E+00,-2.04672E-03,
     * -7.86899E-03,-7.98285E-03, 5.36515E-03,-5.31172E+03, 0.00000E+00,
     * -6.42781E-03,-1.71690E-03, 0.00000E+00,-6.79131E+01,-1.79912E-02,
     * -1.58305E-02,-7.12313E-03, 0.00000E+00, 2.53477E-02, 8.52960E-02,
     *  1.02163E-01, 2.95009E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -6.84625E+03,-6.19098E-03,-2.69289E-03, 0.00000E+00,-5.20231E+02,
     * -6.33463E-03, 0.00000E+00, 0.00000E+00,-6.02428E-03,-4.07077E-03,
     *  5.42264E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  4.07560E-02, 2.82288E-02, 9.08088E-03, 0.00000E+00, 0.00000E+00,
     * -4.05204E-01,-5.97931E-02,-7.31823E+01,-2.06620E-03, 0.00000E+00/
      DATA PA3/
     * -3.72723E-03,-1.88146E-02,-1.01794E-02, 8.04633E-03, 1.01090E-02,
     *  8.73253E-03, 2.38268E-02, 4.80444E-03, 1.71088E-03, 3.96369E-02,
     * -2.13809E-02, 0.00000E+00,-1.02588E-01,-5.91702E-03, 0.00000E+00,
     *  2.70923E-03, 0.00000E+00, 0.00000E+00,-1.75043E+02, 6.03489E-01,
     * -6.17589E-01, 8.38098E-03, 1.83871E-03,-7.05329E-04,-4.06644E+00,
     * -5.09347E-03,-2.84344E-02,-1.24160E-02, 1.33665E-02, 3.93410E-03,
     * -5.03723E-04,-4.57683E+00,-5.29542E-01,-4.25812E-03, 0.00000E+00,
     *  0.00000E+00, 1.91541E+01, 0.00000E+00, 3.23247E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         O DENSITY
      DATA PB1/
     *  9.31113E-01,-1.38721E-01,-1.33457E-01,-5.29542E-02,-4.44983E-03,
     *  1.35264E-02, 5.98075E-02,-3.62880E-02,-3.12798E-02, 3.72068E-01,
     *  2.95974E-02, 1.20509E-02, 5.21995E-02,-7.78888E+00, 0.00000E+00,
     *  1.18634E-01,-2.04495E-02, 1.03280E+02, 9.82432E-02, 4.77694E-04,
     *  0.00000E+00, 2.74372E-03, 0.00000E+00, 0.00000E+00, 7.57809E-02,
     *  1.71403E-01,-1.05205E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-8.73348E+00,-5.81094E-03, 0.00000E+00,-8.14944E-03,
     *  0.00000E+00, 0.00000E+00, 5.17255E-02,-1.53028E+01,-3.48932E-03,
     *  9.61771E-04, 5.57732E-03,-4.54180E-04, 9.88213E-02, 9.40456E-02,
     * -3.18797E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.32122E-03/
      DATA PB2/
     * -6.00220E-03, 2.77654E-05,-3.22019E-03, 0.00000E+00,-3.78551E-03,
     * -3.34809E-03,-1.70668E-03, 0.00000E+00, 6.36184E+03, 0.00000E+00,
     *  1.59986E-03,-3.88204E-03,-1.64825E-03,-7.47955E+01,-1.05360E-02,
     * -9.45723E-03,-1.59824E-03,-7.06730E-04,-1.68513E-02,-1.13023E-01,
     * -6.36637E-02,-1.37709E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -1.52368E+04,-5.86061E-03,-2.53108E-03, 0.00000E+00,-2.54837E+03,
     * -3.28988E-03, 0.00000E+00, 0.00000E+00,-2.76364E-03, 9.67923E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  4.34255E-02, 1.14020E-02,-6.18447E-03, 0.00000E+00, 0.00000E+00,
     * -3.02568E-01,-3.27694E-02,-6.71589E+01,-2.28340E-03, 0.00000E+00/
      DATA PB3/
     *  3.06230E-03,-4.65113E-03,-9.73421E-03, 1.28326E-02, 7.88553E-03,
     *  7.97197E-03,-1.20760E-02,-7.67547E-03,-1.20755E-03,-2.98523E-02,
     * -1.26560E-02, 0.00000E+00,-5.68350E-02,-1.53039E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.42911E-03,-4.01347E-03,-2.19074E-03, 3.11281E+00,
     *  3.23251E-03,-6.39523E-03,-6.63069E-03,-3.04403E-04,-4.01920E-03,
     * -1.18708E-03, 4.15211E+00,-2.01896E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         N2 DENSITY
      DATA PC1/
     *  1.06903E+00, 0.00000E+00, 0.00000E+00, 3.66210E-03, 0.00000E+00,
     *  1.90412E-02,-1.78929E-03, 0.00000E+00,-3.92257E-02,-1.19444E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.45398E+00, 0.00000E+00,
     *  2.08180E-02, 0.00000E+00, 1.39638E+02, 8.98481E-02, 0.00000E+00,
     *  0.00000E+00, 3.77113E-04, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-2.36325E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.43022E-03,
     * -3.99776E-06, 6.32343E-03, 5.48144E-03, 1.13139E-01, 1.69134E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PC2/
     *  0.00000E+00, 2.41470E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PC3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         TLB
      DATA PD1/
     *  9.76619E-01, 0.00000E+00, 0.00000E+00,-2.00200E-02, 0.00000E+00,
     * -9.38391E-03,-1.95833E-03, 0.00000E+00, 1.31480E-02,-1.92414E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.45398E+00, 0.00000E+00,
     *  1.07674E-02, 0.00000E+00, 8.93820E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 5.68478E-04, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 4.66814E-03, 0.00000E+00, 0.00000E+00,
     *  5.11651E-05, 2.55717E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-2.60147E-03,-8.08556E-04, 1.13139E-01, 1.69134E-01,
     *  6.64196E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PD2/
     *  5.82026E-03, 2.41470E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 6.21998E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PD3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         O2 DENSITY
      DATA PE1/
     *  9.31402E-01, 1.37976E-01, 0.00000E+00, 3.23736E-04, 0.00000E+00,
     * -9.10906E-03, 7.07506E-02, 0.00000E+00,-5.16650E-02, 6.89755E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.45398E+00, 0.00000E+00,
     *  2.81140E-02, 0.00000E+00, 7.36009E+01, 5.96604E-02, 0.00000E+00,
     *  0.00000E+00,-1.51792E-03, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 9.48758E+00, 8.84541E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.13139E-01, 1.69134E-01,
     *  1.45192E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PE2/
     *  1.07906E-02, 2.99942E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.48930E-02,
     * -7.87184E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -6.83420E-02,-4.41778E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.29730E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PE3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         AR DENSITY
      DATA PF1/
     *  8.68053E-01, 2.36364E-01, 1.34306E-01, 1.03086E-02, 0.00000E+00,
     * -3.79164E-03,-1.57806E-01, 0.00000E+00,-5.87644E-02,-3.12508E-01,
     *  0.00000E+00, 4.37387E-02,-3.54091E-02,-2.23636E+01, 0.00000E+00,
     * -5.33976E-02, 0.00000E+00, 1.14091E+02, 5.17497E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 3.42702E+02, 1.57033E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-3.66278E-03,
     * -1.16193E-03, 0.00000E+00, 0.00000E+00, 1.13139E-01, 1.69134E-01,
     *  1.78431E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PF2/
     *  1.62864E-02, 3.16963E-05, 1.27968E-02, 0.00000E+00, 0.00000E+00,
     * -7.04599E-03, 2.07921E-03, 6.36660E-03, 2.29940E+04, 0.00000E+00,
     *  1.27833E-02,-2.08036E-03,-4.61820E-03,-6.29391E+01,-1.20745E-02,
     *  1.36675E-02, 1.36011E-02,-5.37162E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.92509E+04, 8.35522E-03, 4.19439E-03, 0.00000E+00, 1.20366E+04,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-1.00034E-02,-2.33267E-03,
     *  9.72374E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -2.65079E-02,-2.09125E-02,-1.09465E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.17252E-02,-7.12385E+01,-1.89428E-03, 0.00000E+00/
      DATA PF3/
     * -6.02006E-03, 1.69058E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.90646E-02,
     *  3.48971E-03, 0.00000E+00, 5.01174E-02, 5.50595E-02, 0.00000E+00,
     * -9.55897E-03, 0.00000E+00, 0.00000E+00,-1.51693E+03, 0.00000E+00,
     *  0.00000E+00, 1.29306E-02, 2.69567E-03, 0.00000E+00, 3.92243E+00,
     * -8.47690E-03, 1.16896E-02, 0.00000E+00, 1.48967E-02, 5.44521E-03,
     *  0.00000E+00, 5.64918E+00, 0.00000E+00,-7.72178E-03, 0.00000E+00,
     *  0.00000E+00,-7.34042E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          H DENSITY
      DATA PG1/
     *  1.27515E+00,-2.10472E-01,-1.77924E-01, 2.18900E-01, 2.88436E-02,
     *  1.90077E-02, 2.91001E-01, 2.17437E-02,-1.05186E-02, 4.36141E-01,
     *  1.07605E-01, 3.30755E-02, 4.00581E-02,-9.58051E+00, 0.00000E+00,
     *  1.54028E-02, 0.00000E+00, 7.34194E+01, 4.96540E-02,-5.95906E-03,
     *  3.84512E-05,-1.36000E-02, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01,-4.16610E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.46276E+02,-1.98408E-02, 0.00000E+00, 1.32530E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.04687E-04,
     * -1.47562E-03, 0.00000E+00, 0.00000E+00, 1.13139E-01, 1.69134E-01,
     * -1.26913E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,-6.08370E-03/
      DATA PG2/
     * -2.57587E-02, 3.19022E-05, 0.00000E+00, 0.00000E+00, 1.56644E-02,
     *  1.03640E-02, 1.05771E-03, 0.00000E+00, 3.57949E+03, 0.00000E+00,
     * -1.25672E-03, 1.52783E-03, 1.30518E-03, 7.55558E+00,-9.20341E-03,
     * -2.09142E-02,-1.34106E-02, 0.00000E+00,-4.83312E-02, 8.30900E-02,
     *  9.88009E-02,-1.41148E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -1.05513E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 6.73442E-03, 2.01691E-03,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  5.98019E-02, 6.33298E-03,-1.12871E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.28604E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PG3/
     * -4.94960E-03,-1.36415E-02,-1.15039E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-5.86860E-03,-1.41732E-03, 2.13697E-03, 2.63845E+00,
     * -8.34186E-03,-1.87336E-02,-1.90870E-02,-8.03810E-03,-2.84279E-03,
     *  2.56722E-03, 1.71429E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          N DENSITY
      DATA PH1/
     *  5.73587E+01,-3.98747E-01, 0.00000E+00,-5.29554E-01,-5.82186E-03,
     *  7.14177E-02,-6.79279E-01,-1.67715E-01,-6.42434E-02,-2.11569E-01,
     * -1.59922E-01,-1.71024E-04,-1.15885E-01, 6.51603E+00, 0.00000E+00,
     * -1.76683E-01, 6.50395E-02, 1.43504E+00, 9.28208E-02, 5.11662E-03,
     *  0.00000E+00, 9.95121E-03, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 1.01451E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 5.67667E+01, 2.38192E-03, 0.00000E+00,-1.88240E-02,
     *  0.00000E+00, 0.00000E+00, 4.76218E-02, 2.35206E+01, 4.75901E-03,
     *  5.76162E-03, 1.51815E-02,-1.92730E-02, 1.13139E-01, 1.69134E-01,
     * -2.88771E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.18418E-03/
      DATA PH2/
     * -3.68927E-03, 3.14704E-05, 8.82198E-03, 0.00000E+00,-1.92562E-02,
     * -2.58674E-03,-2.19913E-02, 0.00000E+00, 4.38655E+03, 0.00000E+00,
     *  7.60126E-03, 2.59438E-03, 1.72310E-03, 7.79204E+01, 7.97786E-04,
     * -7.70510E-03, 1.90982E-03, 2.72707E-03, 1.01016E-02, 1.16537E-01,
     * -3.12236E-03, 1.39783E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -1.30712E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-3.20544E-03,-2.06970E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.59010E-02,-1.91427E-03,-3.42829E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-3.45379E-02, 8.94518E+01, 1.71556E-03, 0.00000E+00/
      DATA PH3/
     * -7.65278E-03,-2.08987E-04,-1.57393E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-8.60673E-03,-1.19922E-02,-6.46356E-03,-3.00107E+00,
     * -9.32511E-03,-1.50205E-02,-8.67835E-03,-7.64801E-03,-1.31495E-02,
     * -6.76720E-03,-1.82396E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C        SPARE
      DATA PI1/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.45398E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.13139E-01, 1.69134E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PI2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PI3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          S PARAM  
      DATA PJ1/
     *  9.51363E-01,-4.67542E-02, 1.20260E-01, 0.00000E+00, 0.00000E+00,
     *  1.91357E-02, 0.00000E+00, 0.00000E+00, 1.25429E-03,-1.33240E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-8.45398E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.52317E-03, 0.00000E+00,-9.73404E-03, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-7.18482E-04, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 7.87683E-03,-2.33698E-03, 1.13139E-01, 1.69134E-01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PJ2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PJ3/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TURBO
      DATA PK1/
     *  9.33804E-01, 5.47446E+00, 1.53263E-01, 9.19303E-01, 1.64109E+01,
     *  4.27083E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.40925E-01,
     *  1.15897E+00, 4.71094E-01, 1.09459E+00, 5.25012E+00, 1.00000E+00,
     *  1.00000E+00, 1.03999E+00, 7.67132E-01, 1.10514E+00, 1.75636E+00,
     *  1.10845E+00, 2.33439E+00, 7.96532E-01, 4.31520E+00, 4.07300E+00,
     *  1.22807E+02, 2.39547E-01, 2.53791E-06, 8.42931E-01, 1.04192E+00,
     *  2.00202E+00, 1.00000E+00, 1.00000E+00, 1.00000E+00, 9.62736E-01/
C         LOWER BOUNDARY
      DATA PTM/
     L  1.04130E+03, 3.86000E+02, 1.95000E+02, 1.66728E+01, 2.13000E+02,
     L  1.20000E+02, 2.40000E+02, 1.87000E+02,-2.00000E+00, 0.00000E+00/
      DATA PDM/
     L  2.45600E+07, 6.71072E-06, 1.00000E+02, 0.00000E+00, 1.10000E+02,
     L  1.00000E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  8.59400E+10, 5.40000E-01, 1.05000E+02,-8.00000E+00, 1.10000E+02,
     L  1.00000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  2.81000E+11, 0.00000E+00, 1.05000E+02, 2.80000E+01, 2.89500E+01,
     L  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  3.30000E+10, 2.68270E-01, 1.05000E+02, 0.00000E+00, 1.10000E+02,
     L  1.00000E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  1.33000E+09, 1.19615E-02, 1.05000E+02, 0.00000E+00, 1.10000E+02,
     L  1.00000E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  1.76100E+05, 1.00000E+00, 9.50000E+01,-8.00000E+00, 1.10000E+02,
     L  1.00000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  1.00000E+07, 1.00000E+00, 1.05000E+02,-8.00000E+00, 1.10000E+02,
     L  1.00000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 0.00000E+00,
C
     L  1.00000E+07, 1.00000E+00, 1.05000E+02,-8.00000E+00, 1.10000E+02,
     L  1.00000E+01, 9.00000E+01, 2.00000E+00, 0.00000E+00, 0.00000E+00/
C         TN1(2)
      DATA PL1/
     *  1.02083E+00, 4.08449E-02,-2.34582E-02, 4.38274E-04,-1.52380E-02,
     * -2.09089E-02, 4.46355E-03,-3.41250E-03,-1.12961E-02,-7.03277E-02,
     * -4.82724E-02, 0.00000E+00, 0.00000E+00,-6.20496E+00, 0.00000E+00,
     * -9.80197E-03,-1.45065E-02,-1.13226E+02, 2.28455E-02, 0.00000E+00,
     *  0.00000E+00, 4.93658E-04, 0.00000E+00, 3.79078E-03, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-8.89051E+03, 2.25900E-03, 1.76142E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.55015E-04,
     *  2.21388E-03,-5.99073E-04,-3.52331E-03, 1.13139E-01, 1.69134E-01,
     *  7.79156E-03,-1.93458E-03,-1.08596E-02,-4.39285E-04, 0.00000E+00/
      DATA PL2/
     *  3.83994E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 6.76608E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         TN1(3)
      DATA PM1/
     *  9.24880E-01, 7.41986E-02,-6.37629E-03, 6.00575E-03, 1.29382E-03,
     *  6.97550E-03,-1.70782E-03, 2.80584E-03,-8.87214E-03,-4.35703E-02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 4.31515E+00, 0.00000E+00,
     * -1.81474E-02,-6.06627E-02,-8.43503E+01, 8.46944E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-2.17081E-02,-2.19500E-03, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.47580E+02, 4.41585E-03, 7.80466E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 6.44155E-04,
     * -2.49166E-03, 2.90482E-03,-3.40501E-04, 1.13139E-01, 1.69134E-01,
     * -6.01460E-03,-1.63368E-03, 0.00000E+00,-4.31340E-03, 0.00000E+00/
      DATA PM2/
     *  4.53979E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-5.43660E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         TN1(4)
      DATA PN1/
     *  9.72669E-01,-4.26748E-02, 1.12876E-02,-8.44951E-03, 7.04114E-03,
     *  1.26036E-02,-3.88164E-03,-5.20509E-04,-6.09710E-04, 1.31603E-01,
     *  1.13804E-01, 0.00000E+00, 0.00000E+00,-6.15970E+00, 0.00000E+00,
     * -2.14214E-02,-6.62913E-02,-2.02884E-01, 2.35350E-02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.13573E-02,-1.84905E-03, 1.32397E-01,
     *  2.13315E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 1.42645E+00,-2.64405E-03,-5.57771E-04, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-2.20621E+01,-1.10313E-03,
     *  3.97063E-05, 5.47632E-05, 3.57577E-03, 1.13139E-01, 1.69134E-01,
     *  0.00000E+00, 1.18897E-03, 0.00000E+00, 7.62305E-04, 0.00000E+00/
      DATA PN2/
     * -3.52015E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-9.52550E-04,
     *  8.56253E-04, 4.33114E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.21223E-03,
     *  2.38694E-04, 9.15245E-04, 1.28385E-03, 8.67668E-04,-5.61425E-06,
     *  1.04445E+00, 3.41112E+01, 0.00000E+00,-8.40704E-01,-2.39639E+02,
     *  7.06668E-01,-2.05873E+01,-3.63696E-01, 2.39245E+01, 1.00000E+01,
     * -1.06657E-03,-7.67292E-04, 1.54534E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         TN1(5) TN2(1)
      DATA PO1/
     *  9.99368E-01, 4.33893E-02,-2.07009E-03, 1.09617E-03, 1.05440E-03,
     *  4.83408E-04, 9.77040E-04, 9.24791E-04, 4.80247E-04, 4.94737E-02,
     *  1.05985E-03, 0.00000E+00, 0.00000E+00, 2.74409E+00, 0.00000E+00,
     * -4.96656E-03,-1.51684E-02, 4.65158E+01,-7.51133E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 6.63808E-04, 1.32397E-01,
     *  2.13315E-01,-2.06652E-03,-6.32046E-03, 0.00000E+00, 0.00000E+00,
     *  5.94545E-03,-1.90958E+02, 0.00000E+00,-4.16892E-03, 0.00000E+00,
     * -1.67499E-02, 0.00000E+00, 2.58987E-03, 5.97781E+02, 0.00000E+00,
     *  0.00000E+00, 4.44890E-04, 4.66444E-04, 1.13139E-01, 1.69134E-01,
     *  0.00000E+00, 7.11360E-04, 1.32186E-02, 2.23948E-03, 0.00000E+00/
      DATA PO2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.60571E-03,
     *  6.28078E-04, 5.05469E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.57829E-03,
     * -4.00855E-04, 5.04077E-05,-1.39001E-03,-2.33406E-03,-4.81197E-04,
     *  1.46758E+00, 6.20332E+00, 0.00000E+00, 3.66476E-01,-6.19760E+01,
     *  3.09198E-01,-1.98999E+01, 0.00000E+00,-3.29933E+02, 0.00000E+00,
     * -1.10080E-03,-9.39310E-05, 1.39638E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN2(2)
      DATA PP1/
     *  9.81637E-01,-1.41317E-03, 3.87323E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-3.58707E-02,
     * -8.63658E-03, 0.00000E+00, 0.00000E+00,-2.02226E+00, 0.00000E+00,
     * -8.69424E-03,-1.91397E-02, 8.76779E+01, 4.52188E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-7.07572E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -4.11210E-03, 3.50060E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  2.23760E-02, 0.00000E+00,-8.36657E-03, 1.61347E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.45130E-02, 0.00000E+00, 0.00000E+00/
      DATA PP2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.24152E-03,
     *  6.43365E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.33255E-03,
     *  2.42657E-03, 1.60666E-03,-1.85728E-03,-1.46874E-03,-4.79163E-06,
     *  1.22464E+00, 3.53510E+01, 0.00000E+00, 4.49223E-01,-4.77466E+01,
     *  4.70681E-01, 8.41861E+00,-2.88198E-01, 1.67854E+02, 0.00000E+00,
     *  7.11493E-04, 6.05601E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN2(3)
      DATA PQ1/
     *  1.00422E+00,-7.11212E-03, 5.24480E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-5.28914E-02,
     * -2.41301E-02, 0.00000E+00, 0.00000E+00,-2.12219E+01, 0.00000E+00,
     * -3.28077E-03, 1.65727E-02, 1.68564E+00,-6.68154E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 8.42365E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-4.34645E-03,-1.03830E-02,-8.08279E-03, 2.16780E-02,
     *  0.00000E+00,-1.38459E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.45155E-02, 0.00000E+00, 7.04573E-03,-4.73204E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.08767E-02, 0.00000E+00, 0.00000E+00/
      DATA PQ2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.21769E-04,
     * -2.27387E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.26769E-03,
     *  3.16901E-03, 4.60316E-04,-1.01431E-04, 1.02131E-03, 9.96601E-04,
     *  1.25707E+00, 2.50114E+01, 0.00000E+00, 4.24472E-01,-2.77655E+01,
     *  3.44625E-01, 2.75412E+01, 0.00000E+00, 7.94251E+02, 0.00000E+00,
     *  2.45835E-03, 1.38871E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN2(4) TN3(1)
      DATA PR1/
     *  1.01890E+00,-2.46603E-02, 1.00078E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-6.70977E-02,
     * -4.02286E-02, 0.00000E+00, 0.00000E+00,-2.29466E+01, 0.00000E+00,
     *  2.26580E-03, 2.63931E-02, 3.72625E+01,-6.39041E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.85291E-03,-7.47019E-03,-7.07265E-03, 0.00000E+00,
     *  0.00000E+00, 1.39717E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  9.58383E-03, 0.00000E+00, 9.19771E-03,-3.69121E+02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.57067E-02, 0.00000E+00, 0.00000E+00/
      DATA PR2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.92953E-03,
     * -2.77739E-03,-4.40092E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.47280E-03,
     *  2.95035E-04,-1.81246E-03, 2.81945E-03, 4.27296E-03, 9.78863E-04,
     *  1.40545E+00,-6.19173E+00, 0.00000E+00, 0.00000E+00,-7.93632E+01,
     *  4.44643E-01,-4.03085E+02, 0.00000E+00, 1.15603E+01, 0.00000E+00,
     *  2.25068E-03, 8.48557E-04,-2.98493E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN3(2)
      DATA PS1/
     *  9.75801E-01, 3.80680E-02,-3.05198E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.85575E-02,
     *  5.04057E-02, 0.00000E+00, 0.00000E+00,-1.76046E+02, 0.00000E+00,
     * -1.48297E-03,-3.68560E-03, 3.02185E+01,-3.23338E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.15558E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 4.89620E-03, 1.44594E-02, 9.91215E-03,-1.00616E-02,
     * -8.21324E-03,-1.57757E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.53569E-02, 0.00000E+00, 6.63564E-03, 4.58410E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-2.51280E-02, 0.00000E+00, 0.00000E+00/
      DATA PS2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-8.73148E-04,
     * -1.29648E-03,-7.32026E-05, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-4.68110E-03,
     * -4.66003E-03,-1.31567E-03,-7.39390E-04, 6.32499E-04,-4.65588E-04,
     * -1.29785E+00,-1.57139E+02, 0.00000E+00, 2.58350E-01,-3.69453E+01,
     *  4.10672E-01, 9.78196E+00,-1.52064E-01,-3.85084E+03, 0.00000E+00,
     * -8.52706E-04,-1.40945E-03,-7.26786E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN3(3)
      DATA PU1/
     *  9.60722E-01, 7.03757E-02,-3.00266E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.22671E-02,
     *  4.10423E-02, 0.00000E+00, 0.00000E+00,-1.63070E+02, 0.00000E+00,
     *  5.40747E-04, 7.79481E-03, 1.44908E+02, 1.51484E-04, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.41844E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 5.77884E-03, 1.06073E-02, 5.36685E-03, 9.74319E-03,
     *  0.00000E+00,-2.88015E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  1.97547E-02, 0.00000E+00,-4.44902E-03,-2.92760E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 2.34419E-02, 0.00000E+00, 0.00000E+00/
      DATA PU2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-4.65325E-04,
     * -5.50628E-04, 3.31465E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.06179E-03,
     * -3.08575E-03,-7.93589E-04,-1.08629E-04, 5.95511E-04,-9.05050E-04,
     *  1.18997E+00, 4.15924E+01, 0.00000E+00,-4.72064E-01,-9.47150E+02,
     *  3.98723E-01, 1.98304E+01, 0.00000E+00, 3.73219E+03, 0.00000E+00,
     * -1.50040E-03,-1.14933E-03,-1.56769E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN3(4)
      DATA PV1/
     *  1.03123E+00,-7.05124E-02, 8.71615E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-3.82621E-02,
     * -9.80975E-03, 0.00000E+00, 0.00000E+00, 2.89286E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 8.66153E+01, 7.91938E-04, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 4.68917E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 7.86638E-03, 9.57341E-03, 5.72268E-03, 9.90827E-03,
     *  0.00000E+00, 6.55573E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00,-4.00200E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 7.07457E-03, 0.00000E+00, 0.00000E+00/
      DATA PV2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.04970E-04,
     *  1.21560E-03,-8.05579E-06, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.49941E-03,
     * -4.57256E-04,-1.59311E-04, 2.96481E-04,-1.77318E-03,-6.37918E-04,
     *  1.02395E+00, 1.28172E+01, 0.00000E+00, 1.49903E-01,-2.63818E+01,
     *  0.00000E+00, 4.70628E+01,-2.22139E-01, 4.82292E-02, 0.00000E+00,
     * -8.67075E-04,-5.86479E-04, 5.32462E-04, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TN3(5) SURFACE TEMP TSL
      DATA PW1/
     *  1.00828E+00,-9.10404E-02,-2.26549E-02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-2.32420E-02,
     * -9.08925E-03, 0.00000E+00, 0.00000E+00, 3.36105E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.24957E+01,-5.87939E-03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.79765E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.01237E+03, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-1.75553E-02, 0.00000E+00, 0.00000E+00/
      DATA PW2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.29699E-03,
     *  1.26659E-03, 2.68402E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.17894E-03,
     *  1.48746E-03, 1.06478E-04, 1.34743E-04,-2.20939E-03,-6.23523E-04,
     *  6.36539E-01, 1.13621E+01, 0.00000E+00,-3.93777E-01, 2.38687E+03,
     *  0.00000E+00, 6.61865E+02,-1.21434E-01, 9.27608E+00, 0.00000E+00,
     *  1.68478E-04, 1.24892E-03, 1.71345E-03, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TGN3(2) SURFACE GRAD TSLG
      DATA PX1/
     *  1.57293E+00,-6.78400E-01, 6.47500E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-7.62974E-02,
     * -3.60423E-01, 0.00000E+00, 0.00000E+00, 1.28358E+02, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 4.68038E+01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.67898E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.90994E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 3.15706E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PX2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TGN2(1) TGN1(2)
      DATA PY1/
     *  8.66492E-01, 3.55807E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.12111E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 1.82458E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 1.01024E+02, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 6.54251E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      DATA PY2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-1.56959E-02,
     *  1.91001E-02, 3.15971E-02, 1.00982E-02,-6.71565E-03, 2.57693E-03,
     *  1.38692E+00, 2.82132E-01, 0.00000E+00, 0.00000E+00, 3.81511E+02,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C          TGN3(1) TGN2(2)
      DATA PZ1/
     *  1.06029E+00,-5.25231E-02, 3.73034E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.31072E-02,
     * -3.88409E-01, 0.00000E+00, 0.00000E+00,-1.65295E+02, 0.00000E+00,
     * -4.38916E-02,-3.22716E-01,-8.82393E+01, 1.18458E-01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00,-1.19782E-01,-2.13801E-01, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 2.62229E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     * -4.35863E-01, 0.00000E+00, 0.00000E+00,-5.37443E+01, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00,-4.55788E-01, 0.00000E+00, 0.00000E+00/
      DATA PZ2/
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 3.84009E-02,
     *  3.96733E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 5.05494E-02,
     *  7.39617E-02, 1.92200E-02,-8.46151E-03,-1.34244E-02, 1.96338E-02,
     *  1.50421E+00, 1.88368E+01, 0.00000E+00, 0.00000E+00,-5.13114E+01,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  5.11923E-02, 3.61225E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     *  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
C         MIDDLE ATMOSPHERE AVERAGES
      DATA PAVGM/
     M  2.61000E+02, 2.64000E+02, 2.29000E+02, 2.17000E+02, 2.17000E+02,
     M  2.23000E+02, 2.86760E+02,-2.93940E+00, 2.50000E+00, 0.00000E+00/
      END

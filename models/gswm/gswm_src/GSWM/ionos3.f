      SUBROUTINE IONOS3(Z,VLAT,COLAT,GMLT,XDIP,
     1   XMW,RHO,V0,V1P,V1M,V2P,V2M)
C
C   CALCULATES MEAN,DIURNAL, AND SEMIDIURNAL COMPONENTS
C   OF ION DRAG; CALLED FROM ABCR1 FOR USE IN INSITU PROGRAM
C
C   This is a simplified version of IONOS2 to save computational time for
C   applications emphasizing the region below 100 km.  Here geomagnetic
C   latitude is the same as geographic latitude and the magnetic dip
C   angle corresponds to a simple (non-tilted) dipole field.

C   INPUT PARAMETERS:
C
C      Z - ALTITUDE(KM)
C      VLAT - GEOGRAPHIC LATITUDE(RAD)
C      COLAT - GEOGRAPHIC COLATITUDE(RAD)
C      XMW - MOLECULAR WEIGHT(A.M.U.)
C      RHO - MEAN MASS DENSITY (KG M-3)
C
C   CALCULATIONS ARE IN CGS UNITS; FOR REFERENCE SEE---
C
C          FORBES AND GARRETT, REV. GEO. SPACE PHYS,17,1979
C          CHIU, JATP, 37, 1975.
C
      DIMENSION QI(3),F(13),A(3),B(3)

      COMPLEX V0,V1P,V1M,V2P,V2M

      COMMON/SEZUN/SWE(3),NSWE
      COMMON /MODEi/NZONAL,MOIS,NSS
      COMMON /MODEr/PERIOD,FREQ,FLUX
C     COMMON/SFLUX/FLUX
C
C   THE FOLLOWING VALUES FOR ZURICH SUNSPOT NUMBER SHOULD BE
C   REEVALUATED FOR FLUX VALUES LESS THAN ABOUT 100.
C  (SEE JACCHIA 1977, P.21 AND NATURE 4-66 P.289)
C
      DATA XNA,OMEGI/6.023E+23,200./

	pi=acos(-1.)

      V0=0.0
      V1P=(0.,0.)
      V1M=(0.,0.)
      V2P=(0.,0.)
      V2M=(0.,0.)

      IF(Z.LT.100.) GO TO 10

      DIP=XDIP
C     write *,"FLUX in ionos3= ", FLUX
      RZUR=1.08*FLUX-62.
C     write *,"RZUR in ionos3= ", RZUR
      XN=.001*RHO*XNA/XMW
      VIN=XN*2.6E-09/SQRT(XMW)
      FACT=(VIN/XN)*(1./(1.+(VIN/OMEGI)**2))
      TMO=(SWE(NSWE)/.39+1.)*3.5

      DL=PI/13.
      DO 1 LST=1,25,2
      HANG=(LST-1)*DL
      RLGM=285.*PI/180.
      RLT=VLAT
      RLTM=VLAT
      CALL IONDEN2(QTOT,QI,Z,RZUR,HANG,TMO,RLT,RLTM,RLGM,DIP)
      LS=(LST-1)/2+1
    1 F(LS)=QTOT

      CALL FORIT(F,6,2,A,B,IER)

      XN10=A(1)
      X1C=-A(2)/2.
      X1S=B(2)/2.
      X2C=A(3)/2.
      X2S=B(3)/2.
      V0=FACT*XN10*1.0E+05
      V1P=CMPLX(X1C,X1S)*FACT*1.0E+5
      V1M=CONJG(V1P)
      V2P=CMPLX(X2C,X2S)*FACT*1.0E+05
      V2M=CONJG(V2P)
  10  RETURN
      END

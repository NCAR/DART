      SUBROUTINE BFUNC(WJ,IJ,IBFUNC)
      COMPLEX WJ(91)

      COMMON/LATVECT/DTHETA,CLATT(91),XLAT(91),SNCLAT(91),CSCLAT(91),
     1    TNCLAT(91),CTNCLAT(91),SIN2I(91),GMLAT(91),DIP(91)

      DO 2 I=1,IJ
      WJ(I)=(0.0,0.0)
2     CONTINUE

      IF(IBFUNC.NE.1) go to 99

      pi=acos(-1.)
      pio2=pi/2.

      DO 1 I=1,IJ
      X=CLATT(I)
      Y=XLAT(I)
      SINY=SIN(Y)
      COSY=COS(Y)
C  Constant with Latitude in Northern Hemisphere; no forcing in Southern Hemis.
      if(clatt(i).lt.pio2) rwj=500.
c     if(clatt(i).gt.pio2) rwj=-300.
C  Approximation to (1,3) shape due to Goretzki and Rose,
C  Beitr. Phys. Atmosph., 58, 74, 1985.
c     rwj=.125*COSY*(112.*SINY**4-33.*SINY**2-3.)
C  N. Hemisphere asymmetric function
C     ARG=-2.*X
C     rwj=X*EXP(ARG)
C  Equatorial symmetric function
C     ARG=-4.*((X-1.5708)**2)
C     rwj=EXP(ARG)
C  Symmetric function peaking at high latitudes
C     rwj=SIN(X)
C     XD=X*180./3.14159
C     IF(XD.GT.60.AND.XD.LT.120.) rwj=0.0

      WJ(I)=CMPLX(rwj,0.0)
1     CONTINUE

99    RETURN
      END


      SUBROUTINE HEATpw(Z,HJ,IJ,IHEAT)
      COMPLEX HJ(91)

      COMMON/LATVECT/DTHETA,CLATT(91),XLAT(91),SNCLAT(91),CSCLAT(91),
     1    TNCLAT(91),CTNCLAT(91),SIN2I(91),GMLAT(91),DIP(91)

      DO 2 I=1,IJ
       HJ(I)=(0.0,0.0)
2     CONTINUE

      IF(IHEAT.NE.1) go to 99

      PI=ACOS(-1.)
      IF(Z.GT.20.) GO TO 99
      FH=1.0*EXP(-Z/10.)

      DO 1 I=1,IJ
      X=CLATT(I)
      Y=XLAT(I)
      SINY=SIN(Y)
      COSY=COS(Y)
      SINx=SIN(x)
      COSx=COS(x)
C  Approximation to (1,3) shape due to Goretzki and Rose,
C  Beitr. Phys. Atmosph., 58, 74, 1985.
c     rhj=FH*.125*COSY*(112.*SINY**4-33.*SINY**2-3.)
c     rhj=5.*rhj 
C  N. Hemisphere asymmetric heating function
C     ARG=-2.*X
C     rhj=FH*X*EXP(ARG)
C  Equatorial symmetric heating function
C     ARG=-4.*((X-1.5708)**2)
C     rhj=FH*EXP(ARG)
C  Symmetric heating function peaking at high latitudes
C     rhj=FH*SIN(X)
C     XD=X*180./3.14159
C     IF(XD.GT.60.AND.XD.LT.120.) rhj=0.0
C  Symmetric shape like diurnal H2O heating
C     rhj=SIN(X)*SIN(X)*FH
C  Asymmetric (1,-1) shape
C     rhj=SIN(X)*COS(X)*FH*10.
C  Approximation to (1,1) mode shape
c     YE=(Y*180./PI)/36.
C     rhj=COS(5.*Y)*EXP(-YE*YE)
C asymmetric step-function for forcing of 2-day wave
C	if(y.lt.0.) hj(i)=.5
C	if(y.gt.0.) hj(i)=-.5
C	if(y.eq.0.) hj(i)=0.
C Approximation to 2-day wave;(3,0) shape from Salby: Part II 1981 Fig 11b
	rhj=.5*sinx**3.*cosx
       HJ(I)=CMPLX(rhj,0.0)
1     CONTINUE

99    RETURN
      END

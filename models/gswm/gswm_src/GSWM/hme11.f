      SUBROUTINE hme11(xh,mois,z,HJ,IJ)

c Artificial heat source to do (1,1) HME Calculation--M. Hagan (12/27/95)
c   hj scaled to be of the order of Groves (1,1) H2O source in J/kg/s...
C			........so multiply by rho in abcr1_2.f (J/m3/s)

      real H(3,5,45),T(45),C(91), H11(91), hj(91)
      real QQR(21),XN(21),S1(4),S2(4),QR(4,21)

C latvect arrays are in colat from the north pole to the south pole
      COMMON/LATVECT/DTHETA,CLATT(91),XLAT(91),SNCLAT(91),CSCLAT(91),
     1    TNCLAT(91),CTNCLAT(91),SIN2I(91),GMLAT(91),DIP(91)

      DATA (H(1,1,J),J=1,45)/-.0007,-.0023,-.0047,-.0087,-.0153,-.0264,
     1 -.0445,-.0730,-.1158,-.1765,-.2568,-.3536,-.4561,-.5429,-.5824,
     2 -.5373,-.3747,-.0806, .3268, .7912,1.2262,1.5372,1.6502,1.5372,
     3 1.2262,.7912,.3268,-.0806,-.3747,-.5373,-.5824,-.5429,-.4561,
     4 -.3536,-.2568,-.1765,-.1158,-.0730,-.0445,-.0264,-.0153,-.0087,
     5 -.0047,-.0023,-.0007/

C Groves' (1982) Water Vapour Forcing:
C     Coefficients are given in units of mW/kg as a function of scale height:
C
C     Coefficients are specified for 0-deg longitude and reckoned to 0GMT
C       must invoke a sign change for consistency in our model where fields
C       are reckoned to noon. No long. differences for migrating modes
C
C MIGRATING DIURNAL COMPONENTS: ONLY REAL COEFFICIENTS
C
C MODE (1,1,1)
C
C  January:
      DATA (QR(1,J),J=1,21)/-2.70,-2.60,-2.68,-3.15,-3.59,-3.83,
     1 -4.36,-5.39,-5.86,-5.90,-5.86,-6.55,-6.49,-4.86,-2.95,-1.59,
     1 -.83,-.44,-.24,-.14,-.09/
C  April:
C     DATA (QR(2,J),J=1,21)/-2.70,-2.65,-2.60,-2.99,-3.35,-3.71,
C    1 -4.44,-5.88,-6.36,-6.29,-6.30,-7.14,-6.98,-5.19,-3.16,-1.83,
C    1 -1.01,-.52,-.27,-.15,-.10/
C  July:
      DATA (QR(3,J),J=1,21)/-1.92,-2.13,-2.33,-2.81,-3.22,-3.49,
     1 -4.00,-4.90,-5.30,-5.33,-5.44,-6.21,-6.05,-4.37,-2.60,-1.40,
     1 -.74,-.40,-.22,-.14,-.09/
C  October:
      DATA (QR(4,J),J=1,21)/-2.06,-2.04,-2.17,-2.81,-3.44,-3.73,
     1 -4.26,-5.34,-6.09,-6.24,-6.32,-7.41,-7.73,-5.90,-3.58,-1.92,
     1 -.99,-.51,-.27,-.16,-.10/

CTEST March 21 CCM3 data.  Interpolated M. Hagan 4/16/99

      data(QR(2,J),J=1,21)/ -2.98051, -2.98051, -2.98051, -2.82990,
     +   -2.52867, -2.22744, -2.41716, -2.60689, -2.97693, -3.52728,
     +   -4.07764, -3.18989, -2.30214, -1.55664,-0.953408,-0.350172,
     +   -0.120518, 0.109135, 0.164191,0.0446508,-0.0748907/
C

CTEST October CCM3 data.  Interpolated J. Hackney 5/6/98

CTEST      data(QR(4,J),J=1,21)/-9.58333,-9.58333,-9.58333,-9.64758,
CTEST     +-9.77608,-9.90458,-10.2066,-10.5087,-10.5876,-10.4433,-10.2989,
CTEST     +-8.32788,-6.35685,-5.01746,-4.30972,-3.60197,-2.68903,-1.77609,
CTEST     +-1.41421,-1.60340,-1.79259/
C

101   FORMAT(1X,6HCOLAT=,1X,7F10.2,/,8X,7F10.2,/,8X,7F10.2,
     +      /,8X,7F10.2,/,8X,7F10.2,/,8X,7F10.2,/,8X,7F10.2,
     +		/,8X,7F10.2,/,8X,3F10.2)
102   FORMAT(/)
100   FORMAT(3X, 4HI11=,1X,7F10.2,/,8X,7F10.2,/,8X,7F10.2,
     +      /,8X,7F10.2,/,8X,7F10.2,/,8X,7F10.2,/,8X,7F10.2,
     +		/,8X,7F10.2,/,8X,3F10.2)
103   FORMAT(3X, 4H HJ=,1X,7e10.3,/,8X,7e10.3,/,8X,7e10.3,
     +      /,8X,7e10.3,/,8X,7e10.3,/,8X,7e10.3,/,8X,7e10.3,
     +		/,8X,7e10.3,/,8X,3e10.3)

	do 99 i=1,91
	hj(i)=0.
	H11(i)=0.
   99 continue
	if(z.gt.20.) go to 12

	pi=acos(-1.)

	iseas=(mois+2)/3

      DO 1 I=1,45
1     T(I)=4.*FLOAT(I-1)+2.

      DO 2 I=1,IJ
2     C(I)=180.*CLATT(I)/pi
 110  CALL HOFFI(H,1,1,IJ,C,T,H11)

      if(z.eq.0.) then 
      PRINT 101, (C(I),I=1,IJ)
      PRINT 102
      PRINT 100, (H11(I),I=1,IJ)
      PRINT 102
      endif

      DO 3 I=1,21
      QQR(I)=QR(iseas,I)
3     XN(I)=FLOAT(I-1)*.1

      CALL ATSM(Xh,XN,QQR,21,1,S1,S2,4)
      CALL ALI(Xh,S1,S2,FH,4,1.0E-03,IER)

C     Convert units here to
C       J/kg/s (W/kg) from mW/kg (multiply by .001)     (M.Hagan --11/4/93)
C       (N.B. mW==>milliWatt and Watt=Newton*meter/sec=Joule/sec)
      do 98 i=1,IJ
      hj(i)=.001*FH*H11(i)
C Introduce sign change due to noon/midnight convention discussed
C       in comments above:
      hj(i)=hj(i)*(-1.)
   98 continue
      if(z.le.10.) PRINT 103, (Hj(I),I=1,IJ)
      if(z.le.10.) PRINT 102
   12 continue

      RETURN
      END

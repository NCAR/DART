	subroutine honglate(z,hlj,ij)
C  adapted from latent.f (M. Hagan 10/25/95)
c----- This subroutine calculates the latitude/height distribution of latent
c----- heating (UNITS = W/Kg or Joules/sec/Kg) corresponding to the semidiurnal 
c----- tide.

c----- The GSWM, due to historical reasons, defines the heating which maximizes
c----- at noon to be a real number.  That is, the forcing is given by

c			A cos w(t-tn)

c----- where w = frequency and t = hrs reckoned from local noon, and tn = time
c----- of maximum forcing reckoned from local noon.  In this notation, the

      REAL HLJ(91)

      COMMON/LATVECT/DTHETA,CLATT(91),XLAT(91),SNCLAT(91),CSCLAT(91),
     1    TNCLAT(91),CTNCLAT(91),SIN2I(91),GMLAT(91),DIP(91)

      COMMON /MODEi/NZONAL,MOIS,NSS
      COMMON /MODEr/PERIOD,FREQ,FLUX

      COMMON/HOFFS/H11(91),H12(91),H13(91),H1M2(91),H1M1(91)
     1,H1M4(91),H22(91),H23(91),H24(91),H25(91),H26(91)
     + ,H33(91),H34(91),H35(91),H36(91)

      DO 2 I=1,IJ
      HLJ(I)=0.0
2     CONTINUE
      IF(Z.GT.20.) GO TO 99
c     go to 99

      PI=ACOS(-1.)
      TPI=2.*PI

c------------------------------------------------------------------------------
c  The following expression for the vertical variation of latent heating due to
c  cumulus convection is from Hong and Wang (1980, Bull. Geophys., 19, 56-84),
c  which is based on the work of Reed and Recker (1971) and Nitta (1972).

      FH=.00854*(EXP(-(Z-6.5)*(Z-6.5)/29.05)-0.23*EXP(-Z/1.31))

c------------------------------------------------------------------------------
      tpart=COS(TPI*2.10/12.)
      DO 1 I=1,IJ
      
      HLJ(I)=FH*tpart*H22(I)

1     CONTINUE

99    RETURN
      END


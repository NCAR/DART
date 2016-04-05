      SUBROUTINE NORTHV(ZZ,THETA,VN)
      COMMON /MODEi/NZONAL,MOIS,NSS
      COMMON /MODEr/PERIOD,FREQ,FLUX
      COMMON/SEZUN/SWE(3),NSWE
C
      PI=acos(-1.)
C
      TT=(PI/2.)-THETA
      T2=SIN(2.*TT)
      T1=SIN(TT)
      F1=(((TT/PI/2.)**3.)*T2*T2/.36)*(5.+(NSS-1)*15.)
      F1=-F1*1.7
      F2=0.0
      IF(NSWE.EQ.3) GO TO 1
      F1=F1*1.5
      F2=COS(TT)*(55.+(NSS-1)*15.)
    1 F3=(1.-EXP(-(ZZ-120.)/125.))/.8
      IF (ZZ.LT.120.) F3=0.0
      VN=(F1+F2)*F3
      RETURN
      END


      SUBROUTINE AMPHZgeo(F,ij)
C Version for non-migrating tide
C     CONVERTS COMPLEX SOLUTION TO (AMP,PHASE) WHERE PHASE IS ***GMT***
      COMPLEX F(*)
      COMMON/MODE/NZONAL,PERIOD,FREQ,MOIS,NSS,FLUX
      
      pi=acos(-1.)
      tpi=2.*pi
      omega=tpi/(24.*period)

      DO 1 I=1,ij
      Y=AIMAG(F(I))
      X=REAL(F(I))+1.0E-08
      PHASE=-ATAN2(Y,X)/omega
	If(phase.lt.0.) phase=phase+(24.*period)
	If(phase.gt.24.) phase=phase-(24.*period)
      AMP=CABS(F(I))
1     F(I)=CMPLX(AMP,PHASE)
      RETURN
      END

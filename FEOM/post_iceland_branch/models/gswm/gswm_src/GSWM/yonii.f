
      FUNCTION YONII(RLTM,RF,R,PHI,TMO,DEC,PI,CLTM,SLTM)
      B=1.3+(.139*(1.+COS(RLTM-PI/4.))+.0517*R)*R*R
      DRF=1./RF
      W1=PI/6.
      W2=W1+W1
      DE=.1778*R*R
      ALTM=ABS(RLTM)
      SNX=SIN(ALTM-.5236)
      AE=.2*(1.-SNX)
      BLTM=ABS(ALTM-PI/9.)
      SX=SIN(BLTM)
      FE=.13-.06*SX
      YM=COS(RLTM+DEC)
      CPHG=COS(PHI)
      XTC=YM**3*(1.-CPHG)**.25
      YTC=-(.15+.3*SIN(ALTM))*XTC
      T1=AE*(1.+.6*COS(W2*(TMO-4.)))*COS(W1*(TMO-1.))
      TRIV=(COS(RLTM-W1))*(COS(W1*(.5*TMO-1.)))**3
      TRIV=TRIV+(COS(RLTM+PI/4.))*(COS(W1*(.5*TMO-4.)))**2
      QQ=1.+.085*TRIV
      T2=.7*(QQ+DE*DRF*COS(W2*(TMO-4.3)))*W(B,RLTM,PHI,DEC)
      T3=FE*COS(W2*(TMO-4.5))+YTC
      YONII=(T1+T3)*DRF+T2
      RETURN
      END

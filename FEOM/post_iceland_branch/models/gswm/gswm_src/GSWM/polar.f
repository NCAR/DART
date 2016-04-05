
      FUNCTION POLAR(RLTM,RLGM,DIP,R,PHI,TMO,DEC,PI,CLTM,SLTM)
      T=PI*TMO/12.
      V=SIN(T)
      U=COS(T+T)
      Y=SIN(RLGM/2.)
      YS=COS(RLGM/2.-PI/20.)
      Z=SIN(RLGM)
      ZA=SQRT(ABS(Z))
      AM=1.+V
      IF(RLTM) 1,2,2
 2    C=-23.5*PI/180.
      POLAR=(2.+1.2*R)*W(1.2,RLTM,PHI,C)*(1.+.3*V)
      GOTO 3
 1    B=V*(.5*Y-.5*Z-Y**8)-AM*U*(Z/ZA)*EXP(-4.*Y*Y)
      POLAR=2.5+2.0*R+U*(0.5+(1.3+.2*R)*YS**4)
      POLAR=POLAR+(1.3+0.5*R)*COS(PHI-PI*(1.+B))
      POLAR=POLAR*(1.+0.4*(1.-V*V))*EXP(-1.0*V*YS**4)
 3    CONTINUE
      RETURN
      END

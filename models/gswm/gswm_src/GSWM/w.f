
      FUNCTION W(B,XI,ETA,DEC)
      P=PSI(XI,ETA,DEC)
      W=EXP(-B*(COS(P)-COS(XI)))
      RETURN
      END

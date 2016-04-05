C---------------------------------------------------------------------- 
C  CONVERT AN 8 DIGIT INTEGER DATE (YYMMDDHH) TO 10 DIGITS (CCYYMMDDHH)
C---------------------------------------------------------------------- 
      FUNCTION I4DY(IDATE)                                              

      IF(IDATE.LT.10**8) THEN
         IY = IDATE/10**6
         IF(IY.GT.20) I4DY = IDATE + 19*10**8 
         IF(IY.LE.20) I4DY = IDATE + 20*10**8 
      ELSE
         I4DY = IDATE
      ENDIF

      RETURN
      END

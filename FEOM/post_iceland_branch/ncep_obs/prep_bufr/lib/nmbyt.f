C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION NMBYT(LUNIT)
 
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      NMBYT = -1
 
C  CHECK THE FILE STATUS
C  ---------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      IF(IM.EQ.0) GOTO 902
      NMBYT = IUPB(MBAY(1,LUN),5,24)


      RETURN
 
C  ERROR EXITS
C  -----------
 
900   CALL BORT('NMBYT - FILE IS CLOSED                   ')
901   CALL BORT('NMBYT - FILE IS OPEN FOR OUTPUT          ')
902   CALL BORT('NMBYT - NO MESSAGE IS OPENED             ')
      END

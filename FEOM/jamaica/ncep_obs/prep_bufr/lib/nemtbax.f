C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NEMTBAX(LUN,NEMO,MTYP,MSBT,INOD)
 
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),
     .                TABA(50,32),TABB(250,32),TABD(250,32)
 
      CHARACTER*(*) NEMO
      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*20  NEMT
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      NEMT = NEMO
      INOD = 0
 
C  LOOK FOR NEMO IN TABLE A
C  ------------------------
 
      DO I=1,NTBA(LUN)
      IF(TABA(I,LUN)(4:11).EQ.NEMO) THEN
         MTYP = IDNA(I,LUN,1)
         MSBT = IDNA(I,LUN,2)
         INOD = MTAB(I,LUN)
         IF(MTYP.LT.0 .OR. MTYP.GT.255) GOTO 900
         IF(MSBT.LT.0 .OR. MSBT.GT.255) GOTO 901
         RETURN
      ENDIF
      ENDDO
 
      RETURN
900   CALL BORT('NEMTBAX - BAD MTYP  '//NEMT)
901   CALL BORT('NEMTBAX - BAD MSBT  '//NEMT)
      END

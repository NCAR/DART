C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READSB(LUNIT,IRET)
 
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
      COMMON /UNPTYP/ MSGUNP(32)
 
C-----------------------------------------------------------------------
      ENTRY READERS(LUNIT,IRET)
C-----------------------------------------------------------------------
 
      IRET = 0
 
C  CHECK THE FILE STATUS
C  ---------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      IF(IM.EQ.0) THEN
         IRET = -1
         RETURN
      ENDIF
 
C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE
C  ---------------------------------------------
 
      IF(NSUB(LUN).EQ.MSUB(LUN)) THEN
         IRET = -1
         RETURN
      ELSE
         NSUB(LUN) = NSUB(LUN) + 1
      ENDIF
 
C  READ THE NEXT SUBSET AND RESET THE POINTERS
C  -------------------------------------------
 
      IF(MSGUNP(LUN).EQ.0) THEN
         IBIT = MBYT(LUN)*8
         CALL UPB(NBYT,16,MBAY(1,LUN),IBIT)
         CALL RDTREE(LUN)
         MBYT(LUN) = MBYT(LUN) + NBYT
      ELSEIF(MSGUNP(LUN).EQ.1) THEN
         IBIT = MBYT(LUN)
         CALL RDTREE(LUN)
         MBYT(LUN) = IBIT
      ELSEIF(MSGUNP(LUN).EQ.2) THEN
         CALL RDCMPS(LUN)
      ELSE
         GOTO 903
      ENDIF
 
      RETURN
900   CALL BORT('READSB - FILE IS CLOSED             ')
901   CALL BORT('READSB - FILE IS OPEN FOR OUTPUT    ')
902   CALL BORT('READSB - NO MESSAGE OPEN            ')
903   CALL BORT('READSB - UNKNOWN MESSAGE UNPACK TYPE')
      END

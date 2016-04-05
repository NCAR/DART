C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READIBM(LUNIT,SUBSET,JDATE,IRET)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
 
      CHARACTER*8 SEC0,SUBSET
      CHARACTER*4 BUFR
      CHARACTER*1 CBAY(8*5000)
      DIMENSION   JBAY(5000)
      EQUIVALENCE (CBAY(1),JBAY(1))
      EQUIVALENCE (CBAY(1),SEC0)
 
C-----------------------------------------------------------------------
      LBMG(SEC0) = IUPM(SEC0(5:7),24)
C-----------------------------------------------------------------------
 
      IRET = 0
 
C  CHECK THE FILE STATUS
C  ---------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      CALL WTSTAT(LUNIT,LUN,IL,1)
 
C  READ A MESSAGE INTO A MESSAGE BUFFER - SKIP DX MESSAGES
C  -------------------------------------------------------
 
1     SEC0 = ' '
      READ(LUNIT,ERR=902,END=100) SEC0,(CBAY(I),I=9,LBMG(SEC0))
      DO I=1,8
      CBAY(I) = SEC0(I:I)
      ENDDO
      DO I=1,LMSG(SEC0)
      MBAY(I,LUN) = JBAY(I)
      ENDDO
 
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,IRET)
      IF(IRET.NE.0) GOTO 1
      RETURN
 
C  EOF ON ATTEMPTED READ
C  ---------------------
 
100   CALL WTSTAT(LUNIT,LUN,IL,0)
      INODE(LUN) = 0
      IDATE(LUN) = 0
      SUBSET = ' '
      JDATE = 0
      IRET = -1
      RETURN
 
C  ERROR EXITS
C  -----------
 
900   CALL BORT('READIBM - FILE IS CLOSED                       ')
901   CALL BORT('READIBM - FILE IS OPEN FOR OUTPUT              ')
902   CALL BORT('READIBM - I/O ERROR READING MESSAGE            ')
      END

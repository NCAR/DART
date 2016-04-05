C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE RDTREE(LUN)
 
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),
     .                JUMP(15000),LINK(15000),JMPB(15000),
     .                IBT(15000),IRF(15000),ISC(15000),
     .                ITP(15000),VALI(15000),KNTI(15000),
     .                ISEQ(15000,2),JSEQ(15000)
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)
      COMMON /USRBIT/ NBIT(15000),MBIT(15000)
 
      CHARACTER*10 TAG
      CHARACTER*8  CVAL
      CHARACTER*3  TYP
      DIMENSION    IVAL(15000)
      EQUIVALENCE  (CVAL,RVAL)
      REAL*8       VAL,RVAL
 
C-----------------------------------------------------------------------
CFPP$ EXPAND (UPBB)
C-----------------------------------------------------------------------
      MPS(NODE) = 2**(IBT(NODE))-1
      UPS(NODE) = float(IVAL(N)+IRF(NODE))*10.**(-ISC(NODE))
C-----------------------------------------------------------------------
 
C  CYCLE THROUGH A SUBSET SETTING UP THE USER ARRAY
C  ------------------------------------------------
 
      MBIT(1) = IBIT
      NBIT(1) = 0
      CALL RCSTPL(LUN)
 
C  UNPACK A SUBSET INTO THE USER ARRAY
C  -----------------------------------
 
      DO N=1,NVAL(LUN)
      CALL UPBB(IVAL(N),NBIT(N),MBIT(N),MBAY(1,LUN))
      ENDDO
 
C  CONVERT THE UNPACKED INTEGERS TO THE PROPER TYPES
C  -------------------------------------------------
 
      DO N=1,NVAL(LUN)
      NODE = INV(N,LUN)
      IF(ITP(NODE).EQ.1) THEN
         VAL(N,LUN) = IVAL(N)
      ELSEIF(ITP(NODE).EQ.2) THEN
         IF(IVAL(N).LT.MPS(NODE)) VAL(N,LUN) = UPS(NODE)
      ENDIF
      ENDDO
 
C  SPECIAL TREATMENT FOR CHARACTERS
C  --------------------------------
 
      DO N=1,NVAL(LUN)
      NODE = INV(N,LUN)
      IF(ITP(NODE).EQ.3) THEN
         CVAL = ' '
         CALL UPC(CVAL,NBIT(N)/8,MBAY(1,LUN),MBIT(N))
         VAL(N,LUN) = RVAL
      ENDIF
      ENDDO
 
      IBIT = NBIT(NVAL(LUN))+MBIT(NVAL(LUN))
 
      RETURN
      END

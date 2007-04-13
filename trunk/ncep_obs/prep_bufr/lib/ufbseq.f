C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE UFBSEQ(LUNIN,USR,I1,I2,IRET,STR)

      PARAMETER(MTAG=10)

      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),
     .                JUMP(15000),LINK(15000),JMPB(15000),
     .                IBT(15000),IRF(15000),ISC(15000),
     .                ITP(15000),VALI(15000),KNTI(15000),
     .                ISEQ(15000,2),JSEQ(15000)
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)

      CHARACTER*(*) STR
      CHARACTER*10  TAG,TAGS(MTAG)
      CHARACTER*3   TYP
      LOGICAL       RPC
      REAL*8        USR(I1,I2),VAL,BMISS

      DATA BMISS /10E10/

C----------------------------------------------------------------------
C----------------------------------------------------------------------

      IRET = 0

C  CHECK THE FILE STATUS AND CLEAR AN INPUT ARRAY
C  ----------------------------------------------

      LUNIT = ABS(LUNIN)
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) CALL BORT('UFBSEQ - FILE IS CLOSED ')
      IF(IM.EQ.0) CALL BORT('UFBSEQ - NO MESSAGE OPEN')

      IO = MIN(MAX(0,IL),1)
      IF(LUNIT.NE.LUNIN) IO = 0

      IF(IO.EQ.0) THEN
         DO I=1,I1*I2
         USR(I,1) = BMISS
         ENDDO
      ENDIF

C  CHECK FOR VALID SEQUENCE AND SEQUENCE LENGTH ARGUMENTS
C  ------------------------------------------------------

      CALL PARSEQ(STR,TAGS,MTAG,NTAG)
      IF(NTAG.LT.1) RETURN
      IF(NTAG.GT.1) CALL BORT('UFBSEQ - MORE THAN ONE STRING ARG!')

      INOD = INODE(LUN)
      IF(INOD.NE.INV(1,LUN)) CALL BORT('UFBSEQ - I-NODE MISMATCH!')

      DO NODE=INOD,ISC(INOD)
      IF(STR.EQ.TAG(NODE)) THEN
         IF(TYP(NODE).EQ.'SEQ' .OR. TYP(NODE).EQ.'RPC') THEN
            NODS = NODE
            DO WHILE(JUMP(NODS).NE.0 .OR.
     .               LINK(NODS).NE.0 .OR.
     .               JMPB(NODS).NE.NODE) 
            NODS = NODS+1     
            IF(NODS.GT.ISC(INOD)) CALL BORT('UFBSEQ - NODS NO END')
            ENDDO
            INS1 = INVWIN(NODE,LUN,   1,NVAL(LUN))
            INS2 = INVWIN(NODS,LUN,INS1,NVAL(LUN))
         ELSEIF(TYP(NODE).EQ.'SUB') THEN
            INS1 = 1             
            INS2 = NVAL(LUN)     
         ELSE
            CALL BORT('UFBSEQ - MNEMONIC NOT A SEQUENCE: '//TAG(NODE))
         ENDIF
         NSEQ = 0
         DO ISQ=INS1,INS2
         ITYP = ITP(INV(ISQ,LUN))
         IF(ITYP.EQ.1) CALL BORT('UFBSEQ - ILLEGAL SEQUENCE!')
         IF(ITYP.GT.1) NSEQ = NSEQ+1 
         ENDDO
         IF(NSEQ.GT.I1) CALL BORT('UFBSEQ - SEQ LENGTH GT I1!')
         GOTO 1
      ENDIF
      ENDDO

      RETURN

C  FRAME A SECTION OF THE BUFFER - RETURN WHEN NO FRAME
C  ----------------------------------------------------

1     INS1 = INVTAG(NODE,LUN,INS1,NVAL(LUN))
      IF(INS1.GT.0.AND.TYP(NODE).EQ.'RPC'.AND.VAL(INS1,LUN).EQ.0.) THEN
         INS1 = INS1+1
         GOTO 1
      ENDIF
      IF(INS1.EQ.0.AND.IO.EQ.1.AND.IRET.LT.I2) THEN
         CALL BORT('UFBSEQ - INCOMPLETE WRITE')
      ELSEIF(INS1.GT.0.AND.IO.EQ.0.AND.IRET+1.GT.I2) THEN
         CALL BORT('UFBSEQ - INCOMPLETE READ ')
      ELSEIF(INS1.EQ.0.OR.IRET.EQ.I2) THEN
         RETURN
      ENDIF

      IRET = IRET+1
      INS1 = INS1+1

C  READ/WRITE USER VALUES
C  ----------------------

      J = INS1
      DO I=1,NSEQ
      DO WHILE(ITP(INV(J,LUN)).LT.2) 
      J = J+1
      ENDDO
      IF(IO.EQ.0) USR(I,IRET) = VAL(J,LUN )
      IF(IO.EQ.1) VAL(J,LUN ) = USR(I,IRET)
      J = J+1
      ENDDO

C  CHECK FOR NEXT FRAME
C  --------------------

      GOTO 1
      END

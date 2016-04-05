      SUBROUTINE TABENT(LUN,NEMO,TAB,ITAB,IREP,IKNT,JUM0)

C************************************************************************
C* TABENT								*
C*									*
C* Given a Table B or Table D mnemonic NEMO, this subroutine builds and	*
C* stores an entry for NEMO within the internal jump/link table.	*
C*									*
C* TABENT  ( LUN, NEMO, TAB, ITAB, IREP, IKNT, JUM0 )			*
C*									*
C* Input parameters:							*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					BUFR table arrays		*
C*	NEMO		CHARACTER*(*)	Mnemonic			*
C*	TAB		CHARACTER	Internal BUFR table array	*
C*					('B' or 'D') in which NEMO	*
C*					is defined			*
C*	ITAB		INTEGER		Positional index of NEMO	*
C*					within TAB			*
C*	IREP		INTEGER		Positional index within		*
C*					COMMON /REPTAB/ arrays, for use	*
C*					when NEMO is replicated:	*
C*					  0 = NEMO is not replicated	*
C*	IKNT		INTEGER		Number of replications, for use	*
C*					when NEMO is replicated using	*
C*					F=1 regular (i.e. non-delayed)	*
C*					replication:			*
C*					  0 = NEMO is not replicated	*
C*					      using F=1 regular (i.e.	*
C*					      non-delayed) replication	*
C*	JUM0		INTEGER		Index value to be stored for	*
C*					NEMO within internal jump/link	*
C*					table array JMPB(*)		*
C*									*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************

C     Note that the values within the COMMON /REPTAB/ arrays were
C     initialized within subroutine BFRINI.

      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)

      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),
     .                JUMP(15000),LINK(15000),JMPB(15000),
     .                IBT(15000),IRF(15000),ISC(15000),
     .                ITP(15000),VALI(15000),KNTI(15000),
     .                ISEQ(15000,2),JSEQ(15000)
      COMMON /TABCCC/ ICDW,ICSC,ICRV

      CHARACTER*24 UNIT
      CHARACTER*10 TAG,RTAG
      CHARACTER*8  NEMO
      CHARACTER*3  TYP,TYPS,TYPT
      CHARACTER*1  REPS,TAB

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  MAKE A JUMP/LINK TABLE ENTRY FOR A REPLICATOR
C  ---------------------------------------------

      IF(IREP.NE.0) THEN
         RTAG = REPS(IREP,1)//NEMO
         DO I=1,10
         IF(RTAG(I:I).EQ.' ') THEN
            RTAG(I:I) = REPS(IREP,2)
            CALL INCTAB(RTAG,TYPS(IREP,1),NODE)
            JUMP(NODE) = NODE+1
            JMPB(NODE) = JUM0
            LINK(NODE) = 0
            IBT (NODE) = LENS(IREP)
            IRF (NODE) = 0
            ISC (NODE) = 0
            IF(IREP.EQ.1) IRF(NODE) = IKNT
            JUM0 = NODE
            GOTO 1
         ENDIF
         ENDDO
         GOTO 900
      ENDIF

C  MAKE AN JUMP/LINK ENTRY FOR AN ELEMENT OR A SEQUENCE
C  ----------------------------------------------------

1     IF(TAB.EQ.'B') THEN
         CALL NEMTBB(LUN,ITAB,UNIT,ISCL,IREF,IBIT)
         IF(UNIT(1:5).EQ.'CCITT') TYPT = 'CHR'
         IF(UNIT(1:5).NE.'CCITT') TYPT = 'NUM'
         CALL INCTAB(NEMO,TYPT,NODE)
         JUMP(NODE) = 0
         JMPB(NODE) = JUM0
         LINK(NODE) = 0
         IBT (NODE) = IBIT
         IRF (NODE) = IREF
         ISC (NODE) = ISCL
         IF(UNIT(1:5).EQ.'CODE') TYPT = 'COD'
         IF(UNIT(1:5).EQ.'FLAG') TYPT = 'FLG'
         IF(TYPT.EQ.'NUM') THEN
            IBT(NODE) = IBT(NODE)+ICDW
            ISC(NODE) = ISC(NODE)+ICSC
         ENDIF
      ELSEIF(TAB.EQ.'D') THEN
         IF(IREP.EQ.0) TYPT = 'SEQ'
         IF(IREP.NE.0) TYPT = TYPS(IREP,2)
         CALL INCTAB(NEMO,TYPT,NODE)
         JUMP(NODE) = NODE+1
         JMPB(NODE) = JUM0
         LINK(NODE) = 0
         IBT (NODE) = 0
         IRF (NODE) = 0
         ISC (NODE) = 0
      ELSE
         GOTO 901
      ENDIF

      RETURN
900   CALL BORT('TABENT - REPLICATOR ERROR: '//RTAG//':'//NEMO)
901   CALL BORT('TABENT - UNDEFINED TAG   : '           //NEMO)
      END

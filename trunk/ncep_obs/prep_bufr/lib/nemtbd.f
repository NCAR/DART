      SUBROUTINE NEMTBD(LUN,ITAB,NSEQ,NEMS,IRPS,KNTS)

C************************************************************************
C* NEMTBD								*
C*									*
C* Given a Table D sequence mnemonic (i.e. a "parent mnemonic"), this	*
C* subroutine returns a list of the mnemonics contained within that	*
C* sequence (i.e. the "child mnemonics").  This information should have	*
C* been packed into the internal BUFR table D entry for the parent	*
C* mnemonic via previous calls to subroutine PKTDD.  Note that this	*
C* subroutine (i.e. NEMTBD) does *not* recursively resolve child	*
C* mnemonics which are themselves Table D sequence mnemonics; rather,	*
C* such resolution must be done via separate subsequent calls to this	*
C* subroutine.								*
C*									*
C* NEMTBD  ( LUN, ITAB, NSEQ, NEMS, IRPS, KNTS )			*
C*									*
C* Input parameters:							*
C*	LUN		INTEGER		I/O stream index of parent	*
C*					mnemonic within internal BUFR	*
C*					Table D array TABD(*,*)		*
C*	ITAB		INTEGER		Positional index of parent	*
C*					mnemonic within internal BUFR	*
C*					Table D array TABD(*,*)		*
C*									*
C* Output parameters:							*
C*	NSEQ		INTEGER		Total number of child mnemonics	*
C*					for the parent mnemonic given	*
C*					by TABD(ITAB,LUN)		*
C*	NEMS(NSEQ)	CHARACTER*8	Child mnemonics			*
C*	IRPS(NSEQ)	INTEGER		Return value:			*
C*				      /               \			*
C*	      ------------------------                 ---------------	*
C* 		The interpretation of the return value IRPS(I) depends	*
C*		upon the type of descriptor corresponding to NEMS(I),	*
C*		as follows:						*
C*									*
C*		IF ( NEMS(I) corresponds to an F=1 regular		*
C*		     (i.e. non-delayed) replication descriptor ) THEN	*
C*		    IRPS(I) = 1						*
C*		ELSE IF ( NEMS(I) corresponds to a delayed replicator	*
C*			  or replication factor descriptor )  THEN	*
C*		    IRPS(I) = positional index of corresponding		*
C*			      descriptor within internal replication	*
C*			      array IDNR(*,*)				*
C*		ELSE							*
C*		    IRPS(I) = 0						*
C*		END IF							*
C*	      --------------------------------------------------------	*
C*									*
C*	KNTS(NSEQ)	INTEGER		Return value:			*
C*				      /               \			*
C*	      ------------------------                 ---------------	*
C* 		The interpretation of the return value KNTS(I) depends	*
C*		upon the type of descriptor corresponding to NEMS(I),	*
C*		as follows:						*
C*									*
C*		IF ( NEMS(I) corresponds to an F=1 regular		*
C*		     (i.e. non-delayed) replication descriptor ) THEN	*
C*		    KNTS(I) = number of replications			*
C*		ELSE							*
C*		    KNTS(I) = 0						*
C*		END IF							*
C*	      --------------------------------------------------------	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation				*
C************************************************************************

      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),
     .                TABA(50,32),TABB(250,32),TABD(250,32)

      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*8   NEMO,NEMS,NEMT,NEMF
      CHARACTER*1   TAB
      DIMENSION     NEMS(250),IRPS(250),KNTS(250)
      LOGICAL       REP

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      IF(ITAB.LE.0 .OR. ITAB.GT.NTBD(LUN)) GOTO 900

      REP  = .FALSE.

C  CLEAR THE RETURN VALUES
C  -----------------------

      NSEQ = 0

      DO I=1,250
      NEMS(I) = ' '
      IRPS(I) = 0
      KNTS(I) = 0
      ENDDO

C  PARSE THE TABLE D ENTRY
C  -----------------------

      NEMO = TABD(ITAB,LUN)(7:14)
      IDSC = IDND(ITAB,LUN)
      CALL UPTDD(ITAB,LUN,0,NDSC)

      IF(IDSC.LT.IFXY('300000')) GOTO 901
      IF(IDSC.GT.IFXY('363255')) GOTO 901
C     IF(NDSC.LE.0             ) GOTO 902


C     Loop through each child mnemonic.

      DO J=1,NDSC
      IF(NSEQ+1.GT.250) GOTO 903
      CALL UPTDD(ITAB,LUN,J,IDSC)
      CALL NUMTAB(LUN,IDSC,NEMT,TAB,IRET)
      IF(TAB.EQ.'R') THEN
         IF(REP) GOTO 904
         REP = .TRUE.
         IF(IRET.LT.0) THEN

C	    F=1 regular (i.e. non-delayed) replication.

            IRPS(NSEQ+1) = 1
            KNTS(NSEQ+1) = ABS(IRET)
         ELSEIF(IRET.GT.0) THEN

C	    Delayed replication.

            IRPS(NSEQ+1) = IRET
         ENDIF
      ELSEIF(TAB.EQ.'F') THEN

C	    Replication factor.

         IF(.NOT.REP) GOTO 904
         IRPS(NSEQ+1) = IRET
         REP = .FALSE.
      ELSEIF(TAB.EQ.'D'.OR.TAB.EQ.'C') THEN
         REP = .FALSE.
         NSEQ = NSEQ+1
         NEMS(NSEQ) = NEMT
      ELSEIF(TAB.EQ.'B') THEN
         REP = .FALSE.
         NSEQ = NSEQ+1
         IF(NEMT(1:1).EQ.'.') THEN

C	    This is a "following value" mnemonic.

            CALL UPTDD(ITAB,LUN,J+1,IDSC)
            CALL NUMTAB(LUN,IDSC,NEMF,TAB,IRET)
            CALL RSVFVM(NEMT,NEMF)
            IF(TAB.NE.'B') GOTO 906
         ENDIF
         NEMS(NSEQ) = NEMT
      ELSE
         GOTO 905
      ENDIF
      ENDDO

      RETURN
900   CALL BORT('NEMTBD - ITAB NOT IN TABLE D   '                )
901   CALL BORT('NEMTBD - BAD DESCRIPTOR VALUE: '          //NEMO)
902   CALL BORT('NEMTBD - ZERO LENGTH SEQUENCE: '          //NEMO)
903   CALL BORT('NEMTBD - TOO MANY DESCRIPTORS IN SEQ: '   //NEMO)
904   CALL BORT('NEMTBD - REPLICATOR OUT OF ORDER IN SEQ: '//NEMO)
905   CALL BORT('NEMTBD - BAD DESCRIPTOR IN SEQUENCE: '    //NEMO)
906   CALL BORT('NEMTBD - FOLLOWING VALUE NOT FROM TABLEB:'//NEMF)
      END

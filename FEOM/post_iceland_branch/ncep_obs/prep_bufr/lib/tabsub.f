      SUBROUTINE TABSUB(LUN,NEMO)

C************************************************************************
C* TABSUB								*
C*									*
C* Given a Table A mnemonic NEMO, this subroutine builds the *entire*	*
C* jump/link tree (i.e. including recursively resolving all child	*
C* mnemonics!) for NEMO within the internal jump/link table.		*
C*									*
C* TABSUB  ( LUN, NEMO )						*
C*									*
C* Input parameters:							*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					BUFR table arrays		*
C*	NEMO		CHARACTER*(*)	Table A mnemonic		*
C*									*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
C*									*
C*	EXAMPLE SHOWING CONTENTS OF INTERNAL JUMP/LINK TABLE		*
C*			(WITHIN COMMON /TABLES/)			*
C*									*
C*  INTEGER MAXTAB =							*
C*	maximum number of jump/link table entries			*
C*									*
C*  INTEGER NTAB =							*
C*	actual number of jump/link table entries currently in use	*
C*									*
C*  FOR I = 1, NTAB:							*
C*									*
C*    CHARACTER*10 TAG(I) =						*
C*	mnemonic							*
C*									*
C*    CHARACTER*3 TYP(I) =						*
C*	mnemonic type indicator:					*
C*	  "SUB" if TAG(I) is a Table A mnemonic				*
C*	  "SEQ" if TAG(I) is a Table D mnemonic using either		*
C*		short (i.e. 1-bit) delayed replication,			*
C*		F=1 regular (i.e. non-delayed) replication, or		*
C*		no replication at all					*
C*	  "RPC" if TAG(I) is a Table D mnemonic using either		*
C*		medium (i.e. 8-bit) delayed replication or		*
C*		long (i.e. 16-bit) delayed replication			*
C*	  "DRB" if TAG(I) denotes the short (i.e. 1-bit) delayed	*
C*		replication of a Table D mnemonic (which would then	*
C*		itself have its own separate entry in the jump/link	*
C*		table with a corresponding TAG value of "SEQ"!)		*
C*	  "DRP" if TAG(I) denotes either the medium (i.e. 8-bit) or	*
C*		long (i.e. 16-bit) delayed replication of a Table D	*
C*		mnemonic (which would then itself have its own separate	*
C*		entry in the jump/link table with a corresponding TAG	*
C*		value of "RPC"!)					*
C*	  "REP" if TAG(I) denotes the F=1 regular (i.e. non-delayed)	*
C*		replication of a Table D mnemonic (which would then	*
C*		itself have its own separate entry in the jump/link	*
C*		table with a corresponding TAG value of "SEQ"!)		*
C*	  "CHR" if TAG(I) is a Table B mnemonic with units "CCITT IA5"	*
C*	  "NUM" if TAG(I) is a Table B mnemonic with any units other	*
C*		than "CCITT IA5"					*
C*									*
C*    INTEGER JMPB(I):							*
C*	IF ( TYP(I) = "SUB" ) THEN					*
C*	  JMPB(I) = 0							*
C*	ELSE IF ( ( TYP(I) = "SEQ" and TAG(I) uses either		*
C*		    short (i.e. 1-bit) delayed replication or		*
C*		    F=1 regular (i.e. non-delayed) replication )	*
C*			OR						*
C*		  ( TYP(I) = "RPC" )  )  THEN				*
C*	  JMPB(I) = the index of the jump/link table entry denoting	*
C*		    the replication of TAG(I)				*
C*	ELSE								*
C*	  JMPB(I) = the index of the jump/link table entry for the	*
C*		    Table A or Table D mnemonic of which TAG(I) is a	*
C*		    child						*
C*	END IF								*
C*									*
C*    INTEGER JUMP(I):							*
C*	IF ( ( TYP(I) = "CHR" ) OR					*
C*	     ( TYP(I) = "NUM" ) ) THEN					*
C*	  JUMP(I) = 0							*
C*	ELSE IF ( ( TYP(I) = "DRB" ) OR					*
C*		  ( TYP(I) = "DRP" ) OR					*
C*		  ( TYP(I) = "REP" ) ) THEN				*
C*	  JUMP(I) = the index of the jump/link table entry for the	*
C*		    Table D mnemonic whose replication is denoted	*
C*		    by TAG(I)						*
C*	ELSE								*
C*	  JUMP(I) = the index of the jump/link table entry for the	*
C*		    Table B or Table D mnemonic which, sequentially,	*
C*		    is the first child of TAG(I)			*
C*	END IF								*
C*									*
C*    INTEGER LINK(I):							*
C*	IF ( ( TYP(I) = "SEQ" and TAG(I) uses either			*
C*		    short (i.e. 1-bit) delayed replication or		*
C*		    F=1 regular (i.e. non-delayed) replication )	*
C*			OR						*
C*	     ( TYP(I) = "SUB" )						*
C*			OR						*
C*	     ( TYP(I) = "RPC" ) )  THEN					*
C*	  LINK(I) = 0							*
C*	ELSE IF ( TAG(I) is, sequentially, the last child Table B or	*
C*		  Table D mnemonic of the parent Table A or Table D	*
C*		  mnemonic indexed by JMPB(I) )  THEN			*
C*	  LINK(I) = 0							*
C*	ELSE								*
C*	  LINK(I) = the index of the jump/link table entry for the	*
C*		    Table B or Table D mnemonic which, sequentially,	*
C*		    is the next (i.e. following TAG(I)) child mnemonic	*
C*		    of the parent Table A or Table D mnemonic indexed	*
C*		    by JMPB(I)						*
C*	END IF								*
C*									*
C*    INTEGER IBT(I):							*
C*	IF ( ( TYP(I) = "CHR" ) OR					*
C*	     ( TYP(I) = "NUM" ) ) THEN					*
C*	  IBT(I) = bit width of Table B mnemonic TAG(I)			*
C*	ELSE IF ( ( TYP(I) = "DRB" ) OR					*
C*		  ( TYP(I) = "DRP" ) )  THEN				*
C*	  IBT(I) = bit width of delayed descriptor replication factor	*
C*		   (i.e. 1, 8, or 16, depending on the replication	*
C*		   scheme denoted by TAG(I))				*
C*	ELSE								*
C*	  IBT(I) = 0							*
C*	END IF								*
C*									*
C*    INTEGER IRF(I):							*
C*	IF ( TYP(I) = "NUM" ) THEN					*
C*	  IRF(I) = reference value of Table B mnemonic TAG(I)		*
C*	ELSE IF ( TYP(I) = "REP" ) THEN					*
C*	  IRF(I) = number of F=1 regular (i.e. non-delayed)		*
C*		   replications of Table D mnemonic TAG(JUMP(I))	*
C*	ELSE								*
C*	  IRF(I) = 0							*
C*	END IF								*
C*									*
C*    INTEGER ISC(I):							*
C*	IF ( TYP(I) = "NUM" ) THEN					*
C*	  ISC(I) = scale factor of Table B mnemonic TAG(I)		*
C*	ELSE IF ( TYP(I) = "SUB" ) THEN					*
C*	  ISC(I) = the index of the jump/link table entry which,	*
C*		   sequentially, constitutes the last element of the	*
C*		   jump/link tree for Table A mnemonic TAG(I)		*
C*	ELSE								*
C*	  ISC(I) = 0							*
C*	END IF								*
C************************************************************************

      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),
     .                JUMP(15000),LINK(15000),JMPB(15000),
     .                IBT(15000),IRF(15000),ISC(15000),
     .                ITP(15000),VALI(15000),KNTI(15000),
     .                ISEQ(15000,2),JSEQ(15000)
      COMMON /TABCCC/ ICDW,ICSC,ICRV

      CHARACTER*10 TAG
      CHARACTER*8  NEMO,NEMS,NEM
      CHARACTER*3  TYP
      CHARACTER*1  TAB
      DIMENSION    NEM(250,10),IRP(250,10),KRP(250,10)
      DIMENSION    DROP(10),JMP0(10),NODL(10),NTAG(10,2)
      LOGICAL      DROP

      DATA MAXLIM /10/

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK THE MNEMONIC
C  ------------------

C     Note that Table A mnemonics, in addition to being stored within
C     internal BUFR Table A array TABA(*,LUN), are also stored as
C     Table D mnemonics within internal BUFR Table D array TABD(*,LUN).
C     Thus, the following test is valid.

      CALL NEMTAB(LUN,NEMO,IDN,TAB,ITAB)
      IF(TAB.NE.'D') GOTO 900

C  STORE A SUBSET NODE AND JUMP/LINK THE TREE
C  ------------------------------------------

      CALL INCTAB(NEMO,'SUB',NODE)
      JUMP(NODE) = NODE+1
      JMPB(NODE) = 0
      LINK(NODE) = 0
      IBT (NODE) = 0
      IRF (NODE) = 0
      ISC (NODE) = 0

      CALL NEMTBD(LUN,ITAB,NSEQ,NEM(1,1),IRP(1,1),KRP(1,1))
      NTAG(1,1) = 1
      NTAG(1,2) = NSEQ
      JMP0(1)   = NODE
      LIMB      = 1

      ICDW = 0
      ICSC = 0
      ICRV = 0

C  THIS LOOP RESOLVES ENTITIES IN A SUBSET BY EMULATING RECURSION
C  --------------------------------------------------------------

1     DO N=NTAG(LIMB,1),NTAG(LIMB,2)

      NTAG(LIMB,1) = N+1
      NODL(LIMB)   = NTAB+1
      DROP(LIMB)   = N.EQ.NTAG(LIMB,2)

      CALL NEMTAB(LUN,NEM(N,LIMB),IDN,TAB,ITAB)
      NEMS = NEM(N,LIMB)

C  SPECIAL TREATMENT FOR CERTAIN OPERATOR DESCRIPTORS (TAB=C)
C  ----------------------------------------------------------

      IF(TAB.EQ.'C') THEN
         NODL(LIMB) = NTAB
         READ(NEMS,'(3X,I3)') IYYY
         IF(ITAB.EQ.1) THEN
            ICDW = IYYY-128
            IF(IYYY.EQ.0) ICDW = 0
         ELSEIF(ITAB.EQ.2) THEN
            ICSC = IYYY-128
            IF(IYYY.EQ.0) ICSC = 0
         ENDIF
      ELSE
         IREP = IRP(N,LIMB)
         IKNT = KRP(N,LIMB)
         JUM0 = JMP0(LIMB)
         CALL TABENT(LUN,NEMS,TAB,ITAB,IREP,IKNT,JUM0)
      ENDIF

      IF(TAB.EQ.'D') THEN

C        Note here how a new tree "LIMB" is created (and is then
C	 *immediately* recursively resolved!) whenever a Table D
C	 mnemonic contains another Table D mnemonic as one of its
C	 children.

         LIMB = LIMB+1
         IF(LIMB.GT.MAXLIM) GOTO 901
         CALL NEMTBD(LUN,ITAB,NSEQ,NEM(1,LIMB),IRP(1,LIMB),KRP(1,LIMB))
         NTAG(LIMB,1) = 1
         NTAG(LIMB,2) = NSEQ
         JMP0(LIMB)   = NTAB
         GOTO 1
      ELSEIF(DROP(LIMB)) THEN
2        LINK(NODL(LIMB)) = 0
         LIMB = LIMB-1
         IF(LIMB.EQ.0 ) THEN
            IF(ICDW.NE.0) GOTO 902
            IF(ICSC.NE.0) GOTO 903
            RETURN
         ENDIF
         IF(DROP(LIMB)) GOTO 2
         LINK(NODL(LIMB)) = NTAB+1
         GOTO 1
      ELSEIF(TAB.NE.'C') THEN
         LINK(NODL(LIMB)) = NTAB+1
      ENDIF

      ENDDO

      CALL BORT('TABSUB - SHOULD NOT GET HERE               ')
900   CALL BORT('TABSUB - SUBSET NODE NOT IN TABLE D: '//NEMO)
901   CALL BORT('TABSUB - TOO MANY LIMBS                    ')
902   CALL BORT('TABSUB - CHANGE DATA WIDTH OPERATOR NOT CANCELED')
903   CALL BORT('TABSUB - CHANGE DATA SCALE OPERATOR NOT CANCELED')
      END

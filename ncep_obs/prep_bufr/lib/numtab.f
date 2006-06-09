      SUBROUTINE NUMTAB(LUN,IDN,NEMO,TAB,IRET)

C************************************************************************
C* NUMTAB								*
C*									*
C* Given an integer IDN containing the bit-wise representation of	*
C* an FXY value, this subroutine searches for IDN within the internal	*
C* BUFR table B and D arrays and, if found, returns the corresponding	*
C* mnemonic and other information from within these arrays.  Otherwise,	*
C* it checks whether IDN is the bit-wise representation of an FXY value	*
C* for a table C operator descriptor, a replication descriptor, or a	*
C* replication factor and, if so, directly computes and returns similar	*
C* information about that descriptor.					*
C*									*
C* NUMTAB  ( LUN, IDN, NEMO, TAB, IRET )				*
C*									*
C* Input parameters:							*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					BUFR table arrays		*
C*	IDN		INTEGER		Bit-wise representation of	*
C*					FXY value			*
C*									*
C* Output parameters:							*
C*	NEMO		CHARACTER*(*)	Mnemonic corresponding to IDN	*
C*	TAB		CHARACTER	Type of FXY value that is	*
C*					bit-wise represented by IDN:	*
C*					 'B' = BUFR Table B descriptor	*
C*					 'C' = BUFR Table C descriptor	*
C*					 'D' = BUFR Table D descriptor	*
C*					 'R' = replication descriptor	*
C*					 'F' = replication factor	*
C*	IRET		INTEGER		Return value:			*
C*				      /               \			*
C*	      ------------------------                 ---------------	*
C* 		The interpretation of the return value IRET depends	*
C*		upon the return value of TAB, as follows:		*
C*									*
C*		IF ( TAB = 'B' ) THEN					*
C*		    IRET = positional index of IDN within internal	*
C*			   BUFR table B array				*
C*		ELSE IF ( TAB = 'C') THEN				*
C*		    IRET = the X portion of the FXY value that is	*
C*			   bit-wise represented by IDN			*
C*		ELSE IF ( TAB = 'D') THEN				*
C*		    IRET = positional index of IDN within internal	*
C*			   BUFR table D array				*
C*		ELSE IF ( TAB = 'R') THEN				*
C*		    IF ( ( the F portion of the FXY value that is	*
C*			   bit-wise represented by IDN ) = 1  (i.e.	*
C*			   regular (non-delayed) replication!) ) THEN	*
C*			IRET = ( (-1) * ( the Y portion of the FXY value*
C*			       that is bit-wise represented by IDN ) )	*
C*			   ( = ( (-1) * ( the number of F=1 regular	*
C*			       (i.e. non-delayed) replications! ) ) )	*
C*		    ELSE						*
C*			IRET = positional index of IDN within internal	*
C*			       replication array IDNR(*,*)		*
C*		    END IF						*
C*		ELSE IF ( TAB = 'F') THEN				*
C*		    IRET = positional index of IDN within internal	*
C*			   replication array IDNR(*,*)			*
C*		ELSE IF ( IRET = 0 ) THEN				*
C*		    IDN was not found in internal BUFR table B or D,	*
C*		    nor does it represent a table C operator descriptor,*
C*		    a replication descriptor, or a replication factor	*  
C*		END IF							*
C*	      --------------------------------------------------------	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************

C     Note that the values within the COMMON /REPTAB/ arrays were
C     initialized within subroutine BFRINI.

      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)

      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),
     .                TABA(50,32),TABB(250,32),TABD(250,32)

      CHARACTER*(*) NEMO
      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*6   ADN30,CID
      CHARACTER*3   TYPS
      CHARACTER*1   REPS,TAB

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      NEMO = ' '
      IRET = 0
      TAB = ' '

C  LOOK FOR A REPLICATOR OR A REPLICATOR FACTOR
C  --------------------------------------------

      IF(IDN.GE.IDNR(1,1) .AND. IDN.LE.IDNR(1,2)) THEN

C	 Note that the above test is checking whether IDN is the bit-wise
C	 representation of an FXY value denoting F=1 regular (i.e. non-delayed)
C	 replication, since, as was initialized within subroutine BFRINI,
C	 IDNR(1,1) = IFXY('101000'), and IDNR(1,2) = IFXY('101255').

         TAB  = 'R'
         IRET = -MOD(IDN,256)
         RETURN
      ENDIF

      DO I=2,5
      IF(IDN.EQ.IDNR(I,1)) THEN
         TAB  = 'R'
         IRET = I
         RETURN
      ELSEIF(IDN.EQ.IDNR(I,2)) THEN
         TAB  = 'F'
         IRET = I
         RETURN
      ENDIF
      ENDDO

C  LOOK FOR IDN IN TABLE D
C  -----------------------

      DO I=1,NTBD(LUN)
      IF(IDN.EQ.IDND(I,LUN)) THEN
         NEMO = TABD(I,LUN)(7:14)
         TAB  = 'D'
         IRET = I
         RETURN
      ENDIF
      ENDDO

C  LOOK FOR IDN IN TABLE B
C  -----------------------

      DO I=1,NTBB(LUN)
      IF(IDN.EQ.IDNB(I,LUN)) THEN
         NEMO = TABB(I,LUN)(7:14)
         TAB  = 'B'
         IRET = I
         RETURN
      ENDIF
      ENDDO

C  LOOK FOR IDN IN TABLE C
C  -----------------------

      CID = ADN30(IDN,6)
      CID = CID(1:3)
      IF(CID.EQ.'201' .OR. CID.EQ.'202' .OR. CID.EQ.'206') THEN
         NEMO = ADN30(IDN,6)
         TAB  = 'C'
         IRET = MOD(ID,10)
         RETURN
      ENDIF

      RETURN
      END

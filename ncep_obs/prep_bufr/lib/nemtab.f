      SUBROUTINE NEMTAB(LUN,NEMO,IDN,TAB,IRET)

C************************************************************************
C* NEMTAB								*
C*									*
C* This subroutine searches for mnemonic NEMO within the internal BUFR	*
C* table B and D arrays and, if found, returns information about that	*
C* mnemonic from within these arrays.  Otherwise, it checks whether	*
C* NEMO is a table C operator descriptor and, if so, directly computes	*
C* and returns similar information about that descriptor.		*
C*									*
C* NEMTAB  ( LUN, NEMO, IDN, TAB, IRET )				*
C*									*
C* Input parameters:							*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					BUFR table arrays		*
C*	NEMO		CHARACTER*(*)	Mnemonic to search for		*
C*									*
C* Output parameters:							*
C*	IDN		INTEGER		Bit-wise representation of FXY	*
C*					value corresponding to NEMO	*
C*					(if NEMO was found!)		*
C*	TAB		CHARACTER	Internal table array ('B', 'C',	*
C*					or 'D') in which NEMO was found	*
C*	IRET		INTEGER		Positional index of NEMO	*
C*					within TAB			*
C*					  0 = NEMO was not found within	*
C*					      any of the table arrays	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************

      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),
     .                TABA(50,32),TABB(250,32),TABD(250,32)

      CHARACTER*(*) NEMO
      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*8   NEMT
      CHARACTER*1   TAB
      LOGICAL       FOLVAL

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      FOLVAL = NEMO(1:1).EQ.'.'
      IRET = 0
      TAB = ' '

C  LOOK FOR NEMO IN TABLE B
C  ------------------------

      DO 1 I=1,NTBB(LUN)
      NEMT = TABB(I,LUN)(7:14)
      IF(NEMT.EQ.NEMO) THEN
         IDN  = IDNB(I,LUN)
         TAB  = 'B'
         IRET = I
         RETURN
      ELSEIF(FOLVAL.AND.NEMT(1:1).EQ.'.') THEN
         DO J=2,LEN(NEMT)
         IF(NEMT(J:J).NE.'.' .AND. NEMT(J:J).NE.NEMO(J:J)) GOTO 1
         ENDDO
         IDN  = IDNB(I,LUN)
         TAB  = 'B'
         IRET = I
         RETURN
      ENDIF
1     ENDDO

C  DON'T LOOK IN TABLE D FOR FOLLOWING VALUE-MNEMONICS
C  ---------------------------------------------------

      IF(FOLVAL) RETURN

C  LOOK IN TABLE D IF WE GOT THIS FAR
C  ----------------------------------

      DO I=1,NTBD(LUN)
      NEMT = TABD(I,LUN)(7:14)
      IF(NEMT.EQ.NEMO) THEN
         IDN  = IDND(I,LUN)
         TAB  = 'D'
         IRET = I
         RETURN
      ENDIF
      ENDDO

C  HERE CHECK FOR TABLE C OPERATOR DESCRIPTORS
C  -------------------------------------------

      IF(NEMO(1:3).EQ.'201' .OR.
     .   NEMO(1:3).EQ.'202' .OR.
     .   NEMO(1:3).EQ.'203' .OR.
     .   NEMO(1:3).EQ.'206' ) THEN
         READ(NEMO,'(I6)') IRET
         IDN = IFXY(NEMO)
         TAB = 'C'
         IRET = MOD(IRET/1000,10)
         RETURN
      ENDIF

      RETURN
      END

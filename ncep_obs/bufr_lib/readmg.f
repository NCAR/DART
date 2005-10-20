      SUBROUTINE READMG(LUNIT,SUBSET,JDATE,IRET)

C************************************************************************
C* READMG								*
C*									*
C* This subroutine reads (into memory) the next BUFR message from	*
C* logical unit number LUNIT, which should itself have already been	*
C* opened for input operations.						*
C*									*
C* READMG  ( LUNIT, SUBSET, JDATE, IRET )				*
C*									*
C* Input parameters:							*
C*	LUNIT		INTEGER		FORTRAN logical unit number	*
C*									*
C* Output parameters:							*
C*	SUBSET		CHARACTER*8	Table A mnemonic for		*
C*					BUFR message			*
C*	JDATE		INTEGER		Date-time from Section 1 of	*
C*					BUFR message, in format of	*
C*					either YYMMDDHH or YYYYMMDDHH,	*
C*					depending on DATELEN() value	*
C*	IRET		INTEGER		Return code:			*
C*					  0 = normal return		*
C*					 -1 = there are no more BUFR	*
C*					      messages in LUNIT		*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
      COMMON /DATELN/ LENDAT
 
      CHARACTER*8 SEC0,SUBSET
      CHARACTER*4 BUFR
      DIMENSION   IEC0(2)
      EQUIVALENCE (SEC0,IEC0)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      IRET = 0
 
C  CHECK THE FILE STATUS
C  ---------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      CALL WTSTAT(LUNIT,LUN,IL,1)
      IMSG = 8/NBYTW+1
 
C  READ A MESSAGE INTO A MESSAGE BUFFER - SKIP DX MESSAGES
C  -------------------------------------------------------
 
1     SEC0 = ' '
      READ(LUNIT,ERR=902,END=100) SEC0,(MBAY(I,LUN),I=IMSG,LMSG(SEC0))

C     Confirm that the first 4 bytes of SEC0 contain 'BUFR' encoded in
C     CCITT IA5 (i.e. ASCII).

      CALL CHRTRNA(BUFR,SEC0,4)
      IF(BUFR.NE.'BUFR') GOTO 100

C     Copy SEC0 into the front of MBAY so that MBAY now contains the
C     entire BUFR message.

      DO I=1,IMSG-1
      MBAY(I,LUN) = IEC0(I)
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

C  ENTRY DATELEN SETS THE LENGTH OF THE DATE INTEGER RETURN FROM READS
C  -------------------------------------------------------------------

      ENTRY DATELEN(LEN)
C************************************************************************
C* DATELEN								*
C*									*
C* This entry point is used to specify the length of date-time values	*
C* that will be output by future calls to any of the subroutines which	*
C* read BUFR messages (e.g. READMG, READERM, READFT, READTJ, etc.).	*
C* Possible values are 8 (which is the default) and 10.			*
C*									*
C* DATELEN  ( LEN )							*
C*									*
C* Input parameters:							*
C*	LEN		INTEGER		Length of date-time values to	*
C*					be output by read subroutines:	*
C*					  8 =   YYMMDDHH (2-digit year)	*
C*					 10 = YYYYMMDDHH (4-digit year)	*
C************************************************************************
      IF(LEN.NE.8 .AND. LEN.NE.10) CALL BORT('DATELEN - BAD LEN')
      LENDAT = LEN
      RETURN
 
C  ERROR EXITS
C  -----------
 
900   CALL BORT('READMG - FILE IS CLOSED                       ')
901   CALL BORT('READMG - FILE IS OPEN FOR OUTPUT              ')
902   CALL BORT('READMG - I/O ERROR READING MESSAGE            ')
      END

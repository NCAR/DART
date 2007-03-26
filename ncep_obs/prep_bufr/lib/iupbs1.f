      FUNCTION IUPBS1  ( MBAY, NBYT )

C************************************************************************
C* IUPBS1								*
C*									*
C* Given a BUFR message	contained within array MBAY, this function	*
C* unpacks and returns the binary integer contained within byte NBYT	*
C* of Section 1 of the BUFR message.  Note that the start of the BUFR	*
C* message (i.e. "BUFR") must be aligned on the first 4 bytes of MBAY.	*
C*									*
C* IUPBS1  ( MBAY, NBYT )						*
C*									*
C* Input parameters:							*
C*	MBAY(*)		INTEGER		Array containing BUFR message	*
C*	NBYT		INTEGER		Byte within Section 1 of BUFR	*
C*					message to be unpacked		*
C*									*
C* Output parameters:							*
C*	IUPBS1		INTEGER		Unpacked integer		*
C**									*
C* Log:									*
C* J. Ator/NCEP		05/01						*
C************************************************************************
 
      DIMENSION MBAY(*)

      DATA IFIRST / 0 /
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C     If this is the first call to this subroutine, then call subroutine
C     WRDLEN to initialize some important information about the local
C     machine, just in case subroutine OPENBF hasn't been called yet!
                                                                        
      IF ( IFIRST .EQ. 0 ) THEN
         CALL WRDLEN                                                    
         IFIRST = 1                                                     
      ENDIF                                                             
 

C     Note that there are 8 bytes within Section 0 that must be skipped.

      IUPBS1 = IUPB ( MBAY, NBYT+8, 8 )


      RETURN
      END

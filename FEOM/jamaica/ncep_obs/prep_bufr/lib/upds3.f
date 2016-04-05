      SUBROUTINE UPDS3  ( MBAY, CDS3, NDS3 )

C************************************************************************
C* UPDS3								*
C*									*
C* Given a BUFR message	contained within array MBAY, this subroutine	*
C* unpacks and returns the descriptors contained within Section 3 of	*
C* the BUFR message.  Note that the start of the BUFR message		*
C* (i.e. "BUFR") must be aligned on the first 4 bytes of MBAY.		*
C* Note also that this subroutine does not recursively resolve sequence	*
C* descriptors that appear within Section 3; rather, what is returned	*
C* is the exact list of descriptors as it appears within Section 3.	*
C*									*
C* UPDS3  ( MBAY, CDS3, NDS3 )						*
C*									*
C* Input parameters:							*
C*	MBAY(*)		INTEGER		Array containing BUFR message	*
C*									*
C* Output parameters:							*
C*	CDS3(*)		CHARACTER*6	Unpacked list of descriptors	*
C*	NDS3		INTEGER		Number of descriptors returned	*
C*					within CDS3(*)			*
C**									*
C* Log:									*
C* J. Ator/NCEP		05/01						*
C************************************************************************

      DIMENSION MBAY(*)
                                                                        
      CHARACTER*6 CDS3(*), ADN30
                                                                        
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
 

C     Skip Section 0.

      IPT = 8


C     Skip Section 1, but check for existence of Section 2.

      ISC2 = IUPB ( MBAY, IPT + 8,  1 )

      IPT = IPT + IUPB ( MBAY, IPT + 1, 24 )
 
 
C     Skip Section 2 if it exists.
 
      IPT = IPT + ( IUPB ( MBAY, IPT + 1, 24 ) * ISC2 )
 

C     Get the length of Section 3.
 
      LEN3 = IUPB ( MBAY, IPT + 1 , 24 )


C     Unpack the Section 3 descriptors.

      NDS3 = 0
      DO JJ = 8, ( LEN3 - 1 ), 2
         NDS3 = NDS3 + 1
         CDS3 (NDS3) = ADN30 ( IUPB ( MBAY, IPT + JJ, 16 ), 6 )
      ENDDO
 

      RETURN
      END

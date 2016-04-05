C-----------------------------------------------------------------------
C  UNPACK UP AN INTEGER FROM A PACKED CHARACTER ARRAY (FUNCTION)
C----------------------------------------------------------------------
 
      FUNCTION IUPM(CBAY,NBITS)
 
C************************************************************************
C* IUPM									*
C*									*
C* This function unpacks and returns a binary integer contained within	*	
C* NBITS bits of CBAY, starting with the first bit of the first byte	*
C* of CBAY.								*
C*									*
C* IUPM  ( CBAY, NBITS )						*
C*									*
C* Input parameters:							*
C*	CBAY		CHARACTER*(*)	Character string containing	*
C*					packed integer			*
C*	NBITS		INTEGER		Number of bits within CBAY to	*
C*					be unpacked			*
C*									*
C* Output parameters:							*
C*	IUPM		INTEGER		Unpacked integer		*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*8 CBAY,CINT
      EQUIVALENCE(CINT,INT)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      IF(NBITS.GT.NBITW) CALL BORT('IUPM - NBITS>WRD LEN')
      CINT = CBAY
      INT  = IREV(INT)
      IUPM = ISHFT(INT,NBITS-NBITW)
      RETURN
      END

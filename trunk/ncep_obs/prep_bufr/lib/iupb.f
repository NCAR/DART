C----------------------------------------------------------------------
C  UNPACK UP AN INTEGER FROM A PACKED INTEGER ARRAY (FUNCTION)
C----------------------------------------------------------------------

      FUNCTION IUPB(MBAY,NBYT,NBIT)

C************************************************************************
C* IUPB									*
C*									*
C* This function unpacks and returns a binary integer contained within	*	
C* NBIT bits of MBAY, starting with the first bit of byte NBYT.		*
C*									*
C* IUPB  ( MBAY, NBYT, NBIT )						*
C*									*
C* Input parameters:							*
C*	MBAY		INTEGER(*)	Packed binary array		*
C*	NBYT		INTEGER		Byte within MBAY at whose first	*
C*					bit to begin unpacking		*
C*	NBIT		INTEGER		Number of bits within MBAY to	*
C*					be unpacked			*
C*									*
C* Output parameters:							*
C*	IUPB		INTEGER		Unpacked integer		*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
 
      DIMENSION MBAY(*)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      MBIT = (NBYT-1)*8
      CALL UPB(IRET,NBIT,MBAY,MBIT)
      IUPB = IRET
      RETURN
      END

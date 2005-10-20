      SUBROUTINE UPC(CHR,NCHR,IBAY,IBIT)

C************************************************************************
C* UPC									*
C*									*
C* This subroutine unpacks and returns a character string contained	*
C* within NCHR bytes of IBAY, starting with bit (IBIT+1).  On output,	*
C* IBIT is updated to point to the last bit that was unpacked.		*
C* Note that the string to be unpacked does not necessarily need to	*
C* be aligned on a byte boundary within IBAY.				* 
C*									*
C* UPC  ( CHR, NCHR, IBAY, IBIT )					*
C*									*
C* Input parameters:							*
C*	NCHR		INTEGER		Number of bytes within IBAY to	*
C*					be unpacked			*
C*	IBAY		INTEGER(*)	Packed binary array		*
C*									*
C* Input and output parameters:						*
C*	IBIT		INTEGER		Pointer within IBAY		*
C*									*
C* Output parameters:							*
C*	CHR		CHARACTER*(*)	Unpacked character string of	*
C*					length NCHR			*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation				*
C************************************************************************
 
      COMMON /CHARAC/ IASCII,IATOE(0:255),IETOA(0:255)
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*(*) CHR
      CHARACTER*8   CVAL
      DIMENSION     IBAY(*)
      EQUIVALENCE   (CVAL,IVAL)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      LB = IORD(NBYTW)
      DO I=1,NCHR
      CALL UPB(IVAL,8,IBAY,IBIT)
      CHR(I:I) = CVAL(LB:LB)
      IF(IASCII.EQ.0) CALL IPKM(CHR(I:I),1,IATOE(IUPM(CHR(I:I),8)))
      ENDDO
 
      RETURN
      END

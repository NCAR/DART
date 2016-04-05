C----------------------------------------------------------------------
C  REVERSE A WEIRD INTEGER
C----------------------------------------------------------------------

      FUNCTION IREV(N)

C************************************************************************
C* IREV									*
C*									*
C* By definition (within WMO Manual 306), a BUFR message is a stream of	*
C* individual octets (i.e. bytes) that is independent of any particular	*
C* machine representation.  However, the BUFRLIB software often needs	*
C* to interpret two or more adjacent bytes as an integer; therefore,	*
C* when doing so, it is critical to know whether the local machine uses	*
C* the "big-endian" (i.e. left->right) or "little-endian" (right->left)	*
C* scheme for numbering the bytes within a machine word!  By default,	*
C* BUFRLIB decodes multi-byte integers according to the "big-endian"	*
C* numbering scheme; thus, if the local machine is "little-endian"	*
C* instead, then this function IREV will return a copy of the input	*
C* integer field with the bytes reversed so that it can be properly	*
C* read or written (depending on whether input or output operations,	*
C* respectively, are being performed!).  If, on the other hand, the	*
C* local machine is already "big-endian", then no such reversal is	*
C* necessary and the function simply returns a copy of the same integer	*
C* that was input.							*
C*									*
C* IREV  ( N )								*
C*									*
C* Input parameters:							*
C*	N		INTEGER		Integer (to possibly reverse)	*
C*									*
C* Output parameters:							*
C*	IREV		INTEGER		(Possibly reversed) copy of N	*
C*									*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*8 CINT,DINT
      EQUIVALENCE(CINT,INT)
      EQUIVALENCE(DINT,JNT)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
C     Note that the value of NREV is set within subroutine WRDLEN.

      IF(NREV.EQ.0) THEN
         IREV = N
      ELSE
         INT = N
         DO I=1,NBYTW
         DINT(I:I) = CINT(IORD(I):IORD(I))
         ENDDO
         IREV = JNT
      ENDIF
 
      RETURN
      END

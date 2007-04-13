C-----------------------------------------------------------------------
C  PACK AN INTEGER INTO A CHARACTER ARRAY
C----------------------------------------------------------------------

      SUBROUTINE IPKM(CBAY,NBYT,N)

C************************************************************************
C* IPKM									*
C*									*
C* This subroutine packs an integer N into a character string of length	*
C* NBYT bytes.								*
C*									*
C* IPKM  ( CBAY, NBYT, N )						*
C*									*
C* Input parameters:							*
C*	NBYT		INTEGER		Number of bytes into which to	*
C*					pack N				*
C*	N		INTEGER		Integer to be packed		*
C*									*
C* Output parameters:							*
C*	CBAY		CHARACTER*(*)	String of length NBYT containing*
C*					packed integer N		*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation				*
C************************************************************************
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*8 CBAY,CINT
      EQUIVALENCE(CINT,INT)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------

      IF(NBYT.GT.NBYTW) CALL BORT('IPKM - NBYT>WRD LEN')
 
C     Note that the widths of input variable N and local variable INT
C     will both be equal to the default size of an integer (= NBYTW),
C     since they aren't specifically declared otherwise!
 
      INT = IREV(ISHFT(N,(NBYTW-NBYT)*8))
      DO I=1,NBYT
      CBAY(I:I) = CINT(I:I)
      ENDDO
      RETURN
      END

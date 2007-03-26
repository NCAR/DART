C----------------------------------------------------------------------
C  COPY CHARACTERS INTO A BIT ARRAY
C----------------------------------------------------------------------

      SUBROUTINE PKC(CHR,NCHR,IBAY,IBIT)

C************************************************************************
C* PKC									*
C*									*
C* This subroutine packs NCHR characters from CHR within NCHR bytes of	*
C* IBAY, starting with bit (IBIT+1).  On output, IBIT is updated to	*
C* point to the last bit that was packed.  Note that there is no	*
C* guarantee that the NCHR characters will be aligned on byte		*
C* boundaries when packed within IBAY.					*
C*									*
C* PKC  ( CHR, NCHR, IBAY, IBIT )					*
C*									*
C* Input parameters:							*
C*	CHR		CHARACTER*(*)	Character string		*
C*	NCHR		INTEGER		Number of characters of CHR to	*
C*				        be packed within IBAY		*
C*									*
C* Input and output parameters:						*
C*	IBIT		INTEGER		Pointer within IBAY		*
C*									*
C* Output parameters:							*
C*	IBAY		INTEGER(*)	Packed binary array		*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
 
      COMMON /CHARAC/ IASCII,IATOE(0:255),IETOA(0:255)
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*(*) CHR
      CHARACTER*1   CVAL(8)
      DIMENSION     IBAY(*)
      EQUIVALENCE   (CVAL,IVAL)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      IF(NCHR.GT.LEN(CHR)) CALL BORT('PKC - CHR < NCHR')
      LB = IORD(NBYTW)

C     LB now points to the "low-order" (i.e. least significant) byte
C     within a machine word.

      IVAL = 0
      NBIT = 8
 
      DO I=1,NCHR
      CVAL(LB) = CHR(I:I)

C     If the machine is EBCDIC, then translate character CVAL(LB) from
C     EBCDIC to ASCII.

      IF(IASCII.EQ.0) CALL IPKM(CVAL(LB),1,IETOA(IUPM(CVAL(LB),8)))

      NWD  = IBIT/NBITW + 1
      NBT  = MOD(IBIT,NBITW)
      INT = ISHFT(IVAL,NBITW-NBIT)
      INT = ISHFT(INT,-NBT)
      MSK = ISHFT(  -1,NBITW-NBIT)
      MSK = ISHFT(MSK,-NBT)
      IBAY(NWD) = IREV(IOR(IAND(IREV(IBAY(NWD)),NOT(MSK)),INT))
      IF(NBT+NBIT.GT.NBITW) THEN

C        This character will not fit within the current word (i.e.
C        array member) of IBAY, because there are less than 8 bits of
C        space left.  Store as many bits as will fit within the current
C        word and then store the remaining bits within the next word.

         INT = ISHFT(IVAL,2*NBITW-(NBT+NBIT))
         MSK = ISHFT(  -1,2*NBITW-(NBT+NBIT))
         IBAY(NWD+1) = IREV(IOR(IAND(IREV(IBAY(NWD+1)),NOT(MSK)),INT))
      ENDIF
      IBIT = IBIT + NBIT
      ENDDO
 
      RETURN
      END

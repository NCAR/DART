C----------------------------------------------------------------------
C  PACK UP A NUMBER ACCORDING TO SPECS
C----------------------------------------------------------------------

      SUBROUTINE PKB(NVAL,NBITS,IBAY,IBIT)

C************************************************************************
C* PKB									*
C*									*
C* This subroutine packs NVAL within NBITS bits of IBAY, starting with	*
C* bit (IBIT+1).  On output, IBIT is updated to point to the last bit	*
C* that was packed.							*
C*									*
C* PKB  ( NVAL, NBITS, IBAY, IBIT )					*
C*									*
C* Input parameters:							*
C*	NVAL		INTEGER		Integer	to be packed		*
C*	NBITS		INTEGER		Number of bits of IBAY within	*
C*					which to pack NVAL		*
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
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      DIMENSION IBAY(*)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      NWD  = IBIT/NBITW + 1
      NBT  = MOD(IBIT,NBITW)
      IVAL = NVAL
      IF(ISHFT(IVAL,-NBITS).GT.0) IVAL = -1
      INT = ISHFT(IVAL,NBITW-NBITS)
      INT = ISHFT(INT,-NBT)
      MSK = ISHFT(  -1,NBITW-NBITS)
      MSK = ISHFT(MSK,-NBT)
      IBAY(NWD) = IREV(IOR(IAND(IREV(IBAY(NWD)),NOT(MSK)),INT))
      IF(NBT+NBITS.GT.NBITW) THEN

C        There are less than NBITS bits remaining within the current word
C        (i.e. array member) of IBAY, so store as many bits as will fit
C        within the current word and then store the remaining bits within
C        the next word.

         INT = ISHFT(IVAL,2*NBITW-(NBT+NBITS))
         MSK = ISHFT(  -1,2*NBITW-(NBT+NBITS))
         IBAY(NWD+1) = IREV(IOR(IAND(IREV(IBAY(NWD+1)),NOT(MSK)),INT))
      ENDIF
 
      IBIT = IBIT + NBITS
 
      RETURN
      END

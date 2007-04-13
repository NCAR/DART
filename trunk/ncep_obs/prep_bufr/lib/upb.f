C----------------------------------------------------------------------
C  UNPACK UP AN INTEGER
C----------------------------------------------------------------------
      SUBROUTINE UPB(NVAL,NBITS,IBAY,IBIT)
C************************************************************************
C* UPB									*
C*									*
C* This subroutine unpacks and returns a binary integer contained	*
C* within NBITS bits of IBAY, starting with bit (IBIT+1).  On output,	*
C* IBIT is updated to point to the last bit that was unpacked.		*
C*									*
C* UPB  ( NVAL, NBITS, IBAY, IBIT )					*
C*									*
C* Input parameters:							*
C*	NBITS		INTEGER		Number of bits within IBAY to	*
C*					be unpacked			*
C*	IBAY		INTEGER(*)	Packed binary array		*
C*									*
C* Input and output parameters:						*
C*	IBIT		INTEGER		Pointer within IBAY		*
C*									*
C* Output parameters:							*
C*	NVAL		INTEGER		Unpacked integer		*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		11/00	Added docblock and check for NBITS=0	*
C************************************************************************
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      DIMENSION IBAY(*)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
C  IF NBITS=0, THEN JUST SET NVAL=0 AND RETURN
C  -------------------------------------------
      IF(NBITS.EQ.0)THEN
        NVAL=0
        RETURN
      ENDIF 

      NWD = (IBIT)/NBITW+1
      NBT = MOD(IBIT,NBITW)
      INT = ISHFT(IREV(IBAY(NWD)),NBT)
      INT = ISHFT(INT,NBITS-NBITW)
      LBT = NBT+NBITS
      IF(LBT.GT.NBITW) JNT = IREV(IBAY(NWD+1))
      IF(LBT.GT.NBITW) INT = IOR(INT,ISHFT(JNT,LBT-2*NBITW))
      IBIT = IBIT+NBITS
      NVAL = INT
      RETURN
      END

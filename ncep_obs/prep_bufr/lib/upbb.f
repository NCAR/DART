C----------------------------------------------------------------------
C  UNPACK UP AN INTEGER FOR SUBROUTINE RDTREE
C----------------------------------------------------------------------
      SUBROUTINE UPBB(INT,NBIT,MBIT,MBAY)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)

      DIMENSION MBAY(*)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      NWD = MBIT/NBITW + 1
      NBT = MOD(MBIT,NBITW)
      LBT = NBT+NBIT
      IBA = IREV(MBAY(NWD))
      INT = ISHFT(ISHFT(IBA,NBT),NBIT-NBITW)
      IF(LBT.GT.NBITW) THEN
         IBA = IREV(MBAY(NWD+1))
         INT = IOR(INT,ISHFT(IBA,LBT-2*NBITW))
      ENDIF
      RETURN
      END

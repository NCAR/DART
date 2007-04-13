C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RDCMPS(LUN)
C************************************************************************
C* RDCMPS								*
C*									*
C* This subroutine uncompresses the next subset from within the current	*
C* compressed BUFR message that is in memory, and it stores the output	*
C* values within the internal VAL(*,LUN) array.				*
C*									*
C* RDCMPS  ( LUN )							*
C*									*
C* Input parameters:							*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					arrays for current BUFR message	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		11/00	Added documentation			*
C************************************************************************
 
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),
     .                JUMP(15000),LINK(15000),JMPB(15000),
     .                IBT(15000),IRF(15000),ISC(15000),
     .                ITP(15000),VALI(15000),KNTI(15000),
     .                ISEQ(15000,2),JSEQ(15000)
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)
 
      CHARACTER*10 TAG
      CHARACTER*8  CREF,CVAL
      CHARACTER*3  TYP
      EQUIVALENCE  (CVAL,RVAL)
      REAL*8       VAL,RVAL
 
C-----------------------------------------------------------------------
C     Statement function to compute BUFR "missing value" for field
C     of length LBIT bits:

      LPS(LBIT) = MAX(2**(LBIT)-1,1)

C     Statement function to decode the encoded BUFR value IVAL according
C     to the scale and reference values that are stored within index NODE
C     of the internal arrays ISC(*) and IRF(*):

      UPS(NODE) = FLOAT(IVAL+IRF(NODE))*10.**(-ISC(NODE))
C-----------------------------------------------------------------------
 
C  SETUP THE SUBSET TEMPLATE
C  -------------------------
 
      CALL USRTPL(LUN,1,1)
 
C  UNCOMPRESS A SUBSET INTO THE VAL ARRAY ACCORDING TO TABLE B
C  -----------------------------------------------------------
 
      NELM = NVAL(LUN)
      NSBS = NSUB(LUN)

C     Note that we are going to unpack the (NSBS)th subset from within
C     the current BUFR message.

      IBIT = MBYT(LUN)
 
C     Loop through each element of the subset.

      DO N=1,NELM
      NODE = INV(N,LUN)
      NBIT = IBT(NODE)
      ITYP = ITP(NODE)

C     In each of the following code blocks, the "local reference value"
C     for the element is determined first, followed by the 6-bit value
C     which indicates how many bits are used to store the increment
C     (i.e. offset) from this "local reference value".  Then, we jump
C     ahead to where this increment is stored for this particular subset,
C     unpack it, and add it to the "local reference value" to determine
C     the final uncompressed value for this element from this subset.

C     Note that, if an element has the same final uncompressed value
C     for each subset in the message, then the encoding rules for BUFR
C     compression dictate that the "local reference value" will be equal
C     to this value, the 6-bit increment length indicator will have
C     a value of zero, and the actual increments themselves will be
C     omitted from the message.

      IF(ITYP.EQ.2) THEN

C        This is a numeric element.

         CALL UPB(LREF,NBIT,MBAY(1,LUN),IBIT)
         CALL UPB(LINC,   6,MBAY(1,LUN),IBIT)
         JBIT = IBIT + LINC*(NSBS-1)
         CALL UPB(NINC,LINC,MBAY(1,LUN),JBIT)
         IF(NINC.EQ.LPS(LINC)) NINC = LPS(NBIT)
         IVAL = LREF+NINC
         IF(IVAL.LT.LPS(NBIT)) VAL(N,LUN) = UPS(NODE)
         IBIT = IBIT + LINC*MSUB(LUN)
      ELSEIF(ITYP.EQ.3) THEN

C        This is a character element.

         CALL UPC(CREF,NBIT/8,MBAY(1,LUN),IBIT)
         CALL UPB(LINC,   6,MBAY(1,LUN),IBIT)
         JBIT = IBIT + LINC*(NSBS-1)*8
         CALL UPC(CVAL,LINC,MBAY(1,LUN),JBIT)
         VAL(N,LUN) = RVAL
         IBIT = IBIT + 8*LINC*MSUB(LUN)
      ENDIF
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CKTABA(LUN,SUBSET,JDATE,IRET)
 
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
      COMMON /PADESC/ IBCT,IPD1,IPD2,IPD3,IPD4                          
      COMMON /UNPTYP/ MSGUNP(32)
      COMMON /DATELN/ LENDAT
 
      CHARACTER*8 SUBSET
      CHARACTER*1 TAB
      LOGICAL     TRYBT
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      TRYBT = .TRUE.
 
C  PARSE SECTION 1
C  ---------------
 
      IAD1 = 8
      LEN1 = IUPB(MBAY(1,LUN),IAD1+ 1,24)
      LEN2 = IUPB(MBAY(1,LUN),IAD1+ 8, 1)
      MTYP = IUPB(MBAY(1,LUN),IAD1+ 9, 8)
      MSBT = IUPB(MBAY(1,LUN),IAD1+10, 8)
      MEAR = MOD(IUPB(MBAY(1,LUN),IAD1+13, 8),100)
      MMON = IUPB(MBAY(1,LUN),IAD1+14, 8)
      MDAY = IUPB(MBAY(1,LUN),IAD1+15, 8)
      MOUR = IUPB(MBAY(1,LUN),IAD1+16, 8)
      MMIN = IUPB(MBAY(1,LUN),IAD1+17, 8)
      MCEN = MAX(0,IUPB(MBAY(1,LUN),IAD1+18, 8)-MIN(MEAR,1))
 
      IF(LENDAT.EQ.10) THEN
         JDATE = MCEN*10**8+MEAR*10**6+MMON*10**4+MDAY*10**2+MOUR
         JDATE = I4DY(JDATE)
      ELSE
         JDATE = MEAR*10**6+MMON*10**4+MDAY*10**2+MOUR
      ENDIF
 
C  DON'T PARSE BUFR TABLE MESSAGES
C  ------------------------------
 
      IF(MTYP.EQ.11) THEN
         IRET = 11
         RETURN
      ENDIF
 
C  PARSE SECTION 2
C  ---------------
 
      IAD2 = IAD1+LEN1
      LEN2 = IUPB(MBAY(1,LUN),IAD2+1,24) * LEN2
 
C  PARSE SECTION 3
C  ---------------
 
      IAD3 = IAD2+LEN2
      LEN3 = IUPB(MBAY(1,LUN),IAD3+1 ,24)
      JSUB = IUPB(MBAY(1,LUN),IAD3+5 ,16)
      NCMP = IUPB(MBAY(1,LUN),IAD3+7 ,8 )
      KSUB = IUPB(MBAY(1,LUN),IAD3+8 ,16)
      ISUB = IUPB(MBAY(1,LUN),IAD3+10,16)
 
C  LOCATE SECTION 4
C  ----------------
 
      IAD4 = IAD3+LEN3
      LEN4 = IUPB(MBAY(1,LUN),IAD4+1,24)
 
C  IF ISUB FROM SECTION 3 DEFINES TABLE A THEN MSGUNP=0
C  ----------------------------------------------------
 
5     CALL NUMTAB(LUN,ISUB,SUBSET,TAB,ITAB)
      CALL NEMTBAX(LUN,SUBSET,MTY1,MSB1,INOD)
      IF(INOD.GT.0) THEN
         MBYT(LUN) = (IAD4+4)
         MSGUNP(LUN) = 0
         GOTO 10
      ENDIF
 
C  IF KSUB FROM SECTION 3 DEFINES TABLE A THEN MSGUNP=1
C  ----------------------------------------------------
 
      CALL NUMTAB(LUN,KSUB,SUBSET,TAB,ITAB)
      CALL NEMTBAX(LUN,SUBSET,MTY1,MSB1,INOD)
      IF(INOD.GT.0) THEN
         MBYT(LUN) = 8*(IAD4+4)
         MSGUNP(LUN) = 1
         GOTO 10
      ENDIF

C  IF MTYP/MSBT DEFINES TABLE A ALSO CHECK FOR A BYTE COUNT DESCRIPTOR
C  -------------------------------------------------------------------
 
      WRITE(SUBSET,'("NC",2I3.3)') MTYP,MSBT
      CALL NEMTBAX(LUN,SUBSET,MTY1,MSB1,INOD)
      IF(INOD.GT.0 .AND. KSUB.EQ.IBCT) THEN
         MBYT(LUN) = (IAD4+4)
         MSGUNP(LUN) = 0
         GOTO 10
      ELSEIF(INOD.GT.0) THEN
         MBYT(LUN) = 8*(IAD4+4)
         MSGUNP(LUN) = 1
         GOTO 10
      ENDIF

C  LAST DESPARATE ATTEMPT - SEE IF A BUFR TABLE IS DEFINED IN OPENBT
C  -----------------------------------------------------------------

      IF(TRYBT) THEN
         TRYBT = .FALSE.
         CALL OPENBT(LUNDX,MTYP)
         IF(LUNDX.GT.0) THEN
            CALL RDUSDX(LUNDX,LUN)                                         
            GOTO 5
         ENDIF
      ENDIF
 
C  IF ALL ATTEMPTS TO DEFINE TABLE A FAIL SKIP GIVE UP
C  ---------------------------------------------------
 
      PRINT*,'UNRECOGNISED TABLE A DESCRIPTOR:',SUBSET
      IRET = -1
      RETURN
 
C  CHECK THE VALIDITY OF THE MTYP/MSBT AND FOR COMPRESSION (MSGUNP=2)
C  ------------------------------------------------------------------
 
10    IF(MTYP.NE.MTY1.OR.MSBT.NE.MSB1) GOTO 900
      IF(IAND(NCMP,64).GT.0) MSGUNP(LUN) = 2
 
C  SET THE OTHER REQUIRED PARAMETERS AND RETURN SUCCESSFULLY
C  ---------------------------------------------------------
 
      IDATE(LUN) = I4DY(JDATE)
      NMSG (LUN) = NMSG(LUN)+1
      INODE(LUN) = INOD
      MSUB (LUN) = JSUB
      NSUB (LUN) = 0
      IRET = 0
      RETURN
 
C  HARD ERROR EXIT
C  ---------------
 
900   CALL BORT('CKTABA - MESG TYP/SUBTYP MISMATCH:'//SUBSET)
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ADDATE(IDATE,JH,JDATE)                                 
                                                                        
      DIMENSION   MON(12)                                               
                                                                        
      DATA MON/31,28,31,30,31,30,31,31,30,31,30,31/                     
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IY = IDATE/1000000
      IM = MOD(IDATE/10000  ,100)                                       
      ID = MOD(IDATE/100    ,100)                                       
      IH = MOD(IDATE        ,100)                                       
      IH = IH+JH                                                        
                                                                        
      IF(MOD(IY,4)  .NE.0) MON(2) = 28                                    
      IF(MOD(IY,4)  .EQ.0) MON(2) = 29                                    
      IF(MOD(IY,100).EQ.0) MON(2) = 28                                    
      IF(MOD(IY,400).EQ.0) MON(2) = 29                                    
                                                                        
1     IF(IH.LT.0) THEN                                                  
         IH = IH+24                                                     
         ID = ID-1                                                      
         IF(ID.EQ.0) THEN                                               
            IM = IM-1                                                   
            IF(IM.EQ.0) THEN                                            
               IM = 12                                                  
               IY = IY-1                                                
               IF(IY.LT.0) IY = 99                                      
            ENDIF                                                       
            ID = MON(IM)                                                
         ENDIF                                                          
         GOTO 1                                                         
      ELSEIF(IH.GE.24) THEN                                             
         IH = IH-24                                                     
         ID = ID+1                                                      
         IF(ID.GT.MON(IM)) THEN                                         
            ID = 1                                                      
            IM = IM+1                                                   
            IF(IM.GT.12) THEN                                           
               IM = 1                                                   
               IY = IY+1
               IF(IY.EQ.100) IY = 00
            ENDIF                                                       
         ENDIF                                                          
         GOTO 1                                                         
      ENDIF                                                             
                                                                        
      JDATE = IY*1000000 + IM*10000 + ID*100 + IH                       
                                                                        
      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C  CONVERT AN INTEGER DESCRIPTOR TO FIVE OR SIX CHARACTER ASCII FORMAT  
C---------------------------------------------------------------------- 
      FUNCTION ADN30(IDN,L30)                                           
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*(*) ADN30                                               
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(LEN(ADN30).LT.L30         ) GOTO 900                           
      IF(IDN.LT.0 .OR. IDN.GT.65535) GOTO 901                           
      IF(L30.EQ.5) THEN                                                 
         WRITE(ADN30,'(I5)') IDN                                        
      ELSEIF(L30.EQ.6) THEN                                             
         IDF = ISHFT(IDN,-14)                                           
         IDX = ISHFT(ISHFT(IDN,NBITW-14),-(NBITW-6))                    
         IDY = ISHFT(ISHFT(IDN,NBITW- 8),-(NBITW-8))                    
         WRITE(ADN30,'(I1,I2,I3)') IDF,IDX,IDY                          
      ELSE                                                              
         GOTO 902                                                       
      ENDIF                                                             
                                                                        
      DO I=1,L30                                                        
      IF(ADN30(I:I).EQ.' ') ADN30(I:I) = '0'                            
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('ADN30 - FUNCTION RETURN STRING TOO SHORT')            
901   CALL BORT('ADN30 - IDN OUT OF RANGE                ')            
902   CALL BORT('ADN30 - CHARACTER LENGTH L30 <> 5 OR 6')              
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE BFRINI                                                 

      PARAMETER (NFILES=32)
                                                                        
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /PADESC/ IBCT,IPD1,IPD2,IPD3,IPD4                          
      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)             
      COMMON /STBFR / IOLUN(32),IOMSG(32)                               
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /BUFRMG/ MSGLEN,MSGTXT(5000)                               
      COMMON /MRGCOM/ NRPL,NMRG,NAMB,NTOT                               
      COMMON /DATELN/ LENDAT
      COMMON /QUIET / IPRT                                              
      common /ACMODE/ iac
                                                                        
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*10  TAG                                                 
      CHARACTER*6   ADSN(5,2),DNDX(25,10)                               
      CHARACTER*3   TYPX(5,2),TYPS,TYP                                  
      CHARACTER*1   REPX(5,2),REPS                                      
      DIMENSION     NDNDX(10),NLDXA(10),NLDXB(10),NLDXD(10),NLD30(10)   
      DIMENSION     LENX(5)                                             
                                                                        
      DATA ADSN   / '101000','360001','360002','360003','360004' ,      
     .              '101255','031002','031001','031001','031000' /      
      DATA TYPX   /    'REP',   'DRP',   'DRP',   'DRS' ,  'DRB' ,      
     .                 'SEQ',   'RPC',   'RPC',   'RPS' ,  'SEQ' /      
      DATA REPX   /      '"',     '(',     '{',     '[' ,    '<' ,      
     .                   '"',     ')',     '}',     ']' ,    '>' /      
      DATA LENX   /       0 ,     16 ,      8 ,      8  ,     1  /      
                                                                        
      DATA (DNDX(I,1),I=1,25)/                                          
     .'102000','031001','000001','000002',                              
     .'110000','031001','000010','000011','000012','000013','000015',   
     .                  '000016','000017','000018','000019','000020',   
     .'107000','031001','000010','000011','000012','000013','101000',   
     .                  '031001','000030'/                              
                                                                        
      DATA (DNDX(I,2),I=1,15)/                                          
     .'103000','031001','000001','000002','000003',                     
     .'101000','031001','300004',                                       
     .'105000','031001','300003','205064','101000','031001','000030'/   
                                                                        
      DATA NDNDX /  25 ,  15 , 8*0 /                                    
      DATA NLDXA /  35 ,  67 , 8*0 /                                    
      DATA NLDXB /  80 , 112 , 8*0 /                                    
      DATA NLDXD /  38 ,  70 , 8*0 /                                    
      DATA NLD30 /   5 ,   6 , 8*0 /                                    

      SAVE
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  INITIALIZE /BITBUF/                                                  
C  -------------------                                                  
                                                                        
      MAXBYT = 9970                                                     
                                                                        
C  INITIALIZE /PADESC/                                                  
C  -------------------                                                  
                                                                        
      IBCT = IFXY('063000')                                             
      IPD1 = IFXY('102000')                                             
      IPD2 = IFXY('031001')                                             
      IPD3 = IFXY('206001')                                             
      IPD4 = IFXY('063255')                                             
                                                                        
C  INITIALIZE /STBFR/                                                   
C  ------------------                                                   
                                                                        
      DO I=1,NFILES                                                     
      IOLUN(I) = 0                                                      
      IOMSG(I) = 0                                                      
      ENDDO                                                             
                                                                        
C  INITIALIZE /REPTAB/                                                  
C  -------------------                                                  
                                                                        
      DO I=1,5                                                          
      LENS(I) = LENX(I)                                                 
      DO J=1,2                                                          
      IDNR(I,J) = IFXY(ADSN(I,J))                                       
      TYPS(I,J) = TYPX(I,J)                                             
      REPS(I,J) = REPX(I,J)                                             
      ENDDO                                                             
      ENDDO                                                             
                                                                        
C  INITIALIZE /TABABD/                                                  
C  -------------------                                                  
                                                                        
      NTBA(0) = 50                                                      
      NTBB(0) = 250                                                     
      NTBD(0) = 250                                                     
                                                                        
C  INITIALIZE /DXTAB/                                                   
C  ------------------                                                   
                                                                        
      MAXDX = MAXBYT                                                    
      IDXV  = 1                                                         
                                                                        
      DO J=1,10                                                         
      LDXA(J)  = NLDXA(J)                                               
      LDXB(J)  = NLDXB(J)                                               
      LDXD(J)  = NLDXD(J)                                               
      LD30(J)  = NLD30(J)                                               
      DXSTR(J) = '      '                                               
      NXSTR(J) = NDNDX(J)*2                                             
      DO I=1,NDNDX(J)                                                   
      I1 = I*2-1                                                        
      CALL IPKM(DXSTR(J)(I1:I1),2,IFXY(DNDX(I,J)))                      
      ENDDO                                                             
      ENDDO                                                             
                                                                        
C  INITIALIZE /TABLES/                                                  
C  -------------------                                                  
                                                                        
      MAXTAB = 15000                                                    
                                                                        
C  INITIALIZE /BUFRMG/                                                  
C  -------------------                                                  
                                                                        
      MSGLEN = 0                                                        
                                                                        
C  INITIALIZE /MRGCOM/                                                  
C  -------------------                                                  
                                                                        
      NRPL = 0                                                          
      NMRG = 0                                                          
      NAMB = 0                                                          
      NTOT = 0                                                          
                                                                        
C  INITIALIZE /DATELN/                                                  
C  -------------------                                                   
                                                                        
      IF(LENDAT.NE.10) LENDAT = 8                                                        
                                                                        
C  INITIALIZE /QUIET/                                                   
C  ------------------                                                   
                                                                        
      IPRT = 0                                                          

C  INITIALIZE /ACMODE/                                                  
C  ------------------                                                   
                                                                        
      IAC = 0                                                          
                                                                        
      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C  ENTRY POINT BORT IS REQUIRED FOR NON-CRAY SYSTEMS                   
C---------------------------------------------------------------------- 
      SUBROUTINE BORT(STR)                                             
      CHARACTER*(*) STR                                                 
      PRINT*
      PRINT*,'************************ABORT**************************'
      PRINT*,STR                                                        
      PRINT*,'************************ABORT**************************'
      PRINT*
      CALLEXIT(49)                                                          
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CHEKSTAB(LUN)                                          
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*24  UNIT                                                
      CHARACTER*8   NEMO,NEMS(250)                                      
      CHARACTER*1   TAB                                                 
      DIMENSION     IRPS(250),KNTS(250)                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  THERE MUST BE ENTRIES IN TABLES A, B, AND D                          
C  -------------------------------------------                          
                                                                        
      IF(NTBA(LUN).EQ.0) GOTO 900                                       
      IF(NTBB(LUN).EQ.0) GOTO 901                                       
      IF(NTBD(LUN).EQ.0) GOTO 902                                       
                                                                        
C  MAKE SURE EACH TABLE A ENTRY DEFINED AS A SEQUENCE                   
C  --------------------------------------------------                   
                                                                        
      DO I=1,NTBA(LUN)                                                  
      NEMO = TABA(I,LUN)(4:11)                                          
      CALL NEMTAB(LUN,NEMO,IDN,TAB,IRET)                                
      IF(TAB.NE.'D') GOTO 903                                           
      ENDDO                                                             
                                                                        
C  CHECK TABLE B CONTENTS                                               
C  ----------------------                                               
                                                                        
      DO ITAB=1,NTBB(LUN)                                               
      CALL NEMTBB(LUN,ITAB,UNIT,ISCL,IREF,IBIT)                         
      ENDDO                                                             
                                                                        
C  CHECK TABLE D CONTNETS                                               
C  ----------------------                                               
                                                                        
      DO ITAB=1,NTBD(LUN)                                               
      CALL NEMTBD(LUN,ITAB,NSEQ,NEMS,IRPS,KNTS)                         
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('CHEKSTAB - EMPTY TABLE A')                            
901   CALL BORT('CHEKSTAB - EMPTY TABLE B')                            
902   CALL BORT('CHEKSTAB - EMPTY TABLE D')                            
903   CALL BORT('CHEKSTAB - NO SEQUENCE DEFINED FOR TABLE A: '//NEMO)  
      END                                                               
C---------------------------------------------------------------------- 
C  CHARACTER TRANSFER TO A STRING                                       
C---------------------------------------------------------------------- 
      SUBROUTINE CHRTRN(STR,CHR,N)                                      
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*1   CHR(N)                                              
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      DO I=1,N                                                          
      STR(I:I) = CHR(I)                                                 
      ENDDO                                                             
      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C  CHARACTER TRANSFER TO A STRING WITH EBCDIC TRANSLATION               
C---------------------------------------------------------------------- 
      SUBROUTINE CHRTRNA(STR,CHR,N)                                     
                                                                        
      COMMON /CHARAC/ IASCII,IATOE(0:255),IETOA(0:255)                  
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*1   CHR(N)                                              
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      DO I=1,N                                                          
      STR(I:I) = CHR(I)                                                 
      IF(IASCII.EQ.0) CALL IPKM(STR(I:I),1,IATOE(IUPM(STR(I:I),8)))     
      ENDDO                                                             
      RETURN                                                            
      END                                                               
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
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CLOSBF(LUNIT)                                          
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.GT.0 .AND. IM.NE.0) CALL CLOSMG(LUNIT)                      
      CALL WTSTAT(LUNIT,LUN,0,0)                                        
      CLOSE(LUNIT)                                                      
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CLOSMG(LUNIT)                                          
                                                                        
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.NE.0) CALL MSGWRT(LUNIT,MBAY(1,LUN),MBYT(LUN))              
      CALL WTSTAT(LUNIT,LUN,IL,0)                                       
                                                                        
      RETURN                                                            
900   CALL BORT('CLOSMG - FILE IS CLOSED            ')                 
901   CALL BORT('CLOSMG - FILE IS OPEN FOR INPUT    ')                 
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE CONWIN(LUN,INC1,INC2,NBMP)                             
                                                                        
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  SPECIAL CASES                                                        
C  -------------                                                        
                                                                        
      IF(NCON.EQ.0) THEN                                                
         INC1 = 1                                                       
         INC2 = NVAL(LUN)                                               
         RETURN                                                         
      ENDIF                                                             
                                                                        
      IF(INC1.GT.1 .AND. KONS(NCON).EQ.5) THEN                          
         CALL NXTWIN(LUN,INC1,INC2)                                     
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  EVALUATE CONDITIONS TO SEE IF ANY MORE CASES                         
C  --------------------------------------------                         
                                                                        
10    DO NC=1,NCON                                                      
      IF(KONS(NC).EQ.5) THEN                                            
         INC1 = INVWIN(NODC(NC),LUN,INC1,NVAL(LUN))                     
         CALL USRTPL(LUN,INC1-1,NBMP)                                   
         CALL NEWWIN(LUN,INC1,INC2)                                     
      ELSE                                                              
15       CALL GETWIN(NODC(NC),LUN,INC1,INC2)                            
         IF(INC1.EQ.0 .AND. NC.EQ.1) RETURN                             
         IF(INC1.EQ.0              ) GOTO10                             
         ICON = INVCON(NC,LUN,INC1,INC2)                                
         IF(ICON.EQ.0) GOTO 15                                          
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE COPYBF(LUNIN,LUNOT)                                    
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 SEC0                                                  
      CHARACTER*1 MOCT(24000)                                           
      DIMENSION   MBAY(5000)                                            
      EQUIVALENCE (MBAY(1),MOCT(1))                                     
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      ISEC0 = 8/NBYTW+1                                                 
      NMSG  = 0                                                         
                                                                        
C  CHECK BUFR FILE STATUSES                                             
C  ------------------------                                             
                                                                        
      CALL STATUS(LUNIN,LUN,IL,IM)                                      
      IF(IL.NE.0) GOTO 900                                              
      CALL STATUS(LUNOT,LUN,IL,IM)                                      
      IF(IL.NE.0) GOTO 901                                              
                                                                        
      REWIND(LUNIN)                                                     
      REWIND(LUNOT)                                                     
                                                                        
C  READ AND COPY A BUFR FILE ON UNIT LUNIN TO UNIT LUNOT                
C  -----------------------------------------------------                
                                                                        
1     READ(LUNIN,END=2,ERR=902) SEC0,(MBAY(I),I=ISEC0,LMSG(SEC0))       
      WRITE(LUNOT     ,ERR=903) SEC0,(MBAY(I),I=ISEC0,LMSG(SEC0))       
      GOTO 1                                                            
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
2     CLOSE(LUNIN)                                                      
      CLOSE(LUNOT)                                                      
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('COPYBF - INPUT  FILE IS CURRENTLY OPEN FOR BUFR')     
901   CALL BORT('COPYBF - OUTPUT FILE IS CURRENTLY OPEN FOR BUFR')     
902   CALL BORT('COPYBF - ERROR READING FILE    ')                     
903   CALL BORT('COPYBF - ERROR WRITING FILE    ')                     
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE COPYMG(LUNIN,LUNOT)                                    
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*8  SUBSET                                               
      CHARACTER*3  TYP                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUSES                                              
C  -----------------------                                              
                                                                        
      CALL STATUS(LUNIN,LIN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      IF(IM.LE.0) GOTO 902                                              
                                                                        
      CALL STATUS(LUNOT,LOT,IL,IM)                                      
      IF(IL.EQ.0) GOTO 903                                              
      IF(IL.LT.0) GOTO 904                                              
      IF(IM.NE.0) GOTO 905                                              
                                                                        
C  MAKE SURE BOTH FILES HAVE THE SAME TABLES                            
C  -----------------------------------------                            
                                                                        
      SUBSET = TAG(INODE(LIN))                                          
      CALL NEMTBA(LOT,SUBSET,MSGT,MSTB,INOD)                            
      IF(INODE(LIN).NE.INOD) GOTO 906                                   
                                                                        
C  EVERYTHING OKAY, COPY A MESSAGE                                      
C  -------------------------------                                      
                                                                        
      MBYM = IUPB(MBAY(1,LIN),5,24)                                     
      CALL MSGWRT(LUNOT,MBAY(1,LIN),MBYM)                               
                                                                        
C  SET THE MESSAGE CONTROL WORLDS FOR LUNOT                             
C  ----------------------------------------                             
                                                                        
      NMSG (LOT) = NMSG(LOT) + 1                                        
      NSUB (LOT) = MSUB(LIN)                                            
      IDATE(LOT) = IDATE(LIN)                                           
      INODE(LOT) = INODE(LIN)                                           
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('COPYMG - INPUT FILE IS CLOSED                    ')   
901   CALL BORT('COPYMG - INPUT FILE IS OPEN FOR OUTPUT           ')   
902   CALL BORT('COPYMG - NO INPUT FILE MESSAGE OPEN              ')   
903   CALL BORT('COPYMG - OUTPUT FILE IS CLOSED                   ')   
904   CALL BORT('COPYMG - OUTPUT FILE IS OPEN FOR OUTPUT          ')   
905   CALL BORT('COPYMG - OUTPUT FILE MESSAGE OPEN                ')   
906   CALL BORT('COPYMG - INPUT/OUTPUT FILES HAVE DIFFERENT TABLES')   
      END                                                               
C-----------------------------------------------------------------------
C  COPY A SUBSET IF LUNOT>0 OTHERWISE JUST ADVANCE THE INPUT SUBSET PTR 
C-----------------------------------------------------------------------
      SUBROUTINE COPYSB(LUNIN,LUNOT,IRET)                               
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND(STATUS,CPYUPD)                                             
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
                                                                        
C  CHECK THE FILE STATUSES                                              
C  -----------------------                                              
                                                                        
      CALL STATUS(LUNIN,LIN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      IF(IM.LE.0) GOTO 902                                              
                                                                        
      IF(LUNOT.GT.0) THEN                                               
         CALL STATUS(LUNOT,LOT,IL,IM)                                   
         IF(IL.EQ.0) GOTO 903                                           
         IF(IL.LT.0) GOTO 904                                           
         IF(IM.EQ.0) GOTO 905                                           
         IF(INODE(LIN).NE.INODE(LOT)) GOTO 906                          
      ENDIF                                                             
                                                                        
C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE                        
C  ---------------------------------------------                        
                                                                        
      IF(NSUB(LIN).EQ.MSUB(LIN)) THEN                                   
         IRET = -1                                                      
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  COPY THE SUBSET TO THE OUTPUT MESSAGE AND/OR RESET THE POINTERS      
C  ---------------------------------------------------------------      
                                                                        
      IBIT = (MBYT(LIN))*8                                              
      CALL UPB(NBYT,16,MBAY(1,LIN),IBIT)                                
      IF(LUNOT.GT.0) CALL CPYUPD(LUNOT,LIN,LOT,NBYT)                    
      MBYT(LIN) = MBYT(LIN) + NBYT                                      
      NSUB(LIN) = NSUB(LIN) + 1                                         
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('COPYSB - INPUT FILE IS CLOSED                    ')   
901   CALL BORT('COPYSB - INPUT FILE IS OPEN FOR OUTPUT           ')   
902   CALL BORT('COPYSB - NO INPUT FILE MESSAGE OPEN              ')   
903   CALL BORT('COPYSB - OUTPUT FILE IS CLOSED                   ')   
904   CALL BORT('COPYSB - OUTPUT FILE IS OPEN FOR INPUT           ')   
905   CALL BORT('COPYSB - NO OUTPUT FILE MESSAGE OPEN             ')   
906   CALL BORT('COPYSB - INPUT/OUTPUT FILES HAVE DIFFERENT TABLES')   
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CPBFDX(LUD,LUN)                                        
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                

      SAVE
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  INITIALIZE THE DX-TABLE PARTITION                                    
C  ---------------------------------                                    
                                                                        
      CALL DXINIT(LUN,0)                                                
                                                                        
C  COPY ONE TABLE PARTITION TO ANOTHER                                  
C  -----------------------------------                                  
                                                                        
      INODE(LUN) = INODE(LUD)                                           
                                                                        
      NTBA(LUN) = NTBA(LUD)                                             
      NTBB(LUN) = NTBB(LUD)                                             
      NTBD(LUN) = NTBD(LUD)                                             
                                                                        
      DO I=1,NTBA(LUD)                                                  
      IDNA(I,LUN,1) = IDNA(I,LUD,1)                                     
      IDNA(I,LUN,2) = IDNA(I,LUD,2)                                     
      TABA(I,LUN) = TABA(I,LUD)                                         
      MTAB(I,LUN) = MTAB(I,LUD)                                         
      ENDDO                                                             
                                                                        
      DO I=1,NTBB(LUD)                                                  
      IDNB(I,LUN) = IDNB(I,LUD)                                         
      TABB(I,LUN) = TABB(I,LUD)                                         
      ENDDO                                                             
                                                                        
      DO I=1,NTBD(LUD)                                                  
      IDND(I,LUN) = IDND(I,LUD)                                         
      TABD(I,LUN) = TABD(I,LUD)                                         
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CPYMEM(LUNOT)                                          
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*8  SUBSET                                               
      CHARACTER*3  TYP                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUSES                                              
C  -----------------------                                              
                                                                        
      CALL STATUS(MUNIT,LIN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      IF(IM.LE.0) GOTO 902                                              
                                                                        
      CALL STATUS(LUNOT,LOT,IL,IM)                                      
      IF(IL.EQ.0) GOTO 903                                              
      IF(IL.LT.0) GOTO 904                                              
      IF(IM.NE.0) GOTO 905                                              
                                                                        
C  MAKE SURE BOTH FILES HAVE THE SAME TABLES                            
C  -----------------------------------------                            
                                                                        
      SUBSET = TAG(INODE(LIN))                                          
      CALL NEMTBA(LOT,SUBSET,MTYP,MSBT,INOD)                            
      IF(INODE(LIN).NE.INOD) GOTO 906                                   
                                                                        
C  EVERYTHING OKAY, COPY A MESSAGE                                      
C  -------------------------------                                      
                                                                        
      MBYM = IUPB(MBAY(1,LIN),5,24)                                     
      CALL MSGWRT(LUNOT,MBAY(1,LIN),MBYM)                               
                                                                        
C  SET THE MESSAGE CONTROL WORLDS FOR LUNOT                             
C  ----------------------------------------                             
                                                                        
      NMSG (LOT) = NMSG(LOT) + 1                                        
      NSUB (LOT) = MSUB(LIN)                                            
      MSUB (LOT) = MSUB(LIN)                                            
      IDATE(LOT) = IDATE(LIN)                                           
      INODE(LOT) = INODE(LIN)                                           
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('CPYMEM - INPUT FILE IS CLOSED                    ')   
901   CALL BORT('CPYMEM - INPUT FILE IS OPEN FOR OUTPUT           ')   
902   CALL BORT('CPYMEM - NO INPUT FILE MESSAGE OPEN              ')   
903   CALL BORT('CPYMEM - OUTPUT FILE IS CLOSED                   ')   
904   CALL BORT('CPYMEM - OUTPUT FILE IS OPEN FOR OUTPUT          ')   
905   CALL BORT('CPYMEM - OUTPUT FILE MESSAGE OPEN                ')   
906   CALL BORT('CPYMEM - INPUT/OUTPUT FILES HAVE DIFFERENT TABLES')   
      END                                                               
C---------------------------------------------------------------------- 
C  UPDATE THE MESSAGE BUFFER WITH NEW SUBSET                            
C---------------------------------------------------------------------- 
      SUBROUTINE CPYUPD(LUNIT,LIN,LUN,IBYT)                             
                                                                        
      COMMON /MSGPTR/ NBY0,NBY1,NBY2,NBY3,NBY4,NBY5                     
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
      CHARACTER*8 CBAY                                                  
      EQUIVALENCE (CBAY,JBAY)                                           
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND(MVB,PKB)                                                   
C-----------------------------------------------------------------------
                                                                        
                                                                        
C  SEE IF THE NEW SUBSET FITS                                           
C  --------------------------                                           
                                                                        
      IF(MBYT(LUN)+IBYT.GT.MAXBYT) THEN                                 
         CALL MSGWRT(LUNIT,MBAY(1,LUN),MBYT(LUN))                       
         CALL MSGINI(LUN)                                               
      ENDIF                                                             
                                                                        
      IF(IBYT.GT.MAXBYT-MBYT(LUN)) GOTO 900                             
                                                                        
C  TRANSFER SUBSET FROM ONE MESSAGE TO THE OTHER                        
C  ---------------------------------------------                        
                                                                        
      CALL MVB(MBAY(1,LIN),MBYT(LIN)+1,MBAY(1,LUN),MBYT(LUN)-3,IBYT)    
                                                                        
C  UPDATE THE SUBSET AND BYTE COUNTERS                                  
C  --------------------------------------                               
                                                                        
      MBYT(LUN)   = MBYT(LUN)   + IBYT                                  
      NSUB(LUN)   = NSUB(LUN)   + 1                                     
                                                                        
      LBIT = (NBY0+NBY1+NBY2+4)*8                                       
      CALL PKB(NSUB(LUN),16,MBAY(1,LUN),LBIT)                           
                                                                        
      LBYT = NBY0+NBY1+NBY2+NBY3                                        
      NBYT = IUPB(MBAY(1,LUN),LBYT+1,24)                                
      LBIT = LBYT*8                                                     
      CALL PKB(NBYT+IBYT,24,MBAY(1,LUN),LBIT)                           
                                                                        
      RETURN                                                            
                                                                        
900   CALL BORT('CPYUPD - SUBSET LONGER THAN ANY POSSIBLE MESSAGE')    
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DATEBF(LUNIT,MEAR,MMON,MDAY,MOUR,IDATE)                        
                                                                        
      COMMON /DATELN/ LENDAT

      CHARACTER*26  MSTR                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IDATE = -1                                                        
                                                                        
C  SEE IF THE FILE IS ALREADY OPEN TO BUFR INTERFACE (A NO-NO)          
C  -----------------------------------------------------------          
                                                                        
      CALL STATUS(LUNIT,LUN,JL,JM)                                      
      IF(JL.NE.0) CALL BORT('DATEBF - FILE ALREADY OPEN')              
                                                                        
C  CHECK FOR NO BUFR DATA OR NO DATA AT ALL                             
C  ----------------------------------------                             
                                                                        
      REWIND LUNIT                                                      
      READ(LUNIT,END=100,ERR=100) MSTR                                  
      IF(MSTR(1:4).NE.'BUFR') GOTO 100                                  
                                                                        
C  READ TO A DATA MESSAGE AND PICK OUT THE DATE                         
C  --------------------------------------------                         
                                                                        
1     READ(LUNIT,END=100,ERR=100) MSTR                                  
      IF(ICHAR(MSTR(17:17)).EQ.11) GOTO 1                               
      MEAR = MOD(ICHAR(MSTR(21:21)),100)                                      
      MMON = ICHAR(MSTR(22:22))                                           
      MDAY = ICHAR(MSTR(23:23))                                           
      MOUR = ICHAR(MSTR(24:24))                                           
      MMIN = ICHAR(MSTR(25:25))                                           
      MCEN = MAX(0,ICHAR(MSTR(26:26))-MIN(MEAR,1))                              
                                                                        
      IF(LENDAT.EQ.10) THEN
         IDATE = MCEN*10**8+MEAR*10**6+MMON*10**4+MDAY*10**2+MOUR               
         IDATE = I4DY(IDATE)
         MEAR  = IDATE/10**6
      ELSE
         IDATE = MEAR*10**6+MMON*10**4+MDAY*10**2+MOUR               
      ENDIF

100   RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      LOGICAL FUNCTION DIGIT(STR)
      CHARACTER*(*) STR
      DIGIT = .FALSE.
      DO I=1,LEN(STR)
      IF(STR(I:I).NE.'0' .AND. STR(I:I).NE.'1' .AND. 
     .   STR(I:I).NE.'2' .AND. STR(I:I).NE.'3' .AND. 
     .   STR(I:I).NE.'4' .AND. STR(I:I).NE.'5' .AND. 
     .   STR(I:I).NE.'6' .AND. STR(I:I).NE.'7' .AND. 
     .   STR(I:I).NE.'8' .AND. STR(I:I).NE.'9') RETURN
      ENDDO
      DIGIT = .TRUE.
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DRSTPL(INOD,LUN,INV1,INV2,INVN)                        
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      REAL*8       VAL                                                  
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND (INVWIN,USRTPL,NEWWIN)                                     
C-----------------------------------------------------------------------
                                                                        
1     NODE = INOD                                                       
2     NODE = JMPB(NODE)                                                 
      IF(NODE.EQ.0) RETURN                                              
      IF(TYP(NODE).EQ.'DRS' .OR. TYP(NODE).EQ.'DRB') THEN               
         INVN = INVWIN(NODE,LUN,INV1,INV2)                              
         IF(INVN.GT.0) THEN                                             
            CALL USRTPL(LUN,INVN,1)                                     
            CALL NEWWIN(LUN,INV1,INV2)                                  
            INVN = INVWIN(INOD,LUN,INVN,INV2)                           
            IF(INVN.GT.0) RETURN                                        
            GOTO 1                                                      
         ENDIF                                                          
      ENDIF                                                             
      GOTO 2                                                            
900   CALL BORT('DRSTPL - CANT FIND NODE:'//TAG(INOD))                 
C     print*,'drstpl:',tag(inod),':',tag(node),inv1,inv2                
C     print'(5a10)',(tag(inv(i,lun)),i=inv1,inv2)                       
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DUMPBF(LUNIT,JDATE,JDUMP)                              
                                                                        
      COMMON /DATELN/ LENDAT

      DIMENSION JDATE(5),JDUMP(5)                                       
                                                                        
      CHARACTER*32  MSTR                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      DO I=1,5                                                          
      JDATE(I) = -1                                                     
      JDUMP(I) = -1                                                     
      ENDDO                                                             
                                                                        
C  SEE IF THE FILE IS ALREADY OPEN TO BUFR INTERFACE (A NO-NO)          
C  -----------------------------------------------------------          
                                                                        
      CALL STATUS(LUNIT,LUN,JL,JM)                                      
      IF(JL.NE.0) CALL BORT('DUMPBF - FILE ALREADY OPEN')              
                                                                        
C  CHECK FOR NO BUFR DATA OR NO DATA AT ALL                             
C  ----------------------------------------                             
                                                                        
      REWIND LUNIT                                                      
1     READ(LUNIT,END=100,ERR=100) MSTR                                  
      IF(MSTR(1:4).NE.'BUFR') GOTO 100                                  
      IF(ICHAR(MSTR(17:17)).EQ.11) GOTO 1                               
                                                                        
C  DUMP CENTER YY,MM,DD,HH,MM IS IN THE FIRST EMPTY MESSAGE             
C  --------------------------------------------------------             
                                                                        
      IF(ICHAR(MSTR(31:31)).EQ.0 .AND. ICHAR(MSTR(32:32)).EQ.0) THEN    
         JDATE(1) = MOD(ICHAR(MSTR(21:21)),100)                             
         JDATE(2) = ICHAR(MSTR(22:22))                                  
         JDATE(3) = ICHAR(MSTR(23:23))                                  
         JDATE(4) = ICHAR(MSTR(24:24))                                  
         JDATE(5) = ICHAR(MSTR(25:25))                                  
         MCEN     = MAX(0,ICHAR(MSTR(26:26))-MIN(JDATE(1),1))                 
      ELSE                                                              
         RETURN                                                         
      ENDIF                                                             

      IF(LENDAT.EQ.10) THEN
         JDATE(1) = I4DY(MCEN*10**8+JDATE(1)*10**6)/10**6
      ENDIF
                                                                        
C  DUMP CLOCK YY,MM,DD,HH,MM IS IN THE SECOND EMPTY MESSAGE             
C  --------------------------------------------------------             
                                                                        
      READ(LUNIT,END=100,ERR=100) MSTR                                  
                                                                        
      IF(ICHAR(MSTR(31:31)).EQ.0 .AND. ICHAR(MSTR(32:32)).EQ.0) THEN    
         JDUMP(1) = MOD(ICHAR(MSTR(21:21)),100)                             
         JDUMP(2) = ICHAR(MSTR(22:22))                                  
         JDUMP(3) = ICHAR(MSTR(23:23))                                  
         JDUMP(4) = ICHAR(MSTR(24:24))                                  
         JDUMP(5) = ICHAR(MSTR(25:25))                                  
         MCEN     = MAX(0,ICHAR(MSTR(26:26))-MIN(JDUMP(1),1))                 
      ELSE                                                              
         RETURN                                                         
      ENDIF                                                             
                                                                        
      IF(LENDAT.EQ.10) THEN
         JDUMP(1) = I4DY(MCEN*10**8+JDUMP(1)*10**6)/10**6
      ENDIF
                                                                        
100   RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DXINIT(LUN,IOI)                                        
                                                                        
      COMMON /PADESC/ IBCT,IPD1,IPD2,IPD3,IPD4                          
      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)             
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*8   INIB(6,5),INID(5)                                   
      CHARACTER*6   ADN30                                               
      CHARACTER*3   TYPS                                                
      CHARACTER*1   REPS                                                
                                                                        
      DATA INIB   /'------','BYTCNT  ','BYTES  ','+0','+0','16',        
     .             '------','BITPAD  ','NONE   ','+0','+0','1 ',        
     .             '031000','DRF1BIT ','NUMERIC','+0','+0','1 ',        
     .             '031001','DRF8BIT ','NUMERIC','+0','+0','8 ',        
     .             '031002','DRF16BIT','NUMERIC','+0','+0','16'/        
      DATA NINIB  /5/                                                   
                                                                        
      DATA INID   /'        ',                                          
     .             'DRP16BIT',                                          
     .             'DRP8BIT ',                                          
     .             'DRPSTAK ',                                          
     .             'DRP1BIT '/                                          
      DATA NINID  /5/                                                   
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CLEAR OUT A MESSAGE CONTROL WORD PARTITION                           
C  ------------------------------------------                           
                                                                        
      NMSG(LUN)  = 0                                                    
      NSUB(LUN)  = 0                                                    
      MSUB(LUN)  = 0                                                    
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
                                                                        
C  CLEAR OUT A TABLE PARTITION                                          
C  ---------------------------                                          
                                                                        
      NTBA(LUN) = 0                                                     
      DO I=1,NTBA(0)                                                    
      TABA(I,LUN) = ' '                                                 
      MTAB(I,LUN) = 0                                                   
      ENDDO                                                             
                                                                        
      NTBB(LUN) = 0                                                     
      DO I=1,NTBB(0)                                                    
      TABB(I,LUN) = ' '                                                 
      ENDDO                                                             
                                                                        
      NTBD(LUN) = 0                                                     
      DO I=1,NTBD(0)                                                    
      TABD(I,LUN) = ' '                                                 
      CALL PKTDD(I,LUN,0,IRET)                                          
      ENDDO                                                             
                                                                        
      IF(IOI.EQ.0) RETURN                                               
                                                                        
C  INITIALIZE TABLE WITH APRIORI TABLE B AND D ENTRIES                  
C  ---------------------------------------------------                  
                                                                        
      INIB(1,1) = ADN30(IBCT,6)                                         
      INIB(1,2) = ADN30(IPD4,6)                                         
                                                                        
      DO I=1,NINIB                                                      
      NTBB(LUN) = NTBB(LUN)+1                                           
      IDNB(I,LUN) = IFXY(INIB(1,I))                                     
      TABB(I,LUN)(  1:  6) = INIB(1,I)                                  
      TABB(I,LUN)(  7: 70) = INIB(2,I)                                  
      TABB(I,LUN)( 71: 94) = INIB(3,I)                                  
      TABB(I,LUN)( 95: 98) = INIB(4,I)                                  
      TABB(I,LUN)( 99:109) = INIB(5,I)                                  
      TABB(I,LUN)(110:112) = INIB(6,I)                                  
      ENDDO                                                             
                                                                        
      DO I=2,NINID                                                      
      N = NTBD(LUN)+1                                                   
      IDND(N,LUN) = IDNR(I,1)                                           
      TABD(N,LUN)(1: 6) = ADN30(IDNR(I,1),6)                            
      TABD(N,LUN)(7:70) = INID(I)                                       
      CALL PKTDD(N,LUN,IDNR(1,1),IRET)                                  
      CALL PKTDD(N,LUN,IDNR(I,2),IRET)                                  
      NTBD(LUN) = N                                                     
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DXMINI(LUN,MBAY,MBYT,MB4,MBA,MBB,MBD)                  
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*56  DXSTR                                               
      DIMENSION     MBAY(5000)                                          
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  ENTRY POINTS FOR SEPARATING TABLES A B D
C  ----------------------------------------

      MSBT = IDXV                                                       
      GOTO 1
      ENTRY DXMINA(LUN,MBAY,MBYT,MB4,MBA,MBB,MBD)                  
      MSBT = 1                                                          
      GOTO 1
      ENTRY DXMINB(LUN,MBAY,MBYT,MB4,MBA,MBB,MBD)                  
      MSBT = 2                                                          
      GOTO 1
      ENTRY DXMIND(LUN,MBAY,MBYT,MB4,MBA,MBB,MBD)                  
      MSBT = 4                                                          
1     CONTINUE
      
                                                                        
C  INITIALIZE THE MESSAGE                                               
C  ----------------------                                               
                                                                        
      MBIT = 0                                                          
      DO I=1,5000                                                       
      MBAY(I) = 0                                                       
      ENDDO                                                             
                                                                        
      IH   = 0                                                          
      ID   = 0                                                          
      IM   = 0                                                          
      IY   = 0                                                          
      MTYP = 11                                                         
      NSUB = 1                                                          
      IDXS = IDXV+1                                                     
      LDXS = NXSTR(IDXS)                                                
                                                                        
      NBY0 = 8                                                          
      NBY1 = 18                                                         
      NBY2 = 0                                                          
      NBY3 = 7 + NXSTR(IDXS) + 1                                        
      NBY4 = 7                                                          
      NBY5 = 4                                                          
      MBYT = NBY0+NBY1+NBY2+NBY3+NBY4+NBY5                              
                                                                        
      IF(MOD(NBY3,2).NE.0) GOTO 900                                     
                                                                        
C  SECTION 0                                                            
C  ---------                                                            
                                                                        
      CALL PKC('BUFR' ,  4 , MBAY,MBIT)                                 
      CALL PKB(  MBYT , 24 , MBAY,MBIT)                                 
      CALL PKB(     3 ,  8 , MBAY,MBIT)                                 
                                                                        
C  SECTION 1                                                            
C  ---------                                                            
                                                                        
      CALL PKB(  NBY1 , 24 , MBAY,MBIT)                                 
      CALL PKB(     0 ,  8 , MBAY,MBIT)                                 
      CALL PKB(     3 ,  8 , MBAY,MBIT)                                 
      CALL PKB(     7 ,  8 , MBAY,MBIT)                                 
      CALL PKB(     0 ,  8 , MBAY,MBIT)                                 
      CALL PKB(     0 ,  8 , MBAY,MBIT)                                 
      CALL PKB(  MTYP ,  8 , MBAY,MBIT)                                 
      CALL PKB(  MSBT ,  8 , MBAY,MBIT)                                 
      CALL PKB(     4 ,  8 , MBAY,MBIT)                                 
      CALL PKB(  IDXV ,  8 , MBAY,MBIT)                                 
      CALL PKB(    IY ,  8 , MBAY,MBIT)                                 
      CALL PKB(    IM ,  8 , MBAY,MBIT)                                 
      CALL PKB(    ID ,  8 , MBAY,MBIT)                                 
      CALL PKB(    IH ,  8 , MBAY,MBIT)                                 
      CALL PKB(     0 ,  8 , MBAY,MBIT)                                 
      CALL PKB(     0 ,  8 , MBAY,MBIT)                                 
                                                                        
C  SECTION 3                                                            
C  ---------                                                            
                                                                        
      CALL PKB(       NBY3 ,   24 , MBAY,MBIT)                          
      CALL PKB(          0 ,    8 , MBAY,MBIT)                          
      CALL PKB(          1 ,   16 , MBAY,MBIT)                          
      CALL PKB(       2**7 ,    8 , MBAY,MBIT)                          
      DO I=1,LDXS                                                       
      CALL PKB(IUPM(DXSTR(IDXS)(I:I),8),8,MBAY,MBIT)                    
      ENDDO                                                             
      CALL PKB(          0 ,    8 , MBAY,MBIT)                          
                                                                        
C  SECTION 4                                                            
C  ---------                                                            
                                                                        
      MB4 = MBIT/8+1                                                    
      CALL PKB(NBY4 , 24 , MBAY,MBIT)                                   
      CALL PKB(   0 ,  8 , MBAY,MBIT)                                   
      MBA = MBIT/8+1                                                    
      CALL PKB(   0 ,  8 , MBAY,MBIT)                                   
      MBB = MBIT/8+1                                                    
      CALL PKB(   0 ,  8 , MBAY,MBIT)                                   
      MBD = MBIT/8+1                                                    
      CALL PKB(   0 ,  8 , MBAY,MBIT)                                   
                                                                        
      IF(MBIT/8+NBY5.NE.MBYT) GOTO 901                                  
                                                                        
      RETURN                                                            
900   CALL BORT('DXMINI - UNEVEN SECTION 3')                           
901   CALL BORT('DXMINI - BYTCNT IS OFF   ')                           
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ELEMDX(CARD,LUN)                                       
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*80  CARD                                                
      CHARACTER*24  UNIT                                                
      CHARACTER*11  REFR                                                
      CHARACTER*8   NEMO                                                
      CHARACTER*4   SCAL                                                
      CHARACTER*3   BITW                                                
      CHARACTER*1   SIGN,TAB                                            
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CAPTURE THE VARIOUS ELEMENTS CHARACTERISTICS                         
C  --------------------------------------------                         
                                                                        
      NEMO = CARD( 3:10)                                                
      SCAL = CARD(14:17)                                                
      REFR = CARD(21:31)                                                
      BITW = CARD(35:37)                                                
      UNIT = CARD(41:64)                                                
                                                                        
C  FIND THE ELEMENT TAG IN TABLE B                                      
C  -------------------------------                                      
                                                                        
      CALL NEMTAB(LUN,NEMO,IDSN,TAB,IELE)                               
      IF(TAB.NE.'B') GOTO 900                                           
                                                                        
C  LEFT JUSTIFY AND STORE CHARACTERISTICS                               
C  --------------------------------------                               
                                                                        
      CALL JSTCHR(UNIT)                                                 
      TABB(IELE,LUN)(71:94) = UNIT                                      
                                                                        
      CALL JSTNUM(SCAL,SIGN,IRET)                                       
      IF(IRET.NE.0) GOTO 901                                            
      TABB(IELE,LUN)(95:95) = SIGN                                      
      TABB(IELE,LUN)(96:98) = SCAL                                      
                                                                        
      CALL JSTNUM(REFR,SIGN,IRET)                                       
      IF(IRET.NE.0) GOTO 902                                            
      TABB(IELE,LUN)( 99: 99) = SIGN                                    
      TABB(IELE,LUN)(100:109) = REFR                                    
                                                                        
      CALL JSTNUM(BITW,SIGN,IRET)                                       
      IF(IRET.NE.0  ) GOTO 903                                          
      IF(SIGN.EQ.'-') GOTO 903                                          
      TABB(IELE,LUN)(110:112) = BITW                                    
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXIT                                                           
C  -----------                                                          
                                                                        
900   CALL BORT('ELEMDX - UNDEFINED ELEMENT: '//CARD)                  
901   CALL BORT('ELEMDX - BAD SCALE VALUE:   '//CARD)                  
902   CALL BORT('ELEMDX - BAD REFERENCE VAL: '//CARD)                  
903   CALL BORT('ELEMDX - BAD BIT WIDTH:     '//CARD)                  
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE GETWIN(NODE,LUN,IWIN,JWIN)                             
                                                                        
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
cfpp$ expand (lstrpc)                                                   
C---------------------------------------------------------------------- 
                                                                        
      IRPC = LSTRPC(NODE,LUN)                                           
                                                                        
      IF(IRPC.EQ.0) THEN                                                
         IWIN = INVWIN(NODE,LUN,JWIN,NVAL(LUN))                         
         IF(IWIN.EQ.0 .and. jwin.gt.1) RETURN                           
         IWIN = 1                                                       
         JWIN = NVAL(LUN)                                               
         RETURN                                                         
      ELSE                                                              
         IWIN = INVWIN(IRPC,LUN,JWIN,NVAL(LUN))                         
         IF(IWIN.EQ.0) THEN                                             
            RETURN                                                      
         ELSEIF(VAL(IWIN,LUN).EQ.0.) THEN                               
            IWIN = 0                                                    
            RETURN                                                      
         ENDIF                                                          
      ENDIF                                                             
                                                                        
      JWIN = INVWIN(IRPC,LUN,IWIN+1,NVAL(LUN))                          
      IF(JWIN.EQ.0) CALL BORT('GETWIN - MISSING BRACKET')              
                                                                        
      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C  CONVERT AN 8 DIGIT INTEGER DATE (YYMMDDHH) TO 10 DIGITS (CCYYMMDDHH)
C---------------------------------------------------------------------- 
      FUNCTION I4DY(IDATE)                                              

      IF(IDATE.LT.10**8) THEN
         IY = IDATE/10**6
         IF(IY.GT.20) I4DY = IDATE + 19*10**8 
         IF(IY.LE.20) I4DY = IDATE + 20*10**8 
      ELSE
         I4DY = IDATE
      ENDIF

      RETURN
      END
C---------------------------------------------------------------------- 
C  CONVERT A FIVE OR SIX CHARACTER ASCII DESCRIPTOR TO AN INTEGER       
C---------------------------------------------------------------------- 
      FUNCTION IDN30(ADN30,L30)                                         
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*(*) ADN30                                               
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(LEN(ADN30).LT.L30) GOTO 900                                    
      IF(L30.EQ.5) THEN                                                 
         READ(ADN30,'(I5)') IDN30                                       
         IF(IDN30.LT.0 .OR. IDN30.GT.65535) GOTO 901                    
      ELSEIF(L30.EQ.6) THEN                                             
         IDN30 = IFXY(ADN30)                                            
      ELSE                                                              
         GOTO 902                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('IDN30 - FUNCTION INPUT STRING TOO SHORT    ')         
901   CALL BORT('IDN30 - IDN OUT OF RANGE, NOT A DESCRIPTOR ')         
902   CALL BORT('IDN30 - CHARACTER LENGTH L30 <> 5 OR 6     ')         
      END                                                               
C---------------------------------------------------------------------- 
C  WILL CHECK TO SEE IF ANY UNREAD SUBSETS IN AN OPEN INPUT MESSAGE     
C---------------------------------------------------------------------- 
      FUNCTION IFBGET(LUNIT)                                            
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  MAKE SURE A FILE/MESSAGE IS OPEN FOR INPUT                           
C  ------------------------------------------                           
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.GE.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
                                                                        
C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE                        
C  ---------------------------------------------                        
                                                                        
      IF(NSUB(LUN).LT.MSUB(LUN)) THEN                                   
         IFBGET = 0                                                     
      ELSE                                                              
         IFBGET = -1                                                    
      ENDIF                                                             
                                                                        
C  EXIT ONE WAY OR ANOTHER                                              
C  -----------------------                                              
                                                                        
      RETURN                                                            
900   CALL BORT('IFBGET - FILE NOT OPEN FOR INPUT')                    
901   CALL BORT('IFBGET - NO MESSAGE OPEN        ')                    
      END                                                               
C---------------------------------------------------------------------- 
C  CONVERT A SIX CHARACTER (FXY) ASCII DESCRIPTOR TO AN INTEGER         
C---------------------------------------------------------------------- 
      FUNCTION IFXY(ADSC)                                               
                                                                        
      CHARACTER*6 ADSC                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      READ(ADSC,'(I1,I2,I3)') IF,IX,IY                                  
      IFXY = IF*2**14 + IX*2**8 + IY                                    
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE INCTAB(ATAG,ATYP,NODE)                                 
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*(*) ATAG,ATYP                                           
      CHARACTER*10  TAG                                                 
      CHARACTER*3   TYP                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NTAB = NTAB+1                                                     
      IF(NTAB.GT.MAXTAB) CALL BORT('INCTAB - TOO MANY ENTRIES')        
                                                                        
      TAG(NTAB) = ATAG                                                  
      TYP(NTAB) = ATYP                                                  
      NODE = NTAB                                                       
                                                                        
      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      FUNCTION INVCON(NC,LUN,INV1,INV2)                                 
                                                                        
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE INVENTORY INTERVAL                                         
C  ----------------------------                                         
                                                                        
      IF(INV1.LE.0 .OR. INV1.GT.NVAL(LUN)) GOTO 99                      
      IF(INV2.LE.0 .OR. INV2.GT.NVAL(LUN)) GOTO 99                      
                                                                        
C  FIND AN OCCURANCE OF NODE IN THE WINDOW MEETING THIS CONDITION       
C  --------------------------------------------------------------       
                                                                        
      DO INVCON=INV1,INV2                                               
      IF(INV(INVCON,LUN).EQ.NODC(NC)) THEN                              
         IF(KONS(NC).EQ.1 .AND. VAL(INVCON,LUN).EQ.IVLS(NC)) RETURN     
         IF(KONS(NC).EQ.2 .AND. VAL(INVCON,LUN).NE.IVLS(NC)) RETURN     
         IF(KONS(NC).EQ.3 .AND. VAL(INVCON,LUN).LT.IVLS(NC)) RETURN     
         IF(KONS(NC).EQ.4 .AND. VAL(INVCON,LUN).GT.IVLS(NC)) RETURN     
      ENDIF                                                             
      ENDDO                                                             
                                                                        
99    INVCON = 0                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE INVMRG(LUBFI,LUBFJ)                                    
                                                                        
      COMMON /MRGCOM/ NRPL,NMRG,NAMB,NTOT                               
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*10  TAG                                                 
      CHARACTER*3   TYP                                                 
      LOGICAL       HEREI,HEREJ,MISSI,MISSJ,SAMEI                       
      REAL*8        VAL,BMISS                                           
                                                                        
      DATA BMISS /10E10/                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IS = 1                                                            
      JS = 1                                                            
                                                                        
C  GET THE UNIT POINTERS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUBFI,LUNI,IL,IM)                                     
      CALL STATUS(LUBFJ,LUNJ,JL,JM)                                     
                                                                        
C  STEP THROUGH THE BUFFERS COMPARING THE INVENTORY AND MERGING DATA    
C  -----------------------------------------------------------------    
                                                                        
      DO WHILE(IS.LE.NVAL(LUNI))                                        
                                                                        
C  CHECK TO SEE WE ARE AT THE SAME NODE IN EACH BUFFER                  
C  ---------------------------------------------------                  
                                                                        
      NODE = INV(IS,LUNI)                                               
      NODJ = INV(JS,LUNJ)                                               
      IF(NODE.NE.NODJ) CALL BORT('TABULAR MISMATCH')                   
                                                                        
      ITYP = ITP(NODE)                                                  
                                                                        
C  FOR TYPE 1 NODES DO AN ENTIRE SEQUENCE REPLACEMENT                   
C  --------------------------------------------------                   
                                                                        
      IF(ITYP.EQ.1) THEN                                                
         IF(TYP(NODE).EQ.'DRB') IOFF = 0                                
         IF(TYP(NODE).NE.'DRB') IOFF = 1                                
         IWRDS = NWORDS(IS,LUNI)+IOFF                                   
         JWRDS = NWORDS(JS,LUNJ)+IOFF                                   
         IF(IWRDS.GT.IOFF .AND. JWRDS.EQ.IOFF) THEN                     
CDIR$       IVDEP                                                       
            DO N=NVAL(LUNJ),JS+1,-1                                     
            INV(N+IWRDS-JWRDS,LUNJ) = INV(N,LUNJ)                       
            VAL(N+IWRDS-JWRDS,LUNJ) = VAL(N,LUNJ)                       
            ENDDO                                                       
            DO N=0,IWRDS                                                
            INV(JS+N,LUNJ) = INV(IS+N,LUNI)                             
            VAL(JS+N,LUNJ) = VAL(IS+N,LUNI)                             
            ENDDO                                                       
            NVAL(LUNJ) = NVAL(LUNJ)+IWRDS-JWRDS                         
            JWRDS = IWRDS                                               
            NRPL = NRPL+1                                               
         ENDIF                                                          
         IS = IS+IWRDS                                                  
         JS = JS+JWRDS                                                  
      ENDIF                                                             
                                                                        
C  FOR TYPES 2 AND 3 FILL MISSINGS                                      
C  -------------------------------                                      
                                                                        
      IF(ITYP.EQ.2) THEN                                                
         HEREI = VAL(IS,LUNI).LT.BMISS                                  
         HEREJ = VAL(JS,LUNJ).LT.BMISS                                  
         MISSI = VAL(IS,LUNI).GE.BMISS                                  
         MISSJ = VAL(JS,LUNJ).GE.BMISS                                  
         SAMEI = VAL(IS,LUNI).EQ.VAL(JS,LUNJ)                           
         IF(HEREI.AND.MISSJ) THEN                                       
            VAL(JS,LUNJ) = VAL(IS,LUNI)                                 
            NMRG = NMRG+1                                               
         ELSEIF(HEREI.AND.HEREJ.AND..NOT.SAMEI) THEN                    
            NAMB = NAMB+1                                               
         ENDIF                                                          
      ENDIF                                                             
                                                                        
      IF(ITYP.EQ.3) THEN                                                
         HEREI = VAL(IS,LUNI).NE.BMISS                                  
         HEREJ = VAL(JS,LUNJ).NE.BMISS                                  
         MISSI = VAL(IS,LUNI).EQ.BMISS                                  
         MISSJ = VAL(JS,LUNJ).EQ.BMISS                                  
         SAMEI = VAL(IS,LUNI).EQ.VAL(JS,LUNJ)                           
         IF(HEREI.AND.MISSJ) THEN                                       
            VAL(JS,LUNJ) = VAL(IS,LUNI)                                 
            NMRG = NMRG+1                                               
         ELSEIF(HEREI.AND.HEREJ.AND..NOT.SAMEI) THEN                    
            NAMB = NAMB+1                                               
         ENDIF                                                          
      ENDIF                                                             
                                                                        
C  BUMP THE COUNTERS AND GO CHECK THE NEXT PAIR                         
C  --------------------------------------------                         
                                                                        
      IS = IS + 1                                                       
      JS = JS + 1                                                       
      ENDDO                                                             
                                                                        
      NTOT = NTOT+1                                                     
                                                                        
      RETURN                                                            
                                                                        
C  ENTRY MRGINV PRINTS A SUMMARY OF MERGE ACTIVITY                      
C  -----------------------------------------------                      
                                                                        
      ENTRY MRGINV                                                      
                                                                        
      PRINT*                                                            
      PRINT*,'-----------------------------'                            
      PRINT*,'INVENTORY FROM INVMRG PROCESS'                            
      PRINT*,'-----------------------------'                            
      PRINT*,'DRB EXPANSION= ',NRPL                                     
      PRINT*,'MERGES       = ',NMRG                                     
      PRINT*,'AMBIGUOUS    = ',NAMB                                     
      PRINT*,'-----------------------------'                            
      PRINT*,'TOTAL VISITS = ',NTOT                                     
      PRINT*,'-----------------------------'                            
      PRINT*                                                            
                                                                        
      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      FUNCTION INVTAG(NODE,LUN,INV1,INV2)                               
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG,tagn                                             
      CHARACTER*3  TYP                                                  
      REAL*8       VAL                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      INVtag = 0                                                        
      IF(NODE.EQ.0) RETURN                                              
      tagn = tag(node)                                                  
                                                                        
C  SEARCH BETWEEN INV1 AND INV2                                         
C  ----------------------------                                         
                                                                        
10    DO INVtag=INV1,INV2                                               
      IF(tag(INV(INVtag,LUN)).EQ.tagn) RETURN                           
      ENDDO                                                             
                                                                        
      INVtag = 0                                                        
      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      FUNCTION INVWIN(NODE,LUN,INV1,INV2)                               
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      REAL*8       VAL                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(NODE.EQ.0) RETURN                                              
                                                                        
C  SEARCH BETWEEN INV1 AND INV2                                         
C  ----------------------------                                         
                                                                        
10    DO INVWIN=INV1,INV2                                               
      IF(INV(INVWIN,LUN).EQ.NODE) RETURN                                
      ENDDO                                                             
                                                                        
      INVWIN = 0                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C  PACK AN INTEGER INTO A CHARACTER ARRAY
C----------------------------------------------------------------------
      SUBROUTINE IPKM(CBAY,NBYT,N)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*8 CBAY,CINT
      EQUIVALENCE(CINT,INT)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      IF(NBYT.GT.NBYTW) CALL BORT('IPKM - NBYT>WRD LEN')
      INT = IREV(ISHFT(N,(NBYTW-NBYT)*8))
      DO I=1,NBYT
      CBAY(I:I) = CINT(I:I)
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION IRDERM(LUNIT,MBAY)                                       
                                                                        
      CHARACTER*4  SEVN                                            
      CHARACTER*1  BAY(5000*8)
      CHARACTER*36 SEC013,SECSAV
      DIMENSION    MBAY(5000),KBAY(5000)                                
      EQUIVALENCE  (BAY(1),KBAY(1),SEC013)                       
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      DO I=1,5000                                                       
      KBAY(I) = 0                                                       
      ENDDO                                                             
      IBSKIP = 0
                                                                        
C  FIND A BUFR MESSAGE                                                  
C  -------------------                                                  
                                                                        
1     READ(LUNIT,END=100,ERR=100) SEC013(1:32)
      LASTRD = 32
2     IBUFR = INDEX(SEC013(1:32),'BUFR')
      IF(IBUFR.EQ.0) THEN                                  
         IBSKIP = IBSKIP + LASTRD
         IF(INDEX(SEC013(1:32),'B').EQ.0) THEN
            READ(LUNIT,END=100,ERR=100) SEC013(1:32)                    
            LASTRD = 32
         ELSE IF(SEC013(32:32).EQ.'B') THEN
            SEC013(1:1) = 'B'
            READ(LUNIT,END=100,ERR=100) SEC013(2:32)
            LASTRD = 31
         ELSE IF(INDEX(SEC013(1:32),'BU').EQ.0) THEN
            READ(LUNIT,END=100,ERR=100) SEC013(1:32)
            LASTRD = 32
         ELSE IF(SEC013(31:32).EQ.'BU') THEN
            SEC013(1:2) = 'BU'
            READ(LUNIT,END=100,ERR=100) SEC013(3:32)
            LASTRD = 30
         ELSE IF(INDEX(SEC013(1:32),'BU').LT.30) THEN
            READ(LUNIT,END=100,ERR=100) SEC013(1:32)                    
            LASTRD = 32
         ELSE IF(SEC013(30:32).EQ.'BUF') THEN
            SEC013(1:3) = 'BUF'
            READ(LUNIT,END=100,ERR=100) SEC013(4:32)
            LASTRD = 29
         ELSE
            READ(LUNIT,END=100,ERR=100) SEC013(1:32)
            LASTRD = 32
         ENDIF
         GOTO 2                                                         
      ELSE IF(IBUFR.GT.1) THEN
         SECSAV = SEC013
         SEC013(1:33-IBUFR) = SECSAV(IBUFR:32)
         READ(LUNIT,END=100,ERR=100) SEC013(34-IBUFR:32)
         IBSKIP = IBSKIP + IBUFR - 1
      ENDIF                                                             
                                                                        
C  IF THIS IS BUFR RELEASE 0, THE FOLLOWING WILL ACCOUNT FOR SECTION-1
C  -------------------------------------------------------------------

      J = IUPM(BAY(5),24)
      I = J + 8

      IF(J.LE.32) THEN

C  DETERMINE WHETHER BYTE COUNT IS FOR SECTION-0 OR SECTION-1
C  ----------------------------------------------------------

         SEVN = SEC013(J-3:J)
         IF(SEVN.EQ.'7777') THEN
            IF(J.LT.32) THEN
               SECSAV = SEC013
               SEC013(1:32-J) = SECSAV(J+1:32)
            ENDIF
            READ(LUNIT,END=100,ERR=100) SEC013(33-J:32)
            PRINT*,'SHORT RECORD SKIPPED'                    
            GOTO 2
         ENDIF

C  IF THIS IS BUFR RELEASE 0, SHIFT BYTES 5 AND UP 4 BYTES TO THE RIGHT 
C  --------------------------------------------------------------------

         SECSAV = SEC013
         SEC013(9:36) = SECSAV(5:32)

C  IF SECTION-2 IS ABSENT, BEGIN WITH SECTION-3
C  --------------------------------------------

         IF(IUPM(BAY(16),1).GT.0) THEN         
            KCONT = 2
         ELSE
            KCONT = 3
         ENDIF

C  DETERMINE WHETHER SECTION-2 AND SECTION-3 ARE IN BYTES 9-36
C  -----------------------------------------------------------

         KREST = 0
         DO K=KCONT,3
         IF(I.LE.33) THEN
            I = I + IUPM(BAY(I+1),24)                                
            IF(I.GT.36) THEN
               READ(LUNIT,END=100,ERR=100) (BAY(J),J=37,I)
            ENDIF
            KREST = K + 1
         ELSE IF(I.LE.35) THEN
            READ(LUNIT,END=100,ERR=100) (BAY(J),J=37,I+3),
     .                  (BAY(J),J=I+4,I+IUPM(BAY(I+1),24))           
            I = I + IUPM(BAY(I+1),24)                                
            KREST = K + 1
         ENDIF                                                       
         ENDDO                                                          
         IF(KREST.NE.0) KCONT = KREST
      ELSE
         READ(LUNIT,END=100,ERR=100) (BAY(K),K=33,J)                

C  DETERMINE WHETHER BYTE COUNT IS FOR SECTION-0 OR SECTION-1
C  ----------------------------------------------------------

         SEVN = BAY(J-3)//BAY(J-2)//BAY(J-1)//BAY(J)
         IF(SEVN.EQ.'7777') GOTO 50

C  IF THIS IS BUFR RELEASE 0, SHIFT BYTES 5 AND UP 4 BYTES TO THE RIGHT 
C  --------------------------------------------------------------------

         READ(LUNIT,END=100,ERR=100) (BAY(K),K=J+5,J+8)                
         DO K=J,33,-1
         BAY(K+4) = BAY(K)
         ENDDO
         SECSAV = SEC013
         SEC013(9:36) = SECSAV(5:32)

C  IF SECTION-2 IS ABSENT, BEGIN WITH SECTION-3
C  --------------------------------------------

         IF(IUPM(BAY(16),1).GT.0) THEN         
            KCONT = 2
         ELSE
            KCONT = 3
         ENDIF

      ENDIF

C  FOR REMAINING SECTIONS (UP TO SECTION-4) READ BYTE COUNT AND BYTES
C  ------------------------------------------------------------------

      DO K=KCONT,4                                                      
      READ(LUNIT,END=100,ERR=100) (BAY(J),J=I+1,I+3),                
     .            (BAY(J),J=I+4,I+IUPM(BAY(I+1),24))                 
      I = I + IUPM(BAY(I+1),24)                                      
      ENDDO                                                             
                                                                        
C  CHECK ON SECTION 5 FOR BAD RECORD INDICATOR                          
C  -------------------------------------------                          
                                                                        
      READ(LUNIT,END=100,ERR=100) SEVN                                  
      IF(SEVN.NE.'7777') THEN
         PRINT*,'BAD RECORD SKIPPED'                    
         GOTO 1                                         
      ENDIF
      I = I+4                                                           
                                                                        
C  FILL IN THE ARRAY TO RETURN                                          
C  ---------------------------                                          
                                                                        
50    DO I=1,5000                                                       
      MBAY(I) = KBAY(I)                                                 
      ENDDO                                                             
                                                                        
      IRDERM =  0                                                       
      RETURN                                                            
                                                                        
100   IRDERM = -1                                                       
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION IREADERM(LUNIT,SUBSET,IDATE)                             
      CHARACTER*8 SUBSET                                                
      CALL READERM(LUNIT,SUBSET,IDATE,IRET)                         
      IREADERM = IRET
      RETURN                                                            
                                                                        
      ENTRY IREADERS(LUNIT)                                             
      CALL READERS(LUNIT,IRET)                                      
      IREADERS = IRET
      RETURN                                                            
                                                                        
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION IREADMG(LUNIT,SUBSET,IDATE)                              
      CHARACTER*8 SUBSET                                                
      CALL READMG(LUNIT,SUBSET,IDATE,IRET)                           
      IREADMG = IRET
      RETURN                                                            
                                                                        
      ENTRY IREADMM(IMSG,SUBSET,IDATE)                                 
      CALL READMM(IMSG,SUBSET,IDATE,IRET)                         
      IREADMM = IRET
      RETURN                                                            

      ENTRY IREADFT(LUNIT,SUBSET,IDATE)                                 
      CALL READFT(LUNIT,SUBSET,IDATE,IRET)                           
      IREADFT = IRET
      RETURN                                                            
                                                                        
      ENTRY IREADNS(LUNIT,SUBSET,IDATE)                                 
      CALL READNS(LUNIT,SUBSET,IDATE,IRET)                           
      IREADNS = IRET
      RETURN                                                            
                                                                        
      ENTRY IREADSB(LUNIT)                                              
      CALL READSB(LUNIT,IRET)                                        
      IREADSB = IRET
      RETURN                                                            
                                                                        
      ENTRY ICOPYSB(LUNIN,LUNOT)                                        
      CALL COPYSB(LUNIN,LUNOT,IRET)                                  
      ICOPYSB = IRET
      RETURN                                                            
                                                                        
      ENTRY IREADIBM(LUNIT,SUBSET,IDATE)                              
      CALL READIBM(LUNIT,SUBSET,IDATE,IRET)                           
      IREADIBM = IRET
      RETURN                                                            

      END                                                               
C----------------------------------------------------------------------
C  REVERSE A WIERD INTEGER
C----------------------------------------------------------------------
      FUNCTION IREV(N)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*8 CINT,DINT
      EQUIVALENCE(CINT,INT)
      EQUIVALENCE(DINT,JNT)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
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
C----------------------------------------------------------------------
C  UNPACK UP AN INTEGER FROM A PACKED INTEGER ARRAY (FUNCTION)
C----------------------------------------------------------------------
      FUNCTION IUPB(MBAY,NBYT,NBIT)
 
      DIMENSION MBAY(*)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      MBIT = (NBYT-1)*8
      CALL UPB(IRET,NBIT,MBAY,MBIT)
      IUPB = IRET
      RETURN
      END
C-----------------------------------------------------------------------
C  UNPACK UP AN INTEGER FROM A PACKED CHARACTER ARRAY (FUNCTION)
C----------------------------------------------------------------------
      FUNCTION IUPM(CBAY,NBITS)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*8 CBAY,CINT
      EQUIVALENCE(CINT,INT)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      IF(NBITS.GT.NBITW) CALL BORT('IUPM - NBITS>WRD LEN')
      CINT = CBAY
      INT  = IREV(INT)
      IUPM = ISHFT(INT,NBITS-NBITW)
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE JSTIFY                                                 
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*1  SIGN                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      ENTRY JSTCHR(STR)                                                 
                                                                        
      LSTR = LEN(STR)                                                   
                                                                        
      IF(STR.EQ.' ') GOTO 900                                           
1     IF(STR(1:1).EQ.' ') THEN                                          
         STR  = STR(2:LSTR)                                             
         GOTO 1                                                         
      ENDIF                                                             
      RETURN                                                            
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      ENTRY JSTNUM(STR,SIGN,IRET)                                       
                                                                        
      IRET = 0                                                          
      LSTR = LEN(STR)                                                   
                                                                        
      IF(STR.EQ.' ') GOTO 900                                           
2     IF(STR(1:1).EQ.' ') THEN                                          
         STR  = STR(2:LSTR)                                             
         GOTO 2                                                         
      ENDIF                                                             
      IF(STR(1:1).EQ.'+') THEN                                          
         STR  = STR(2:LSTR)                                             
         SIGN = '+'                                                     
      ELSEIF(STR(1:1).EQ.'-') THEN                                      
         STR  = STR(2:LSTR)                                             
         SIGN = '-'                                                     
      ELSE                                                              
         SIGN = '+'                                                     
      ENDIF                                                             
                                                                        
      CALL STRNUM(STR,NUM)                                              
      IF(NUM.LT.0) IRET = -1                                            
      RETURN                                                            
                                                                        
900   CALL BORT('JSTIFY - BLANK STRING NOT ALLOWED')                   
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION LJUST(STR)                                               
      CHARACTER*(*) STR                                                 
      LJUST = 0                                                         
      IF(STR.EQ.' ') RETURN                                             
      LSTR = LEN(STR)                                                   
      DO I=1,LSTR
      DO WHILE(STR(I:I).EQ.' ' .AND. STR(I+1:LSTR).NE.' ')
         STR(I:LSTR) = STR(I+1:LSTR)                                      
      ENDDO                                                             
      ENDDO                                                             
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C  UNPACK BUFR MESSAGE LENGTH                                           
C-----------------------------------------------------------------------
      FUNCTION LMSG(SEC0)                                               
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 SEC0                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IMSG = 8/NBYTW                                                    
      LMSG = IUPM(SEC0(5:7),24)/8                                       
      IF(LMSG.EQ.0) RETURN                                              
      LMSG = (LMSG+1)*IMSG                                              
      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      FUNCTION LSTJPB(NODE,LUN,JBTYP)                                   
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*(*) JBTYP                                               
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  MAKE SURE WE ARE ALL ON THE SAME PAGE                                
C  -------------------------------------                                
                                                                        
      IF(NODE.LT.INODE(LUN) .OR. NODE.GT.ISC(INODE(LUN))) THEN          
         PRINT*,INODE(LUN),':',NODE,':',TAG(NODE)                       
         CALL BORT('LSTJPB - TABLE NODE IS OUT OF BOUNDS')             
      ENDIF                                                             
                                                                        
C  FIND THIS OR THE PREVIOUS RPC NODE                                   
C  ----------------------------------                                   
                                                                        
      LSTJPB = NODE                                                     
                                                                        
10    IF(TYP(LSTJPB).NE.JBTYP) THEN                                     
         LSTJPB = JMPB(LSTJPB)                                          
         IF(LSTJPB.NE.0) GOTO 10                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      FUNCTION LSTRPC(NODE,LUN)                                         
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(NODE.LT.INODE(LUN) .OR. NODE.GT.ISC(INODE(LUN))) GOTO 900      
                                                                        
      NOD = NODE                                                        
                                                                        
C  FIND THIS OR THE PREVIOUS RPC NODE                                   
C  ----------------------------------                                   
                                                                        
10    IF(TYP(NOD).NE.'RPC') THEN                                        
         NOD = JMPB(NOD)                                                
         IF(NOD.NE.0) GOTO 10                                           
      ENDIF                                                             
                                                                        
      LSTRPC = NOD                                                      
                                                                        
      RETURN                                                            
900   PRINT*,INODE(LUN),':',NODE                                        
      CALL BORT('LSTRPC - TABLE NODE IS OUT OF BOUNDS')                
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      FUNCTION LSTRPS(NODE,LUN)                                         
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(NODE.LT.INODE(LUN) .OR. NODE.GT.ISC(INODE(LUN))) GOTO 900      
                                                                        
      NOD = NODE                                                        
                                                                        
C  FIND THIS OR THE PREVIOUS RPS NODE                                   
C  ----------------------------------                                   
                                                                        
10    IF(TYP(NOD).NE.'RPS') THEN                                        
         NOD = JMPB(NOD)                                                
         IF(NOD.NE.0) GOTO 10                                           
      ENDIF                                                             
                                                                        
      LSTRPS = NOD                                                      
                                                                        
      RETURN                                                            
900   CALL BORT('LSTRPS - TABLE NODE IS OUT OF BOUNDS')                
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MAKESTAB                                               
                                                                        
      COMMON /QUIET/  IPRT                                              
      COMMON /STBFR/  IOLUN(32),IOMSG(32)                               
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*10  TAG                                                 
      CHARACTER*8   NEMO                                                
      CHARACTER*3   TYP                                                 
      DIMENSION     LUS(32)                                             
      LOGICAL       EXPAND,PRTTAB                                       
      REAL*8        VAL,BMISS                                           
                                                                        
      DATA PRTTAB /.FALSE./                                            
      DATA BMISS  /10E10/                                               
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      PRTTAB = IPRT.GE.2                                                
                                                                        
C  RESET POINTER TABLE AND STRING CACHE                                 
C  ------------------------------------                                 
                                                                        
      NTAB = 0                                                          
      CALL STRCLN                                                       
                                                                        
C  FIGURE OUT WHICH UNITS SHARE TABLES                                  
C  -----------------------------------                                  
                                                                        
      DO LUN=1,32                                                       
      LUS(LUN) = 0                                                      
      IF(IOLUN(LUN).NE.0) THEN                                          
         DO LUM=1,LUN-1                                                 
         IF(MTAB(1,LUN).EQ.MTAB(1,LUM)) LUS(LUN) = LUM                  
         ENDDO                                                          
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  INITIALIZE JUMP-LINK TABLES WITH SUBSETS/SEQUENCES/ELEMENTS          
C  ----------------------------------------------------------           
                                                                        
      DO LUN=1,32                                                       
                                                                        
      IF(IOLUN(LUN).NE.0) THEN                                          
                                                                        
C  RESET ANY EXISTING INVENTORY POINTERS                                
C  -------------------------------------                                
                                                                        
         IF(IOMSG(LUN).NE.0) THEN                                       
            IF(LUS(LUN).EQ.0) INC = (NTAB+1)-MTAB(1,LUN)                
            IF(LUS(LUN).NE.0) INC = MTAB(1,LUS(LUN))-MTAB(1,LUN)        
            DO N=1,NVAL(LUN)                                            
            INV(N,LUN) = INV(N,LUN)+INC                                 
            ENDDO                                                       
         ENDIF                                                          

         print *, 'LUN=', LUN
         print *, 'IOLUN(LUN)=', IOLUN(LUN)
         print *, 'IOMSG(LUN)=', IOMSG(LUN)

         print *, 'INC=', INC
                                                                        
C  CREATE NEW TABLE ENTRIES IF THIS UNIT DOESN'T SHARE EXISTING ONES    
C  -----------------------------------------------------------------    
                                                                        
         print *, 'NTBA(LUN)=', NTBA(LUN)
         print *, 'NTBB(LUN)=', NTBB(LUN)
         print *, 'NTBD(LUN)=', NTBD(LUN)
         NTBA(LUN)= 15
         NTBB(LUN)= 181
         NTBD(LUN)= 113
         print *, 'NTBA(LUN)=', NTBA(LUN)
         print *, 'NTBB(LUN)=', NTBB(LUN)
         print *, 'NTBD(LUN)=', NTBD(LUN)

         IF(LUS(LUN).EQ.0) THEN                                         
            CALL CHEKSTAB(LUN)                                          
            DO ITBA=1,NTBA(LUN)                                         
            INOD = NTAB+1                                               
            NEMO = TABA(ITBA,LUN)(4:11)                                 
            CALL TABSUB(LUN,NEMO)                                       
            MTAB(ITBA,LUN) = INOD                                       
            ISC(INOD)      = NTAB                                       
C           DO N1=INOD,ISC(INOD)-1                                      
C           DO N2=N1+1,ISC(INOD)                                        
C           IF(TAG(N1).EQ.TAG(N2)) GOTO 900                             
C           ENDDO                                                       
C           ENDDO                                                       
            ENDDO                                                       
         ENDIF                                                          
                                                                        
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  STORE TYPES AND INITIAL VALUES AND COUNTS                            
C  -----------------------------------------                            
                                                                        
      DO NODE=1,NTAB                                                    
      IF(TYP(NODE).EQ.'SUB') THEN                                       
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 0                                                 
      ELSEIF(TYP(NODE).EQ.'SEQ') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 0                                                 
      ELSEIF(TYP(NODE).EQ.'RPC') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 0                                                 
         ITP (NODE) = 0                                                 
      ELSEIF(TYP(NODE).EQ.'RPS') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 0                                                 
         ITP (NODE) = 0                                                 
      ELSEIF(TYP(NODE).EQ.'REP') THEN                                   
         VALI(NODE) = BMISS                                             
         KNTI(NODE) = IRF(NODE)                                         
         ITP (NODE) = 0                                                 
      ELSEIF(TYP(NODE).EQ.'DRS') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 1                                                 
      ELSEIF(TYP(NODE).EQ.'DRP') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 1                                                 
      ELSEIF(TYP(NODE).EQ.'DRB') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 0                                                 
         ITP (NODE) = 1                                                 
      ELSEIF(TYP(NODE).EQ.'NUM') THEN                                   
         VALI(NODE) = BMISS                                             
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 2                                                 
      ELSEIF(TYP(NODE).EQ.'CHR') THEN                                   
         VALI(NODE) = BMISS                                             
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 3                                                 
      ELSE                                                              
         GOTO 901                                                       
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  SET UP EXPANSION SEGMENTS FOR TYPE 'SUB', 'DRP', AND 'DRS' NODES     
C  ----------------------------------------------------------------     
                                                                        
      NEWN = 0                                                          
                                                                        
      DO N=1,NTAB                                                       
      ISEQ(N,1) = 0                                                     
      ISEQ(N,2) = 0                                                     
      EXPAND = TYP(N).EQ.'SUB' .OR. TYP(N).EQ.'DRP' .OR. TYP(N).EQ.'DRS'
     .                         .OR. TYP(N).EQ.'REP' .OR. TYP(N).EQ.'DRB'
      IF(EXPAND) THEN                                                   
         ISEQ(N,1) = NEWN+1                                             
         NODA = N                                                       
         NODE = N+1                                                     
         DO K=1,15000                                                   
         KNT(K) = 0                                                     
         ENDDO                                                          
         IF(TYP(NODA).EQ.'REP') KNT(NODE) = KNTI(NODA)                  
         IF(TYP(NODA).NE.'REP') KNT(NODE) = 1                           
                                                                        
1        NEWN = NEWN+1                                                  
         IF(NEWN.GT.15000) GOTO 902                                     
         JSEQ(NEWN) = NODE                                              
         KNT(NODE) = MAX(KNTI(NODE),KNT(NODE))                          
2        IF(JUMP(NODE)*KNT(NODE).GT.0) THEN                             
            NODE = JUMP(NODE)                                           
            GOTO 1                                                      
         ELSE IF(LINK(NODE).GT.0) THEN                                  
            NODE = LINK(NODE)                                           
            GOTO 1                                                      
         ELSE                                                           
            NODE = JMPB(NODE)                                           
            IF(NODE.EQ.NODA) GOTO 3                                     
            IF(NODE.EQ.0   ) GOTO 903                                   
            KNT(NODE) = MAX(KNT(NODE)-1,0)                              
            GOTO 2                                                      
         ENDIF                                                          
3        ISEQ(N,2) = NEWN                                               
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  PRINT THE SEQUENCE TABLES                                            
C  ------------------------                                             
                                                                        
      IF(PRTTAB) THEN                                                   
         PRINT*                                                         
         DO I=1,NTAB                                                    
         PRINT99,I,                                                     
     .   TAG(I),TYP(I),JMPB(I),JUMP(I),LINK(I),IBT(I),IRF(I),ISC(I)     
         ENDDO                                                          
         PRINT*                                                         
99       FORMAT(I5,2X,A10,A5,6I8)                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('MAKESTAB - DUP IN SUBSET: '//TAG(N1)//':'//NEMO)      
901   CALL BORT('MAKESTAB - UNKNOWN TYPE : '         //TYP(NODE))      
902   CALL BORT('MAKESTAB - JSEQ OVERFLOW       : '  //TAG(N   ))      
903   CALL BORT('MAKESTAB - FAILED TO CIRCULATE : '  //TAG(N   ))      
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MESGBF(LUNIT,MESGTYP)                                  
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*20 MESG                                                 
      CHARACTER*8  SEC0                                                 
      DIMENSION    MBAY(5000)                                           
      EQUIVALENCE  (MESG,MBAY(1))                                       
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      MESGTYP = -1                                                      
                                                                        
C  READ PAST ANY BUFR TABLES AND RETURN THE FIRST MESSAGE TYPE FOUND    
C  -----------------------------------------------------------------    
                                                                        
      CALL WRDLEN                                                       
      IMSG = 8/NBYTW+1                                                  
                                                                        
      REWIND LUNIT                                                      
1     READ(LUNIT,ERR=900,END=900) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))      
      MESGTYP = IUPM(MESG(17:17),8)                                      
      IF(MESGTYP.EQ.11) GOTO 1                                          
      REWIND LUNIT                                                      
                                                                        
900   RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MSGINI(LUN)                                            
                                                                        
      COMMON /PADESC/ IBCT,IPD1,IPD2,IPD3,IPD4                          
      COMMON /MSGPTR/ NBY0,NBY1,NBY2,NBY3,NBY4,NBY5                     
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*8  SUBTAG                                               
      CHARACTER*4  BUFR,SEVN                                            
      CHARACTER*3  TYP                                                  
      CHARACTER*1  TAB                                                  
                                                                        
      DATA BUFR/'BUFR'/                                                 
      DATA SEVN/'7777'/                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  GET THE MESSAGE TAG AND TYPE, AND BREAK UP THE DATE                  
C  ---------------------------------------------------                  
                                                                        
      SUBTAG = TAG(INODE(LUN))                                          
      CALL NEMTBA(LUN,SUBTAG,MTYP,MSBT,INOD)                            
      IF(INODE(LUN).NE.INOD) GOTO 900                                   
      CALL NEMTAB(LUN,SUBTAG,ISUB,TAB,IRET)                             
      IF(IRET.EQ.0) GOTO 901                                            
                                                                        
C  DATE CAN BE YYMMDDHH OR YYYYMMDDHH
C  ----------------------------------                                   
                                                                        
      MCEN = MOD(IDATE(LUN)/10**8,100)+1
      MEAR = MOD(IDATE(LUN)/10**6,100)                             
      MMON = MOD(IDATE(LUN)/10**4,100)                             
      MDAY = MOD(IDATE(LUN)/10**2,100)                             
      MOUR = MOD(IDATE(LUN)      ,100)                             
      MMIN = 0                                                       

      IF(MCEN.EQ.1) CALL BORT('MSGINI - IDATE(LUN) = 0!')
      IF(MEAR.EQ.0) MCEN = MCEN-1
      IF(MEAR.EQ.0) MEAR = 100
                                                                        
C  INITIALIZE THE MESSAGE                                               
C  ----------------------                                               
                                                                        
      MBIT = 0                                                          
      NBY0 = 8                                                          
      NBY1 = 18                                                         
      NBY2 = 0                                                          
      NBY3 = 20                                                         
      NBY4 = 4                                                          
      NBY5 = 4                                                          
      NBYT = NBY0+NBY1+NBY2+NBY3+NBY4+NBY5                              
                                                                        
C  SECTION 0                                                            
C  ---------                                                            
                                                                        
      CALL PKC(BUFR ,  4 , MBAY(1,LUN),MBIT)                            
      CALL PKB(NBYT , 24 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   3 ,  8 , MBAY(1,LUN),MBIT)                            
                                                                        
C  SECTION 1                                                            
C  ---------                                                            
                                                                        
      CALL PKB(NBY1 , 24 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   3 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   7 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MTYP ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MSBT ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   4 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MEAR ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MMON ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MDAY ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MOUR ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MMIN ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MCEN ,  8 , MBAY(1,LUN),MBIT)                            
                                                                        
C  SECTION 3                                                            
C  ---------                                                            
                                                                        
      CALL PKB(NBY3 , 24 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(2**7 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(IBCT , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(ISUB , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(IPD1 , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(IPD2 , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(IPD3 , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(IPD4 , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
                                                                        
C  SECTION 4                                                            
C  ---------                                                            
                                                                        
      CALL PKB(NBY4 , 24 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
                                                                        
C  SECTION 5                                                            
C  ---------                                                            
                                                                        
      CALL PKC(SEVN ,  4 , MBAY(1,LUN),MBIT)                            
                                                                        
C  DOUBLE CHECK INITIAL MESSAGE LENGTH                                  
C  -----------------------------------                                  
                                                                        
      IF(MOD(MBIT,8).NE.0) GOTO 902                                     
      IF(MBIT/8.NE.NBYT  ) GOTO 903                                     
                                                                        
      NMSG(LUN) = NMSG(LUN)+1                                           
      NSUB(LUN) = 0                                                     
      MBYT(LUN) = NBYT                                                  
                                                                        
      RETURN                                                            

C  ENTRY MINIMG WILL WRITE MINUTES INTO A MESSAGE SECTION ONE 
C  ----------------------------------------------------------

      ENTRY MINIMG(LUNIT,MINI)

      CALL STATUS(LUNIT,LUT,IL,IM)                                      
      IF(IL.LE.0) GOTO 904                                              
      IF(IM.EQ.0) GOTO 905                                              
      MBIT = 24*8
      CALL PKB(MINI,8,MBAY(1,LUT),MBIT)                            

      RETURN

C  ERROR EXITS
C  -----------

900   CALL BORT('MSGINI - INODE MISMATCH                      ')       
901   CALL BORT('MSGINI - CANT FIND SUBSET IDN: '//SUBTAG      )       
902   CALL BORT('MSGINI - INITIAL MESSAGE OUT OF BYTE BOUNDARY')       
903   CALL BORT('MSGINI - INITIAL MESSAGE COUNT FAILS CHECK   ')       
904   CALL BORT('MINIMG - FILE IS NOT OPEN FOR OUTPUT         ')       
905   CALL BORT('MINIMG - NO MESSGAGE OPEN                    ')       
      END                                                               
C---------------------------------------------------------------------- 
C  UPDATE THE MESSAGE BUFFER WITH NEW SUBSET                            
C---------------------------------------------------------------------- 
      SUBROUTINE MSGUPD(LUNIT,LUN)                                      
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /MSGPTR/ NBY0,NBY1,NBY2,NBY3,NBY4,NBY5                     
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
      CHARACTER*8 CBAY                                                  
      EQUIVALENCE (CBAY,JBAY)                                           
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  PAD THE SUBSET BUFFER                                                
C  ---------------------                                                
                                                                        
      CALL PAD(IBAY,IBIT,IBYT,8)                                        
      GOTO 1                                                            
                                                                        
C  SPECIAL ENTRY POINT FOR COPYSB                                       
C  ------------------------------                                       
                                                                        
      ENTRY SUBUPD(LUNIT,LUN,JBYT)                                      
      IBYT = JBYT                                                       
                                                                        
C  SEE IF THE NEW SUBSET FITS                                           
C  --------------------------                                           
                                                                        
1     IF(MBYT(LUN)+IBYT.GT.MAXBYT) THEN                                 
         CALL MSGWRT(LUNIT,MBAY(1,LUN),MBYT(LUN))                       
         CALL MSGINI(LUN)                                               
      ENDIF                                                             
                                                                        
      IF(IBYT.GT.MAXBYT-MBYT(LUN)) GOTO 900                             
                                                                        
C  SET A BYTE COUNT AND TRANSFER THE SUBSET BUFFER INTO THE MESSAGE     
C  ----------------------------------------------------------------     
                                                                        
      LBIT = 0                                                          
      CALL PKB(IBYT,16,IBAY,LBIT)                                       
      CALL MVB(IBAY,1,MBAY(1,LUN),MBYT(LUN)-3,IBYT)                     
                                                                        
C  UPDATE THE SUBSET AND BYTE COUNTERS                                  
C  --------------------------------------                               
                                                                        
      MBYT(LUN)   = MBYT(LUN)   + IBYT                                  
      NSUB(LUN)   = NSUB(LUN)   + 1                                     
                                                                        
      LBIT = (NBY0+NBY1+NBY2+4)*8                                       
      CALL PKB(NSUB(LUN),16,MBAY(1,LUN),LBIT)                           
                                                                        
      LBYT = NBY0+NBY1+NBY2+NBY3                                        
      NBYT = IUPB(MBAY(1,LUN),LBYT+1,24)                                
      LBIT = LBYT*8                                                     
      CALL PKB(NBYT+IBYT,24,MBAY(1,LUN),LBIT)                           
                                                                        
C  RESET THE USER ARRAYS AND EXIT NORMALLY                              
C  ---------------------------------------
                                                                        
      CALL USRTPL(LUN,1,1)                                              
      RETURN                                                            

C  ON ENCOUTERING OVERLARGE REPORTS RESET THE USER ARRAYS AND EXIT 'GRACEFULLY'
C  ----------------------------------------------------------------------------

900   PRINT*,'MSGUPD - SUBSET LONGER THAN ANY POSSIBLE MESSAGE   '
      PRINT*,'>>>>>>>OVERLARGE SUBSET DISCARDED FROM FILE<<<<<<<<'
      PRINT*,'MSGUPD - SUBSET LONGER THAN ANY POSSIBLE MESSAGE   '
      RETURN                                                            

      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MSGWRT(LUNIT,MBAY,MBYT)                                
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /BUFRMG/ MSGLEN,MSGTXT(5000)                               
                                                                        
      CHARACTER*4 BUFR,SEVN                                             
      DIMENSION   MBAY(*)                                               
                                                                        
      DATA BUFR/'BUFR'/                                                 
      DATA SEVN/'7777'/                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  MAKE SURE ALL SECTIONS HAVE EVEN NUMBER OF BYTES                     
C  ------------------------------------------------                     
                                                                        
      IAD1 = 8                                                          
      LEN1 = IUPB(MBAY,IAD1+1,24)                                       
      LEN2 = IUPB(MBAY,IAD1+8, 1)                                       
      MTYP = IUPB(MBAY,IAD1+9, 8)                                       
      IAD2 = IAD1+LEN1                                                  
      LEN2 = IUPB(MBAY,IAD2+1,24)*LEN2                                  
      IAD3 = IAD2+LEN2                                                  
      LEN3 = IUPB(MBAY,IAD3+1,24)                                       
      IAD4 = IAD3+LEN3                                                  
      LEN4 = IUPB(MBAY,IAD4+1,24)                                       
                                                                        
      IF(MOD(LEN1,2).NE.0) GOTO 901                                     
      IF(MOD(LEN2,2).NE.0) GOTO 902                                     
      IF(MOD(LEN3,2).NE.0) GOTO 903                                     
      IF(MOD(LEN4,2).NE.0) THEN                                         
         IAD5 = IAD4+LEN4                                               
         IBIT = IAD4*8                                                  
         LEN4 = LEN4+1                                                  
         CALL PKB(LEN4,24,MBAY,IBIT)                                    
         IBIT = IAD5*8                                                  
         CALL PKB(0,8,MBAY,IBIT)                                        
         MBYX = MBYT+1                                                  
      ELSE                                                              
         MBYX = MBYT                                                    
      ENDIF                                                             
                                                                        
C  WRITE SECTION 0 BYTE COUNT AND SECTION 5                             
C  ----------------------------------------                             
                                                                        
      IBIT = 0                                                          
      KBIT = (MBYX-4)*8                                                 
                                                                        
      CALL PKC(BUFR, 4,MBAY,IBIT)                                       
      CALL PKB(MBYX,24,MBAY,IBIT)                                       
      CALL PKC(SEVN, 4,MBAY,KBIT)                                       
                                                                        
C  ZERO OUT THE EXTRA BYTES WHICH WILL BE WRITTEN                       
C  ----------------------------------------------                       
                                                                        
      IMSG = 8/NBYTW                                                    
      MWRD = (MBYX/8+1)*IMSG                                            
      MBZZ = MWRD*NBYTW-MBYX
      DO I=1,MBZZ
      CALL PKB(0,8,MBAY,KBIT)                                       
      ENDDO
                                                                        
C  WRITE THE MESSAGE PLUS PADDING TO A WORD BOUNDARY                    
C  -------------------------------------------------                    
                                                                        
      IMSG = 8/NBYTW                                                    
      MWRD = (MBYX/8+1)*IMSG                                            
      WRITE(LUNIT) (MBAY(I),I=1,MWRD)                                   
C     PRINT*,'MSGWRT - LUNIT=',LUNIT,' BYTES=',MBYX                     
                                                                        
C  save a memory copy of this message - no bufr tables though           
C  ----------------------------------------------------------           
                                                                        
      IF(MTYP.NE.11) THEN                                               
         MSGLEN = MWRD                                                  
         DO I=1,MSGLEN                                                  
         MSGTXT(I) = MBAY(I)                                            
         ENDDO                                                          
      ENDIF                                                             
                                                                        
      RETURN                                                            
901   CALL BORT('MSGWRT - UNEVEN SECTION 1')                           
902   CALL BORT('MSGWRT - UNEVEN SECTION 2')                           
903   CALL BORT('MSGWRT - UNEVEN SECTION 3')                           
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MVB(IB1,NB1,IB2,NB2,NBM)                               
                                                                        
      DIMENSION IB1(*),IB2(*),NVAL(24000)
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND(UPB,PKB)                                                   
C-----------------------------------------------------------------------
                                                                        
      IF(NBM.GT.24000) CALL BORT('MVB - NBM>24000')
      JB1 = 8*(NB1-1)                                                   
      JB2 = 8*(NB2-1)                                                   
                                                                        
      DO N=1,NBM                                                        
      CALL UPB(NVAL(N),8,IB1,JB1)                                          
      ENDDO

      DO N=1,NBM
      CALL PKB(NVAL(N),8,IB2,JB2)                                          
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION NEMOCK(NEMO)                                             
                                                                        
      CHARACTER*(*) NEMO                                                
      CHARACTER*38  CHRSET                                              
                                                                        
      DATA CHRSET /'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.'/            
      DATA NCHR   /38/                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  GET THE LENGTH OF NEMO                                               
C  ----------------------                                               
                                                                        
      LNEMO = 0                                                         
                                                                        
      DO I=LEN(NEMO),1,-1                                               
      IF(NEMO(I:I).NE.' ') THEN                                         
         LNEMO = I                                                      
         GOTO 1                                                         
      ENDIF                                                             
      ENDDO                                                             
                                                                        
1     IF(LNEMO.LT.1 .OR. LNEMO.GT.8) THEN                               
         NEMOCK = -1                                                    
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  SCAN NEMO FOR ALLOWABLE CHARACTERS                                   
C  ----------------------------------                                   
                                                                        
      DO 10 I=1,LNEMO                                                   
      DO J=1,NCHR                                                       
      IF(NEMO(I:I).EQ.CHRSET(J:J)) GOTO 10                              
      ENDDO                                                             
      NEMOCK = -1                                                       
      RETURN                                                            
10    ENDDO                                                             
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      NEMOCK = 0                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NEMTAB(LUN,NEMO,IDN,TAB,IRET)

      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),
     .                TABA(50,32),TABB(250,32),TABD(250,32)

      CHARACTER*(*) NEMO
      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*8   NEMT
      CHARACTER*1   TAB
      LOGICAL       FOLVAL

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      FOLVAL = NEMO(1:1).EQ.'.'
      IRET = 0
      TAB = ' '

C  LOOK FOR NEMO IN TABLE B
C  ------------------------

      DO 1 I=1,NTBB(LUN)
      NEMT = TABB(I,LUN)(7:14)
      IF(NEMT.EQ.NEMO) THEN
         IDN  = IDNB(I,LUN)
         TAB  = 'B'
         IRET = I
         RETURN
      ELSEIF(FOLVAL.AND.NEMT(1:1).EQ.'.') THEN
         DO J=2,LEN(NEMT)
         IF(NEMT(J:J).NE.'.' .AND. NEMT(J:J).NE.NEMO(J:J)) GOTO 1
         ENDDO
         IDN  = IDNB(I,LUN)
         TAB  = 'B'
         IRET = I
         RETURN
      ENDIF
1     ENDDO

C  DON'T LOOK IN TABLE D FOR FOLLOWING VALUE-MNEMONICS
C  ---------------------------------------------------

      IF(FOLVAL) RETURN

C  LOOK IN TABLE D IF WE GOT THIS FAR
C  ----------------------------------

      DO I=1,NTBD(LUN)
      NEMT = TABD(I,LUN)(7:14)
      IF(NEMT.EQ.NEMO) THEN
         IDN  = IDND(I,LUN)
         TAB  = 'D'
         IRET = I
         RETURN
      ENDIF
      ENDDO

C  HERE CHECK FOR TABLE C OPERATOR DESCRIPTORS
C  -------------------------------------------

      IF(NEMO(1:3).EQ.'201' .OR.
     .   NEMO(1:3).EQ.'202' .OR.
     .   NEMO(1:3).EQ.'203' .OR.
     .   NEMO(1:3).EQ.'206' ) THEN
         READ(NEMO,'(I6)') IRET
         IDN = IFXY(NEMO)
         TAB = 'C'
         IRET = MOD(IRET/1000,10)
         RETURN
      ENDIF

      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NEMTBA(LUN,NEMO,MTYP,MSBT,INOD)                        
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*(*) NEMO                                                
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*20  NEMT                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NEMT = NEMO                                                       
      IRET = 0                                                          
                                                                        
C  LOOK FOR NEMO IN TABLE A                                             
C  ------------------------                                             
                                                                        
      DO I=1,NTBA(LUN)                                                  
      IF(TABA(I,LUN)(4:11).EQ.NEMO) THEN                                
         MTYP = IDNA(I,LUN,1)                                           
         MSBT = IDNA(I,LUN,2)                                           
         INOD = MTAB(I,LUN)                                             
         IF(MTYP.LT.0 .OR. MTYP.GT.255) GOTO 900                        
         IF(MSBT.LT.0 .OR. MSBT.GT.255) GOTO 901                        
         RETURN                                                         
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      CALL BORT('NEMTBA - CANT FIND '//NEMT)                           
900   CALL BORT('NEMTBA - BAD MTYP  '//NEMT)                           
901   CALL BORT('NEMTBA - BAD MSBT  '//NEMT)                           
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NEMTBAX(LUN,NEMO,MTYP,MSBT,INOD)
 
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),
     .                TABA(50,32),TABB(250,32),TABD(250,32)
 
      CHARACTER*(*) NEMO
      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*20  NEMT
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      NEMT = NEMO
      INOD = 0
 
C  LOOK FOR NEMO IN TABLE A
C  ------------------------
 
      DO I=1,NTBA(LUN)
      IF(TABA(I,LUN)(4:11).EQ.NEMO) THEN
         MTYP = IDNA(I,LUN,1)
         MSBT = IDNA(I,LUN,2)
         INOD = MTAB(I,LUN)
         IF(MTYP.LT.0 .OR. MTYP.GT.255) GOTO 900
         IF(MSBT.LT.0 .OR. MSBT.GT.255) GOTO 901
         RETURN
      ENDIF
      ENDDO
 
      RETURN
900   CALL BORT('NEMTBAX - BAD MTYP  '//NEMT)
901   CALL BORT('NEMTBAX - BAD MSBT  '//NEMT)
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NEMTBB(LUN,ITAB,UNIT,ISCL,IREF,IBIT)                   
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*24  UNIT                                                
      CHARACTER*8   NEMO                                                
      REAL*8        MXR                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      MXR = 1E11-1                                                      
                                                                        
      IF(ITAB.LE.0 .OR. ITAB.GT.NTBB(LUN)) GOTO 900                     
                                                                        
C  PULL OUT TABLE B INFORMATION                                         
C  ----------------------------                                         
                                                                        
      IDN  = IDNB(ITAB,LUN)                                             
      NEMO = TABB(ITAB,LUN)( 7:14)                                      
      UNIT = TABB(ITAB,LUN)(71:94)                                      
      ISCL = VALX(TABB(ITAB,LUN)( 95: 98))                              
      IREF = VALX(TABB(ITAB,LUN)( 99:109))                              
      IBIT = VALX(TABB(ITAB,LUN)(110:112))                              
                                                                        
C  CHECK TABLE B CONTENTS                                               
C  ----------------------                                               
                                                                        
      IF(IDN.LT.IFXY('000000')) GOTO 901                                
      IF(IDN.GT.IFXY('063255')) GOTO 901                                
                                                                        
      IF(ISCL.LT.-999 .OR. ISCL.GT.999) GOTO 902                        
      IF(IREF.LE.-MXR .OR. IREF.GE.MXR) GOTO 903                        
      IF(IBIT.LE.0) GOTO 904 
      IF(UNIT(1:5).NE.'CCITT' .AND. IBIT.GT.64      ) GOTO 904 
      IF(UNIT(1:5).EQ.'CCITT' .AND. MOD(IBIT,8).NE.0) GOTO 905
                                                                        
      RETURN                                                            
900   CALL BORT('NEMTBB - ITAB NOT IN TABLE B'         )               
901   CALL BORT('NEMTBB - BAD DESCRIPTOR VALUE: '//NEMO)               
902   CALL BORT('NEMTBB - BAD SCALE VALUE     : '//NEMO)               
903   CALL BORT('NEMTBB - BAD REFERENCE VALUE : '//NEMO)               
904   CALL BORT('NEMTBB - BAD BIT WIDTH       : '//NEMO)               
905   CALL BORT('NEMTBB - BAD CHAR BIT WIDTH  : '//NEMO)               
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NEMTBD(LUN,ITAB,NSEQ,NEMS,IRPS,KNTS)

      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),
     .                TABA(50,32),TABB(250,32),TABD(250,32)

      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*8   NEMO,NEMS,NEMT,NEMF
      CHARACTER*1   TAB
      DIMENSION     NEMS(250),IRPS(250),KNTS(250)
      LOGICAL       REP

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      IF(ITAB.LE.0 .OR. ITAB.GT.NTBD(LUN)) GOTO 900

      REP  = .FALSE.

C  CLEAR THE RETURN VALUES
C  -----------------------

      NSEQ = 0

      DO I=1,250
      NEMS(I) = ' '
      IRPS(I) = 0
      KNTS(I) = 0
      ENDDO

C  PARSE THE TABLE D ENTRY
C  -----------------------

      NEMO = TABD(ITAB,LUN)(7:14)
      IDSC = IDND(ITAB,LUN)
      CALL UPTDD(ITAB,LUN,0,NDSC)

      IF(IDSC.LT.IFXY('300000')) GOTO 901
      IF(IDSC.GT.IFXY('363255')) GOTO 901
C     IF(NDSC.LE.0             ) GOTO 902

      DO J=1,NDSC
      IF(NSEQ+1.GT.250) GOTO 903
      CALL UPTDD(ITAB,LUN,J,IDSC)
      CALL NUMTAB(LUN,IDSC,NEMT,TAB,IRET)
      IF(TAB.EQ.'R') THEN
         IF(REP) GOTO 904
         REP = .TRUE.
         IF(IRET.LT.0) THEN
            IRPS(NSEQ+1) = 1
            KNTS(NSEQ+1) = ABS(IRET)
         ELSEIF(IRET.GT.0) THEN
            IRPS(NSEQ+1) = IRET
         ENDIF
      ELSEIF(TAB.EQ.'F') THEN
         IF(.NOT.REP) GOTO 904
         IRPS(NSEQ+1) = IRET
         REP = .FALSE.
      ELSEIF(TAB.EQ.'D'.OR.TAB.EQ.'C') THEN
         REP = .FALSE.
         NSEQ = NSEQ+1
         NEMS(NSEQ) = NEMT
      ELSEIF(TAB.EQ.'B') THEN
         REP = .FALSE.
         NSEQ = NSEQ+1
         IF(NEMT(1:1).EQ.'.') THEN
            CALL UPTDD(ITAB,LUN,J+1,IDSC)
            CALL NUMTAB(LUN,IDSC,NEMF,TAB,IRET)
            CALL RSVFVM(NEMT,NEMF)
            IF(TAB.NE.'B') GOTO 906
         ENDIF
         NEMS(NSEQ) = NEMT
      ELSE
         GOTO 905
      ENDIF
      ENDDO

      RETURN
900   CALL BORT('NEMTBD - ITAB NOT IN TABLE D   '                )
901   CALL BORT('NEMTBD - BAD DESCRIPTOR VALUE: '          //NEMO)
902   CALL BORT('NEMTBD - ZERO LENGTH SEQUENCE: '          //NEMO)
903   CALL BORT('NEMTBD - TOO MANY DESCRIPTORS IN SEQ: '   //NEMO)
904   CALL BORT('NEMTBD - REPLICATOR OUT OF ORDER IN SEQ: '//NEMO)
905   CALL BORT('NEMTBD - BAD DESCRIPTOR IN SEQUENCE: '    //NEMO)
906   CALL BORT('NEMTBD - FOLLOWING VALUE NOT FROM TABLEB:'//NEMF)
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NENUCK(NEMO,NUMB,LUN)                                  
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*8   NEMO                                                
      CHARACTER*6   NUMB                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK TABLE A                                                        
C  -------------                                                        
                                                                        
      ENTRY NENUAA(NEMO,NUMB,LUN)                                       
                                                                        
      DO N=1,NTBA(LUN)                                                  
      IF(NUMB(4:6).EQ.TABA(N,LUN)(1: 3)) GOTO 900                       
      IF(NEMO     .EQ.TABA(N,LUN)(4:11)) GOTO 900                       
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  CHECK TABLE B AND D                                                  
C  -------------------                                                  
                                                                        
      ENTRY NENUBD(NEMO,NUMB,LUN)                                       
                                                                        
      DO N=1,NTBB(LUN)                                                  
      IF(NUMB.EQ.TABB(N,LUN)(1: 6)) GOTO 900                            
      IF(NEMO.EQ.TABB(N,LUN)(7:14)) GOTO 900                            
      ENDDO                                                             
                                                                        
      DO N=1,NTBD(LUN)                                                  
      IF(NUMB.EQ.TABD(N,LUN)(1: 6)) GOTO 900                            
      IF(NEMO.EQ.TABD(N,LUN)(7:14)) GOTO 900                            
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXIT                                                           
C  ----------                                                           
                                                                        
900   CALL BORT('NENUCK - DUPLICATE NEM/NUM '//NEMO//' '//NUMB)        
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE NEWWIN(LUN,IWIN,JWIN)                                  
                                                                        
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
cfpp$ expand (lstrpc)                                                   
C---------------------------------------------------------------------- 
                                                                        
      IF(IWIN.EQ.1) THEN                                                
         JWIN = NVAL(LUN)                                               
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  REFIND THE JWIN BOUNDARY FROM IWIN                                   
C  ----------------------------------                                   
                                                                        
      NODE = INV(IWIN,LUN)                                              
      IF(LSTRPC(NODE,LUN).NE.NODE) CALL BORT('NEWWIN - NOT RPC')       
      JWIN = IWIN+VAL(IWIN,LUN)                                         
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION NMBYT(LUNIT)
 
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      NMBYT = -1
 
C  CHECK THE FILE STATUS
C  ---------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      IF(IM.EQ.0) GOTO 902
      NMBYT = IUPB(MBAY(1,LUN),5,24)


      RETURN
 
C  ERROR EXITS
C  -----------
 
900   CALL BORT('NMBYT - FILE IS CLOSED                   ')
901   CALL BORT('NMBYT - FILE IS OPEN FOR OUTPUT          ')
902   CALL BORT('NMBYT - NO MESSAGE IS OPENED             ')
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION NMSUB(LUNIT)                                             
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NMSUB = 0                                                         
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      IF(IM.EQ.0) GOTO 902                                              
                                                                        
      NMSUB = MSUB(LUN)                                                 
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('NMSUB - FILE IS CLOSED                   ')           
901   CALL BORT('NMSUB - FILE IS OPEN FOR OUTPUT          ')           
902   CALL BORT('NMSUB - NO MESSAGE IS OPENED             ')           
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION NUMBCK(NUMB)                                             
                                                                        
      CHARACTER*10 CHRSET                                               
      CHARACTER*6  NUMB                                                 
      CHARACTER*1  FC                                                   
      LOGICAL      DIGIT
                                                                        
      DATA CHRSET /'0123456789'/                                        
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NUMBCK = 0                                                        
      LNUMB  = 0                                                        
      FC     = NUMB(1:1)                                                
                                                                        
C  CHECK THE FIRST CHARACTER OF NUMB                                    
C  ---------------------------------                                    
                                                                        
      IF(.NOT.(FC.EQ.'A' .OR. FC.EQ.'0' .OR. FC.EQ.'3')) GOTO 900       
                                                                        
C  CHECK THE REST OF NUMB                                               
C  ----------------------                                               
                                                                        
      DO 10 I=2,6                                                       
      DO J=1,10                                                         
      IF(NUMB(I:I).EQ.CHRSET(J:J)) GOTO 10                              
      ENDDO                                                             
      GOTO 900                                                          
10    ENDDO                                                             
                                                                        
C  CHECK FOR A VALID DESCRIPTOR                                         
C  ----------------------------                                         
                                                                        
      IF(DIGIT(NUMB(2:6))) THEN
         READ(NUMB,'(1X,I2,I3)') IX,IY                             
      ELSE
         GOTO 900
      ENDIF
      IF(IX.LT.0 .OR. IX.GT. 63) GOTO 900                               
      IF(IY.LT.0 .OR. IY.GT.255) GOTO 900                               
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      NUMBCK = 0                                                        
      RETURN                                                            
                                                                        
C  ERROR EXIT                                                           
C  ----------                                                           
                                                                        
900   NUMBCK = -1                                                       
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NUMTAB(LUN,IDN,NEMO,TAB,IRET)

      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),
     .                TABA(50,32),TABB(250,32),TABD(250,32)

      CHARACTER*(*) NEMO
      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*6   ADN30,CID
      CHARACTER*3   TYPS
      CHARACTER*1   REPS,TAB

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      NEMO = ' '
      IRET = 0
      TAB = ' '

C  LOOK FOR A REPLICATOR OR A REPLICATOR FACTOR
C  --------------------------------------------

      IF(IDN.GE.IDNR(1,1) .AND. IDN.LE.IDNR(1,2)) THEN
         TAB  = 'R'
         IRET = -MOD(IDN,256)
         RETURN
      ENDIF

      DO I=2,5
      IF(IDN.EQ.IDNR(I,1)) THEN
         TAB  = 'R'
         IRET = I
         RETURN
      ELSEIF(IDN.EQ.IDNR(I,2)) THEN
         TAB  = 'F'
         IRET = I
         RETURN
      ENDIF
      ENDDO

C  LOOK FOR IDN IN TABLE D
C  -----------------------

      DO I=1,NTBD(LUN)
      IF(IDN.EQ.IDND(I,LUN)) THEN
         NEMO = TABD(I,LUN)(7:14)
         TAB  = 'D'
         IRET = I
         RETURN
      ENDIF
      ENDDO

C  LOOK FOR IDN IN TABLE B
C  -----------------------

      DO I=1,NTBB(LUN)
      IF(IDN.EQ.IDNB(I,LUN)) THEN
         NEMO = TABB(I,LUN)(7:14)
         TAB  = 'B'
         IRET = I
         RETURN
      ENDIF
      ENDDO

C  LOOK FOR IDN IN TABLE C
C  -----------------------

      CID = ADN30(IDN,6)
      CID = CID(1:3)
      IF(CID.EQ.'201' .OR. CID.EQ.'202' .OR. CID.EQ.'206') THEN
         NEMO = ADN30(IDN,6)
         TAB  = 'C'
         IRET = MOD(IDN,10)
         RETURN
      ENDIF

      RETURN
      END
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      FUNCTION NVNWIN(NODE,LUN,INV1,INV2,INVN,NMAX)                     
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      DIMENSION    INVN(NMAX)                                           
      REAL*8       VAL,BMISS                                            
                                                                        
      DATA BMISS/10E10/                                                 
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(NODE.EQ.0) RETURN                                              
      NVNWIN = 0                                                        
                                                                        
      DO I=1,NMAX                                                       
      INVN(I) = BMISS                                                   
      ENDDO                                                             
                                                                        
C  SEARCH BETWEEN INV1 AND INV2                                         
C  ----------------------------                                         
                                                                        
      DO N=INV1,INV2                                                    
      IF(INV(N,LUN).EQ.NODE) THEN                                       
         NVNWIN = NVNWIN+1                                              
         INVN(NVNWIN) = N                                               
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      IF(NVNWIN.GT.NMAX) CALL BORT('NVNWIN - TOO MANY EVENTS')         
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION NWORDS(N,LUN)                                            
                                                                        
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      REAL*8 VAL                                                        
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NWORDS = 0                                                        
                                                                        
      DO K=1,NINT(VAL(N,LUN))                                           
      NWORDS = NWORDS + NINT(VAL(NWORDS+N+1,LUN))                       
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE NXTWIN(LUN,IWIN,JWIN)                                  
                                                                        
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
cfpp$ expand (lstrpc)                                                   
C---------------------------------------------------------------------- 
                                                                        
      IF(JWIN.EQ.NVAL(LUN)) THEN                                        
         IWIN = 0                                                       
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  FIND THE NEXT SEQUENTIAL WINDOW                                      
C  -------------------------------                                      
                                                                        
      NODE = INV(IWIN,LUN)                                              
      IF(LSTRPC(NODE,LUN).NE.NODE) print*,'bad node=',node,iwin         
      IF(LSTRPC(NODE,LUN).NE.NODE) CALL BORT('NXTWIN - NOT RPC')       
      IF(VAL(JWIN,LUN).EQ.0) THEN                                       
         IWIN = 0                                                       
      ELSE                                                              
         IWIN = JWIN                                                    
         JWIN = IWIN+VAL(IWIN,LUN)                                      
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE OPENBF(LUNIT,IO,LUNDX)                                 
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /STBFR / IOLUN(32),IOMSG(32)                               
      COMMON /QUIET / IPRT                                              
                                                                        
      CHARACTER*(*) IO                                                  
      CHARACTER*4   BUFR,MSTR                                           
      LOGICAL       SKIPDX,APPEND                                       
                                                                        
      DATA IFIRST/0/                                                    
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      write(unit=*, fmt='(a, i4)') 'LUNIT=', LUNIT
      write(unit=*, fmt='(a, i4)') 'LUNDX=', LUNDX

      write(unit=*, fmt='(2a)') ' IO =', trim(IO)
                                                                        
      IF(IFIRST.EQ.0) THEN                                              
         CALL WRDLEN                                                    
         CALL BFRINI                                                    
         IFIRST = 1                                                     
      ENDIF                                                             
                                                                        
      IF(IO.EQ.'QUIET') THEN                                            
         IPRT = LUNDX                                                   
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  SEE IF A FILE CAN BE OPENED                                          
C  ---------------------------                                          
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      

      write(unit=*, fmt='(a, i4)') 'LUN=', LUN
      write(unit=*, fmt='(a, i4)') ' IL=',  IL
      write(unit=*, fmt='(a, i4)') ' IM=',  IM

      IF(LUN.EQ.0) GOTO 900                                             
      IF(IL .NE.0) GOTO 901                                             
                                                                        
C  CHECK FOR NO BUFR DATA OR NO DATA AT ALL IN AN "IN" FILE             
C  --------------------------------------------------------             
                                                                        
      IF(IO.EQ.'IN' .AND. LUNIT.EQ.LUNDX) THEN                          
         REWIND LUNIT                                                   
         READ(LUNIT,END=100,ERR=902) MSTR                               

         write(unit=*, fmt='(2a)') ' MSTR =', trim(MSTR)

         IBIT = 0                                                       
         CALL UPC(BUFR,4,MSTR,IBIT)                                     
         IF(BUFR.NE.'BUFR') GOTO 902                                    
      ENDIF                                                             
                                                                        
C  SET INITIAL OPEN DEFAULTS                                            
C  -------------------------                                            
                                                                        
      REWIND LUNIT                                                      
      NMSG (LUN) = 0                                                    
      NSUB (LUN) = 0                                                    
      MSUB (LUN) = 0                                                    
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
      SKIPDX = .FALSE.                                                  
      APPEND = .FALSE.                                                  
                                                                        
C  DECIDE HOW TO SETUP THE DICTIONARY                                   
C  ----------------------------------                                   
                                                                        
      IF(IO.EQ.'IN') THEN                                               
         CALL WTSTAT(LUNIT,LUN,-1,0)                                    
         CALL READDX(LUNIT,LUN,LUNDX)                                   
      ELSE IF(IO.EQ.'OUT') THEN                                         
         CALL WTSTAT(LUNIT,LUN, 1,0)                                    
         CALL WRITDX(LUNIT,LUN,LUNDX)                                   
      ELSE IF(IO.EQ.'APN' .OR. IO.EQ.'APX') THEN                        
         CALL WTSTAT(LUNIT,LUN, 1,0)                                    
         CALL READDX(LUNIT,LUN,LUNDX)                                   
         IF(IO.EQ.'APN') CALL POSAPN(LUNIT)                             
         IF(IO.EQ.'APX') CALL POSAPX(LUNIT)                             
      ELSE                                                              
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
                                                                        
C  FILE OPENED FOR INPUT IS EMPTY - LET READMG GIVE THE BAD NEWS        
C  -------------------------------------------------------------        
                                                                        
100   REWIND LUNIT                                                      
      CALL WTSTAT(LUNIT,LUN,-1,0)                                       
      CALL DXINIT(LUN,0)                                                
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('OPENBF - TOO MANY FILES OPENED ALREADY       ')       
901   CALL BORT('OPENBF - FILE ALREADY OPEN                   ')       
902   CALL BORT('OPENBF - INPUT FILE HAS NON-BUFR DATA        ')       
903   CALL BORT('OPENBF - IO MUST BE ONE OF "IN" "OUT" "APN"  ')       
      END                                                               
C-----------------------------------------------------------------------
C  DEFAULT OPENBT RETURNS LUNDX=0 TO INDICATE NO LINKS TO BUFRTABLES
C-----------------------------------------------------------------------
      SUBROUTINE OPENBT(LUNDX,MTYP)
      LUNDX = 0
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE OPENMB(LUNIT,SUBSET,JDATE)                             
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
                                                                        
      CHARACTER*(*) SUBSET                                              
      LOGICAL       OPEN                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
                                                                        
C  GET SOME SUBSET PARTICULARS                                          
C  ---------------------------                                          
                                                                         
      CALL NEMTBA(LUN,SUBSET,MTYP,MSTB,INOD)                            
      OPEN = IM.EQ.0.OR.INOD.NE.INODE(LUN).OR.I4DY(JDATE).NE.IDATE(LUN)   
                                                                        
C  MAYBE OPEN A NEW OR DIFFERENT TYPE OF MESSAGE                        
C  ---------------------------------------------                        
                                                                        
      IF(OPEN) THEN                                                     
         CALL CLOSMG(LUNIT)                                             
         CALL WTSTAT(LUNIT,LUN,IL, 1)                                   
         INODE(LUN) = INOD                                              
         IDATE(LUN) = I4DY(JDATE)
         CALL MSGINI(LUN)                                               
         CALL USRTPL(LUN,1,1)                                           
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('OPENMB - FILE IS CLOSED            ')                 
901   CALL BORT('OPENMB - FILE IS OPEN FOR INPUT    ')                 
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE OPENMG(LUNIT,SUBSET,JDATE)                             
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
                                                                        
      CHARACTER*(*) SUBSET                                              
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.NE.0) CALL CLOSMG(LUNIT)                                    
      CALL WTSTAT(LUNIT,LUN,IL, 1)                                      
                                                                        
C  GET SOME SUBSET PARTICULARS                                          
C  ---------------------------                                          
                                                                        
      CALL NEMTBA(LUN,SUBSET,MTYP,MSTB,INOD)                            
      INODE(LUN) = INOD                                                 
      IDATE(LUN) = I4DY(JDATE)
                                                                        
C  INITIALIZE THE OPEN MESSAGE                                          
C  ---------------------------                                          
                                                                        
      CALL MSGINI(LUN)                                                  
      CALL USRTPL(LUN,1,1)                                              
                                                                        
      RETURN                                                            
900   CALL BORT('OPENMG - FILE IS CLOSED            ')                 
901   CALL BORT('OPENMG - FILE IS OPEN FOR INPUT    ')                 
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE PAD(IBAY,IBIT,IBYT,IPADB)                              
                                                                        
      DIMENSION IBAY(*)                                                 
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  PAD THE SUBSET TO AN IPADB BIT BOUNDARY                              
C  ----------------------------------------                             
                                                                        
      IPAD = IPADB - MOD(IBIT+8,IPADB)                                  
      CALL PKB(IPAD,8,IBAY,IBIT)                                        
      CALL PKB(   0,IPAD,IBAY,IBIT)                                     
      IBYT = IBIT/8                                                     
                                                                        
      IF(MOD(IBIT,IPADB).NE.0) GOTO 900                                 
      IF(MOD(IBIT,8    ).NE.0) GOTO 900                                 
                                                                        
      RETURN                                                            
900   CALL BORT('PAD - BIT PAD FAILURE              ')                 
      END                                                               
C-----------------------------------------------------------------------
C  PARSE SEPARATE WORDS FROM A STRING SEQUENCE                          
C-----------------------------------------------------------------------
      SUBROUTINE PARSEQ(STR,TAGS,MTAG,NTAG)                             
                                                                        
      CHARACTER*(*) STR,TAGS(MTAG)                                      
      CHARACTER*80  ASTR                                                
      LOGICAL       WORD                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      ASTR = STR                                                        
      LSTR = LEN(STR)                                                   
      LTAG = LEN(TAGS(1))                                               
      IF(LSTR.GT.80) GOTO 900                                           
      NTAG = 0                                                          
      NCHR = 0                                                          
      WORD = .FALSE.                                                    
                                                                        
      DO 10 I=1,LSTR                                                    
                                                                        
      IF(.NOT.WORD .AND. STR(I:I).NE.' ') THEN                          
         NTAG = NTAG+1                                                  
         IF(NTAG.GT.MTAG) GOTO 901                                      
         TAGS(NTAG) = ' '                                               
      ENDIF                                                             
                                                                        
      IF(WORD .AND. STR(I:I).EQ.' ') NCHR = 0                           
      WORD = STR(I:I).NE.' '                                            
                                                                        
      IF(WORD) THEN                                                     
         NCHR = NCHR+1                                                  
         IF(NCHR.GT.LTAG) GOTO 902                                      
         TAGS(NTAG)(NCHR:NCHR) = STR(I:I)                               
      ENDIF                                                             
                                                                        
10    CONTINUE                                                          
                                                                        
      RETURN                                                            
900   CALL BORT('PARSEQ - STRING TOO LONG  '//ASTR)                    
901   CALL BORT('PARSEQ - TOO MANY TAGS    '//ASTR)                    
902   CALL BORT('PARSEQ - TAG IS TOO LONG  '//ASTR)                    
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE PARUSR(STR,LUN,I1,IO)                                  
                                                                        
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
      common /ACMODE/ iac
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*80  UST                                                 
      CHARACTER*20  UTG(30)                                             
      LOGICAL       BUMP                                                
                                                                        
      DATA MAXUSR /30/                                                  
      DATA MAXNOD /20/                                                  
      DATA MAXCON /10/                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      UST  = STR                                                        
      IF(LEN(STR).GT.80) GOTO 900                                       
                                                                        
      NCON = 0                                                          
      NNOD = 0                                                          
                                                                        
C  PROCESS STRING PIECES(S) INTO COND NODES AND STORE NODES             
C  --------------------------------------------------------             
                                                                        
      CALL PARSEQ(UST,UTG,MAXUSR,NTOT)                                  
                                                                        
      DO N=1,NTOT                                                       
      CALL PARUTG(LUN,IO,UTG(N),NOD,KON,VAL,*908)                       
      IF(KON.NE.0) THEN                                                 
         NCON = NCON+1                                                  
         IF(NCON.GT.MAXCON) GOTO 901                                    
         NODC(NCON) = NOD                                               
         KONS(NCON) = KON                                               
         IVLS(NCON) = NINT(VAL)                                              
      ELSE                                                              
         NNOD = NNOD+1                                                  
         IF(NNOD.GT.MAXNOD) GOTO 902                                    
         NODS(NNOD) = NOD                                               
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  SORT COND NODES IN JUMP/LINK TABLE ORDER                             
C  ----------------------------------------                             
                                                                        
      DO I=1,NCON                                                       
      DO J=I+1,NCON                                                     
      IF(NODC(I).GT.NODC(J)) THEN                                       
         NOD     = NODC(I)                                              
         NODC(I) = NODC(J)                                              
         NODC(J) = NOD                                                  
                                                                        
         KON     = KONS(I)                                              
         KONS(I) = KONS(J)                                              
         KONS(J) = KON                                                  
                                                                        
         VAL     = IVLS(I)                                              
         IVLS(I) = IVLS(J)                                              
         IVLS(J) = VAL                                                  
      ENDIF                                                             
      ENDDO                                                             
      ENDDO                                                             
                                                                        
C  CHECK ON SPECIAL RULES FOR BUMP NODES                                
C  -------------------------------------                                
                                                                        
      BUMP = .FALSE.                                                    
                                                                        
      DO N=1,NCON                                                       
      IF(KONS(N).EQ.5) THEN                                             
         IF(IO.EQ.0)   GOTO 903                                         
         IF(N.NE.NCON) GOTO 904                                         
         BUMP = .TRUE.                                                  
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  CHECK STORE NODE COUNT AND ALIGNMENT                                 
C  ------------------------------------                                 
                                                                        
      IF(.NOT.BUMP .AND. NNOD.EQ.0) GOTO 905                            
      IF(NNOD.GT.I1)                GOTO 906                            
                                                                        
      IRPC = -1                                                         
      DO I=1,NNOD                                                       
      IF(NODS(I).GT.0) THEN                                             
         IF(IRPC.LT.0) IRPC = LSTRPC(NODS(I),LUN)                       
         IF(IRPC.NE.LSTRPC(NODS(I),LUN).and.iac.eq.0) GOTO 907  
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('PARUSR - USER STRING > 80 CHARS         :'//UST)      
901   CALL BORT('PARUSR - TOO MANY COND NODES            :'//UST)      
902   CALL BORT('PARUSR - TOO MANY STOR NODES            :'//UST)      
903   CALL BORT('PARUSR - BUMP ON INPUT NOT ALLOWED      :'//UST)      
904   CALL BORT('PARUSR - BUMP MUST BE ON INNER NODE     :'//UST)      
905   CALL BORT('PARUSR - USER STRING HAS NO STORE NODES :'//UST)      
906   CALL BORT('PARUSR - MUST BE AT LEAST I1 STORE NODES:'//UST)      
907   CALL BORT('PARUSR - STORE NODES MUST IN ONE REP GRP:'//UST)      
908   CALL BORT('PARUSR - PARUTG:'                         //UST)      
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE PARUTG(LUN,IO,UTG,NOD,KON,VAL,*)                       
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /UTGPRM/ PICKEY                                            
                                                                        
      CHARACTER*20  UTG,ATAG                                            
      CHARACTER*10  TAG                                                 
      CHARACTER*3   TYP,ATYP,BTYP                                       
      CHARACTER*1   COND(5)                                             
      DIMENSION     BTYP(8),IOK(8)                                      
      LOGICAL       PICKEY                                              
                                                                        
      DATA NCHK   / 8/                                                  
      DATA BTYP   /'SUB','SEQ','REP','RPC','RPS','DRB','DRP','DRS'/     
      DATA IOK    /  -1 ,  -1 ,  -1 ,  -1 ,  -1 ,  -1 ,   0 ,   0 /     
      DATA LTG    /20/                                                  
                                                                        
C---------------------------------------------------------------------- 
      PICKEY = .FALSE.                                                  
      COND(1) = '='                                                     
      COND(2) = '!'                                                     
      COND(3) = '<'                                                     
      COND(4) = '>'                                                     
      COND(5) = '^'                                                     
      NCOND   = 5                                                       
C---------------------------------------------------------------------- 
                                                                        
      ATAG  = ' '                                                       
      ATYP  = ' '                                                       
      KON   = 0                                                         
      NOD   = 0                                                         
      VAL   = 0                                                         
                                                                        
C  PARSE THE TAG                                                        
C  -------------                                                        
                                                                        
      DO I=1,LTG                                                        
      IF(UTG(I:I).EQ.' ') GOTO 1                                        
      DO J=1,NCOND                                                      
      IF(UTG(I:I).EQ.COND(J)) THEN                                      
         KON = J                                                        
         ICV = I+1                                                      
         GOTO 1                                                         
      ENDIF                                                             
      ENDDO                                                             
      ATAG(I:I) = UTG(I:I)                                              
      ENDDO                                                             
                                                                        
C  FIND THE TAG IN THE SUBSET TABLE                                     
C  --------------------------------                                     
                                                                        
1     INOD = INODE(LUN)                                                 
      DO NOD=INOD,ISC(INOD)                                             
      IF(ATAG.EQ.TAG(NOD)) GOTO 2                                       
      ENDDO                                                             
                                                                        
      IF(KON.EQ.0 .AND. (IO.EQ.0.OR.ATAG.EQ.'NUL'.OR..NOT.PICKEY)) THEN 
C     IF(KON.EQ.0) THEN                                                 
         NOD = 0                                                        
         RETURN                                                         
      ELSE                                                              
         PRINT*,'TRYING TO WRITE A NON-EXISTANT MNEMONIC:'//ATAG  
         RETURN 1                                                       
      ENDIF                                                             
                                                                        
C  CHECK FOR A VALID NODE TYP                                           
C  --------------------------                                           
                                                                        
2     IF(KON.EQ.5) THEN                                                 
         IF(TYP(NOD-1).NE.'DRP' .AND. TYP(NOD-1).NE.'DRS') GOTO 901     
      ELSE                                                              
         ATYP = TYP(NOD)                                                
         DO I=1,NCHK                                                    
         IF(ATYP.EQ.BTYP(I) .AND. IO.NE.IOK(I)) GOTO 902                
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  IF A COND NODE GET THE COND VALUE                                    
C  ---------------------------------                                    
                                                                        
      IF(KON.NE.0) THEN                                                 
         CALL STRNUM(UTG(ICV:LTG),NUM)                                  
         IF(NUM.LT.0) GOTO 903                                          
         VAL = NUM                                                      
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('PARUTG - NO VALID TAG FOUND IN              :'//UTG)  
901   CALL BORT('PARUTG - BUMP NODE MUST BE TYPE RPC(DRP)    :'//UTG)  
902   CALL BORT('PARUTG - ILLEGAL NODE TYPE:'//ATYP//       ':'//UTG)  
903   CALL BORT('PARUTG - BAD OR MISSING COND VALUE IN       :'//UTG)  
      END                                                               
C----------------------------------------------------------------------
C  PACK UP A NUMBER ACCORDING TO SPECS
C----------------------------------------------------------------------
      SUBROUTINE PKB(NVAL,NBITS,IBAY,IBIT)
 
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
         INT = ISHFT(IVAL,2*NBITW-(NBT+NBITS))
         MSK = ISHFT(  -1,2*NBITW-(NBT+NBITS))
         IBAY(NWD+1) = IREV(IOR(IAND(IREV(IBAY(NWD+1)),NOT(MSK)),INT))
      ENDIF
 
      IBIT = IBIT + NBITS
 
      RETURN
      END
C----------------------------------------------------------------------
C  COPY CHARACTERS INTO A BIT ARRAY
C----------------------------------------------------------------------
      SUBROUTINE PKC(CHR,NCHR,IBAY,IBIT)
 
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
      IVAL = 0
      NBIT = 8
 
      DO I=1,NCHR
      CVAL(LB) = CHR(I:I)
      IF(IASCII.EQ.0) CALL IPKM(CVAL(LB),1,IETOA(IUPM(CVAL(LB),8)))
      NWD  = IBIT/NBITW + 1
      NBT  = MOD(IBIT,NBITW)
      INT = ISHFT(IVAL,NBITW-NBIT)
      INT = ISHFT(INT,-NBT)
      MSK = ISHFT(  -1,NBITW-NBIT)
      MSK = ISHFT(MSK,-NBT)
      IBAY(NWD) = IREV(IOR(IAND(IREV(IBAY(NWD)),NOT(MSK)),INT))
      IF(NBT+NBIT.GT.NBITW) THEN
         INT = ISHFT(IVAL,2*NBITW-(NBT+NBIT))
         MSK = ISHFT(  -1,2*NBITW-(NBT+NBIT))
         IBAY(NWD+1) = IREV(IOR(IAND(IREV(IBAY(NWD+1)),NOT(MSK)),INT))
      ENDIF
      IBIT = IBIT + NBIT
      ENDDO
 
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PKTDD(ID,LUN,IDN,IRET)                                 
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*56  DXSTR                                               
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      LDD = LDXD(IDXV+1)+1                                              
                                                                        
C  ZERO THE COUNTER IF IDN IS ZERO                                      
C  -------------------------------                                      
                                                                        
      IF(IDN.EQ.0) THEN                                                 
         CALL IPKM(TABD(ID,LUN)(LDD:LDD),1,0)                           
         IRET = 0                                                       
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  UPDATE THE STORED DESCRIPTOR COUNT FOR THIS TABLE D ENTRY            
C  ---------------------------------------------------------            
                                                                        
      ND = IUPM(TABD(ID,LUN)(LDD:LDD),8)                                
                                                                        
      IF(ND.LT.0 .OR. ND.EQ.250) THEN                                   
         IRET = -1                                                      
         RETURN                                                         
      ELSE                                                              
         ND = ND+1                                                      
         CALL IPKM(TABD(ID,LUN)(LDD:LDD),1,ND)                          
         IRET = ND                                                      
      ENDIF                                                             
                                                                        
C  PACK AND STORE THE DESCRIPTOR                                        
C  -----------------------------                                        
                                                                        
      IDM = LDD+1 + (ND-1)*2                                            
      CALL IPKM(TABD(ID,LUN)(IDM:IDM),2,IDN)                            
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE POSAPN(LUNIT)                                          
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 SEC0                                                  
      DIMENSION   MBAY(5000)                                            
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  READ AND COUNT MESSAGES IN ORDER TO POSTION FOR APPEND               
C  ------------------------------------------------------               
                                                                        
      IMSG = 8/NBYTW+1                                                  
      REWIND LUNIT                                                      
      IREC = 0                                                          
                                                                        
1     READ(LUNIT,ERR=2,END=2) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))          
      IREC = IREC+1                                                     
      GOTO 1                                                            
                                                                        
2     REWIND LUNIT                                                      
      DO J=1,IREC                                                       
      READ(LUNIT,ERR=900,END=901) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))      
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('POSAPN - IO ERR READING A MESSAGE')                   
901   CALL BORT('POSAPN - FAILURE TO READ TO EOFLE')                   
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE POSAPX(LUNIT)                                          
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 SEC0                                                  
      DIMENSION   MBAY(5000)                                            
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  TRY TO READ TO THE END OF THE FILE AND BACKSPACE FOR APPEND          
C  ----------------------------------------------------                 
                                                                        
      REWIND LUNIT                                                      
      IMSG = 8/NBYTW+1                                                  
      IREC = 0                                                          
                                                                        
1     READ(LUNIT,END=2,ERR=3) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))          
      IREC = IREC+1                                                     
      GOTO 1                                                            
                                                                        
C  IF SUCCESSFUL BACKSPACE FOR APPENDING AND RETURN                     
C  ------------------------------------------------                     
                                                                        
2     BACKSPACE LUNIT                                                   
      RETURN                                                            
                                                                        
C  IF AN I/O ERROR IS ENCOUNTERED REREAD THE GOOD RECORDS AND RETURN    
C  -----------------------------------------------------------------    
                                                                        
3     REWIND LUNIT                                                      
      DO J=1,IREC                                                       
      READ(LUNIT,END=2,ERR=900) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))        
      ENDDO                                                             
      RETURN                                                            
                                                                        
C  IF ALL ELSE FAILS JUST GIVE UP AT THIS POINT                         
C  --------------------------------------------                         
                                                                        
900   CALL BORT('POSAPX - IO ERR REREADING A MESSAGE')                 
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE RCSTPL(LUN)                                            
                                                                        
      PARAMETER (MAXINV=15000)                                          
      PARAMETER (MAXRCR=100 )                                           
                                                                        
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRBIT/ NBIT(15000),MBIT(15000)                           
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      DIMENSION ITMP(MAXINV,MAXRCR),VTMP(MAXINV,MAXRCR)                 
      DIMENSION NBMP(2,MAXRCR),NEWN(2,MAXRCR)                           
      DIMENSION KNX(MAXRCR)                                             
      REAL*8    VAL,VTMP                                                
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND (UPBB)                                                     
C-----------------------------------------------------------------------
                                                                        
C  SET THE INITIAL VALUES FOR THE TEMPLATE                              
C  ---------------------------------------                              
                                                                        
      INV(1,LUN) = INODE(LUN)                                           
      VAL(1,LUN) = 0                                                    
      NBMP(1,1) = 1                                                     
      NBMP(2,1) = 1                                                     
      NODI = INODE(LUN)                                                 
      NODE = INODE(LUN)                                                 
      MBMP = 1                                                          
      KNVN = 1                                                          
      NR   = 0                                                          
                                                                        
      DO I=1,MAXRCR                                                     
      KNX(I) = 0                                                        
      ENDDO                                                             
                                                                        
C  SET UP THE PARAMETRES FOR A LEVEL OF RECURSION                       
C  ----------------------------------------------                       
                                                                        
10    CONTINUE                                                          
                                                                        
      NR = NR+1                                                         
      NBMP(1,NR) = 1                                                    
      NBMP(2,NR) = MBMP                                                 
                                                                        
      N1 = ISEQ(NODE,1)                                                 
      N2 = ISEQ(NODE,2)                                                 
      IF(N1.EQ.0          ) GOTO 905                                    
      IF(N2-N1+1.GT.MAXINV) GOTO 906                                    
      NEWN(1,NR) = 1                                                    
      NEWN(2,NR) = N2-N1+1                                              
                                                                        
      DO N=1,NEWN(2,NR)                                                 
      NN = JSEQ(N+N1-1)                                                 
      ITMP(N,NR) = NN                                                   
      VTMP(N,NR) = VALI(NN)                                             
      if(vtmp(n,nr).gt.10e9) vtmp(n,nr) = 10e10                         
      ENDDO                                                             
                                                                        
C  STORE NODES AT SOME RECURSION LEVEL                                  
C  -----------------------------------                                  
                                                                        
20    DO I=NBMP(1,NR),NBMP(2,NR)                                        
      IF(KNX(NR).EQ.0000) KNX(NR) = KNVN                                
      IF(I.GT.NBMP(1,NR)) NEWN(1,NR) = 1                                
      DO J=NEWN(1,NR),NEWN(2,NR)                                        
      KNVN = KNVN+1                                                     
      NODE = ITMP(J,NR)                                                 
      INV(KNVN,LUN) = NODE                                              
      VAL(KNVN,LUN) = VTMP(J,NR)                                        
      MBIT(KNVN) = MBIT(KNVN-1)+NBIT(KNVN-1)                            
      NBIT(KNVN) = IBT(NODE)                                            
      IF(ITP(NODE).EQ.1) THEN                                           
         CALL UPBB(MBMP,NBIT(KNVN),MBIT(KNVN),MBAY(1,LUN))              
         NEWN(1,NR) = J+1                                               
         NBMP(1,NR) = I                                                 
         GOTO 10                                                        
      ENDIF                                                             
      ENDDO                                                             
      NEW = KNVN-KNX(NR)                                                
      VAL(KNX(NR)+1,LUN) = VAL(KNX(NR)+1,LUN) + NEW                     
      KNX(NR) = 0                                                       
      ENDDO                                                             
                                                                        
C  CONTINUE AT ONE RECUSION LEVEL BACK                                  
C  -----------------------------------                                  
                                                                        
      IF(NR-1.NE.0) THEN                                                
         NR = NR-1                                                      
         GOTO 20                                                        
      ENDIF                                                             
                                                                        
C  FINALLY STORE THE LENGTH OF THE SUBSET TEMPLATE                      
C  -----------------------------------------------                      
                                                                        
      NVAL(LUN) = KNVN                                                  
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('RCSTPL - NBMP <> 1 FOR        : '//TAG(NODI))         
901   CALL BORT('RCSTPL - NODE NOT SUB,DRP,DRS : '//TAG(NODI))         
902   CALL BORT('RCSTPL - NEGATIVE REP FACTOR  : '//TAG(NODI))         
903   CALL BORT('RCSTPL - REP FACTOR OVERFLOW  : '//TAG(NODI))         
904   CALL BORT('RCSTPL - INVENTORY INDEX OUT OF BOUNDS     ')         
905   CALL BORT('RCSTPL - UNSET EXPANSION SEG  : '//TAG(NODI))         
906   CALL BORT('RCSTPL - TEMP ARRAY OVERFLOW  : '//TAG(NODI))         
907   CALL BORT('RCSTPL - INVENTORY OVERFLOW   : '//TAG(NODI))         
908   CALL BORT('RCSTPL - TPL CACHE OVERFLOW   : '//TAG(NODI))         
909   CALL BORT('RCSTPL - BAD BACKUP STRATEGY  : '//TAG(NODI))         
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RDBFDX(LUNIT,LUN)                                      
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB,TABB1,TABB2                                    
      CHARACTER*128 TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*50  DXCMP                                               
      CHARACTER*8   SEC0,NEMO
      CHARACTER*6   NUMB,CIDN                                           
      CHARACTER*1   MOCT(24000)                                         
      DIMENSION     MBAY(5000),LDXBD(10),LDXBE(10)                      
      EQUIVALENCE   (MBAY(1),MOCT(1))                                   
      LOGICAL       DIGIT
                                                                        
      DATA LDXBD /38,70,8*0/                                            
      DATA LDXBE /42,42,8*0/                                            
                                                                        
C-----------------------------------------------------------------------
      JA(I) = IA+1+LDA*(I-1)                                            
      JB(I) = IB+1+LDB*(I-1)                                            
C-----------------------------------------------------------------------
                                                                        
C  INITIALIZE THE DX-TABLE PARTITION AND SOURCE FILE                    
C  -------------------------------------------------                    
                                                                        
      CALL DXINIT(LUN,0)                                                
      REWIND LUNIT                                                      
      IDX = 0                                                           
                                                                        
C  CLEAR THE BUFFER AND READ A MESSAGE                                  
C  -----------------------------------                                  
                                                                        
1     DO I=1,5000                                                       
      MBAY(I) = 0                                                       
      ENDDO                                                             
                                                                        
      IMSG = 8/NBYTW+1                                                  
      READ(LUNIT,ERR=900,END=2) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))        
      IDX = IDX+1                                                       
                                                                        
C  GET THE SECTION START OCTETS AND LENGTHS                             
C  ----------------------------------------                             
                                                                        
2     I1 = 8                                                            
      L1 = IUPM(MOCT(I1+1),24)                                          
      I2 = I1+L1                                                        
      L2 = IUPM(MOCT(I2+1),24)*IUPM(MOCT(I1+8),1)                       
      I3 = I2+L2                                                        
      L3 = IUPM(MOCT(I3+1),24)                                          
      I4 = I3+L3                                                        
      L4 = IUPM(MOCT(I4+1),24)                                          
                                                                        
C  SEE IF THIS IS A BUFR DX MESSAGE - CHECK FOR RECOGNISABLE DX VERSION 
C  -------------------------------------------------------------------- 
                                                                        
      IF(IUPM(MOCT(I1+9),8).NE.11) THEN                                 
         REWIND LUNIT                                                   
         DO NDX=1,IDX-1                                                 
         READ(LUNIT,ERR=910,END=910) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))   
         ENDDO                                                          
         CALL MAKESTAB                                                  
         RETURN                                                         
      ENDIF                                                             
                                                                        
      IDXS = IUPM(MOCT(I1+10),8)+1
      IF(IDXS.GT.IDXV+1) IDXS = IUPM(MOCT(I1+12),8)+1
      IF(LDXA(IDXS).EQ.0) GOTO 902                                      
      IF(LDXB(IDXS).EQ.0) GOTO 902                                      
      IF(LDXD(IDXS).EQ.0) GOTO 902                                      
      L30 = LD30(IDXS)                                                  
                                                                        
      DXCMP = ' '                                                       
      CALL CHRTRN(DXCMP,MOCT(I3+8),NXSTR(IDXS))                         
      IF(DXCMP.NE.DXSTR(IDXS)) GOTO 902                                 
                                                                        
C  SECTION 4 - READ DEFINITIONS FOR TABLES A B AND D                    
C  -------------------------------------------------                    
                                                                        
      LDA  = LDXA (IDXS)                                                
      LDB  = LDXB (IDXS)                                                
      LDD  = LDXD (IDXS)                                                
      LDBD = LDXBD(IDXS)                                                
      LDBE = LDXBE(IDXS)                                                
                                                                        
      IA = I4+5                                                         
      LA = IUPM(MOCT(IA),8)                                             
      IB = JA(LA+1)                                                     
      LB = IUPM(MOCT(IB),8)                                             
      ID = JB(LB+1)                                                     
      LD = IUPM(MOCT(ID),8)                                             
                                                                        
C  TABLE A - MESSAGE TYPE/SUBTYPE FROM THE NEMONIC OR THE SEQ NUMBER    
C  -----------------------------------------------------------------    
                                                                        
      DO I=1,LA                                                         
      N = NTBA(LUN)+1                                                   
      IF(N.GT.NTBA(0)) GOTO 903                                         
      CALL CHRTRNA(TABA(N,LUN),MOCT(JA(I)),LDA)                         
      NUMB = '   '//TABA(N,LUN)(1:3)                                           
      NEMO = TABA(N,LUN)(4:11)                                          
      CALL NENUAA(NEMO,NUMB,LUN)                                        
      NTBA(LUN) = N                                                     
                                                                        
      IF(DIGIT(NEMO(3:8))) THEN
         READ(NEMO,'(2X,2I3)') MTYP,MSBT                            
         IDNA(N,LUN,1) = MTYP                                              
         IDNA(N,LUN,2) = MSBT                                              
      ELSE
         READ(NUMB(4:6),'(I3)') IDNA(N,LUN,1)                           
         IDNA(N,LUN,2) = 0                                              
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  TABLE B                                                              
C  -------                                                              
                                                                        
      DO I=1,LB                                                         
      N = NTBB(LUN)+1                                                   
      IF(N.GT.NTBB(0)) GOTO 904                                         
      CALL CHRTRNA(TABB1,MOCT(JB(I)     ),LDBD)                         
      CALL CHRTRNA(TABB2,MOCT(JB(I)+LDBD),LDBE)                         
      TABB(N,LUN) = TABB1(1:LDXBD(IDXV+1))//TABB2(1:LDXBE(IDXV+1))      
      NUMB = TABB(N,LUN)(1:6)                                           
      NEMO = TABB(N,LUN)(7:14)                                          
      CALL NENUBD(NEMO,NUMB,LUN)                                        
      IDNB(N,LUN) = IFXY(NUMB)                                          
      NTBB(LUN) = N                                                     
      ENDDO                                                             
                                                                        
C  TABLE D                                                              
C  -------                                                              
                                                                        
      DO I=1,LD                                                         
      N = NTBD(LUN)+1                                                   
      IF(N.GT.NTBD(0)) GOTO 905                                         
      CALL CHRTRNA(TABD(N,LUN),MOCT(ID+1),LDD)                          
      NUMB = TABD(N,LUN)(1:6)                                           
      NEMO = TABD(N,LUN)(7:14)                                          
      CALL NENUBD(NEMO,NUMB,LUN)                                        
      IDND(N,LUN) = IFXY(NUMB)                                          
      ND = IUPM(MOCT(ID+LDD+1),8)                                       
      IF(ND.GT.250) GOTO 906                                            
      DO J=1,ND                                                         
      NDD = ID+LDD+2 + (J-1)*L30                                        
      CALL CHRTRNA(CIDN,MOCT(NDD),L30)                                  
      IDN = IDN30(CIDN,L30)                                             
      CALL PKTDD(N,LUN,IDN,IRET)                                        
      IF(IRET.LT.0) GOTO 908                                            
      ENDDO                                                             
      ID = ID+LDD+1 + ND*L30                                            
      IF(IUPM(MOCT(ID+1),8).EQ.0) ID = ID+1                             
      NTBD(LUN) = N                                                     
      ENDDO                                                             
                                                                        
C  GOTO READ THE NEXT MESSAGE                                           
C  --------------------------                                           
                                                                        
      GOTO 1                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('RDBFDX - I/O ERROR READING DX MESSAGE          ')     
901   CALL BORT('RDBFDX - EOF >>>>> READING DX MESSAGE          ')     
902   CALL BORT('RDBFDX - UNEXPECTED DX MESSAGE TYPE OR CONTENTS')     
903   CALL BORT('RDBFDX - TOO MANY TABLE A ENTRIES              ')     
904   CALL BORT('RDBFDX - TOO MANY TABLE B ENTRIES              ')     
905   CALL BORT('RDBFDX - TOO MANY TABLE D ENTRIES              ')     
906   CALL BORT('RDBFDX - TOO MANY DESCRIPTORS IN TABLE D ENTRY ')     
907   CALL BORT('RDBFDX - ERROR READING IDN SEQ FROM MOCT       ')     
908   CALL BORT('RDBFDX - BAD RETURN FROM PKTDD                 ')     
909   CALL BORT('RDBFDX - DESC COUNT IN TABD <> MOCT            ')     
910   CALL BORT('RDBFDX - ERR/EOF POSITIONING AFTER DX MESSAGES ')     
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RDCMPS(LUN)
 
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
      LPS(LBIT) = MAX(2**(LBIT)-1,1)
      UPS(NODE) = FLOAT(IVAL+IRF(NODE))*10.**(-ISC(NODE))
C-----------------------------------------------------------------------
 
C  SETUP THE SUBSET TEMPLATE
C  -------------------------
 
      CALL USRTPL(LUN,1,1)
 
C  UNCOMPRESS A SUBSET INTO THE VAL ARRAY ACCORDING TO TABLE B
C  -----------------------------------------------------------
 
      NELM = NVAL(LUN)
      NSBS = NSUB(LUN)
      IBIT = MBYT(LUN)
 
      DO N=1,NELM
      NODE = INV(N,LUN)
      NBIT = IBT(NODE)
      ITYP = ITP(NODE)
      IF(ITYP.EQ.2) THEN
         CALL UPB(LREF,NBIT,MBAY(1,LUN),IBIT)
         CALL UPB(LINC,   6,MBAY(1,LUN),IBIT)
         JBIT = IBIT + LINC*(NSBS-1)
         CALL UPB(NINC,LINC,MBAY(1,LUN),JBIT)
         IF(NINC.EQ.LPS(LINC)) NINC = LPS(NBIT)
         IVAL = LREF+NINC
         IF(IVAL.LT.LPS(NBIT)) VAL(N,LUN) = UPS(NODE)
         IBIT = IBIT + LINC*MSUB(LUN)
      ELSEIF(ITYP.EQ.3) THEN
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
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RDMEMM(IMSG,SUBSET,JDATE,IRET)                         
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
                                                                        
      CHARACTER*8 SUBSET                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE MESSAGE REQUEST AND FILE STATUS                            
C  -----------------------------------------                            
                                                                        
      CALL STATUS(MUNIT,LUN,IL,IM)                                      
      CALL WTSTAT(MUNIT,LUN,IL, 1)                                      
      IF(IL.GE.0) GOTO 900                                              
      IRET = 0                                                          
                                                                        
      IF(IMSG.EQ.0 .OR.IMSG.GT.MSGP(0)) THEN                            
         CALL WTSTAT(MUNIT,LUN,IL,0)                                    
         IRET = -1                                                      
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  READ MESSAGE# IMSG INTO A MESSAGE BUFFER                             
C  ----------------------------------------                             
                                                                        
      IPTR = MSGP(IMSG)                                                 
      IF(IMSG.LT.MSGP(0)) LPTR = MSGP(IMSG+1)-IPTR                      
      IF(IMSG.EQ.MSGP(0)) LPTR = MLAST-IPTR+1                           
      IPTR = IPTR-1                                                     
                                                                        
      DO I=1,LPTR                                                       
      MBAY(I,LUN) = MSGS(IPTR+I)                                        
      ENDDO                                                             
                                                                        
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,JRET)
      NMSG(LUN) = IMSG
      RETURN
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('RDMEMM - FILE IS NOT OPEN FOR INPUT       ')          
901   CALL BORT('RDMEMM - BAD IPTR IN MEMORY BUFFER        ')          
902   CALL BORT('RDMEMM - BAD NPTR IN MEMORY BUFFER        ')          
903   CALL BORT('RDMEMM - MSGTYPE MISMATCH FOR '//SUBSET    )          
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RDMEMS(ISUB,IRET)                                      
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /UNPTYP/ MSGUNP(32)
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      CALL STATUS(MUNIT,LUN,IL,IM)                                      
      IF(IL       .GE.0) GOTO 900                                       
      IF(IM       .EQ.0) GOTO 901                                       
      IF(NSUB(LUN).NE.0) GOTO 902                                       
                                                                        
      IF(ISUB.GT.MSUB(LUN)) IRET = -1                                   
      IF(ISUB.GT.MSUB(LUN)) RETURN                                      
                                                                        
      MBYM = MBYT(LUN)                                                  
      NBYT = 0                                                          
                                                                        
C  POSITION TO READ SUBSET #ISUB FROM MEMORY MESSAGE #IMSG            
C  -------------------------------------------------------                          
      IF(MSGUNP(LUN).EQ.0) THEN
         NSUB(LUN) = ISUB-1                                                
         DO I=1,ISUB-1                                                   
         MBYT(LUN) = MBYT(LUN) + IUPB(MBAY(1,LUN),MBYT(LUN)+1,16)
         ENDDO                                                             
      ELSEIF(MSGUNP(LUN).EQ.1) THEN
         DO I=1,ISUB-1
         CALL READSB(MUNIT,IRET)
         ENDDO
      ELSEIF(MSGUNP(LUN).EQ.2) THEN
         NSUB(LUN) = ISUB-1                                                
      ENDIF

C  READ SUBSET #ISUB FROM MEMORY MESSAGE #IMSG            
C  -------------------------------------------                          

      CALL READSB(MUNIT,IRET)                                           
      IF(IRET.NE.0) GOTO 903                                            

C  RESET SUBSET POINTER BACK TO ZERO AND RETURN
C  --------------------------------------------

      MBYT(LUN) = MBYM                                                  
      NSUB(LUN) = 0                                                     
      RETURN                                                            

C  ERROR EXITS
C  -----------
                                                                        
900   CALL BORT('RDMEMS - MEMORY FILE NOT OPEN FOR INPUT')             
901   CALL BORT('RDMEMS - MEMORY MESSAGE NOT OPEN       ')             
902   CALL BORT('RDMEMS - MEMORY MESSAGE NOT AT BEGINING')             
903   CALL BORT('RDMEMS - ERROR CALLING READSB          ')             
      END                                                               
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE RDTREE(LUN)
 
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),
     .                JUMP(15000),LINK(15000),JMPB(15000),
     .                IBT(15000),IRF(15000),ISC(15000),
     .                ITP(15000),VALI(15000),KNTI(15000),
     .                ISEQ(15000,2),JSEQ(15000)
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)
      COMMON /USRBIT/ NBIT(15000),MBIT(15000)
 
      CHARACTER*10 TAG
      CHARACTER*8  CVAL
      CHARACTER*3  TYP
      DIMENSION    IVAL(15000)
      EQUIVALENCE  (CVAL,RVAL)
      REAL*8       VAL,RVAL
 
C-----------------------------------------------------------------------
CFPP$ EXPAND (UPBB)
C-----------------------------------------------------------------------
      MPS(NODE) = 2**(IBT(NODE))-1
      UPS(NODE) = float(IVAL(N)+IRF(NODE))*10.**(-ISC(NODE))
C-----------------------------------------------------------------------
 
C  CYCLE THROUGH A SUBSET SETTING UP THE USER ARRAY
C  ------------------------------------------------
 
      MBIT(1) = IBIT
      NBIT(1) = 0
      CALL RCSTPL(LUN)
 
C  UNPACK A SUBSET INTO THE USER ARRAY
C  -----------------------------------
 
      DO N=1,NVAL(LUN)
      CALL UPBB(IVAL(N),NBIT(N),MBIT(N),MBAY(1,LUN))
      ENDDO
 
C  CONVERT THE UNPACKED INTEGERS TO THE PROPER TYPES
C  -------------------------------------------------
 
      DO N=1,NVAL(LUN)
      NODE = INV(N,LUN)
      IF(ITP(NODE).EQ.1) THEN
         VAL(N,LUN) = IVAL(N)
      ELSEIF(ITP(NODE).EQ.2) THEN
         IF(IVAL(N).LT.MPS(NODE)) VAL(N,LUN) = UPS(NODE)
      ENDIF
      ENDDO
 
C  SPECIAL TREATMENT FOR CHARACTERS
C  --------------------------------
 
      DO N=1,NVAL(LUN)
      NODE = INV(N,LUN)
      IF(ITP(NODE).EQ.3) THEN
         CVAL = ' '
         CALL UPC(CVAL,NBIT(N)/8,MBAY(1,LUN),MBIT(N))
         VAL(N,LUN) = RVAL
      ENDIF
      ENDDO
 
      IBIT = NBIT(NVAL(LUN))+MBIT(NVAL(LUN))
 
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RDUSDX(LUNDX,LUN)                                      
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*80  CARD                                                
      CHARACTER*8   NEMO                                  
      CHARACTER*6   NUMB                                                
      LOGICAL       DIGIT
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  INITIALIZE THE DX-TABLE PARTITION AND SOURCE FILE                    
C  -------------------------------------------------                    
                                                                        
      CALL DXINIT(LUN,1)                                                
      REWIND LUNDX                                                      
                                                                        
C  READ USER CARDS UNTIL THERE ARE NO MORE                              
C  ---------------------------------------                              
                                                                        
1     READ(LUNDX,'(A80)',END=100) CARD                                  
                                                                        
C  REREAD IF NOT A DEFINITION CARD                                      
C  -------------------------------                                      
                                                                        
      IF(CARD(1: 1).EQ.       '*') GOTO 1                               
      IF(CARD(3:10).EQ.'--------') GOTO 1                               
      IF(CARD(3:10).EQ.'        ') GOTO 1                               
      IF(CARD(3:10).EQ.'MNEMONIC') GOTO 1                               
      IF(CARD(3:10).EQ.'TABLE  D') GOTO 1                               
      IF(CARD(3:10).EQ.'TABLE  B') GOTO 1                               
                                                                        
C  PARSE A DESCRIPTOR DEFINITION CARD                                   
C  ----------------------------------                                   
                                                                        
      IF(CARD(12:12).EQ.'|' .AND. CARD(21:21).EQ.'|') THEN              
                                                                        
         NEMO = CARD(3:10)                                              
         NUMB = CARD(14:19)                                             
         IF(NEMOCK(NEMO).NE.0) GOTO 900                                 
         IF(NUMBCK(NUMB).NE.0) GOTO 900                                 
                                                                        
         IF(NUMB(1:1).EQ.'A') THEN                                      
            N = NTBA(LUN)+1                                             
            IF(N.GT.NTBA(0)) GOTO 901                                   
            CALL NENUAA(NEMO,NUMB,LUN)                             
            TABA(N,LUN)( 1: 3) = NUMB(4:6)                              
            TABA(N,LUN)( 4:11) = NEMO                                   
            TABA(N,LUN)(13:67) = CARD(23:77)                            
            NTBA(LUN) = N                                               
                                                                        
            IF(DIGIT(NEMO(3:8))) THEN
               READ(NEMO,'(2X,2I3)') MTYP,MSBT                            
               IDNA(N,LUN,1) = MTYP                                     
               IDNA(N,LUN,2) = MSBT              
            ELSE
               READ(NUMB(4:6),'(I3)') IDNA(N,LUN,1)                           
               IDNA(N,LUN,2) = 0                                              
            ENDIF                                                             
                                                                        
            NUMB(1:1) = '3'                                             
         ENDIF                                                          
                                                                        
         IF(NUMB(1:1).EQ.'0') THEN                                      
            N = NTBB(LUN)+1                                             
            IF(N.GT.NTBB(0)) GOTO 902                                   
            CALL NENUBD(NEMO,NUMB,LUN)                                  
            IDNB(N,LUN) = IFXY(NUMB)                                    
            TABB(N,LUN)( 1: 6) = NUMB                                   
            TABB(N,LUN)( 7:14) = NEMO                                   
            TABB(N,LUN)(16:70) = CARD(23:77)                            
            NTBB(LUN) = N                                               
            GOTO 1                                                      
         ENDIF                                                          
                                                                        
         IF(NUMB(1:1).EQ.'3') THEN                                      
            N = NTBD(LUN)+1                                             
            IF(N.GT.NTBD(0)) GOTO 903                                   
            CALL NENUBD(NEMO,NUMB,LUN)                                  
            IDND(N,LUN) = IFXY(NUMB)                                    
            TABD(N,LUN)( 1: 6) = NUMB                                   
            TABD(N,LUN)( 7:14) = NEMO                                   
            TABD(N,LUN)(16:70) = CARD(23:77)                            
            NTBD(LUN) = N                                               
            GOTO 1                                                      
         ENDIF                                                          
                                                                        
         GOTO 904                                                       
      ENDIF                                                             
                                                                        
C  PARSE A SEQUENCE DEFINITION CARD                                     
C  --------------------------------                                     
                                                                        
      IF(CARD(12:12).EQ.'|' .AND. CARD(19:19).NE.'|') THEN              
         CALL SEQSDX(CARD,LUN)                                          
         GOTO 1                                                         
      ENDIF                                                             
                                                                        
C  PARSE AN ELEMENT DEFINITION CARD                                     
C  --------------------------------                                     
                                                                        
      IF(CARD(12:12).EQ.'|' .AND. CARD(19:19).EQ.'|') THEN              
         CALL ELEMDX(CARD,LUN)                                          
         GOTO 1                                                         
      ENDIF                                                             
                                                                        
C  CAN'T FIGURE OUT WHAT KIND OF CARD IT IS                             
C  ----------------------------------------                             
                                                                        
      GOTO 905                                                          
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
100   CALL MAKESTAB                                                     
      RETURN                                                            
                                                                        
C  ERROR EXIT                                                           
C  ----------                                                           
                                                                        
900   PRINT*,CARD                                                       
      CALL BORT('RDUSDX - NEMO OR NUMB ERROR             '//CARD)      
901   CALL BORT('RDUSDX - TOO MANY TABLE A ENTRIES       '//CARD)      
902   CALL BORT('RDUSDX - TOO MANY TABLE B ENTRIES       '//CARD)      
903   CALL BORT('RDUSDX - TOO MANY TABLE D ENTRIES       '//CARD)      
904   CALL BORT('RDUSDX - BAD DESCRIPTOR NUMBER          '//CARD)      
905   CALL BORT('RDUSDX - BAD CARD FORMAT                '//CARD)      
      END                                                               
C-----------------------------------------------------------------------
C  BUFR TABLE INFORMATION CONTAINED IN A FILE CONNECTED TO UNIT LUNDX   
C  IS USED TO INITIALIZE PROCESSING TABLES FOR A BUFR FILE CONNECTED    
C  TO UNIT LUNIT. LUNIT AND LUNDX MAY BE THE SAME ONLY IF THE UNIT IS   
C  CONNECTED TO A BUFR FILE, CURRENTLY OPEN FOR INPUT PROCESSING,       
C  POSITIONED AT A DX-TABLE MESSAGE (ANYWHERE IN THE FILE). OTHERWISE,  
C  LUNDX MAY BE CONNECTED TO ANOTHER CURRENTLY OPEN AND DEFINED BUFR    
C  FILE, OR TO A USER SUPPLIED, CHARACTER FORMAT, DX-TABLE FILE.        
C                                                                       
C  NOTE: READDX IS USED TO INITIALIZE INTERNAL BUFR DX-TABLES ONLY.     
C        IF A BUFR OUTPUT FILE IS BEING OPENED, SUBROUTINE WRITDX       
C        CALLS READDX TO INITIALIZE THE INTERNAL DX-TABLES, AND THEN    
C        WRITES BUFR DX-TABLE MESSAGES INTO THE OUTPUT FILE.            
C                                                                       
C                                                                       
C  INPUT ARGUMENTS:                                                     
C     LUNIT    - UNIT CONNECTED TO BUFR FILE TO BE INITIALIZED/UPDATED  
C     LUN      - INTERNAL BUFR UNIT ASSOCIATED WITH FORTRAN UNIT LUNIT  
C     LUNDX    - UNIT CONTAINING DX-TABLES                              
C                                                                       
C-----------------------------------------------------------------------
      SUBROUTINE READDX(LUNIT,LUN,LUNDX)                                
                                                                        
      COMMON /QUIET/ IPRT                                               
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  GET THE BUFR STATUS OF UNIT LUNDX                                    
C  ---------------------------------                                    
                                                                        
      CALL STATUS(LUNDX,LUD,ILDX,IMDX)                                  
                                                                        
C  READ A DX-TABLE FROM THE INDICATED SOURCE                            
C  -----------------------------------------                            
                                                                        
      IF (LUNIT.EQ.LUNDX) THEN                                          
         IF(IPRT.GE.1) PRINT100,LUNDX,LUNIT                             
         REWIND LUNIT                                                   
         CALL RDBFDX(LUNIT,LUN)                                         
      ELSEIF(ILDX.NE.0) THEN                                            
         IF(IPRT.GE.1) PRINT101,LUNDX,LUNIT                             
         CALL CPBFDX(LUD,LUN)                                           
      ELSEIF(ILDX.EQ.0) THEN                                            
         IF(IPRT.GE.1) PRINT102,LUNDX,LUNIT                             
         REWIND LUNDX                                                   
         CALL RDUSDX(LUNDX,LUN)                                         
      ELSE                                                              
         CALL BORT('READDX - SCREWUP')                                 
      ENDIF                                                             
                                                                        
100   FORMAT(' READING BUFR DX-TABLES FROM ',I2,' TO ',I2)              
101   FORMAT(' COPYING BUFR DX-TABLES FROM ',I2,' TO ',I2)              
102   FORMAT(' READING USER DX-TABLES FROM ',I2,' TO ',I2)              
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READERM(LUNIT,SUBSET,JDATE,IRET)                       
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
      CHARACTER*8 SUBSET                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      CALL WTSTAT(LUNIT,LUN,IL, 1)                                      
                                                                        
C  READ A MESSAGE INTO A MESSAGE BUFFER                                 
C  ------------------------------------                                 
                                                                        
1     IF(IRDERM(LUNIT,MBAY(1,LUN)).NE.0) GOTO100                        
                                                                        
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,IRET)
      IF(IRET.NE.0) GOTO 1
      RETURN
 
C  EOF ON ATTEMPTED READ                                                
C  ---------------------                                                
                                                                        
100   CALL WTSTAT(LUNIT,LUN,IL,0)                                      
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
      SUBSET = ' '                                                      
      JDATE = 0                                                         
      IRET = -1                                                         
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('READERM - FILE IS CLOSED                   ')         
901   CALL BORT('READERM - FILE IS OPEN FOR OUTPUT          ')         
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READERME(MESG,LUNIT,SUBSET,JDATE,IRET)                 
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
      CHARACTER*8 SUBSET,SEC0                                           
      DIMENSION   MESG(*),IEC0(2)                                       
      EQUIVALENCE (SEC0,IEC0)                                           
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      CALL WTSTAT(LUNIT,LUN,IL, 1)                                      
                                                                        
C  READ A MESSAGE INTO A MESSAGE BUFFER                                 
C  ------------------------------------                                 
                                                                        
      IEC0(1) = MESG(1)                                                 
      IEC0(2) = MESG(2)                                                 
      DO I=1,LMSG(SEC0)                                                 
      MBAY(I,LUN) = MESG(I)                                             
      ENDDO                                                             
      IF(SEC0(1:4).NE.'BUFR') GOTO 902                                  
                                                                        
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,IRET)
      RETURN
                                                                        
C  EOF ON ATTEMPTED READ                                                
C  ---------------------                                                
                                                                        
100   CALL WTSTAT(LUNIT,LUN,IL,0)                                      
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
      SUBSET = ' '                                                      
      JDATE = 0                                                         
      IRET = -1                                                         
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('READERME - FILE IS CLOSED                   ')        
901   CALL BORT('READERME - FILE IS OPEN FOR OUTPUT          ')        
902   CALL BORT('READERME - NOT A BUFR FILE                  ')        
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READFT(LUNIT,SUBSET,JDATE,IRET)                        
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
      CHARACTER*8 SEC0,SUBSET                                           
      CHARACTER*4 BUFR                                                  
      DIMENSION   IEC0(2)
      EQUIVALENCE (SEC0,IEC0)
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      CALL WTSTAT(LUNIT,LUN,IL, 1)                                      
                                                                        
C  READ A MESSAGE INTO A MESSAGE BUFFER - SKIP DX MESSAGES              
C  -------------------------------------------------------              
                                                                        
1     MBIT = 0                                                          
      SEC0 = ' '                                                        
      IMSG = 8/NBYTW+1                                                  
      READ(LUNIT,ERR=100,END=100) SEC0,(MBAY(I,LUN),I=IMSG,LMSG(SEC0))  
      CALL CHRTRNA(BUFR,SEC0,4)                                         
      IF(BUFR.NE.'BUFR') GOTO 100                                       
      DO I=1,IMSG-1
      MBAY(I,LUN) = IEC0(I)
      ENDDO
                                                                        
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,IRET)
      IF(IRET.NE.0) GOTO 1
      RETURN
 
C  EOF ON ATTEMPTED READ                                                
C  ---------------------                                                
                                                                        
100   CALL WTSTAT(LUNIT,LUN,IL,0)                                      
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
      SUBSET = ' '                                                      
      JDATE = 0                                                         
      IRET = -1                                                         
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('READFT - FILE IS CLOSED                       ')      
901   CALL BORT('READFT - FILE IS OPEN FOR OUTPUT              ')      
903   CALL BORT('READFT - MSGTYPE MISMATCH  FOR '//SUBSET       )      
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READIBM(LUNIT,SUBSET,JDATE,IRET)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
 
      CHARACTER*8 SEC0,SUBSET
      CHARACTER*4 BUFR
      CHARACTER*1 CBAY(8*5000)
      DIMENSION   JBAY(5000)
      EQUIVALENCE (CBAY(1),JBAY(1))
      EQUIVALENCE (CBAY(1),SEC0)
 
C-----------------------------------------------------------------------
      LBMG(SEC0) = IUPM(SEC0(5:7),24)
C-----------------------------------------------------------------------
 
      IRET = 0
 
C  CHECK THE FILE STATUS
C  ---------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      CALL WTSTAT(LUNIT,LUN,IL,1)
 
C  READ A MESSAGE INTO A MESSAGE BUFFER - SKIP DX MESSAGES
C  -------------------------------------------------------
 
1     SEC0 = ' '
      READ(LUNIT,ERR=902,END=100) SEC0,(CBAY(I),I=9,LBMG(SEC0))
      DO I=1,8
      CBAY(I) = SEC0(I:I)
      ENDDO
      DO I=1,LMSG(SEC0)
      MBAY(I,LUN) = JBAY(I)
      ENDDO
 
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,IRET)
      IF(IRET.NE.0) GOTO 1
      RETURN
 
C  EOF ON ATTEMPTED READ
C  ---------------------
 
100   CALL WTSTAT(LUNIT,LUN,IL,0)
      INODE(LUN) = 0
      IDATE(LUN) = 0
      SUBSET = ' '
      JDATE = 0
      IRET = -1
      RETURN
 
C  ERROR EXITS
C  -----------
 
900   CALL BORT('READIBM - FILE IS CLOSED                       ')
901   CALL BORT('READIBM - FILE IS OPEN FOR OUTPUT              ')
902   CALL BORT('READIBM - I/O ERROR READING MESSAGE            ')
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READMG(LUNIT,SUBSET,JDATE,IRET)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
      COMMON /DATELN/ LENDAT
 
      CHARACTER*8 SEC0,SUBSET
      CHARACTER*4 BUFR
      DIMENSION   IEC0(2)
      EQUIVALENCE (SEC0,IEC0)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      IRET = 0
 
C  CHECK THE FILE STATUS
C  ---------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      CALL WTSTAT(LUNIT,LUN,IL,1)
      IMSG = 8/NBYTW+1
 
C  READ A MESSAGE INTO A MESSAGE BUFFER - SKIP DX MESSAGES
C  -------------------------------------------------------
 
1     SEC0 = ' '
      READ(LUNIT,ERR=902,END=100) SEC0,(MBAY(I,LUN),I=IMSG,LMSG(SEC0))
      CALL CHRTRNA(BUFR,SEC0,4)
      IF(BUFR.NE.'BUFR') GOTO 100
      DO I=1,IMSG-1
      MBAY(I,LUN) = IEC0(I)
      ENDDO
 
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,IRET)
      IF(IRET.NE.0) GOTO 1
      RETURN
 
C  EOF ON ATTEMPTED READ
C  ---------------------
 
100   CALL WTSTAT(LUNIT,LUN,IL,0)
      INODE(LUN) = 0
      IDATE(LUN) = 0
      SUBSET = ' '
      JDATE = 0
      IRET = -1
      RETURN
 
C  ENTRY DATELEN SETS THE LENGTH OF THE DATE INTEGER RETURN FROM READS
C  -------------------------------------------------------------------
 
      ENTRY DATELEN(LEN)
      IF(LEN.NE.8 .AND. LEN.NE.10) CALL BORT('DATELEN - BAD LEN')
      LENDAT = LEN
      RETURN
 
C  ERROR EXITS
C  -----------
 
900   CALL BORT('READMG - FILE IS CLOSED                       ')
901   CALL BORT('READMG - FILE IS OPEN FOR OUTPUT              ')
902   CALL BORT('READMG - I/O ERROR READING MESSAGE            ')
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READMM(IMSG,SUBSET,JDATE,IRET)                         
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
                                                                        
      CHARACTER*8 SUBSET                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE MESSAGE REQUEST AND FILE STATUS                            
C  -----------------------------------------                            
                                                                        
      CALL STATUS(MUNIT,LUN,IL,IM)                                      
      CALL WTSTAT(MUNIT,LUN,IL, 1)                                      
      IF(IL.GE.0) GOTO 900                                              
      IRET = 0                                                          
                                                                        
      IF(IMSG.EQ.0 .OR.IMSG.GT.MSGP(0)) THEN                            
         CALL WTSTAT(MUNIT,LUN,IL,0)                                    
         IRET = -1                                                      
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  READ MESSAGE# IMSG INTO A MESSAGE BUFFER                             
C  ----------------------------------------                             
                                                                        
      IPTR = MSGP(IMSG)                                                 
      IF(IMSG.LT.MSGP(0)) LPTR = MSGP(IMSG+1)-IPTR                      
      IF(IMSG.EQ.MSGP(0)) LPTR = MLAST-IPTR+1                           
      IPTR = IPTR-1                                                     
                                                                        
      DO I=1,LPTR                                                       
      MBAY(I,LUN) = MSGS(IPTR+I)                                        
      ENDDO                                                             
                                                                        
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,JRET)
      NMSG(LUN) = IMSG 
      IMSG = IMSG+1
      RETURN
 
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('READMM - FILE IS NOT OPEN FOR INPUT       ')          
901   CALL BORT('READMM - BAD IPTR IN MEMORY BUFFER        ')          
902   CALL BORT('READMM - BAD NPTR IN MEMORY BUFFER        ')          
903   CALL BORT('READMM - MSGTYPE MISMATCH FOR '//SUBSET    )          
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READNS(LUNIT,SUBSET,JDATE,IRET)                        
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*8  SUBSET                                               
      CHARACTER*3  TYP                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  REFRESH THE SUBSET AND IDATE PARAMETERS                              
C  ---------------------------------------                              
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.GE.0) CALL BORT('READNS - FILE NOT OPEN FOR INPUT')        
      SUBSET = TAG(INODE(LUN))                                          
      JDATE  = IDATE(LUN)                                               
                                                                        
C  READ THE NEXT SUBSET IN THE BUFR FILE                                
C  -------------------------------------                                
                                                                        
1     CALL READSB(LUNIT,IRET)                                           
      IF(IRET.NE.0) THEN                                                
         CALL READMG(LUNIT,SUBSET,JDATE,IRET)                           
         IF(IRET.NE.0) RETURN                                           
         GOTO 1                                                         
      ELSE                                                              
         RETURN                                                         
      ENDIF                                                             
                                                                        
      END                                                               
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READSB(LUNIT,IRET)
 
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)
      COMMON /UNPTYP/ MSGUNP(32)
 
C-----------------------------------------------------------------------
      ENTRY READERS(LUNIT,IRET)
C-----------------------------------------------------------------------
 
      IRET = 0
 
C  CHECK THE FILE STATUS
C  ---------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      IF(IM.EQ.0) THEN
         IRET = -1
         RETURN
      ENDIF
 
C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE
C  ---------------------------------------------
 
      IF(NSUB(LUN).EQ.MSUB(LUN)) THEN
         IRET = -1
         RETURN
      ELSE
         NSUB(LUN) = NSUB(LUN) + 1
      ENDIF
 
C  READ THE NEXT SUBSET AND RESET THE POINTERS
C  -------------------------------------------
 
      IF(MSGUNP(LUN).EQ.0) THEN
         IBIT = MBYT(LUN)*8
         CALL UPB(NBYT,16,MBAY(1,LUN),IBIT)
         CALL RDTREE(LUN)
         MBYT(LUN) = MBYT(LUN) + NBYT
      ELSEIF(MSGUNP(LUN).EQ.1) THEN
         IBIT = MBYT(LUN)
         CALL RDTREE(LUN)
         MBYT(LUN) = IBIT
      ELSEIF(MSGUNP(LUN).EQ.2) THEN
         CALL RDCMPS(LUN)
      ELSE
         GOTO 903
      ENDIF
 
      RETURN
900   CALL BORT('READSB - FILE IS CLOSED             ')
901   CALL BORT('READSB - FILE IS OPEN FOR OUTPUT    ')
902   CALL BORT('READSB - NO MESSAGE OPEN            ')
903   CALL BORT('READSB - UNKNOWN MESSAGE UNPACK TYPE')
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READTJ(LUNIT,SUBSET,JDATE,IRET)                        
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
      CHARACTER*8 SEC0,SUBSET                                           
      CHARACTER*4 BUFR                                                  
      DIMENSION   IEC0(2)
      EQUIVALENCE (SEC0,IEC0)
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      CALL WTSTAT(LUNIT,LUN,IL, 1)                                      
                                                                        
C  READ A MESSAGE INTO A MESSAGE BUFFER - SKIP DX MESSAGES              
C  -------------------------------------------------------              
                                                                        
1     MBIT = 0                                                          
      SEC0 = ' '                                                        
      IMSG = 8/NBYTW+1                                                  
      READ(LUNIT,ERR=902,END=100) SEC0,(MBAY(I,LUN),I=IMSG,LMSG(SEC0))  
      CALL CHRTRNA(BUFR,SEC0,4)                                         
      IF(BUFR.NE.'BUFR') GOTO 100                                       
      DO I=1,IMSG-1
      MBAY(I,LUN) = IEC0(I)
      ENDDO
                                                                        
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,IRET)
      IF(IRET.NE.0) GOTO 1
      RETURN
                                                                        
C  EOF ON ATTEMPTED READ                                                
C  ---------------------                                                
                                                                        
100   CALL WTSTAT(LUNIT,LUN,IL,0)                                      
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
      SUBSET = ' '                                                      
      JDATE = 0                                                         
      IRET = -1                                                         
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('READTJ - FILE IS CLOSED                       ')      
901   CALL BORT('READTJ - FILE IS OPEN FOR OUTPUT              ')      
902   CALL BORT('READTJ - I/O ERROR READING MESSAGE            ')      
903   CALL BORT('READTJ - MSGTYPE MISMATCH  FOR '//SUBSET       )      
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION RJUST(STR)                                               
      CHARACTER*(*) STR                                                 
      RJUST = 0                                                         
      IF(STR.EQ.' ') RETURN                                             
      LSTR = LEN(STR)                                                   
      DO WHILE(STR(LSTR:LSTR).EQ.' ')                                   
         DO I=LSTR,2,-1
         STR(I:I) = STR(I-1:I-1)                                        
         ENDDO
         STR(1:1) = ' '
      ENDDO                                                             
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RSVFVM(NEM1,NEM2)                                      
                                                                        
      CHARACTER*8 NEM1,NEM2                                             
                                                                        
      DO I=1,LEN(NEM1)                                                  
      IF(I.EQ.1) THEN                                                   
         J = 1                                                          
      ELSE                                                              
         IF(NEM1(I:I).EQ.'.') THEN                                      
            NEM1(I:I) = NEM2(J:J)                                       
            J = J+1                                                     
         ENDIF                                                          
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SEQSDX(CARD,LUN)                                       
                                                                        
      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)             
                                                                        
      CHARACTER*80  CARD,SEQS                                           
      CHARACTER*12  ATAG,TAGS(250)                                      
      CHARACTER*8   NEMO,NEMA,NEMB                                      
      CHARACTER*3   TYPS                                                
      CHARACTER*1   REPS,TAB                                            
                                                                        
      DATA MAXTGS /250/                                                 
      DATA MAXTAG /12/                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  FIND THE SEQUENCE TAG IN TABLE D AND PARSE THE SEQUENCE STRING       
C  --------------------------------------------------------------       
                                                                        
      NEMO = CARD( 3:10)                                                
      SEQS = CARD(14:78)                                                
                                                                        
      CALL NEMTAB(LUN,NEMO,IDN,TAB,ISEQ)                                
      CALL PARSEQ(SEQS,TAGS,MAXTGS,NTAG)                                
      IF(TAB.NE.'D') GOTO 900                                           
      IF(NTAG.EQ.0 ) GOTO 900                                           
                                                                        
      DO N=1,NTAG                                                       
      ATAG = TAGS(N)                                                    
      IREP = 0                                                          
                                                                        
C  CHECK FOR REPLICATOR                                                 
C  --------------------                                                 
                                                                        
      DO I=1,5                                                          
      IF(ATAG(1:1).EQ.REPS(I,1)) THEN                                   
         DO J=2,MAXTAG                                                  
         IF(ATAG(J:J).EQ.REPS(I,2)) THEN                                
            IF(J.EQ.MAXTAG) GOTO 901                                    
            CALL STRNUM(ATAG(J+1:MAXTAG),NUMR)                          
            IF(I.EQ.1 .AND. NUMR.LE.0  ) GOTO 901                       
            IF(I.EQ.1 .AND. NUMR.GT.255) GOTO 901                       
            IF(I.NE.1 .AND. NUMR.NE.0  ) GOTO 901                       
            ATAG = ATAG(2:J-1)                                          
            IREP = I                                                    
            GOTO 1                                                      
         ENDIF                                                          
         ENDDO                                                          
         GOTO 901                                                       
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  CHECK FOR VALID TAG                                                  
C  -------------------                                                  
                                                                        
1     IF(NEMOCK(ATAG).NE.0) GOTO 901                                    
      CALL NEMTAB(LUN,ATAG,IDN,TAB,IRET)                                
      IF(IRET.GT.0) THEN                                                
         IF(TAB.EQ.'B' .AND. IREP.NE.0) GOTO 902                        
         IF(ATAG(1:1).EQ.'.') THEN                                      
            NEMB = TAGS(N+1)                                            
            CALL NUMTAB(LUN,IDN,NEMA,TAB,ITAB)                          
            CALL NEMTAB(LUN,NEMB,JDN,TAB,IRET)                          
            CALL RSVFVM(NEMA,NEMB)                                      
            IF(NEMA.NE.ATAG) GOTO 903                                   
            IF(N.GT.NTAG ) GOTO 905                                     
            IF(TAB.NE.'B') GOTO 906                                     
         ENDIF                                                          
      ELSE                                                              
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
C  WRITE THE DESCRIPTOR STRING INTO TABD ARRAY                          
C  -------------------------------------------                          
                                                                        
10    IF(IREP.GT.0) CALL PKTDD(ISEQ,LUN,IDNR(IREP,1)+NUMR,IRET)         
      IF(IRET.LT.0) GOTO 904                                            
      CALL PKTDD(ISEQ,LUN,IDN,IRET)                                     
      IF(IRET.LT.0) GOTO 904                                            
                                                                        
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('SEQSDX - UNDEFINED SEQUENCE: '             //   NEMO) 
901   CALL BORT('SEQSDX - BAD TAG IN SEQUENCE: '            //TAGS(N)) 
902   CALL BORT('SEQSDX - REPLICATED ELEMENTS NOT ALLOWED:' //TAGS(N)) 
903   CALL BORT('SEQSDX - UNDEFINED TAG: '                  //TAGS(N)) 
904   CALL BORT('SEQSDX - TOO MANY DESCRIPTORS IN STRING:'  //   NEMO) 
905   CALL BORT('SEQSDX - FOLLOWING-VALUE LAST IN STRING:'  //   NEMA) 
906   CALL BORT('SEQSDX - FOLLOWING VALUE NOT FROM TABLEB:' //   NEMB) 
      END                                                               
C-----------------------------------------------------------------------
C  SUBROUTINE STANDARD WILL REWRITE A MESSAGE WRITTEN BY THE NCEP BUFR  
C  INTERFACE PROGRAM INTO A MORE STANDARD BUFR FORM. SECTION THREE IS   
C  REWRITTEN AS THE EXPANSION (1 LEVEL DEEP) OF THE NCEP SUBSET SEQUENCE
C  DESCRIPTOR. SECTION FOUR IS REWRITTEN TO CONFORM TO THE NEW SECTION  
C  THREE DESCRIPTOR LIST. ALL SUBSET BYTE COUNTERS AND BIT PADS ARE REMO
C  FROM THE DATA. IF THE 1 LEVEL EXPANSION OF THE SUBSET SEQUENCE DESCRI
C  CONTAINS ONLY STANDARD DESCRIPTORS, THEN THE NEW MESSAGE IS ENTIRLY  
C  AND STRICTLY STANDARD.                                               
C                                                                       
C  THE SUBROUTINE ARGUMENTS ARE:                                        
C                                                                       
C  INPUT:  LUNIT - UNIT OPENED WITH BUFR TABLES USING OPENBF            
C          MSGIN - ARRAY CONTAINING AN NCEP BUFR MESSAGE                
C                                                                       
C  OUTPUT: MSGOT - ARRAY CONTAINING STANDARDIZED FORM OF THE INPUT MESSA
C                                                                       
C  ----------------------------------------------                       
C  NOTE: MSGIN AND MSGOT MUST BE SEPARATE ARRAYS.                       
C  ----------------------------------------------                       
C                                                                       
C-----------------------------------------------------------------------
      SUBROUTINE STANDARD(LUNIT,MSGIN,MSGOT)                            
                                                                        
      DIMENSION MSGIN(*),MSGOT(*)                                       
                                                                        
      CHARACTER*8 SUBSET                                                
      CHARACTER*4 SEVN                                                  
      CHARACTER*1 TAB                                                   
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  LUNIT MUST POINT TO AN OPEN BUFR FILE                                
C  -------------------------------------                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) CALL BORT('STANDARD - NO OPEN BUFR FILE!!')          
                                                                        
C  IDENTIFY THE SECTION LENGTHS AND ADDRESSES IN MSGIN                  
C  ---------------------------------------------------                  
                                                                        
      IAD0 = 0                                                          
      LEN0 = 8                                                          
      LENN = IUPB(MSGIN,IAD0+5,24)                                      
                                                                        
      IAD1 = IAD0+LEN0                                                  
      LEN1 = IUPB(MSGIN,IAD1+1,24)                                      
      LEN2 = IUPB(MSGIN,IAD1+8,1)                                       
                                                                        
      IAD2 = IAD1+LEN1                                                  
      LEN2 = IUPB(MSGIN,IAD2+1,24)*LEN2                                 
                                                                        
      IAD3 = IAD2+LEN2                                                  
      LEN3 = IUPB(MSGIN,IAD3+1,24)                                      
                                                                        
      IAD4 = IAD3+LEN3                                                  
      LEN4 = IUPB(MSGIN,IAD4+1,24)                                      
                                                                        
      LENM = LEN0+LEN1+LEN2+LEN3+LEN4+4                                 
                                                                        
      IF(LENN.NE.LENM) CALL BORT('STANDARD - BAD INPUT BYTE COUNTS')   
                                                                        
      MBIT = (LENN-4)*8                                                 
      CALL UPC(SEVN,4,MSGIN,MBIT)                                       
      IF(SEVN.NE.'7777') CALL BORT('STANDARD - CANT FIND 7777')        
                                                                        
C  COPY SECTIONS 0 THROUGH PART OF SECTION 3 INTO MSGOT                 
C  ----------------------------------------------------                 
                                                                        
      CALL MVB(MSGIN,1,MSGOT,1,LEN0+LEN1+LEN2+7)                        
                                                                        
C  REWRITE NEW SECTION 3 IN A "STANDARD" FORM                           
C  ------------------------------------------                           
                                                                        
      NSUB = IUPB(MSGIN,IAD3+ 5,16)                                     
      ISUB = IUPB(MSGIN,IAD3+10,16)                                     
      IBIT = (IAD3+7)*8                                                 
                                                                        
C  LOOK UP THE SUBSET DESCRIPTOR AND ITS LENGTH IN DESCRIPTORS          
C  -----------------------------------------------------------          
                                                                        
      CALL NUMTAB(LUN,ISUB,SUBSET,TAB,ITAB)                             
      IF(ITAB.EQ.0) CALL BORT('STANDARD - UNKNOWN SUBSET DESCRIPTOR')  
      CALL UPTDD(ITAB,LUN,0,NSEQ)                                       
                                                                        
C  COPY EACH DESCRIPTOR IN THE SUBSET SEQUENCE INTO THE NEW SECTION 3   
C  ------------------------------------------------------------------   
                                                                        
      DO N=1,NSEQ                                                       
      CALL UPTDD(ITAB,LUN,N,IDSC)                                       
      CALL PKB(IDSC,16,MSGOT,IBIT)                                      
      IF(N.EQ.NSEQ) CALL PKB(0,8,MSGOT,IBIT)                            
      ENDDO                                                             
                                                                        
      IBIT = IAD3*8                                                     
      LEN3 = 8+NSEQ*2                                                   
      NAD4 = IAD3+LEN3                                                  
      CALL PKB(LEN3,24,MSGOT,IBIT)                                      
                                                                        
C  NOW THE TRICKY PART - NEW SECTION 4                                  
C  -----------------------------------                                  
                                                                        
      IBIT = (IAD4+4)*8                                                 
      JBIT = (NAD4+4)*8                                                 
                                                                        
C  COPY THE SUBSETS, MINUS THE BYTE COUNTER AND PAD, INTO THE NEW SECTIO
C  ---------------------------------------------------------------------
                                                                        
      DO 10 I=1,NSUB                                                    
      CALL UPB(LSUB,16,MSGIN,IBIT)                                      
                                                                        
      DO L=1,LSUB-2                                                     
      CALL UPB(NVAL,8,MSGIN,IBIT)                                       
      CALL PKB(NVAL,8,MSGOT,JBIT)                                       
      ENDDO                                                             
                                                                        
      DO K=1,8                                                          
      KBIT = IBIT-K-8                                                   
      CALL UPB(KVAL,8,MSGIN,KBIT)                                       
      IF(KVAL.EQ.K) THEN                                                
         JBIT = JBIT-K-8                                                
         GOTO 10                                                        
      ENDIF                                                             
      ENDDO                                                             
      CALL BORT('STANDARD - KBIT ERROR')                               
                                                                        
10    ENDDO                                                             
                                                                        
C  MAKE SURE NEW SECTION 4 HAS AN EVEN NUMBER OF BYTES AND ENTER THE COU
C  ---------------------------------------------------------------------
                                                                        
      DO WHILE(.NOT.(MOD(JBIT,8).EQ.0 .AND. MOD(JBIT/8,2).EQ.0))        
         CALL PKB(0,1,MSGOT,JBIT)                                       
      ENDDO                                                             
                                                                        
      IBIT = NAD4*8                                                     
      LEN4 = JBIT/8 - NAD4                                              
      CALL PKB(LEN4,24,MSGOT,IBIT)                                      
                                                                        
C  FINISH THE NEW MESSAGE WITH AN UPDATED SECTION-0 BYTE COUNT          
C  -----------------------------------------------------------          
                                                                        
      IBIT = 32                                                         
      LENM = LEN0+LEN1+LEN2+LEN3+LEN4+4                                 
      CALL PKB(LENM,24,MSGOT,IBIT)                                      
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      CALL PKC('7777', 4,MSGOT,JBIT)                                    
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE STATUS(LUNIT,LUN,IL,IM)                                

      PARAMETER (NFILES=32)
                                                                        
      COMMON /STBFR/ IOLUN(NFILES),IOMSG(NFILES)                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IF(LUNIT.LE.0 .OR. LUNIT.GT.99) GOTO 900                          
                                                                        
C  CLEAR THE STATUS INDICATORS                                          
C  ---------------------------                                          
                                                                        
      LUN = 0                                                           
      IL  = 0                                                           
      IM  = 0                                                           
                                                                        
C  SEE IF THE UNIT IS DEFINED                                           
C  --------------------------                                           
                                                                        
      DO I=1,NFILES                                                     
      IF(ABS(IOLUN(I)).EQ.LUNIT) LUN = I                                
      ENDDO                                                             
                                                                        
C  IF NOT, CHECK FOR FILE SPACE - RETURN LUN=0 IF NO FILE SPACE         
C  ------------------------------------------------------------         
                                                                        
      IF(LUN.EQ.0) THEN                                                 
         DO I=1,NFILES                                                  
         IF(IOLUN(I).EQ.0) LUN = I                                      
         IF(IOLUN(I).EQ.0) RETURN                                       
         ENDDO                                                          
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  IF FILE DEFINED RETURN STATUSES                                      
C  -------------------------------                                      
                                                                        
      IL = SIGN(1,IOLUN(LUN))                                           
      IM = IOMSG(LUN)                                                   
                                                                        
      RETURN                                                            
900   CALL BORT('STATUS - ILLEGAL UNIT GIVEN')                         
      END                                                               
C---------------------------------------------------------------------- 
C  ENTRY TO RESET STRING CACHE                                          
C---------------------------------------------------------------------- 
      SUBROUTINE STRCLN                                                 
      PARAMETER(MXS=1000)                                               
      COMMON /STCACH/ MSTR,NSTR,LSTR,LUNS(MXS,2),USRS(MXS),ICON(52,MXS)    
      CHARACTER*80 USRS                                                 
                                                                        
      MSTR = MXS                                                        
      NSTR = 0                                                          
      LSTR = 0                                                          
      RETURN                                                            
      END                                                               
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE STRING(STR,LUN,I1,IO)
 
      PARAMETER (MXS=1000,JCONS=52)
 
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /STCACH/ MSTR,NSTR,LSTR,LUX(MXS,2),USR(MXS),ICON(52,MXS)
      COMMON /USRSTR/ JCON(JCONS)
      COMMON /STORDS/ IORD(MXS),IORX(MXS)
 
      CHARACTER*(*) STR
      CHARACTER*80  USR,UST
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------

      NXT = 0
      UST = STR
      IND = INODE(LUN)
      IF(LEN(STR).GT.80) GOTO 900
 
C  SEE IF STRING IS IN THE CACHE
C  -----------------------------

      DO N=1,NSTR
      IF(LUX(IORD(N),2).EQ.IND) THEN
         IORX(NXT+1) = IORD(N)
         NXT = NXT+1
      ENDIF
      ENDDO
      DO N=1,NXT 
      IF(UST.EQ.USR(IORX(N)))GOTO1
      ENDDO
      GOTO2
 
C  IF IT IS COPY PARAMETERS FROM THE CACHE
C  ---------------------------------------
 
1     DO J=1,JCONS
      JCON(J) = ICON(J,IORX(N))
      ENDDO
      GOTO 100
 
C  IF NOT PARSE IT AND PUT IT THERE
C  --------------------------------
 
2     CALL PARUSR(STR,LUN,I1,IO)
      LSTR = MAX(MOD(LSTR+1,MSTR+1),1)
      NSTR = MIN(NSTR+1,MSTR)
      LUX(LSTR,1) = LUN
      LUX(LSTR,2) = IND
      USR(LSTR) = STR
      DO J=1,JCONS
      ICON(J,LSTR) = JCON(J)
      ENDDO
 
C  REARRANGE THE CACHE ORDER AFTER AN UPDATE
C  -----------------------------------------
 
      DO N=NSTR,2,-1
      IORD(N) = IORD(N-1)
      ENDDO
      IORD(1) = LSTR
 
C  NORMAL AND ERROR EXITS
C  ----------------------
 
100   IF(JCON(1).GT.I1) GOTO 901
      RETURN
900   CALL BORT('STRING - USER STRING > 80 CHARS         :'//UST)
901   CALL BORT('STRING - MUST BE AT LEAST I1 STORE NODES:'//UST)
      END
C-----------------------------------------------------------------------
C  INTEGER FROM A STRING                                                
C-----------------------------------------------------------------------
      SUBROUTINE STRNUM(STR,NUM)                                        
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*20  STR2                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NUM = 0                                                           
      K = 0                                                             
                                                                        
      CALL STRSUC(STR,STR2,NUM)                                         
                                                                        
      DO I=1,NUM                                                        
      READ(STR(I:I),'(I1)',ERR=99) J                                    
      IF(J.EQ.0 .AND. STR(I:I).NE.'0') GOTO 99                          
      K = K*10+J                                                        
      ENDDO                                                             
                                                                        
      NUM = K                                                           
      RETURN                                                            
                                                                        
99    NUM = -1                                                          
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C  DEAD SPACE FORM A STRING                                             
C-----------------------------------------------------------------------
      SUBROUTINE STRSUC(STR1,STR2,LENS)                                 
                                                                        
      CHARACTER*(*) STR1,STR2                                           
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      LENS = 0                                                          
      LSTR = LEN(STR1)                                                  
                                                                        
      DO I=1,LSTR                                                       
      IF(STR1(I:I).NE.' ') GOTO 2                                       
      ENDDO                                                             
      RETURN                                                            
                                                                        
2     DO J=I,LSTR                                                       
      IF(STR1(J:J).EQ.' ') GOTO 3                                       
      LENS = LENS+1                                                     
      STR2(LENS:LENS) = STR1(J:J)                                       
      ENDDO                                                             
      RETURN                                                            
                                                                        
3     DO I=J,LSTR                                                       
      IF(STR1(I:I).NE.' ') LENS = -1                                    
      ENDDO                                                             
      RETURN                                                            
                                                                        
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE TABENT(LUN,NEMO,TAB,ITAB,IREP,IKNT,JUM0)

      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),
     .                JUMP(15000),LINK(15000),JMPB(15000),
     .                IBT(15000),IRF(15000),ISC(15000),
     .                ITP(15000),VALI(15000),KNTI(15000),
     .                ISEQ(15000,2),JSEQ(15000)
      COMMON /TABCCC/ ICDW,ICSC,ICRV

      CHARACTER*24 UNIT
      CHARACTER*10 TAG,RTAG
      CHARACTER*8  NEMO
      CHARACTER*3  TYP,TYPS,TYPT
      CHARACTER*1  REPS,TAB

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  MAKE A JUMP/LINK TABLE ENTRY FOR A REPLICATOR
C  ---------------------------------------------

      IF(IREP.NE.0) THEN
         RTAG = REPS(IREP,1)//NEMO
         DO I=1,10
         IF(RTAG(I:I).EQ.' ') THEN
            RTAG(I:I) = REPS(IREP,2)
            CALL INCTAB(RTAG,TYPS(IREP,1),NODE)
            JUMP(NODE) = NODE+1
            JMPB(NODE) = JUM0
            LINK(NODE) = 0
            IBT (NODE) = LENS(IREP)
            IRF (NODE) = 0
            ISC (NODE) = 0
            IF(IREP.EQ.1) IRF(NODE) = IKNT
            JUM0 = NODE
            GOTO 1
         ENDIF
         ENDDO
         GOTO 900
      ENDIF

C  MAKE AN JUMP/LINK ENTRY FOR AN ELEMENT OR A SEQUENCE
C  ----------------------------------------------------

1     IF(TAB.EQ.'B') THEN
         CALL NEMTBB(LUN,ITAB,UNIT,ISCL,IREF,IBIT)
         IF(UNIT(1:5).EQ.'CCITT') TYPT = 'CHR'
         IF(UNIT(1:5).NE.'CCITT') TYPT = 'NUM'
         CALL INCTAB(NEMO,TYPT,NODE)
         JUMP(NODE) = 0
         JMPB(NODE) = JUM0
         LINK(NODE) = 0
         IBT (NODE) = IBIT
         IRF (NODE) = IREF
         ISC (NODE) = ISCL
         IF(UNIT(1:5).EQ.'CODE') TYPT = 'COD'
         IF(UNIT(1:5).EQ.'FLAG') TYPT = 'FLG'
         IF(TYPT.EQ.'NUM') THEN
            IBT(NODE) = IBT(NODE)+ICDW
            ISC(NODE) = ISC(NODE)+ICSC
         ENDIF
      ELSEIF(TAB.EQ.'D') THEN
         IF(IREP.EQ.0) TYPT = 'SEQ'
         IF(IREP.NE.0) TYPT = TYPS(IREP,2)
         CALL INCTAB(NEMO,TYPT,NODE)
         JUMP(NODE) = NODE+1
         JMPB(NODE) = JUM0
         LINK(NODE) = 0
         IBT (NODE) = 0
         IRF (NODE) = 0
         ISC (NODE) = 0
      ELSE
         GOTO 901
      ENDIF

      RETURN
900   CALL BORT('TABENT - REPLICATOR ERROR: '//RTAG//':'//NEMO)
901   CALL BORT('TABENT - UNDEFINED TAG   : '           //NEMO)
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE TABSUB(LUN,NEMO)

      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),
     .                JUMP(15000),LINK(15000),JMPB(15000),
     .                IBT(15000),IRF(15000),ISC(15000),
     .                ITP(15000),VALI(15000),KNTI(15000),
     .                ISEQ(15000,2),JSEQ(15000)
      COMMON /TABCCC/ ICDW,ICSC,ICRV

      CHARACTER*10 TAG
      CHARACTER*8  NEMO,NEMS,NEM
      CHARACTER*3  TYP
      CHARACTER*1  TAB
      DIMENSION    NEM(250,10),IRP(250,10),KRP(250,10)
      DIMENSION    DROP(10),JMP0(10),NODL(10),NTAG(10,2)
      LOGICAL      DROP

      DATA MAXLIM /10/

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK THE MNEMONIC
C  ------------------

      CALL NEMTAB(LUN,NEMO,IDN,TAB,ITAB)
      IF(TAB.NE.'D') GOTO 900

C  STORE A SUBSET NODE AND JUMP/LINK THE TREE
C  ------------------------------------------

      CALL INCTAB(NEMO,'SUB',NODE)
      JUMP(NODE) = NODE+1
      JMPB(NODE) = 0
      LINK(NODE) = 0
      IBT (NODE) = 0
      IRF (NODE) = 0
      ISC (NODE) = 0

      CALL NEMTBD(LUN,ITAB,NSEQ,NEM(1,1),IRP(1,1),KRP(1,1))
      NTAG(1,1) = 1
      NTAG(1,2) = NSEQ
      JMP0(1)   = NODE
      LIMB      = 1

      ICDW = 0
      ICSC = 0
      ICRV = 0

C  THIS LOOP RESOLVES ENTITIES IN A SUBSET BY EMULATING RECURSION
C  --------------------------------------------------------------

1     DO N=NTAG(LIMB,1),NTAG(LIMB,2)

      NTAG(LIMB,1) = N+1
      NODL(LIMB)   = NTAB+1
      DROP(LIMB)   = N.EQ.NTAG(LIMB,2)

      CALL NEMTAB(LUN,NEM(N,LIMB),IDN,TAB,ITAB)
      NEMS = NEM(N,LIMB)

C  SPECIAL TREATMENT FOR CERTAIN OPERATOR DESCRIPTORS (TAB=C)
C  ----------------------------------------------------------

      IF(TAB.EQ.'C') THEN
         NODL(LIMB) = NTAB
         READ(NEMS,'(3X,I3)') IYYY
         IF(ITAB.EQ.1) THEN
            ICDW = IYYY-128
            IF(IYYY.EQ.0) ICDW = 0
         ELSEIF(ITAB.EQ.2) THEN
            ICSC = IYYY-128
            IF(IYYY.EQ.0) ICSC = 0
         ENDIF
      ELSE
         IREP = IRP(N,LIMB)
         IKNT = KRP(N,LIMB)
         JUM0 = JMP0(LIMB)
         CALL TABENT(LUN,NEMS,TAB,ITAB,IREP,IKNT,JUM0)
      ENDIF

      IF(TAB.EQ.'D') THEN
         LIMB = LIMB+1
         IF(LIMB.GT.MAXLIM) GOTO 901
         CALL NEMTBD(LUN,ITAB,NSEQ,NEM(1,LIMB),IRP(1,LIMB),KRP(1,LIMB))
         NTAG(LIMB,1) = 1
         NTAG(LIMB,2) = NSEQ
         JMP0(LIMB)   = NTAB
         GOTO 1
      ELSEIF(DROP(LIMB)) THEN
2        LINK(NODL(LIMB)) = 0
         LIMB = LIMB-1
         IF(LIMB.EQ.0 ) THEN
            IF(ICDW.NE.0) GOTO 902
            IF(ICSC.NE.0) GOTO 903
            RETURN
         ENDIF
         IF(DROP(LIMB)) GOTO 2
         LINK(NODL(LIMB)) = NTAB+1
         GOTO 1
      ELSEIF(TAB.NE.'C') THEN
         LINK(NODL(LIMB)) = NTAB+1
      ENDIF

      ENDDO

      CALL BORT('TABSUB - SHOULD NOT GET HERE               ')
900   CALL BORT('TABSUB - SUBSET NODE NOT IN TABLE D: '//NEMO)
901   CALL BORT('TABSUB - TOO MANY LIMBS                    ')
902   CALL BORT('TABSUB - CHANGE DATA WIDTH OPERATOR NOT CANCELED')
903   CALL BORT('TABSUB - CHANGE DATA SCALE OPERATOR NOT CANCELED')
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE TRYBUMP(LUNIT,LUN,USR,I1,I2,IO,IRET)                   
                                                                        
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
                                                                        
      REAL*8 USR(I1,I2),VAL                                             
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  SEE IF THERE IS A DRP GROUP INVOLVED                                 
C  ------------------------------------                                 
                                                                        
      NDRP = LSTJPB(NODS(1),LUN,'DRP')                                  
      IF(NDRP.LE.0) RETURN                                              
                                                                        
C  IF SO, CLEAN IT OUT, BUMP IT TO I2, AND TRY UFBRW AGAIN              
C  -------------------------------------------------------              
                                                                        
      INVN = INVWIN(NDRP,LUN,1,NVAL(LUN))                               
      VAL(INVN,LUN) = 0                                                 
      JNVN = INVN+1                                                     
      DO WHILE(NINT(VAL(JNVN,LUN)).GT.0)                                
         JNVN = JNVN+NINT(VAL(JNVN,LUN))                                
      ENDDO                                                             
      DO KNVN=1,NVAL(LUN)-JNVN+1                                        
      INV(INVN+KNVN,LUN) = INV(JNVN+KNVN-1,LUN)                         
      VAL(INVN+KNVN,LUN) = VAL(JNVN+KNVN-1,LUN)                         
      ENDDO                                                             
      NVAL(LUN) = NVAL(LUN)-(JNVN-INVN-1)                               
      CALL USRTPL(LUN,INVN,I2)                                          
      CALL UFBRW(LUN,USR,I1,I2,IO,IRET)                                 
                                                                        
      RETURN                                                            
900   CALL BORT('TRYBUMP - ATTEMPT TO BUMP NON-ZERO REP FACTOR')       
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBCNT(LUNIT,KMSG,KSUB)                                
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUS - RETURN THE MESSAGE AND SUBSET COUNTERS       
C  --------------------------------------------------------------       
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) CALL BORT('UFBCNT - FILE IS CLOSED')                 
      KMSG = NMSG(LUN)                                                  
      KSUB = NSUB(LUN)                                                  
      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBCPY(LUBIN,LUBOT)                                    
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE FILE STATUSES AND I-NODE                                   
C  ----------------------------------                                   
                                                                        
      CALL STATUS(LUBIN,LUI,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUI).NE.INV(1,LUI)) GOTO 902                             
                                                                        
      CALL STATUS(LUBOT,LUO,IL,IM)                                      
      IF(IL.EQ.0) GOTO 903                                              
      IF(IM.EQ.0) GOTO 904                                              
      IF(INODE(LUI).NE.INODE(LUO)) GOTO 905                             
                                                                        
C  EVERYTHING OKAY COPY USER ARRAY FROM LUI TO LUO                      
C  -----------------------------------------------                      
                                                                        
      NVAL(LUO) = NVAL(LUI)                                             
                                                                        
      DO N=1,NVAL(LUI)                                                  
      INV(N,LUO) = INV(N,LUI)                                           
      VAL(N,LUO) = VAL(N,LUI)                                           
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('UFBCPY - INPUT  FILE IS NOT OPEN             ')       
901   CALL BORT('UFBCPY - INPUT  MESG IS NOT OPEN             ')       
902   CALL BORT('UFBCPY - INPUT  I-NODE  MISMATCH             ')       
903   CALL BORT('UFBCPY - OUTPUT FILE IS NOT OPEN             ')       
904   CALL BORT('UFBCPY - OUTPUT MESG IS NOT OPEN             ')       
905   CALL BORT('UFBCPY - IN/OUT I-NODE  MISMATCH             ')       
      END                                                               
C---------------------------------------------------------------------- 
C  WILL MAKE ONE COPY OF EACH UNIQUE ELEMENT IN AN INPUT SUBSET BUFFER  
C  INTO IDENTICAL MNEMONIC SLOT IN THE OUTPUT SUBSET BUFFER             
C---------------------------------------------------------------------- 
      SUBROUTINE UFBCUP(LUBIN,LUBOT)                                    
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG,TAGI(15000),TAGO                                 
      CHARACTER*3  TYP                                                  
      DIMENSION    NINI(15000)                                          
      REAL*8       VAL                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE FILE STATUSES AND I-NODE                                   
C  ----------------------------------                                   
                                                                        
      CALL STATUS(LUBIN,LUI,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUI).NE.INV(1,LUI)) GOTO 902                             
                                                                        
      CALL STATUS(LUBOT,LUO,IL,IM)                                      
      IF(IL.EQ.0) GOTO 903                                              
      IF(IM.EQ.0) GOTO 904                                              
                                                                        
C  MAKE A LIST OF UNIQUE TAGS IN INPUT BUFFER                           
C  ------------------------------------------                           
                                                                        
      NTAG = 0                                                          
                                                                        
      DO 5 NI=1,NVAL(LUI)                                               
      NIN = INV(NI,LUI)                                                 
      IF(ITP(NIN).GE.2) THEN                                            
         DO NV=1,NTAG                                                   
         IF(TAGI(NV).EQ.TAG(NIN)) GOTO 5                                
         ENDDO                                                          
         NTAG = NTAG+1                                                  
         NINI(NTAG) = NI                                                
         TAGI(NTAG) = TAG(NIN)                                          
      ENDIF                                                             
5     ENDDO                                                             
                                                                        
      IF(NTAG.EQ.0) GOTO 905                                            
                                                                        
C  GIVEN A list MAKE ONE COPY OF COMMON ELEMENTS TO OUTPUT BUFFER       
C  --------------------------------------------------------------       
                                                                        
      DO 10 NV=1,NTAG                                                   
      NI = NINI(NV)                                                     
      DO NO=1,NVAL(LUO)                                                 
      TAGO = TAG(INV(NO,LUO))                                           
      IF(TAGI(NV).EQ.TAGO) THEN                                         
         VAL(NO,LUO) = VAL(NI,LUI)                                      
         GOTO 10                                                        
      ENDIF                                                             
      ENDDO                                                             
10    ENDDO                                                             
                                                                        
C  ALL EXITS HERE                                                       
C  --------------                                                       
                                                                        
      RETURN                                                            
900   CALL BORT('UFBCUP - INPUT  FILE IS NOT OPEN             ')       
901   CALL BORT('UFBCUP - INPUT  MESG IS NOT OPEN             ')       
902   CALL BORT('UFBCUP - INPUT  I-NODE  MISMATCH             ')       
903   CALL BORT('UFBCUP - OUTPUT FILE IS NOT OPEN             ')       
904   CALL BORT('UFBCUP - OUTPUT MESG IS NOT OPEN             ')       
905   CALL BORT('UFBCUP - NO TAGS IN INPUT BUFFER             ')       
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBDMP(LUNIT,LUPRT)                                    
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG,TG                                               
      CHARACTER*8  VC                                                   
      CHARACTER*3  TYP,TP                                               
      CHARACTER*1  YOU                                                  
      EQUIVALENCE  (VL,VC)                                              
      REAL*8       VAL,VL,BMISS
                                                                        
      DATA BMISS /10E10/                                                
                                                                        
C---------------------------------------------------------------------- 
CFFP$ EXPAND (STATUS)                                                   
C---------------------------------------------------------------------- 
                                                                        
      if(luprt.eq.0) luout = 6                                          
      if(luprt.ne.0) luout = luprt                                      
                                                                        
C  CHECK THE FILE STATUS AND I-NODE                                     
C  --------------------------------                                     
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
                                                                        
C  DUMP THE CONTENTS OF COMMON /USRINT/ FOR UNIT LUNIT                  
C  ---------------------------------------------------                  
                                                                        
      DO NV=1,NVAL(LUN)                                                 
      if(luprt.eq.0 .and. mod(nv,20).eq.0) then                         
         print*,'(MORE)'                                                
         read(5,'(a1)') you                                             
         if(you.eq.'q') return                                          
      endif                                                             
      ND = INV (NV,LUN)                                                 
      VL = VAL (NV,LUN)                                                 
      TG = TAG (ND)                                                     
      TP = TYP (ND)                                                     
      IT = ITP (ND)                                                     
      IB = IBT (ND)                                                     
      IS = ISC (ND)                                                     
      IR = IRF (ND)                                                     
      JP = JUMP(ND)                                                     
      LK = LINK(ND)                                                     
      JB = JMPB(ND)                                                     
      RJ = RJUST(TG)                                                    
      IF(TP.NE.'CHR' .AND. VL.LT.BMISS) THEN                      
         WRITE(LUOUT,1) NV,TP,IT,TG,VL,IB,IS,IR,ND,JP,LK,JB             
      ELSE                                                              
         IF(VL.EQ.BMISS) THEN 
            VC = 'MISSING'                                 
         ELSE
            VC = VC(1:IB/8)                                
         ENDIF
         RJ = RJUST(VC)                                                 
         WRITE(LUOUT,2) NV,TP,IT,TG,VC,IB,IS,IR,ND,JP,LK,JB             
      ENDIF                                                             
      ENDDO                                                             
                                                                        
1     FORMAT(I5,1X,A3,'-',I1,1X,A10,1X,F10.1,7(1X,I5))                  
2     FORMAT(I5,1X,A3,'-',I1,1X,A10,1X,A10  ,7(1X,I5))                  
                                                                        
                                                                        
C  EXITS                                                                
C  -----                                                                
                                                                        
      RETURN                                                            
900   CALL BORT('UFBDMP - FILE IS CLOSED                     ')        
901   CALL BORT('UFBDMP - NO MESSAGE OPEN                    ')        
902   CALL BORT('UFBDMP - I-NODE MISMATCH                    ')        
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBEVN(LUNIT,USR,I1,I2,I3,IRET,STR)                    
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*10  TAG                                                 
      CHARACTER*3   TYP                                                 
      DIMENSION     USR(I1,I2,I3),INVN(255)                             
      REAL*8        VAL,USR,BMISS                                       
                                                                        
      DATA BMISS /10E10/                                                
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE FILE STATUS AND I-NODE                                     
C  --------------------------------                                     
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
                                                                        
C  PARSE THE INPUT STRING                                               
C  ----------------------                                               
                                                                        
      CALL STRING(STR,LUN,I1,0)                                         
                                                                        
C  SET INITIAL VALUES FOR RETURNING ARGUMENTS                           
C  ------------------------------------------                           
                                                                        
      DO I=1,I1*I2*I3                                                   
      USR(I,1,1) = BMISS                                                
      ENDDO                                                             
                                                                        
      IRET = 0                                                          
                                                                        
C  LOOP OVER COND WINDOWS                                               
C  ----------------------                                               
                                                                        
      INC1 = 1                                                          
      INC2 = 1                                                          
                                                                        
1     CALL CONWIN(LUN,INC1,INC2,I2)                                     
      IF(NNOD.EQ.0) THEN                                                
         IRET = I2                                                      
         RETURN                                                         
      ELSEIF(INC1.EQ.0) THEN                                            
         RETURN                                                         
      ELSE                                                              
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            INS2 = INC1                                                 
            CALL GETWIN(NODS(I),LUN,INS1,INS2)                          
            IF(INS1.EQ.0) RETURN                                        
            GOTO 2                                                      
         ENDIF                                                          
         ENDDO                                                          
         INS1 = INC1                                                    
         INS2 = INC2                                                    
      ENDIF                                                             
                                                                        
C  READ PUSH DOWN STACK DATA INTO 3D ARRAYS                             
C  ----------------------------------------                             
                                                                        
2     IRET = IRET+1                                                     
      IF(IRET.LE.I2) THEN                                               
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            NNVN = NVNWIN(NODS(I),LUN,INS1,INS2,INVN,I3)                
            DO N=1,NNVN                                                 
            USR(I,IRET,N) = VAL(INVN(N),LUN)                            
            ENDDO                                                       
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  DECIDE WHAT TO DO NEXT                                               
C  ----------------------                                               
                                                                        
      CALL NXTWIN(LUN,INS1,INS2)                                        
      IF(INS1.GT.0 .AND. INS1.LT.INC2) GOTO 2                           
      IF(NCON.GT.0) GOTO 1                                              
                                                                        
      RETURN                                                            
900   CALL BORT('UFBEVN - FILE IS CLOSED                     ')        
901   CALL BORT('UFBEVN - NO MESSAGE OPEN                    ')        
902   CALL BORT('UFBEVN - I-NODE MISMATCH                    ')        
      END                                                               
C---------------------------------------------------------------------- 
C  WILL GET (UNPACK/RETURN) 1-D DESCRIPTORS IN THE INPUT STRING WITHOUT 
C  ADVANCING THE SUBSET POINTER                                         
C---------------------------------------------------------------------- 
      SUBROUTINE UFBGET(LUNIT,TAB,I1,IRET,STR)                          
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRBIT/ NBIT(15000),MBIT(15000)                           
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*10  TAG                                                 
      CHARACTER*8   CVAL                                                
      CHARACTER*3   TYP                                                 
      DIMENSION     TAB(I1)                                             
      EQUIVALENCE   (CVAL,RVAL)
      REAL*8        VAL,RVAL,TAB,BMISS
                                                                        
      DATA BMISS /10E10/                                                
      DATA MAXTG /100/                                                  
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND (UPBB,USRTPL,INVWIN)                                       
C-----------------------------------------------------------------------
      MPS(NODE) = 2**(IBT(NODE))-1                                      
      UPS(NODE) = (IVAL+IRF(NODE))*10.**(-ISC(NODE))                    
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          

      DO I=1,I1
      TAB(I) = BMISS
      ENDDO
                                                                        
C  MAKE SURE A FILE/MESSAGE IS OPEN FOR INPUT                           
C  ------------------------------------------                           
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.GE.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
                                                                        
C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE                        
C  ---------------------------------------------                        
                                                                        
      IF(NSUB(LUN).EQ.MSUB(LUN)) THEN                                   
         IRET = -1                                                      
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  PARSE THE STRING                                                     
C  ----------------                                                     
                                                                        
      CALL STRING(STR,LUN,I1,0)                                         
                                                                        
C  EXPAND THE TEMPLATE FOR THIS SUBSET AS LITTLE AS POSSIBLE            
C  ---------------------------------------------------------            
                                                                        
      N = 1                                                             
      NBIT(N) = 0                                                       
      MBIT(N) = MBYT(LUN)*8 + 16                                        
      CALL USRTPL(LUN,N,N)                                              
                                                                        
10    DO N=N+1,NVAL(LUN)                                                
      NODE = INV(N,LUN)                                                 
      NBIT(N) = IBT(NODE)                                               
      MBIT(N) = MBIT(N-1)+NBIT(N-1)                                     
      IF(NODE.EQ.NODS(NNOD)) THEN                                       
         NVAL(LUN) = N                                                  
         GOTO 20                                                        
      ELSEIF(ITP(NODE).EQ.1) THEN                                       
         CALL UPBB(IVAL,NBIT(N),MBIT(N),MBAY(1,LUN))                    
         CALL USRTPL(LUN,N,IVAL)                                        
         GOTO 10                                                        
      ENDIF                                                             
      ENDDO                                                             
20    CONTINUE                                                          
                                                                        
C  UNPACK ONLY THE NODES FOUND IN THE STRING                            
C  -----------------------------------------                            
                                                                        
      DO I=1,NNOD                                                       
      NODE = NODS(I)                                                    
      INVN = INVWIN(NODE,LUN,1,NVAL(LUN))                               
      IF(INVN.GT.0) THEN                                                
         CALL UPBB(IVAL,NBIT(INVN),MBIT(INVN),MBAY(1,LUN))              
         IF(ITP(NODE).EQ.1) THEN                                        
            TAB(I) = IVAL                                               
         ELSEIF(ITP(NODE).EQ.2) THEN                                    
            IF(IVAL.LT.MPS(NODE)) TAB(I) = UPS(NODE)                    
         ELSEIF(ITP(NODE).EQ.3) THEN                                 
            CVAL = ' '                                               
            KBIT = MBIT(INVN)                                        
            CALL UPC(CVAL,NBIT(INVN)/8,MBAY(1,LUN),KBIT)                   
            TAB(I) = RVAL                                       
         ENDIF                                                          
      ELSE                                                              
         TAB(I) = BMISS                                                 
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('FILE NOT OPEN FOR INPUT')                             
901   CALL BORT('NO MESSAGE OPEN        ')                             
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBINT(LUNIN,USR,I1,I2,IRET,STR)                       
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*(*) STR                                                 
      DIMENSION     USR(I1,I2)                                          
      REAL*8        USR,VAL,BMISS
                                                                        
      DATA BMISS /10E10/                                                
                                                                        
C---------------------------------------------------------------------- 
CFPP$ EXPAND (STATUS,UFBRW)                                             
C---------------------------------------------------------------------- 
                                                                        
      IRET = 0                                                          
      IF(I1.LE.0) RETURN                                                
      IF(I2.LE.0) RETURN                                                
                                                                        
C  CHECK THE FILE STATUS AND I-NODE AND PARSE OR RECALL THE STRING      
C  ---------------------------------------------------------------      
                                                                        
      LUNIT = ABS(LUNIN)                                                
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
                                                                        
      IO = MIN(MAX(0,IL),1)                                             
      IF(LUNIT.NE.LUNIN) IO = 0                                         
                                                                        
      CALL STRING(STR,LUN,I1,IO)                                        
                                                                        
C  INITIALIZE USR ARRAY PRECEEDING AN INPUT OPERATION                   
C  --------------------------------------------------                   
                                                                        
      IF(IO.EQ.0) THEN                                                  
         DO I=1,I1*I2                                                   
         USR(I,1) = BMISS                                               
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  CALL THE MNEMONIC READER/WRITER                                      
C  -------------------------------                                      
                                                                        
      CALL UFBRW(LUN,USR,I1,I2,IO,IRET)                                 
                                                                        
C  IF INCOMPLETE WRITE TRY TO INITIALIZE REPLICATION SEQUENCE OR RETURN 
C  ---------------------------------------------------------------------
                                                                        
      IF(IO.EQ.1 .AND. IRET.NE.I2 .AND. IRET.GE.0) THEN                 
         CALL TRYBUMP(LUNIT,LUN,USR,I1,I2,IO,IRET)                      
         IF(IRET.NE.I2) PRINT*,STR                                      
         IF(IRET.NE.I2) GOTO 903                                        
      ELSEIF(IRET.EQ.-1) THEN                                           
         IRET = 0                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('UFBINT - FILE IS CLOSED                     ')        
901   CALL BORT('UFBINT - NO MESSAGE OPEN                    ')        
902   CALL BORT('UFBINT - I-NODE MISMATCH                    ')        
903   CALL BORT('UFBINT - INCOMPLETE WRITE                   ')        
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBMEM(LUNIT,INEW,IRET,IUNIT)                          
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 SEC0,CUNIT                                            
      DIMENSION   MBAY(5000)                                            
      EQUIVALENCE (SEC0,MBAY(1))                                        
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      WRITE(CUNIT,'(I4)') LUNIT                                         
                                                                        
C  TRY TO OPEN BUFR FILE AND SET TO INITIALIZE OR CONCATENATE           
C  ----------------------------------------------------------           
                                                                        
      CALL OPENBF(LUNIT,'IN',LUNIT)                                     
                                                                        
      IF(INEW.EQ.0) THEN                                                
         MSGP(0) = 0                                                    
         MUNIT = 0                                                      
         MLAST = 0                                                      
      ENDIF                                                             
                                                                        
      IMSG = 8/NBYTW+1                                                  
      NMSG = MSGP(0)                                                    
      IRET = 0                                                          
                                                                        
C  TRANSFER MESSAGES FROM FILE TO MEMORY - SET MESSAGE POINTERS         
C  ------------------------------------------------------------         
                                                                        
1     READ(LUNIT,ERR=900,END=100) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))      
      IRET = IRET+1                                                     
      NMSG = NMSG+1                                                     
      LMEM = LMSG(SEC0)+IMSG-1                                          
ccccc print *, '^^^^^ NMSG,MAXMSG = ',nmsg,maxmsg
      IF(NMSG      .GT.MAXMSG) CALL BORT('UFBMEM - MSG BUFFER')        
ccccc print *, '^^^^^ LMEM+MLAST,MAXMEM = ',lmem+mlast,maxmem
      IF(LMEM+MLAST.GT.MAXMEM) CALL BORT('UFBMEM - MEM BUFFER')        
                                                                        
      DO I=1,LMEM                                                       
      MSGS(MLAST+I) = MBAY(I)                                           
      ENDDO                                                             
                                                                        
      MSGP(0000) = NMSG                                                 
      MSGP(NMSG) = MLAST+1                                              
      MLAST = MLAST+LMEM                                                
      GOTO 1                                                            
                                                                        
C  EXITS                                                                
C  -----                                                                
                                                                        
100   IF(IRET.EQ.0) THEN                                                
         CALL CLOSBF(LUNIT)                                             
      ELSE                                                              
         IF(MUNIT.NE.0) CALL CLOSBF(LUNIT)                              
         IF(MUNIT.EQ.0) MUNIT = LUNIT                                   
      ENDIF                                                             
      IUNIT = MUNIT                                                     
      RETURN                                                            
                                                                        
900   CALL BORT('UFBMEM - ERROR READING MESSAGE ON UNIT:'//CUNIT)      
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBMMS(IMSG,ISUB,SUBSET,IDATE)                         
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
                                                                        
      CHARACTER*8   SUBSET                                              
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  UFBINT SUBSET #ISUB FROM MEMORY MESSAGE #IMSG                        
C  ---------------------------------------------                        
                                                                        
      CALL RDMEMM(IMSG,SUBSET,IDATE,IRET)                               
      IF(IRET.NE.0) GOTO 900                                            
      CALL RDMEMS(ISUB,IRET)                                            
      IF(IRET.NE.0) GOTO 901                                            
                                                                        
      RETURN                                                            
900   CALL BORT('UFBMMS - BAD RETURN FROM RDMEMM')                     
901   CALL BORT('UFBMMS - BAD RETURN FROM RDMEMS')                     
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBMNS(IREP,SUBSET,IDATE)                              
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
                                                                        
      CHARACTER*8   SUBSET                                              
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
      JREP = 0                                                          
      IMSG = 1                                                          
                                                                        
C  READ SUBSET #ISUB FROM MEMORY MESSAGE #IMSG                          
C  -------------------------------------------                          
                                                                        
      DO WHILE(IRET.EQ.0)                                               
      CALL RDMEMM(IMSG,SUBSET,IDATE,IRET)                               
      IF(IRET.NE.0) GOTO 900                                            
      IF(JREP+NMSUB(MUNIT).GE.IREP) THEN                                
         CALL RDMEMS(IREP-JREP,IRET)                                    
         IF(IRET.NE.0) GOTO 901                                         
         RETURN                                                         
      ELSE                                                              
         JREP = JREP+NMSUB(MUNIT)                                       
         IMSG = IMSG+1                                                  
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      CALL BORT('UFBMNS - REPORT NUMBER OUT OF RANGE')                 
900   CALL BORT('UFBMNS - BAD RETURN FROM RDMEMM')                     
901   CALL BORT('UFBMNS - BAD RETURN FROM RDMEMS')                     
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBOVR(LUNIT,USR,I1,I2,IRET,STR)                       
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*(*) STR                                                 
      DIMENSION     USR(I1,I2)                                          
      REAL*8        USR,VAL                                             
                                                                        
C---------------------------------------------------------------------- 
CFPP$ EXPAND (STATUS,UFBRW)                                             
C---------------------------------------------------------------------- 
                                                                        
      IRET = 0                                                          
      IF(I2.LE.0) RETURN                                                
                                                                        
C  CHECK THE FILE STATUS AND I-NODE                                     
C  --------------------------------                                     
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.NE.1) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
                                                                        
      IO = MIN(MAX(0,IL),1)                                             
                                                                        
C  PARSE THE INPUT STRING - READ/WRITE VALUES                           
C  ------------------------------------------                           
                                                                        
      CALL STRING(STR,LUN,I1,IO)                                        
      CALL TRYBUMP(LUNIT,LUN,USR,I1,I2,IO,IRET)                         
                                                                        
      IF(IO.EQ.1 .AND. IRET.NE.I2) THEN                                 
         IF(IRET.NE.I2) PRINT*,STR                                      
         IF(IRET.NE.I2) GOTO 903                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('UFBOVR - FILE IS NOT OPEN FOR OUTPUT        ')        
901   CALL BORT('UFBOVR - NO MESSAGE OPEN                    ')        
902   CALL BORT('UFBOVR - I-NODE MISMATCH                    ')        
903   CALL BORT('UFBOVR - INCOMPLETE WRITE                   ')        
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBQCD(LUNIT,NEMO,QCD)                                 
                                                                        
      CHARACTER*(*) NEMO                                                
      CHARACTER*6  FXY,ADN30                                            
      CHARACTER*1  TAB                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
                                                                        
      CALL NEMTAB(LUN,NEMO,IDN,TAB,IRET)                                
      IF(TAB.NE.'D') GOTO 901                                           
                                                                        
      FXY = ADN30(IDN,6)                                                
      IF(FXY(2:3).NE.'63') GOTO 902                                     
      READ(FXY(4:6),'(F3.0)',ERR=903) QCD                               
                                                                        
      RETURN                                                            
900   CALL BORT('UFBQCD - FILE IS CLOSED                       ')      
901   CALL BORT('UFBQCD - MISSING OR INVALID TABLE D QC CODE   ')      
902   CALL BORT('UFBQCD - TABLE D QC CODE DESCRIPTOR NOT 363YYY')      
903   CALL BORT('UFBQCD - ERROR READING YYY FROM QC CODE DESCRP')      
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBQCP(LUNIT,QCP,NEMO)                                 
                                                                        
      CHARACTER*(*) NEMO                                                
      CHARACTER*1  TAB                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
                                                                        
      IDN = IFXY('363000')+IFIX(QCP)                                    
      CALL NUMTAB(LUN,IDN,NEMO,TAB,IRET)                                
                                                                        
      RETURN                                                            
900   CALL BORT('UFBQCP - FILE IS CLOSED                       ')      
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBREP(LUNIO,USR,I1,I2,IRET,STR)                       
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*(*) STR                                                 
      DIMENSION     USR(I1,I2)                                          
      REAL*8        USR,VAL,BMISS

      DATA BMISS /10E10/
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE FILE STATUS AND I-NODE                                     
C  --------------------------------                                     
                                                                        
      LUNIT = ABS(LUNIO)                                                
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
      IO = MIN(MAX(0,IL),1)                                             
      IF(LUNIO.NE.LUNIT) IO = 0                                         
                                                                        
C  INITIALIZE USR ARRAY PRECEEDING AN INPUT OPERATION                   
C  --------------------------------------------------                   
                                                                        
      IF(IO.EQ.0) THEN                                                  
         DO I=1,I1*I2                                                   
         USR(I,1) = BMISS                                               
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  PARSE THE INPUT STRING - READ/WRITE VALUES                           
C  ------------------------------------------                           
                                                                        
      CALL STRING(STR,LUN,I1,IO)                                        
      CALL UFBRP(LUN,USR,I1,I2,IO,IRET)                                 
                                                                        
      IF(IO.EQ.1 .AND. IRET.LT.I2) THEN                                 
         PRINT*,STR                                                     
         PRINT*,' UFBREP STORED ',IRET,' OF ',I2,' LEVELS '
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('UFBREP - FILE IS CLOSED                     ')        
901   CALL BORT('UFBREP - NO MESSAGE OPEN                    ')        
902   CALL BORT('UFBREP - I-NODE MISMATCH                    ')        
903   CALL BORT('UFBREP - INCOMPLETE WRITE                   ')        
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBRMS(IMSG,ISUB,USR,I1,I2,IRET,STR)                   
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*8   SUBSET                                              
      DIMENSION     USR(I1,I2)                                          
      REAL*8        USR                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  UFBINT SUBSET #ISUB FROM MEMORY MESSAGE #IMSG                        
C  ---------------------------------------------                        
                                                                        
      CALL RDMEMM(IMSG,SUBSET,IDATE,IRET)                               
      IF(IRET.NE.0) GOTO 900                                            
      CALL RDMEMS(ISUB,IRET)                                            
      IF(IRET.NE.0) GOTO 901                                            
                                                                        
      CALL UFBINT(MUNIT,USR,I1,I2,IRET,STR)                             
                                                                        
      RETURN                                                            
900   CALL BORT('UFBRMS - BAD RETURN FROM RDMEMM')                     
901   CALL BORT('UFBRMS - BAD RETURN FROM RDMEMS')                     
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBRP(LUN,USR,I1,I2,IO,IRET)                           
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      REAL*8       USR(I1,I2),VAL                                       
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IRET = 0                                                          
      INS1 = 0                                                          
      INS2 = 0                                                          

c  find first non-zero node in string
c  ----------------------------------

      do nz=1,nnod
      if(nods(nz).gt.0) goto 1
      enddo
      return
                                                                        
C  FRAME A SECTION OF THE BUFFER - RETURN WHEN NO FRAME                 
C  ----------------------------------------------------                 
                                                                        
1     IF(INS1+1.GT.NVAL(LUN)) RETURN
      IF(IO.EQ.1 .AND. IRET.EQ.I2) RETURN                               
      INS1 = INVtag(NODS(nz),LUN,INS1+1,NVAL(LUN))                       
      IF(INS1.EQ.0) RETURN                                              
                                                                        
      INS2 = INVtag(NODS(nz),LUN,INS1+1,NVAL(LUN))                       
      IF(INS2.EQ.0) INS2 = NVAL(LUN)                                    
      IRET = IRET+1                                                     
                                                                        
C  READ USER VALUES                                                     
C  ----------------                                                     
                                                                        
      IF(IO.EQ.0 .AND. IRET.LE.I2) THEN                                 
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            INVN = INVTAG(NODS(I),LUN,INS1,INS2)                        
            IF(INVN.GT.0) USR(I,IRET) = VAL(INVN,LUN)                   
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  WRITE USER VALUES                                                    
C  -----------------                                                    
                                                                        
      IF(IO.EQ.1 .AND. IRET.LE.I2) THEN                                 
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            INVN = INVtag(NODS(I),LUN,INS1,INS2)                        
            IF(INVN.GT.0) VAL(INVN,LUN) = USR(I,IRET)                   
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  GO FOR NEXT FRAME                                                    
C  -----------------                                                    
                                                                        
      GOTO 1                                                            
                                                                        
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBRW(LUN,USR,I1,I2,IO,IRET)                           
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      REAL*8       USR(I1,I2),VAL,BMISS                                 

      DATA BMISS/10E10/
                                                                        
C---------------------------------------------------------------------- 
CFPP$ EXPAND (CONWIN,DRSTPL,GETWIN,INVWIN,LSTRPS,NEWWIN,NXTWIN)          
C---------------------------------------------------------------------- 
                                                                        
      IRET = 0                                                          
                                                                        
C  LOOP OVER COND WINDOWS                                               
C  ----------------------                                               
                                                                        
      INC1 = 1                                                          
      INC2 = 1                                                          
                                                                        
1     CALL CONWIN(LUN,INC1,INC2,I2)                                     
      IF(NNOD.EQ.0) THEN                                                
         IRET = I2                                                      
         RETURN                                                         
      ELSEIF(INC1.EQ.0) THEN                                            
         RETURN                                                         
      ELSE                                                              
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            INS2 = INC1                                                 
            CALL GETWIN(NODS(I),LUN,INS1,INS2)                          
            IF(INS1.EQ.0) RETURN                                        
            GOTO 2                                                      
         ENDIF                                                          
         ENDDO                                                          
         IRET = -1                                                      
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  LOOP OVER STORE NODES                                                
C  ---------------------                                                
                                                                        
2     IRET = IRET+1                                                     
                                                                        
C     print*,'ufbrw:',iret,':',ins1,':',ins2,':',inc1,':',inc2          
C     print'(5a10)',(tag(inv(i,lun)),i=ins1,ins2)                       
                                                                        
C  WRITE USER VALUES                                                    
C  -----------------                                                    
                                                                        
      IF(IO.EQ.1 .AND. IRET.LE.I2) THEN                                 
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            IF(USR(I,IRET).NE.BMISS) THEN                               
               INVN = INVWIN(NODS(I),LUN,INS1,INS2)                     
               IF(INVN.EQ.0) THEN                                       
                  CALL DRSTPL(NODS(I),LUN,INS1,INS2,INVN)               
                  if(invn.eq.0) then                                    
                     iret = 0                                           
                     return                                             
                  endif                                                 
                  CALL NEWWIN(LUN,INC1,INC2)                            
                  VAL(INVN,LUN) = USR(I,IRET)                           
               ELSEIF(LSTRPS(NODS(I),LUN).EQ.0) THEN                    
                  VAL(INVN,LUN) = USR(I,IRET)                           
               ELSEIF(VAL(INVN,LUN).EQ.BMISS) THEN                      
                  VAL(INVN,LUN) = USR(I,IRET)                           
               ELSE                                                     
                  CALL DRSTPL(NODS(I),LUN,INS1,INS2,INVN)               
                  if(invn.eq.0) then                                    
                     iret = 0                                           
                     return                                             
                  endif                                                 
                  CALL NEWWIN(LUN,INC1,INC2)                            
                  VAL(INVN,LUN) = USR(I,IRET)                           
               ENDIF                                                    
            ENDIF                                                       
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  READ USER VALUES                                                     
C  ----------------                                                     
                                                                        
      IF(IO.EQ.0 .AND. IRET.LE.I2) THEN                                 
         DO I=1,NNOD                                                    
         USR(I,IRET) = BMISS                                            
         IF(NODS(I).GT.0) THEN                                          
            INVN = INVWIN(NODS(I),LUN,INS1,INS2)                        
            IF(INVN.GT.0) USR(I,IRET) = VAL(INVN,LUN)                   
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  DECIDE WHAT TO DO NEXT                                               
C  ----------------------                                               
                                                                        
      IF(IO.EQ.1.AND.IRET.EQ.I2) RETURN                                 
      CALL NXTWIN(LUN,INS1,INS2)                                        
      IF(INS1.GT.0 .AND. INS1.LT.INC2) GOTO 2                           
      IF(NCON.GT.0) GOTO 1                                              
                                                                        
      RETURN                                                            
      END                                                               
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE UFBSEQ(LUNIN,USR,I1,I2,IRET,STR)

      PARAMETER(MTAG=10)

      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),
     .                JUMP(15000),LINK(15000),JMPB(15000),
     .                IBT(15000),IRF(15000),ISC(15000),
     .                ITP(15000),VALI(15000),KNTI(15000),
     .                ISEQ(15000,2),JSEQ(15000)
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)

      CHARACTER*(*) STR
      CHARACTER*10  TAG,TAGS(MTAG)
      CHARACTER*3   TYP
      REAL*8        USR(I1,I2),VAL,BMISS

      DATA BMISS /10E10/

C----------------------------------------------------------------------
C----------------------------------------------------------------------

      IRET = 0

C  CHECK THE FILE STATUS AND CLEAR AN INPUT ARRAY
C  ----------------------------------------------

      LUNIT = ABS(LUNIN)
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) CALL BORT('UFBSEQ - FILE IS CLOSED ')
      IF(IM.EQ.0) CALL BORT('UFBSEQ - NO MESSAGE OPEN')

      IO = MIN(MAX(0,IL),1)
      IF(LUNIT.NE.LUNIN) IO = 0

      IF(IO.EQ.0) THEN
         DO I=1,I1*I2
         USR(I,1) = BMISS
         ENDDO
      ENDIF

C  CHECK FOR VALID SEQUENCE AND SEQUENCE LENGTH ARGUMENTS
C  ------------------------------------------------------

      CALL PARSEQ(STR,TAGS,MTAG,NTAG)
      IF(NTAG.LT.1) RETURN
      IF(NTAG.GT.1) CALL BORT('UFBSEQ - MORE THAN ONE STRING ARG!')

      INOD = INODE(LUN)
      IF(INOD.NE.INV(1,LUN)) CALL BORT('UFBSEQ - I-NODE MISMATCH!')

      DO NODE=INOD,ISC(INOD)
      IF(STR.EQ.TAG(NODE)) THEN
         IF(TYP(NODE).EQ.'SEQ') THEN
            INS1 = INVTAG(NODE,LUN,1,NVAL(LUN))
            NODS = NODE
            DO WHILE(LINK(NODS).EQ.0.AND.JMPB(NODS).GT.0)
            NODS = JMPB(NODS)
            ENDDO
            IF(LINK(NODS).EQ.0) THEN
               INS2 = NVAL(LUN)
            ELSEIF(LINK(NODS).GT.0) THEN
               INS2 = INVWIN(LINK(NODS)-1,LUN,INS1+1,NVAL(LUN))
            ENDIF
         ELSEIF(TYP(NODE).EQ.'SUB') THEN
            INS1 = 1             
            INS2 = NVAL(LUN)     
         ELSE
            CALL BORT('UFBSEQ - MNEMONIC NOT A SEQUENCE!')
         ENDIF
         NSEQ = 0
         DO ISQ=INS1,INS2
         ITYP = ITP(INV(ISQ,LUN))
         IF(ITYP.EQ.1) CALL BORT('UFBSEQ - ILLEGAL SEQUENCE!')
         IF(ITYP.GT.1) NSEQ = NSEQ+1 
         ENDDO
         IF(NSEQ.GT.I1) CALL BORT('UFBSEQ - SEQ LENGTH GT I1!')
         GOTO 1
      ENDIF
      ENDDO

      RETURN

C  FRAME A SECTION OF THE BUFFER - RETURN WHEN NO FRAME
C  ----------------------------------------------------

1     INS1 = INVTAG(NODE,LUN,INS1,NVAL(LUN))
      IF(INS1.EQ.0.AND.IO.EQ.1.AND.IRET.LT.I2) THEN
         CALL BORT('UFBSEQ - INCOMPLETE WRITE')
      ELSEIF(INS1.GT.0.AND.IO.EQ.0.AND.IRET+1.GT.I2) THEN
         CALL BORT('UFBSEQ - INCOMPLETE READ ')
      ELSEIF(INS1.EQ.0.OR.IRET.EQ.I2) THEN
         RETURN
      ENDIF

      IRET = IRET+1
      INS1 = INS1+1

C  READ/WRITE USER VALUES
C  ----------------------

      J = INS1
      DO I=1,NSEQ
      DO WHILE(ITP(INV(J,LUN)).LT.2) 
      J = J+1
      ENDDO
      IF(IO.EQ.0) USR(I,IRET) = VAL(J,LUN )
      IF(IO.EQ.1) VAL(J,LUN ) = USR(I,IRET)
      J = J+1
      ENDDO

C  CHECK FOR NEXT FRAME
C  --------------------

      GOTO 1
      END
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBSP(LUN,USR,I1,I2,IO,IRET)                           
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      DIMENSION    USR(I1,I2)
      REAL*8       USR,VAL,BMISS                                 

      DATA BMISS /10E10/
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IRET = 0                                                          
      INS1 = 0                                                          
      INS2 = 0                                                          
                                                                        
C  FRAME A SECTION OF THE BUFFER - RETURN WHEN NO FRAME                 
C  ----------------------------------------------------                 
                                                                        
1     IF(INS1+1.GT.NVAL(LUN)) RETURN                                    
      INS1 = INVtag(NODS(1),LUN,INS1+1,NVAL(LUN))                       
      IF(INS1.EQ.0) RETURN                                              
                                                                        
      INS2 = INVtag(NODS(1),LUN,INS1+1,NVAL(LUN))                       
      IF(INS2.EQ.0) INS2 = NVAL(LUN)                                    
      IRET = IRET+1                                                     
                                                                        
C  READ USER VALUES                                                     
C  ----------------                                                     
                                                                        
      IF(IO.EQ.0 .AND. IRET.LE.I2) THEN                                 
         INVM = INS1
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            INVN = INVTAG(NODS(I),LUN,INVM,INS2)                        
            IF(INVN.GT.0) USR(I,IRET) = VAL(INVN,LUN)                   
            INVM = MAX(INVN,INVM)
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  WRITE USER VALUES                                                    
C  -----------------                                                    
                                                                        
      IF(IO.EQ.1 .AND. IRET.LE.I2) THEN                                 
         INVM = INS1
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            INVN = INVtag(NODS(I),LUN,INVM,INS2)                        
            IF(INVN.GT.0) VAL(INVN,LUN) = USR(I,IRET)                   
            INVM = MAX(INVN,INVM)
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  GO FOR NEXT FRAME                                                    
C  -----------------                                                    
                                                                        
      GOTO 1                                                            
                                                                        
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBSTP(LUNIO,USR,I1,I2,IRET,STR)                       
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*(*) STR                                                 
      DIMENSION     USR(I1,I2)                                          
      REAL*8        USR,VAL,BMISS                                 

      DATA BMISS /10E10/
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE FILE STATUS AND I-NODE                                     
C  --------------------------------                                     
                                                                        
      LUNIT = ABS(LUNIO)                                                
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
      IO = MIN(MAX(0,IL),1)                                             
      IF(LUNIO.NE.LUNIT) IO = 0                                         

C  INITIALIZE USR ARRAY PRECEEDING AN INPUT OPERATION                   
C  --------------------------------------------------                   
                                                                        
      IF(IO.EQ.0) THEN                                                  
         DO I=1,I1*I2                                                   
         USR(I,1) = BMISS                                               
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  PARSE THE INPUT STRING - READ/WRITE VALUES                           
C  ------------------------------------------                           
                                                                        
      CALL STRING(STR,LUN,I1,IO)                                        
      CALL UFBSP(LUN,USR,I1,I2,IO,IRET)                                 
                                                                        
      IF(IO.EQ.1 .AND. IRET.NE.I2) THEN                                 
         PRINT*,STR,' i2=',i2,' iret=',iret                             
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('UFBSTP - FILE IS CLOSED                     ')        
901   CALL BORT('UFBSTP - NO MESSAGE OPEN                    ')        
902   CALL BORT('UFBSTP - I-NODE MISMATCH                    ')        
903   CALL BORT('UFBSTP - INCOMPLETE WRITE                   ')        
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBTAB(LUNIT,TAB,I1,I2,IRET,STR)                       
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      common /acmode/ iac
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*10  TAG,TGS(100)                                        
      CHARACTER*8   SUBSET,CVAL                                         
      CHARACTER*3   TYP                                                 
      DIMENSION     TAB(I1,I2)                                          
      EQUIVALENCE   (CVAL,RVAL)                                         
      LOGICAL       OPENIT                                              
      REAL*8        VAL,TAB,RVAL,BMISS
                                                                        
      DATA BMISS /10E10/                                                
      DATA MAXTG /100/                                                  
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND (UPBB,USRTPL)                                              
C-----------------------------------------------------------------------
      MPS(NODE) = 2**(IBT(NODE))-1                                      
      UPS(NODE) = (IVAL+IRF(NODE))*10.**(-ISC(NODE))                    
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
      IREC = 0                                                          
      ISUB = 0                                                          
      iacc = iac
                                                                        
      DO I=1,I1*I2                                                      
      TAB(I,1) = BMISS                                                  
      ENDDO                                                             
                                                                        
C  SEE IF WE NEED TO OPEN A FILE                                        
C  -----------------------------                                        
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      OPENIT = IL.EQ.0                                                  
      CALL OPENBF(LUNIT,'IN',LUNIT)                                  
                                                                        
      iac = 1

C  CHECK FOR SPECIAL TAGS IN STRING                                     
C  --------------------------------                                     
                                                                        
      CALL PARSEQ(STR,TGS,MAXTG,NTG)                                    
      DO I=1,NTG                                                        
      IF(TGS(I).EQ.'IREC') IREC = I                                     
      IF(TGS(I).EQ.'ISUB') ISUB = I                                     
      ENDDO                                                             
                                                                        
C  READ A MESSAGE AND PARSE A STRING                                    
C  ---------------------------------                                    
                                                                        
10    CALL READMG(LUNIT,SUBSET,JDATE,MRET)                              
      IF(MRET.NE.0) GOTO 25                                             
      CALL STRING(STR,LUN,I1,0)                                         
      IF(IREC.GT.0) NODS(IREC) = 0                                      
      IF(ISUB.GT.0) NODS(ISUB) = 0                                      
                                                                        
C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE                        
C  ---------------------------------------------                        
                                                                        
15    IF(NSUB(LUN).EQ.MSUB(LUN)) GOTO 10                                
      IF(IRET+1.GT.I2) CALL BORT('UFBTAB - TAB TOO SMALL')             
      IRET = IRET+1                                                     
                                                                        
      DO I=1,NNOD                                                       
      NODS(I) = ABS(NODS(I))                                            
      ENDDO                                                             
                                                                        
C  PARSE THE STRING NODES FROM A SUBSET                                 
C  ------------------------------------                                 
                                                                        
      MBIT = MBYT(LUN)*8 + 16                                           
      NBIT = 0                                                          
      N = 1                                                             
      CALL USRTPL(LUN,N,N)                                              
20    IF(N+1.LE.NVAL(LUN)) THEN                                         
         N = N+1                                                        
         NODE = INV(N,LUN)                                              
         MBIT = MBIT+NBIT                                               
         NBIT = IBT(NODE)                                               
         IF(ITP(NODE).EQ.1) THEN                                        
            CALL UPBB(IVAL,NBIT,MBIT,MBAY(1,LUN))                       
            CALL USRTPL(LUN,N,IVAL)                                     
         ENDIF                                                          
         DO I=1,NNOD                                                    
         IF(NODS(I).EQ.NODE) THEN                                       
            IF(ITP(NODE).EQ.1) THEN                                     
               CALL UPBB(IVAL,NBIT,MBIT,MBAY(1,LUN))                    
               TAB(I,IRET) = IVAL                                       
            ELSEIF(ITP(NODE).EQ.2) THEN                                 
               CALL UPBB(IVAL,NBIT,MBIT,MBAY(1,LUN))                    
               IF(IVAL.LT.MPS(NODE)) TAB(I,IRET) = UPS(NODE)            
            ELSEIF(ITP(NODE).EQ.3) THEN                                 
               CVAL = ' '                                               
               KBIT = MBIT                                              
               CALL UPC(CVAL,NBIT/8,MBAY(1,LUN),KBIT)                   
               TAB(I,IRET) = RVAL                                       
            ENDIF                                                       
            NODS(I) = -NODS(I)                                          
            GOTO 20                                                     
         ENDIF                                                          
         ENDDO                                                          
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) GOTO 20                                       
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  UPDATE THE SUBSET POINTERS BEFORE NEXT READ                          
C  -------------------------------------------                          
                                                                        
      IBIT = MBYT(LUN)*8                                                
      CALL UPB(NBYT,16,MBAY(1,LUN),IBIT)                                
      MBYT(LUN) = MBYT(LUN) + NBYT                                      
      NSUB(LUN) = NSUB(LUN) + 1                                         
      IF(IREC.GT.0) TAB(IREC,IRET) = NMSG(LUN)                          
      IF(ISUB.GT.0) TAB(ISUB,IRET) = NSUB(LUN)                          
      GOTO 15                                                           
                                                                        
C  LEAVE THE FILE AS IT WAS BEFORE                                      
C  -------------------------------                                      
                                                                        
25    CALL CLOSBF(LUNIT)                                             
      iac = iacc

      RETURN                                                            
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBTAM(TAB,I1,I2,IRET,STR)                             
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),VALS(10),KONS(10)     
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*10  TAG,TGS(100)                                        
      CHARACTER*8   SUBSET,CVAL                                         
      CHARACTER*3   TYP                                                 
      DIMENSION     TAB(I1,I2)                                          
      EQUIVALENCE   (CVAL,RVAL)                                         
      REAL*8        TAB,VAL,RVAL,BMISS
                                                                        
      DATA BMISS /10E10/
      DATA MAXTG /100/                                                  
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND (UPBB,UPC,USRTPL)                                          
C-----------------------------------------------------------------------
      MPS(NODE) = 2**(IBT(NODE))-1                                      
      UPS(NODE) = (IVAL+IRF(NODE))*10.**(-ISC(NODE))                    
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
                                                                        
      IF(MSGP(0).EQ.0) RETURN                                           
                                                                        
      DO I=1,I1*I2                                                      
      TAB(I,1) = 10E10                                                  
      ENDDO                                                             
                                                                        
C  CHECK FOR SPECIAL TAGS IN STRING                                     
C  --------------------------------                                     
                                                                        
      CALL PARSEQ(STR,TGS,MAXTG,NTG)                                    
      IREC = 0                                                          
      ISUB = 0                                                          
      DO I=1,NTG                                                        
      IF(TGS(I).EQ.'IREC') IREC = I                                     
      IF(TGS(I).EQ.'ISUB') ISUB = I                                     
      ENDDO                                                             
                                                                        
C  READ A MESSAGE AND PARSE A STRING                                    
C  ---------------------------------                                    
                                                                        
      CALL STATUS(MUNIT,LUN,IL,IM)                                      
                                                                        
      DO IMSG=1,MSGP(0)                                                 
      CALL RDMEMM(IMSG,SUBSET,JDATE,MRET)                               
      IF(MRET.NE.0) GOTO 900                                            
                                                                        
      CALL STRING(STR,LUN,I1,0)                                         
      IF(IREC.GT.0) NODS(IREC) = 0                                      
      IF(ISUB.GT.0) NODS(ISUB) = 0                                      
                                                                        
C  PROCESS ALL THE SUBSETS IN THE MEMORY MESSAGE                        
C  ---------------------------------------------                        
                                                                        
      DO WHILE (NSUB(LUN).LT.MSUB(LUN))                                 
         IF(IRET+1.GT.I2) GOTO 99
         IRET = IRET+1                                                  
                                                                        
         DO I=1,NNOD                                                    
         NODS(I) = ABS(NODS(I))                                         
         ENDDO                                                          
                                                                        
         CALL USRTPL(LUN,1,1)                                           
         MBIT = MBYT(LUN)*8+16                                          
         NBIT = 0                                                       
         N = 1                                                          
                                                                        
20       IF(N+1.LE.NVAL(LUN)) THEN                                      
            N = N+1                                                     
            NODE = INV(N,LUN)                                           
            MBIT = MBIT+NBIT                                            
            NBIT = IBT(NODE)                                            
            IF(ITP(NODE).EQ.1) THEN                                     
               CALL UPBB(IVAL,NBIT,MBIT,MBAY(1,LUN))                    
               CALL USRTPL(LUN,N,IVAL)                                  
            ENDIF                                                       
            DO I=1,NNOD                                                 
            IF(NODS(I).EQ.NODE) THEN                                    
               IF(ITP(NODE).EQ.1) THEN                                  
                  CALL UPBB(IVAL,NBIT,MBIT,MBAY(1,LUN))                 
                  TAB(I,IRET) = IVAL                                    
               ELSEIF(ITP(NODE).EQ.2) THEN                              
                  CALL UPBB(IVAL,NBIT,MBIT,MBAY(1,LUN))                 
                  IF(IVAL.LT.MPS(NODE)) TAB(I,IRET) = UPS(NODE)         
               ELSEIF(ITP(NODE).EQ.3) THEN                              
                  CVAL = ' '                                            
                  KBIT = MBIT                                           
                  CALL UPC(CVAL,NBIT/8,MBAY(1,LUN),KBIT)                
                  TAB(I,IRET) = RVAL                                    
               ENDIF                                                    
               NODS(I) = -NODS(I)                                       
               GOTO 20                                                  
            ENDIF                                                       
            ENDDO                                                       
            DO I=1,NNOD                                                 
            IF(NODS(I).GT.0) GOTO 20                                    
            ENDDO                                                       
         ENDIF                                                          
                                                                        
C  UPDATE THE SUBSET POINTERS BEFORE NEXT READ                          
C  -------------------------------------------                          
                                                                        
         IBIT = MBYT(LUN)*8                                             
         CALL UPB(NBYT,16,MBAY(1,LUN),IBIT)                             
         MBYT(LUN) = MBYT(LUN) + NBYT                                   
         NSUB(LUN) = NSUB(LUN) + 1                                      
         IF(IREC.GT.0) TAB(IREC,IRET) = NMSG(LUN)                       
         IF(ISUB.GT.0) TAB(ISUB,IRET) = NSUB(LUN)                       
      ENDDO                                                             
                                                                        
      ENDDO                                                             
                                                                        
C  RESET THE MEMORY FILE AND EXIT                                       
C  ------------------------------                                       
                                                                        
      CALL RDMEMM(0,SUBSET,JDATE,MRET)                                  
      RETURN                                                            

C  VARIOUS ERROR EXITS
C  -------------------
                                                                        
99    CALL RDMEMM(0,SUBSET,JDATE,MRET)                                  
      NREP = 0
      DO IMSG=1,MSGP(0)                                                 
      CALL RDMEMM(IMSG,SUBSET,JDATE,MRET)                               
      IF(MRET.NE.0) GOTO 900                                            
      NREP = NREP+NMSUB(MUNIT)
      ENDDO
      PRINT*
      PRINT*,'>>>UFBTAM STORED ',IRET,' REPORTS OUT OF ',NREP,'<<<'
      PRINT*
      CALL RDMEMM(0,SUBSET,JDATE,MRET)                                  
      RETURN                                                            

900   CALL BORT('UFBTAM - EOF ON MEMORY MESSAGES')                     
      END                                                               
C----------------------------------------------------------------------
C  UNPACK UP AN INTEGER
C----------------------------------------------------------------------
      SUBROUTINE UPB(NVAL,NBITS,IBAY,IBIT)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      DIMENSION IBAY(*)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
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
C----------------------------------------------------------------------
C  COPY CHARACTERS FROM A BIT ARRAY
C----------------------------------------------------------------------
      SUBROUTINE UPC(CHR,NCHR,IBAY,IBIT)
 
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
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UPTDD(ID,LUN,IENT,IRET)                                
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*56  DXSTR                                               
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      LDD = LDXD(IDXV+1)+1                                              
                                                                        
C  CHECK IF IENT IS IN BOUNDS                                           
C  --------------------------                                           
                                                                        
      NDSC = IUPM(TABD(ID,LUN)(LDD:LDD),8)                              
                                                                        
      IF(IENT.EQ.0) THEN                                                
         IRET = NDSC                                                    
         RETURN                                                         
      ELSEIF(IENT.LT.0 .OR. IENT.GT.NDSC) THEN                          
         CALL BORT('UPTDD - IENT OUT OF RANGE')                        
      ENDIF                                                             
                                                                        
C  RETURN THE DESCRIPTOR INDICATED BY IENT                              
C  ---------------------------------------                              
                                                                        
      IDSC = LDD+1 + (IENT-1)*2                                         
      IRET = IUPM(TABD(ID,LUN)(IDSC:IDSC),16)                           
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE USRTPL(LUN,INVN,NBMP)                                  
                                                                        
      PARAMETER (MAXINV=15000)                                          
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      DIMENSION    ITMP(MAXINV),VTMP(MAXINV)                            
      LOGICAL      DRP,DRS,DRB,DRX                                      
      REAL*8       VAL,VTMP                                             
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     PRINT*,'USRTPL:',LUN,':',INVN,':',NBMP,':',tag(inode(lun))        
                                                                        
      IF(NBMP.LE.0) RETURN                                              
                                                                        
      DRP = .FALSE.                                                     
      DRS = .FALSE.                                                     
      DRX = .FALSE.                                                     
                                                                        
C  SET UP A NODE EXPANSION                                              
C  -----------------------                                              
                                                                        
      IF(INVN.EQ.1) THEN                                                
         NODI = INODE(LUN)                                              
         INV(1,LUN) = NODI                                              
         NVAL(LUN)  = 1                                                 
         IF(NBMP.NE.1) GOTO 900                                         
      ELSEIF(INVN.GT.0 .AND. INVN.LE.NVAL(LUN)) THEN                    
         NODI = INV(INVN,LUN)                                           
         DRP  = TYP(NODI) .EQ. 'DRP'                                    
         DRS  = TYP(NODI) .EQ. 'DRS'                                    
         DRB  = TYP(NODI) .EQ. 'DRB'                                    
         DRX  = DRP .OR. DRS .OR. DRB                                   
         IVAL = VAL(INVN,LUN)                                           
         JVAL = 2**IBT(NODI)-1                                          
         VAL(INVN,LUN) = IVAL+NBMP                                      
         IF(DRB.AND.NBMP.NE.1) GOTO 900                                 
         IF(.NOT.DRX         ) GOTO 901                                 
         IF(IVAL.LT.0.       ) GOTO 902                                 
         IF(IVAL+NBMP.GT.JVAL) GOTO 903                                 
      ELSE                                                              
         GOTO 904                                                       
      ENDIF                                                             
                                                                        
C  RECALL A PRE-FAB NODE EXPANSION SEGMENT                              
C  ---------------------------------------                              
                                                                        
      NEWN = 0                                                          
      N1 = ISEQ(NODI,1)                                                 
      N2 = ISEQ(NODI,2)                                                 
                                                                        
      IF(N1.EQ.0          ) GOTO 905                                    
      IF(N2-N1+1.GT.MAXINV) GOTO 906                                    
                                                                        
      DO N=N1,N2                                                        
      NEWN = NEWN+1                                                     
      ITMP(NEWN) = JSEQ(N)                                              
      VTMP(NEWN) = VALI(JSEQ(N))                                        
      if(vtmp(newn).gt.10e9) vtmp(newn) = 10e10                         
      ENDDO                                                             
                                                                        
C  MOVE OLD NODES - STORE NEW ONES                                      
C  -------------------------------                                      
                                                                        
      IF(NVAL(LUN)+NEWN*NBMP.GT.MAXINV) print*,'@:',nval(lun)+newn*nbmp 
      IF(NVAL(LUN)+NEWN*NBMP.GT.MAXINV) GOTO 907                        
                                                                        
CDIR$ IVDEP                                                             
      DO J=NVAL(LUN),INVN+1,-1                                          
      INV(J+NEWN*NBMP,LUN) = INV(J,LUN)                                 
      VAL(J+NEWN*NBMP,LUN) = VAL(J,LUN)                                 
      ENDDO                                                             
                                                                        
      IF(DRP.OR.DRS) VTMP(1) = NEWN                                     
      KNVN = INVN                                                       
                                                                        
      DO I=1,NBMP                                                       
      DO J=1,NEWN                                                       
      KNVN = KNVN+1                                                     
      INV(KNVN,LUN) = ITMP(J)                                           
      VAL(KNVN,LUN) = VTMP(J)                                           
      ENDDO                                                             
      ENDDO                                                             
                                                                        
C  RESET POINTERS AND COUNTERS                                          
C  ---------------------------                                          
                                                                        
      NVAL(LUN) = NVAL(LUN) + NEWN*NBMP                                 
                                                                        
C     print*,tag(inv(invn,lun)),' ',newn,' ',nbmp,' ',nval(lun)         
C     DO I=1,NEWN                                                       
C     PRINT*,TAG(ITMP(I))                                               
C     ENDDO                                                             
                                                                        
                                                                        
      IF(DRX) THEN                                                      
         NODE = NODI                                                    
         INVR = INVN                                                    
4        NODE = JMPB(NODE)                                              
         IF(NODE.GT.0) THEN                                             
            IF(ITP(NODE).EQ.0) THEN                                     
               DO INVR=INVR-1,1,-1                                      
               IF(INV(INVR,LUN).EQ.NODE) THEN                           
                  VAL(INVR,LUN) = VAL(INVR,LUN)+NEWN*NBMP               
                  GOTO 4                                                
               ENDIF                                                    
               ENDDO                                                    
               GOTO 909                                                 
            ELSE                                                        
               GOTO 4                                                   
            ENDIF                                                       
         ENDIF                                                          
      ENDIF                                                             
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('USRTPL - NBMP <> 1 FOR        : '//TAG(NODI))         
901   CALL BORT('USRTPL - NODE NOT SUB,DRP,DRS : '//TAG(NODI))         
902   CALL BORT('USRTPL - NEGATIVE REP FACTOR  : '//TAG(NODI))         
903   CALL BORT('USRTPL - REP FACTOR OVERFLOW  : '//TAG(NODI))         
904   CALL BORT('USRTPL - INVENTORY INDEX OUT OF BOUNDS     ')         
905   CALL BORT('USRTPL - UNSET EXPANSION SEG  : '//TAG(NODI))         
906   CALL BORT('USRTPL - TEMP ARRAY OVERFLOW  : '//TAG(NODI))         
907   CALL BORT('USRTPL - INVENTORY OVERFLOW   : '//TAG(NODI))         
908   CALL BORT('USRTPL - TPL CACHE OVERFLOW   : '//TAG(NODI))         
909   CALL BORT('USRTPL - BAD BACKUP STRATEGY  : '//TAG(NODI))         
      END                                                               
C---------------------------------------------------------------------- 
C  REAL NUMBER FROM A STRING                                            
C---------------------------------------------------------------------- 
      FUNCTION VALX(STR)                                                
                                                                        
      CHARACTER*(*) STR
      CHARACTER*99  BSTR
      CHARACTER*8   FMT
      REAL*8        BMISS                                           
                                                                        
      DATA BMISS /10E10/
      data noinline /0/                                                 
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      LENS = LEN(STR)
      IF(LENS.GT.99) CALL BORT('VALX - ARG TOO LONG')
      BSTR(1:LENS) = STR            
      RJ = RJUST(BSTR(1:LENS))
      WRITE(FMT,'(''(F'',I2,''.0)'')') LENS                             
      READ(BSTR,FMT,ERR=900) VAL                                         
      VALX = VAL                                                        
      RETURN                                                            
900   VALX = BMISS                                                      
      RETURN                                                            
      END                                                               
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE WRDLEN
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
      COMMON /CHARAC/ IASCII,IATOE(0:255),IETOA(0:255)
      COMMON /QUIET / IPRT
 
      CHARACTER*8 CINT,DINT
      EQUIVALENCE (CINT,INT)
      EQUIVALENCE (DINT,JNT)
      LOGICAL     PRINT
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      PRINT = NBYTW.EQ.0 .AND. IPRT.EQ.1
 
C  COUNT THE BITS IN A WORD - MAX 64 ALLOWED
C  -----------------------------------------
 
      INT = 1
      DO I=1,65
      INT = ISHFT(INT,1)
      IF(INT.EQ.0) GOTO 10
      ENDDO
10    IF(I.GE.65)       GOTO 900
      IF(MOD(I,8).NE.0) GOTO 901
      NBITW = I
      NBYTW = I/8
 
C  INDEX THE BYTE STORAGE ORDER -  HIGH BYTE TO LOW BYTE
C  -----------------------------------------------------
 
      JNT = 0
      DO I=1,NBYTW
      INT = ISHFT(1,(NBYTW-I)*8)
      DO J=1,NBYTW
      IF(CINT(J:J).NE.DINT(J:J)) GOTO 20
      ENDDO
20    IF(J.GT.NBYTW) GOTO 902
      IORD(I) = J
      ENDDO
 
C  SET THE NOREVERSE FLAG - 0=NOREVERSE;1=REVERSE
C  ----------------------------------------------
 
      NREV = 0
      DO I=1,NBYTW
      IF(IORD(I).NE.I) NREV = 1
      ENDDO
 
C  SETUP AN ASCII/EBCDIC TRANSALTOR AND DETERMINE WHICH IS NATIVE
C  --------------------------------------------------------------
 
 
      IF(IUPM('A',8).EQ. 65) THEN
         IASCII = 1
      ELSEIF(IUPM('A',8).EQ.193) THEN
         IASCII = 0
      ELSE
         CALL BORT('WRDLEN - CANT DETERMINE NATIVE LANGUAGE')
      ENDIF
 
      DO I=0,255
      IETOA(I) = 0
      IATOE(I) = 0
      ENDDO
 
      IETOA(  1) =   1
      IATOE(  1) =   1
      IETOA(  2) =   2
      IATOE(  2) =   2
      IETOA(  3) =   3
      IATOE(  3) =   3
      IETOA(  5) =   9
      IATOE(  9) =   5
      IETOA(  7) = 127
      IATOE(127) =   7
      IETOA( 11) =  11
      IATOE( 11) =  11
      IETOA( 12) =  12
      IATOE( 12) =  12
      IETOA( 13) =  13
      IATOE( 13) =  13
      IETOA( 14) =  14
      IATOE( 14) =  14
      IETOA( 15) =  15
      IATOE( 15) =  15
      IETOA( 16) =  16
      IATOE( 16) =  16
      IETOA( 17) =  17
      IATOE( 17) =  17
      IETOA( 18) =  18
      IATOE( 18) =  18
      IETOA( 19) =  19
      IATOE( 19) =  19
      IETOA( 22) =   8
      IATOE(  8) =  22
      IETOA( 24) =  24
      IATOE( 24) =  24
      IETOA( 25) =  25
      IATOE( 25) =  25
      IETOA( 29) =  29
      IATOE( 29) =  29
      IETOA( 31) =  31
      IATOE( 31) =  31
      IETOA( 34) =  28
      IATOE( 28) =  34
      IETOA( 37) =  10
      IATOE( 10) =  37
      IETOA( 38) =  23
      IATOE( 23) =  38
      IETOA( 39) =  27
      IATOE( 27) =  39
      IETOA( 45) =   5
      IATOE(  5) =  45
      IETOA( 46) =   6
      IATOE(  6) =  46
      IETOA( 47) =   7
      IATOE(  7) =  47
      IETOA( 50) =  22
      IATOE( 22) =  50
      IETOA( 53) =  30
      IATOE( 30) =  53
      IETOA( 55) =   4
      IATOE(  4) =  55
      IETOA( 60) =  20
      IATOE( 20) =  60
      IETOA( 61) =  21
      IATOE( 21) =  61
      IETOA( 63) =  26
      IATOE( 26) =  63
      IETOA( 64) =  32
      IATOE( 32) =  64
      IETOA( 74) =  91
      IATOE( 91) =  74
      IETOA( 75) =  46
      IATOE( 46) =  75
      IETOA( 76) =  60
      IATOE( 60) =  76
      IETOA( 77) =  40
      IATOE( 40) =  77
      IETOA( 78) =  43
      IATOE( 43) =  78
      IETOA( 79) =  33
      IATOE( 33) =  79
      IETOA( 80) =  38
      IATOE( 38) =  80
      IETOA( 90) =  93
      IATOE( 93) =  90
      IETOA( 91) =  36
      IATOE( 36) =  91
      IETOA( 92) =  42
      IATOE( 42) =  92
      IETOA( 93) =  41
      IATOE( 41) =  93
      IETOA( 94) =  59
      IATOE( 59) =  94
      IETOA( 95) =  94
      IATOE( 94) =  95
      IETOA( 96) =  45
      IATOE( 45) =  96
      IETOA( 97) =  47
      IATOE( 47) =  97
      IETOA(106) = 124
      IATOE(124) = 106
      IETOA(107) =  44
      IATOE( 44) = 107
      IETOA(108) =  37
      IATOE( 37) = 108
      IETOA(109) =  95
      IATOE( 95) = 109
      IETOA(110) =  62
      IATOE( 62) = 110
      IETOA(111) =  63
      IATOE( 63) = 111
      IETOA(121) =  96
      IATOE( 96) = 121
      IETOA(122) =  58
      IATOE( 58) = 122
      IETOA(123) =  35
      IATOE( 35) = 123
      IETOA(124) =  64
      IATOE( 64) = 124
      IETOA(125) =  39
      IATOE( 39) = 125
      IETOA(126) =  61
      IATOE( 61) = 126
      IETOA(127) =  34
      IATOE( 34) = 127
      IETOA(129) =  97
      IATOE( 97) = 129
      IETOA(130) =  98
      IATOE( 98) = 130
      IETOA(131) =  99
      IATOE( 99) = 131
      IETOA(132) = 100
      IATOE(100) = 132
      IETOA(133) = 101
      IATOE(101) = 133
      IETOA(134) = 102
      IATOE(102) = 134
      IETOA(135) = 103
      IATOE(103) = 135
      IETOA(136) = 104
      IATOE(104) = 136
      IETOA(137) = 105
      IATOE(105) = 137
      IETOA(145) = 106
      IATOE(106) = 145
      IETOA(146) = 107
      IATOE(107) = 146
      IETOA(147) = 108
      IATOE(108) = 147
      IETOA(148) = 109
      IATOE(109) = 148
      IETOA(149) = 110
      IATOE(110) = 149
      IETOA(150) = 111
      IATOE(111) = 150
      IETOA(151) = 112
      IATOE(112) = 151
      IETOA(152) = 113
      IATOE(113) = 152
      IETOA(153) = 114
      IATOE(114) = 153
      IETOA(161) = 126
      IATOE(126) = 161
      IETOA(162) = 115
      IATOE(115) = 162
      IETOA(163) = 116
      IATOE(116) = 163
      IETOA(164) = 117
      IATOE(117) = 164
      IETOA(165) = 118
      IATOE(118) = 165
      IETOA(166) = 119
      IATOE(119) = 166
      IETOA(167) = 120
      IATOE(120) = 167
      IETOA(168) = 121
      IATOE(121) = 168
      IETOA(169) = 122
      IATOE(122) = 169
      IETOA(173) =  91
      IATOE( 91) = 173
      IETOA(176) =  48
      IATOE( 48) = 176
      IETOA(177) =  49
      IATOE( 49) = 177
      IETOA(178) =  50
      IATOE( 50) = 178
      IETOA(179) =  51
      IATOE( 51) = 179
      IETOA(180) =  52
      IATOE( 52) = 180
      IETOA(181) =  53
      IATOE( 53) = 181
      IETOA(182) =  54
      IATOE( 54) = 182
      IETOA(183) =  55
      IATOE( 55) = 183
      IETOA(184) =  56
      IATOE( 56) = 184
      IETOA(185) =  57
      IATOE( 57) = 185
      IETOA(189) =  93
      IATOE( 93) = 189
      IETOA(192) = 123
      IATOE(123) = 192
      IETOA(193) =  65
      IATOE( 65) = 193
      IETOA(194) =  66
      IATOE( 66) = 194
      IETOA(195) =  67
      IATOE( 67) = 195
      IETOA(196) =  68
      IATOE( 68) = 196
      IETOA(197) =  69
      IATOE( 69) = 197
      IETOA(198) =  70
      IATOE( 70) = 198
      IETOA(199) =  71
      IATOE( 71) = 199
      IETOA(200) =  72
      IATOE( 72) = 200
      IETOA(201) =  73
      IATOE( 73) = 201
      IETOA(208) = 125
      IATOE(125) = 208
      IETOA(209) =  74
      IATOE( 74) = 209
      IETOA(210) =  75
      IATOE( 75) = 210
      IETOA(211) =  76
      IATOE( 76) = 211
      IETOA(212) =  77
      IATOE( 77) = 212
      IETOA(213) =  78
      IATOE( 78) = 213
      IETOA(214) =  79
      IATOE( 79) = 214
      IETOA(215) =  80
      IATOE( 80) = 215
      IETOA(216) =  81
      IATOE( 81) = 216
      IETOA(217) =  82
      IATOE( 82) = 217
      IETOA(224) =  92
      IATOE( 92) = 224
      IETOA(226) =  83
      IATOE( 83) = 226
      IETOA(227) =  84
      IATOE( 84) = 227
      IETOA(228) =  85
      IATOE( 85) = 228
      IETOA(229) =  86
      IATOE( 86) = 229
      IETOA(230) =  87
      IATOE( 87) = 230
      IETOA(231) =  88
      IATOE( 88) = 231
      IETOA(232) =  89
      IATOE( 89) = 232
      IETOA(233) =  90
      IATOE( 90) = 233
      IETOA(240) =  48
      IATOE( 48) = 240
      IETOA(241) =  49
      IATOE( 49) = 241
      IETOA(242) =  50
      IATOE( 50) = 242
      IETOA(243) =  51
      IATOE( 51) = 243
      IETOA(244) =  52
      IATOE( 52) = 244
      IETOA(245) =  53
      IATOE( 53) = 245
      IETOA(246) =  54
      IATOE( 54) = 246
      IETOA(247) =  55
      IATOE( 55) = 247
      IETOA(248) =  56
      IATOE( 56) = 248
      IETOA(249) =  57
      IATOE( 57) = 249
 
C  SHOW SOME RESULTS
C  -----------------
 
      IF(PRINT) THEN
         PRINT100,NBYTW,NBITW,NREV,(IORD(I),I=1,NBYTW)
         IF(IASCII.EQ.0) PRINT*,'EBCDIC IS NATIVE LANGUAGE'
         IF(IASCII.EQ.1) PRINT*,'ASCII  IS NATIVE LANGUAGE'
      ENDIF
100   FORMAT(' WRDLEN:NBYTW=',I1,' NBITW=',I2,' IREV=',I1,' IORD=',8I1)
 
      RETURN
900   CALL BORT('WRDLEN - A WORD IS MORE THAN 64 BITS')
901   CALL BORT('WRDLEN - A WORD IS NOT MADE OF BYTES')
902   CALL BORT('WRDLEN - BYTE ORDER CHECKING MISTAKE')
      END
C-----------------------------------------------------------------------
C  WRITDX IS CALLED TO INITIALIZE DESCRIPTOR PROCESSING TABLES FOR A    
C  BUFR OUTPUT FILE CONNECTED TO UNIT LUNIT FROM A DX TABLE SOURCE      
C  CONNECTED TO UNIT LUNDX. IN THIS CASE, SINCE LUNIT IS PRESUMABLY     
C  CONNECTED TO AN OUTPUT FILE, LUNDX CANNOT BE THE SAME AS LUNIT.      
C  SUBROUTINE READDX IS USED TO READ THE DX TABLES FROM UNIT LUNDX,     
C  AND BUFR DX (DICTIONARY) MESSAGES CONTAINING THIS INFORMATION ARE    
C  WRITTEN AS THE INITIAL CONTENTS OF THE FILE CONNECTED TO UNIT LUNIT. 
C                                                                       
C  INPUT ARGUMENTS:                                                     
C     LUNIT    - UNIT CONNECTED TO BUFR FILE TO BE INITIALIZED/UPDATED  
C     LUN      - INTERNAL BUFR UNIT ASSOCIATED WITH FORTRAN UNIT LUNIT  
C     LUNDX    - UNIT CONTAINING DX-TABLES                              
C                                                                       
C-----------------------------------------------------------------------
      SUBROUTINE WRITDX(LUNIT,LUN,LUNDX)                                
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*6   ADN30                                               
      CHARACTER*1   MOCT(24000)                                         
      DIMENSION     MBAY(5000)                                          
      EQUIVALENCE   (MOCT(1),MBAY(1))                                   
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK UNITS, GET A DX TABLE, AND START A DX MESSAGE                  
C  ---------------------------------------------------                  
                                                                        
      IF(LUNIT.EQ.LUNDX) GOTO 900                                       
      CALL READDX(LUNIT,LUN,LUNDX)                                      
      CALL DXMINI(LUN,MBAY,MBYT,MBY4,MBYA,MBYB,MBYD)                    
                                                                        
      LDA = LDXA(IDXV+1)                                                
      LDB = LDXB(IDXV+1)                                                
      LDD = LDXD(IDXV+1)                                                
      L30 = LD30(IDXV+1)                                                
                                                                        
C  COPY TABLE A CONTENTS TO A BUFR DX MESSAGE                           
C  ------------------------------------------                           
                                                                        
      DO I=1,NTBA(LUN)                                                  
      IF(MBYT+LDA.GT.MAXDX) THEN                                        
         CALL MSGWRT(LUNIT,MBAY,MBYT)                                   
         CALL DXMINI(LUN,MBAY,MBYT,MBY4,MBYA,MBYB,MBYD)                 
      ENDIF                                                             
      CALL IPKM(MOCT(MBY4),3,IUPM(MOCT(MBY4),24)+LDA)                   
      CALL IPKM(MOCT(MBYA),1,IUPM(MOCT(MBYA), 8)+  1)                   
      MBIT = 8*(MBYB-1)                                                 
      CALL PKC(TABA(I,LUN),LDA,MBAY,MBIT)                               
      CALL PKB(          0,  8,MBAY,MBIT)                               
      CALL PKB(          0,  8,MBAY,MBIT)                               
      MBYT = MBYT+LDA                                                   
      MBYB = MBYB+LDA                                                   
      MBYD = MBYD+LDA                                                   
      ENDDO                                                             
                                                                        
C  COPY TABLE B CONTENTS TO A BUFR DX MESSAGE                           
C  ------------------------------------------                           
                                                                        
      DO I=1,NTBB(LUN)                                                  
      IF(MBYT+LDB.GT.MAXDX) THEN                                        
         CALL MSGWRT(LUNIT,MBAY,MBYT)                                   
         CALL DXMINI(LUN,MBAY,MBYT,MBY4,MBYA,MBYB,MBYD)                 
      ENDIF                                                             
      CALL IPKM(MOCT(MBY4),3,IUPM(MOCT(MBY4),24)+LDB)                   
      CALL IPKM(MOCT(MBYB),1,IUPM(MOCT(MBYB), 8)+  1)                   
      MBIT = 8*(MBYD-1)                                                 
      CALL PKC(TABB(I,LUN),LDB,MBAY,MBIT)                               
      CALL PKB(          0,  8,MBAY,MBIT)                               
      MBYT = MBYT+LDB                                                   
      MBYD = MBYD+LDB                                                   
      ENDDO                                                             
                                                                        
C  COPY TABLE D CONTENTS TO A BUFR DX MESSAGE                           
C  ------------------------------------------                           
                                                                        
      DO I=1,NTBD(LUN)                                                  
      NSEQ = IUPM(TABD(I,LUN)(LDD+1:LDD+1),8)                           
      LEND = LDD+1 + L30*NSEQ                                           
      IF(MBYT+LEND.GT.MAXDX) THEN                                       
         CALL MSGWRT(LUNIT,MBAY,MBYT)                                   
         CALL DXMINI(LUN,MBAY,MBYT,MBY4,MBYA,MBYB,MBYD)                 
      ENDIF                                                             
      CALL IPKM(MOCT(MBY4),3,IUPM(MOCT(MBY4),24)+LEND)                  
      CALL IPKM(MOCT(MBYD),1,IUPM(MOCT(MBYD), 8)+   1)                  
      MBIT = 8*(MBYT-4)                                                 
      CALL PKC(TABD(I,LUN),LDD,MBAY,MBIT)                               
      CALL PKB(       NSEQ,  8,MBAY,MBIT)                               
         DO J=1,NSEQ                                                    
         JJ  = LDD+2 + (J-1)*2                                          
         IDN = IUPM(TABD(I,LUN)(JJ:JJ),16)                              
         CALL PKC(ADN30(IDN,L30),L30,MBAY,MBIT)                         
         ENDDO                                                          
      MBYT = MBYT+LEND                                                  
      ENDDO                                                             
                                                                        
C  WRITE THE UNWRITTEN MESSAGE                                          
C  ---------------------------                                          
                                                                        
      CALL MSGWRT(LUNIT,MBAY,MBYT)                                      
                                                                        
      RETURN                                                            
900   CALL BORT('WRITDX - INPUR AND OUTPUT UNIT MUST NOT BE THE SAME') 
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE WRITSA(LUNXX,MSGT,MSGL)                                
                                                                        
      COMMON /BUFRMG/ MSGLEN,MSGTXT(5000)                               
                                                                        
      DIMENSION MSGT(*)                                                 
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      LUNIT = ABS(LUNXX)                                                
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.EQ.0) GOTO 902                                              
                                                                        
C  SEE IF A MEMORY MESSAGE IS WAITING OR FORCED                         
C  --------------------------------------------                         
                                                                        
      IF(LUNXX.LT.0) CALL CLOSMG(LUNIT)                                 
                                                                        
      IF(MSGLEN.GT.0) THEN                                              
         MSGL = MSGLEN                                                  
         DO N=1,MSGL                                                    
         MSGT(N) = MSGTXT(N)                                            
         ENDDO                                                          
         MSGLEN = 0                                                     
      ELSE                                                              
         MSGL = 0                                                       
      ENDIF                                                             
                                                                        
      IF(LUNXX.LT.0) RETURN                                             
                                                                        
C  PACK UP THE SUBSET AND PUT IT INTO THE MESSAGE                       
C  ----------------------------------------------                       
                                                                        
      CALL WRTREE(LUN)                                                  
      CALL MSGUPD(LUNIT,LUN)                                            
                                                                        
      RETURN                                                            
900   CALL BORT('WRITSA - FILE IS CLOSED                     ')        
901   CALL BORT('WRITSA - FILE IS OPEN FOR INPUT             ')        
902   CALL BORT('WRITSA - NO MESSAGE OPEN                    ')        
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE WRITSB(LUNIT)                                          
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.EQ.0) GOTO 902                                              
                                                                        
C  PACK UP THE SUBSET AND PUT IT INTO THE MESSAGE                       
C  ----------------------------------------------                       
                                                                        
      CALL WRTREE(LUN)                                                  
      CALL MSGUPD(LUNIT,LUN)                                            
                                                                        
      RETURN                                                            
900   CALL BORT('WRITSB - FILE IS CLOSED                     ')        
901   CALL BORT('WRITSB - FILE IS OPEN FOR INPUT             ')        
902   CALL BORT('WRITSB - NO MESSAGE OPEN                    ')        
      END                                                               
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE WRTREE(LUN)                                            
                                                                        
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*8  CVAL                                                 
      CHARACTER*3  TYP                                                  
      DIMENSION    IVAL(15000)                                          
      EQUIVALENCE  (CVAL,RVAL)                                          
      REAL*8       VAL,RVAL,BMISS                                       
      DATA         BMISS/10E10/
                                                                        
C-----------------------------------------------------------------------
      PKS(NODE) = VAL(N,LUN)*10.**ISC(NODE)-IRF(NODE)                   
C-----------------------------------------------------------------------
                                                                        
C  CONVERT USER NUMBERS INTO SCALED INTEGERS                            
C  -----------------------------------------                            
                                                                        
      DO N=1,NVAL(LUN)                                                  
      NODE = INV(N,LUN)                                                 
      IF(ITP(NODE).EQ.1) THEN                                           
         IVAL(N) = VAL(N,LUN)                                           
      ELSEIF(TYP(NODE).EQ.'NUM') THEN                                       
         IF(VAL(N,LUN).NE.BMISS) THEN
            IVAL(N) = ANINT(PKS(NODE))
         ELSE
            IVAL(N) = -1
         ENDIF
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  PACK THE USER ARRAY INTO THE SUBSET BUFFER                           
C  ------------------------------------------                           
                                                                        
      IBIT = 16                                                         
                                                                        
      DO N=1,NVAL(LUN)                                                  
      NODE = INV(N,LUN)                                                 
      IF(ITP(NODE).LT.3) THEN                                           
         CALL PKB(IVAL(N),IBT(NODE),IBAY,IBIT)                          
      ELSE                                                              
         RVAL = VAL(N,LUN)                                              
         CALL PKC(CVAL,IBT(NODE)/8,IBAY,IBIT)                           
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE WTSTAT(LUNIT,LUN,IL,IM)                                
                                                                        
      COMMON /STBFR/ IOLUN(32),IOMSG(32)                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK ON THE ARGUMENTS                                               
C  ----------------------                                               
                                                                        
      IF(LUNIT.LE.0)            GOTO 900                                
      IF(LUN  .LE.0)            GOTO 901                                
      IF(IL.LT.-1 .OR. IL.GT.1) GOTO 902                                
      IF(IM.LT. 0 .OR. IL.GT.1) GOTO 903                                
                                                                        
C  CHECK ON LUNIT-LUN COMBINATION                                       
C  ------------------------------                                       
                                                                        
      IF(ABS(IOLUN(LUN)).NE.LUNIT) THEN                                 
         IF(IOLUN(LUN).NE.0) GOTO 905                                   
      ENDIF                                                             
                                                                        
C  RESET THE FILE STATUSES                                              
C  -----------------------                                              
                                                                        
      IF(IL.NE.0) THEN                                                  
         IOLUN(LUN) = SIGN(LUNIT,IL)                                    
         IOMSG(LUN) = IM                                                
      ELSE                                                              
         IOLUN(LUN) = 0                                                 
         IOMSG(LUN) = 0                                                 
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('WTSTAT - BAD LUNIT                               ')   
901   CALL BORT('WTSTAT - BAD LUN                                 ')   
902   CALL BORT('WTSTAT - BAD IL                                  ')   
903   CALL BORT('WTSTAT - BAD IM                                  ')   
905   CALL BORT('WTSTAT - ATTEMPT TO REDEFINE EXISITING FILE UNIT ')   
      END                                                               


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

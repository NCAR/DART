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

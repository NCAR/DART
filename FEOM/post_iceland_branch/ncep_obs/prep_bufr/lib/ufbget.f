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

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

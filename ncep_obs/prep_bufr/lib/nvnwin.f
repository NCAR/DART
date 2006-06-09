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

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

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

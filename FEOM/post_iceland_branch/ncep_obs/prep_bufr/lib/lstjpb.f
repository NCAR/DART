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

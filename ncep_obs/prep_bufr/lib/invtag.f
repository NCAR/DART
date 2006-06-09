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

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

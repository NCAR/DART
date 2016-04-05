C---------------------------------------------------------------------- 
C  WILL MAKE ONE COPY OF EACH UNIQUE ELEMENT IN AN INPUT SUBSET BUFFER  
C  INTO IDENTICAL MNEMONIC SLOT IN THE OUTPUT SUBSET BUFFER             
C---------------------------------------------------------------------- 
      SUBROUTINE UFBCUP(LUBIN,LUBOT)                                    
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG,TAGI(15000),TAGO                                 
      CHARACTER*3  TYP                                                  
      DIMENSION    NINI(15000)                                          
      REAL*8       VAL                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE FILE STATUSES AND I-NODE                                   
C  ----------------------------------                                   
                                                                        
      CALL STATUS(LUBIN,LUI,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUI).NE.INV(1,LUI)) GOTO 902                             
                                                                        
      CALL STATUS(LUBOT,LUO,IL,IM)                                      
      IF(IL.EQ.0) GOTO 903                                              
      IF(IM.EQ.0) GOTO 904                                              
                                                                        
C  MAKE A LIST OF UNIQUE TAGS IN INPUT BUFFER                           
C  ------------------------------------------                           
                                                                        
      NTAG = 0                                                          
                                                                        
      DO 5 NI=1,NVAL(LUI)                                               
      NIN = INV(NI,LUI)                                                 
      IF(ITP(NIN).GE.2) THEN                                            
         DO NV=1,NTAG                                                   
         IF(TAGI(NV).EQ.TAG(NIN)) GOTO 5                                
         ENDDO                                                          
         NTAG = NTAG+1                                                  
         NINI(NTAG) = NI                                                
         TAGI(NTAG) = TAG(NIN)                                          
      ENDIF                                                             
5     ENDDO                                                             
                                                                        
      IF(NTAG.EQ.0) GOTO 905                                            
                                                                        
C  GIVEN A list MAKE ONE COPY OF COMMON ELEMENTS TO OUTPUT BUFFER       
C  --------------------------------------------------------------       
                                                                        
      DO 10 NV=1,NTAG                                                   
      NI = NINI(NV)                                                     
      DO NO=1,NVAL(LUO)                                                 
      TAGO = TAG(INV(NO,LUO))                                           
      IF(TAGI(NV).EQ.TAGO) THEN                                         
         VAL(NO,LUO) = VAL(NI,LUI)                                      
         GOTO 10                                                        
      ENDIF                                                             
      ENDDO                                                             
10    ENDDO                                                             
                                                                        
C  ALL EXITS HERE                                                       
C  --------------                                                       
                                                                        
      RETURN                                                            
900   CALL BORT('UFBCUP - INPUT  FILE IS NOT OPEN             ')       
901   CALL BORT('UFBCUP - INPUT  MESG IS NOT OPEN             ')       
902   CALL BORT('UFBCUP - INPUT  I-NODE  MISMATCH             ')       
903   CALL BORT('UFBCUP - OUTPUT FILE IS NOT OPEN             ')       
904   CALL BORT('UFBCUP - OUTPUT MESG IS NOT OPEN             ')       
905   CALL BORT('UFBCUP - NO TAGS IN INPUT BUFFER             ')       
      END                                                               

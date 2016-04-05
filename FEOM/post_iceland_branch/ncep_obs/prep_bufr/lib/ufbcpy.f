C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBCPY(LUBIN,LUBOT)                                    
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      REAL*8 VAL                                                        
                                                                        
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
      IF(INODE(LUI).NE.INODE(LUO)) GOTO 905                             
                                                                        
C  EVERYTHING OKAY COPY USER ARRAY FROM LUI TO LUO                      
C  -----------------------------------------------                      
                                                                        
      NVAL(LUO) = NVAL(LUI)                                             
                                                                        
      DO N=1,NVAL(LUI)                                                  
      INV(N,LUO) = INV(N,LUI)                                           
      VAL(N,LUO) = VAL(N,LUI)                                           
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('UFBCPY - INPUT  FILE IS NOT OPEN             ')       
901   CALL BORT('UFBCPY - INPUT  MESG IS NOT OPEN             ')       
902   CALL BORT('UFBCPY - INPUT  I-NODE  MISMATCH             ')       
903   CALL BORT('UFBCPY - OUTPUT FILE IS NOT OPEN             ')       
904   CALL BORT('UFBCPY - OUTPUT MESG IS NOT OPEN             ')       
905   CALL BORT('UFBCPY - IN/OUT I-NODE  MISMATCH             ')       
      END                                                               

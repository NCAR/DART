C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION NMSUB(LUNIT)                                             
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NMSUB = 0                                                         
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      IF(IM.EQ.0) GOTO 902                                              
                                                                        
      NMSUB = MSUB(LUN)                                                 
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('NMSUB - FILE IS CLOSED                   ')           
901   CALL BORT('NMSUB - FILE IS OPEN FOR OUTPUT          ')           
902   CALL BORT('NMSUB - NO MESSAGE IS OPENED             ')           
      END                                                               

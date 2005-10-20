C---------------------------------------------------------------------- 
C  WILL CHECK TO SEE IF ANY UNREAD SUBSETS IN AN OPEN INPUT MESSAGE     
C---------------------------------------------------------------------- 
      FUNCTION IFBGET(LUNIT)                                            
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  MAKE SURE A FILE/MESSAGE IS OPEN FOR INPUT                           
C  ------------------------------------------                           
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.GE.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
                                                                        
C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE                        
C  ---------------------------------------------                        
                                                                        
      IF(NSUB(LUN).LT.MSUB(LUN)) THEN                                   
         IFBGET = 0                                                     
      ELSE                                                              
         IFBGET = -1                                                    
      ENDIF                                                             
                                                                        
C  EXIT ONE WAY OR ANOTHER                                              
C  -----------------------                                              
                                                                        
      RETURN                                                            
900   CALL BORT('IFBGET - FILE NOT OPEN FOR INPUT')                    
901   CALL BORT('IFBGET - NO MESSAGE OPEN        ')                    
      END                                                               

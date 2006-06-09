C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE NXTWIN(LUN,IWIN,JWIN)                                  
                                                                        
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
cfpp$ expand (lstrpc)                                                   
C---------------------------------------------------------------------- 
                                                                        
      IF(JWIN.EQ.NVAL(LUN)) THEN                                        
         IWIN = 0                                                       
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  FIND THE NEXT SEQUENTIAL WINDOW                                      
C  -------------------------------                                      
                                                                        
      NODE = INV(IWIN,LUN)                                              
      IF(LSTRPC(NODE,LUN).NE.NODE) print*,'bad node=',node,iwin         
      IF(LSTRPC(NODE,LUN).NE.NODE) CALL BORT('NXTWIN - NOT RPC')       
      IF(VAL(JWIN,LUN).EQ.0) THEN                                       
         IWIN = 0                                                       
      ELSE                                                              
         IWIN = JWIN                                                    
         JWIN = IWIN+VAL(IWIN,LUN)                                      
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               

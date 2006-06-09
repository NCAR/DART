C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE NEWWIN(LUN,IWIN,JWIN)                                  
                                                                        
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
cfpp$ expand (lstrpc)                                                   
C---------------------------------------------------------------------- 
                                                                        
      IF(IWIN.EQ.1) THEN                                                
         JWIN = NVAL(LUN)                                               
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  REFIND THE JWIN BOUNDARY FROM IWIN                                   
C  ----------------------------------                                   
                                                                        
      NODE = INV(IWIN,LUN)                                              
      IF(LSTRPC(NODE,LUN).NE.NODE) CALL BORT('NEWWIN - NOT RPC')       
      JWIN = IWIN+VAL(IWIN,LUN)                                         
                                                                        
      RETURN                                                            
      END                                                               

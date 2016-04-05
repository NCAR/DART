C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE GETWIN(NODE,LUN,IWIN,JWIN)                             
                                                                        
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
cfpp$ expand (lstrpc)                                                   
C---------------------------------------------------------------------- 
                                                                        
      IRPC = LSTRPC(NODE,LUN)                                           
                                                                        
      IF(IRPC.EQ.0) THEN                                                
         IWIN = INVWIN(NODE,LUN,JWIN,NVAL(LUN))                         
         IF(IWIN.EQ.0 .and. jwin.gt.1) RETURN                           
         IWIN = 1                                                       
         JWIN = NVAL(LUN)                                               
         RETURN                                                         
      ELSE                                                              
         IWIN = INVWIN(IRPC,LUN,JWIN,NVAL(LUN))                         
         IF(IWIN.EQ.0) THEN                                             
            RETURN                                                      
         ELSEIF(VAL(IWIN,LUN).EQ.0.) THEN                               
            IWIN = 0                                                    
            RETURN                                                      
         ENDIF                                                          
      ENDIF                                                             
                                                                        
      JWIN = INVWIN(IRPC,LUN,IWIN+1,NVAL(LUN))                          
      IF(JWIN.EQ.0) CALL BORT('GETWIN - MISSING BRACKET')              
                                                                        
      RETURN                                                            
      END                                                               

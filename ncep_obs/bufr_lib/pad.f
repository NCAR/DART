C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE PAD(IBAY,IBIT,IBYT,IPADB)                              
                                                                        
      DIMENSION IBAY(*)                                                 
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  PAD THE SUBSET TO AN IPADB BIT BOUNDARY                              
C  ----------------------------------------                             
                                                                        
      IPAD = IPADB - MOD(IBIT+8,IPADB)                                  
      CALL PKB(IPAD,8,IBAY,IBIT)                                        
      CALL PKB(   0,IPAD,IBAY,IBIT)                                     
      IBYT = IBIT/8                                                     
                                                                        
      IF(MOD(IBIT,IPADB).NE.0) GOTO 900                                 
      IF(MOD(IBIT,8    ).NE.0) GOTO 900                                 
                                                                        
      RETURN                                                            
900   CALL BORT('PAD - BIT PAD FAILURE              ')                 
      END                                                               

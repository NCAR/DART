C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION IREADERM(LUNIT,SUBSET,IDATE)                             
      CHARACTER*8 SUBSET                                                
      CALL READERM(LUNIT,SUBSET,IDATE,IRET)                         
      IREADERM = IRET
      RETURN                                                            
                                                                        
      ENTRY IREADERS(LUNIT)                                             
      CALL READERS(LUNIT,IRET)                                      
      IREADERS = IRET
      RETURN                                                            
                                                                        
      END                                                               

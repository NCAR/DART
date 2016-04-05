C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION IREADMG(LUNIT,SUBSET,IDATE)                              
      CHARACTER*8 SUBSET                                                
      CALL READMG(LUNIT,SUBSET,IDATE,IRET)                           
      IREADMG = IRET
      RETURN                                                            
                                                                        
      ENTRY IREADMM(IMSG,SUBSET,IDATE)                                 
      CALL READMM(IMSG,SUBSET,IDATE,IRET)                         
      IREADMM = IRET
      RETURN                                                            

      ENTRY IREADFT(LUNIT,SUBSET,IDATE)                                 
      CALL READFT(LUNIT,SUBSET,IDATE,IRET)                           
      IREADFT = IRET
      RETURN                                                            
                                                                        
      ENTRY IREADNS(LUNIT,SUBSET,IDATE)                                 
      CALL READNS(LUNIT,SUBSET,IDATE,IRET)                           
      IREADNS = IRET
      RETURN                                                            
                                                                        
      ENTRY IREADSB(LUNIT)                                              
      CALL READSB(LUNIT,IRET)                                        
      IREADSB = IRET
      RETURN                                                            
                                                                        
      ENTRY ICOPYSB(LUNIN,LUNOT)                                        
      CALL COPYSB(LUNIN,LUNOT,IRET)                                  
      ICOPYSB = IRET
      RETURN                                                            
                                                                        
      ENTRY IREADIBM(LUNIT,SUBSET,IDATE)                              
      CALL READIBM(LUNIT,SUBSET,IDATE,IRET)                           
      IREADIBM = IRET
      RETURN                                                            

      END                                                               

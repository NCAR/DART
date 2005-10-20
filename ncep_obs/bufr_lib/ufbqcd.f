C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBQCD(LUNIT,NEMO,QCD)                                 
                                                                        
      CHARACTER*(*) NEMO                                                
      CHARACTER*6  FXY,ADN30                                            
      CHARACTER*1  TAB                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
                                                                        
      CALL NEMTAB(LUN,NEMO,IDN,TAB,IRET)                                
      IF(TAB.NE.'D') GOTO 901                                           
                                                                        
      FXY = ADN30(IDN,6)                                                
      IF(FXY(2:3).NE.'63') GOTO 902                                     
      READ(FXY(4:6),'(F3.0)',ERR=903) QCD                               
                                                                        
      RETURN                                                            
900   CALL BORT('UFBQCD - FILE IS CLOSED                       ')      
901   CALL BORT('UFBQCD - MISSING OR INVALID TABLE D QC CODE   ')      
902   CALL BORT('UFBQCD - TABLE D QC CODE DESCRIPTOR NOT 363YYY')      
903   CALL BORT('UFBQCD - ERROR READING YYY FROM QC CODE DESCRP')      
      END                                                               

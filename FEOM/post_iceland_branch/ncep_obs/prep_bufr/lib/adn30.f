C---------------------------------------------------------------------- 
C  CONVERT AN INTEGER DESCRIPTOR TO FIVE OR SIX CHARACTER ASCII FORMAT  
C---------------------------------------------------------------------- 
      FUNCTION ADN30(IDN,L30)                                           
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*(*) ADN30                                               
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(LEN(ADN30).LT.L30         ) GOTO 900                           
      IF(IDN.LT.0 .OR. IDN.GT.65535) GOTO 901                           
      IF(L30.EQ.5) THEN                                                 
         WRITE(ADN30,'(I5)') IDN                                        
      ELSEIF(L30.EQ.6) THEN                                             
         IDF = ISHFT(IDN,-14)                                           
         IDX = ISHFT(ISHFT(IDN,NBITW-14),-(NBITW-6))                    
         IDY = ISHFT(ISHFT(IDN,NBITW- 8),-(NBITW-8))                    
         WRITE(ADN30,'(I1,I2,I3)') IDF,IDX,IDY                          
      ELSE                                                              
         GOTO 902                                                       
      ENDIF                                                             
                                                                        
      DO I=1,L30                                                        
      IF(ADN30(I:I).EQ.' ') ADN30(I:I) = '0'                            
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('ADN30 - FUNCTION RETURN STRING TOO SHORT')            
901   CALL BORT('ADN30 - IDN OUT OF RANGE                ')            
902   CALL BORT('ADN30 - CHARACTER LENGTH L30 <> 5 OR 6')              
      END                                                               

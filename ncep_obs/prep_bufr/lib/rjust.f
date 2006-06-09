C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION RJUST(STR)                                               
      CHARACTER*(*) STR                                                 
      RJUST = 0                                                         
      IF(STR.EQ.' ') RETURN                                             
      LSTR = LEN(STR)                                                   
      DO WHILE(STR(LSTR:LSTR).EQ.' ')                                   
         DO I=LSTR,2,-1
         STR(I:I) = STR(I-1:I-1)                                        
         ENDDO
         STR(1:1) = ' '
      ENDDO                                                             
      RETURN                                                            
      END                                                               

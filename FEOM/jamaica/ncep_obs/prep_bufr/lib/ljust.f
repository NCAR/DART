C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION LJUST(STR)                                               
      CHARACTER*(*) STR                                                 
      LJUST = 0                                                         
      IF(STR.EQ.' ') RETURN                                             
      LSTR = LEN(STR)                                                   
      DO I=1,LSTR
      DO WHILE(STR(I:I).EQ.' ' .AND. STR(I+1:LSTR).NE.' ')
         STR(I:LSTR) = STR(I+1:LSTR)
      ENDDO                                                             
      ENDDO                                                             
      RETURN                                                            
      END                                                               

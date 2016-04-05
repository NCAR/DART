C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CLOSBF(LUNIT)                                          

      COMMON /NULBFR/ NULL(32)

      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.GT.0 .AND. IM.NE.0) CALL CLOSMG(LUNIT)                      
      CALL WTSTAT(LUNIT,LUN,0,0)                                        
      IF(NULL(LUN).EQ.0) THEN
        CLOSE(LUNIT)                                                      
      ENDIF

      RETURN                                                            
      END                                                               

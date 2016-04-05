C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      LOGICAL FUNCTION DIGIT(STR)
      CHARACTER*(*) STR
      DIGIT = .FALSE.
      DO I=1,LEN(STR)
      IF(STR(I:I).NE.'0' .AND. STR(I:I).NE.'1' .AND. 
     .   STR(I:I).NE.'2' .AND. STR(I:I).NE.'3' .AND. 
     .   STR(I:I).NE.'4' .AND. STR(I:I).NE.'5' .AND. 
     .   STR(I:I).NE.'6' .AND. STR(I:I).NE.'7' .AND. 
     .   STR(I:I).NE.'8' .AND. STR(I:I).NE.'9') RETURN
      ENDDO
      DIGIT = .TRUE.
      RETURN
      END

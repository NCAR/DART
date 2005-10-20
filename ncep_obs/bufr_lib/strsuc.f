      SUBROUTINE STRSUC(STR1,STR2,LENS)                                 

C************************************************************************
C* STRSUC								*
C*									*
C* This subroutine removes leading and trailing blanks from a string.	*
C*									*
C* STRSUC  ( STR1, STR2, LENS )						*
C*									*
C* Input parameters:							*
C*	STR1		CHARACTER*(*)	String				*
C*									*
C* Output parameters:							*
C*	STR2		CHARACTER*(*)	Copy of STR1 with leading and	*
C*					trailing blanks removed		*
C*	LENS		INTEGER		Length of STR2:			*
C*					  -1 = STR1 contained embedded	*
C*					       blanks			*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      CHARACTER*(*) STR1,STR2                                           
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      LENS = 0                                                          
      LSTR = LEN(STR1)                                                  

C     Find the first non-blank in the input string.
                                                                        
      DO I=1,LSTR                                                       
      IF(STR1(I:I).NE.' ') GOTO 2                                       
      ENDDO                                                             
      RETURN                                                            

C     Now, starting with the first non-blank in the input string,
C     copy characters from the input string into the output string
C     until reaching the next blank in the input string.
                                                                        
2     DO J=I,LSTR                                                       
      IF(STR1(J:J).EQ.' ') GOTO 3                                       
      LENS = LENS+1                                                     
      STR2(LENS:LENS) = STR1(J:J)                                       
      ENDDO                                                             
      RETURN                                                            

C     Now, continuing on within the input string, make sure that
C     there are no more non-blank characters.  If there are, then
C     the blank at which we stopped copying from the input string
C     into the output string was an embedded blank.
                                                                        
3     DO I=J,LSTR                                                       
      IF(STR1(I:I).NE.' ') LENS = -1                                    
      ENDDO                                                             
      RETURN                                                            
                                                                        
      END                                                               

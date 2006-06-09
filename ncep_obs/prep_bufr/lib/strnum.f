      SUBROUTINE STRNUM(STR,NUM)                                        

C************************************************************************
C* STRNUM								*
C*									*
C* This subroutine decodes an integer from a character string.		*
C* The input string should contain only digits and (optional) trailing	*
C* blanks and should *not* contain any sign characters (e.g. '+', '-')	*
C* nor leading blanks nor embedded blanks.				*
C*									*
C* STRNUM  ( STR, NUM )							*
C*									*
C* Input parameters:							*
C*	STR		CHARACTER*(*)	String containing encoded	*
C*					integer value			*
C*									*
C* Output parameters:							*
C*	NUM		INTEGER		Decoded integer			*
C*					  -1 = decode was unsuccessful	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*20  STR2                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NUM = 0                                                           
      K = 0                                                             
                                                                        
C     Note that, in the following call to subroutine STRSUC, the output
C     string STR2 is not used anywhere else in this routine.  In fact,
C     the only reason that subroutine STRSUC is being called here is to
C     determine NUM, which, owing to the fact that the input string STR
C     cannot contain any leading blanks, is equal to the number of
C     digits to be decoded from the beginning of STR.

      CALL STRSUC(STR,STR2,NUM)                                         
                                                                        
      DO I=1,NUM                                                        
      READ(STR(I:I),'(I1)',ERR=99) J                                    
      IF(J.EQ.0 .AND. STR(I:I).NE.'0') GOTO 99                          
      K = K*10+J                                                        
      ENDDO                                                             
                                                                        
      NUM = K                                                           
      RETURN                                                            
                                                                        
C     Note that NUM = -1 unambiguously indicates a bad decode since
C     the input string cannot contain sign characters; thus, NUM is
C     always positive if the decode is successful.

99    NUM = -1                                                          
      RETURN                                                            
      END                                                               

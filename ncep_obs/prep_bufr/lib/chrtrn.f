      SUBROUTINE CHRTRN(STR,CHR,N)                                      

C************************************************************************
C* CHRTRN								*
C*									*
C* This subroutine copies a specified number of characters from a	*
C* character array into a character string.				*
C*									*
C* CHRTRN  ( STR, CHR, N )						*
C*									*
C* Input parameters:							*
C*	CHR		CHARACTER(N)	Character array			*
C*	N		INTEGER		Number of characters to copy	*
C*									*
C* Output parameters:							*
C*	STR		CHARACTER*(N)	Character string		*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation				*
C************************************************************************
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*1   CHR(N)                                              
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      DO I=1,N                                                          
      STR(I:I) = CHR(I)                                                 
      ENDDO                                                             
      RETURN                                                            
      END                                                               

C---------------------------------------------------------------------- 
C  REAL NUMBER FROM A STRING                                            
C---------------------------------------------------------------------- 

      FUNCTION VALX(STR)                                                

C************************************************************************
C* VALX									*
C*									*
C* This function decodes a real number from a string.  If the decode	*
C* fails, then the value BMISS (=10E10) is returned.  Note that, unlike	*
C* for subroutine STRNUM, the input string may contain a leading sign	*
C* character (e.g. '+', '-').						*
C*									*
C* VALX  ( STR )							*
C*									*
C* Input parameters:							*
C*	STR		CHARACTER*(*)	String containing encoded	*
C*					real value			*
C*									*
C* Output parameters:							*
C*	VALX		REAL		Decoded real value		*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      CHARACTER*(*) STR
      CHARACTER*99  BSTR
      CHARACTER*8   FMT
      REAL*8        BMISS                                           
                                                                        
      DATA BMISS /10E10/
      data noinline /0/                                                 
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      LENS = LEN(STR)
      IF(LENS.GT.99) CALL BORT('VALX - ARG TOO LONG')
      BSTR(1:LENS) = STR            
      RJ = RJUST(BSTR(1:LENS))
      WRITE(FMT,'(''(F'',I2,''.0)'')') LENS                             
      READ(BSTR,FMT,ERR=900) VAL
      VALX = VAL                                                        
      RETURN                                                            
900   VALX = BMISS                                                      
      RETURN                                                            
      END                                                               

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE JSTIFY                                                 
C************************************************************************
C* JSTIFY								*
C*									*
C* This subroutine consists solely of two separate entry points for	*
C* left-justifying strings containing integer and non-integer values.	*
C*									*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*1  SIGN                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      ENTRY JSTCHR(STR)                                                 
C************************************************************************
C* JSTCHR								*
C*									*
C* This entry point removes all leading blanks from a string.		*
C*									*
C* JSTCHR  ( STR )							*
C*									*
C* Input parameters:							*
C*	STR		CHARACTER*(*)	String				*
C*									*
C* Output parameters:							*
C*	STR		CHARACTER*(*)	Copy of input STR with leading	*
C*					blanks removed			*
C************************************************************************
                                                                        
      LSTR = LEN(STR)                                                   
                                                                        
      IF(STR.EQ.' ') GOTO 900                                           
1     IF(STR(1:1).EQ.' ') THEN                                          
         STR  = STR(2:LSTR)                                             
         GOTO 1                                                         
      ENDIF                                                             
      RETURN                                                            
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      ENTRY JSTNUM(STR,SIGN,IRET)                                       
C************************************************************************
C* JSTNUM								*
C*									*
C* This entry point removes all leading blanks from a string containing	*
C* an encoded integer value.  If the value has a leading sign character	*
C* ('+' or '-'), then this character is also removed and is returned	*
C* separately within SIGN.						*
C*									*
C* JSTNUM  ( STR )							*
C*									*
C* Input parameters:							*
C*	STR		CHARACTER*(*)	String containing encoded	*
C*					integer value			*
C*									*
C* Output parameters:							*
C*	STR		CHARACTER*(*)	Copy of input STR with leading	*
C*					blanks and sign character	*
C*					removed				*
C*	SIGN		CHARACTER	Sign of encoded integer value:	*
C*					  '+' = positive value		*
C*					  '-' = negative value		*
C*	IRET		INTEGER		Return code:			*
C*					  0 = normal return		*
C*					 -1 = encoded value within STR	*
C*					      was not an integer	*
C************************************************************************
                                                                        
      IRET = 0                                                          
      LSTR = LEN(STR)                                                   
                                                                        
      IF(STR.EQ.' ') GOTO 900                                           
2     IF(STR(1:1).EQ.' ') THEN                                          
         STR  = STR(2:LSTR)                                             
         GOTO 2                                                         
      ENDIF                                                             
      IF(STR(1:1).EQ.'+') THEN                                          
         STR  = STR(2:LSTR)                                             
         SIGN = '+'                                                     
      ELSEIF(STR(1:1).EQ.'-') THEN                                      
         STR  = STR(2:LSTR)                                             
         SIGN = '-'                                                     
      ELSE                                                              
         SIGN = '+'                                                     
      ENDIF                                                             
                                                                        
      CALL STRNUM(STR,NUM)                                              
      IF(NUM.LT.0) IRET = -1                                            
      RETURN                                                            
                                                                        
900   CALL BORT('JSTIFY - BLANK STRING NOT ALLOWED')                   
      END                                                               

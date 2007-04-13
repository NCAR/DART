C-----------------------------------------------------------------------
C  PARSE SEPARATE WORDS FROM A STRING SEQUENCE                          
C-----------------------------------------------------------------------

      SUBROUTINE PARSEQ(STR,TAGS,MTAG,NTAG)                             

C************************************************************************
C* PARSEQ								*
C*									*
C* This subroutine parses a string containing one or more mnemonics	*
C* into an array of mnemonics.  The mnemonics within the string must	*
C* be separated by one or more blank characters.			*
C*									*
C* PARSEQ  ( STR, TAGS, MTAG, NTAG )					*
C*									*
C* Input parameters:							*
C*	STR		CHARACTER*(*)	String				*
C*	MTAG		INTEGER		Maximum number of mnemonics to	*
C*					be parsed from string		*
C*									*
C* Output parameters:							*
C*	TAGS (NTAG)	CHARACTER*(*)	Array of mnemonics		*
C*	NTAG		INTEGER		Number of mnemonics returned	*
C*									*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      CHARACTER*(*) STR,TAGS(MTAG)                                      
      CHARACTER*80  ASTR                                                
      LOGICAL       WORD                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      ASTR = STR                                                        
      LSTR = LEN(STR)                                                   
      LTAG = LEN(TAGS(1))                                               
      IF(LSTR.GT.80) GOTO 900                                           
      NTAG = 0                                                          
      NCHR = 0                                                          
      WORD = .FALSE.                                                    
                                                                        
      DO 10 I=1,LSTR                                                    
                                                                        
      IF(.NOT.WORD .AND. STR(I:I).NE.' ') THEN                          
         NTAG = NTAG+1                                                  
         IF(NTAG.GT.MTAG) GOTO 901                                      
         TAGS(NTAG) = ' '                                               
      ENDIF                                                             
                                                                        
      IF(WORD .AND. STR(I:I).EQ.' ') NCHR = 0                           
      WORD = STR(I:I).NE.' '                                            
                                                                        
      IF(WORD) THEN                                                     
         NCHR = NCHR+1                                                  
         IF(NCHR.GT.LTAG) GOTO 902                                      
         TAGS(NTAG)(NCHR:NCHR) = STR(I:I)                               
      ENDIF                                                             
                                                                        
10    CONTINUE                                                          
                                                                        
      RETURN                                                            
900   CALL BORT('PARSEQ - STRING TOO LONG  '//ASTR)                    
901   CALL BORT('PARSEQ - TOO MANY TAGS    '//ASTR)                    
902   CALL BORT('PARSEQ - TAG IS TOO LONG  '//ASTR)                    
      END                                                               

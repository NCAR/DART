      SUBROUTINE CHRTRNA(STR,CHR,N)                                     

C************************************************************************
C* CHRTRNA								*
C*									*
C* This subroutine copies a specified number of characters from a	*
C* character array into a character string.  The difference between	*
C* this subroutine and subroutine CHRTRN is that, in this subroutine,	*
C* the input character array is assumed to be in ASCII; thus, for cases	*
C* where the native machine is EBCDIC, an ASCII -> EBCDIC translation	*
C* is done on the final string before it is output.			*
C*									*
C* CHRTRNA  ( STR, CHR, N )						*
C*									*
C* Input parameters:							*
C*	CHR		CHARACTER(N)	Character array	in ASCII	*
C*	N		INTEGER		Number of characters to copy	*
C*									*
C* Output parameters:							*
C*	STR		CHARACTER*(N)	Character string in ASCII or	*
C*					EBCDIC, depending on native	*
C*					machine				*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /CHARAC/ IASCII,IATOE(0:255),IETOA(0:255)                  
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*1   CHR(N)                                              
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 

C     Loop on N characters of CHR

      DO I=1,N                                                          
      STR(I:I) = CHR(I)                                                 

C     If this is an EBCDIC machine, then translate the character
C     from ASCII -> EBCDIC.

      IF(IASCII.EQ.0) CALL IPKM(STR(I:I),1,IATOE(IUPM(STR(I:I),8)))     
      ENDDO                                                             
      RETURN                                                            
      END                                                               

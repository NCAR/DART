      FUNCTION NEMOCK(NEMO)                                             

C************************************************************************
C* NEMOCK								*
C*									*
C* This function checks a mnemonic that was read via subroutine RDUSDX	*
C* (i.e. a mnemonic that was read from a user-supplied DX table) to	*
C* verify that it has a length of between 1 and 8 characters and that	*
C* it only contains characters from the allowable character set.	*
C*									*
C* NEMOCK  ( NEMO )							*
C*									*
C* Input parameters:							*
C*	NEMO		CHARACTER*8	Mnemonic to be checked		*
C*									*
C* Output parameters:							*
C*	NEMOCK		INTEGER		Indicator as to whether		*
C*					NEMO is valid:			*
C*					  0 = yes			*
C*					 -1 = no			*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      CHARACTER*(*) NEMO                                                
      CHARACTER*38  CHRSET                                              
                                                                        
      DATA CHRSET /'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.'/            
      DATA NCHR   /38/                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  GET THE LENGTH OF NEMO                                               
C  ----------------------                                               
                                                                        
      LNEMO = 0                                                         
                                                                        
      DO I=LEN(NEMO),1,-1                                               
      IF(NEMO(I:I).NE.' ') THEN                                         
         LNEMO = I                                                      
         GOTO 1                                                         
      ENDIF                                                             
      ENDDO                                                             
                                                                        
1     IF(LNEMO.LT.1 .OR. LNEMO.GT.8) THEN                               
         NEMOCK = -1                                                    
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  SCAN NEMO FOR ALLOWABLE CHARACTERS                                   
C  ----------------------------------                                   
                                                                        
      DO 10 I=1,LNEMO                                                   
      DO J=1,NCHR                                                       
      IF(NEMO(I:I).EQ.CHRSET(J:J)) GOTO 10                              
      ENDDO                                                             
      NEMOCK = -1                                                       
      RETURN                                                            
10    ENDDO                                                             
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      NEMOCK = 0                                                        
      RETURN                                                            
      END                                                               

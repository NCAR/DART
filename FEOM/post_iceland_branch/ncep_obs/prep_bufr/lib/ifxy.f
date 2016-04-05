      FUNCTION IFXY(ADSC)                                               

C************************************************************************
C* IFXY									*
C*									*
C* This function returns the integer corresponding to the bitwise	*
C* representation of an input FXY value.				*
C*									*
C* IFXY  ( ADSC )							*
C*									*
C* Input parameters:							*
C*	ADSC		CHARACTER*6	FXY value			*
C*									*
C* Output parameters:							*
C*	IFXY		INTEGER		Integer corresponding to the	*
C*					bitwise representation of ADSC	*
C*									*
C************************************************************************
C*									*
C*	EXAMPLE:							*
C*									*
C*	If ADSC = '063022', then IFXY = 16150 since:			*
C*									*
C*	0	63	     22						*
C*	  |	      |							*
C*	F |	X     |	      Y						*
C*	  |	      |							*
C*     0 0 1 1 1 1 1 1 0 0 0 1 0 1 1 0	=				*
C*									*
C*	( 2**13 + 2**12 + 2**11 + 2**10 +				*
C*		2**9 + 2**8 + 2**4 + 2**2 + 2**1 )  = 16150		*
C*									*
C************************************************************************
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      CHARACTER*6 ADSC                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      READ(ADSC,'(I1,I2,I3)') IF,IX,IY                                  
      IFXY = IF*2**14 + IX*2**8 + IY                                    
      RETURN                                                            
      END                                                               

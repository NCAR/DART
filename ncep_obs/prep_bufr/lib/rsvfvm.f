      SUBROUTINE RSVFVM(NEM1,NEM2)                                      

C************************************************************************
C* RSVFVM								*
C*									*
C* This subroutine steps through the "following value" mnemonic NEM1	*
C* and, for each "." character encountered (except for the initial one),*
C* overwrites it with the next corresponding character from NEM2.	*
C*									*
C************************************************************************
C* For example, if, on input:	NEM1 = ".DTH...."			*
C*				NEM2 = "MXTM    "			*
C*									*
C*	     then, on output:	NEM1 = ".DTHMXTM"			*
C************************************************************************
C*									*
C* RSVFVM  ( NEM1, NEM2 )						*
C*									*
C* Input parameters:							*
C*	NEM1		CHARACTER*8	"Following value" mnemonic	*
C*	NEM2		CHARACTER*8	Mnemonic immediately following	*
C*					NEM1 within user DX table	*
C*									*
C* Output parameters:							*
C*	NEM1		CHARACTER*8	Copy of input NEM1 with all "."	*
C*					characters (except initial one)	*
C*					overwritten with corresponding	*
C*					characters from NEM2		*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      CHARACTER*8 NEM1,NEM2                                             
                                                                        
      DO I=1,LEN(NEM1)                                                  
      IF(I.EQ.1) THEN                                                   

C	 Skip initial "." and initialize J.

         J = 1                                                          
      ELSE                                                              
         IF(NEM1(I:I).EQ.'.') THEN                                      
            NEM1(I:I) = NEM2(J:J)                                       
            J = J+1                                                     
         ENDIF                                                          
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               

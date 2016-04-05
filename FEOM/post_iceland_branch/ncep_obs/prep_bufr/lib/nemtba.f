      SUBROUTINE NEMTBA(LUN,NEMO,MTYP,MSBT,INOD)                        

C************************************************************************
C* NEMTBA								*
C*									*
C* This subroutine searches for mnemonic NEMO within the internal BUFR	*
C* table A arrays and, if found, returns information about that	mnemonic*
C* from within these arrays.						*
C*									*
C* NEMTBA  ( LUN, NEMO, MTYP, MSBT, INOD )				*
C*									*
C* Input parameters:							*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					BUFR table arrays		*
C*	NEMO		CHARACTER*(*)	Table A mnemonic to search for	*
C*									*
C* Output parameters:							*
C*	MTYP		INTEGER		Message type corresponding	*
C*					to NEMO				*
C*	MSBT		INTEGER		Message subtype corresponding	*
C*					to NEMO				*
C*	INOD		INTEGER		Positional index of NEMO within	*
C*					internal jump/link table	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*(*) NEMO                                                
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*20  NEMT                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NEMT = NEMO                                                       
      IRET = 0                                                          
                                                                        
C  LOOK FOR NEMO IN TABLE A                                             
C  ------------------------                                             
                                                                        
      DO I=1,NTBA(LUN)                                                  
      IF(TABA(I,LUN)(4:11).EQ.NEMO) THEN                                
         MTYP = IDNA(I,LUN,1)                                           
         MSBT = IDNA(I,LUN,2)                                           
         INOD = MTAB(I,LUN)                                             
         IF(MTYP.LT.0 .OR. MTYP.GT.255) GOTO 900                        
         IF(MSBT.LT.0 .OR. MSBT.GT.255) GOTO 901                        
         RETURN                                                         
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      CALL BORT('NEMTBA - CANT FIND '//NEMT)                           
900   CALL BORT('NEMTBA - BAD MTYP  '//NEMT)                           
901   CALL BORT('NEMTBA - BAD MSBT  '//NEMT)                           
      END                                                               

      SUBROUTINE UPTDD(ID,LUN,IENT,IRET)                                

C************************************************************************
C* UPTDD								*
C*									*
C* Given a Table D sequence mnemonic (i.e. a "parent mnemonic") and an	*
C* integer IENT, this subroutine returns the bit-wise representation of	*
C* the FXY value corresponding to, sequentially, the (IENT)th child	*
C* mnemonic of that parent mnemonic.					*
C*									*
C* UPTDD  ( ID, LUN, IENT, IRET )					*
C*									*
C* Input parameters:							*
C*	ID		INTEGER		Positional index of parent	*
C*					mnemonic within internal BUFR	*
C*					Table D array TABD(*,*)		*
C*	LUN		INTEGER		I/O stream index of parent	*
C*					mnemonic within internal BUFR	*
C*					Table D array TABD(*,*)		*
C*	IENT		INTEGER		Ordinal indicator of child	*
C*					mnemonic to return from within	*
C*					TABD(ID,LUN) sequence:		*
C*					  0 = return a count of the	*
C*					      total number of child	*
C*					      mnemonics within		*
C*					      TABD(ID,LUN)		*
C* Output parameters:							*
C*	IRET		INTEGER		Return value:			*
C*				      /               \			*
C*	      ------------------------                 ---------------	*
C* 		The interpretation of the return value IRET depends	*
C*		upon the input value IENT, as follows:			*
C*									*
C*		IF ( IENT = 0 ) THEN					*
C*		    IRET = a count of the total number of child		*
C*			   mnemonics within TABD(ID,LUN)		*
C*		ELSE							*
C*		    IRET = the bit-wise representation of the FXY value	*
C*			   corresponding to the (IENT)th child mnemonic	*
C*			   of TABD(ID,LUN)				*
C*		END IF							*
C*	      --------------------------------------------------------	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*56  DXSTR                                               
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      LDD = LDXD(IDXV+1)+1                                              
                                                                        
C  CHECK IF IENT IS IN BOUNDS                                           
C  --------------------------                                           
                                                                        
      NDSC = IUPM(TABD(ID,LUN)(LDD:LDD),8)                              
                                                                        
      IF(IENT.EQ.0) THEN                                                
         IRET = NDSC                                                    
         RETURN                                                         
      ELSEIF(IENT.LT.0 .OR. IENT.GT.NDSC) THEN                          
         CALL BORT('UPTDD - IENT OUT OF RANGE')                        
      ENDIF                                                             
                                                                        
C  RETURN THE DESCRIPTOR INDICATED BY IENT                              
C  ---------------------------------------                              
                                                                        
      IDSC = LDD+1 + (IENT-1)*2                                         
      IRET = IUPM(TABD(ID,LUN)(IDSC:IDSC),16)                           
                                                                        
      RETURN                                                            
      END                                                               

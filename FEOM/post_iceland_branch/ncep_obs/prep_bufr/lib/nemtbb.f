      SUBROUTINE NEMTBB(LUN,ITAB,UNIT,ISCL,IREF,IBIT)                   

C************************************************************************
C* NEMTBB								*
C*									*
C* This subroutine checks all of the properties (e.g. FXY value, units,	*
C* scale factor, reference value, etc.) of a specified mnemonic within	*
C* the internal BUFR Table B arrays in order to verify that the values	*
C* of those properties are all legal and well-defined.  If any errors	*
C* are found, then an appropriate call is made to subroutine BORT.	*
C*									*
C* NEMTBB  ( LUN, ITAB, UNIT, ISCL, IREF, IBIT )			*
C*									*
C* Input parameters:							*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					BUFR Table B arrays for mnemonic*
C*					to be checked			*
C*	ITAB		INTEGER		Positional index into internal	*
C*					BUFR Table B arrays for mnemonic*
C*					to be checked			*
C*									*
C* Output parameters:							*
C*	UNIT		CHARACTER*(*)	Units of mnemonic		*
C*	ISCL		INTEGER		Scale factor of mnemonic	*
C*	IREF		INTEGER		Reference value of mnemonic	*
C*	IBIT		INTEGER		Bit width of mnemonic		*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*24  UNIT                                                
      CHARACTER*8   NEMO                                                
      REAL*8        MXR                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      MXR = 1E11-1                                                      
                                                                        
      IF(ITAB.LE.0 .OR. ITAB.GT.NTBB(LUN)) GOTO 900                     
                                                                        
C  PULL OUT TABLE B INFORMATION                                         
C  ----------------------------                                         
                                                                        
      IDN  = IDNB(ITAB,LUN)                                             
      NEMO = TABB(ITAB,LUN)( 7:14)                                      
      UNIT = TABB(ITAB,LUN)(71:94)                                      
      ISCL = VALX(TABB(ITAB,LUN)( 95: 98))                              
      IREF = VALX(TABB(ITAB,LUN)( 99:109))                              
      IBIT = VALX(TABB(ITAB,LUN)(110:112))                              
                                                                        
C  CHECK TABLE B CONTENTS                                               
C  ----------------------                                               
                                                                        
      IF(IDN.LT.IFXY('000000')) GOTO 901                                
      IF(IDN.GT.IFXY('063255')) GOTO 901                                
                                                                        
      IF(ISCL.LT.-999 .OR. ISCL.GT.999) GOTO 902                        
      IF(IREF.LE.-MXR .OR. IREF.GE.MXR) GOTO 903                        
      IF(IBIT.LE.0) GOTO 904 
      IF(UNIT(1:5).NE.'CCITT' .AND. IBIT.GT.64      ) GOTO 904 
      IF(UNIT(1:5).EQ.'CCITT' .AND. MOD(IBIT,8).NE.0) GOTO 905
                                                                        
      RETURN                                                            
900   CALL BORT('NEMTBB - ITAB NOT IN TABLE B'         )               
901   CALL BORT('NEMTBB - BAD DESCRIPTOR VALUE: '//NEMO)               
902   CALL BORT('NEMTBB - BAD SCALE VALUE     : '//NEMO)               
903   CALL BORT('NEMTBB - BAD REFERENCE VALUE : '//NEMO)               
904   CALL BORT('NEMTBB - BAD BIT WIDTH       : '//NEMO)               
905   CALL BORT('NEMTBB - BAD CHAR BIT WIDTH  : '//NEMO)               
      END                                                               

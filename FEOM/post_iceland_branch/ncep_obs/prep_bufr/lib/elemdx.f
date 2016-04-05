      SUBROUTINE ELEMDX(CARD,LUN)                                       

C************************************************************************
C* ELEMDX								*
C*									*
C* This subroutine decodes the scale factor, reference value, bit width,*
C* and units from a mnemonic definition card that was previously read	*
C* from a user DX table by subroutine RDUSDX.  These decoded values are	*
C* then added to the already-existing entry for that mnemonic within the*
C* internal BUFR table B array TABB(*,LUN).				*
C*									*
C* ELEMDX  ( CARD, LUN )						*
C*									*
C* Input parameters:							*
C*	CARD		CHARACTER*80	Mnemonic definition card that	*
C*					was read from a user DX table	*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					arrays for this user DX table	*
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
      CHARACTER*80  CARD                                                
      CHARACTER*24  UNIT                                                
      CHARACTER*11  REFR                                                
      CHARACTER*8   NEMO                                                
      CHARACTER*4   SCAL                                                
      CHARACTER*3   BITW                                                
      CHARACTER*1   SIGN,TAB                                            
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CAPTURE THE VARIOUS ELEMENTS CHARACTERISTICS                         
C  --------------------------------------------                         
                                                                        
      NEMO = CARD( 3:10)                                                
      SCAL = CARD(14:17)                                                
      REFR = CARD(21:31)                                                
      BITW = CARD(35:37)                                                
      UNIT = CARD(41:64)                                                
                                                                        
C  FIND THE ELEMENT TAG IN TABLE B                                      
C  -------------------------------                                      

C     Note that an entry for this mnemonic should already exist within
C     the internal BUFR table B array TABB(*,LUN); this entry should
C     have been created by subroutine RDUSDX when the mnemonic and its
C     associated FXY value and description were initially defined within
C     a card read from near the top of the user DX table.  Now, we need
C     to retrieve the positional index for that entry within TABB(*,LUN)
C     so that we can access the entry and then add the scale factor,
C     reference value, bit width, and units to it.
                                                                        
      CALL NEMTAB(LUN,NEMO,IDSN,TAB,IELE)                               
      IF(TAB.NE.'B') GOTO 900                                           
                                                                        
C  LEFT JUSTIFY AND STORE CHARACTERISTICS                               
C  --------------------------------------                               
                                                                        
      CALL JSTCHR(UNIT)                                                 
      TABB(IELE,LUN)(71:94) = UNIT                                      
                                                                        
      CALL JSTNUM(SCAL,SIGN,IRET)                                       
      IF(IRET.NE.0) GOTO 901                                            
      TABB(IELE,LUN)(95:95) = SIGN                                      
      TABB(IELE,LUN)(96:98) = SCAL                                      
                                                                        
      CALL JSTNUM(REFR,SIGN,IRET)                                       
      IF(IRET.NE.0) GOTO 902                                            
      TABB(IELE,LUN)( 99: 99) = SIGN                                    
      TABB(IELE,LUN)(100:109) = REFR                                    
                                                                        
      CALL JSTNUM(BITW,SIGN,IRET)                                       
      IF(IRET.NE.0  ) GOTO 903                                          
      IF(SIGN.EQ.'-') GOTO 903                                          
      TABB(IELE,LUN)(110:112) = BITW                                    
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXIT                                                           
C  -----------                                                          
                                                                        
900   CALL BORT('ELEMDX - UNDEFINED ELEMENT: '//CARD)                  
901   CALL BORT('ELEMDX - BAD SCALE VALUE:   '//CARD)                  
902   CALL BORT('ELEMDX - BAD REFERENCE VAL: '//CARD)                  
903   CALL BORT('ELEMDX - BAD BIT WIDTH:     '//CARD)                  
      END                                                               

      SUBROUTINE RDUSDX(LUNDX,LUN)                                      

C************************************************************************
C* RDUSDX								*
C*									*
C* This subroutine reads and parses a file containing a user DX table	*
C* and then stores this information into the following internal arrays:	*
C*									*
C* For Table A entries:							*
C*	NTBA(LUN)	INTEGER		Number of Table A entries	*
C*					(note that NTBA(0) contains the	*
C*					maximum number of such entries	*
C*					as set within subroutine BFRINI)*
C*	TABA(N,LUN)	CHARACTER*128	Table A entries, where 		*
C*					N=1,2,3,...,NTBA(LUN)		*
C*	IDNA(N,LUN,1)	INTEGER		Message type corresponding	*
C*					to TABA(N,LUN)			*
C*	IDNA(N,LUN,2)	INTEGER		Message subtype corresponding	*
C*					to TABA(N,LUN)			*
C*									*
C* For Table B entries:							*
C*	NTBB(LUN)	INTEGER		Number of Table B entries	*
C*					(note that NTBB(0) contains the	*
C*					maximum number of such entries	*
C*					as set within subroutine BFRINI)*
C*	TABB(N,LUN)	CHARACTER*128	Table B entries, where 		*
C*					N=1,2,3,...,NTBB(LUN)		*
C*	IDNB(N,LUN)	INTEGER		Bit-wise representation of the	*
C*					FXY value corresponding to	*
C*					TABB(N,LUN)			*
C*									*
C* For Table D entries:							*
C*	NTBD(LUN)	INTEGER		Number of Table D entries	*
C*					(note that NTBD(0) contains the	*
C*					maximum number of such entries	*
C*					as set within subroutine BFRINI)*
C*	TABD(N,LUN)	CHARACTER*600	Table D entries, where 		*
C*					N=1,2,3,...,NTBD(LUN)		*
C*	IDND(N,LUN)	INTEGER		Bit-wise representation of the	*
C*					FXY value corresponding to	*
C*					TABD(N,LUN)			*
C*									*
C* RDUSDX  ( LUNDX, LUN )						*
C*									*
C* Input parameters:							*
C*	LUNDX		INTEGER		FORTRAN logical unit # of file	*
C*					containing user DX table	*
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
      CHARACTER*8   NEMO                                  
      CHARACTER*6   NUMB                                                
      LOGICAL       DIGIT
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  INITIALIZE THE DX-TABLE PARTITION AND SOURCE FILE                    
C  -------------------------------------------------                    
                                                                        
      CALL DXINIT(LUN,1)                                                
      REWIND LUNDX                                                      
                                                                        
C  READ USER CARDS UNTIL THERE ARE NO MORE                              
C  ---------------------------------------                              
                                                                        
1     READ(LUNDX,'(A80)',END=100) CARD                                  
                                                                        
C  REREAD IF NOT A DEFINITION CARD                                      
C  -------------------------------                                      
                                                                        
      IF(CARD(1: 1).EQ.       '*') GOTO 1                               
      IF(CARD(3:10).EQ.'--------') GOTO 1                               
      IF(CARD(3:10).EQ.'        ') GOTO 1                               
      IF(CARD(3:10).EQ.'MNEMONIC') GOTO 1                               
      IF(CARD(3:10).EQ.'TABLE  D') GOTO 1                               
      IF(CARD(3:10).EQ.'TABLE  B') GOTO 1                               
                                                                        
C  PARSE A DESCRIPTOR DEFINITION CARD                                   
C  ----------------------------------                                   
                                                                        
      IF(CARD(12:12).EQ.'|' .AND. CARD(21:21).EQ.'|') THEN              
                                                                        
         NEMO = CARD(3:10)                                              

C	 NEMO is the CHARACTER*8 mnemonic name.

         NUMB = CARD(14:19)                                             

C	 NUMB is the CHARACTER*6 FXY value corresponding to NEMO.

         IF(NEMOCK(NEMO).NE.0) GOTO 900                                 
         IF(NUMBCK(NUMB).NE.0) GOTO 900                                 
                                                                        
                                                                        
         IF(NUMB(1:1).EQ.'A') THEN                                      

C	    This is a Table A mnemonic.

            N = NTBA(LUN)+1                                             
            IF(N.GT.NTBA(0)) GOTO 901                                   
            CALL NENUAA(NEMO,NUMB,LUN)                             
            TABA(N,LUN)( 1: 3) = NUMB(4:6)                              
            TABA(N,LUN)( 4:11) = NEMO                                   
            TABA(N,LUN)(13:67) = CARD(23:77)                            
            NTBA(LUN) = N                                               
                                                                        
            IF(DIGIT(NEMO(3:8))) THEN

C	       Set the message type and subtype.

               READ(NEMO,'(2X,2I3)') MTYP,MSBT                            
               IDNA(N,LUN,1) = MTYP                                     
               IDNA(N,LUN,2) = MSBT              
            ELSE

C	       If either the message type or subtype could not be determined,
C	       then set default values for both.

               READ(NUMB(4:6),'(I3)') IDNA(N,LUN,1)                           
               IDNA(N,LUN,2) = 0                                              
            ENDIF                                                             
                                                                        
            NUMB(1:1) = '3'                                             
         ENDIF                                                          
                                                                        
                                                                        
         IF(NUMB(1:1).EQ.'0') THEN                                      

C	    This is a Table B mnemonic.

            N = NTBB(LUN)+1                                             
            IF(N.GT.NTBB(0)) GOTO 902                                   
            CALL NENUBD(NEMO,NUMB,LUN)                                  
            IDNB(N,LUN) = IFXY(NUMB)                                    
            TABB(N,LUN)( 1: 6) = NUMB                                   
            TABB(N,LUN)( 7:14) = NEMO                                   
            TABB(N,LUN)(16:70) = CARD(23:77)                            
            NTBB(LUN) = N                                               
            GOTO 1                                                      
         ENDIF                                                          
                                                                        
                                                                        
         IF(NUMB(1:1).EQ.'3') THEN                                      

C	    This is a Table D mnemonic.

            N = NTBD(LUN)+1                                             
            IF(N.GT.NTBD(0)) GOTO 903                                   
            CALL NENUBD(NEMO,NUMB,LUN)                                  
            IDND(N,LUN) = IFXY(NUMB)                                    
            TABD(N,LUN)( 1: 6) = NUMB                                   
            TABD(N,LUN)( 7:14) = NEMO                                   
            TABD(N,LUN)(16:70) = CARD(23:77)                            
            NTBD(LUN) = N                                               
            GOTO 1                                                      
         ENDIF                                                          
                                                                        
         GOTO 904                                                       
      ENDIF                                                             
                                                                        
C  PARSE A SEQUENCE DEFINITION CARD                                     
C  --------------------------------                                     
                                                                        
      IF(CARD(12:12).EQ.'|' .AND. CARD(19:19).NE.'|') THEN              
         CALL SEQSDX(CARD,LUN)                                          
         GOTO 1                                                         
      ENDIF                                                             
                                                                        
C  PARSE AN ELEMENT DEFINITION CARD                                     
C  --------------------------------                                     
                                                                        
      IF(CARD(12:12).EQ.'|' .AND. CARD(19:19).EQ.'|') THEN              
         CALL ELEMDX(CARD,LUN)                                          
         GOTO 1                                                         
      ENDIF                                                             
                                                                        
C  CAN'T FIGURE OUT WHAT KIND OF CARD IT IS                             
C  ----------------------------------------                             
                                                                        
      GOTO 905                                                          
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
100   CALL MAKESTAB                                                     
      RETURN                                                            
                                                                        
C  ERROR EXIT                                                           
C  ----------                                                           
                                                                        
900   PRINT*,CARD                                                       
      CALL BORT('RDUSDX - NEMO OR NUMB ERROR             '//CARD)      
901   CALL BORT('RDUSDX - TOO MANY TABLE A ENTRIES       '//CARD)      
902   CALL BORT('RDUSDX - TOO MANY TABLE B ENTRIES       '//CARD)      
903   CALL BORT('RDUSDX - TOO MANY TABLE D ENTRIES       '//CARD)      
904   CALL BORT('RDUSDX - BAD DESCRIPTOR NUMBER          '//CARD)      
905   CALL BORT('RDUSDX - BAD CARD FORMAT                '//CARD)      
      END                                                               

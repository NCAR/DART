      SUBROUTINE NENUCK(NEMO,NUMB,LUN)                                  

C************************************************************************
C* NENUCK								*
C*									*
C* This subroutine consists solely of two separate entry points for	*
C* checking a mnemonic and FXY value pair that were read from a user	*
C* DX table, in order to make sure that neither value has already been	*
C* defined within the internal BUFR table arrays for the given LUN.	*
C*									*
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
      CHARACTER*8   NEMO                                                
      CHARACTER*6   NUMB                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK TABLE A                                                        
C  -------------                                                        
                                                                        
      ENTRY NENUAA(NEMO,NUMB,LUN)                                       
C************************************************************************
C* NENUAA								*
C*									*
C* This entry point checks a mnemonic and FXY value pair that were read	*
C* from a user DX table, in order to make sure that neither value has	*
C* already been defined within internal BUFR table A for the given LUN.	*
C* If either value has already been defined for this LUN, then an	*
C* appropriate call is made to subroutine BORT.				*
C*									*
C* NENUAA  ( NEMO, NUMB, LUN )						*
C*									*
C* Input parameters:							*
C*	NEMO		CHARACTER*(*)	Mnemonic 			*
C*	NUMB		INTEGER		FXY value associated with NEMO	*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					arrays for this user DX table	*
C************************************************************************
C*									*
C*	EXAMPLE SHOWING LAYOUT OF INTERNAL BUFR TABLE A			*
C* (FROM A DBX DEBUG SESSION USING "bufrtab.002", AND WHERE LUN = 1)	*
C*									*
C* (dbx) print NTBA[1]							*
C* 8									*
C*									*
C* (dbx) print TABA[1,1]						*
C* 0x1002c764 = "218NC002001 MESSAGE TYPE 002-001  RAWINSONDE - FIXED LAND                                                                       "	*
C* (dbx) print TABA[2,1]						*
C* 0x1002c7e4 = "219NC002002 MESSAGE TYPE 002-002  RAWINSONDE - MOBIL LAND                                                                       "	* 
C* (dbx) print TABA[3,1]						*
C* 0x1002c864 = "220NC002003 MESSAGE TYPE 002-003  RAWINSONDE - SHIP                                                                             "	* 
C*									*
C*	and so on, up through TABA[8,1] ( = TABA[NTBA[LUN],LUN] )	*
C*									*
C************************************************************************
                                                                        
      DO N=1,NTBA(LUN)                                                  
      IF(NUMB(4:6).EQ.TABA(N,LUN)(1: 3)) GOTO 900                       
      IF(NEMO     .EQ.TABA(N,LUN)(4:11)) GOTO 900                       
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  CHECK TABLE B AND D                                                  
C  -------------------                                                  
                                                                        
      ENTRY NENUBD(NEMO,NUMB,LUN)                                       
C************************************************************************
C* NENUBD								*
C*									*
C* This entry point checks a mnemonic and FXY value pair that were read	*
C* from a user DX table, in order to make sure that neither value has	*
C* already been defined within internal BUFR table B or D for the given	*
C* LUN.  If either value has already been defined for this LUN, then an	*
C* appropriate call is issued to subroutine BORT.			*
C*									*
C* NENUBD  ( NEMO, NUMB, LUN )						*
C*									*
C* Input parameters:							*
C*	NEMO		CHARACTER*(*)	Mnemonic 			*
C*	NUMB		INTEGER		FXY value associated with NEMO	*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					arrays for this user DX table	*
C************************************************************************
C*									*
C*	EXAMPLE SHOWING LAYOUT OF INTERNAL BUFR TABLE B			*
C* (FROM A DBX DEBUG SESSION USING "bufrtab.002", AND WHERE LUN = 1)	*
C*									*
C* (dbx) print NTBB[1]							*
C* 95									* 
C*									*
C* (dbx) print TABB[1,1]						*
C* 0x1003c164 = "063000BYTCNT                                                          BYTES                   +0  +0         16                 " 	*
C*									*
C* (dbx) print TABB[2,1]						*
C* 0x1003c1e4 = "063255BITPAD                                                          NONE                    +0  +0         1                  "	* 
C*									*
C* (dbx) print TABB[3,1]						*
C* 0x1003c264 = "031000DRF1BIT                                                         NUMERIC                 +0  +0         1                  "	* 
C*									*
C* (dbx) print TABB[8,1]						*
C* 0x1003c4e4 = "001003WMOR     WMO REGION NUMBER                                      CODE TABLE              +0  +0         3                  "	* 
C*									*
C* (dbx) print TABB[11,1]						*
C* 0x1003c664 = "001194BUHD     BULLETIN HEADER                                        CCITT IA5               +0  +0         64                 "	* 
C*									*
C* (dbx) print TABB[21,1]						*
C* 0x1003cb64 = "004003DAYS     DAY                                                    DAY                     +0  +0         6                  "	* 
C*									*
C* (dbx) print TABB[33,1]						*
C* 0x1003d164 = "005002CLAT     LATITUDE (COARSE ACCURACY)                             DEGREES                 +2  -9000      15                 "	* 
C*									*
C*									*
C*	and so on, up through TABB[95,1] ( = TABB[NTBB[LUN],LUN] )	*
C*									*
C************************************************************************
C*									*
C*	EXAMPLE SHOWING LAYOUT OF INTERNAL BUFR TABLE D			*
C* (FROM A DBX DEBUG SESSION USING "bufrtab.002", AND WHERE LUN = 1)	*
C*									*
C* (dbx) print NTBD[1]							*
C* 43									*
C*									*
C* (dbx) &TABD[1,1]/14c							*
C* 1008a364:  '3' '6' '0' '0' '0' '1' 'D' 'R' 'P' '1' '6' 'B' 'I' 'T'	*
C*									*
C* (dbx) &TABD[2,1]/14c							*
C* 1008a5bc:  '3' '6' '0' '0' '0' '2' 'D' 'R' 'P' '8' 'B' 'I' 'T' ' '	*
C*									*
C* (dbx) &TABD[3,1]/14c							*
C* 1008a814:  '3' '6' '0' '0' '0' '3' 'D' 'R' 'P' 'S' 'T' 'A' 'K' ' '	*
C*									*
C* (dbx) &TABD[4,1]/14c							*
C* 1008aa6c:  '3' '6' '0' '0' '0' '4' 'D' 'R' 'P' '1' 'B' 'I' 'T' ' '	*
C*									*
C* (dbx) &TABD[5,1]/14c							*
C* 1008acc4:  '3' '6' '3' '2' '1' '8' 'N' 'C' '0' '0' '2' '0' '0' '1'	*
C*									*
C* (dbx) &TABD[6,1]/14c							*
C* 1008af1c:  '3' '6' '3' '2' '1' '9' 'N' 'C' '0' '0' '2' '0' '0' '2'	*
C*									*
C* (dbx) &TABD[24,1]/14c						*
C* 1008d94c:  '3' '6' '1' '1' '3' '0' 'U' 'A' 'A' 'D' 'F' ' ' ' ' ' '	*
C*									*
C*	and so on, up through TABD[43,1] ( = TABD[NTBD[LUN],LUN] )	*
C*									*
C************************************************************************
                                                                        
      DO N=1,NTBB(LUN)                                                  
      IF(NUMB.EQ.TABB(N,LUN)(1: 6)) GOTO 900                            
      IF(NEMO.EQ.TABB(N,LUN)(7:14)) GOTO 900                            
      ENDDO                                                             
                                                                        
      DO N=1,NTBD(LUN)                                                  
      IF(NUMB.EQ.TABD(N,LUN)(1: 6)) GOTO 900                            
      IF(NEMO.EQ.TABD(N,LUN)(7:14)) GOTO 900                            
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXIT                                                           
C  ----------                                                           
                                                                        
900   CALL BORT('NENUCK - DUPLICATE NEM/NUM '//NEMO//' '//NUMB)        
      END                                                               

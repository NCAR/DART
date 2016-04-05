      SUBROUTINE OPENBF(LUNIT,IO,LUNDX)                                 

C************************************************************************
C* OPENBF								*
C*									*
C* This subroutine is called in order to identify a new logical unit to	*
C* the BUFRLIB software for input or output operations.  Up to 32 such	*
C* logical units can be connected to the software at any one time.	*
C*									*
C* OPENBF  ( LUNIT, IO, LUNDX )						*
C*									*
C* Input parameters:							*
C*	LUNIT		INTEGER		FORTRAN logical unit number	*
C*	IO		CHARACTER*(*)	Flag indicating how LUNIT is	*
C*					to be used by the software:	*
C*					  'IN' = input operations	*
C*					  'OUT' = output operations	*
C*					  'APN' = same as 'OUT', except	*
C*						  begin writing at end	*
C*						  of file ("append")	*
C*					  'NUL' = same as 'OUT', except	*
C*						  don't actually write	*
C*						  to LUNIT (for use with*
C*						  subroutine WRITSA)	*
C*	LUNDX		INTEGER		FORTRAN logical unit number	*
C*					containing DX table information	*
C*					to be used in reading/writing	*
C*					from/to LUNIT (depending on the	*
C*					case); may be set equal to	*
C*					LUNIT if DX table information	*
C*					is already embedded in LUNIT	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C* J. Ator/NCEP		06/01	Added IO='NUL' option			*
C************************************************************************
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /STBFR / IOLUN(32),IOMSG(32)                               
      COMMON /NULBFR/ NULL(32)
      COMMON /QUIET / IPRT                                              
                                                                        
      CHARACTER*(*) IO                                                  
      CHARACTER*4   BUFR,MSTR                                           
      LOGICAL       SKIPDX,APPEND                                       
                                                                        
      DATA IFIRST/0/                                                    
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IF(IFIRST.EQ.0) THEN                                              
         CALL WRDLEN                                                    
         CALL BFRINI                                                    
         IFIRST = 1                                                     
      ENDIF                                                             
                                                                        
      IF(IO.EQ.'QUIET') THEN                                            
         IPRT = LUNDX                                                   
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  SEE IF A FILE CAN BE OPENED                                          
C  ---------------------------                                          
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(LUN.EQ.0) GOTO 900                                             
      IF(IL .NE.0) GOTO 901                                             
      NULL(LUN) = 0
                                                                        
C  CHECK FOR NO BUFR DATA OR NO DATA AT ALL IN AN "IN" FILE             
C  --------------------------------------------------------             

C  Note that we only want to do this check if LUNIT=LUNDX, in order 
C  to allow for the possibility that input BUFR messages may be passed
C  to the BUFRLIB software via an alternative method (e.g. a future
C  call to subroutine READERME) rather than read directly from LUNIT,
C  which is the usual method.
                                                                        
      IF(IO.EQ.'IN' .AND. LUNIT.EQ.LUNDX) THEN                          
         REWIND LUNIT                                                   
         READ(LUNIT,END=100,ERR=902) MSTR                               
C
C        Note that we need to pass MSTR through subroutine UPC before
C        comparing it to the string 'BUFR' because the native machine
C        might be EBCDIC rather than ASCII; in this case MSTR, which
C        contains the first 4 bytes of the first BUFR message in LUNIT
C        (and which thus must always be CCITT IA5 (i.e. ASCII) by
C        definition, regardless of the native machine!), will be
C        converted to EBCDIC, and thus the subsequent comparison
C        with the string 'BUFR' will always be valid.
C
         IBIT = 0                                                       
         CALL UPC(BUFR,4,MSTR,IBIT)                                     
         IF(BUFR.NE.'BUFR') GOTO 902                                    
      ENDIF                                                             
                                                                        
C  SET INITIAL OPEN DEFAULTS                                            
C  -------------------------                                            
                                                                        
      IF(IO.NE.'NUL') THEN
        REWIND LUNIT                                                      
      ENDIF
      NMSG (LUN) = 0                                                    
      NSUB (LUN) = 0                                                    
      MSUB (LUN) = 0                                                    
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
      SKIPDX = .FALSE.                                                  
      APPEND = .FALSE.                                                  
                                                                        
C  DECIDE HOW TO SETUP THE DICTIONARY                                   
C  ----------------------------------                                   
                                                                        
      IF(IO.EQ.'IN') THEN                                               
         CALL WTSTAT(LUNIT,LUN,-1,0)                                    
         CALL READDX(LUNIT,LUN,LUNDX)                                   
      ELSE IF(IO.EQ.'OUT') THEN                                         
         CALL WTSTAT(LUNIT,LUN, 1,0)                                    
         CALL WRITDX(LUNIT,LUN,LUNDX)                                   
      ELSE IF(IO.EQ.'NUL') THEN
         CALL WTSTAT(LUNIT,LUN, 1,0)
         CALL READDX(LUNIT,LUN,LUNDX)                                   
         NULL(LUN) = 1
      ELSE IF(IO.EQ.'APN' .OR. IO.EQ.'APX') THEN                        
         CALL WTSTAT(LUNIT,LUN, 1,0)                                    
         CALL READDX(LUNIT,LUN,LUNDX)                                   
         IF(IO.EQ.'APN') CALL POSAPN(LUNIT)                             
         IF(IO.EQ.'APX') CALL POSAPX(LUNIT)                             
      ELSE                                                              
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
                                                                        
C  FILE OPENED FOR INPUT IS EMPTY - LET READMG GIVE THE BAD NEWS        
C  -------------------------------------------------------------        
                                                                        
100   REWIND LUNIT                                                      
      CALL WTSTAT(LUNIT,LUN,-1,0)                                       
      CALL DXINIT(LUN,0)                                                
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('OPENBF - TOO MANY FILES OPENED ALREADY       ')       
901   CALL BORT('OPENBF - FILE ALREADY OPEN                   ')       
902   CALL BORT('OPENBF - INPUT FILE HAS NON-BUFR DATA        ')       
903   CALL BORT('OPENBF - IO MUST BE ONE OF "IN" "OUT" "APN"  ')       
      END                                                               

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DATEBF(LUNIT,MEAR,MMON,MDAY,MOUR,IDATE)                        
                                                                        
      COMMON /DATELN/ LENDAT

      CHARACTER*26  MSTR                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IDATE = -1                                                        
                                                                        
C  SEE IF THE FILE IS ALREADY OPEN TO BUFR INTERFACE (A NO-NO)          
C  -----------------------------------------------------------          
                                                                        
      CALL STATUS(LUNIT,LUN,JL,JM)                                      
      IF(JL.NE.0) CALL BORT('DATEBF - FILE ALREADY OPEN')              
                                                                        
C  CHECK FOR NO BUFR DATA OR NO DATA AT ALL                             
C  ----------------------------------------                             
                                                                        
      REWIND LUNIT                                                      
      READ(LUNIT,END=100,ERR=100) MSTR                                  
      IF(MSTR(1:4).NE.'BUFR') GOTO 100                                  
                                                                        
C  READ TO A DATA MESSAGE AND PICK OUT THE DATE                         
C  --------------------------------------------                         
                                                                        
1     READ(LUNIT,END=100,ERR=100) MSTR                                  
      IF(ICHAR(MSTR(17:17)).EQ.11) GOTO 1                               
      MEAR = MOD(ICHAR(MSTR(21:21)),100)                                      
      MMON = ICHAR(MSTR(22:22))                                           
      MDAY = ICHAR(MSTR(23:23))                                           
      MOUR = ICHAR(MSTR(24:24))                                           
      MMIN = ICHAR(MSTR(25:25))                                           
      MCEN = MAX(0,ICHAR(MSTR(26:26))-MIN(MEAR,1))                              
                                                                        
      IF(LENDAT.EQ.10) THEN
         IDATE = MCEN*10**8+MEAR*10**6+MMON*10**4+MDAY*10**2+MOUR               
         IDATE = I4DY(IDATE)
         MEAR  = IDATE/10**6
      ELSE
         IDATE = MEAR*10**6+MMON*10**4+MDAY*10**2+MOUR               
      ENDIF

100   RETURN                                                            
      END                                                               

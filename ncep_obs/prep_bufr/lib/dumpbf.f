C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DUMPBF(LUNIT,JDATE,JDUMP)                              
                                                                        
      COMMON /DATELN/ LENDAT

      DIMENSION JDATE(5),JDUMP(5)                                       
                                                                        
      CHARACTER*32  MSTR                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      DO I=1,5                                                          
      JDATE(I) = -1                                                     
      JDUMP(I) = -1                                                     
      ENDDO                                                             
                                                                        
C  SEE IF THE FILE IS ALREADY OPEN TO BUFR INTERFACE (A NO-NO)          
C  -----------------------------------------------------------          
                                                                        
      CALL STATUS(LUNIT,LUN,JL,JM)                                      
      IF(JL.NE.0) CALL BORT('DUMPBF - FILE ALREADY OPEN')              
                                                                        
C  CHECK FOR NO BUFR DATA OR NO DATA AT ALL                             
C  ----------------------------------------                             
                                                                        
      REWIND LUNIT                                                      
1     READ(LUNIT,END=100,ERR=100) MSTR                                  
      IF(MSTR(1:4).NE.'BUFR') GOTO 100                                  
      IF(ICHAR(MSTR(17:17)).EQ.11) GOTO 1                               
                                                                        
C  DUMP CENTER YY,MM,DD,HH,MM IS IN THE FIRST EMPTY MESSAGE             
C  --------------------------------------------------------             
                                                                        
      IF(ICHAR(MSTR(31:31)).EQ.0 .AND. ICHAR(MSTR(32:32)).EQ.0) THEN    
         JDATE(1) = MOD(ICHAR(MSTR(21:21)),100)
         JDATE(2) = ICHAR(MSTR(22:22))                                  
         JDATE(3) = ICHAR(MSTR(23:23))                                  
         JDATE(4) = ICHAR(MSTR(24:24))                                  
         JDATE(5) = ICHAR(MSTR(25:25))                                  
         MCEN     = MAX(0,ICHAR(MSTR(26:26))-MIN(JDATE(1),1))
      ELSE                                                              
         RETURN                                                         
      ENDIF                                                             

      IF(LENDAT.EQ.10) THEN
         JDATE(1) = I4DY(MCEN*10**8+JDATE(1)*10**6)/10**6
      ENDIF
                                                                        
C  DUMP CLOCK YY,MM,DD,HH,MM IS IN THE SECOND EMPTY MESSAGE             
C  --------------------------------------------------------             
                                                                        
      READ(LUNIT,END=100,ERR=100) MSTR                                  
                                                                        
      IF(ICHAR(MSTR(31:31)).EQ.0 .AND. ICHAR(MSTR(32:32)).EQ.0) THEN    
         JDUMP(1) = MOD(ICHAR(MSTR(21:21)),100)
         JDUMP(2) = ICHAR(MSTR(22:22))                                  
         JDUMP(3) = ICHAR(MSTR(23:23))                                  
         JDUMP(4) = ICHAR(MSTR(24:24))                                  
         JDUMP(5) = ICHAR(MSTR(25:25))                                  
         MCEN     = MAX(0,ICHAR(MSTR(26:26))-MIN(JDUMP(1),1))
      ELSE                                                              
         RETURN                                                         
      ENDIF                                                             
                                                                        
      IF(LENDAT.EQ.10) THEN
         JDUMP(1) = I4DY(MCEN*10**8+JDUMP(1)*10**6)/10**6
      ENDIF
                                                                        
100   RETURN                                                            
      END

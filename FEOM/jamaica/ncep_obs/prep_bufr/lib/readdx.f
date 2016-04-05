C-----------------------------------------------------------------------
C  BUFR TABLE INFORMATION CONTAINED IN A FILE CONNECTED TO UNIT LUNDX   
C  IS USED TO INITIALIZE PROCESSING TABLES FOR A BUFR FILE CONNECTED    
C  TO UNIT LUNIT. LUNIT AND LUNDX MAY BE THE SAME ONLY IF THE UNIT IS   
C  CONNECTED TO A BUFR FILE, CURRENTLY OPEN FOR INPUT PROCESSING,       
C  POSITIONED AT A DX-TABLE MESSAGE (ANYWHERE IN THE FILE). OTHERWISE,  
C  LUNDX MAY BE CONNECTED TO ANOTHER CURRENTLY OPEN AND DEFINED BUFR    
C  FILE, OR TO A USER SUPPLIED, CHARACTER FORMAT, DX-TABLE FILE.        
C                                                                       
C  NOTE: READDX IS USED TO INITIALIZE INTERNAL BUFR DX-TABLES ONLY.     
C        IF A BUFR OUTPUT FILE IS BEING OPENED, SUBROUTINE WRITDX       
C        CALLS READDX TO INITIALIZE THE INTERNAL DX-TABLES, AND THEN    
C        WRITES BUFR DX-TABLE MESSAGES INTO THE OUTPUT FILE.            
C                                                                       
C                                                                       
C  INPUT ARGUMENTS:                                                     
C     LUNIT    - UNIT CONNECTED TO BUFR FILE TO BE INITIALIZED/UPDATED  
C     LUN      - INTERNAL BUFR UNIT ASSOCIATED WITH FORTRAN UNIT LUNIT  
C     LUNDX    - UNIT CONTAINING DX-TABLES                              
C                                                                       
C-----------------------------------------------------------------------
      SUBROUTINE READDX(LUNIT,LUN,LUNDX)                                
                                                                        
      COMMON /QUIET/ IPRT                                               
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  GET THE BUFR STATUS OF UNIT LUNDX                                    
C  ---------------------------------                                    
                                                                        
      CALL STATUS(LUNDX,LUD,ILDX,IMDX)                                  
                                                                        
C  READ A DX-TABLE FROM THE INDICATED SOURCE                            
C  -----------------------------------------                            
                                                                        
      IF (LUNIT.EQ.LUNDX) THEN                                          
         IF(IPRT.GE.1) PRINT100,LUNDX,LUNIT                             
         REWIND LUNIT                                                   
         CALL RDBFDX(LUNIT,LUN)                                         
      ELSEIF(ILDX.NE.0) THEN                                            
         IF(IPRT.GE.1) PRINT101,LUNDX,LUNIT                             
         CALL CPBFDX(LUD,LUN)                                           
      ELSEIF(ILDX.EQ.0) THEN                                            
         IF(IPRT.GE.1) PRINT102,LUNDX,LUNIT                             
         REWIND LUNDX                                                   
         CALL RDUSDX(LUNDX,LUN)                                         
      ELSE                                                              
         CALL BORT('READDX - SCREWUP')                                 
      ENDIF                                                             
                                                                        
100   FORMAT(' READING BUFR DX-TABLES FROM ',I2,' TO ',I2)              
101   FORMAT(' COPYING BUFR DX-TABLES FROM ',I2,' TO ',I2)              
102   FORMAT(' READING USER DX-TABLES FROM ',I2,' TO ',I2)              
                                                                        
      RETURN                                                            
      END                                                               

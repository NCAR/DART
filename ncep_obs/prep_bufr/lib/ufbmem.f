C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBMEM(LUNIT,INEW,IRET,IUNIT)                          
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)                           
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 SEC0,CUNIT                                            
      DIMENSION   MBAY(5000)                                            
      EQUIVALENCE (SEC0,MBAY(1))                                        
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      WRITE(CUNIT,'(I4)') LUNIT                                         
                                                                        
C  TRY TO OPEN BUFR FILE AND SET TO INITIALIZE OR CONCATENATE           
C  ----------------------------------------------------------           
                                                                        
      CALL OPENBF(LUNIT,'IN',LUNIT)                                     
                                                                        
      IF(INEW.EQ.0) THEN                                                
         MSGP(0) = 0                                                    
         MUNIT = 0                                                      
         MLAST = 0                                                      
      ENDIF                                                             
                                                                        
      IMSG = 8/NBYTW+1                                                  
      NMSG = MSGP(0)                                                    
      IRET = 0                                                          
                                                                        
C  TRANSFER MESSAGES FROM FILE TO MEMORY - SET MESSAGE POINTERS         
C  ------------------------------------------------------------         
                                                                        
1     READ(LUNIT,ERR=900,END=100) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))      
      IRET = IRET+1                                                     
      NMSG = NMSG+1                                                     
      LMEM = LMSG(SEC0)+IMSG-1                                          
      IF(NMSG      .GT.MAXMSG) CALL BORT('UFBMEM - MSG BUFFER')        
      IF(LMEM+MLAST.GT.MAXMEM) CALL BORT('UFBMEM - MEM BUFFER')        
                                                                        
      DO I=1,LMEM                                                       
      MSGS(MLAST+I) = MBAY(I)                                           
      ENDDO                                                             
                                                                        
      MSGP(0000) = NMSG                                                 
      MSGP(NMSG) = MLAST+1                                              
      MLAST = MLAST+LMEM                                                
      GOTO 1                                                            
                                                                        
C  EXITS                                                                
C  -----                                                                
                                                                        
100   IF(IRET.EQ.0) THEN                                                
         CALL CLOSBF(LUNIT)                                             
      ELSE                                                              
         IF(MUNIT.NE.0) CALL CLOSBF(LUNIT)                              
         IF(MUNIT.EQ.0) MUNIT = LUNIT                                   
      ENDIF                                                             
      IUNIT = MUNIT                                                     
      RETURN                                                            
                                                                        
900   CALL BORT('UFBMEM - ERROR READING MESSAGE ON UNIT:'//CUNIT)      
      END                                                               

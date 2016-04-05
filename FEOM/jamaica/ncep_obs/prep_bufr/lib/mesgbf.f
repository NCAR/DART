C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MESGBF(LUNIT,MESGTYP)                                  
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*20 MESG                                                 
      CHARACTER*8  SEC0                                                 
      DIMENSION    MBAY(5000)                                           
      EQUIVALENCE  (MESG,MBAY(1))                                       
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      MESGTYP = -1                                                      
                                                                        
C  READ PAST ANY BUFR TABLES AND RETURN THE FIRST MESSAGE TYPE FOUND    
C  -----------------------------------------------------------------    
                                                                        
      CALL WRDLEN                                                       
      IMSG = 8/NBYTW+1                                                  
                                                                        
      REWIND LUNIT                                                      
1     READ(LUNIT,ERR=900,END=900) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))      
      MESGTYP = IUPM(MESG(17:17),8)                                      
      IF(MESGTYP.EQ.11) GOTO 1                                          
      REWIND LUNIT                                                      
                                                                        
900   RETURN                                                            
      END                                                               

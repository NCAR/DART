C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE COPYBF(LUNIN,LUNOT)                                    
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 SEC0                                                  
      CHARACTER*1 MOCT(24000)                                           
      DIMENSION   MBAY(5000)                                            
      EQUIVALENCE (MBAY(1),MOCT(1))                                     
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      ISEC0 = 8/NBYTW+1                                                 
      NMSG  = 0                                                         
                                                                        
C  CHECK BUFR FILE STATUSES                                             
C  ------------------------                                             
                                                                        
      CALL STATUS(LUNIN,LUN,IL,IM)                                      
      IF(IL.NE.0) GOTO 900                                              
      CALL STATUS(LUNOT,LUN,IL,IM)                                      
      IF(IL.NE.0) GOTO 901                                              
                                                                        
      REWIND(LUNIN)                                                     
      REWIND(LUNOT)                                                     
                                                                        
C  READ AND COPY A BUFR FILE ON UNIT LUNIN TO UNIT LUNOT                
C  -----------------------------------------------------                
                                                                        
1     READ(LUNIN,END=2,ERR=902) SEC0,(MBAY(I),I=ISEC0,LMSG(SEC0))       
      WRITE(LUNOT     ,ERR=903) SEC0,(MBAY(I),I=ISEC0,LMSG(SEC0))       
      GOTO 1                                                            
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
2     CLOSE(LUNIN)                                                      
      CLOSE(LUNOT)                                                      
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('COPYBF - INPUT  FILE IS CURRENTLY OPEN FOR BUFR')     
901   CALL BORT('COPYBF - OUTPUT FILE IS CURRENTLY OPEN FOR BUFR')     
902   CALL BORT('COPYBF - ERROR READING FILE    ')                     
903   CALL BORT('COPYBF - ERROR WRITING FILE    ')                     
      END                                                               

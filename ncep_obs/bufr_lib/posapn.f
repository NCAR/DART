C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE POSAPN(LUNIT)                                          
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 SEC0                                                  
      DIMENSION   MBAY(5000)                                            
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  READ AND COUNT MESSAGES IN ORDER TO POSTION FOR APPEND               
C  ------------------------------------------------------               
                                                                        
      IMSG = 8/NBYTW+1                                                  
      REWIND LUNIT                                                      
      IREC = 0                                                          
                                                                        
1     READ(LUNIT,ERR=2,END=2) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))          
      IREC = IREC+1                                                     
      GOTO 1                                                            
                                                                        
2     REWIND LUNIT                                                      
      DO J=1,IREC                                                       
      READ(LUNIT,ERR=900,END=901) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))      
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('POSAPN - IO ERR READING A MESSAGE')                   
901   CALL BORT('POSAPN - FAILURE TO READ TO EOFLE')                   
      END                                                               

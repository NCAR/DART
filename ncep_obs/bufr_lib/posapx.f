C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE POSAPX(LUNIT)                                          
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 SEC0                                                  
      DIMENSION   MBAY(5000)                                            
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  TRY TO READ TO THE END OF THE FILE AND BACKSPACE FOR APPEND          
C  ----------------------------------------------------                 
                                                                        
      REWIND LUNIT                                                      
      IMSG = 8/NBYTW+1                                                  
      IREC = 0                                                          
                                                                        
1     READ(LUNIT,END=2,ERR=3) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))          
      IREC = IREC+1                                                     
      GOTO 1                                                            
                                                                        
C  IF SUCCESSFUL BACKSPACE FOR APPENDING AND RETURN                     
C  ------------------------------------------------                     
                                                                        
2     BACKSPACE LUNIT                                                   
      RETURN                                                            
                                                                        
C  IF AN I/O ERROR IS ENCOUNTERED REREAD THE GOOD RECORDS AND RETURN    
C  -----------------------------------------------------------------    
                                                                        
3     REWIND LUNIT                                                      
      DO J=1,IREC                                                       
      READ(LUNIT,END=2,ERR=900) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))        
      ENDDO                                                             
      RETURN                                                            
                                                                        
C  IF ALL ELSE FAILS JUST GIVE UP AT THIS POINT                         
C  --------------------------------------------                         
                                                                        
900   CALL BORT('POSAPX - IO ERR REREADING A MESSAGE')                 
      END                                                               

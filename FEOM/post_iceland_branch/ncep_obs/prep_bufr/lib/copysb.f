C-----------------------------------------------------------------------
C  COPY A SUBSET IF LUNOT>0 OTHERWISE JUST ADVANCE THE INPUT SUBSET PTR 
C-----------------------------------------------------------------------
      SUBROUTINE COPYSB(LUNIN,LUNOT,IRET)                               
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND(STATUS,CPYUPD)                                             
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
                                                                        
C  CHECK THE FILE STATUSES                                              
C  -----------------------                                              
                                                                        
      CALL STATUS(LUNIN,LIN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      IF(IM.LE.0) GOTO 902                                              
                                                                        
      IF(LUNOT.GT.0) THEN                                               
         CALL STATUS(LUNOT,LOT,IL,IM)                                   
         IF(IL.EQ.0) GOTO 903                                           
         IF(IL.LT.0) GOTO 904                                           
         IF(IM.EQ.0) GOTO 905                                           
         IF(INODE(LIN).NE.INODE(LOT)) GOTO 906                          
      ENDIF                                                             
                                                                        
C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE                        
C  ---------------------------------------------                        
                                                                        
      IF(NSUB(LIN).EQ.MSUB(LIN)) THEN                                   
         IRET = -1                                                      
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  COPY THE SUBSET TO THE OUTPUT MESSAGE AND/OR RESET THE POINTERS      
C  ---------------------------------------------------------------      
                                                                        
      IBIT = (MBYT(LIN))*8                                              
      CALL UPB(NBYT,16,MBAY(1,LIN),IBIT)                                
      IF(LUNOT.GT.0) CALL CPYUPD(LUNOT,LIN,LOT,NBYT)                    
      MBYT(LIN) = MBYT(LIN) + NBYT                                      
      NSUB(LIN) = NSUB(LIN) + 1                                         
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('COPYSB - INPUT FILE IS CLOSED                    ')   
901   CALL BORT('COPYSB - INPUT FILE IS OPEN FOR OUTPUT           ')   
902   CALL BORT('COPYSB - NO INPUT FILE MESSAGE OPEN              ')   
903   CALL BORT('COPYSB - OUTPUT FILE IS CLOSED                   ')   
904   CALL BORT('COPYSB - OUTPUT FILE IS OPEN FOR INPUT           ')   
905   CALL BORT('COPYSB - NO OUTPUT FILE MESSAGE OPEN             ')   
906   CALL BORT('COPYSB - INPUT/OUTPUT FILES HAVE DIFFERENT TABLES')   
      END                                                               

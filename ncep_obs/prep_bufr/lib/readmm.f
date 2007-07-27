C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READMM(IMSG,SUBSET,JDATE,IRET)                         
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)                           
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
                                                                        
      CHARACTER*8 SUBSET                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE MESSAGE REQUEST AND FILE STATUS                            
C  -----------------------------------------                            
                                                                        
      CALL STATUS(MUNIT,LUN,IL,IM)                                      
      CALL WTSTAT(MUNIT,LUN,IL, 1)                                      
      IF(IL.GE.0) GOTO 900                                              
      IRET = 0                                                          
                                                                        
      IF(IMSG.EQ.0 .OR.IMSG.GT.MSGP(0)) THEN                            
         CALL WTSTAT(MUNIT,LUN,IL,0)                                    
         IRET = -1                                                      
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  READ MESSAGE# IMSG INTO A MESSAGE BUFFER                             
C  ----------------------------------------                             
                                                                        
      IPTR = MSGP(IMSG)                                                 
      IF(IMSG.LT.MSGP(0)) LPTR = MSGP(IMSG+1)-IPTR                      
      IF(IMSG.EQ.MSGP(0)) LPTR = MLAST-IPTR+1                           
      IPTR = IPTR-1                                                     
                                                                        
      DO I=1,LPTR                                                       
      MBAY(I,LUN) = MSGS(IPTR+I)                                        
      ENDDO                                                             
                                                                        
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,JRET)
      NMSG(LUN) = IMSG 
      IMSG = IMSG+1
      RETURN
 
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('READMM - FILE IS NOT OPEN FOR INPUT       ')          
901   CALL BORT('READMM - BAD IPTR IN MEMORY BUFFER        ')          
902   CALL BORT('READMM - BAD NPTR IN MEMORY BUFFER        ')          
903   CALL BORT('READMM - MSGTYPE MISMATCH FOR '//SUBSET    )          
      END                                                               

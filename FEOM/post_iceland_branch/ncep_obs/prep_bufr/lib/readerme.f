C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READERME(MESG,LUNIT,SUBSET,JDATE,IRET)                 
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
      CHARACTER*8 SUBSET,SEC0                                           
      DIMENSION   MESG(*),IEC0(2)                                       
      EQUIVALENCE (SEC0,IEC0)                                           
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      CALL WTSTAT(LUNIT,LUN,IL, 1)                                      
                                                                        
C  READ A MESSAGE INTO A MESSAGE BUFFER                                 
C  ------------------------------------                                 
                                                                        
      IEC0(1) = MESG(1)                                                 
      IEC0(2) = MESG(2)                                                 
      DO I=1,LMSG(SEC0)                                                 
      MBAY(I,LUN) = MESG(I)                                             
      ENDDO                                                             
      IF(SEC0(1:4).NE.'BUFR') GOTO 902                                  
                                                                        
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,IRET)
      RETURN
                                                                        
C  EOF ON ATTEMPTED READ                                                
C  ---------------------                                                
                                                                        
100   CALL WTSTAT(LUNIT,LUN,IL,0)                                      
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
      SUBSET = ' '                                                      
      JDATE = 0                                                         
      IRET = -1                                                         
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('READERME - FILE IS CLOSED                   ')        
901   CALL BORT('READERME - FILE IS OPEN FOR OUTPUT          ')        
902   CALL BORT('READERME - NOT A BUFR FILE                  ')        
      END                                                               

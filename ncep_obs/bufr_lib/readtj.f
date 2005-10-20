C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READTJ(LUNIT,SUBSET,JDATE,IRET)                        
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
      CHARACTER*8 SEC0,SUBSET                                           
      CHARACTER*4 BUFR                                                  
      DIMENSION   IEC0(2)
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
                                                                        
C  READ A MESSAGE INTO A MESSAGE BUFFER - SKIP DX MESSAGES              
C  -------------------------------------------------------              
                                                                        
1     MBIT = 0                                                          
      SEC0 = ' '                                                        
      IMSG = 8/NBYTW+1                                                  
      READ(LUNIT,ERR=902,END=100) SEC0,(MBAY(I,LUN),I=IMSG,LMSG(SEC0))  
      CALL CHRTRNA(BUFR,SEC0,4)                                         
      IF(BUFR.NE.'BUFR') GOTO 100                                       
      DO I=1,IMSG-1
      MBAY(I,LUN) = IEC0(I)
      ENDDO
                                                                        
C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------
 
      CALL CKTABA(LUN,SUBSET,JDATE,IRET)
      IF(IRET.NE.0) GOTO 1
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
                                                                        
900   CALL BORT('READTJ - FILE IS CLOSED                       ')      
901   CALL BORT('READTJ - FILE IS OPEN FOR OUTPUT              ')      
902   CALL BORT('READTJ - I/O ERROR READING MESSAGE            ')      
903   CALL BORT('READTJ - MSGTYPE MISMATCH  FOR '//SUBSET       )      
      END                                                               

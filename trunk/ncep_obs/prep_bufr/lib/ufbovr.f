C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBOVR(LUNIT,USR,I1,I2,IRET,STR)                       
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*(*) STR                                                 
      DIMENSION     USR(I1,I2)                                          
      REAL*8        USR,VAL                                             
                                                                        
C---------------------------------------------------------------------- 
CFPP$ EXPAND (STATUS,UFBRW)                                             
C---------------------------------------------------------------------- 
                                                                        
      IRET = 0                                                          
      IF(I2.LE.0) RETURN                                                
                                                                        
C  CHECK THE FILE STATUS AND I-NODE                                     
C  --------------------------------                                     
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.NE.1) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
                                                                        
      IO = MIN(MAX(0,IL),1)                                             
                                                                        
C  PARSE THE INPUT STRING - READ/WRITE VALUES                           
C  ------------------------------------------                           
                                                                        
      CALL STRING(STR,LUN,I1,IO)                                        
      CALL TRYBUMP(LUNIT,LUN,USR,I1,I2,IO,IRET)                         
                                                                        
      IF(IO.EQ.1 .AND. IRET.NE.I2) THEN                                 
         IF(IRET.NE.I2) PRINT*,STR                                      
         IF(IRET.NE.I2) GOTO 903                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('UFBOVR - FILE IS NOT OPEN FOR OUTPUT        ')        
901   CALL BORT('UFBOVR - NO MESSAGE OPEN                    ')        
902   CALL BORT('UFBOVR - I-NODE MISMATCH                    ')        
903   CALL BORT('UFBOVR - INCOMPLETE WRITE                   ')        
      END                                                               

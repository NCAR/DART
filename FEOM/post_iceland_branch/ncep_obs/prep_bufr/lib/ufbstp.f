C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBSTP(LUNIO,USR,I1,I2,IRET,STR)                       
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*(*) STR                                                 
      DIMENSION     USR(I1,I2)                                          
      REAL*8        USR,VAL,BMISS                                 

      DATA BMISS /10E10/
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE FILE STATUS AND I-NODE                                     
C  --------------------------------                                     
                                                                        
      LUNIT = ABS(LUNIO)                                                
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
      IO = MIN(MAX(0,IL),1)                                             
      IF(LUNIO.NE.LUNIT) IO = 0                                         

C  INITIALIZE USR ARRAY PRECEEDING AN INPUT OPERATION                   
C  --------------------------------------------------                   
                                                                        
      IF(IO.EQ.0) THEN                                                  
         DO I=1,I1*I2                                                   
         USR(I,1) = BMISS                                               
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  PARSE THE INPUT STRING - READ/WRITE VALUES                           
C  ------------------------------------------                           
                                                                        
      CALL STRING(STR,LUN,I1,IO)                                        
      CALL UFBSP(LUN,USR,I1,I2,IO,IRET)                                 
                                                                        
      IF(IO.EQ.1 .AND. IRET.NE.I2) THEN                                 
         PRINT*,STR,' i2=',i2,' iret=',iret                             
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('UFBSTP - FILE IS CLOSED                     ')        
901   CALL BORT('UFBSTP - NO MESSAGE OPEN                    ')        
902   CALL BORT('UFBSTP - I-NODE MISMATCH                    ')        
903   CALL BORT('UFBSTP - INCOMPLETE WRITE                   ')        
      END                                                               

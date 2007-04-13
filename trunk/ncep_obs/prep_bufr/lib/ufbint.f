C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBINT(LUNIN,USR,I1,I2,IRET,STR)                       
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*(*) STR                                                 
      DIMENSION     USR(I1,I2)                                          
      REAL*8        USR,VAL,BMISS
                                                                        
      DATA BMISS /10E10/                                                
                                                                        
C---------------------------------------------------------------------- 
CFPP$ EXPAND (STATUS,UFBRW)                                             
C---------------------------------------------------------------------- 
                                                                        
      IRET = 0                                                          
      IF(I1.LE.0) RETURN                                                
      IF(I2.LE.0) RETURN                                                
                                                                        
C  CHECK THE FILE STATUS AND I-NODE AND PARSE OR RECALL THE STRING      
C  ---------------------------------------------------------------      
                                                                        
      LUNIT = ABS(LUNIN)                                                
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
                                                                        
      IO = MIN(MAX(0,IL),1)                                             
      IF(LUNIT.NE.LUNIN) IO = 0                                         
                                                                        
      CALL STRING(STR,LUN,I1,IO)                                        
                                                                        
C  INITIALIZE USR ARRAY PRECEEDING AN INPUT OPERATION                   
C  --------------------------------------------------                   
                                                                        
      IF(IO.EQ.0) THEN                                                  
         DO I=1,I1*I2                                                   
         USR(I,1) = BMISS                                               
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  CALL THE MNEMONIC READER/WRITER                                      
C  -------------------------------                                      
                                                                        
      CALL UFBRW(LUN,USR,I1,I2,IO,IRET)                                 
                                                                        
C  IF INCOMPLETE WRITE TRY TO INITIALIZE REPLICATION SEQUENCE OR RETURN 
C  ---------------------------------------------------------------------
                                                                        
      IF(IO.EQ.1 .AND. IRET.NE.I2 .AND. IRET.GE.0) THEN                 
         CALL TRYBUMP(LUNIT,LUN,USR,I1,I2,IO,IRET)                      
         IF(IRET.NE.I2) PRINT*,STR                                      
         IF(IRET.NE.I2) GOTO 903                                        
      ELSEIF(IRET.EQ.-1) THEN                                           
         IRET = 0                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('UFBINT - FILE IS CLOSED                     ')        
901   CALL BORT('UFBINT - NO MESSAGE OPEN                    ')        
902   CALL BORT('UFBINT - I-NODE MISMATCH                    ')        
903   CALL BORT('UFBINT - INCOMPLETE WRITE                   ')        
      END                                                               

C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBEVN(LUNIT,USR,I1,I2,I3,IRET,STR)                    
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*10  TAG                                                 
      CHARACTER*3   TYP                                                 
      DIMENSION     USR(I1,I2,I3),INVN(255)                             
      REAL*8        VAL,USR,BMISS                                       
                                                                        
      DATA BMISS /10E10/                                                
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE FILE STATUS AND I-NODE                                     
C  --------------------------------                                     
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
                                                                        
C  PARSE THE INPUT STRING                                               
C  ----------------------                                               
                                                                        
      CALL STRING(STR,LUN,I1,0)                                         
                                                                        
C  SET INITIAL VALUES FOR RETURNING ARGUMENTS                           
C  ------------------------------------------                           
                                                                        
      DO I=1,I1*I2*I3                                                   
      USR(I,1,1) = BMISS                                                
      ENDDO                                                             
                                                                        
      IRET = 0                                                          
                                                                        
C  LOOP OVER COND WINDOWS                                               
C  ----------------------                                               
                                                                        
      INC1 = 1                                                          
      INC2 = 1                                                          
                                                                        
1     CALL CONWIN(LUN,INC1,INC2,I2)                                     
      IF(NNOD.EQ.0) THEN                                                
         IRET = I2                                                      
         RETURN                                                         
      ELSEIF(INC1.EQ.0) THEN                                            
         RETURN                                                         
      ELSE                                                              
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            INS2 = INC1                                                 
            CALL GETWIN(NODS(I),LUN,INS1,INS2)                          
            IF(INS1.EQ.0) RETURN                                        
            GOTO 2                                                      
         ENDIF                                                          
         ENDDO                                                          
         INS1 = INC1                                                    
         INS2 = INC2                                                    
      ENDIF                                                             
                                                                        
C  READ PUSH DOWN STACK DATA INTO 3D ARRAYS                             
C  ----------------------------------------                             
                                                                        
2     IRET = IRET+1                                                     
      IF(IRET.LE.I2) THEN                                               
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            NNVN = NVNWIN(NODS(I),LUN,INS1,INS2,INVN,I3)                
            DO N=1,NNVN                                                 
            USR(I,IRET,N) = VAL(INVN(N),LUN)                            
            ENDDO                                                       
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  DECIDE WHAT TO DO NEXT                                               
C  ----------------------                                               
                                                                        
      CALL NXTWIN(LUN,INS1,INS2)                                        
      IF(INS1.GT.0 .AND. INS1.LT.INC2) GOTO 2                           
      IF(NCON.GT.0) GOTO 1                                              
                                                                        
      RETURN                                                            
900   CALL BORT('UFBEVN - FILE IS CLOSED                     ')        
901   CALL BORT('UFBEVN - NO MESSAGE OPEN                    ')        
902   CALL BORT('UFBEVN - I-NODE MISMATCH                    ')        
      END                                                               

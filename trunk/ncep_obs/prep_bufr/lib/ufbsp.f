C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBSP(LUN,USR,I1,I2,IO,IRET)                           
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      DIMENSION    USR(I1,I2)
      REAL*8       USR,VAL,BMISS                                 

      DATA BMISS /10E10/
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IRET = 0                                                          
      INS1 = 0                                                          
      INS2 = 0                                                          
                                                                        
C  FRAME A SECTION OF THE BUFFER - RETURN WHEN NO FRAME                 
C  ----------------------------------------------------                 
                                                                        
1     IF(INS1+1.GT.NVAL(LUN)) RETURN                                    
      INS1 = INVtag(NODS(1),LUN,INS1+1,NVAL(LUN))                       
      IF(INS1.EQ.0) RETURN                                              
                                                                        
      INS2 = INVtag(NODS(1),LUN,INS1+1,NVAL(LUN))                       
      IF(INS2.EQ.0) INS2 = NVAL(LUN)                                    
      IRET = IRET+1                                                     
                                                                        
C  READ USER VALUES                                                     
C  ----------------                                                     
                                                                        
      IF(IO.EQ.0 .AND. IRET.LE.I2) THEN                                 
         INVM = INS1
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            INVN = INVTAG(NODS(I),LUN,INVM,INS2)                        
            IF(INVN.GT.0) USR(I,IRET) = VAL(INVN,LUN)                   
            INVM = MAX(INVN,INVM)
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  WRITE USER VALUES                                                    
C  -----------------                                                    
                                                                        
      IF(IO.EQ.1 .AND. IRET.LE.I2) THEN                                 
         INVM = INS1
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            INVN = INVtag(NODS(I),LUN,INVM,INS2)                        
            IF(INVN.GT.0) VAL(INVN,LUN) = USR(I,IRET)                   
            INVM = MAX(INVN,INVM)
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  GO FOR NEXT FRAME                                                    
C  -----------------                                                    
                                                                        
      GOTO 1                                                            
                                                                        
      END                                                               

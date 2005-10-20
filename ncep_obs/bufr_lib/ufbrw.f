C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBRW(LUN,USR,I1,I2,IO,IRET)                           
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      REAL*8       USR(I1,I2),VAL,BMISS                                 

      DATA BMISS/10E10/
                                                                        
C---------------------------------------------------------------------- 
CFPP$ EXPAND (CONWIN,DRSTPL,GETWIN,INVWIN,LSTRPS,NEWWIN,NXTWIN)          
C---------------------------------------------------------------------- 
                                                                        
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
         IRET = -1                                                      
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  LOOP OVER STORE NODES                                                
C  ---------------------                                                
                                                                        
2     IRET = IRET+1                                                     
                                                                        
C     print*,'ufbrw:',iret,':',ins1,':',ins2,':',inc1,':',inc2          
C     print'(5a10)',(tag(inv(i,lun)),i=ins1,ins2)                       
                                                                        
C  WRITE USER VALUES                                                    
C  -----------------                                                    
                                                                        
      IF(IO.EQ.1 .AND. IRET.LE.I2) THEN                                 
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            IF(USR(I,IRET).NE.BMISS) THEN                               
               INVN = INVWIN(NODS(I),LUN,INS1,INS2)                     
               IF(INVN.EQ.0) THEN                                       
                  CALL DRSTPL(NODS(I),LUN,INS1,INS2,INVN)               
                  if(invn.eq.0) then                                    
                     iret = 0                                           
                     return                                             
                  endif                                                 
                  CALL NEWWIN(LUN,INC1,INC2)                            
                  VAL(INVN,LUN) = USR(I,IRET)                           
               ELSEIF(LSTRPS(NODS(I),LUN).EQ.0) THEN                    
                  VAL(INVN,LUN) = USR(I,IRET)                           
               ELSEIF(VAL(INVN,LUN).EQ.BMISS) THEN                      
                  VAL(INVN,LUN) = USR(I,IRET)                           
               ELSE                                                     
                  CALL DRSTPL(NODS(I),LUN,INS1,INS2,INVN)               
                  if(invn.eq.0) then                                    
                     iret = 0                                           
                     return                                             
                  endif                                                 
                  CALL NEWWIN(LUN,INC1,INC2)                            
                  VAL(INVN,LUN) = USR(I,IRET)                           
               ENDIF                                                    
            ENDIF                                                       
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  READ USER VALUES                                                     
C  ----------------                                                     
                                                                        
      IF(IO.EQ.0 .AND. IRET.LE.I2) THEN                                 
         DO I=1,NNOD                                                    
         USR(I,IRET) = BMISS                                            
         IF(NODS(I).GT.0) THEN                                          
            INVN = INVWIN(NODS(I),LUN,INS1,INS2)                        
            IF(INVN.GT.0) USR(I,IRET) = VAL(INVN,LUN)                   
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  DECIDE WHAT TO DO NEXT                                               
C  ----------------------                                               
                                                                        
      IF(IO.EQ.1.AND.IRET.EQ.I2) RETURN                                 
      CALL NXTWIN(LUN,INS1,INS2)                                        
      IF(INS1.GT.0 .AND. INS1.LT.INC2) GOTO 2                           
      IF(NCON.GT.0) GOTO 1                                              
                                                                        
      RETURN                                                            
      END                                                               

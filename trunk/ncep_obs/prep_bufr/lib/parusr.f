C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE PARUSR(STR,LUN,I1,IO)                                  
                                                                        
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
      common /ACMODE/ iac
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*80  UST                                                 
      CHARACTER*20  UTG(30)                                             
      LOGICAL       BUMP                                                
                                                                        
      DATA MAXUSR /30/                                                  
      DATA MAXNOD /20/                                                  
      DATA MAXCON /10/                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      UST  = STR                                                        
      IF(LEN(STR).GT.80) GOTO 900                                       
                                                                        
      NCON = 0                                                          
      NNOD = 0                                                          
                                                                        
C  PROCESS STRING PIECES(S) INTO COND NODES AND STORE NODES             
C  --------------------------------------------------------             
                                                                        
      CALL PARSEQ(UST,UTG,MAXUSR,NTOT)                                  
                                                                        
      DO N=1,NTOT                                                       
      CALL PARUTG(LUN,IO,UTG(N),NOD,KON,VAL,*908)                       
      IF(KON.NE.0) THEN                                                 
         NCON = NCON+1                                                  
         IF(NCON.GT.MAXCON) GOTO 901                                    
         NODC(NCON) = NOD                                               
         KONS(NCON) = KON                                               
         IVLS(NCON) = NINT(VAL)
      ELSE                                                              
         NNOD = NNOD+1                                                  
         IF(NNOD.GT.MAXNOD) GOTO 902                                    
         NODS(NNOD) = NOD                                               
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  SORT COND NODES IN JUMP/LINK TABLE ORDER                             
C  ----------------------------------------                             
                                                                        
      DO I=1,NCON                                                       
      DO J=I+1,NCON                                                     
      IF(NODC(I).GT.NODC(J)) THEN                                       
         NOD     = NODC(I)                                              
         NODC(I) = NODC(J)                                              
         NODC(J) = NOD                                                  
                                                                        
         KON     = KONS(I)                                              
         KONS(I) = KONS(J)                                              
         KONS(J) = KON                                                  
                                                                        
         VAL     = IVLS(I)                                              
         IVLS(I) = IVLS(J)                                              
         IVLS(J) = VAL                                                  
      ENDIF                                                             
      ENDDO                                                             
      ENDDO                                                             
                                                                        
C  CHECK ON SPECIAL RULES FOR BUMP NODES                                
C  -------------------------------------                                
                                                                        
      BUMP = .FALSE.                                                    
                                                                        
      DO N=1,NCON                                                       
      IF(KONS(N).EQ.5) THEN                                             
         IF(IO.EQ.0)   GOTO 903                                         
         IF(N.NE.NCON) GOTO 904                                         
         BUMP = .TRUE.                                                  
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  CHECK STORE NODE COUNT AND ALIGNMENT                                 
C  ------------------------------------                                 
                                                                        
      IF(.NOT.BUMP .AND. NNOD.EQ.0) GOTO 905                            
      IF(NNOD.GT.I1)                GOTO 906                            
                                                                        
      IRPC = -1                                                         
      DO I=1,NNOD                                                       
      IF(NODS(I).GT.0) THEN                                             
         IF(IRPC.LT.0) IRPC = LSTRPC(NODS(I),LUN)                       
         IF(IRPC.NE.LSTRPC(NODS(I),LUN).and.iac.eq.0) GOTO 907  
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('PARUSR - USER STRING > 80 CHARS         :'//UST)      
901   CALL BORT('PARUSR - TOO MANY COND NODES            :'//UST)      
902   CALL BORT('PARUSR - TOO MANY STOR NODES            :'//UST)      
903   CALL BORT('PARUSR - BUMP ON INPUT NOT ALLOWED      :'//UST)      
904   CALL BORT('PARUSR - BUMP MUST BE ON INNER NODE     :'//UST)      
905   CALL BORT('PARUSR - USER STRING HAS NO STORE NODES :'//UST)      
906   CALL BORT('PARUSR - MUST BE AT LEAST I1 STORE NODES:'//UST)      
907   CALL BORT('PARUSR - STORE NODES MUST IN ONE REP GRP:'//UST)      
908   CALL BORT('PARUSR - PARUTG:'                         //UST)      
      END                                                               

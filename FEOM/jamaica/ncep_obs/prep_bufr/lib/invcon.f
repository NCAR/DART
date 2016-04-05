C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      FUNCTION INVCON(NC,LUN,INV1,INV2)                                 
                                                                        
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),IVLS(10),KONS(10)     
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE INVENTORY INTERVAL                                         
C  ----------------------------                                         
                                                                        
      IF(INV1.LE.0 .OR. INV1.GT.NVAL(LUN)) GOTO 99                      
      IF(INV2.LE.0 .OR. INV2.GT.NVAL(LUN)) GOTO 99                      
                                                                        
C  FIND AN OCCURANCE OF NODE IN THE WINDOW MEETING THIS CONDITION       
C  --------------------------------------------------------------       
                                                                        
      DO INVCON=INV1,INV2                                               
      IF(INV(INVCON,LUN).EQ.NODC(NC)) THEN                              
         IF(KONS(NC).EQ.1 .AND. VAL(INVCON,LUN).EQ.IVLS(NC)) RETURN     
         IF(KONS(NC).EQ.2 .AND. VAL(INVCON,LUN).NE.IVLS(NC)) RETURN     
         IF(KONS(NC).EQ.3 .AND. VAL(INVCON,LUN).LT.IVLS(NC)) RETURN     
         IF(KONS(NC).EQ.4 .AND. VAL(INVCON,LUN).GT.IVLS(NC)) RETURN     
      ENDIF                                                             
      ENDDO                                                             
                                                                        
99    INVCON = 0                                                        
      RETURN                                                            
      END                                                               

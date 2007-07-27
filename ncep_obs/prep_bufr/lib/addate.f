C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ADDATE(IDATE,JH,JDATE)                                 
                                                                        
      DIMENSION   MON(12)                                               
                                                                        
      DATA MON/31,28,31,30,31,30,31,31,30,31,30,31/                     
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IY = IDATE/1000000
      IM = MOD(IDATE/10000  ,100)                                       
      ID = MOD(IDATE/100    ,100)                                       
      IH = MOD(IDATE        ,100)                                       
      IH = IH+JH                                                        
                                                                        
      IF(MOD(IY,4)  .NE.0) MON(2) = 28                                    
      IF(MOD(IY,4)  .EQ.0) MON(2) = 29                                    
      IF(MOD(IY,100).EQ.0) MON(2) = 28                                    
      IF(MOD(IY,400).EQ.0) MON(2) = 29                                    
                                                                        
1     IF(IH.LT.0) THEN                                                  
         IH = IH+24                                                     
         ID = ID-1                                                      
         IF(ID.EQ.0) THEN                                               
            IM = IM-1                                                   
            IF(IM.EQ.0) THEN                                            
               IM = 12                                                  
               IY = IY-1                                                
               IF(IY.LT.0) IY = 99                                      
            ENDIF                                                       
            ID = MON(IM)                                                
         ENDIF                                                          
         GOTO 1                                                         
      ELSEIF(IH.GE.24) THEN                                             
         IH = IH-24                                                     
         ID = ID+1                                                      
         IF(ID.GT.MON(IM)) THEN                                         
            ID = 1                                                      
            IM = IM+1                                                   
            IF(IM.GT.12) THEN                                           
               IM = 1                                                   
               IY = IY+1
               IF(IY.EQ.100) IY = 00
            ENDIF                                                       
         ENDIF                                                          
         GOTO 1                                                         
      ENDIF                                                             
                                                                        
      JDATE = IY*1000000 + IM*10000 + ID*100 + IH                       
                                                                        
      RETURN                                                            
      END                                                               

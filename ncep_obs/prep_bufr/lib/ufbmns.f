C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBMNS(IREP,SUBSET,IDATE)                              
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)                           
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
                                                                        
      CHARACTER*8   SUBSET                                              
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
      JREP = 0                                                          
      IMSG = 1                                                          
                                                                        
C  READ SUBSET #ISUB FROM MEMORY MESSAGE #IMSG                          
C  -------------------------------------------                          
                                                                        
      DO WHILE(IRET.EQ.0)                                               
      CALL RDMEMM(IMSG,SUBSET,IDATE,IRET)                               
      IF(IRET.NE.0) GOTO 900                                            
      IF(JREP+NMSUB(MUNIT).GE.IREP) THEN                                
         CALL RDMEMS(IREP-JREP,IRET)                                    
         IF(IRET.NE.0) GOTO 901                                         
         RETURN                                                         
      ELSE                                                              
         JREP = JREP+NMSUB(MUNIT)                                       
         IMSG = IMSG+1                                                  
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      CALL BORT('UFBMNS - REPORT NUMBER OUT OF RANGE')                 
900   CALL BORT('UFBMNS - BAD RETURN FROM RDMEMM')                     
901   CALL BORT('UFBMNS - BAD RETURN FROM RDMEMS')                     
      END                                                               

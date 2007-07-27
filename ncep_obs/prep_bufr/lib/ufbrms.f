C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBRMS(IMSG,ISUB,USR,I1,I2,IRET,STR)                   
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)                           
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*8   SUBSET                                              
      DIMENSION     USR(I1,I2)                                          
      REAL*8        USR                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  UFBINT SUBSET #ISUB FROM MEMORY MESSAGE #IMSG                        
C  ---------------------------------------------                        
                                                                        
      CALL RDMEMM(IMSG,SUBSET,IDATE,IRET)                               
      IF(IRET.NE.0) GOTO 900                                            
      CALL RDMEMS(ISUB,IRET)                                            
      IF(IRET.NE.0) GOTO 901                                            
                                                                        
      CALL UFBINT(MUNIT,USR,I1,I2,IRET,STR)                             
                                                                        
      RETURN                                                            
900   CALL BORT('UFBRMS - BAD RETURN FROM RDMEMM')                     
901   CALL BORT('UFBRMS - BAD RETURN FROM RDMEMS')                     
      END                                                               

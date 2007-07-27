C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBMMS(IMSG,ISUB,SUBSET,IDATE)                         
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)                           
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
                                                                        
      CHARACTER*8   SUBSET                                              
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  UFBINT SUBSET #ISUB FROM MEMORY MESSAGE #IMSG                        
C  ---------------------------------------------                        
                                                                        
      CALL RDMEMM(IMSG,SUBSET,IDATE,IRET)                               
      IF(IRET.NE.0) GOTO 900                                            
      CALL RDMEMS(ISUB,IRET)                                            
      IF(IRET.NE.0) GOTO 901                                            
                                                                        
      RETURN                                                            
900   CALL BORT('UFBMMS - BAD RETURN FROM RDMEMM')                     
901   CALL BORT('UFBMMS - BAD RETURN FROM RDMEMS')                     
      END                                                               

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CLOSMG(LUNIT)                                          
                                                                        
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.NE.0) CALL MSGWRT(LUNIT,MBAY(1,LUN),MBYT(LUN))              
      CALL WTSTAT(LUNIT,LUN,IL,0)                                       
                                                                        
      RETURN                                                            
900   CALL BORT('CLOSMG - FILE IS CLOSED            ')                 
901   CALL BORT('CLOSMG - FILE IS OPEN FOR INPUT    ')                 
      END                                                               

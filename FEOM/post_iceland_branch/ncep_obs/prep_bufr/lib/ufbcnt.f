C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBCNT(LUNIT,KMSG,KSUB)                                
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUS - RETURN THE MESSAGE AND SUBSET COUNTERS       
C  --------------------------------------------------------------       
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) CALL BORT('UFBCNT - FILE IS CLOSED')                 
      KMSG = NMSG(LUN)                                                  
      KSUB = NSUB(LUN)                                                  
      RETURN                                                            
      END                                                               

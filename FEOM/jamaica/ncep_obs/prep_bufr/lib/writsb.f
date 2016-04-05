      SUBROUTINE WRITSB(LUNIT)                                          

C************************************************************************
C* WRITSB								*
C*									*
C* This subroutine should only be called when logical unit LUNIT has	*
C* been opened for output operations.  It packs up the current subset	*
C* within memory and then tries to add it to the BUFR message that is	*
C* currently open within memory for this LUNIT.  If the subset will not	*
C* fit into the currently open message (as determined via a call to	*
C* subroutine MSGUPD), then that message is flushed to LUNIT and a new	*
C* one is created in order to hold the current subset.			*
C*									*
C* WRITSB  ( LUNIT )							*
C*									*
C* Input parameters:							*
C*	LUNIT		INTEGER		FORTRAN logical unit number	*
C*									*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************

C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.EQ.0) GOTO 902                                              
                                                                        
C  PACK UP THE SUBSET AND PUT IT INTO THE MESSAGE                       
C  ----------------------------------------------                       
                                                                        
      CALL WRTREE(LUN)                                                  
      CALL MSGUPD(LUNIT,LUN)                                            
                                                                        
      RETURN                                                            
900   CALL BORT('WRITSB - FILE IS CLOSED                     ')        
901   CALL BORT('WRITSB - FILE IS OPEN FOR INPUT             ')        
902   CALL BORT('WRITSB - NO MESSAGE OPEN                    ')        
      END                                                               

      SUBROUTINE OPENMG(LUNIT,SUBSET,JDATE)                             

C************************************************************************
C* OPENMG								*
C*									*
C* This subroutine should only be called when logical unit LUNIT has	*
C* been opened for output operations.  Like subroutine OPENMB, it opens	*
C* and initializes a new BUFR message within memory; however, unlike	*
C* subroutine OPENMB, it does so regardless of the values of SUBSET and	*
C* JDATE.  If there is already a BUFR message open within memory for	*
C* this	LUNIT, then that message will be closed and flushed to LUNIT	*
C* before opening the new one.						*
C*									*
C* OPENMG  ( LUNIT, SUBSET, JDATE )					*
C*									*
C* Input parameters:							*
C*	LUNIT		INTEGER		FORTRAN logical unit number	*
C*	SUBSET		CHARACTER*8	Table A mnemonic for type of	*
C*					BUFR message to be opened	*
C*	JDATE		INTEGER		Date-time to be stored within	*
C*					Section 1 of BUFR message,	*
C*					in format of either YYMMDDHH	*
C*					or YYYYMMDDHH			*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
                                                                        
      CHARACTER*(*) SUBSET                                              
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.NE.0) CALL CLOSMG(LUNIT)                                    
      CALL WTSTAT(LUNIT,LUN,IL, 1)                                      
                                                                        
C  GET SOME SUBSET PARTICULARS                                          
C  ---------------------------                                          
                                                                        
      CALL NEMTBA(LUN,SUBSET,MTYP,MSTB,INOD)                            
      INODE(LUN) = INOD                                                 
      IDATE(LUN) = I4DY(JDATE)
                                                                        
C  INITIALIZE THE OPEN MESSAGE                                          
C  ---------------------------                                          
                                                                        
      CALL MSGINI(LUN)                                                  
      CALL USRTPL(LUN,1,1)                                              
                                                                        
      RETURN                                                            
900   CALL BORT('OPENMG - FILE IS CLOSED            ')                 
901   CALL BORT('OPENMG - FILE IS OPEN FOR INPUT    ')                 
      END                                                               

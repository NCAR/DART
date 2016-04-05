      SUBROUTINE OPENMB(LUNIT,SUBSET,JDATE)                             

C************************************************************************
C* OPENMB								*
C*									*
C* This subroutine should only be called when logical unit LUNIT has	*
C* been opened for output operations.  It opens and initializes a new	*
C* BUFR message within memory, unless there is already a BUFR message	*
C* open within memory for this LUNIT which has the same SUBSET and	*
C* JDATE values, in which case it does nothing.  Otherwise, if there is	*
C* already a BUFR message open within memory for this LUNIT but which	*
C* has a different SUBSET or JDATE value, then that message will be	*
C* closed and flushed to LUNIT before opening the new one.		*
C*									*
C* OPENMB  ( LUNIT, SUBSET, JDATE )					*
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
      LOGICAL       OPEN                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
                                                                        
C  GET SOME SUBSET PARTICULARS                                          
C  ---------------------------                                          

      CALL NEMTBA(LUN,SUBSET,MTYP,MSTB,INOD)
      OPEN = IM.EQ.0.OR.INOD.NE.INODE(LUN).OR.I4DY(JDATE).NE.IDATE(LUN)
                                                                        
C  MAYBE OPEN A NEW OR DIFFERENT TYPE OF MESSAGE                        
C  ---------------------------------------------                        
                                                                        
      IF(OPEN) THEN                                                     
         CALL CLOSMG(LUNIT)                                             
         CALL WTSTAT(LUNIT,LUN,IL, 1)                                   
         INODE(LUN) = INOD                                              
         IDATE(LUN) = I4DY(JDATE)
         CALL MSGINI(LUN)                                               
         CALL USRTPL(LUN,1,1)                                           
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('OPENMB - FILE IS CLOSED            ')                 
901   CALL BORT('OPENMB - FILE IS OPEN FOR INPUT    ')                 
      END                                                               

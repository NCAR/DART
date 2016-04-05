      SUBROUTINE WRITSA(LUNXX,MSGT,MSGL)                                

C************************************************************************
C* WRITSA								*
C*									*
C* This subroutine should only be called when logical unit ABS(LUNXX)	*
C* has been opened for output operations.  Like subroutine WRITSB, it	*
C* packs up the current subset within memory and then tries to add it	*
C* to the BUFR message that is currently open within memory for		*
C* ABS(LUNXX), automatically closing and flushing (to ABS(LUNXX)) this	*
C* current message and opening a new one when necessary in order to	*
C* hold the current subset.  However, unlike subroutine WRITSB, this	*
C* subroutine also returns a copy of each completed BUFR message to the	*
C* calling routine (in addition to flushing it to ABS(LUNXX)!).		*
C* Further, this subroutine can even be called in an alternate context	*
C* with no current subset in memory but where LUNXX contains a negative	*
C* value, in which case this is a signal to force the return (and	*
C* corresponding flush to ABS(LUNXX)!) of the current message within	*
C* memory (if there is one!).						*
C*									*
C* WRITSA  ( LUNXX, MSGT, MSGL )					*
C*									*
C* Input parameters:							*
C*	LUNXX		INTEGER		FORTRAN logical unit number	*
C*					and "force" flag:		*
C*					* ABS(LUNXX) always denotes the	*
C*					  FORTRAN logical unit number of*
C*					  the output file.		*
C*					* If LUNXX < 0, then no subset	*
C*					  is in memory, but any current	*
C*					  message in memory will be	*
C*					  forcibly flushed to ABS(LUNXX)*
C*					  *and* to array MSGT.		*
C*									*
C* Output parameters:							*
C*	MSGT(MSGL)	INTEGER		BUFR message			*
C*	MSGL		INTEGER		Length of MSGT:			*
C*					  0 = No message was returned	*
C*									*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /BUFRMG/ MSGLEN,MSGTXT(5000)                               
                                                                        
      DIMENSION MSGT(*)                                                 
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      LUNIT = ABS(LUNXX)                                                
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.EQ.0) GOTO 902                                              
                                                                        
C  SEE IF A MEMORY MESSAGE IS WAITING OR FORCED                         
C  --------------------------------------------                         
                                                                        
      IF(LUNXX.LT.0) CALL CLOSMG(LUNIT)                                 
                                                                        
      IF(MSGLEN.GT.0) THEN                                              
         MSGL = MSGLEN                                                  
         DO N=1,MSGL                                                    
         MSGT(N) = MSGTXT(N)                                            
         ENDDO                                                          
         MSGLEN = 0                                                     
      ELSE                                                              
         MSGL = 0                                                       
      ENDIF                                                             
                                                                        
      IF(LUNXX.LT.0) RETURN                                             
                                                                        
C  PACK UP THE SUBSET AND PUT IT INTO THE MESSAGE                       
C  ----------------------------------------------                       
                                                                        
      CALL WRTREE(LUN)                                                  
      CALL MSGUPD(LUNIT,LUN)                                            
                                                                        
      RETURN                                                            
900   CALL BORT('WRITSA - FILE IS CLOSED                     ')        
901   CALL BORT('WRITSA - FILE IS OPEN FOR INPUT             ')        
902   CALL BORT('WRITSA - NO MESSAGE OPEN                    ')        
      END                                                               

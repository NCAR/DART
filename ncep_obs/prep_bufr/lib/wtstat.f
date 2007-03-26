      SUBROUTINE WTSTAT(LUNIT,LUN,IL,IM)                                

C************************************************************************
C* WTSTAT								*
C*									*
C* Depending on the value of IL, this subroutine either stores or	*
C* deletes the values of LUNIT, IL, and IM within the internal arrays	*
C* IOLUN(LUN) and IOMSG(LUN).						*
C*									*
C* WTSTAT  ( LUNIT, LUN, IL, IM )					*
C*									*
C* Input parameters:							*
C*	LUNIT		INTEGER		FORTRAN logical unit number	*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					arrays for logical unit LUNIT	*
C*	IL		INTEGER		Logical unit status indicator:	*
C*					  0 = delete information about	*
C*					      LUNIT from within		*
C*					      internal arrays		*
C*					  1 = LUNIT is an output file	*
C*					 -1 = LUNIT is an input file	*
C*	IM		INTEGER		Indicator as to whether there is*
C*					a BUFR message currently open	*
C*					within memory for this LUNIT:	*
C*					  0 = no			*
C*					  1 = yes			*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /STBFR/ IOLUN(32),IOMSG(32)                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK ON THE ARGUMENTS                                               
C  ----------------------                                               
                                                                        
      IF(LUNIT.LE.0)            GOTO 900                                
      IF(LUN  .LE.0)            GOTO 901                                
      IF(IL.LT.-1 .OR. IL.GT.1) GOTO 902                                
      IF(IM.LT. 0 .OR. IM.GT.1) GOTO 903                                
                                                                        
C  CHECK ON LUNIT-LUN COMBINATION                                       
C  ------------------------------                                       
                                                                        
      IF(ABS(IOLUN(LUN)).NE.LUNIT) THEN                                 
         IF(IOLUN(LUN).NE.0) GOTO 905                                   
      ENDIF                                                             
                                                                        
C  RESET THE FILE STATUSES                                              
C  -----------------------                                              
                                                                        
      IF(IL.NE.0) THEN                                                  
         IOLUN(LUN) = SIGN(LUNIT,IL)                                    
         IOMSG(LUN) = IM                                                
      ELSE                                                              
         IOLUN(LUN) = 0                                                 
         IOMSG(LUN) = 0                                                 
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('WTSTAT - BAD LUNIT                               ')   
901   CALL BORT('WTSTAT - BAD LUN                                 ')   
902   CALL BORT('WTSTAT - BAD IL                                  ')   
903   CALL BORT('WTSTAT - BAD IM                                  ')   
905   CALL BORT('WTSTAT - ATTEMPT TO REDEFINE EXISITING FILE UNIT ')   
      END                                                               

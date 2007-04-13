      SUBROUTINE STATUS(LUNIT,LUN,IL,IM)                                

C************************************************************************
C* STATUS								*
C*									*
C* This subroutine checks whether logical unit number LUNIT is currently*
C* defined to the BUFRLIB software and, if so, returns LUN, IL, and IM	*
C* for that I/O stream.  Otherwise, it checks whether there is space	*
C* within the internal arrays for a new I/O stream and, if so, returns	*
C* within LUN the next available index into those internal arrays.	*
C*									*
C* STATUS  ( LUNIT, LUN, IL, IM )					*
C*									*
C* Input parameters:							*
C*	LUNIT		INTEGER		FORTRAN logical unit number	*
C*									*
C* Output parameters:							*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					arrays for logical unit LUNIT:	*
C*					  0 = LUNIT is not currently	*
C*					      defined *and* could not	*
C*					      subsequently allocate	*
C*					      space for it within	*
C*					      internal arrays		*
C*	IL		INTEGER		Logical unit status indicator:	*
C*					  0 = LUNIT is not currently	*
C*					      defined			*
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

      PARAMETER (NFILES=32)
                                                                        
      COMMON /STBFR/ IOLUN(NFILES),IOMSG(NFILES)
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IF(LUNIT.LE.0 .OR. LUNIT.GT.99) GOTO 900                          
                                                                        
C  CLEAR THE STATUS INDICATORS                                          
C  ---------------------------                                          
                                                                        
      LUN = 0                                                           
      IL  = 0                                                           
      IM  = 0                                                           
                                                                        
C  SEE IF THE UNIT IS DEFINED                                           
C  --------------------------                                           
                                                                        
      DO I=1,NFILES                                                     
      IF(ABS(IOLUN(I)).EQ.LUNIT) LUN = I                                
      ENDDO                                                             
                                                                        
C  IF NOT, CHECK FOR FILE SPACE - RETURN LUN=0 IF NO FILE SPACE         
C  ------------------------------------------------------------         
                                                                        
      IF(LUN.EQ.0) THEN                                                 
         DO I=1,NFILES                                                  
         IF(IOLUN(I).EQ.0) LUN = I                                      
         IF(IOLUN(I).EQ.0) RETURN                                       
         ENDDO                                                          
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  IF FILE DEFINED RETURN STATUSES                                      
C  -------------------------------                                      
                                                                        
      IL = SIGN(1,IOLUN(LUN))                                           
      IM = IOMSG(LUN)                                                   
                                                                        
      RETURN                                                            
900   CALL BORT('STATUS - ILLEGAL UNIT GIVEN')                         
      END                                                               

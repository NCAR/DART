      SUBROUTINE BORT(STR)                                             

C************************************************************************
C* BORT									*
C*									*
C* This subroutine prints (to STDOUT) a given error string and then	*
C* aborts the main program which called the BUFRLIB software.		*
C*									*
C* BORT  ( STR )							*
C*									*
C* Input parameters:							*
C*	STR		CHARACTER*(*)	Error message to be printed to	*
C*					standard output			*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************

      CHARACTER*(*) STR                                                 
      PRINT*
      PRINT*,'************************ABORT**************************'
      PRINT*,STR                                                        
      PRINT*,'************************ABORT**************************'
      PRINT*
      CALLEXIT(49)                                                          
      END                                                               

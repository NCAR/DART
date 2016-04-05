MODULE  module_getsm
!
! DART $Id$
!
CONTAINS

  SUBROUTINE getsmlsm (MMINLU, isltyp, ivgtyp, smcmin, smcmax)
!-----------------------------------------------------------------

    IMPLICIT NONE                                                        
    
    INTEGER :: isltyp, ivgtyp
    REAL :: smcmin,smcmax
    
    INTEGER :: LUMATCH, IINDEX, LC
    
    CHARACTER*4 :: MMINLU, MMINSL
    CHARACTER*4 SLTYPE
    INTEGER :: SLCATS
    INTEGER, PARAMETER :: NSLTYPE=30

    REAL, DIMENSION (1:NSLTYPE) :: BB,DRYSMC,F11, &
        MAXSMC, REFSMC,SATPSI,SATDK,SATDW, WLTSMC,QTZ


    
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
!                                                                       
    OPEN(19, FILE='SOILPARM.TBL',FORM='FORMATTED',STATUS='OLD') 
    
    MMINSL='STAS'                       !oct2
    
    LUMATCH=0                                                       
2000 FORMAT (A4)    

    READ (19,*)
    READ (19,2000,END=2003)SLTYPE
    READ (19,*)SLCATS,IINDEX                                        
    IF(SLTYPE.EQ.MMINSL)THEN                                        
       LUMATCH=1                                                     
    ENDIF
    IF(SLTYPE.EQ.MMINSL)THEN                                    
       DO LC=1,SLCATS                                                
          READ (19,*) IINDEX,BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
               REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
               WLTSMC(LC), QTZ(LC)
       ENDDO
    ENDIF
    
2003 CONTINUE                                                        
    
    CLOSE (19)                                                        
    
    smcmin=drysmc(isltyp)
    smcmax=maxsmc(isltyp)
    
    IF(ivgtyp==1) smcmax = 0.45
    
    IF(LUMATCH.EQ.0)THEN                                            
       PRINT *,'SOIl TEXTURE IN INPUT FILE DOES NOT '
       PRINT *, 'MATCH SOILPARM TABLE'        
       CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
    ENDIF
    
  END SUBROUTINE getsmlsm

  SUBROUTINE getsmruclsm (MMINLU, isltyp, ivgtyp, smcmin, smcmax)

    IMPLICIT NONE
    INTEGER :: isltyp, ivgtyp
    REAL :: smcmin, smcmax
    CHARACTER(len=4) :: mminlu

    INTEGER,   PARAMETER      ::      nsoilclas=19

    REAL  :: LQMA(nsoilclas), LQMI(nsoilclas)

!-- Clapp et al. [1978]
    LQMA =(/0.395, 0.410, 0.435, 0.485, 0.485, 0.451, 0.420,      &
          0.477, 0.476, 0.426, 0.492, 0.482, 0.451, 1.0,        &
          0.20,  0.435, 0.468, 0.200, 0.339/)
     
!-- Carsel and Parrish [1988]
    LQMI=(/0.045, 0.057, 0.065, 0.067, 0.034, 0.078, 0.10,     &
         0.089, 0.095, 0.10,  0.070, 0.068, 0.078, 0.0,      &
         0.004, 0.065, 0.020, 0.004, 0.008/)

    smcmin=lqmi(isltyp)
    smcmax=lqma(isltyp)

  END SUBROUTINE getsmruclsm
    
END MODULE module_getsm

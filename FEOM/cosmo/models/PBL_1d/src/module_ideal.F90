MODULE module_ideal
!
! DART $Id$
!
REAL,PARAMETER :: PI=3.1415926

CONTAINS
  
  FUNCTION glwfunc(time)

    IMPLICIT NONE

    REAL :: glwfunc,time

    IF (MODULO(NINT(time)/3600,24) < 12) THEN
       glwfunc=400.
    ELSE
       glwfunc=400.
    ENDIF

  END FUNCTION glwfunc

  FUNCTION gswfunc(time,albedo)

    IMPLICIT NONE

    REAL :: gswfunc,time,albedo

    gswfunc=MAX(0.,(1.-albedo)*900.*&
         &SIN(MODULO(time,86400.)/86400.*PI*2-.5*PI))

  END FUNCTION gswfunc

  SUBROUTINE mysoil(xland,ts_ref,ps_ref,dtamplitude,qsavail,time,&
       &Tsk,Qsfc)
    
    IMPLICIT NONE
    
    REAL :: xland,ps_ref,ts_ref,dtamplitude,qsavail,&
         &time,Tsk,Qsfc
    
    IF (xland < 1.5) THEN
       tsk=ts_ref+dtamplitude*&
            &SIN(MODULO(time,86400.)/86400.*3.1415926*2.-.5*PI)
       qsfc=qsavail*qsat(ps_ref,tsk)
    ELSE
       tsk=ts_ref
       qsfc=qsat(ps_ref,tsk)
    ENDIF

  END SUBROUTINE mysoil

  FUNCTION qsat(p,t)
    
    IMPLICIT NONE
    
    REAL :: qsat,p,t
    
    REAL :: esat
    
    esat=611.2*EXP(17.67*(t-273.15)/(t-29.65))
    qsat=.622*esat/(p-0.378*esat)
    
  END FUNCTION qsat
  
END MODULE module_ideal

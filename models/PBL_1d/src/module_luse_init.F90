MODULE module_luse_init
!
! DART $Id$
!

CONTAINS

  LOGICAL FUNCTION wrf_dm_on_monitor()
    wrf_dm_on_monitor = .TRUE.
  END FUNCTION wrf_dm_on_monitor

  SUBROUTINE landuse_init(lu_index, snowc, albedo, albbck, mavail, emiss,  &
       znt,Z0,thc,cs,xland, julday, cen_lat, iswater, mminlu, &
       ids, ide, jds, jde, kds, kde,                       &
       ims, ime, jms, jme, kms, kme,                       &
       its, ite, jts, jte, kts, kte                       )

    USE module_wrf_error
    IMPLICIT NONE

!---------------------------------------------------------------------
    INTEGER , INTENT(IN)           :: ids, ide, jds, jde, kds, kde,   &
         ims, ime, jms, jme, kms, kme,   &
         its, ite, jts, jte, kts, kte

    INTEGER , INTENT(IN)           :: iswater, julday
    REAL    , INTENT(IN)           :: cen_lat
    CHARACTER*4, INTENT(IN)        :: mminlu
    REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: lu_index, snowc
    REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(OUT  ) :: albedo, albbck, mavail, emiss, &
         znt, Z0, thc, cs, xland

!---------------------------------------------------------------------
! Local
    CHARACTER*4 LUTYPE
    INTEGER  :: ISICE, LUCATS, LUSEAS
    INTEGER  :: landuse_unit, LS, LC, LI, LUN, NSN=0
    INTEGER  :: i, j, itf, jtf, is, isn
    INTEGER , PARAMETER :: max_cats = 100 , max_seas = 12 
    REAL, DIMENSION( max_cats, max_seas ) :: ALBD, SLMO, SFEM, SFZ0, THERIN, SFHC
    REAL, DIMENSION( max_cats )     :: SCFX
    LOGICAL :: FOUND_LU
!    LOGICAL, EXTERNAL :: wrf_dm_on_monitor

!---------------------------------------------------------------------
    CALL wrf_debug( 100 , 'top of landuse_init' )
    landuse_unit = 29
    IF ( wrf_dm_on_monitor() ) THEN
       OPEN(landuse_unit, FILE='LANDUSE.TBL',FORM='FORMATTED',STATUS='OLD')
    ENDIF
! Determine season (summer=1, winter=2)
    ISN=1                                                            
    IF(JULDAY.LT.105.OR.JULDAY.GT.288)ISN=2                         
    IF(CEN_LAT.LT.0.0)ISN=3-ISN                                   

! Read info from file LANDUSE.TBL
    IF(MMINLU.EQ.'OLD ')THEN
!       ISWATER=7
       ISICE=11 
    ELSE IF(MMINLU.EQ.'USGS')THEN
!       ISWATER=16
       ISICE=24
    ELSE IF(MMINLU.EQ.'SiB ')THEN
!       ISWATER=15
       ISICE=16
    ELSE IF(MMINLU.EQ.'LW12')THEN
!       ISWATER=15
       ISICE=3
    ENDIF
    PRINT *, 'INPUT LANDUSE = ',MMINLU
    FOUND_LU = .FALSE.
1999 CONTINUE                                                      
    if ( wrf_dm_on_monitor() ) then
       READ (landuse_unit,2000,END=2001)LUTYPE                                
       READ (landuse_unit,*)LUCATS,LUSEAS                                    
       FOUND_LU = LUTYPE.EQ.MMINLU
    endif

2000 FORMAT (A4)                                                


    IF(FOUND_LU)THEN                                  
       LUN=LUCATS                                             
       NSN=LUSEAS                                            
       PRINT *, 'LANDUSE TYPE = ',LUTYPE,' FOUND',        &
            LUCATS,' CATEGORIES',LUSEAS,' SEASONS',     &
            ' WATER CATEGORY = ',ISWATER,               &
            ' SNOW CATEGORY = ',ISICE                
    ENDIF
    DO LS=1,LUSEAS                                   
       if ( wrf_dm_on_monitor() ) then
          READ (landuse_unit,*)                                   
       endif
       DO LC=1,LUCATS                               
          IF(FOUND_LU)THEN                  
             IF ( wrf_dm_on_monitor() ) THEN
                READ (landuse_unit,*)LI,ALBD(LC,LS),SLMO(LC,LS),SFEM(LC,LS),        &       
                     SFZ0(LC,LS),THERIN(LC,LS),SCFX(LC),SFHC(LC,LS)       
             ENDIF
             IF(LC.NE.LI)CALL wrf_error_fatal ( 'module_start: MISSING LANDUSE UNIT ' )
          ELSE                                                            
             IF ( wrf_dm_on_monitor() ) THEN
                READ (landuse_unit,*)                                                  
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    IF(NSN.EQ.1) THEN
       ISN = 1
    END IF

    IF(.NOT. FOUND_LU) GOTO 1999


2001 CONTINUE                                                      
    IF(.NOT. FOUND_LU)THEN                                         
       CALL wrf_message ( 'LANDUSE IN INPUT FILE DOES NOT MATCH LUTABLE: TABLE NOT USED' )
    ENDIF

    IF(FOUND_LU)THEN
! Set arrays according to lu_index
       itf = min0(ite, ide)
       jtf = min0(jte, jde)
       DO j = jts, jtf
          DO i = its, itf
             IS=nint(lu_index(i,j))
             IF(IS.LT.0.OR.IS.GT.LUN)THEN                                        
                WRITE ( wrf_err_message , * ) 'module_luse_init: landuse_init: ERROR: LANDUSE OUTSIDE RANGE =',IS,' AT ',I,J
                CALL wrf_error_fatal ( TRIM ( wrf_err_message ) )
             ENDIF
!   SET NO-DATA POINTS (IS=0) TO WATER                                    
             IF(IS.EQ.0)THEN                                                
                IS=ISWATER                                                  
             ENDIF
             ALBBCK(I,J)=ALBD(IS,ISN)/100.                                  
             ALBEDO(I,J)=ALBBCK(I,J)
             THC(I,J)=THERIN(IS,ISN)/100.                               
             cs(i,j)= sfhc(IS,ISN)
             Z0(I,J)=SFZ0(IS,ISN)/100.                                
             ZNT(I,J)=Z0(I,J)
             EMISS(I,J)=SFEM(IS,ISN)                                  
             MAVAIL(I,J)=SLMO(IS,ISN)                                
             IF(IS.NE.ISWATER)THEN                                  
                XLAND(I,J)=1.0                                      
             ELSE                                                 
                XLAND(I,J)=2.0                                    
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    if ( wrf_dm_on_monitor() ) then
       CLOSE (landuse_unit)
    endif
    CALL wrf_debug( 100 , 'returning from of landuse_init' )
    RETURN

  END SUBROUTINE landuse_init

  SUBROUTINE bucket_init(mminlu,lu_index,Maxm,Minm,Erate,&
       ids, ide, jds, jde, kds, kde,   &
       ims, ime, jms, jme, kms, kme,   &
       its, ite, jts, jte, kts, kte)


    IMPLICIT NONE

    INTEGER , INTENT(IN)           :: ids, ide, jds, jde, kds, kde,   &
         ims, ime, jms, jme, kms, kme,   &
         its, ite, jts, jte, kts, kte
    
    
    REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: lu_index
    
    REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(OUT  ) :: &
         &Maxm,Minm,Erate
    
    INTEGER             :: bucket_tbl = 18
    INTEGER , PARAMETER :: max_cats = 100

    REAL, DIMENSION( max_cats )     ::  Maxmm,Minmm,Eratem
    
    CHARACTER(len=4), INTENT(IN)        :: mminlu

    INTEGER  :: LUCATS,LUSEAS
    INTEGER  :: LC, LI, LUN, NSN=0,LUMATCH
    CHARACTER(len=4) ::lutype
    INTEGER :: i,j,itf,jtf,is

    OPEN(bucket_tbl, FILE='BUCKET.TBL',FORM='FORMATTED',STATUS='OLD',ERR=9004)   

    LUMATCH=0
2010 CONTINUE
    READ (bucket_tbl,2000,END=2020)LUTYPE 
    READ (bucket_tbl,*)LUCATS,LUSEAS      
    IF(LUTYPE.EQ.MMINLU)THEN      
       LUN=LUCATS                 
       LUMATCH=1                  
    ENDIF

    DO LC=1,LUCATS                
       IF(LUTYPE.EQ.MMINLU)THEN   
          READ (bucket_tbl,*)LI,maxmm(LC),minmm(LC),eratem(LC) 
          IF(LC.NE.LI) STOP 'MISSING BUCKET.TBL: UNIT 18'
       ELSE                                              
          READ (bucket_tbl,*)                                    
       ENDIF
    ENDDO
    
    GOTO 2010 

2020 CONTINUE 

    IF(LUMATCH.EQ.0)THEN 
       PRINT *,&
            &'LANDUSE IN INPUT FILE DOES NOT MATCH DATA IN BUCKET TBL'
       print*,LUTYPE,MMINLU
       STOP 'INCONSISTENT OR MISSING BUCKET.TBL FILE'                
    ELSE
       itf = min0(ite, ide)
       jtf = min0(jte, jde)
       DO j = jts, jtf
          DO i = its, itf
             is=NINT(lu_index(i,j))
             maxm(i,j)=maxmm(is)
             minm(i,j)=minmm(is)
             erate(i,j)=eratem(is)
          ENDDO
       ENDDO
    ENDIF

    IF ( wrf_dm_on_monitor() ) THEN
       CLOSE(bucket_tbl)
    ENDIF

    RETURN

2000 FORMAT (A4)

9004 PRINT *, 'ERRPR OPENING BUCKET MODEL DATA FILE FROM UNIT 18'
    STOP '9004 IN  module_luse_init'                 

  END SUBROUTINE bucket_init


END MODULE module_luse_init

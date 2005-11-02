PROGRAM driver

! program to drive wrf1d for compatibility with DART

  USE module_wrf

  IMPLICIT NONE

  CHARACTER(len=120) :: &
       &namelistfile='wrf1d_namelist.input',&
       &outlogfile='wrf1d_log.out'
  INTEGER  :: unit_nml=151,unit_log=152
  LOGICAL :: is_it_there = .FALSE.
  INTEGER :: dart_days, dart_seconds

  INTEGER :: wrf_rnd_seed,itime ! equivalent to itimestep in module_wrf.F

  LOGICAL :: allocate_wrf = .TRUE.

  INQUIRE ( FILE = namelistfile , EXIST = is_it_there )

  IF ( is_it_there ) THEN
     OPEN ( FILE   = namelistfile     , UNIT   =  unit_nml        ,&
          & STATUS = 'OLD'            , FORM   = 'FORMATTED'      ,&
          & ACTION = 'READ'           , ACCESS = 'SEQUENTIAL'     )
  ELSE
     PRINT '(A)','Could not find the namelist: ',namelistfile
     STOP 'No_namelist_found'
  ENDIF

  OPEN ( FILE   = outlogfile        , UNIT   =  unit_log        ,&
       & STATUS = 'REPLACE'         , FORM   = 'FORMATTED'      ,&
       & ACTION = 'WRITE'           , ACCESS = 'SEQUENTIAL'     )

  CALL do_namelist_wrf1d(unit_nml,unit_log)

  CLOSE(unit_log)
  
  CALL static_init_wrf()
  
  wrf_rnd_seed = rnd_seed_val

  CALL init_wrf(allocate_wrf,wrf_rnd_seed)
  
  allocate_wrf=.FALSE.

  IF (init_gen) THEN 
     OPEN(ncunit,file=out_gen_file)
  ENDIF

  DO itime=1,ntime+1
     IF (init_gen) THEN
        dart_days=INT(REAL(REAL(itime-1)*dt+start_seconds)/86400.)
        dart_seconds=NINT(REAL(itime-1)*dt-REAL(dart_days*86400))+&
             &start_seconds
     ELSE
        dart_days=0
        dart_seconds=(itime-1)*dt
     ENDIF
     CALL wrf(dart_seconds,dart_days)
     CALL output_wrf_profiles()
  ENDDO
  
END PROGRAM driver

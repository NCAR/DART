! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

PROGRAM driver

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! 1d WRF incorporates original WRF PBL routines for idealized runs.
! It can be initialized with idealized profiles of wind, temperature,
! and mixing ratio or with profiles derived from 3d WRF runs.
! All input is driven by parameters specified in wrf1d_namelist.input
! file. This code is compatible with DART software used for Data
! Assimilatiuon studies.

! The code currently has a problem when compiled with ifc
! as it crashes when calculating MINVAL of an array for unknown reason.
! It runs OK when compiled with pgf90


! The code is available as is but had been checked quite thoroughlly.
! It is of course free to use and modify but I would appreciate
! if you let me know your interest, acknowledge the authorship,
! and, even better if you contribute
! to the development and make the modifications available.

! File README gives more details on input parameters in the namelist.
! File PBLcolumn.pdf gives some general information on the model
! and its coupling with DART.

! 1d WRF was originally developed by Mariusz Pagowski 
! NOAA/GSD/CIRA, Boulder, CO, Mariusz.Pagowski@noaa.gov
! Josh Hacker,NCAR/RAL, hacker@ucar.edu
! contributed reading routines to real data and 3dWRF output and 
! interface to DART and made the code more flexible

! Based on WRFV2.1.2

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
  
  CALL static_init_wrf(allocate_wrf)
  
  wrf_rnd_seed = rnd_seed_val
  allocate_wrf=.FALSE.

  CALL init_wrf(wrf_rnd_seed)
  

! now this is done in init_wrf
!  IF (init_f) THEN 
!     OPEN(ncunit,file=out_f_file)
!  ENDIF

  DO itime=1,ntime+1
     IF (init_f) THEN
        dart_days=INT(REAL(REAL(itime-1)*dt+start_seconds)/86400.)
        dart_seconds=NINT(REAL(itime-1)*dt-REAL(dart_days*86400))+&
             &start_seconds
     ELSE
        dart_days=0
        dart_seconds=(itime-1)*dt
     ENDIF
     if ( forecast_length == 0 ) stop '0 time steps'
     CALL wrf(dart_seconds,dart_days)
     CALL output_wrf_profiles()
  ENDDO
  
END PROGRAM driver

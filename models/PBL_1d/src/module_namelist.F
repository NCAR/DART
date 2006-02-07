MODULE module_namelist

  USE module_model_constants

  IMPLICIT NONE

! record 1
  LOGICAL :: init_f, rnd_init
  CHARACTER(len=3)  :: init_f_type
  CHARACTER(len=15) :: pbltype,sfctype,surface_physics,bucket_model
  REAL :: dt,deep_soil_moisture
  INTEGER :: nz
  INTEGER :: n_moist,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG
  INTEGER :: ifsnow,isfflx
  REAL :: pblh_ref
  CHARACTER(len=120) :: indir,outdir

! record 2
  CHARACTER(len=120) :: init_f_file,out_f_file,init_soil_file,&
                        init_flux_file, init_smos_file
  LOGICAL :: output_state_vector
  INTEGER :: start_year_f, start_month_f, start_day_f, &
             start_hour_f, interval_f, interval_flux, &
             interval_soil, interval_smos, start_minute_f, &
             start_forecast, forecast_length
  INTEGER :: splineinterval, splineinterval_flux, splineinterval_smos
  REAL :: outfinterval
  INTEGER :: rnd_seed_val
  REAL :: z_g

!record 3

  REAL :: odayfraction,totdayfraction
  INTEGER :: outinterval
  REAL :: dtamplitude_ref,dtdz_ref,u_g_ref,v_g_ref,qsrat_ref,&
       rland_ref,mavail_ref

  INTEGER :: ivgtyp_ref,isltyp_ref,lu_index_ref,julday_ref
  CHARACTER(len=4) :: mminlu_ref

  REAL :: vegfra_ref,zo_ref,emiss_ref,thc_ref,cs_ref,albedo_ref,&
       &maxm_ref,minm_ref,erate_ref,&
       ts_ref,tmn_ref,ps_ref,prate_ref,lat_ref,lon_ref

! record 4 - parameter estimation
  INTEGER                :: num_est_params
  INTEGER, dimension(10) :: est_param_types
  REAL   , dimension(10) :: pert_init_sd
  REAL   , dimension(10) :: pert_param_sd
  REAL   , dimension(10) :: pert_param_min
  REAL   , dimension(10) :: pert_param_max
  LOGICAL                :: maintain_initial_spread
  CHARACTER(len=4), dimension(10) :: dist_shape

! additional variables

  REAL :: timeo,timetot

CONTAINS

  SUBROUTINE do_namelist_wrf1d(unit_nml,logfileunit)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)  :: logfileunit, unit_nml
    
    NAMELIST /RECORD1/ init_f,init_f_type,pbltype,sfctype,surface_physics,&
         bucket_model,&
         dt,deep_soil_moisture,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG,&
         ifsnow,isfflx,pblh_ref,indir,outdir, nz, rnd_init
    
    NAMELIST /RECORD2/ init_f_file,out_f_file,init_soil_file,&
                        init_flux_file, init_smos_file, &
                        output_state_vector,&
                        start_year_f, start_month_f, start_day_f, &
                        start_hour_f, start_minute_f, &
                        interval_f, interval_flux, interval_soil, &
                        interval_smos, &
                        start_forecast,forecast_length, &
                        splineinterval,splineinterval_flux, &
                        splineinterval_smos,&
                        outfinterval,&
                        rnd_seed_val,z_g

    NAMELIST /RECORD3/ odayfraction,totdayfraction,outinterval,&
         dtamplitude_ref,dtdz_ref,u_g_ref,v_g_ref,qsrat_ref,rland_ref,&
         mminlu_ref,julday_ref,lu_index_ref,mavail_ref,ivgtyp_ref,&
         isltyp_ref,vegfra_ref,zo_ref,pblh_ref,&
         emiss_ref,thc_ref,cs_ref,albedo_ref,&
         &maxm_ref,minm_ref,erate_ref,&
         &ts_ref,tmn_ref,ps_ref,prate_ref,&
         lat_ref,lon_ref
    
    NAMELIST /RECORD4/ num_est_params, est_param_types, pert_param_sd, &
         pert_init_sd, maintain_initial_spread, dist_shape, pert_param_min, &
         pert_param_max

! Local variables.
    
    LOGICAL :: is_it_there = .FALSE.
    
    
! Does the file exist?
    
       
!  The file exists, get a unit number.
       
!       OPEN ( FILE   = namelistfile     , UNIT   =  unit_nml        ,&
!            & STATUS = 'OLD'            , FORM   = 'FORMATTED'      ,&
!            & ACTION = 'READ'           , ACCESS = 'SEQUENTIAL'     )
       
!  A couple of defaults

       rnd_init = .false.
       num_est_params = 0
       est_param_types(:) = 0
       pert_init_sd(:) = 0.0
       pert_param_sd(:) = 0.0
       pert_param_min(:) = 0.0
       pert_param_max(:) = 0.0
       maintain_initial_spread = .false.
       dist_shape(:) = 'norm'

!  File is opened, so read it.
       
       READ (unit_nml , RECORD1 )
       READ (unit_nml , RECORD2 )
       READ (unit_nml , RECORD3 )
       READ (unit_nml , RECORD4 )
       
       n_moist=MAX(P_QV,P_QC,P_QR,P_QI,P_QS,P_QG)

       IF (init_f) THEN
          timeo=start_forecast
       ELSE
          timeo=86400.*odayfraction
          timetot=86400.*totdayfraction
       ENDIF

       IF (dtdz_ref < -999. ) dtdz_ref=g/cp

       prate_ref=prate_ref/3600. ! convert hours to secs

! Silly error checking - don't know where else to put it
    IF ( init_f_type == 'OBS' .and. start_forecast > 0 ) THEN
       print*,'OBS init must be with start_forecast = 0'
       stop 'module_namelist'
    ENDIF
    
! Record the namelist to the logfile
    write(logfileunit,nml=record1)
    write(logfileunit,nml=record2)
    write(logfileunit,nml=record3)
    write(logfileunit,nml=record4)

    init_f_file=TRIM(indir)//'/'//init_f_file
    out_f_file=TRIM(outdir)//'/'//out_f_file
    init_soil_file=TRIM(indir)//'/'//init_soil_file
    init_flux_file=TRIM(indir)//'/'//init_flux_file
    init_smos_file=TRIM(indir)//'/'//init_smos_file

  END SUBROUTINE do_namelist_wrf1d

END MODULE module_namelist


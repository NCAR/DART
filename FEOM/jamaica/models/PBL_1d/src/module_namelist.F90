MODULE module_namelist

  USE module_model_constants

  IMPLICIT NONE

! record 1
  LOGICAL :: init_f 
  CHARACTER(len=3)  :: init_f_type, force_f_type
  CHARACTER(len=15) :: bl_pbl_physics,sf_sfclay_physics,sf_surface_physics,&
       &bucket_model
  INTEGER :: rnd_init, rnd_force     
  REAL :: dt,deep_soil_moisture
  INTEGER :: nz
  INTEGER :: n_moist,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG
  INTEGER :: ifsnow,isfflx
  REAL :: pblh_ref
  CHARACTER(len=120) :: indir,outdir

! record 2
  CHARACTER(len=120) :: init_f_file,out_f_file,init_soil_file,&
                        init_flux_file, init_smos_file, uvg_file,&
                        &sfc_file,&
                        eofs_file_init, eofs_file_forc
  LOGICAL :: output_state_vector, force_uvg, force_flux, calc_sfc, &
             calc_ust, t_advection, qv_advection 
  INTEGER :: start_year_f, start_month_f, start_day_f, &
             start_hour_f, interval_f, interval_flux, &
             interval_soil, interval_smos, interval_sfc,&
             &start_minute_f, &
             start_forecast, forecast_length, interval_uvg,n_eo
  INTEGER :: splineinterval, splineinterval_flux, splineinterval_smos,&
       &splineinterval_sfc, splineinterval_advection
  REAL :: outfinterval, scales
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
       ts_ref,tmn_ref,ps_ref,prate_ref,lat_ref,lon_ref,&
       &hflux_ref,qvflux_ref

! additional variables

  REAL :: timeo,timetot

CONTAINS

  SUBROUTINE do_namelist_wrf1d(unit_nml,logfileunit)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)  :: logfileunit, unit_nml
    
    NAMELIST /RECORD1/ init_f,init_f_type,force_f_type,bl_pbl_physics,&
         &sf_sfclay_physics,sf_surface_physics,&
         bucket_model,&
         dt,deep_soil_moisture,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG,&
         ifsnow,isfflx,pblh_ref,indir,outdir, nz, rnd_init, rnd_force
    
    NAMELIST /RECORD2/ init_f_file,out_f_file,init_soil_file,&
                        init_flux_file, init_smos_file, uvg_file, &
                        eofs_file_init, eofs_file_forc,&
                        sfc_file,output_state_vector, force_uvg, &
                        &calc_sfc, calc_ust, t_advection, qv_advection, &
                        start_year_f, start_month_f, start_day_f, &
                        start_hour_f, start_minute_f, &
                        interval_f, interval_flux, interval_soil, &
                        interval_smos, interval_uvg, &
                        &interval_sfc,&
                        start_forecast,forecast_length, &
                        splineinterval,splineinterval_flux, &
                        splineinterval_smos,splineinterval_sfc,&
                        splineinterval_advection, &
                        outfinterval,&
                        rnd_seed_val,z_g,n_eo,scales

    NAMELIST /RECORD3/ odayfraction,totdayfraction,outinterval,&
         dtamplitude_ref,dtdz_ref,u_g_ref,v_g_ref,qsrat_ref,rland_ref,&
         mminlu_ref,julday_ref,lu_index_ref,mavail_ref,ivgtyp_ref,&
         isltyp_ref,vegfra_ref,zo_ref,pblh_ref,&
         emiss_ref,thc_ref,cs_ref,albedo_ref,&
         &maxm_ref,minm_ref,erate_ref,&
         &ts_ref,tmn_ref,ps_ref,prate_ref,&
         hflux_ref,qvflux_ref,lat_ref,lon_ref
    
! Local variables.
    
    LOGICAL :: is_it_there = .FALSE.
    
! Any defaults?  Should be a long list but starting here...
       force_uvg = .true.
       t_advection = .false.
       qv_advection = .false.
       calc_sfc = .false.
       calc_ust = .false.
       splineinterval = 3600.0
       splineinterval_advection = 3600.0
       splineinterval_flux = 3600.0
       splineinterval_smos = 3600.0
       splineinterval_sfc = 3600.0
       
!  File is opened, so read it.
       
       READ (unit_nml , RECORD1 )
       READ (unit_nml , RECORD2 )
       READ (unit_nml , RECORD3 )

       n_moist=MAX(P_QV,P_QC,P_QR,P_QI,P_QS,P_QG)

       IF  (sf_surface_physics=='FLUXFORCE') THEN
          force_flux=.TRUE.
       ELSE
          force_flux=.FALSE.
       ENDIF

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

    IF ( splineinterval_smos > interval_smos ) THEN
       PRINT *,'Missing precip accumulated at ',interval_smos,&
            &'intervals.'
       PRINT *,'Reset splineinterval_smos =< interval_smos ',&
            &'in the namelist. Stopping'
       STOP
    ENDIF

    IF ( init_f_type /= 'WRF' .and. t_advection ) then
       PRINT*, 'If including T advection, only WRF initialization valid'
       STOP
    ENDIF
    IF ( init_f_type /= 'WRF' .and. qv_advection ) then
       PRINT*, 'If including QV advection, only WRF initialization valid'
       STOP
    ENDIF
    
! Record the namelist to the logfile
    write(logfileunit,nml=record1)
    write(logfileunit,nml=record2)
    write(logfileunit,nml=record3)

    init_f_file=TRIM(indir)//'/'//init_f_file
    out_f_file=TRIM(outdir)//'/'//out_f_file
    init_soil_file=TRIM(indir)//'/'//init_soil_file
    init_flux_file=TRIM(indir)//'/'//init_flux_file
    init_smos_file=TRIM(indir)//'/'//init_smos_file
    uvg_file=TRIM(indir)//'/'//uvg_file
    eofs_file_init=TRIM(indir)//'/'//eofs_file_init
    eofs_file_forc=TRIM(indir)//'/'//eofs_file_forc
    sfc_file=TRIM(indir)//'/'//sfc_file

  END SUBROUTINE do_namelist_wrf1d

END MODULE module_namelist


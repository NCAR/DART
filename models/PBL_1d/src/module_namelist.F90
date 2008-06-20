MODULE module_namelist

  USE module_model_constants

  IMPLICIT NONE

! record 1
  LOGICAL :: init_f 
  CHARACTER(len=3)  :: init_f_type, force_f_type, ra_type
  CHARACTER(len=4)  :: nudge_f_type
  CHARACTER(len=15) :: bl_pbl_physics,sf_sfclay_physics,sf_surface_physics,&
       &bucket_model,ra_lw_physics,ra_sw_physics, mp_physics
  INTEGER :: rnd_init, rnd_force     
  REAL :: dt,deep_soil_moisture, radt, nudge_f_coeff, nudge_f_hwpress
  INTEGER :: nz
  INTEGER :: n_moist,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG
  INTEGER :: ifsnow,isfflx
  REAL :: pblh_ref
  CHARACTER(len=120) :: indir,outdir
  CHARACTER(len=120) :: gridfile !...RT

! record 2
  CHARACTER(len=120) :: init_f_file,out_f_file,init_soil_file,&
                        init_flux_file, init_smos_file, uvg_file,&
                        &sfc_file,&
                        eofs_file_init, eofs_file_forc
  LOGICAL :: output_state_vector, force_uvg, force_flux, calc_sfc, &
             calc_ust, t_advection, qv_advection, u_advection, &
             qr_advection, qi_advection, qg_advection, qc_advection, &
             fixed_timeofday, rotate_sfc_winds, force_soil
  INTEGER :: start_year_f, start_month_f, start_day_f, &
             start_hour_f, interval_f, interval_flux, &
             interval_soil, interval_smos, interval_sfc,&
             &start_minute_f, &
             start_forecast, forecast_length, interval_uvg,n_eo
  INTEGER :: splineinterval, splineinterval_flux, splineinterval_smos,&
       &splineinterval_sfc, splineinterval_advection, splineinterval_soil
  REAL :: outfinterval, scales
  INTEGER :: rnd_seed_val
  REAL :: z_g

!record 3

  REAL :: odayfraction,totdayfraction
  INTEGER :: outinterval
  REAL :: dtamplitude_ref,dtdz_ref,u_g_ref,v_g_ref,qsrat_ref,&
       rland_ref,mavail_ref,timeofday_ref

  INTEGER :: ivgtyp_ref,isltyp_ref,lu_index_ref,julday_ref
  CHARACTER(len=4) :: mminlu_ref

  REAL :: vegfra_ref,zo_ref,emiss_ref,thc_ref,cs_ref,albedo_ref,&
       &maxm_ref,minm_ref,erate_ref,&
       ts_ref,tmn_ref,ps_ref,prate_ref,lat_ref,lon_ref,&
       &hflux_ref,qvflux_ref

!record 4 !! ...added by RT

  LOGICAL :: forc_stochastic_cloud 
  CHARACTER(len=120) :: cld_file
  CHARACTER(len=3) :: cld_sequence
  INTEGER :: nbrealizations
  REAL :: spin_up_period, cloud_seq_start, cloud_realization_length

! additional variables

  REAL :: timeo,timetot

CONTAINS

  SUBROUTINE do_namelist_wrf1d(unit_nml,logfileunit)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)  :: logfileunit, unit_nml
    
    NAMELIST /RECORD1/ init_f,init_f_type,force_f_type,&
        &ra_type,radt,ra_lw_physics,ra_sw_physics,&
        &bl_pbl_physics,sf_sfclay_physics,sf_surface_physics,&
         mp_physics, bucket_model,&
         dt,deep_soil_moisture,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG,&
         ifsnow,isfflx,pblh_ref,indir,outdir, nz, rnd_init, rnd_force, &
         nudge_f_type, nudge_f_coeff, nudge_f_hwpress, &
         gridfile, fixed_timeofday
    
    NAMELIST /RECORD2/ init_f_file,out_f_file,init_soil_file,&
                        init_flux_file, init_smos_file, uvg_file, &
                        eofs_file_init, eofs_file_forc,&
                        sfc_file,output_state_vector, force_uvg, &
                        calc_sfc, calc_ust, t_advection, qv_advection, &
                        qc_advection, qr_advection,  &
                        qi_advection, qg_advection,  &
                        u_advection, &
                        start_year_f, start_month_f, start_day_f, &
                        start_hour_f, start_minute_f, &
                        interval_f, interval_flux, interval_soil, &
                        interval_smos, interval_uvg, &
                        interval_sfc,rotate_sfc_winds,&
                        start_forecast,forecast_length, &
                        splineinterval,splineinterval_flux, &
                        splineinterval_smos,splineinterval_sfc,&
                        splineinterval_advection, splineinterval_soil, &
                        outfinterval, force_soil, &
                        rnd_seed_val,z_g,n_eo,scales

    NAMELIST /RECORD3/ odayfraction,totdayfraction,outinterval,&
         dtamplitude_ref,dtdz_ref,u_g_ref,v_g_ref,qsrat_ref,rland_ref,&
         mminlu_ref,julday_ref,lu_index_ref,mavail_ref,ivgtyp_ref,&
         isltyp_ref,vegfra_ref,zo_ref,pblh_ref,&
         emiss_ref,thc_ref,cs_ref,albedo_ref,&
         maxm_ref,minm_ref,erate_ref,timeofday_ref,&
         ts_ref,tmn_ref,ps_ref,prate_ref,&
         hflux_ref,qvflux_ref,lat_ref,lon_ref
    
    NAMELIST /RECORD4/ forc_stochastic_cloud,cld_file,cld_sequence,&
      nbrealizations,spin_up_period, cloud_seq_start, cloud_realization_length


! Local variables.
    
    LOGICAL :: is_it_there = .FALSE.
    
! Any defaults?  Should be a long list but starting here...
       mp_physics = 'NONE'
       nudge_f_type = 'NONE'
       rotate_sfc_winds = .true.
       nudge_f_coeff = 0.0
       P_QV = 1
       P_QC = 1
       P_QR = 1
       P_QI = 1
       P_QS = 1
       P_QG = 1
       force_uvg = .true.
       force_soil = .false.
       t_advection = .false.
       qv_advection = .false.
       qc_advection = .false.
       qr_advection = .false.
       qi_advection = .false.
       qg_advection = .false.
       u_advection = .false.
       calc_sfc = .false.
       calc_ust = .false.
       splineinterval = 3600.0
       splineinterval_advection = 3600.0
       splineinterval_flux = 3600.0
       splineinterval_smos = 3600.0
       splineinterval_sfc = 3600.0
       splineinterval_soil = 3600.0
       forc_stochastic_cloud = .false. !...added by RT
       spin_up_period = 0.0 !...added by RT
       cloud_realization_length = 86400.0
       gridfile = 'grid_wrf1d.ascii' !...added by RT
       cloud_seq_start = 43200.0
       fixed_timeofday = .false.
       timeofday_ref = 0
       
!  File is opened, so read it.
       
       READ (unit_nml , RECORD1 )
       READ (unit_nml , RECORD2 )
       READ (unit_nml , RECORD3 )
       READ (unit_nml , RECORD4 ) !...added by RT

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

! random seed should be negative
       if ( rnd_seed_val > 0 ) rnd_seed_val = -rnd_seed_val

! Silly error checking - don't know where else to put it
    IF ( init_f_type == 'OBS' .and. start_forecast > 0 ) THEN
       PRINT*,'OBS init must be with start_forecast = 0'
       STOP 'module_namelist'
    ENDIF

    IF ( trim(nudge_f_type) /= 'NONE' .and. &
         trim(nudge_f_type) .ne. trim(force_f_type) ) then
       PRINT*,'nudge_f_type and force_f_type must be the same if using nudging:'
       PRINT*, nudge_f_type, force_f_type 
       STOP 'module_namelist'
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
    IF ( init_f_type /= 'WRF' .and. u_advection ) then
       PRINT*, 'If including U,V advection, only WRF initialization valid'
       STOP
    ENDIF

! simple check of options for stochastic cloud and interactive radiation...RT
    IF (forc_stochastic_cloud .and. ra_type /= 'INT') then
       PRINT*, 'Stochastic cloud forcing activated only when interactive radiation'
       STOP
    ELSEIF ( forc_stochastic_cloud ) then
       PRINT*, 'Stochastic cloud forcing activated'
    ENDIF
      
! for now, need advection of species to run Lin scheme
    IF ( trim(mp_physics) /= 'NONE' ) then
       IF ( P_QC < 2 ) stop 'Need QC allocations -- activate species'
       IF ( P_QR < 2 ) stop 'Need QR allocations -- activate species'
       IF ( P_QI < 2 ) stop 'Need QI allocations -- activate species'
       IF ( P_QG < 2 ) stop 'Need QG allocations -- activate species'
    ENDIF

! Record the namelist to the logfile
    write(logfileunit,nml=record1)
    write(logfileunit,nml=record2)
    write(logfileunit,nml=record3)
    write(logfileunit,nml=record4) !...added by RT

    init_f_file=TRIM(indir)//'/'//init_f_file
    out_f_file=TRIM(outdir)//'/'//out_f_file
    init_soil_file=TRIM(indir)//'/'//init_soil_file
    init_flux_file=TRIM(indir)//'/'//init_flux_file
    init_smos_file=TRIM(indir)//'/'//init_smos_file
    uvg_file=TRIM(indir)//'/'//uvg_file
    eofs_file_init=TRIM(indir)//'/'//eofs_file_init
    eofs_file_forc=TRIM(indir)//'/'//eofs_file_forc
    sfc_file=TRIM(indir)//'/'//sfc_file

    cld_file=TRIM(indir)//'/'//cld_file !...added by RT

  END SUBROUTINE do_namelist_wrf1d

END MODULE module_namelist


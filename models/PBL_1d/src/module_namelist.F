MODULE module_namelist

  USE module_model_constants

  IMPLICIT NONE

!  INTEGER, PARAMETER :: nz=40,nz_gen=27,ns_gen=5,ntimes_gen=16
!  INTEGER, PARAMETER :: nz=33,nz_gen=27,ns_gen=5,ntimes_gen=16 !WRFRT
  INTEGER, PARAMETER :: nz=33,nz_gen=34,ns_gen=4,ntimes_gen=33 !BAMEX
! the above is hardcoded for now - ie. programs require recompilation

! record 1
  LOGICAL :: init_gen
  CHARACTER(len=15) :: pbltype,sfctype,surface_physics
  REAL :: dt,deep_soil_moisture
  INTEGER :: n_moist,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG
  INTEGER :: ifsnow,isfflx
  REAL :: pblh_ref
  CHARACTER(len=120) :: indir,outdir

! record 2
  LOGICAL :: bamex
  CHARACTER(len=120) :: init_gen_file,out_gen_file
  LOGICAL :: output_state_vector
  INTEGER :: start_hour_wrf, hour_interval_wrf, start_forecast,&
       &forecast_length
  INTEGER :: splineinterval
  REAL :: outgeninterval
  INTEGER :: rnd_seed_val
  REAL :: z_g

!record 3

  REAL :: odayfraction,totdayfraction
  INTEGER :: outinterval
  REAL :: dtamplitude_ref,dtdz_ref,u_g_ref,v_g_ref,qsrat_ref,&
       &rland_ref,mavail_ref

  INTEGER :: ivgtyp_ref,isltyp_ref,lu_index_ref,julday_ref
  CHARACTER(len=4) :: mminlu_ref

  REAL :: vegfra_ref,zo_ref,emiss_ref,thc_ref,albedo_ref,&
       &ts_ref,tmn_ref,ps_ref,lat_ref,lon_ref

! record 4 - parameter estimation
  INTEGER                :: num_est_params
  INTEGER, dimension(10) :: est_param_types
  REAL   , dimension(10) :: pert_init_sd
  REAL   , dimension(10) :: pert_param_sd
  LOGICAL                :: maintain_initial_spread
  CHARACTER(len=4), dimension(10) :: dist_shape

! additional variables

  REAL :: timeo,timetot

CONTAINS

  SUBROUTINE do_namelist_wrf1d(unit_nml,logfileunit)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)  :: logfileunit, unit_nml
    
    NAMELIST /RECORD1/ init_gen,pbltype,sfctype,surface_physics,&
         &dt,deep_soil_moisture,P_QV,P_QC,P_QR,P_QI,P_QS,P_QG,&
         &ifsnow,isfflx,pblh_ref,indir,outdir
    
    NAMELIST /RECORD2/ bamex,init_gen_file,out_gen_file,&
         &output_state_vector,&
         &start_hour_wrf,hour_interval_wrf,&
         &start_forecast,forecast_length,splineinterval,outgeninterval,&
         &rnd_seed_val,z_g

    NAMELIST /RECORD3/ odayfraction,totdayfraction,outinterval,&
         &dtamplitude_ref,dtdz_ref,u_g_ref,v_g_ref,qsrat_ref,rland_ref,&
         &mminlu_ref,julday_ref,lu_index_ref,mavail_ref,ivgtyp_ref,&
         &isltyp_ref,vegfra_ref,zo_ref,pblh_ref,&
         &emiss_ref,thc_ref,albedo_ref,ts_ref,tmn_ref,ps_ref,&
         &lat_ref,lon_ref
    
    NAMELIST /RECORD4/ num_est_params, est_param_types, pert_param_sd, &
         &pert_init_sd, maintain_initial_spread, dist_shape

!  Local variables.
    
    LOGICAL :: is_it_there = .FALSE.
    
    
! Does the file exist?
    
       
!  The file exists, get a unit number.
       
!       OPEN ( FILE   = namelistfile     , UNIT   =  unit_nml        ,&
!            & STATUS = 'OLD'            , FORM   = 'FORMATTED'      ,&
!            & ACTION = 'READ'           , ACCESS = 'SEQUENTIAL'     )
       
!  A couple of defaults

       num_est_params = 0
       est_param_types(:) = 0
       pert_init_sd(:) = 0.0
       pert_param_sd(:) = 0.0
       maintain_initial_spread = .false.
       dist_shape(:) = 'norm'

!  File is opened, so read it.
       
       READ (unit_nml , RECORD1 )
       READ (unit_nml , RECORD2 )
       READ (unit_nml , RECORD3 )
       READ (unit_nml , RECORD4 )
       
       n_moist=MAX(P_QV,P_QC,P_QR,P_QI,P_QS,P_QG)

       IF (init_gen) THEN
          timeo=start_forecast*3600
          outgeninterval=3600.*outgeninterval
       ELSE
          timeo=86400.*odayfraction
          timetot=86400.*totdayfraction
       ENDIF

       IF (dtdz_ref < -999. ) dtdz_ref=g/cp

    
! Record the namelist to the logfile
    write(logfileunit,nml=record1)
    write(logfileunit,nml=record2)
    write(logfileunit,nml=record3)
    write(logfileunit,nml=record4)

    init_gen_file=TRIM(indir)//'/'//init_gen_file
    out_gen_file=TRIM(outdir)//'/'//out_gen_file

  END SUBROUTINE do_namelist_wrf1d

END MODULE module_namelist


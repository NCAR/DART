MODULE module_init_soil_real_fluxforce
!
! DART $Id$
!

! Module takes care of driving the soil initialization for the ideal
! case.

USE module_ideal,           ONLY: qsat

USE module_luse_init,       ONLY: landuse_init,bucket_init

USE module_model_constants, only: DEGRAD, eomeg

USE module_namelist,        ONLY: lat_ref,mminlu_ref, &
                                  lu_index_ref, ivgtyp_ref, isltyp_ref, &
                                  vegfra_ref, rland_ref, albedo_ref, &
                                  zo_ref, emiss_ref,&
                                  ps_ref,  mavail_ref, &
                                  ts_ref,tmn_ref, &
                                  deep_soil_moisture

USE module_soil_pre,        ONLY: init_soil_depth_0,init_soil_0_real

IMPLICIT NONE

private

PUBLIC::   init_soil_real_fluxforce

CONTAINS

  SUBROUTINE init_soil_real_fluxforce(julday, lu_index, ivgtyp, &
                           isltyp, vegfra, iswater, &
                           cent_lat, cor, mminlu, &
                           snowc, albedo, albbck,  mavail, emiss, &
                           maxm,minm,erate,&
                           znt, z0, thc, cs, xland, xice, seamask, &
                           snow, snowh, canwat, smstav, smstot, &
                           sfcrunoff, udrunoff, acsnow, acsnom, &
                           fndsoilw, fndsnowh, smfr3d,  &
                           landmask_input, time, &
                           tsk, qsfc, tmn, tslb, smois, sh2o, &
                           sst_input, flag_sst, &
                           num_soil_layers,  num_st_levels_input , &
                           num_sm_levels_input, &
                           st_levels_input, sm_levels_input, &
                           zs, dzs, st_input, sm_input, keepfr3dflag, &
                           smcmin, smcmax, restart, &
                           ns_f, zs_f, &
                           itime_f, vegfra_f, ivgtyp_f, isltyp_f, &
                           lu_index_f, &
                           tsk_init_f, tmn_init_f, qsfc_init_f, &
                           smois_init_f, tslb_init_f, &
                           ts_init_f, qvs_init_f, &
                           ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte )

   IMPLICIT NONE

   INTEGER, INTENT(in)         :: ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  num_soil_layers,  num_st_levels_input, &
                                  num_sm_levels_input, itime_f, ns_f,&
                                  julday
   INTEGER, DIMENSION(:), INTENT( in) :: lu_index_f, isltyp_f,ivgtyp_f

   INTEGER, INTENT(out)        :: iswater
   INTEGER, DIMENSION(:,:), INTENT(out) ::  ivgtyp, isltyp
   INTEGER, DIMENSION(:), INTENT(out) :: st_levels_input, sm_levels_input

   REAL,    INTENT( in)            :: time
   REAL, DIMENSION(:), INTENT( in) :: tsk_init_f, tmn_init_f, qsfc_init_f, &
                                      vegfra_f

   REAL, DIMENSION(:), INTENT( in) :: ts_init_f,qvs_init_f

   REAL, DIMENSION(:,:), INTENT( in) :: smois_init_f, tslb_init_f

   REAL,    INTENT(out)        :: cent_lat, cor
   REAL,    INTENT(out)        :: smcmin, smcmax
   REAL, DIMENSION(:), INTENT(in) :: zs_f
   REAL, DIMENSION(:,:), INTENT(out) :: lu_index, vegfra, snowc, &
                                        albedo, albbck, mavail, emiss, &
                                        &maxm,minm,erate,&
                                        znt, z0, thc, cs,xland, seamask, &
                                        landmask_input, tsk, qsfc, tmn, &
                                        sst_input, snow, snowh, canwat, &
                                        smstav, smstot, sfcrunoff, udrunoff, &
                                        acsnow, acsnom, xice 
   REAL, DIMENSION(:,:,:), INTENT(out) :: tslb, smois, sh2o, st_input, &
                                          sm_input, smfr3d, keepfr3dflag
   REAL, DIMENSION(:),     INTENT(out) :: zs, dzs
   CHARACTER(len=4), INTENT(inout) :: mminlu
   LOGICAL,  INTENT(in)         :: flag_sst, fndsoilw, fndsnowh, restart

! local
   INTEGER                      :: k ,kns_f_t,kns_f_m

   cent_lat=lat_ref*DEGRAD
   cor=2.*eomeg*SIN(cent_lat)

   lu_index=  REAL(lu_index_f(itime_f))
   vegfra =   vegfra_f(itime_f)
   ivgtyp =   ivgtyp_f(itime_f)
   isltyp =   isltyp_f(itime_f)

   ! if missing, replace with *_ref values
   PRINT *,'WARNING: Using defaults for landuse, vegetation,&
        &and soil type from _ref namelist - reset to the actual values'
   PRINT *,'WARNING: Using USGS table'
   
   IF ( MAXVAL(lu_index) < 0.0 ) lu_index = lu_index_ref
   IF ( MAXVAL(vegfra) < 0.0 )   vegfra = vegfra_ref
   IF ( MAXVAL(ivgtyp) < 0 )     ivgtyp = ivgtyp_ref
   IF ( MAXVAL(isltyp) < 0 )     isltyp = isltyp_ref
   
   IF (mminlu=='USGS') THEN
      iswater=16
   ELSE IF (mminlu(1:3)=='OLD') THEN
      iswater=7
   ELSE IF (mminlu(1:3)=='SiB') THEN
      iswater=15
   ELSE IF (mminlu=='LW12') THEN
      iswater=2
   ELSE
      PRINT *,'Unknown Landuse Table ',mminlu
      PRINT *,'Stopping'
      STOP
   ENDIF

   CALL landuse_init(lu_index, snowc, Albedo, Albbck, Mavail, &
        Emiss,  &
        Znt,Z0,Thc,Cs, Xland, julday, cent_lat, iswater, mminlu, &
        ids, ide, jds, jde, kds, kde,                       &
        ims, ime, jms, jme, kms, kme,                       &
        its, ite, jts, jte, kts, kte                       )

   CALL bucket_init(mminlu,lu_index,Maxm,Minm,Erate,&
        &ids, ide, jds, jde, kds, kde,   &
        ims, ime, jms, jme, kms, kme,   &
        its, ite, jts, jte, kts, kte)


   seamask=xland-1.
   landmask_input=1.-seamask
   tsk=ts_init_f(itime_f)
   qsfc=qvs_init_f(itime_f)


   IF (tsk(1,1) < 0.) THEN
      tsk=ts_ref
      PRINT *,'WARNING: No soil temperature found. Using reference',&
           &'for skin temperature'
   ENDIF
   
   tmn=tsk
   tslb=tsk(1,1)

   smois=qsfc(1,1)/qsat(ps_ref,tsk(1,1))
   
   sst_input=tsk

   CALL init_soil_depth_0 ( zs , dzs , num_soil_layers,thc,cs)

   
 END SUBROUTINE init_soil_real_fluxforce


END MODULE module_init_soil_real_fluxforce


MODULE module_init_soil_ideal

! Module takes care of driving the soil initialization for the ideal
! case.

USE module_getsm,           only: getsmlsm, getsmruclsm

USE module_ideal,           ONLY: mysoil,qsat

USE module_luse_init,       ONLY: landuse_init,bucket_init

USE module_model_constants, only: DEGRAD, eomeg

USE module_namelist,        only: lat_ref, julday_ref, mminlu_ref, &
                                  lu_index_ref, ivgtyp_ref, isltyp_ref, &
                                  vegfra_ref, rland_ref, albedo_ref, &
                                  zo_ref, thc_ref, cs_ref, emiss_ref, ts_ref, &
                                  ps_ref, dtamplitude_ref, mavail_ref, &
                                  sf_surface_physics, tmn_ref, deep_soil_moisture,&
                                  &Maxm_ref,Minm_ref,Erate_ref, bucket_model

USE module_soil_pre,        only: init_soil_depth_1, init_soil_depth_2, &
                                  init_soil_depth_3, init_soil_1_real, &
                                  init_soil_2_real, init_soil_3_real,&
                                  &init_soil_depth_0

USE module_sf_noahlsm,      only: lsminit

USE module_sf_ruclsm,       only: lsmrucinit


IMPLICIT NONE

private

public::   init_soil_ideal

CONTAINS

SUBROUTINE init_soil_ideal(julday, lu_index, ivgtyp, &
                           isltyp, vegfra, iswater, &
                           cent_lat, cor, mminlu, rland_ref, &
                           snowc, albedo, albbck,  mavail, emiss, &
                           maxm,minm,erate,&
                           znt, z0, thc, cs, xland, xice, seamask, &
                           snow, snowh, canwat, smstav, smstot, &
                           sfcrunoff, udrunoff, acsnow, acsnom, &
                           fndsoilw, fndsnowh, smfr3d,  &
                           landmask_input, time, tsoil, qsoil, &
                           tsk, qsfc, tmn, tslb, smois, sh2o, &
                           sst_input, flag_sst, &
                           num_soil_layers,  num_st_levels_input , &
                           num_sm_levels_input, &
                           st_levels_input, sm_levels_input, &
                           zs, dzs, st_input, sm_input, keepfr3dflag, &
                           smcmin, smcmax, restart, &
                           ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte )

   INTEGER, INTENT(in)         :: ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  num_soil_layers,  num_st_levels_input, &
                                  num_sm_levels_input
   INTEGER, INTENT(out)        :: julday, iswater
   INTEGER, DIMENSION(:,:), INTENT(out) ::  ivgtyp, isltyp
   INTEGER, DIMENSION(:), INTENT(out) :: st_levels_input, sm_levels_input
   REAL,    INTENT( in)        :: time
   REAL,    INTENT(out)        :: cent_lat, cor, rland_ref
   REAL,    INTENT(out)        :: tsoil, qsoil, smcmin, smcmax
   REAL, DIMENSION(:,:), INTENT(out) :: lu_index, vegfra, snowc, &
                                        albedo, albbck, mavail, emiss, &
                                        &maxm,minm,erate,&
                                        znt, z0, thc, cs, xland, seamask, &
                                        landmask_input, tsk, qsfc, tmn, &
                                        sst_input, snow, snowh, canwat, &
                                        smstav, smstot, sfcrunoff, udrunoff, &
                                        acsnow, acsnom, xice 
   REAL, DIMENSION(:,:,:), INTENT(out) :: tslb, smois, sh2o, st_input, &
                                          sm_input, smfr3d, keepfr3dflag
   REAL, DIMENSION(:),     INTENT(out) :: zs, dzs
   CHARACTER(len=4), INTENT(out) :: mminlu
   LOGICAL,  INTENT(in)         :: flag_sst, fndsoilw, fndsnowh, restart

   cent_lat=lat_ref*DEGRAD
   cor=2.*eomeg*SIN(cent_lat)

   julday =   julday_ref
   mminlu =   mminlu_ref
   lu_index=  lu_index_ref
   ivgtyp =   ivgtyp_ref
   isltyp =   isltyp_ref
   vegfra =   vegfra_ref

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

   IF (rland_ref < .5) rland_ref=2.

   IF (rland_ref > 1.5) THEN
      lu_index=16.
      CALL landuse_init(lu_index, snowc, Albedo, Albbck, Mavail, &
           &Emiss, Znt,Z0,Thc,Cs, Xland, julday, cent_lat, iswater, &
           &mminlu, &
            ids, ide, jds, jde, kds, kde,                       &
            ims, ime, jms, jme, kms, kme,                       &
            its, ite, jts, jte, kts, kte                       )
      CALL bucket_init(mminlu,lu_index,Maxm,Minm,Erate,&
           &ids, ide, jds, jde, kds, kde,   &
           ims, ime, jms, jme, kms, kme,   &
           its, ite, jts, jte, kts, kte)
      
  
   ELSE
      xland=rland_ref
      z0=zo_ref
      znt=z0
      thc=thc_ref
      cs=cs_ref
      emiss=emiss_ref
      albedo=albedo_ref
      albbck=albedo_ref
      maxm=maxm_ref
      minm=minm_ref
      erate=erate_ref
   ENDIF

   seamask=xland-1.
   landmask_input=1.-seamask

   CALL mysoil(rland_ref,ts_ref,ps_ref,dtamplitude_ref,&
           &mavail_ref,time,tsoil,qsoil)

   tsk=tsoil
   qsfc=qsoil
   tmn=tmn_ref
   sst_input=tsk
 
   IF (sf_surface_physics=='SIMPLESCHEME') THEN 
      
      tslb=tsoil
      smois=0.
      sh2o = 0.

   ELSEIF (sf_surface_physics=='FLUXFORCE') THEN

      tslb=tsoil
      smois=qsfc(1,1)/qsat(ps_ref,tsk(1,1))
      
      CALL init_soil_depth_0 ( zs , dzs , num_soil_layers,thc,cs)
      

   ELSEIF (sf_surface_physics=='FRSCHEME') THEN

      CALL init_soil_depth_0 ( zs, dzs, num_soil_layers, thc, cs )

      tslb=tsoil

      IF (bucket_model=='BUCKET2') THEN
! need fixing still
      ELSE
         smois=0.
         sh2o = 0.
      ENDIF
      
   ELSEIF (sf_surface_physics=='SLABSCHEME') THEN

      CALL init_soil_depth_1 ( zs , dzs , num_soil_layers )
      CALL init_soil_1_real ( tsk , tmn , tslb , zs , dzs , &
           num_soil_layers , 0 , &
           landmask_input , sst_input , flag_sst , &
           ids , ide , jds , jde , kds , kde , &
           ims , ime , jms , jme , kms , kme , &
           its , ite , jts , jte , kts , kte )
       smois=0.
       sh2o = 0.

   ELSEIF (sf_surface_physics=='LSMSCHEME') THEN

      CALL init_soil_depth_2 ( zs , dzs , num_soil_layers )

      st_levels_input(1)=5
      st_levels_input(2)=250
      sm_levels_input(1)=5
      sm_levels_input(2)=250

      st_input(:,:,2)=tsk
      st_input(:,:,3)=tmn

      CALL getsmlsm (mminlu, isltyp(1,1), ivgtyp(1,1), &
           Smcmin, Smcmax)

      sm_input(:,:,2)=MAX(MIN(mavail_ref,smcmax),smcmin)
      sm_input(:,:,3)=MAX(deep_soil_moisture*smcmax,smcmin)
      CALL init_soil_2_real ( tsk , tmn , smois , tslb , &
              st_input , sm_input , landmask_input , sst_input , &
              zs , dzs , &
              st_levels_input , sm_levels_input , &
              num_soil_layers , num_st_levels_input , &
              &num_sm_levels_input , &
              flag_sst , &
              ids , ide , jds , jde , kds , kde , &
              ims , ime , jms , jme , kms , kme , &
              its , ite , jts , jte , kts , kte )

      CALL LSMINIT(VEGFRA,SNOW,SNOWC,SNOWH,CANWAT,SMSTAV,  &
           SMSTOT, SFCRUNOFF,UDRUNOFF,ACSNOW,        &
           ACSNOM,IVGTYP,ISLTYP,TSLB,SMOIS,SH2O,ZS,DZS, &
           FNDSOILW, FNDSNOWH,                       &
           num_soil_layers, .FALSE.,                 &
           .TRUE., &
           ids,ide, jds,jde, kds,kde,                &
           ims,ime, jms,jme, kms,kme,                &
           its,ite, jts,jte, kts,kte                 )

      ELSEIF (sf_surface_physics=='RUCLSMSCHEME') THEN

         CALL init_soil_depth_3 ( zs , dzs , num_soil_layers )

         st_levels_input(1)=0
         st_levels_input(2)=300
         sm_levels_input(1)=0
         sm_levels_input(2)=300
         st_input(:,:,1)=tsk
         st_input(:,:,2)=tmn

         CALL getsmruclsm (mminlu,isltyp(1,1), ivgtyp(1,1), &
              &Smcmin, Smcmax)

         sm_input(:,:,1)=MAX(MIN(mavail_ref,smcmax),smcmin)-smcmin
         sm_input(:,:,2)=MAX(deep_soil_moisture*smcmax,smcmin)-smcmin

        CALL init_soil_3_real ( tsk , tmn , smois , tslb , &
              st_input , sm_input , landmask_input , sst_input , &
              zs , dzs , &
              st_levels_input , sm_levels_input , &
              num_soil_layers , num_st_levels_input , &
              &num_sm_levels_input , &
              flag_sst , &
              ids , ide , jds , jde , kds , kde , &
              ims , ime , jms , jme , kms , kme , &
              its , ite , jts , jte , kts , kte )

        mavail=mavail_ref

         CALL LSMRUCINIT(SMFR3D,TSLB,SMOIS,ISLTYP,mavail,&
              num_soil_layers, .FALSE.,.TRUE.,                 &
              ids,ide, jds,jde, kds,kde,                &
              ims,ime, jms,jme, kms,kme,                &
              its,ite, jts,jte, kts,kte                 )

      ENDIF

      xice = 0.
      snow = 0.
      snowc = 0.
      canwat = 0.
      KEEPFR3DFLAG=0.
      SMFR3D=0.
      sh2o = 0.

END SUBROUTINE init_soil_ideal

END MODULE module_init_soil_ideal


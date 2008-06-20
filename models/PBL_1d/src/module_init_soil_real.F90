MODULE module_init_soil_real

! Module takes care of driving the soil initialization for the ideal
! case.

USE module_getsm,           only: getsmlsm, getsmruclsm

USE module_ideal,           ONLY: mysoil,qsat

USE module_luse_init,       ONLY: landuse_init,bucket_init

USE module_model_constants, only: DEGRAD, eomeg

USE module_namelist,        ONLY: lat_ref,mminlu_ref, &
                                  lu_index_ref, ivgtyp_ref, isltyp_ref, &
                                  vegfra_ref, rland_ref, albedo_ref, &
                                  zo_ref, thc_ref, cs_ref,emiss_ref, ts_ref, &
                                  ps_ref, dtamplitude_ref, mavail_ref, &
                                  sf_surface_physics, tmn_ref, &
                                  deep_soil_moisture,&
                                  &Maxm_ref,Minm_ref,Erate_ref, bucket_model

USE module_soil_pre,        only: init_soil_depth_1, init_soil_depth_2, &
                                  init_soil_depth_3, init_soil_1_real, &
                                  init_soil_2_real, init_soil_3_real,&
                                  &init_soil_depth_0,init_soil_0_real

USE module_sf_noahlsm,      only: lsminit, maxsmc, wltsmc

USE module_sf_ruclsm,       only: lsmrucinit


IMPLICIT NONE

private

PUBLIC::   init_soil_real_wrf,init_soil_real_snd

CONTAINS

SUBROUTINE init_soil_real_wrf(julday, lu_index, ivgtyp, &
                           isltyp, vegfra, iswater, &
                           cent_lat, cor, mminlu, &
                           snowc, albedo, albbck,  mavail, emiss, &
                           maxm,minm,erate,&
                           znt, z0, thc, cs,xland, xice, seamask, &
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
                           ns_f, zs_f, dzs_f, &
                           itime_f, vegfra_f, ivgtyp_f, isltyp_f, &
                           lu_index_f, &
                           tsk_init_f, tmn_init_f, qsfc_init_f, &
                           smois_init_f, tslb_init_f, &
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
   REAL, DIMENSION(:,:), INTENT( in) :: smois_init_f, tslb_init_f

   REAL,    INTENT(out)        :: cent_lat, cor
   REAL,    INTENT(out)        :: smcmin, smcmax
   REAL, DIMENSION(:), INTENT(out) :: zs_f, dzs_f
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
   REAL, DIMENSION(ims:ime,jms:jme) :: smois_frac

! local
   INTEGER                      :: k 

   cent_lat=lat_ref*DEGRAD
   cor=2.*eomeg*SIN(cent_lat)

   lu_index=  REAL(lu_index_f(itime_f))
   vegfra =   vegfra_f(itime_f)
   ivgtyp =   ivgtyp_f(itime_f)
   isltyp =   isltyp_f(itime_f)

   ! if missing, replace with *_ref values
   IF ( maxval(lu_index) < 0.0 ) lu_index = lu_index_ref
   IF ( maxval(vegfra) < 0.0 )   vegfra = vegfra_ref
   IF ( maxval(ivgtyp) < 0 )     ivgtyp = ivgtyp_ref
   IF ( maxval(isltyp) < 0 )     isltyp = isltyp_ref

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
        Znt,Z0,Thc,cs,Xland, julday, cent_lat, iswater, mminlu, &
        ids, ide, jds, jde, kds, kde,                       &
        ims, ime, jms, jme, kms, kme,                       &
        its, ite, jts, jte, kts, kte                       )

   CALL bucket_init(mminlu,lu_index,Maxm,Minm,Erate,&
        &ids, ide, jds, jde, kds, kde,   &
        ims, ime, jms, jme, kms, kme,   &
        its, ite, jts, jte, kts, kte)
   
   seamask=xland-1.
   landmask_input=1.-seamask


   tsk=tsk_init_f(itime_f)
   tmn=tmn_init_f(itime_f)
   qsfc=qsfc_init_f(itime_f)
   sst_input=tsk_init_f(itime_f)


   SELECT CASE (ns_f)
      CASE (4)
         CALL init_soil_depth_2 ( zs_f , dzs_f , ns_f )
      CASE (5)
         CALL init_soil_depth_1 ( zs_f , dzs_f , ns_f )
      CASE DEFAULT
         print*,'Do not know how to initialize ',ns_f,' INPUT soil layers'
         stop 'module_init_soil_real'
   END SELECT

   IF (sf_surface_physics=='SIMPLESCHEME') THEN
      PRINT *,'SIMPLESCHEME cannot be run with init_f = .TRUE.'
      PRINT *, 'Stopping'
      STOP

   ELSEIF (sf_surface_physics=='FRSCHEME') THEN
      
      CALL init_soil_depth_0 ( zs , dzs , num_soil_layers,thc,cs)
      
      tslb=tsk(1,1)
      
      CALL init_soil_0_real ( tsk , tmn , tslb_init_f(:,itime_f),&
           &zs , dzs , zs_f, &
           num_soil_layers , ns_f, 1 , &
           landmask_input , sst_input , flag_sst , &
           ids , ide , jds , jde , kds , kde , &
           ims , ime , jms , jme , kms , kme , &
           its , ite , jts , jte , kts , kte )

      IF (bucket_model=='BUCKET2') THEN
      ELSE
         smois=0.
         sh2o = 0.
      ENDIF

      
   ELSEIF (sf_surface_physics=='SLABSCHEME') THEN

      CALL init_soil_depth_1 ( zs , dzs , num_soil_layers )

      SELECT CASE (ns_f)
         CASE (4)
            CALL init_soil_1_real ( tsk , tmn , tslb , zs , dzs , &
                 num_soil_layers , 1 , &
                 landmask_input , sst_input , flag_sst , &
                 ids , ide , jds , jde , kds , kde , &
                 ims , ime , jms , jme , kms , kme , &
                 its , ite , jts , jte , kts , kte )
            smois=0.

         CASE (5)
            DO k=1,num_soil_layers
               tslb(:,k,:)=tslb_init_f(k,itime_f)
            ENDDO
            smois=0.

         CASE DEFAULT
            print*,'Do not know how to initialize ',ns_f,' INPUT soil layers'
            print*,'for SLABSCHEME'
            stop 'module_init_soil_real'

      END SELECT

      IF  (qsfc(1,1) <= 0.) THEN
         qsfc(1,1)=qsat(ps_ref,tsk(1,1))*mavail(1,1)
      ENDIF

      keepfr3dflag=0.
      smfr3d=0.
      xice = 0.
      snow = 0.
      snowc = 0.
      canwat = 0.
      sh2o = 0.

   ELSEIF (sf_surface_physics=='LSMSCHEME') THEN

      CALL init_soil_depth_2 ( zs , dzs , num_soil_layers )

      SELECT CASE (ns_f)
         CASE (4)
            DO k=1,num_soil_layers
               tslb(:,k,:)=tslb_init_f(k,itime_f)
               smois(:,k,:)=smois_init_f(k,itime_f)
            ENDDO

         CASE (5)
            st_levels_input(1)=NINT(zs_f(ns_f-2)*100.)
            st_levels_input(2)=NINT(zs_f(ns_f)*100.)
            sm_levels_input(1)=NINT(zs_f(ns_f-2)*100.)
            sm_levels_input(2)=NINT(zs_f(ns_f)*100.)
   
            st_input(1,1,2)=tslb_init_f(ns_f-2,itime_f)
            st_input(1,1,3)=tslb_init_f(ns_f,itime_f)

            CALL getsmlsm (mminlu, isltyp(1,1), ivgtyp(1,1), &
                 Smcmin, Smcmax)

            sm_input(:,:,2)=MAX(MIN(mavail,smcmax),smcmin)
            sm_input(:,:,3)=MAX(deep_soil_moisture*smcmax,smcmin)

            CALL init_soil_2_real ( tsk , tmn , smois , tslb , &
                 st_input , sm_input , landmask_input , sst_input , &
                 zs , dzs , &
                 st_levels_input , sm_levels_input , &
                 num_soil_layers , num_st_levels_input , &
                 num_sm_levels_input , &
                 flag_sst , &
                 ids , ide , jds , jde , kds , kde , &
                 ims , ime , jms , jme , kms , kme , &
                 its , ite , jts , jte , kts , kte )

         CASE DEFAULT
            print*,'Do not know how to initialize ',ns_f,' INPUT soil layers'
            print*,'for LSMSCHEME'
            stop 'module_init_soil_real'

      END SELECT

      keepfr3dflag=0.
      smfr3d=0.
      xice = 0.
      snow = 0.
      snowc = 0.
      canwat = 0.
      sh2o = 0.

      CALL LSMINIT(VEGFRA,SNOW,SNOWC,SNOWH,CANWAT,SMSTAV,  &
           SMSTOT, SFCRUNOFF,UDRUNOFF,ACSNOW,        &
           ACSNOM,IVGTYP,ISLTYP,TSLB,SMOIS,SH2O,ZS,DZS, &
           FNDSOILW, FNDSNOWH,                       &
           num_soil_layers,                  &
           .FALSE.,.TRUE., &
           ids,ide, jds,jde, kds,kde,                &
           ims,ime, jms,jme, kms,kme,                &
           its,ite, jts,jte, kts,kte                 )

      smois_frac = (smois(:,1,:)-wltsmc(isltyp(1,1)))/(maxsmc(isltyp(1,1))-wltsmc(isltyp(1,1)))
      IF  (qsfc(1,1) <= 0.) THEN
         qsfc(1,1)=qsat(ps_ref,tsk(1,1))*smois_frac(1,1) !smois(1,1,1)
!         qsfc(1,1)=qsat(ps_ref,tsk(1,1))*smois(1,1,1)
      ENDIF
      
   ELSEIF (sf_surface_physics=='RUCLSMSCHEME') THEN

      CALL init_soil_depth_3 ( zs , dzs , num_soil_layers )

      SELECT CASE (ns_f)
         CASE (4)
            st_levels_input(1)=0
            st_levels_input(2)=300
            st_input(:,:,1)=tsk
            st_input(:,:,2)=tmn

            CALL getsmruclsm (mminlu,isltyp(1,1), ivgtyp(1,1), &
                 Smcmin, Smcmax)

            sm_levels_input(1)=0.
            sm_levels_input(2)=300.

            sm_input(:,:,1)= &
               MAX(MIN(smois_init_f(1,itime_f),smcmax),smcmin)-smcmin
            sm_input(:,:,2)= &
               MAX(MIN(smois_init_f(ns_f,itime_f),smcmax),smcmin)-smcmin

         CASE (5)
            st_levels_input(1)=0
            st_levels_input(2)=300
            sm_levels_input(1)=0
            sm_levels_input(2)=300
            st_input(:,:,1)=tsk
            st_input(:,:,2)=tmn

            CALL getsmruclsm (mminlu,isltyp(1,1), ivgtyp(1,1), &
                 Smcmin, Smcmax)

             sm_input(:,:,1)=MAX(MIN(mavail,smcmax),smcmin)-smcmin
             sm_input(:,:,2)=MAX(deep_soil_moisture*smcmax,smcmin)-smcmin

         CASE DEFAULT
            print*,'Do not know how to initialize ',ns_f,' INPUT soil layers'
            print*,'for RUCLSMSCHEME'
            stop 'module_init_soil_real'

      END SELECT

      CALL init_soil_3_real ( tsk , tmn , smois , tslb , &
           st_input , sm_input , landmask_input , sst_input , &
           zs , dzs , &
           st_levels_input , sm_levels_input , &
           num_soil_layers , num_st_levels_input , &
           num_sm_levels_input , &
           flag_sst , &
           ids , ide , jds , jde , kds , kde , &
           ims , ime , jms , jme , kms , kme , &
           its , ite , jts , jte , kts , kte )

      IF ( ns_f == 4 ) THEN
         smois(:,1,:)=(smois(:,2,:)-smois(:,3,:)*zs_f(2)/zs_f(3))/&
              (1.-zs_f(2)/zs_f(3))
      ENDIF

      keepfr3dflag=0.
      smfr3d=0.
      xice = 0.
      snow = 0.
      snowc = 0.
      canwat = 0.
      sh2o = 0.

      CALL LSMRUCINIT(SMFR3D,TSLB,SMOIS,ISLTYP,mavail,&
           num_soil_layers, .FALSE.,.TRUE., &
           ids,ide, jds,jde, kds,kde,                &
           ims,ime, jms,jme, kms,kme,                &
           its,ite, jts,jte, kts,kte                 )

      IF  (qsfc(1,1) <= 0.) THEN
         qsfc(1,1)=qsat(ps_ref,tsk(1,1))*(smois(1,1,1)+smcmin)
      ENDIF

   ENDIF

 END SUBROUTINE init_soil_real_wrf

SUBROUTINE init_soil_real_snd(julday, lu_index, ivgtyp, &
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
   tsk=tsk_init_f(itime_f)

   IF (tsk(1,1) < 0.) THEN
      IF (tslb_init_f(1,itime_f) > 0. .AND. &
           &tslb_init_f(2,itime_f) > 0.) THEN
         tsk=tslb_init_f(1,itime_f)-&
              &(tslb_init_f(2,itime_f)-tslb_init_f(1,itime_f))/&
              &(zs_f(2)-zs_f(1))*zs_f(1)
      ELSEIF (tslb_init_f(1,itime_f) > 0. .AND. &
           &tslb_init_f(2,itime_f) < 0.) THEN
         tsk=tslb_init_f(1,itime_f)
      ELSE
         tsk=ts_ref
         PRINT *,'WARNING: No soil temperature found. Using reference',&
              &'for skin temperature'
      ENDIF
   ENDIF

   tmn=tmn_init_f(itime_f)

   IF (tmn(1,1) < 0.) THEN
      DO k=ns_f,1,-1
         IF (tslb_init_f(k,itime_f) > 0.) THEN
            tmn=tslb_init_f(k,itime_f)
            kns_f_t=k
            EXIT
         ENDIF
      ENDDO

      IF (tmn(1,1) < 0.) THEN
         PRINT *,'WARNING: No soil temperature found using reference ',&
              &'for deep soil temperature'
         tmn = tmn_ref
      ENDIF
   ENDIF

   DO k=ns_f,1,-1
      IF (smois_init_f(k,itime_f) > 0.) THEN
         smois=smois_init_f(k,itime_f)
         kns_f_m=k
         EXIT
      ENDIF
   ENDDO
   
   sst_input=tsk_init_f(itime_f)

   IF (sf_surface_physics=='SIMPLESCHEME') THEN
      PRINT *,'SIMPLESCHEME cannot be run with init_f = .TRUE.'
      PRINT *, 'Stopping'
      STOP

  ELSEIF (sf_surface_physics=='FRSCHEME') THEN
      
      CALL init_soil_depth_0 ( zs , dzs , num_soil_layers,thc,cs)
      
      tslb=tsk(1,1)
      
      CALL init_soil_0_real ( tsk , tmn , tslb_init_f(:,itime_f),&
           &zs , dzs , zs_f, &
           num_soil_layers , ns_f, 1 , &

           landmask_input , sst_input , flag_sst , &
           ids , ide , jds , jde , kds , kde , &
           ims , ime , jms , jme , kms , kme , &
           its , ite , jts , jte , kts , kte )

      IF (bucket_model=='BUCKET2') THEN
      ELSE
         smois=0.
         sh2o = 0.
      ENDIF

   ELSEIF (sf_surface_physics=='SLABSCHEME') THEN
      
      CALL init_soil_depth_1 ( zs , dzs , num_soil_layers )

      IF (tsk(1,1) < 0.) THEN
         PRINT *,'WARNING: No soil temperature found using reference ',&
              &'for skin and deep soil temperature'
      ELSE         
         CALL init_soil_1_real ( tsk , tmn , tslb , zs , dzs , &
              num_soil_layers , 1 , &
           landmask_input , sst_input , flag_sst , &
           ids , ide , jds , jde , kds , kde , &
           ims , ime , jms , jme , kms , kme , &
           its , ite , jts , jte , kts , kte )
      ENDIF

      smois=0.
      keepfr3dflag=0.
      smfr3d=0.
      xice = 0.
      snow = 0.
      snowc = 0.
      canwat = 0.
      sh2o = 0.

      IF  (qsfc(1,1) <= 0.) THEN
         qsfc(1,1)=qsat(ps_ref,tsk(1,1))*mavail(1,1)
      ENDIF

   ELSEIF (sf_surface_physics=='LSMSCHEME') THEN

      CALL init_soil_depth_2 ( zs , dzs , num_soil_layers )
print*,tslb_init_f(:,itime_f)

      IF (tslb_init_f(1,itime_f) < 0.) THEN
         PRINT *,'WARNING: No soil temperature found using reference',&
              &'for skin and deep soil temperature'
         st_levels_input(1)=NINT(zs_f(1))
         st_levels_input(2)=NINT(zs_f(kns_f_t))
         st_input(1,1,2)=ts_ref
         st_input(1,1,3)=tmn_ref
      ELSE
         st_levels_input(1)=NINT(zs_f(1))
         st_levels_input(2)=NINT(zs_f(kns_f_t))
         st_input(1,1,2)=tslb_init_f(1,itime_f)
         st_input(1,1,3)=tslb_init_f(kns_f_t,itime_f)
      ENDIF

      sm_levels_input(1)=NINT(zs_f(1))
      sm_levels_input(2)=NINT(zs_f(kns_f_m))

      CALL getsmlsm (mminlu, isltyp(1,1), ivgtyp(1,1), &
           Smcmin, Smcmax)

      IF (smois_init_f(kns_f_m,itime_f) < 0.)  THEN
         PRINT *,'WARNING: No soil moisture found using reference',&
              &'for surface moisture availability and &
              &deep soil moisture'
         sm_input(:,:,2)=MAX(MIN(mavail,smcmax),smcmin)
         sm_input(:,:,3)=MAX(deep_soil_moisture*smcmax,smcmin)
      ELSE
         sm_input(1,1,2)=smois_init_f(1,itime_f)
         IF (sm_input(1,1,2) < 0.) THEN
            PRINT *,'WARNING: No surface soil moisture found',&
                 &'Using surface moisture availability'
            sm_input(:,:,2)=MAX(MIN(mavail,smcmax),smcmin)
         ENDIF
         sm_input(1,1,3)=smois_init_f(kns_f_m,itime_f)
      ENDIF

      CALL init_soil_2_real ( tsk , tmn , smois , tslb , &
           st_input , sm_input , landmask_input , sst_input , &
           zs , dzs , &
           st_levels_input , sm_levels_input , &
           num_soil_layers , num_st_levels_input , &
           num_sm_levels_input , &
           flag_sst , &
           ids , ide , jds , jde , kds , kde , &
           ims , ime , jms , jme , kms , kme , &
           its , ite , jts , jte , kts , kte )

      keepfr3dflag=0.
      smfr3d=0.
      xice = 0.
      snow = 0.
      snowc = 0.
      canwat = 0.
      sh2o = 0.
      
      CALL LSMINIT(VEGFRA,SNOW,SNOWC,SNOWH,CANWAT,SMSTAV,  &
           SMSTOT, SFCRUNOFF,UDRUNOFF,ACSNOW,        &
           ACSNOM,IVGTYP,ISLTYP,TSLB,SMOIS,SH2O,ZS,DZS, &
           FNDSOILW, FNDSNOWH,                       &
           num_soil_layers, .FALSE.,                 &
           .TRUE., &
           ids,ide, jds,jde, kds,kde,                &
           ims,ime, jms,jme, kms,kme,                &
           its,ite, jts,jte, kts,kte                 )

      IF  (qsfc(1,1) <= 0.) THEN
         qsfc(1,1)=qsat(ps_ref,tsk(1,1))*smois(1,1,1)
      ENDIF

   ELSEIF (sf_surface_physics=='RUCLSMSCHEME') THEN

      CALL init_soil_depth_3 ( zs , dzs , num_soil_layers )

      IF (tslb_init_f(1,itime_f) < 0.) THEN
         PRINT *,'WARNING: No soil temperature found using reference',&
              &'for skin and deep soil temperature'
         st_input(1,1,1)=ts_ref
         st_input(1,1,2)=tmn_ref
      ELSE
         st_input(:,:,1)=tsk
         st_input(:,:,2)=tmn
      ENDIF

      st_levels_input(1)=0.
      st_levels_input(2)=300.
      

      sm_levels_input(1)=zs_f(1)
      sm_levels_input(2)=300.

      CALL getsmruclsm (mminlu,isltyp(1,1), ivgtyp(1,1), &
           Smcmin, Smcmax)
      
      IF (smois_init_f(kns_f_m,itime_f) < 0.)  THEN
         PRINT *,'WARNING: No soil moisture found using reference',&
              &'for surface moisture availability and &
              &deep soil moisture'
         sm_input(:,:,1)=MAX(MIN(mavail,smcmax),smcmin)
         sm_input(:,:,2)=MAX(deep_soil_moisture*smcmax,smcmin)
      ELSE
         sm_input(1,1,1)=smois_init_f(1,itime_f)
         IF (sm_input(1,1,1) < 0.) THEN
            PRINT *,'WARNING: No surface soil moisture found',&
                 &'Using surface moisture availability'
            sm_input(:,:,1)=MAX(MIN(mavail,smcmax),smcmin)
         ENDIF
         sm_input(1,1,2)=smois_init_f(kns_f_m,itime_f)
      ENDIF
      
      sm_input(:,:,1)= &
           MAX(MIN(sm_input(:,:,1),smcmax),smcmin)-smcmin
      sm_input(:,:,2)= &
           MAX(MIN(sm_input(:,:,2),smcmax),smcmin)-smcmin
      
      
      CALL init_soil_3_real ( tsk , tmn , smois , tslb , &
           st_input , sm_input , landmask_input , sst_input , &
           zs , dzs , &
           st_levels_input , sm_levels_input , &
           num_soil_layers , num_st_levels_input , &
           num_sm_levels_input , &
           flag_sst , &
           ids , ide , jds , jde , kds , kde , &
           ims , ime , jms , jme , kms , kme , &
           its , ite , jts , jte , kts , kte )

      smois(:,1,:)=(smois(:,2,:)-smois(:,3,:)*zs_f(2)/zs_f(3))/&
           (1.-zs_f(2)/zs_f(3))
      
      keepfr3dflag=0.
      smfr3d=0.
      xice = 0.
      snow = 0.
      snowc = 0.
      canwat = 0.
      sh2o = 0.

      CALL LSMRUCINIT(SMFR3D,TSLB,SMOIS,ISLTYP,mavail,&
           num_soil_layers,  .FALSE.,.TRUE.,                 &
           ids,ide, jds,jde, kds,kde,                &
           ims,ime, jms,jme, kms,kme,                &
           its,ite, jts,jte, kts,kte                 )


      IF  (qsfc(1,1) <= 0.) THEN
         qsfc(1,1)=qsat(ps_ref,tsk(1,1))*(smois(1,1,1)+smcmin)
      ENDIF


   ENDIF

 END SUBROUTINE init_soil_real_snd


END MODULE module_init_soil_real


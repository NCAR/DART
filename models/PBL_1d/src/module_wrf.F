MODULE module_wrf

  USE module_model_constants
  USE module_namelist
  USE module_ideal
  USE module_initialize
  USE module_soil_pre
  USE module_surface_driver
  USE module_pbl_driver
  USE module_sf_myjsfc
  USE module_bl_myjpbl
  USE module_sf_sfclay
  USE module_bl_mrf
  USE module_bl_ysu
  USE module_sf_noahlsm
  USE module_sf_ruclsm
  USE module_jh_procedures
  USE module_luse_init
  USE module_getsm
  
  IMPLICIT NONE

  INTEGER :: isfcunit=55,iprofunit=56,isimil=57,ncunit=58

  INTEGER,PARAMETER :: ids=1,ide=1, jds=1,jde=1, kds=1,kde=nz, &
        &ims=1,ime=1, jms=1,jme=1, kms=1,kme=nz, &
        &its=1,ite=1,jts=1,jte=1,&
        &kts=1,kte=nz-1, num_tiles=1

  INTEGER, PARAMETER :: k_simil=5,PARAM_FIRST_SCALAR = 2
  INTEGER :: niter,iter

  REAL, DIMENSION(1:k_simil) :: u_simil,v_simil,q_simil,th_simil
  REAL :: hvfx,thvstar,qstar,dlmonin,z_level,zo,zol,x,psimzl,y,psihzl,&
       &rhosfc

  REAL, PARAMETER :: a_stable=1.,b_stable=2./3.,c_stable=5.,d_stable=.35

  INTEGER :: ntime,num_soil_layers,ninputtimes,nsplinetimes,ihour,&
       &imin,imax

  INTEGER, DIMENSION(num_tiles):: &
        &     i_start=1,i_end=1,j_start=1 ,j_end=1
   
   INTEGER       ::     itimestep,STEPBL,itime_gen=1
   INTEGER,    DIMENSION( ims:ime , jms:jme )                    &
                   ::                        LOWLYR
!
   LOGICAL          ::   warm_rain
!
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme )            &
            ::                         p_phy, &
        pi_phy, &
        p8w, &
        rho, &
        t_phy, &
        u_phy, &
        v_phy, &
        dz8w, &
        z, &
        th_phy

   REAL, DIMENSION( ims:ime, jms:jme ) ::   PSFC

   REAL,       DIMENSION( ims:ime, kms:kme+1, jms:jme ) :: z8w
!
   REAL, ALLOCATABLE, DIMENSION(:,:,:,: )         &
       ::                                    moist
!
!

   REAL,       DIMENSION( ims:ime , jms:jme )  ::        &
        XLAND, &
        HT, &
        PSIM, &
        PSIH, &
        GZ1OZ0, &
        BR, &
        CHKLOWQ

!mp@ begin
   REAL,       DIMENSION( ims:ime , jms:jme )  ::        &
        &psimfac,psihfac,psiqfac,zl2
!mp@ begin

   REAL,       DIMENSION( ims:ime, jms:jme )  &
       ::                           TSK, &
        UST, &
        HOL, &
        MOL, &
        PBLH, &
        HFX, &
        QFX, &
        REGIME, &
        ZNT, &
        QSFC, &
        AKHS, &
        AKMS, &
        QZ0, &
        THZ0, &
        UZ0, &
        VZ0, &
        GRDFLX  , &
        WSPD
   
!
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme )            &
            ::                       RUBLTEN, &
        RVBLTEN, &
        RTHBLTEN, &
        RQVBLTEN, &
        RQCBLTEN, &
        RQIBLTEN, &
        TKE_MYJ
   
   REAL                                ::  u_frame=0.,v_frame=0.

   
   INTEGER,    DIMENSION( ims:ime , jms:jme )                    &
         ::                             KPBL
 
   REAL,    DIMENSION( ims:ime , jms:jme ) :: landmask_input,sst_input,&
        &cpm,chs
   LOGICAL :: flag_sst=.TRUE.

   REAL,       DIMENSION( ims:ime , jms:jme )                    &
         :: XICE, SEAMASK, CT, SNOW, LH
   

   REAL    :: DTMIN,DTBL

   INTEGER :: i,J,K,NK,jj,ij,iitime

   REAL :: dx,time,cor,fract
   CHARACTER(len=4) :: mminlu
   INTEGER :: julday

   INTEGER, DIMENSION( ims:ime , jms:jme ) :: ISLTYP,IVGTYP
   INTEGER :: iswater   
   REAL :: smcmin,smcmax

   REAL, DIMENSION( ims:ime, jms:jme ) ::   &
        &ACSNOW,ACSNOM,ALBEDO,CANWAT,CAPG,EMISS,GLW,GSW,MAVAIL,Q10,Q2,&
        &T2,TH10,TH2,SNOWC,U10,V10,z0,PSHLTR,QSHLTR,RAINBL,&
        &RAINCV,RAINNCV,SFCEVP,SFCEXC,THC,TMN,SFCRUNOFF,SNOWH,POTEVP,&
        &SNOPCX,SOILTB,SOILT1,TSNAV,QSG,QVG,QCG,FLHC,FLQC,SNOALB,&
        &SMSTAV,SMSTOT,TSHLTR,UDRUNOFF,VEGFRA,ALBBCK,SHDMAX,SHDMIN
   
   REAL, DIMENSION( ims:ime, jms:jme ) ::  lu_index

   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: &
        &TSLB,SMOIS,KEEPFR3DFLAG,SMFR3D,SH2O

   REAL,  ALLOCATABLE, DIMENSION(:) :: DZS,ZS

   INTEGER :: num_st_levels_input,num_sm_levels_input
   
   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: &
        &st_input,sm_input
   INTEGER, ALLOCATABLE, DIMENSION(:) :: &
           &st_levels_input,sm_levels_input

   REAL, DIMENSION(1:nz) :: z_init,U_init,V_init,T_init,&
        &Th_init,Exn_init,Q_init,P_init,rho_init

   REAL :: Zo_init,tsoil,thsoil,qsoil

   REAL, DIMENSION(1:nz+1) :: z8w_init,p8w_init

   REAL, ALLOCATABLE, DIMENSION(:) :: u_mid,v_mid

   REAL :: cent_lat,cent_lon,truelat1,truelat2,&
        &sw_corner_lon,sw_corner_lat,lat,lon

   INTEGER :: projcode

   REAL,     DIMENSION(1:nz) :: uflux,vflux,hflux,qflux, k_t,k_m
   REAL,     DIMENSION(1:nz) :: u_g,v_g
   
   REAL :: a,b,c,d,e
   REAL, PARAMETER  :: deltat = 2.*3600

   LOGICAL, PARAMETER          :: FNDSOILW=.TRUE., FNDSNOWH=.TRUE.,&
        &restart=.FALSE.


   INTEGER :: idum , start_seconds, control_index1, control_index2
   REAL :: lenfor, control_w
   REAL, DIMENSION(ntimes_gen) :: times_gen
   REAL, DIMENSION(nz_gen,ntimes_gen) :: U_init_gen,V_init_gen,&
        &T_init_gen,&
        &Th_init_gen,Exn_init_gen,Q_init_gen,P_init_gen

! results of temporal and spatial interpolations
   REAL, ALLOCATABLE, DIMENSION(:) :: splinetimes
   REAL, ALLOCATABLE, DIMENSION(:,:) :: &
        &u_g_gen,v_g_gen,p_gen,p8w_gen,&
        &t_gen_uadv,t_gen_vadv,t_gen_wadv,&
        &q_gen_uadv,q_gen_vadv,q_gen_wadv,&
        &u_gen_uadv,u_gen_vadv,u_gen_wadv,&
        &v_gen_uadv,v_gen_vadv,v_gen_wadv
   

   REAL, ALLOCATABLE, DIMENSION(:) :: &
        glw_gen,gsw_gen

   REAL, DIMENSION(ntimes_gen) :: Th2_init_gen,Tsk_init_gen,&
        &U10_init_gen,V10_init_gen,     &
        &Q2_init_gen,Glw_init_gen,Gsw_init_gen,Qsfc_init_gen,&
        &Tmn_init_gen,vegfra_gen
   INTEGER, DIMENSION(ntimes_gen) :: isltyp_gen,ivgtyp_gen,lu_index_gen

   REAL, DIMENSION(ns_gen,ntimes_gen) :: Tslb_init_gen,Smois_init_gen
   REAL, DIMENSION(ns_gen) :: zs_gen,dzs_gen
   REAL, DIMENSION(nz_gen,ntimes_gen) :: z_gen
   REAL, DIMENSION(nz_gen+1,ntimes_gen) :: z_gen_stag
   
contains

SUBROUTINE STATIC_INIT_WRF(allocate_wrf)

  LOGICAL,INTENT(INOUT)           :: allocate_wrf

  IF (.NOT.init_gen) allocate_wrf=.TRUE.

  IF ( allocate_wrf ) THEN
    ALLOCATE(moist( ims:ime, kms:kme, jms:jme, n_moist))
    moist = 0.0
  ENDIF

  IF (surface_physics=='SIMPLESCHEME') THEN
     num_soil_layers=1
     IF ( allocate_wrf ) THEN         
        ALLOCATE(&
             &TSLB( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SMOIS( ims:ime , 1:num_soil_layers, jms:jme ),&
             &KEEPFR3DFLAG( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SMFR3D( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SH2O( ims:ime , 1:num_soil_layers, jms:jme ),&
             &DZS(1:num_soil_layers),&
             &ZS(1:num_soil_layers))
     ENDIF

  ELSEIF (surface_physics=='SLABSCHEME') THEN
     num_soil_layers=5
     IF ( allocate_wrf ) THEN
        ALLOCATE(&
             &TSLB( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SMOIS( ims:ime , 1:num_soil_layers, jms:jme ),&
             &KEEPFR3DFLAG( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SMFR3D( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SH2O( ims:ime , 1:num_soil_layers, jms:jme ),&
             &DZS(1:num_soil_layers),&
             &ZS(1:num_soil_layers))
     ENDIF
   ELSEIF (surface_physics=='LSMSCHEME') THEN
     num_soil_layers=4
     num_st_levels_input=2 
     num_sm_levels_input=2 
     IF ( allocate_wrf ) THEN
        ALLOCATE(&
             &st_levels_input(num_st_levels_input),&
             &sm_levels_input(num_sm_levels_input),&
             &st_input(ims:ime, jms:jme,num_st_levels_input+2),&
             &sm_input(ims:ime, jms:jme,num_sm_levels_input+2))
        ALLOCATE(&
             &TSLB( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SMOIS( ims:ime , 1:num_soil_layers, jms:jme ),&
             &KEEPFR3DFLAG( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SMFR3D( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SH2O( ims:ime , 1:num_soil_layers, jms:jme ),&
             &DZS(1:num_soil_layers),&
             &ZS(1:num_soil_layers))
     ENDIF

   ELSEIF (surface_physics=='RUCLSMSCHEME') THEN
     num_soil_layers=6
     num_st_levels_input=2
     num_sm_levels_input=2
     IF ( allocate_wrf ) THEN
        ALLOCATE(&
             &st_levels_input(num_st_levels_input),&
             &sm_levels_input(num_sm_levels_input),&
             &st_input(ims:ime, jms:jme,num_st_levels_input),&
             &sm_input(ims:ime, jms:jme,num_sm_levels_input))
        ALLOCATE(&
             &TSLB( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SMOIS( ims:ime , 1:num_soil_layers, jms:jme ),&
             &KEEPFR3DFLAG( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SMFR3D( ims:ime , 1:num_soil_layers, jms:jme ),&
             &SH2O( ims:ime , 1:num_soil_layers, jms:jme ),&
             &DZS(1:num_soil_layers),&
             &ZS(1:num_soil_layers))
     ENDIF

  ENDIF
  
! some other stuff we'll need
  nsplinetimes=INT((REAL(forecast_length)/REAL(splineinterval))+&
       &.999)+1
  
   IF (init_gen) THEN
      ntime=NINT(REAL(forecast_length)*3600./dt)
   ELSE
      ntime=NINT(timetot/dt)
   ENDIF

   IF ( allocate_wrf ) THEN
      ALLOCATE(splinetimes(nsplinetimes),&
           &u_g_gen(nz,nsplinetimes),v_g_gen(nz,nsplinetimes),&
           &p_gen(nz,nsplinetimes),p8w_gen(nz,nsplinetimes),&
           &t_gen_uadv(nz,nsplinetimes),&
           &t_gen_vadv(nz,nsplinetimes),&
           &q_gen_uadv(nz,nsplinetimes),&
           &q_gen_vadv(nz,nsplinetimes),&
           &u_gen_uadv(nz,nsplinetimes),&
           &u_gen_vadv(nz,nsplinetimes),&
           &v_gen_uadv(nz,nsplinetimes),&
           &v_gen_vadv(nz,nsplinetimes),&
           &u_gen_wadv(nz,nsplinetimes),&
           &v_gen_wadv(nz,nsplinetimes),&
           &t_gen_wadv(nz,nsplinetimes),&
           &q_gen_wadv(nz,nsplinetimes),&
           &glw_gen(nsplinetimes),gsw_gen(nsplinetimes))
           t_gen_vadv = 0.0
           u_gen_vadv = 0.0
           v_gen_vadv = 0.0
           q_gen_vadv = 0.0
           t_gen_uadv = 0.0
           u_gen_uadv = 0.0
           v_gen_uadv = 0.0
           q_gen_uadv = 0.0
           t_gen_wadv = 0.0
           u_gen_wadv = 0.0
           v_gen_wadv = 0.0
           q_gen_wadv = 0.0
   ENDIF

 
END SUBROUTINE STATIC_INIT_WRF


SUBROUTINE INIT_WRF(wrf_rnd_seed)

   INTEGER,INTENT(INOUT)        :: wrf_rnd_seed

!!! basic stuff

   a=1.-.25*dt**2*cor**2
   b=dt*cor
   c=.5*dt**2*cor**2
   d=.5*dt**2*cor
   e=1.+.25*dt**2*cor**2


   IF (.NOT.init_gen) THEN

!! allocate and assign soil, surface and atmosphere for ideal begin

      cent_lat=lat_ref*DEGRAD
      cor=2.*eomeg*SIN(cent_lat)

! initialize soil begin
      
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
         PRINT *,'Unknown Landuse Category',mminlu
         PRINT *,'Stopping'
         STOP
      ENDIF
      
      IF (rland_ref < .5) rland_ref=2.

      IF (rland_ref > 1.5) THEN
         lu_index=16.
         CALL landuse_init(lu_index, snowc, Albedo, Albbck, Mavail, &
              &Emiss, Znt,Z0,Thc,Xland, julday, cent_lat, iswater, &
              &mminlu, &
              ids, ide, jds, jde, kds, kde,                       &
              ims, ime, jms, jme, kms, kme,                       &
              its, ite, jts, jte, kts, kte                       )
      ELSE
         xland=rland_ref
         z0=zo_ref 
         znt=z0
         thc=thc_ref
         emiss=emiss_ref
         albedo=albedo_ref      
         albbck=albedo_ref
      ENDIF

      seamask=xland-1.
      landmask_input=1.-seamask

      time=timeo
      
      CALL mysoil(rland_ref,ts_ref,ps_ref,dtamplitude_ref,&
           &mavail_ref,time,Tsoil,Qsoil)
      
      tsk=tsoil
      qsfc=qsoil
      tmn=tmn_ref
      sst_input=tsk
      
      IF (surface_physics=='SIMPLESCHEME') THEN

         num_soil_layers=1
         tslb=tsoil
         smois=0.
         sh2o = 0.
         
      ELSEIF (surface_physics=='SLABSCHEME') THEN

         num_soil_layers=5
        CALL init_soil_depth_1 ( zs , dzs , num_soil_layers )
         CALL init_soil_1_real ( tsk , tmn , tslb , zs , dzs , &
           &num_soil_layers , 0 , &
           landmask_input , sst_input , flag_sst , &
           ids , ide , jds , jde , kds , kde , &
           ims , ime , jms , jme , kms , kme , &
           its , ite , jts , jte , kts , kte )
         smois=0.
         sh2o = 0.
         
      ELSEIF (surface_physics=='LSMSCHEME') THEN
         
         num_soil_layers=4
         
         CALL init_soil_depth_2 ( zs , dzs , num_soil_layers )
         
         st_levels_input(1)=5
         st_levels_input(2)=250
         sm_levels_input(1)=5
         sm_levels_input(2)=250
         
         st_input(:,:,2)=tsk
         st_input(:,:,3)=tmn

         CALL getsmlsm (mminlu, isltyp(1,1), ivgtyp(1,1), &
              &Smcmin, Smcmax)
         
         sm_input(:,:,2)=MAX(MIN(mavail_ref,smcmax),smcmin)
         sm_input(:,:,3)=&
              &MAX(deep_soil_moisture*smcmax,smcmin)

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
              num_soil_layers, restart,                 &
              ids,ide, jds,jde, kds,kde,                &
              ims,ime, jms,jme, kms,kme,                &
              its,ite, jts,jte, kts,kte                 )

      ELSEIF (surface_physics=='RUCLSMSCHEME') THEN
         
         num_soil_layers=6
         
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
         sm_input(:,:,2)=&
              &MAX(deep_soil_moisture*smcmax,smcmin)-smcmin


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
         
         CALL LSMRUCINIT(SMFR3D,TSLB,SMOIS,ISLTYP,&
              num_soil_layers, restart,                 &
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
      
! initialize soil end

! initialize atmosphere begin
      
      CALL initgrid(U_init,U_g,V_init,V_g,T_init,Th_init,&
           &Exn_init,Q_init,P_init,P8w_init,Rho_init,&
           &Z_init,Z8w_init,Zo_init,tsk(1,1))

      zl2=z_init(2)
      CALL initvar(z_init,z8w_init,u_init,v_init,t_init,th_init,&
           &exn_init,&
           &q_init,p_init,p8w_init,rho_init,tsoil,qsoil,&
           &zo_init,itimestep,stepbl,&
           &lowlyr,ht,znt,pblh,tsk,qsfc,mavail,&
           &p_phy,p8w,th_phy,t_phy,moist,u_phy,v_phy,pi_phy,rho,tke_myj,&
           &z,z8w,dz8w,num_soil_layers,&
           &ims,ime,jms,jme,kms,kme)

! initialize atmosphere end

      start_seconds = 0

!! allocate and assign soil, surface and atmosphere for ideal end
   
   ELSE  ! will use WRF input instead
      
!! allocate and assign soil, surface and atmosphere for init_gen begin

      idum = wrf_rnd_seed
      
      CALL genctrl(bamex,nz_gen,ns_gen,ntimes_gen,Ninputtimes, &
           &Dx,Cent_lat,Cent_lon,Truelat1,Truelat2,Mminlu,Julday,&
           &Sw_corner_lon,Sw_corner_lat,Projcode,&
           &Z_gen,Z_gen_stag,T_init_gen,U_init_gen,V_init_gen,&
           &Q_init_gen,P_init_gen,&
           &Th2_init_gen,Tsk_init_gen,&
           &U10_init_gen,V10_init_gen,     &
           &Q2_init_gen,Glw_init_gen,Gsw_init_gen,Qsfc_init_gen,&
           &Tslb_init_gen,Smois_init_gen,Tmn_init_gen,&
           &Vegfra_gen,Isltyp_gen,Lu_index_gen,Ivgtyp_gen,&
           &Lat,Lon,Cor,Times_gen,idum, &
           &start_hour_wrf, hour_interval_wrf, start_forecast,&
           &forecast_length,&
           &Control_index1,Control_index2,Control_w,init_gen_file)

      wrf_rnd_seed = idum
      start_seconds = start_forecast * 3600
  
! initialize soil begin

      lu_index=  REAL(lu_index_gen(itime_gen))
      vegfra =   vegfra_gen(itime_gen)
      ivgtyp =   ivgtyp_gen(itime_gen)
      isltyp =   isltyp_gen(itime_gen)

      IF (mminlu=='USGS') THEN
         iswater=16
      ELSE IF (mminlu(1:3)=='OLD') THEN
         iswater=7
      ELSE IF (mminlu(1:3)=='SiB') THEN
         iswater=15
      ELSE IF (mminlu=='LW12') THEN
         iswater=2
      ELSE
         PRINT *,'Unknown Landuse Category',mminlu
         PRINT *,'Stopping'
         STOP
      ENDIF

      CALL landuse_init(lu_index, snowc, Albedo, Albbck, Mavail, &
           &Emiss,  &
           Znt,Z0,Thc,Xland, julday, cent_lat, iswater, mminlu, &
           ids, ide, jds, jde, kds, kde,                       &
           ims, ime, jms, jme, kms, kme,                       &
           its, ite, jts, jte, kts, kte                       )

      seamask=xland-1.
      landmask_input=1.-seamask

      tsk=tsk_init_gen(itime_gen)
      tmn=tmn_init_gen(itime_gen)
      qsfc=qsfc_init_gen(itime_gen)
      sst_input=tsk_init_gen(itime_gen)

      IF (bamex) THEN
         CALL init_soil_depth_2 ( zs_gen , dzs_gen , ns_gen )
      ELSE
         CALL init_soil_depth_1 ( zs_gen , dzs_gen , ns_gen )
      ENDIF

      IF (surface_physics=='SIMPLESCHEME') THEN

         PRINT *,'SIMPLESCHEME cannot be run with init_gen = .TRUE.'
         PRINT *, 'Stopping'
         STOP

      ELSEIF (surface_physics=='SLABSCHEME') THEN

         num_soil_layers=5

         CALL init_soil_depth_1 ( zs , dzs , num_soil_layers )

         IF (bamex) THEN
            CALL init_soil_1_real ( tsk , tmn , tslb , zs , dzs , &
                 &num_soil_layers , 1 , &
                 landmask_input , sst_input , flag_sst , &
                 ids , ide , jds , jde , kds , kde , &
                 ims , ime , jms , jme , kms , kme , &
                 its , ite , jts , jte , kts , kte )
            smois=0.
         ELSE
            DO k=1,num_soil_layers
               tslb(:,k,:)=tslb_init_gen(k,itime_gen)
            ENDDO
            smois=0.
         ENDIF


         keepfr3dflag=0.
         smfr3d=0.
         xice = 0.
         snow = 0.
         snowc = 0.
         canwat = 0.
         sh2o = 0.

      ELSEIF (surface_physics=='LSMSCHEME') THEN
         
         num_soil_layers=4
         
         CALL init_soil_depth_2 ( zs , dzs , num_soil_layers )


         IF (bamex) THEN

            DO k=1,num_soil_layers
               tslb(:,k,:)=tslb_init_gen(k,itime_gen)
               smois(:,k,:)=smois_init_gen(k,itime_gen)
            ENDDO

         ELSE

            st_levels_input(1)=NINT(zs_gen(ns_gen-2)*100.)
            st_levels_input(2)=NINT(zs_gen(ns_gen)*100.)
            sm_levels_input(1)=NINT(zs_gen(ns_gen-2)*100.)
            sm_levels_input(2)=NINT(zs_gen(ns_gen)*100.)
            
            st_input(1,1,2)=tslb_init_gen(ns_gen-2,itime_gen)
            st_input(1,1,3)=tslb_init_gen(ns_gen,itime_gen)
            
            CALL getsmlsm (mminlu, isltyp(1,1), ivgtyp(1,1), &
                 &Smcmin, Smcmax)
            
            sm_input(:,:,2)=MAX(MIN(mavail,smcmax),smcmin)
            sm_input(:,:,3)=&
                 &MAX(deep_soil_moisture*smcmax,smcmin)
            
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
         
         ENDIF

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
              num_soil_layers, restart,                 &
              ids,ide, jds,jde, kds,kde,                &
              ims,ime, jms,jme, kms,kme,                &
              its,ite, jts,jte, kts,kte                 )

      ELSEIF (surface_physics=='RUCLSMSCHEME') THEN
         
         num_soil_layers=6
         
         CALL init_soil_depth_3 ( zs , dzs , num_soil_layers )
         
         IF (bamex) THEN

            st_levels_input(1)=0
            st_levels_input(2)=300
            st_input(:,:,1)=tsk
            st_input(:,:,2)=tmn

            CALL getsmruclsm (mminlu,isltyp(1,1), ivgtyp(1,1), &
                 &Smcmin, Smcmax)

            sm_levels_input(1)=zs_gen(1)*100.
            sm_levels_input(2)=300.

            sm_input(:,:,1)=MAX(MIN(&
                 &smois_init_gen(1,itime_gen),smcmax),smcmin)-smcmin
            sm_input(:,:,2)=MAX(MIN(&
                 &smois_init_gen(ns_gen,itime_gen),smcmax),smcmin)-smcmin

         ELSE

            st_levels_input(1)=0
            st_levels_input(2)=300
            sm_levels_input(1)=0
            sm_levels_input(2)=300
            st_input(:,:,1)=tsk
            st_input(:,:,2)=tmn

            CALL getsmruclsm (mminlu,isltyp(1,1), ivgtyp(1,1), &
                 &Smcmin, Smcmax)
            
            sm_input(:,:,1)=MAX(MIN(mavail,smcmax),smcmin)-smcmin
            sm_input(:,:,2)=&
                 &MAX(deep_soil_moisture*smcmax,smcmin)-smcmin

         ENDIF

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

         IF (bamex) THEN
            smois(:,1,:)=(smois(:,2,:)-smois(:,3,:)*zs_gen(2)/zs_gen(3))/&
                 &(1.-zs_gen(2)/zs_gen(3))
         ENDIF

         keepfr3dflag=0.
         smfr3d=0.
         xice = 0.
         snow = 0.
         snowc = 0.
         canwat = 0.
         sh2o = 0.

         CALL LSMRUCINIT(SMFR3D,TSLB,SMOIS,ISLTYP,&
              num_soil_layers, restart,                 &
              ids,ide, jds,jde, kds,kde,                &
              ims,ime, jms,jme, kms,kme,                &
              its,ite, jts,jte, kts,kte                 )
         
      ENDIF

! initialize soil end
      
     DO i=1,nsplinetimes
         splinetimes(i)=times_gen(1)+(i-1)*splineinterval
      ENDDO

      CALL initgen(bamex,&
           &u_init_gen(:,itime_gen),v_init_gen(:,itime_gen),&
           &t_init_gen(:,itime_gen),q_init_gen(:,itime_gen),&
           &u_init_gen(:,1:ninputtimes),v_init_gen(:,1:ninputtimes),&
           &glw_init_gen(1:ninputtimes),gsw_init_gen(1:ninputtimes),&
           &p_init_gen(:,1:ninputtimes),&
           &th2_init_gen(itime_gen),q2_init_gen(itime_gen),&
           &u10_init_gen(itime_gen),v10_init_gen(itime_gen),&
           &tsk_init_gen(itime_gen),qsfc_init_gen(itime_gen),&
           &z_gen(:,1:ninputtimes),nz_gen,z_g,ninputtimes,&
           &times_gen(1:ninputtimes),&
           &nsplinetimes,splinetimes,&
           &pblh_ref,Stepbl,Lowlyr,Ht,Pblh,&
           &Th_phy,T_phy,Moist,U_phy,V_phy,Tke_myj,&
           &U_g_gen,V_g_gen,Glw_gen,Gsw_gen,P_gen,P8w_gen,&
           &Z,Z8w,Dz8w,&
           &n_moist,P_QV,ims,ime,jms,jme,kms,kme)

      zl2=z(:,2,:)

!! allocate and assign soil, surface and atmosphere for init_gen end
   
   ENDIF

   k_t=0.
   k_m=0.
   uflux=0.
   vflux=0.
   hflux=0.
   qflux=0.

!! initialize pbl begin

   IF (pbltype=='MYJPBLSCHEME') THEN

      OPEN(isfcunit,file=TRIM(outdir)//'/myjsfc.txt')
      OPEN(iprofunit,file=TRIM(outdir)//'/myjprof.txt')
      OPEN(isimil,file=TRIM(outdir)//'/myjsimil.txt')

      CALL MYJSFCINIT(LOWLYR,UST,Znt            &
           ,SEAMASK,XICE,IVGTYP,.FALSE.                &
           ,.TRUE.                           &
           ,IDS,IDE,JDS,JDE,KDS,KDE                    &
           ,IMS,IME,JMS,JME,KMS,KME                    &
           ,ITS,ITE,JTS,JTE,KTS,KTE)
      
      CALL myjpblinit(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN, &
           TKE_MYJ,.FALSE.,.TRUE.,                      &
           ids, ide, jds, jde, kds, kde,         &
           ims, ime, jms, jme, kms, kme,         &
           its, ite, jts, jte, kts, kte          )
      
   ELSEIF (pbltype=='YSUPBLSCHEME') THEN

      OPEN(isfcunit,file=TRIM(outdir)//'/ysusfc.txt')
      OPEN(iprofunit,file=TRIM(outdir)//'/ysuprof.txt')
      OPEN(isimil,file=TRIM(outdir)//'/ysusimil.txt')
      
      CALL sfclayinit(.TRUE.)
      CALL ysuinit(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,    &
           RQCBLTEN,RQIBLTEN,P_QI,               &
           PARAM_FIRST_SCALAR,                   &
           .FALSE.,.TRUE.,                              &
           ids, ide, jds, jde, kds, kde,         &
           ims, ime, jms, jme, kms, kme,         &
           its, ite, jts, jte, kts, kte          )

   ELSEIF (pbltype=='MRFPBLSCHEME') THEN
      OPEN(isfcunit,file=TRIM(outdir)//'/mrfsfc.txt')
      OPEN(iprofunit,file=TRIM(outdir)//'/mrfprof.txt')
      OPEN(isimil,file=TRIM(outdir)//'/mrfsimil.txt')

      CALL sfclayinit(.TRUE.)
      CALL mrfinit(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,    &
           RQCBLTEN,RQIBLTEN,P_QI,               &
           PARAM_FIRST_SCALAR,                   &
           .FALSE.,.TRUE.,                       &
           ids, ide, jds, jde, kds, kde,         &
           ims, ime, jms, jme, kms, kme,         &
           its, ite, jts, jte, kts, kte          )
   ENDIF

END SUBROUTINE INIT_WRF

!------------------------------------------------------------------
SUBROUTINE WRF(dart_seconds,dart_days)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: dart_days, dart_seconds
   LOGICAL             :: debug1, debug2, debug3, debug4, debug5

   debug1 = .false.  ! before sfc
   debug2 = .false.  ! before pbl
   debug3 = .false.  ! after sfc
   debug4 = .false.  ! after pbl
   debug5 = .false.  ! after pbl and tendencies applied

!!! basic stuff

   a=1.-.25*dt**2*cor**2
   b=dt*cor
   c=.5*dt**2*cor**2
   d=.5*dt**2*cor
   e=1.+.25*dt**2*cor**2

!  here we compute the time step to keep pace with dart
   itimestep = (dart_days * 86400 + dart_seconds - start_seconds) / dt + 1

!   DO itimestep=1,ntime+1

      time=timeo+REAL(itimestep-1)*dt
     
      IF (init_gen) THEN

         ihour=INT(time/3600.+epsilon)
         imin=ihour+1-start_forecast

         imax=MIN(ihour+2-start_forecast,nsplinetimes)
         fract=MOD(REAL(itimestep-1)*dt,3600.)/3600.
         glw=glw_gen(imin)+(glw_gen(imax)-glw_gen(imin))*fract
         gsw=gsw_gen(imin)+(gsw_gen(imax)-gsw_gen(imin))*fract
         DO i=ims,ime
            DO j=jms,jme
               p_phy(i,:,j)=p_gen(:,imin)+&
                    &(p_gen(:,imax)-p_gen(:,imin))*fract
               p8w(i,:,j)=p8w_gen(:,imin)+&
                    &(p8w_gen(:,imax)-p8w_gen(:,imin))*fract
               PSFC(I,J)=p8w(I,kts,J)
               pi_phy(i,:,j)=(p_phy(i,:,j)/p1000mb)**rcp
            ENDDO
         ENDDO

         u_g(:)=u_g_gen(:,imin)+&
              &(u_g_gen(:,imax)-u_g_gen(:,imin))*fract
         v_g(:)=v_g_gen(:,imin)+&
              &(v_g_gen(:,imax)-v_g_gen(:,imin))*fract

      ELSE
         glw=glwfunc(time)
         gsw=gswfunc(time,albedo(1,1))
      ENDIF

      IF ( debug1 ) then
         write(65,*)'BEFORE SFC'
         write(65,*)'rad ',gsw,glw
         write(65,*)'u,v ',u_g(1),v_g(1)
         write(65,*)'p,p8 ',p_phy(1,1,1),p8w(1,1,1)
         write(65,*)'l1 ',sfctype,surface_physics,pbltype
         write(65,*)'l2 ',ts_ref,ps_ref,dtamplitude_ref,mavail_ref,time,p_qc,p_qv
         write(65,*)'l3 ',ACSNOM,ACSNOW,AKHS,AKMS,ALBEDO,BR,CANWAT,CAPG
         write(65,*)'l4 ',CHKLOWQ,DT,DX,DZ8W(1,1,1),DZS,EMISS,GLW
         write(65,*)'l5 ',GRDFLX,GSW,GZ1OZ0,HFX,HOL,HT,IFSNOW,ISFFLX
         write(65,*)'l6 ',ISLTYP,ITIMESTEP,IVGTYP,LOWLYR,MAVAIL,MOL    
         write(65,*)'MOIST ',moist(1,1,1,:)
         write(65,*)'l7 ',NUM_SOIL_LAYERS,n_moist
         write(65,*)'l8 ',P8W(1,1,1),PBLH,PI_PHY(1,1,1),PSHLTR,PSIH 
         write(65,*)'l9 ',PSIM,P_PHY(1,1,1),Q10,Q2,QFX,QSFC,QSHLTR,QZ0,RAINBL         
         write(65,*)'l10 ',RAINCV,RAINNCV,REGIME,RHO(1,1,1),SFCEVP,SFCEXC,SFCRUNOFF    
         write(65,*)'l11 ',SMOIS(1,1,1),SMSTAV,SMSTOT,SNOALB,SNOW,SNOWC,SNOWH,STEPBL   
         write(65,*)'l12 ',T2,TH10,TH2,THC,THZ0,TH_PHY(1,1,1),TMN,TSHLTR,TSK,TSLB(1,1,1) 
         WRITE(65,*)'l13 ',T_PHY(1,1,1),U10,UDRUNOFF,UST,UZ0,U_FRAME,U_PHY(1,1,1),V10,VEGFRA  
         WRITE(65,*)'l14 ',VZ0,V_FRAME,V_PHY(1,1,1),WARM_RAIN,WSPD,XICE,XLAND,Z(1,1,1),ZNT,ZS(1) 
         write(65,*)'l15 ',CT,TKE_MYJ(1,1,1)             
!         write(65,*)'l16',ALBBCK,LH,SH2O,SHDMAX,SHDMIN,Z0                      
         write(65,*)'l17',flqc,flhc,qsg,qvg,qcg,soilt1,tsnav                   
         write(65,*)'l18',SMFR3D(1,1,1),KEEPFR3DFLAG(1,1,1)                    
         write(65,*)'l19',PSFC                                                 
         write(65,*)'l20',POTEVP,SNOPCX,SOILTB
      ENDIF
      rho=p_phy/(r_d*t_phy*(1.+SVP1*moist(:,:,:,P_QV)))

      CALL surface_driver(sfctype,surface_physics,pbltype&
           &,ts_ref,ps_ref,dtamplitude_ref,mavail_ref,time,p_qc,p_qv&
           &,ACSNOM,ACSNOW,AKHS,AKMS,ALBEDO,BR,CANWAT,CAPG        &
           &,CHKLOWQ,DT,DX,DZ8W,DZS,EMISS,GLW        &
           &,GRDFLX,GSW,GZ1OZ0,HFX,HOL,HT,IFSNOW,ISFFLX           &
           &,ISLTYP,ITIMESTEP,IVGTYP,LOWLYR,MAVAIL,MOIST,MOL    &
           &,NUM_SOIL_LAYERS,n_moist&
           &,P8W,PBLH,PI_PHY,PSHLTR,PSIH &
           &,PSIM,&
!mp@ begin
                     &psimfac,psihfac,psiqfac,zl2,&
!mp@ end
           &P_PHY,Q10,Q2,QFX,QSFC,QSHLTR,QZ0,RAINBL         &
           &,RAINCV,RAINNCV,REGIME,RHO,SFCEVP,SFCEXC,SFCRUNOFF    &
           &,SMOIS,SMSTAV,SMSTOT,SNOALB,SNOW,SNOWC,SNOWH,STEPBL   &
           &,T2,TH10,TH2,THC,THZ0,TH_PHY,TMN,TSHLTR,TSK,TSLB      &
           &,T_PHY,U10,UDRUNOFF,UST,UZ0,U_FRAME,U_PHY,V10,VEGFRA  &
           &,VZ0,V_FRAME,V_PHY,WARM_RAIN,WSPD,XICE,XLAND,Z,ZNT,ZS &
           &,CT,TKE_MYJ                                           &
           &,ALBBCK,LH,SH2O,SHDMAX,SHDMIN,Z0                      &
           &,flqc,flhc,qsg,qvg,qcg,soilt1,tsnav                   &
           &,SMFR3D,KEEPFR3DFLAG                                  &
           &,PSFC                                                 &
           &,ids,ide,jds,jde,kds,kde                              &
           &,ims,ime,jms,jme,kms,kme                              &
           &,i_start,i_end,j_start,j_end,kts,kte,num_tiles        &
           &,POTEVP,SNOPCX,SOILTB)

      IF ( debug2 ) then
         write(75,*)'AFTER SFC'
         write(75,*)'rad ',gsw,glw
         write(75,*)'u,v ',u_g(1),v_g(1)
         write(75,*)'p,p8 ',p_phy(1,1,1),p8w(1,1,1)
         write(75,*)'l1 ',sfctype,surface_physics,pbltype
         write(75,*)'l2 ',ts_ref,ps_ref,dtamplitude_ref,mavail_ref,time,p_qc,p_qv
         write(75,*)'l3 ',ACSNOM,ACSNOW,AKHS,AKMS,ALBEDO,BR,CANWAT,CAPG
         write(75,*)'l4 ',CHKLOWQ,DT,DX,DZ8W(1,1,1),DZS,EMISS,GLW
         write(75,*)'l5 ',GRDFLX,GSW,GZ1OZ0,HFX,HOL,HT,IFSNOW,ISFFLX
         write(75,*)'l6 ',ISLTYP,ITIMESTEP,IVGTYP,LOWLYR,MAVAIL,MOL    
         write(75,*)'MOIST ',moist(1,1,1,:)
         write(75,*)'l7 ',NUM_SOIL_LAYERS,n_moist
         write(75,*)'l8 ',P8W(1,1,1),PBLH,PI_PHY(1,1,1),PSHLTR,PSIH 
         write(75,*)'l9 ',PSIM,P_PHY(1,1,1),Q10,Q2,QFX,QSFC,QSHLTR,QZ0,RAINBL         
         write(75,*)'l10 ',RAINCV,RAINNCV,REGIME,RHO(1,1,1),SFCEVP,SFCEXC,SFCRUNOFF    
         write(75,*)'l11 ',SMOIS(1,1,1),SMSTAV,SMSTOT,SNOALB,SNOW,SNOWC,SNOWH,STEPBL   
         write(75,*)'l12 ',T2,TH10,TH2,THC,THZ0,TH_PHY(1,1,1),TMN,TSHLTR,TSK,TSLB(1,1,1) 
         WRITE(75,*)'l13 ',T_PHY(1,1,1),U10,UDRUNOFF,UST,UZ0,U_FRAME,U_PHY(1,1,1),V10,VEGFRA  
         WRITE(75,*)'l14 ',VZ0,V_FRAME,V_PHY(1,1,1),WARM_RAIN,WSPD,XICE,XLAND,Z(1,1,1),ZNT,ZS(1) 
         write(75,*)'l15 ',CT,TKE_MYJ(1,1,1)             
!         write(75,*)'l16',ALBBCK,LH,SH2O,SHDMAX,SHDMIN,Z0                      
         write(75,*)'l17',flqc,flhc,qsg,qvg,qcg,soilt1,tsnav                   
         write(75,*)'l18',SMFR3D(1,1,1),KEEPFR3DFLAG(1,1,1)                    
         write(75,*)'l19',PSFC                                                 
         write(75,*)'l20',POTEVP,SNOPCX,SOILTB
      ENDIF
      

      IF (debug3 ) then
         write(66,*)'BEFORE PBL'
         write(66,*)'l1', pbltype,surface_physics
         write(66,*)'l2',itimestep,dt,u_frame,v_frame
         write(66,*)'l3',RUBLTEN(1,1,1),RVBLTEN(1,1,1),RTHBLTEN(1,1,1)
         write(66,*)'l4',RQVBLTEN(1,1,1),RQCBLTEN(1,1,1),RQIBLTEN(1,1,1)
         write(66,*)'l5',TSK,XLAND,ZNT,HT
         write(66,*)'l6',UST,HOL,MOL,PBLH
         write(66,*)'l7',HFX,QFX,REGIME,GRDFLX
         write(66,*)'l8',u_phy(1,1,1),v_phy(1,1,1),th_phy(1,1,1),rho(1,1,1)
         write(66,*)'MOIST ',moist(1,1,1,:)
         write(66,*)'l9',p_phy(1,1,1),pi_phy(1,1,1),p8w(1,1,1),t_phy(1,1,1),dz8w(1,1,1),z(1,1,1)
         write(66,*)'l10',TKE_MYJ(1,1,1),AKHS,AKMS
         write(66,*)'l11',THZ0,QZ0,UZ0,VZ0,QSFC,LOWLYR
         write(66,*)'l12',PSIM, PSIH, GZ1OZ0, WSPD, BR, CHKLOWQ
         write(66,*)'l13',DX,n_moist
         write(66,*)'l14',STEPBL,warm_rain
         write(66,*)'l15',KPBL,CT,LH,SNOW,XICE
         write(66,*)'l16',P_QI,P_QV,P_QC,PARAM_FIRST_SCALAR
         write(66,*)'l17',uflux(1),vflux(1),hflux(1),qflux(1), k_t(1),k_m(1)
      ENDIF

      IF (surface_physics=='SIMPLESCHEME') THEN
         qvg=qsfc
      ELSEIF (surface_physics=='SLABSCHEME') THEN
         qvg=qsfc/(1-qsfc)
      ELSEIF (surface_physics=='LSMSCHEME') THEN
         qvg=qsfc/(1-qsfc)
      ENDIF

      CALL  pbl_driver(pbltype,surface_physics&
           &,itimestep,dt,u_frame,v_frame                    &
           &,RUBLTEN,RVBLTEN,RTHBLTEN                            &
           &,RQVBLTEN,RQCBLTEN,RQIBLTEN                          &
           &,TSK,XLAND,ZNT,HT                                    &
           &,UST,HOL,MOL,PBLH                                    &
           &,HFX,QFX,REGIME,GRDFLX                               &
           &,u_phy,v_phy,th_phy,rho,moist                      &
           &,p_phy,pi_phy,p8w,t_phy,dz8w,z                       &
           &,TKE_MYJ,AKHS,AKMS                                   &
           &,THZ0,QZ0,UZ0,VZ0,QSFC,LOWLYR                        &
           &,PSIM, PSIH, &
!mp@ begin
                  &psimfac,psihfac,psiqfac,          &
!mp@ end
           &GZ1OZ0, WSPD, BR, CHKLOWQ               &
           &,DX,n_moist                                        &
           &,STEPBL,warm_rain                                    &
           &,KPBL,CT,LH,SNOW,XICE                                &
           &,P_QI,P_QV,P_QC,PARAM_FIRST_SCALAR    &
           &,uflux,vflux,hflux,qflux, k_t,k_m&
           &,ids,ide, jds,jde, kds,kde                      &
           &,ims,ime, jms,jme, kms,kme                      &
           &,i_start,i_end, j_start,j_end, kts,kte, num_tiles)

      IF (debug4 ) then
         write(76,*)'AFTER PBL'
         write(76,*)'l1', pbltype,surface_physics
         write(76,*)'l2',itimestep,dt,u_frame,v_frame
         write(76,*)'l3',RUBLTEN(1,1,1),RVBLTEN(1,1,1),RTHBLTEN(1,1,1)
         write(76,*)'l4',RQVBLTEN(1,1,1),RQCBLTEN(1,1,1),RQIBLTEN(1,1,1)
         write(76,*)'l5',TSK,XLAND,ZNT,HT
         write(76,*)'l6',UST,HOL,MOL,PBLH
         write(76,*)'l7',HFX,QFX,REGIME,GRDFLX
         write(76,*)'l8',u_phy(1,1,1),v_phy(1,1,1),th_phy(1,1,1),rho(1,1,1)
         write(76,*)'MOIST ',moist(1,1,1,:)
         write(76,*)'l9',p_phy(1,1,1),pi_phy(1,1,1),p8w(1,1,1),t_phy(1,1,1),dz8w(1,1,1),z(1,1,1)
         write(76,*)'l10',TKE_MYJ(1,1,1),AKHS,AKMS
         write(76,*)'l11',THZ0,QZ0,UZ0,VZ0,QSFC,LOWLYR
         write(76,*)'l12',PSIM, PSIH, GZ1OZ0, WSPD, BR, CHKLOWQ
         write(76,*)'l13',DX,n_moist
         write(76,*)'l14',STEPBL,warm_rain
         write(76,*)'l15',KPBL,CT,LH,SNOW,XICE
         write(76,*)'l16',P_QI,P_QV,P_QC,PARAM_FIRST_SCALAR
         write(76,*)'l17',uflux(1),vflux(1),hflux(1),qflux(1), k_t(1),k_m(1)
      ENDIF
      DO i=ims,ime
         DO j=jms,jme
            DO k=kts,kte
               u_phy(i,k,j)=(&
                    &a*u_phy(i,k,j)+&
                    &b*(v_phy(i,k,j)-v_g(k))+&
                    &RUBLTEN(i,k,j)*dt+&
                    &c*u_g(k)+&
                    &d*RVBLTEN(i,k,j)&
                    &)/e
               v_phy(i,k,j)=(&
                    &a*v_phy(i,k,j)-&
                    &b*(u_phy(i,k,j)-u_g(k))+&
                    &RVBLTEN(i,k,j)*dt-&
                    &c*v_g(k)-&
                    &d*RUBLTEN(i,k,j)&
                    &)/e
               th_phy(i,k,j)=th_phy(i,k,j)+RTHBLTEN(i,k,j)*dt
               moist(i,k,j,P_QV)=moist(i,k,j,P_QV)+RQVBLTEN(i,k,j)*dt
               t_phy(i,k,j)=th_phy(i,k,j)*pi_phy(i,k,j)
               rho(i,k,j)=p_phy(i,k,j)/(r_d*&
                    &t_phy(i,k,j)*(1.+SVP1*moist(i,k,j,P_QV)))
            ENDDO
         ENDDO
      ENDDO

      IF (debug5 ) then
         write(85,*)'AFTER PBL and TENDS'
         write(85,*)'COEFFS ',a,b,c,d
         write(85,*)'l1', pbltype,surface_physics
         write(85,*)'l2',itimestep,dt,u_frame,v_frame
         write(85,*)'l3',RUBLTEN(1,1,1),RVBLTEN(1,1,1),RTHBLTEN(1,1,1)
         write(85,*)'l4',RQVBLTEN(1,1,1),RQCBLTEN(1,1,1),RQIBLTEN(1,1,1)
         write(85,*)'l5',TSK,XLAND,ZNT,HT
         write(85,*)'l6',UST,HOL,MOL,PBLH
         write(85,*)'l7',HFX,QFX,REGIME,GRDFLX
         write(85,*)'l8',u_phy(1,1,1),v_phy(1,1,1),th_phy(1,1,1),rho(1,1,1)
         write(85,*)'MOIST ',moist(1,1,1,:)
         write(85,*)'l9',p_phy(1,1,1),pi_phy(1,1,1),p8w(1,1,1),t_phy(1,1,1),dz8w(1,1,1),z(1,1,1)
         write(85,*)'l10',TKE_MYJ(1,1,1),AKHS,AKMS
         write(85,*)'l11',THZ0,QZ0,UZ0,VZ0,QSFC,LOWLYR
         write(85,*)'l12',PSIM, PSIH, GZ1OZ0, WSPD, BR, CHKLOWQ
         write(85,*)'l13',DX,n_moist
         write(85,*)'l14',STEPBL,warm_rain
         write(85,*)'l15',KPBL,CT,LH,SNOW,XICE
         write(85,*)'l16',P_QI,P_QV,P_QC,PARAM_FIRST_SCALAR
         write(85,*)'l17',uflux(1),vflux(1),hflux(1),qflux(1), k_t(1),k_m(1)
      ENDIF
!   ENDDO

END SUBROUTINE wrf

SUBROUTINE output_wrf_profiles()

  IMPLICIT NONE
  
  IF (init_gen) THEN
     IF (MOD(time,outgeninterval) < epsilon) THEN
        WRITE(ncunit,'(20f16.7)')time/3600.
        DO k=kts,kte-1
           WRITE(ncunit,'(20f16.7)')z(1,k,1),u_phy(1,k,1),&
                &v_phy(1,k,1),th_phy(1,k,1),moist(1,k,1,P_QV)
        ENDDO
     ENDIF
  ENDIF

  IF (MOD(NINT(time-timeo),outinterval)==0) THEN 
     WRITE(isfcunit,'(30f16.7)')time/3600.,tsk(1,1),tslb(1,:,1),&
          &tmn(1,1),qsfc(1,1),SMOIS(1,:,1),gsw(1,1),ust,hfx,xlv*qfx
  ENDIF


  IF (itimestep==ntime+1) THEN

     CLOSE(isfcunit)
     CLOSE(ncunit)

     rhosfc=p_phy(1,1,1)/(r_d*t_phy(1,1,1)*(1.+svp1*moist(1,1,1,P_QV)))

     DO k=kts,kte-1
        WRITE(iprofunit,'(20f16.7)')z(1,k,1),z8w(1,k,1),u_phy(1,k,1),&
             &v_phy(1,k,1),SQRT(u_phy(1,k,1)**2+v_phy(1,k,1)**2),&
             &th_phy(1,k,1),moist(1,k,1,P_QV),tke_myj(1,k,1),&
             &k_t(k),k_m(k),z8w(1,k,1)*cor/ust(1,1),&
             &-uflux(k)/ust**2,vflux(k)/ust**2,hflux(k)/hfx*cp*rhosfc,&
             &SQRT(uflux(k)**2+vflux(k)**2)/ust**2,qflux(k)/qfx*rhosfc
     ENDDO

     CLOSE(iprofunit)
  
     zo=znt(1,1)

     IF(ABS(hfx(1,1)) < epsilon) hfx=epsilon
     IF(ABS(qfx(1,1)) < epsilon) qfx=epsilon
     
     hvfx=hfx(1,1)/(cp*rhosfc)+SVP1*t_phy(1,1,1)*qfx(1,1)
     
     thvstar=-hvfx/(ust(1,1)*rhosfc)
     
     qstar=-qfx(1,1)/(ust(1,1)*rhosfc)
     
     
     IF (ABS(hvfx).LT.1.e-3) THEN
        dlmonin=1.e8
     ELSE
        dlmonin=-ust(1,1)**3/(karman*g/th_phy(1,1,1)*hvfx)
     ENDIF
     
     DO k=kts,k_simil
        
        z_level=z(1,k,1)
        
        zol=z_level/dlmonin
        
        IF (dlmonin.LT.0) THEN
           x=(1.-16.*zol)**0.25
           psimzl=2*LOG(0.5*(1+x))+LOG(0.5*(1+x*x)) &
                &        -2.*ATAN(x)+2.*ATAN(1.)
           
           y=(1.-16.*zol)**0.5
           psihzl=2.*LOG(0.5*(1+y))
        ELSE
           psimzl=-b_stable*(zol-c_stable/d_stable)*EXP(-d_stable*zol)-&
                &a_stable*zol-b_stable*c_stable/d_stable
           psihzl=-b_stable*(zol-c_stable/d_stable)*EXP(-d_stable*zol)-&
                &(SQRT(1.+a_stable*b_stable*zol))**3-&
                &b_stable*c_stable/d_stable+1.
!         psimzl=-5.*zol
!         psihzl=-5.*zol
        ENDIF
        
        u_simil(k)=ust(1,1)/karman*(LOG(z_level/zo)-psimzl)*&
             &u_phy(1,1,1)/&
             &SQRT(u_phy(1,1,1)**2+v_phy(1,1,1)**2)
        v_simil(k)=ust(1,1)/karman*(LOG(z_level/zo)-psimzl)*&
             &v_phy(1,1,1)/&
             &SQRT(u_phy(1,1,1)**2+v_phy(1,1,1)**2)
        q_simil(k)=qstar/karman*(LOG(z_level/zo)-LOG(2.)-psihzl)+qsfc(1,1)
        th_simil(k)=(thvstar/karman*(LOG(z_level/zo)-LOG(2.)-psihzl)+&
             &tsk(1,1)*(1.+SVP1*qsfc(1,1)))/(1.+SVP1*moist(1,k,1,P_QV))
        
        WRITE(isimil,'(20f16.7)')&
             
             &LOG(z_level/zo),&
             
             &LOG(z_level/zo)-LOG(2.),&
             
             &u_simil(k)/&
             &(ust(1,1)*u_simil(1)/&
             &SQRT(u_simil(1)**2+v_simil(1)**2)),&
             
             &u_phy(1,k,1)/&
             &(ust(1,1)*u_phy(1,1,1)/&
             &SQRT(u_phy(1,1,1)**2+v_phy(1,1,1)**2)),&
             
             &v_simil(k)/&
             &(ust(1,1)*v_simil(1)/&
             &SQRT(u_simil(1)**2+v_simil(1)**2)),&
             
             &v_phy(1,k,1)/&
             &(ust(1,1)*v_phy(1,1,1)/&
             &SQRT(u_phy(1,1,1)**2+v_phy(1,1,1)**2)),&
             
             &SQRT(u_simil(k)**2+v_simil(k)**2)/ust,&
             
             &SQRT(u_phy(1,k,1)**2+v_phy(1,k,1)**2)/ust,&
             
             &(th_simil(k)*(1.+SVP1*moist(1,k,1,P_QV))-&
             &tsk(1,1)*(1.+SVP1*qsfc(1,1)))/thvstar,&
             
             &(th_phy(1,k,1)*(1.+SVP1*moist(1,k,1,P_QV))-&
             &tsk(1,1)*(1.+SVP1*qsfc(1,1)))/thvstar,&
             
!           &(th_simil(k)*(1.+SVP1*moist(1,k,1,P_QV))-&
!           &tsk(1,1)*(1.+SVP1*qsfc(1,1)))/&
!           &(th_phy(1,k,1)*(1.+SVP1*moist(1,k,1,P_QV))-&
!           &tsk(1,1)*(1.+SVP1*qsfc(1,1))),&
          
             &(q_simil(k)-qsfc(1,1))/qstar,&
             &(moist(1,k,1,P_QV)-qsfc(1,1))/qstar !,&
        
!           &(q_simil(k)-qsfc(1,1))/(moist(1,k,1,P_QV)-qsfc(1,1))
        
        
     ENDDO

     CLOSE(isimil)
     
  ENDIF
  
END SUBROUTINE output_wrf_profiles


END MODULE module_wrf

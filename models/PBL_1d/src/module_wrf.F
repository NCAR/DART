MODULE module_wrf

! Primary WRF column model driver.  All memory allocations are done in here.

  USE module_model_constants
  USE module_namelist
  USE module_ideal
  USE module_initialize
  USE module_init_soil_ideal
  USE module_init_soil_real
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
  USE module_wrf_init_and_bc
  USE module_snd_init_and_bc
  USE module_luse_init
  USE module_getsm
  
  IMPLICIT NONE
 
   INTEGER                             :: isfcunit=55,&
                                          iprofunit=56,&
                                          isimil=57,&
                                          ncunit=58

   INTEGER, PARAMETER                  :: ids=1,&
                                          ide=1,&
                                          jds=1,&
                                          jde=1,&
                                          kds=1,&
                                          ims=1,&
                                          ime=1,&
                                          jms=1,&
                                          jme=1,&
                                          kms=1,&
                                          its=1,&
                                          ite=1,&
                                          jts=1,&
                                          jte=1,&
                                          kts=1,&
                                          num_tiles=1

   INTEGER                             :: kde, kme, kte

   INTEGER, PARAMETER                  :: k_simil=5,&
                                          PARAM_FIRST_SCALAR = 2

   INTEGER                             :: niter,iter

   REAL, DIMENSION(1:k_simil)          :: u_simil,v_simil,&
                                          q_simil,th_simil

   REAL                                :: hvfx,thvstar,qstar,&
                                          dlmonin,z_level,zo,&
                                          zol,x,psimzl,y,psihzl,&
                                          rhosfc

   REAL, PARAMETER                     :: a_stable=1.,&
                                          b_stable=2./3.,&
                                          c_stable=5.,&
                                          d_stable=.35

   INTEGER                             :: ntime,num_soil_layers,&
                                          nsplinetimes, &
                                          nsplinetimes_flux, &
                                          nsplinetimes_smos, &
                                          n1dsplines, &
                                          iitime,&
                                          imin,imax, &
                                          imin_flux,imax_flux,&
                                          &imin_smos,imax_smos

   INTEGER, DIMENSION(num_tiles)       :: i_start=1,i_end=1, &
                                          j_start=1 ,j_end=1 
! forcing netCDF file IDs
   INTEGER                             :: ncid_f, ncid_soil, &
                                          ncid_flux, ncid_smos

! forcing number for soil, atmos, and times
   INTEGER                             :: ns_f, nz_f, nt_f, nt_f_flux, &
                                          nt_f_soil, nt_f_smos 

   INTEGER                             :: itimestep,STEPBL,itime_f=1

   INTEGER, DIMENSION( ims:ime , jms:jme )  :: LOWLYR

   LOGICAL                             ::   warm_rain

   REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  p_phy, &
                                           pi_phy, &
                                           p8w, &
                                           rho, &
                                           t_phy, &
                                           u_phy, &
                                           v_phy, &
                                           dz8w, &
                                           z, &
                                           th_phy

   REAL, DIMENSION( ims:ime, jms:jme )  ::   PSFC

   REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: z8w

   REAL, ALLOCATABLE, DIMENSION(:,:,:,:):: moist

   REAL, DIMENSION( ims:ime , jms:jme ):: XLAND, &
                                          HT, &
                                          PSIM, &
                                          PSIH, &
                                          GZ1OZ0, &
                                          BR, &
                                          CHKLOWQ

   REAL, DIMENSION( ims:ime , jms:jme ):: psimfac,psihfac,psiqfac,zl2

   REAL, DIMENSION( ims:ime, jms:jme ) :: TSK, &
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
   
   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: RUBLTEN, &
                                          RVBLTEN, &
                                          RTHBLTEN, &
                                          RQVBLTEN, &
                                          RQCBLTEN, &
                                          RQIBLTEN, &
                                          TKE_MYJ
  
   REAL                                ::  u_frame=0.,v_frame=0.
   
   INTEGER, DIMENSION( ims:ime , jms:jme ) :: KPBL
 
   REAL, DIMENSION( ims:ime , jms:jme ) :: landmask_input,sst_input,&
                                           cpm,chs

   LOGICAL                              :: flag_sst=.TRUE.

   REAL, DIMENSION( ims:ime , jms:jme ) :: XICE, SEAMASK, CT, SNOW, LH
   

   REAL                                 :: DTMIN,DTBL

   INTEGER                              :: i,J,K,NK,jj,ij

   REAL                                 :: dx,time,cor,fract,fract_flux

   CHARACTER(len=4)                     :: mminlu

   INTEGER                              :: julday

   INTEGER, DIMENSION( ims:ime , jms:jme ) :: ISLTYP,IVGTYP

   INTEGER                              :: iswater   

   REAL                                 :: smcmin,smcmax

   REAL, DIMENSION( ims:ime, jms:jme )  :: ACSNOW,ACSNOM,ALBEDO, &
                                           CANWAT,CAPG,EMISS, &
                                           GLW,GSW,MAVAIL,Q10,Q2,&
                                           T2,TH10,TH2,SNOWC,U10, &
                                           V10,z0,PSHLTR,QSHLTR, &
                                           RAINBL,RAINCV,RAINNCV,&
                                           SFCEVP,SFCEXC,THC,cs,TMN,&
                                           SFCRUNOFF,SNOWH,POTEVP,&
                                           SNOPCX,SOILTB,SOILT1,&
                                           TSNAV,QSG,QVG,QCG,FLHC,&
                                           FLQC,SNOALB,SMSTAV,&
                                           SMSTOT,TSHLTR,UDRUNOFF,&
                                           VEGFRA,ALBBCK,SHDMAX,SHDMIN,&
                                           Maxm,Minm,Erate
   
   REAL, DIMENSION( ims:ime, jms:jme ) ::  lu_index

   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: TSLB,SMOIS,KEEPFR3DFLAG, &
                                          SMFR3D,SH2O

   REAL,  ALLOCATABLE, DIMENSION(:)    :: DZS,ZS

   INTEGER                             :: num_st_levels_input, &
                                          num_sm_levels_input
   
   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: st_input,sm_input
   INTEGER, ALLOCATABLE, DIMENSION(:)  :: st_levels_input,sm_levels_input

   REAL, ALLOCATABLE, DIMENSION(:)     :: z_init,u_init,v_init,t_init,&
                                          th_init,exn_init,q_init, &
                                          p_init,rho_init

   REAL                                :: zo_init,tsoil,thsoil,qsoil

   REAL, ALLOCATABLE, DIMENSION(:)     :: z8w_init, p8w_init

   REAL, ALLOCATABLE, DIMENSION(:)     :: u_mid,v_mid

! input map information if horizontal interpolation is necessary
   REAL                                :: cent_lat,cent_lon, &
                                          truelat1,truelat2,&
                                          sw_corner_lon,sw_corner_lat, &
                                          lat,lon

   INTEGER                             :: projcode

   REAL, ALLOCATABLE, DIMENSION(:)     :: uflux,vflux,hflux,qflux,&
                                          k_t,k_m

   REAL, ALLOCATABLE, DIMENSION(:)     :: u_g,v_g
   
   REAL                                :: a,b,c,d,e
   REAL, PARAMETER                     :: deltat = 2.*3600

   LOGICAL, PARAMETER                  :: FNDSOILW=.TRUE.,&
                                          FNDSNOWH=.TRUE.,&
                                          restart=.FALSE.

! variables controlling the random selection of profiles, etc
   INTEGER                             :: idum, start_seconds 

! number of available forcing times from initialization
   REAL, DIMENSION(:), ALLOCATABLE     :: times_f, times_f_flux, &
                                          times_f_soil, times_f_smos

! time series of profiles for forcing
   REAL, DIMENSION(:,:), ALLOCATABLE   :: u_init_f,v_init_f,&
                                          t_init_f,th_init_f,&
                                          exn_init_f,q_init_f,&
                                          p_init_f

! results of temporal and spatial interpolations
   REAL, ALLOCATABLE, DIMENSION(:)     :: splinetimes, splinetimes_flux,&
        &splinetimes_smos
   REAL, ALLOCATABLE, DIMENSION(:,:)   :: u_g_f,v_g_f,p_f,p8w_f,&
                                          t_f_uadv,t_f_vadv,t_f_wadv,&
                                          q_f_uadv,q_f_vadv,q_f_wadv,&
                                          u_f_uadv,u_f_vadv,u_f_wadv,&
                                          v_f_uadv,v_f_vadv,v_f_wadv
   

   REAL, ALLOCATABLE, DIMENSION(:)     :: glw_f,gsw_f,precip_f

   REAL, DIMENSION(:), ALLOCATABLE     :: th2_init_f,tsk_init_f,&
                                          t2_init_f, &
                                          u10_init_f,v10_init_f,&
                                          q2_init_f,glw_init_f, &
                                          precip_init_f, &
                                          glw_up_init_f, &
                                          gsw_init_f,qsfc_init_f,&
                                          tmn_init_f,vegfra_f

   INTEGER, DIMENSION(:), ALLOCATABLE  :: isltyp_f,ivgtyp_f,lu_index_f

   REAL, DIMENSION(:,:), ALLOCATABLE   :: tslb_init_f,smois_init_f
   REAL, DIMENSION(:), ALLOCATABLE     :: zs_f,dzs_f
   REAL, DIMENSION(:,:), ALLOCATABLE   :: z_f
   REAL, DIMENSION(:,:), ALLOCATABLE   :: z_f_stag
   
contains

SUBROUTINE STATIC_INIT_WRF(allocate_wrf)

  LOGICAL,INTENT(INOUT)           :: allocate_wrf

  kde = nz
  kme = nz
  kte = nz - 1

  IF (.NOT.init_f) allocate_wrf=.TRUE.

  IF ( allocate_wrf ) THEN
    ALLOCATE(moist( ims:ime, kms:kme, jms:jme, n_moist))
    ALLOCATE(p_phy(ims:ime, kms:kme, jms:jme),            &
             pi_phy(ims:ime, kms:kme, jms:jme), &
             p8w(ims:ime, kms:kme, jms:jme), &
             rho(ims:ime, kms:kme, jms:jme), &
             t_phy(ims:ime, kms:kme, jms:jme), &
             u_phy(ims:ime, kms:kme, jms:jme), &
             v_phy(ims:ime, kms:kme, jms:jme), &
             dz8w(ims:ime, kms:kme, jms:jme), &
             z(ims:ime, kms:kme, jms:jme), &
             th_phy(ims:ime, kms:kme, jms:jme) )
    ALLOCATE(z8w(ims:ime, kms:kme+1, jms:jme ))
    ALLOCATE(RUBLTEN(ims:ime, kms:kme, jms:jme),            &
             RVBLTEN(ims:ime, kms:kme, jms:jme), &
             RTHBLTEN(ims:ime, kms:kme, jms:jme), &
             RQVBLTEN(ims:ime, kms:kme, jms:jme), &
             RQCBLTEN(ims:ime, kms:kme, jms:jme), &
             RQIBLTEN(ims:ime, kms:kme, jms:jme), &
             TKE_MYJ(ims:ime, kms:kme, jms:jme) )
    ALLOCATE(z_init(1:nz), u_init(1:nz), v_init(1:nz), t_init(1:nz), &
             th_init(1:nz), exn_init(1:nz), q_init(1:nz), p_init(1:nz), &
             rho_init(1:nz))
    ALLOCATE(z8w_init(1:nz+1), p8w_init(1:nz+1))
    ALLOCATE(uflux(1:nz), vflux(1:nz), hflux(1:nz), &
             qflux(1:nz), k_t(1:nz),k_m(1:nz))
    ALLOCATE(u_g(1:nz), v_g(1:nz))

!   Explicit zeroing is unfortunately necessary for PGF 
    moist = 0.0   
    p_phy = 0.0
    pi_phy = 0.0
    p8w = 0.0
    rho = 0.0
    t_phy = 0.0
    u_phy = 0.0
    v_phy = 0.0
    dz8w = 0.0
    z = 0.0
    th_phy = 0.0
    z8w = 0.0
    RUBLTEN = 0.0
    RVBLTEN = 0.0
    RTHBLTEN = 0.0
    RQVBLTEN = 0.0
    RQCBLTEN = 0.0
    RQIBLTEN = 0.0
    TKE_MYJ = 0.0
    z_init = 0.0
    u_init = 0.0
    v_init = 0.0
    t_init = 0.0
    th_init = 0.0
    exn_init = 0.0
    q_init = 0.0
    p_init = 0.0
    rho_init = 0.0
    z8w_init = 0.0
    p8w_init = 0.0
    uflux = 0.0
    vflux = 0.0
    hflux = 0.0
    qflux = 0.0
    k_t = 0.0
    k_m = 0.0
    u_g = 0.0
    v_g = 0.0
  ENDIF

  IF (surface_physics=='SIMPLESCHEME') THEN
     num_soil_layers=1
     PRINT *,'BUCKET MODEL INACTIVE'
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
  ELSEIF (surface_physics=='FRSCHEME') THEN

     num_soil_layers=1
     PRINT *,'BUCKET MODEL ACTIVE: ',bucket_model
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

     PRINT *,'BUCKET MODEL INACTIVE'
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
     PRINT *,'BUCKET MODEL INACTIVE'
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
     PRINT *,'BUCKET MODEL INACTIVE'
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
  nsplinetimes_flux=INT((REAL(forecast_length)/REAL(splineinterval_flux))+&
       &.999)+1
  nsplinetimes_smos=INT((REAL(forecast_length)/REAL(splineinterval_smos))+&
       &.999)+1
! this is needed only for the DART interface
  n1dsplines = 2*nsplinetimes_flux + nsplinetimes_smos

   IF (init_f) THEN
      ntime=NINT(REAL(forecast_length)/dt)
   ELSE
      ntime=NINT(timetot/dt)
   ENDIF

! atmospheric profiles, forcing
! get the necessary dimensions from the input forcing file
! times are computed from the namelist
   SELECT CASE (init_f_type)

      CASE('WRF')
         CALL wrf_f_dims(ncid_f, nz_f,ns_f, nt_f, &
                        dx,cent_lat,cent_lon, &
                        truelat1,truelat2,mminlu,julday,&
                        sw_corner_lon,sw_corner_lat,projcode,&
                        lat, lon, cor)
         nt_f_flux = nt_f
         nt_f_soil = nt_f
         nt_f_smos = nt_f

      CASE('OBS')
         CALL snd_f_dims(ncid_f,ncid_soil, ncid_flux, ncid_smos, &
                        nz_f, nt_f, nt_f_flux, nt_f_smos, &
                        ns_f, nt_f_soil, &
                        mminlu,julday,&
                        lat, lon, cor)

      CASE DEFAULT
         print*,'Do not know how to initialize from type ',init_f_type
         stop 'module_wrf'
   END SELECT

   IF ( allocate_wrf ) THEN
      ALLOCATE(splinetimes(nsplinetimes),&
           &splinetimes_flux(nsplinetimes_flux),&
           &splinetimes_smos(nsplinetimes_smos),&
           u_g_f(nz,nsplinetimes),v_g_f(nz,nsplinetimes),&
         p_f(nz,nsplinetimes),p8w_f(nz,nsplinetimes),&
         t_f_uadv(nz,nsplinetimes),&
         t_f_vadv(nz,nsplinetimes),&
         q_f_uadv(nz,nsplinetimes),&
         q_f_vadv(nz,nsplinetimes),&
         u_f_uadv(nz,nsplinetimes),&
         u_f_vadv(nz,nsplinetimes),&
         v_f_uadv(nz,nsplinetimes),&
         v_f_vadv(nz,nsplinetimes),&
         u_f_wadv(nz,nsplinetimes),&
         v_f_wadv(nz,nsplinetimes),&
         t_f_wadv(nz,nsplinetimes),&
         q_f_wadv(nz,nsplinetimes),&
         glw_f(nsplinetimes_flux),gsw_f(nsplinetimes_flux),&
         precip_f(nsplinetimes_smos))
         t_f_vadv = 0.0
         u_f_vadv = 0.0
         v_f_vadv = 0.0
         q_f_vadv = 0.0
         t_f_uadv = 0.0
         u_f_uadv = 0.0
         v_f_uadv = 0.0
         q_f_uadv = 0.0
         t_f_wadv = 0.0
         u_f_wadv = 0.0
         v_f_wadv = 0.0
         q_f_wadv = 0.0
      ALLOCATE(isltyp_f(nt_f), ivgtyp_f(nt_f), lu_index_f(nt_f))
      ALLOCATE(times_f(nt_f), times_f_flux(nt_f_flux), &
            times_f_soil(nt_f_soil), times_f_smos(nt_f_smos), &
            th2_init_f(nt_f_smos), t2_init_f(nt_f_smos), &
            tsk_init_f(nt_f_smos), u10_init_f(nt_f_smos), &
            v10_init_f(nt_f_smos), q2_init_f(nt_f_smos), &
            precip_init_f(nt_f_smos), &
            glw_init_f(nt_f_flux), gsw_init_f(nt_f_flux), &
            glw_up_init_f(nt_f_flux), &
            qsfc_init_f(nt_f), tmn_init_f(nt_f), &
            vegfra_f(nt_f))
      ALLOCATE(u_init_f(nz_f,nt_f), v_init_f(nz_f,nt_f),&
            t_init_f(nz_f,nt_f), th_init_f(nz_f,nt_f), &
            exn_init_f(nz_f,nt_f),q_init_f(nz_f,nt_f), &
            p_init_f(nz_f,nt_f), z_f(nz_f,nt_f))
      ALLOCATE(z_f_stag(nz_f+1,nt_f))
      ALLOCATE(zs_f(ns_f), dzs_f(ns_f))
      ALLOCATE(tslb_init_f(ns_f,nt_f_soil), smois_init_f(ns_f,nt_f_soil))
   ENDIF

END SUBROUTINE STATIC_INIT_WRF

!*****************************************************

SUBROUTINE INIT_WRF(wrf_rnd_seed)

   INTEGER,INTENT(INOUT)        :: wrf_rnd_seed

!!! basic stuff

   a=1.-.25*dt**2*cor**2
   b=dt*cor
   c=.5*dt**2*cor**2
   d=.5*dt**2*cor
   e=1.+.25*dt**2*cor**2

   IF (.NOT.init_f) THEN

!! allocate and assign soil, surface and atmosphere for ideal begin

      time=timeo
      
      CALL init_soil_ideal(julday, lu_index, ivgtyp, &
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

!      start_seconds = 0

!! allocate and assign soil, surface and atmosphere for ideal end
   
   ELSE  ! will use WRF input instead
      
!! allocate and assign soil, surface and atmosphere for init_f begin

      idum = wrf_rnd_seed
      
!! get all we need from the input files (model)
!! these are already selected for time, but no interpolation
      SELECT CASE (init_f_type)

         CASE('WRF')
            CALL wrf_init_and_bc(ncid_f,nz_f, ns_f, nt_f, &
                    z_f,z_f_stag,t_init_f,u_init_f,v_init_f,&
                    q_init_f,p_init_f,&
                    th2_init_f,t2_init_f,tsk_init_f,&
                    u10_init_f,v10_init_f,     &
                    q2_init_f,glw_init_f,gsw_init_f,qsfc_init_f,&
                    tslb_init_f,smois_init_f,tmn_init_f,&
                    precip_init_f, &
                    vegfra_f,isltyp_f,lu_index_f,ivgtyp_f,&
                    times_f,times_f_flux,times_f_soil,times_f_smos,idum)

         CASE('OBS')
            CALL snd_init_and_bc(ncid_f, ncid_flux, ncid_soil, ncid_smos, &
                    nz_f, ns_f, nt_f, nt_f_flux, nt_f_soil, nt_f_smos,&
                    z_f,t_init_f,u_init_f,v_init_f,&
                    q_init_f,p_init_f,&
                    th2_init_f,t2_init_f,tsk_init_f,&
                    u10_init_f,v10_init_f,     &
                    q2_init_f, precip_init_f, &
                    glw_init_f,glw_up_init_f,gsw_init_f,qsfc_init_f,&
                    zs_f, &
                    tslb_init_f,smois_init_f,tmn_init_f,&
                    vegfra_f,isltyp_f,lu_index_f,ivgtyp_f,&
                    times_f,times_f_flux,times_f_soil,times_f_smos,idum)

            WHERE ( gsw_init_f < 0.0 ) gsw_init_f = 0.0
            WHERE ( glw_init_f < 0.0 ) glw_init_f = 0.0

         CASE DEFAULT
            print*,'Do not know how to initialize from type ',init_f_type
            stop 'module_wrf'

      END SELECT

      wrf_rnd_seed = idum
!      start_seconds = start_forecast
  
! initialize soil 
! this is in a separate block for now because it will undoubtedly change
      SELECT CASE (init_f_type)

         CASE('WRF')
             CALL init_soil_real_wrf(julday, lu_index, ivgtyp, &
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

         CASE('OBS')

            CALL init_soil_real_snd(julday, lu_index, ivgtyp, &
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

         CASE DEFAULT
            print*,'Do not know how to initialize from type ',init_f_type
            stop 'module_wrf'

      END SELECT

! initialize soil end
      
     DO i=1,nsplinetimes
         splinetimes(i)=times_f(1)+(i-1)*splineinterval
     ENDDO
     DO i=1,nsplinetimes_flux
         splinetimes_flux(i)=times_f(1)+(i-1)*splineinterval_flux
     ENDDO
     DO i=1,nsplinetimes_smos
         splinetimes_smos(i)=times_f(1)+(i-1)*splineinterval_smos
     ENDDO

     CALL initf(u_init_f(:,itime_f),v_init_f(:,itime_f),&
           t_init_f(:,itime_f),q_init_f(:,itime_f),&
           u_init_f(:,1:nt_f),v_init_f(:,1:nt_f),&
           glw_init_f(1:nt_f_flux),gsw_init_f(1:nt_f_flux),&
           precip_init_f(1:nt_f_smos),&
           p_init_f(:,1:nt_f),&
           th2_init_f(itime_f),q2_init_f(itime_f),&
           u10_init_f(itime_f),v10_init_f(itime_f),&
           tsk_init_f(itime_f),qsfc_init_f(itime_f),&
           z_f(:,1:nt_f),nz_f,z_g,nt_f,nt_f_flux,nt_f_smos,&
           times_f(1:nt_f),times_f_flux(1:nt_f_flux),&
           times_f_smos(1:nt_f_smos),&
           nsplinetimes,splinetimes,&
           nsplinetimes_flux,splinetimes_flux,&
           nsplinetimes_smos,splinetimes_smos,&
           pblh_ref,stepbl,lowlyr,ht,pblh,&
           th_phy,t_phy,moist,u_phy,v_phy,tke_myj,&
           u_g_f,v_g_f,glw_f,gsw_f,precip_f,p_f,p8w_f,&
           z,z8w,dz8w,&
           n_moist,P_QV,ims,ime,jms,jme,kms,kme)

     zl2=z(:,2,:)

   ENDIF ! end if init ideal or real

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
   itimestep = (dart_days * 86400 + dart_seconds) / dt + 1

!   DO itimestep=1,ntime+1

      time=start_forecast+REAL(itimestep-1)*dt
     
!   error checking to make sure DART does not over shoot the forcing
      IF ( time > start_forecast + forecast_length ) THEN
         print*,"You do not have a long enough forcing series to "
         print*,"and make sure the invervals are correct."
         stop "module_wrf"
      ENDIF
      IF (init_f) THEN

         imin = 1+INT(REAL(time-start_forecast)/REAL(splineinterval))
         imin = MIN(imin,nsplinetimes)
         imax = MIN(imin+1,nsplinetimes)
         fract = (time-splinetimes(imin))/REAL(splineinterval)

         imin_flux=1+INT(REAL(time-start_forecast)/REAL(splineinterval_flux))
         imin_flux = MIN(imin_flux,nsplinetimes_flux)
         imax_flux = MIN(imin_flux+1,nsplinetimes_flux)
         fract_flux = (time-splinetimes_flux(imin_flux))/REAL(splineinterval_flux)

         glw=glw_f(imin_flux)+(glw_f(imax_flux)-glw_f(imin_flux))*fract_flux
         gsw=gsw_f(imin_flux)+(gsw_f(imax_flux)-gsw_f(imin_flux))*fract_flux

         imin_smos=1+INT(REAL(time-start_forecast)/REAL(splineinterval_smos))
         imin_smos = min(imin_smos,nsplinetimes_smos)
         imax_smos = min(imin_smos+1,nsplinetimes_smos)

         raincv=0.
         rainncv=precip_f(imax_smos)*dt/REAL(splineinterval_smos)

         DO i=ims,ime
            DO j=jms,jme
               p_phy(i,:,j)=p_f(:,imin)+&
                    &(p_f(:,imax)-p_f(:,imin))*fract
               p8w(i,:,j)=p8w_f(:,imin)+&
                    &(p8w_f(:,imax)-p8w_f(:,imin))*fract
               PSFC(I,J)=p8w(I,kts,J)
               pi_phy(i,:,j)=(p_phy(i,:,j)/p1000mb)**rcp
            ENDDO
         ENDDO

         u_g(:)=u_g_f(:,imin)+&
              &(u_g_f(:,imax)-u_g_f(:,imin))*fract
         v_g(:)=v_g_f(:,imin)+&
              &(v_g_f(:,imax)-v_g_f(:,imin))*fract

      ELSE
         glw=glwfunc(time)
         gsw=gswfunc(time,albedo(1,1))
         raincv=0.
         rainncv=prate_ref*dt
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
         call flush(65)
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
           &,bucket_model,maxm,minm,erate&
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
         call flush(75)
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
  
  IF (init_f) THEN
     IF (MOD(time-start_forecast,outfinterval) < epsilon) THEN
        WRITE(ncunit,'(20f16.7)')time/3600.
        print*,'Writing output ',time
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

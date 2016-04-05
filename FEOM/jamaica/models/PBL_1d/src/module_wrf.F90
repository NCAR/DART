MODULE module_wrf

! Primary WRF column model driver.  All memory allocations are done in here.

  USE module_model_constants
  USE module_namelist
  USE module_ideal
  USE module_initialize
  USE module_init_soil_ideal
  USE module_init_soil_real
  USE module_init_soil_real_fluxforce
  USE module_soil_pre
  USE module_surface_driver
  USE module_pbl_driver
  USE module_sf_myjsfc
  USE module_sf_myjsfc_fluxforce
  USE module_bl_myjpbl
  USE module_sf_gfs
  USE module_sf_gfs_fluxforce
  USE module_bl_gfs
  USE module_sf_sfclay
  USE module_sf_sfclay_fluxforce
  USE module_bl_mrf
  USE module_bl_ysu
  USE module_sf_noahlsm
  USE module_sf_ruclsm
  USE module_wrf_init_and_bc
  USE module_snd_init_and_bc
  USE module_sfc_init_and_bc
  USE module_uvg_force
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
                                          nsplinetimes_advection, &
                                          nsplinetimes_flux, &
                                          nsplinetimes_smos, &
                                          nsplinetimes_sfc, &
                                          n1dsplines, &
                                          iitime,&
                                          imin,imax, &
                                          imin_adv,imax_adv,&
                                          imin_flux,imax_flux,&
                                          &imin_smos,imax_smos,&
                                          &imin_sfc,imax_sfc

   INTEGER, DIMENSION(num_tiles)       :: i_start=1,i_end=1, &
                                          j_start=1 ,j_end=1,&
                                          k_start=1,k_end
 
! forcing netCDF file IDs
   INTEGER                             :: ncid_f, ncid_soil, &
                                          ncid_flux, ncid_smos, &
                                          ncid_uvg, ncid_sfc, &
                                          ncid_eofs_forc
                                          
! base grid heights
  REAL, dimension(:), ALLOCATABLE      :: z_grid_stag, z_grid

! forcing number for soil, atmos, and times
   INTEGER                             :: ns_f, nz_f, nz_uvg, &
                                          nt_f, nt_f_flux, &
                                          nt_f_soil, nt_f_smos, nt_uvg,&
                                          &nt_f_sfc 
! DO. I added here the following
   INTEGER                             :: nz_stag_uvg,nrecords_uvg

   INTEGER                             :: itimestep,STEPBL,itime_f=1

   REAL                                :: terrain_f

   INTEGER, DIMENSION( ims:ime , jms:jme )  :: LOWLYR

   LOGICAL                             ::   warm_rain

   INTEGER :: ucmcall=0 ! all urban variables set to zero

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

   REAL, DIMENSION( ims:ime , jms:jme ):: psimfac,psihfac,psiqfac,zl2,zll

   REAL, DIMENSION( ims:ime, jms:jme ) :: TSK, &
                                          UST, uust,&
                                          HOL, &
                                          MOL, &
                                          RMOL, &
                                          PBLH, &
                                          HFX, hhfx,&
                                          QFX, &
                                          &qqfx,&
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
   
   REAL :: z_o,z_t,z_q

   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: RUBLTEN, &
                                          RVBLTEN, &
                                          RTHBLTEN, &
                                          RQVBLTEN, &
                                          RQCBLTEN, &
                                          RQIBLTEN, &
                                          TKE_MYJ,&
                                          &EXCH_H,&
                                          &EL_MYJ
  
   REAL                                ::  u_frame=0.,v_frame=0.
   
   INTEGER, DIMENSION( ims:ime , jms:jme ) :: KPBL
 
   REAL, DIMENSION( ims:ime , jms:jme ) :: landmask_input,sst_input,&
                                           cpm,chs

   LOGICAL                              :: flag_sst=.TRUE.

   REAL, DIMENSION( ims:ime , jms:jme ) :: XICE, SEAMASK, CT, SNOW, LH
   

   REAL                                 :: DTMIN,DTBL

   INTEGER                              :: i,J,K,NK,jj,ij

   REAL                                 :: dx,time,cor,&
        &fract,fract_flux,fract_sfc, fract_adv

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
                                          SMFR3D,SH2O,TRL_URB3D, &
                                          TBL_URB3D, TGL_URB3D

   REAL,  ALLOCATABLE, DIMENSION(:)    :: DZS,ZS,DZR,DZG,DZB

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

   REAL, ALLOCATABLE, DIMENSION(:)     :: u_g,v_g, &
                                          tau_u,tau_v, &
                                          th_upstream_x, th_upstream_y, &
                                          qv_upstream_x, qv_upstream_y
   
   REAL                                :: a,b,c,d,e
   REAL, PARAMETER                     :: deltat = 2.*3600

   LOGICAL, PARAMETER                  :: FNDSOILW=.TRUE.,&
                                          FNDSNOWH=.TRUE.,&
                                          restart=.FALSE.

   LOGICAL ::                             &
                                                      f_qv=.FALSE.      &
                                                     ,f_qc=.FALSE.      &
                                                     ,f_qr=.FALSE.      &
                                                     ,f_qi=.FALSE.      &
                                                     ,f_qs=.FALSE.      &
                                                     ,f_qg=.FALSE.



! variables controlling the random selection of profiles, etc
   INTEGER                             :: idum, start_seconds 

! number of available forcing times from initialization
   REAL, DIMENSION(:), ALLOCATABLE     :: times_f, times_f_flux, &
                                          times_f_soil, times_f_smos,&
                                          &times_f_sfc

! time series of profiles for forcing
   REAL, DIMENSION(:,:), ALLOCATABLE   :: u_init_f,v_init_f,&
                                          t_init_f,th_init_f,&
                                          exn_init_f,q_init_f,&
                                          p_init_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: tau_init_u, tau_init_v
   REAL, DIMENSION(:,:), ALLOCATABLE :: t_init_upstream_x, t_init_upstream_y
   REAL, DIMENSION(:,:), ALLOCATABLE :: qv_init_upstream_x, qv_init_upstream_y

! results of temporal and spatial interpolations
   REAL, ALLOCATABLE, DIMENSION(:)     :: splinetimes, splinetimes_flux,&
        &splinetimes_smos, splinetimes_sfc, splinetimes_advection
   REAL, ALLOCATABLE, DIMENSION(:,:)   :: u_g_f,v_g_f,p_f,p8w_f,&
                                          t_f_uadv,t_f_vadv,t_f_wadv,&
                                          q_f_uadv,q_f_vadv,q_f_wadv,&
                                          u_f_uadv,u_f_vadv,u_f_wadv,&
                                          v_f_uadv,v_f_vadv,v_f_wadv
   REAL, DIMENSION(:,:), ALLOCATABLE :: tau_u_f, tau_v_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: th_upstream_x_f, th_upstream_y_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: qv_upstream_x_f, qv_upstream_y_f
   
   REAL :: h_gamma
   REAL, ALLOCATABLE, DIMENSION(:)     :: glw_f,gsw_f,precip_f

   REAL, ALLOCATABLE, DIMENSION(:)     :: ts_f,qvs_f,&
        &ustar_f,hflux_f,qvflux_f


   REAL, ALLOCATABLE, DIMENSION(:)     :: ts_init_f,qvs_init_f,&
        &ustar_init_f,hflux_init_f,qvflux_init_f


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

   INTEGER :: ra_lw_physics=1,sst_update=0

! DO eofs related stuff *** Maybe this is not needed at all
  INTEGER  :: ncid_eofs_init, nz_eo, ns_eo, nt_eo

! DO I add the following to pass information from init to uvg
  INTEGER  :: itran1, itran2
  REAL     :: control_w
  REAL, DIMENSION(:), ALLOCATABLE :: gasdo
   
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
             TKE_MYJ(ims:ime, kms:kme, jms:jme),&
             &EXCH_H(ims:ime, kms:kme, jms:jme),&
             &EL_MYJ(ims:ime, kms:kme, jms:jme) )
    ALLOCATE(z_init(1:nz), u_init(1:nz), v_init(1:nz), t_init(1:nz), &
             th_init(1:nz), exn_init(1:nz), q_init(1:nz), p_init(1:nz), &
             rho_init(1:nz))
    ALLOCATE(z8w_init(1:nz+1), p8w_init(1:nz+1))
    ALLOCATE(z_grid_stag(1:nz+1), z_grid(1:nz))
    ALLOCATE(uflux(1:nz), vflux(1:nz), hflux(1:nz), &
             qflux(1:nz), k_t(1:nz),k_m(1:nz))
    ALLOCATE(u_g(1:nz), v_g(1:nz))
    
    PRINT*,'n_eo',n_eo
    ALLOCATE(gasdo(n_eo))

!   go ahead and read the grid here - redundant
    OPEN(55,file='grid_wrf1d.ascii')

    z_grid_stag=0.
    DO k=1,nz
       READ(55,*)i,z_grid_stag(k+1)
    ENDDO

    CLOSE(55)

    DO k=1,nz
       z_grid(k)=.5*(z_grid_stag(k)+z_grid_stag(k+1))
    ENDDO

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

  IF (sf_surface_physics=='FLUXFORCE') THEN
     mavail=1
     num_soil_layers=1
     PRINT *,'No sf_sfclay_physics and sf_surface_physics drivers employed'
     PRINT *,'Surface forcing provided by observed/idealized fluxes'
     ALLOCATE(&
          &TSLB( ims:ime , 1:num_soil_layers, jms:jme ),&
          &SMOIS( ims:ime , 1:num_soil_layers, jms:jme ),&
          &KEEPFR3DFLAG( ims:ime , 1:num_soil_layers, jms:jme ),&
          &SMFR3D( ims:ime , 1:num_soil_layers, jms:jme ),&
          &SH2O( ims:ime , 1:num_soil_layers, jms:jme ),&
          &DZS(1:num_soil_layers),&
          &ZS(1:num_soil_layers))
  ELSEIF (sf_surface_physics=='SIMPLESCHEME') THEN
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
  ELSEIF (sf_surface_physics=='FRSCHEME') THEN

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
  ELSEIF (sf_surface_physics=='SLABSCHEME') THEN
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
   ELSEIF (sf_surface_physics=='LSMSCHEME') THEN
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
             &ZS(1:num_soil_layers), &
             &DZR(1:num_soil_layers), &
             &DZB(1:num_soil_layers), &
             &DZG(1:num_soil_layers), &
             &TRL_URB3D( ims:ime , 1:num_soil_layers, jms:jme ),&
             &TBL_URB3D( ims:ime , 1:num_soil_layers, jms:jme ),&
             &TGL_URB3D( ims:ime , 1:num_soil_layers, jms:jme ))
     ENDIF

   ELSEIF (sf_surface_physics=='RUCLSMSCHEME') THEN
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
  
! more initialization to accommodate pgf
  TSLB = 0.0
  SMOIS = 0.0
  KEEPFR3DFLAG = 0.0
  SMFR3D = 0.0
  SH2O = 0.0
  DZS = 0.0
  ZS = 0.0

! some other stuff we'll need
  nsplinetimes=INT((REAL(forecast_length)/REAL(splineinterval))+&
       &.999)+1
  nsplinetimes_advection=INT((REAL(forecast_length)/ &
       REAL(splineinterval_advection))+.999)+1
  nsplinetimes_flux=INT((REAL(forecast_length)/REAL(splineinterval_flux))+&
       &.999)+1
  nsplinetimes_smos=INT((REAL(forecast_length)/REAL(splineinterval_smos))+&
       &.999)+1
  nsplinetimes_sfc=INT((REAL(forecast_length)/REAL(splineinterval_sfc))+&
       &.999)+1


! this is needed only for the DART interface
  n1dsplines = 2*nsplinetimes_flux + nsplinetimes_smos

   IF (init_f) THEN
      ntime=NINT(REAL(forecast_length)/dt)
   ELSE
      ntime=NINT(timetot/dt)
   ENDIF

! atmospheric profiles, forcing
! forcing will come from same files, but go and get different information
! if force_uvg = .true.
! get the necessary dimensions from the input forcing file
! times are computed from the namelist

   IF (init_f .EQV. .TRUE.) THEN
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
         nt_uvg    = nt_f
         nz_uvg    = nz_f

      CASE('OBS')
         CALL snd_f_dims(ncid_f,ncid_soil, ncid_flux, ncid_smos, &
                        nz_f, nt_f, nt_f_flux, nt_f_smos, &
                        ns_f, nt_f_soil, &
                        mminlu,julday,&
                        lat, lon, cor, terrain_f)
         nt_uvg    = nt_f
         nz_uvg    = nz_f

         IF (force_flux) THEN
            CALL sfc_f_dims(ncid_sfc, nt_f_sfc)
         ENDIF
         
      CASE DEFAULT
         print*,'Do not know how to initialize from type ',init_f_type
         stop 'module_wrf'
   END SELECT

   IF ( allocate_wrf ) THEN
      ALLOCATE(splinetimes(nsplinetimes),&
           &splinetimes_advection(nsplinetimes_advection),&
           &splinetimes_flux(nsplinetimes_flux),&
           &splinetimes_smos(nsplinetimes_smos),&
           u_g_f(nz,nsplinetimes),v_g_f(nz,nsplinetimes),&
         p_f(nz,nsplinetimes),p8w_f(nz,nsplinetimes),&
        glw_f(nsplinetimes_flux),gsw_f(nsplinetimes_flux),&
         precip_f(nsplinetimes_smos))
         p_f = 0.0
         p8w_f = 0.0
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
         ALLOCATE(t_init_upstream_x(nz_f,nt_f), &
                     t_init_upstream_y(nz_f,nt_f), &
                     qv_init_upstream_x(nz_f,nt_f), &
                     qv_init_upstream_y(nz_f,nt_f), &
                     tau_init_u(nz_f, nt_f), &
                     tau_init_v(nz_f, nt_f), &
                     th_upstream_x_f(nz,nsplinetimes_advection), &
                     th_upstream_y_f(nz,nsplinetimes_advection), &
                     qv_upstream_x_f(nz,nsplinetimes_advection), &
                     qv_upstream_y_f(nz,nsplinetimes_advection), &
                     tau_u_f(nz, nsplinetimes_advection), &
                     tau_v_f(nz, nsplinetimes_advection)) 
         ALLOCATE(tau_u(1:nz), tau_v(1:nz), &
                     th_upstream_x(1:nz),th_upstream_y(1:nz), &
                     qv_upstream_x(1:nz),qv_upstream_y(1:nz))
            tau_init_u = 0.0
            tau_init_v = 0.0
            tau_u_f = 0.0
            tau_v_f = 0.0
            tau_u = 0.0
            tau_v = 0.0
            t_init_upstream_x = 0.0
            t_init_upstream_y = 0.0
            qv_init_upstream_x = 0.0
            qv_init_upstream_y = 0.0
            th_upstream_x_f = 0.0
            th_upstream_y_f = 0.0
            th_upstream_x = 0.0
            th_upstream_y = 0.0
            qv_upstream_x_f = 0.0
            qv_upstream_y_f = 0.0
            qv_upstream_x = 0.0
            qv_upstream_y = 0.0

         ALLOCATE(z_f_stag(nz_f+1,nt_f))
         ALLOCATE(zs_f(ns_f), dzs_f(ns_f))
         ALLOCATE(tslb_init_f(ns_f,nt_f_soil), &
              &smois_init_f(ns_f,nt_f_soil))
         
         IF (force_flux) THEN
            ALLOCATE(&
                 &ts_init_f(nsplinetimes_sfc),&
                 &qvs_init_f(nsplinetimes_sfc),&
                 &ustar_init_f(nsplinetimes_sfc),&
                 &hflux_init_f(nsplinetimes_sfc),&
                 &qvflux_init_f(nsplinetimes_sfc))
            
            ALLOCATE(times_f_sfc(nt_f_sfc),&
                 &ts_f(nt_f_sfc),qvs_f(nt_f_sfc),&
                 &ustar_f(nt_f_sfc),hflux_f(nt_f_sfc),qvflux_f(nt_f_sfc))

            ALLOCATE(splinetimes_sfc(nsplinetimes_sfc))

         ELSE

! otherwise would need optional parameters in surface_driver

            ALLOCATE(&
                 &ts_init_f(1),&
                 &qvs_init_f(1),&
                 &ustar_init_f(1),&
                 &hflux_init_f(1),&
                 &qvflux_init_f(1))
            
            ALLOCATE(times_f_sfc(1),&
                 &ts_f(1),qvs_f(1),&
                 &ustar_f(1),hflux_f(1),qvflux_f(1))

            ALLOCATE(splinetimes_sfc(1))

         ENDIF
         
      ENDIF

  ENDIF
! need to get another set of time series for the geostrophic wind if using
      IF ( force_uvg ) THEN

! DO -  Add the possibility of forcing geostrophic winds with WRF or with
! RUC. RUC files structure is 'soundings like', while WRF's is 'WRF model like'
    SELECT CASE (force_f_type) 
    
      CASE('RUC')
      CALL uvg_ruc_f_dims(ncid_uvg, nz_uvg, nt_uvg)!the same subroutine that was
!  before just changed the name to be RUC specific      

      CASE('WRF')
      CALL uvg_wrf_f_dims(ncid_uvg,nz_uvg,nz_stag_uvg,nrecords_uvg,nt_uvg)!I consider that nz_uvg, nt_uvg are the same
! as for wrf, previous select case, further check *****

   END SELECT

   ENDIF

END SUBROUTINE STATIC_INIT_WRF

!*****************************************************

SUBROUTINE INIT_WRF(wrf_rnd_seed)

   INTEGER,INTENT(INOUT)        :: wrf_rnd_seed

   CHARACTER(len=10) :: sfscheme,blscheme

!!! basic stuff

   a=1.-.25*dt**2*cor**2
   b=dt*cor
   c=.5*dt**2*cor**2
   d=.5*dt**2*cor
   e=1.+.25*dt**2*cor**2

   IF ( init_f ) THEN
     OPEN(ncunit,file=out_f_file)
   ENDIF

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
   
   ELSE  ! will use WRF/OBS input instead
      
!! allocate and assign soil, surface and atmosphere for init_f begin

      idum = wrf_rnd_seed
      
!! get all we need from the input files (model)
!! these are already selected for time, but no interpolation
      SELECT CASE (init_f_type)

         CASE('WRF')
            CALL wrf_init_and_bc(ncid_f,nz_f, ns_f, nt_f, &
                    z_f,z_f_stag,t_init_f,u_init_f,v_init_f,&
                    q_init_f,p_init_f,&
                    t_init_upstream_x, t_init_upstream_y, &
                    qv_init_upstream_x, qv_init_upstream_y, &
                    tau_init_u, tau_init_v, &
                    th2_init_f,t2_init_f,tsk_init_f,&
                    u10_init_f,v10_init_f,     &
                    q2_init_f,glw_init_f,gsw_init_f,qsfc_init_f,&
                    tslb_init_f,smois_init_f,tmn_init_f,&
                    precip_init_f, &
                    vegfra_f,isltyp_f,lu_index_f,ivgtyp_f, terrain_f, dx, &
                    times_f,times_f_flux,times_f_soil,times_f_smos,idum,&
                    itran1,itran2,control_w,gasdo)

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

            IF (force_flux) THEN
               CALL sfc_init_and_bc(ncid_sfc, nt_f_sfc,&
                    &ts_init_f,qvs_init_f,&
                    &ustar_init_f,hflux_init_f,qvflux_init_f,&
                    &times_f_sfc,idum)
            ENDIF


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

            IF (.NOT. force_flux) THEN

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
               
            ELSE
               

               CALL init_soil_real_fluxforce(julday, lu_index, ivgtyp, &
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

               
            ENDIF

         CASE DEFAULT
            print*,'Do not know how to initialize from type ',init_f_type
            stop 'module_wrf'

         END SELECT

! initialize soil end
      
     DO i=1,nsplinetimes
         splinetimes(i)=times_f(1)+(i-1)*splineinterval
     ENDDO
     DO i=1,nsplinetimes_advection
         splinetimes_advection(i)=times_f(1)+(i-1)*splineinterval_advection
     ENDDO
     DO i=1,nsplinetimes_flux
         splinetimes_flux(i)=times_f(1)+(i-1)*splineinterval_flux
     ENDDO
     DO i=1,nsplinetimes_smos
         splinetimes_smos(i)=times_f(1)+(i-1)*splineinterval_smos
     ENDDO

     IF (force_flux) THEN
        DO i=1,nsplinetimes_smos
           splinetimes_sfc(i)=times_f(1)+(i-1)*splineinterval_sfc
        ENDDO
     ENDIF
     CALL initf(u_init_f(:,itime_f),v_init_f(:,itime_f),&
           t_init_f(:,itime_f),q_init_f(:,itime_f),&
           u_init_f(:,1:nt_f),v_init_f(:,1:nt_f),&
           t_init_upstream_x(:,1:nt_f), t_init_upstream_y(:,1:nt_f), &
           qv_init_upstream_x(:,1:nt_f), qv_init_upstream_y(:,1:nt_f), &
           tau_init_u(:,1:nt_f), tau_init_v(:,1:nt_f), &
           glw_init_f(1:nt_f_flux),gsw_init_f(1:nt_f_flux),&
           precip_init_f(1:nt_f_smos),&
           ts_init_f(1:nt_f_sfc),&
           qvs_init_f(1:nt_f_sfc),&
           ustar_init_f(1:nt_f_sfc),&
           hflux_init_f(1:nt_f_sfc),&
           qvflux_init_f(1:nt_f_sfc),&
           p_init_f(:,1:nt_f),&
           th2_init_f(itime_f),q2_init_f(itime_f),&
           u10_init_f(itime_f),v10_init_f(itime_f),&
           tsk_init_f(itime_f),qsfc_init_f(itime_f),&
           z_f(:,1:nt_f),nz_f,z_g,nt_f,nt_f_flux,nt_f_smos,nt_f_sfc,&
           times_f(1:nt_f),times_f_flux(1:nt_f_flux),&
           times_f_smos(1:nt_f_smos),&
           times_f_sfc(1:nt_f_sfc),&
           nsplinetimes,splinetimes,&
           nsplinetimes_advection,splinetimes_advection,&
           nsplinetimes_flux,splinetimes_flux,&
           nsplinetimes_smos,splinetimes_smos,&
           nsplinetimes_sfc,splinetimes_sfc,&
           pblh_ref,stepbl,lowlyr,ht,pblh,&
           th_phy,t_phy,moist,u_phy,v_phy,p_phy, p8w, tke_myj,&
           u_g_f,v_g_f,glw_f,gsw_f,precip_f,p_f,p8w_f,&
           th_upstream_x_f, th_upstream_y_f, &
           qv_upstream_x_f, qv_upstream_y_f, &
           tau_u_f, tau_v_f, &
           ts_f,qvs_f,ustar_f,hflux_f,qvflux_f,&
           z,z8w,dz8w,&
           ims,ime,jms,jme,kms,kme)

     zl2=z(:,2,:)

! if we are forcing with an independent UVG source, get the data and replace the
! forcing u_g_f and v_g_f

     IF ( force_uvg ) THEN

        idum = wrf_rnd_seed

      SELECT CASE (force_f_type)
       
       CASE('RUC')

        CALL uvg_ruc_f_bc(ncid_uvg, nz_uvg, nt_uvg, terrain_f, nz, z,  &
                      u_g_f, v_g_f, nsplinetimes, splinetimes, init_f_type, &
                      idum, ims,ime,jms,jme,kms,kme)

       CASE('WRF')

         CALL uvg_wrf_f_bc(ncid_uvg, nz_uvg, nt_uvg, terrain_f, nz, z,  &
                      u_g_f, v_g_f, nsplinetimes, splinetimes, &
                      glw_init_f, gsw_init_f, &
                      glw_f, gsw_f, &
                      nt_f_flux, nsplinetimes_flux, splinetimes_flux, &
                      init_f_type, &
                      idum, ims,ime,jms,jme,kms,kme,ncid_eofs_forc, &
                      itran1,itran2,control_w,gasdo,nz_stag_uvg,nrecords_uvg)

       CASE DEFAULT
        print*, 'Do not know how to force from type ',force_f_type
        stop 'module_wrf'
       
       END SELECT
        
        wrf_rnd_seed = idum

     ENDIF


  ENDIF ! end if init ideal or real

   k_t=0.
   k_m=0.
   uflux=0.
   vflux=0.
   hflux=0.
   qflux=0.


!! initialize pbl begin

   IF (sf_surface_physics=='FLUXFORCE') THEN
      sfscheme='_flux'
   ELSEIF (sf_surface_physics=='FRSCHEME') THEN
      sfscheme='_frb'
   ELSEIF (sf_surface_physics=='SLABSCHEME') THEN
      sfscheme='_slab'
   ELSEIF (sf_surface_physics=='LSMSCHEME') THEN
      sfscheme='_noah'
   ELSEIF (sf_surface_physics=='RUCLSMSCHEME') THEN
      sfscheme='_ruc'
   ELSE
      PRINT *,'Unknown surface physics option ',sfscheme
      PRINT *,'Stopping'
      STOP
   ENDIF

   IF (bl_pbl_physics=='MYJPBLSCHEME') THEN
      blscheme='myj'
   ELSEIF (bl_pbl_physics=='MRFPBLSCHEME') THEN
      blscheme='mrf'
   ELSEIF (bl_pbl_physics=='YSUPBLSCHEME') THEN
      blscheme='ysu'
   ELSEIF (bl_pbl_physics=='GFSPBLSCHEME') THEN
      blscheme='gfs'
   ELSE
      PRINT *,'Unknown bouldary layer option ',blscheme
      PRINT *,'Stopping'
      STOP
   ENDIF
   
   OPEN(isfcunit ,file=TRIM(outdir)//'/'//&
        &TRIM(blscheme)//TRIM(sfscheme)//'_sfc.txt')
   OPEN(iprofunit,file=TRIM(outdir)//'/'//&
        &TRIM(blscheme)//TRIM(sfscheme)//'_prof.txt')
   OPEN(isimil,  file=TRIM(outdir)//'/'//&
        &TRIM(blscheme)//TRIM(sfscheme)//'_simil.txt')

   IF (bl_pbl_physics=='MYJPBLSCHEME') THEN


      CALL MYJSFCINIT(LOWLYR,UST,Znt            &
           ,SEAMASK,XICE,IVGTYP,.FALSE.                &
           ,.TRUE.                           &
           ,IDS,IDE,JDS,JDE,KDS,KDE                    &
           ,IMS,IME,JMS,JME,KMS,KME                    &
           ,ITS,ITE,JTS,JTE,KTS,KTE)
      
      CALL myjpblinit(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN, &
           TKE_MYJ,EXCH_H,.FALSE.,.TRUE.,                      &
           ids, ide, jds, jde, kds, kde,         &
           ims, ime, jms, jme, kms, kme,         &
           its, ite, jts, jte, kts, kte          )
     
   ELSEIF (bl_pbl_physics=='YSUPBLSCHEME') THEN

      CALL sfclayinit(.TRUE.)
      CALL ysuinit(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,    &
           RQCBLTEN,RQIBLTEN,P_QI,               &
           PARAM_FIRST_SCALAR,                   &
           .FALSE.,.TRUE.,                              &
           ids, ide, jds, jde, kds, kde,         &
           ims, ime, jms, jme, kms, kme,         &
           its, ite, jts, jte, kts, kte          )

   ELSEIF (bl_pbl_physics=='MRFPBLSCHEME') THEN

      CALL sfclayinit(.TRUE.)
      CALL mrfinit(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,    &
           RQCBLTEN,RQIBLTEN,P_QI,               &
           PARAM_FIRST_SCALAR,                   &
           .FALSE.,.TRUE.,                       &
           ids, ide, jds, jde, kds, kde,         &
           ims, ime, jms, jme, kms, kme,         &
           its, ite, jts, jte, kts, kte          )

   ELSEIF (bl_pbl_physics=='GFSPBLSCHEME') THEN

      CALL gfsinit(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,    &
           RQCBLTEN,RQIBLTEN,P_QI,               &
           PARAM_FIRST_SCALAR,                   &
           .FALSE.,.TRUE.,                       &
           ids, ide, jds, jde, kds, kde,         &
           ims, ime, jms, jme, kms, kme,         &
           its, ite, jts, jte, kts, kte          )
   ENDIF


END SUBROUTINE INIT_WRF

!------------------------------------------------------------------
SUBROUTINE WRF(seconds,days)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: days, seconds
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
   zll=z(:,1,:)


!  here we compute the time step to keep pace with dart
   itimestep = (days * 86400 + seconds) / dt + 1

!   DO itimestep=1,ntime+1

      time=start_forecast+REAL(itimestep-1)*dt
     
      IF (init_f) THEN

!   error checking to make sure DART does not over shoot the forcing
         IF ( time > start_forecast + forecast_length ) THEN
            PRINT*,"You do not have a long enough forcing series to "
            PRINT*,"and make sure the invervals are correct."
            STOP "module_wrf"
         ENDIF
         
         imin = 1+INT(REAL(time-start_forecast)/REAL(splineinterval))
         imin = MIN(imin,nsplinetimes)
         imax = MIN(imin+1,nsplinetimes)
         fract = (time-splinetimes(imin))/REAL(splineinterval)

         IF ( t_advection ) THEN
            imin_adv=1+INT(REAL(time-start_forecast)/REAL(splineinterval_advection))
            imin_adv = MIN(imin_adv,nsplinetimes_advection)
            imax_adv = MIN(imin_adv+1,nsplinetimes_advection)
            fract_adv = (time-splinetimes_advection(imin_adv))/REAL(splineinterval_advection)
         ENDIF

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
!               p_phy(i,:,j)=p_f(:,imin)+&
!                    &(p_f(:,imax)-p_f(:,imin))*fract
!               p8w(i,:,j)=p8w_f(:,imin)+&
!                    &(p8w_f(:,imax)-p8w_f(:,imin))*fract
               PSFC(I,J)=p8w(I,kts,J)
               pi_phy(i,:,j)=(p_phy(i,:,j)/p1000mb)**rcp
            ENDDO
         ENDDO

         u_g(:)=u_g_f(:,imin)+&
              &(u_g_f(:,imax)-u_g_f(:,imin))*fract
         v_g(:)=v_g_f(:,imin)+&
              &(v_g_f(:,imax)-v_g_f(:,imin))*fract

         IF ( t_advection .or. qv_advection ) THEN
            tau_u(:)=tau_u_f(:,imin_adv)+&
                 &(tau_u_f(:,imax_adv)-tau_u_f(:,imin_adv))*fract_adv
            tau_v(:)=tau_v_f(:,imin_adv)+&
                 &(tau_v_f(:,imax_adv)-tau_v_f(:,imin_adv))*fract_adv
            IF ( t_advection ) THEN
               th_upstream_x(:)=th_upstream_x_f(:,imin_adv)+&
                 &(th_upstream_x_f(:,imax_adv)-th_upstream_x_f(:,imin_adv))*fract_adv
               th_upstream_y(:)=th_upstream_y_f(:,imin_adv)+&
                 &(th_upstream_y_f(:,imax_adv)-th_upstream_y_f(:,imin_adv))*fract_adv

            ENDIF
            IF ( qv_advection ) THEN
               qv_upstream_x(:)=qv_upstream_x_f(:,imin_adv)+&
                 &(qv_upstream_x_f(:,imax_adv)-qv_upstream_x_f(:,imin_adv))*fract_adv
               qv_upstream_y(:)=qv_upstream_y_f(:,imin_adv)+&
                 &(qv_upstream_y_f(:,imax_adv)-qv_upstream_y_f(:,imin_adv))*fract_adv
            ENDIF
         ENDIF

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
         WRITE(65,*)'l1 ',sf_sfclay_physics,sf_surface_physics,bl_pbl_physics
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

      IF (force_flux) THEN

         uust=ust
         hhfx=hfx
         qqfx=qfx

         IF (init_f) THEN
            imin_sfc=1+INT(REAL(time-start_forecast)/&
                 &REAL(splineinterval_sfc))
            
            imin_sfc = MIN(imin_sfc,nsplinetimes_sfc)
            imax_sfc = MIN(imin_sfc+1,nsplinetimes_sfc)
            
            fract_sfc = (time-splinetimes_sfc(imin_sfc))/&
                 &REAL(splineinterval_sfc)

            IF (itimestep==1) THEN
               hfx=hflux_f(imin_sfc)
               qfx=qvflux_f(imin_sfc)
               ust=ustar_f(imin_sfc)
               tsk=ts_f(imin_sfc)
               qsfc=qvs_f(imin_sfc)
            ELSE

               IF (.NOT. calc_ust) THEN
                  ust=ustar_f(imin_sfc)+&
                       &(ustar_f(imax_sfc)-ustar_f(imin_sfc))*fract_sfc
               ENDIF

               IF (calc_sfc) THEN
                  IF (xland(1,1) < 1.5) THEN
                     CALL calcsfc(hfx(1,1),flhc(1,1),qfx(1,1),flqc(1,1),&
                          &th_phy(1,1,1),moist(1,1,1,P_QV),&
                          &psfc(1,1),&
                          &tsk(1,1),qsfc(1,1),bl_pbl_physics)
                  ENDIF
               ELSE
                  tsk=ts_f(imin_sfc)+(ts_f(imax_sfc)-ts_f(imin_sfc))*&
                       &fract_sfc
                  
                  qsfc=qvs_f(imin_sfc)+(qvs_f(imax_sfc)-qvs_f(imin_sfc))*&
                       &fract_sfc
               ENDIF
         

               hfx=hflux_f(imin_sfc)+&
                    &(hflux_f(imax_sfc)-hflux_f(imin_sfc))*&
                    &fract_sfc
               qfx=qvflux_f(imin_sfc)+&
                    &(qvflux_f(imax_sfc)-qvflux_f(imin_sfc))*fract_sfc


            ENDIF

            
         ELSE !ideal

            calc_ust=.TRUE.

            IF (calc_sfc) THEN
! for idealized case only constant hflux=hfx/(rho*cp) and 
! qvflux=qfx/(rho) possible
               
               psfc=p8w(:,kts,:)
               rhosfc=psfc(1,1)/(r_d*tsk(1,1)*(1.+svp1*qsfc(1,1)))
               hfx(1,1)=hflux_ref/rhosfc
               qfx(1,1)=qvflux_ref/(rhosfc*xlv)

               IF (xland(1,1) < 1.5) THEN
                  CALL calcsfc(hfx(1,1),flhc(1,1),qfx(1,1),flqc(1,1),&
                       &th_phy(1,1,1),moist(1,1,1,P_QV),&
                       &psfc(1,1),&
                       &tsk(1,1),qsfc(1,1),bl_pbl_physics)
               ENDIF

            ELSE
               PRINT*,'For the idealized case tsk and qsfc need to be &
                    &recalculated to avoid inconsistencies'
               PRINT *,'In namelist change calc_sfc -> .TRUE.'
               PRINT *,'Stopping in module_wrf.F'
               STOP
            ENDIF

         ENDIF ! split real/idealized fluxforce
   
      ENDIF ! end fluxforce split


      CALL surface_driver(                                                &
     &         ACSNOM=acsnom      ,ACSNOW=acsnow      ,AKHS=akhs          &
     &        ,AKMS=akms          ,ALBBCK=albbck      ,ALBEDO=albedo      &
     &        ,BR=br              ,CANWAT=canwat      ,CHKLOWQ=chklowq    &
     &        ,CT=ct              ,DT=dt              ,DX=dx              &
     &        ,DZ8W=dz8w          ,DZS=dzs            ,FLHC=flhc          &
     &        ,FLQC=flqc          ,GLW=glw            ,GRDFLX=grdflx      &
     &        ,GSW=gsw            ,GZ1OZ0=gz1oz0      ,HFX=hfx            &
     &        ,HT=ht              ,IFSNOW=ifsnow      ,ISFFLX=isfflx      &
     &        ,ISLTYP=isltyp      ,ITIMESTEP=itimestep                    &
     &        ,IVGTYP=ivgtyp      ,LH=lh              ,LOWLYR=lowlyr      &
     &        ,MAVAIL=mavail      ,NUM_SOIL_LAYERS=num_soil_layers        &
     &        ,P8W=p8w            ,PBLH=pblh          ,PI_PHY=pi_phy      &
     &        ,PSFC=psfc          ,PSHLTR=pshltr      ,PSIH=psih          &
     &        ,PSIM=psim          ,P_PHY=p_phy        ,Q10=q10            &
     &        ,Q2=q2              ,QFX=qfx            ,QSFC=qsfc          &
     &        ,QSHLTR=qshltr      ,QZ0=qz0            ,RAINCV=raincv      &
     &        ,RA_LW_PHYSICS=ra_lw_physics            ,RHO=rho            &
     &        ,RMOL=rmol          ,SFCEVP=sfcevp      ,SFCEXC=sfcexc      &
     &        ,SFCRUNOFF=sfcrunoff                                        &
     &         ,BL_PBL_PHYSICS=bl_pbl_physics&
     &        ,SF_SFCLAY_PHYSICS=sf_sfclay_physics                        &
     &        ,SF_SURFACE_PHYSICS=sf_surface_physics  ,SH2O=sh2o          &
     &        ,SHDMAX=shdmax      ,SHDMIN=shdmin      ,SMOIS=smois        &
     &        ,SMSTAV=smstav      ,SMSTOT=smstot      ,SNOALB=snoalb      &
     &        ,SNOW=snow          ,SNOWC=snowc        ,SNOWH=snowh        &
     &        ,SST=sst_input         ,SST_UPDATE=sst_update                  &
     &        ,STEPBL=stepbl      ,TH10=th10          ,TH2=th2            &
     &        ,THZ0=thz0          ,TH_PHY=th_phy      ,TKE_MYJ=tke_myj    &
     &        ,TMN=tmn            ,TSHLTR=tshltr      ,TSK=tsk            &
     &        ,TSLB=tslb          ,T_PHY=t_phy        ,U10=u10            &
     &        ,UDRUNOFF=udrunoff  ,UST=ust            ,UZ0=uz0            &
     &        ,U_FRAME=u_frame    ,U_PHY=u_phy        ,V10=v10            &
     &        ,VEGFRA=vegfra      ,VZ0=vz0            ,V_FRAME=v_frame    &
     &        ,V_PHY=v_phy        ,WARM_RAIN=warm_rain                    &
     &        ,WSPD=wspd          ,XICE=xice          ,XLAND=xland        &
     &        ,Z0=z0              ,Z=z                ,ZNT=znt            &
     &        ,ZS=zs                                                      &
     &        ,DECLIN_URB=0.  ,COSZ_URB2D=(/0.,0./) & !I urban
     &        ,OMG_URB2D=(/0.,0./)    ,xlat_urb2d=(/0.,0/)        & !I urban
     &        ,NUM_ROOF_LAYERS=0                          & !I urban
     &        ,NUM_WALL_LAYERS=0                          & !I urban
     &        ,NUM_ROAD_LAYERS=0                          &
     &        ,DZR=DZR ,DZB=DZB ,DZG=DZG                 & !I urban
     &        ,TR_URB2D=(/0.,0./) ,TB_URB2D=(/0.,0./)           &
     &        ,TG_URB2D=(/0.,0./)                                   & !H urban
     &        ,TC_URB2D=(/0.,0./) ,QC_URB2D=(/0.,0./)           & !H urban
     &        ,UC_URB2D=(/0.,0./)                                 & !H urban
     &        ,XXXR_URB2D=(/0.,0./)                               &
     &        ,XXXB_URB2D=(/0.,0./)                               & !H urban
     &        ,XXXG_URB2D=(/0.,0./)                               &
     &        ,XXXC_URB2D=(/0.,0./)                               & !H urban
     &        ,TRL_URB3D=TRL_URB3D   ,TBL_URB3D=TBL_URB3D     & !H urban
     &        ,TGL_URB3D=TBL_URB3D                                 & !H urban
     &        ,SH_URB2D=(/0.,0./)     ,LH_URB2D=(/0.,0./)       &
     &        ,G_URB2D=(/0.,0./)                                     & !H urban
     &        ,RN_URB2D=(/0.,0./)     , TS_URB2D=(/0.,0./)      & !H urban
     &        ,FRC_URB2D=(/0.,0./)                                 & !H urban
     &        ,UTYPE_URB2D=(/0,0/)                             & !H urban
     &        ,ucmcall=ucmcall                                     & !H urban

           ! Indexes
     &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde          &
     &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme          &
     &        , I_START=i_start,I_END=min(i_end, ide)         &
     &        , J_START=j_start,J_END=min(j_end, jde)         &
     &        , KTS=k_start(1), KTE=kte                         &
     &        , NUM_TILES=num_tiles                                  &
           ! Optional
     &        ,QV_CURR=moist(ims,kms,jms,P_QV) ,F_QV=F_QV                 &
     &        ,QC_CURR=moist(ims,kms,jms,P_QC), F_QC=F_QC                 &
     &        ,QR_CURR=moist(ims,kms,jms,P_QR), F_QR=F_QR                 &
     &        ,QI_CURR=moist(ims,kms,jms,P_QI), F_QI=F_QI                 &
     &        ,QS_CURR=moist(ims,kms,jms,P_QS), F_QS=F_QS                 &
     &        ,QG_CURR=moist(ims,kms,jms,P_QG), F_QG=F_QG                 &
     &        ,CAPG=capg, EMISS=emiss, HOL=hol,MOL=mol                    &
     &        ,RAINBL=rainbl                                              &
     &        ,RAINNCV=rainncv,REGIME=regime,T2=t2,THC=thc                &
     &        ,QSG=qsg,QVG=qvg,QCG=qcg,SOILT1=soilt1,TSNAV=tsnav          & ! ruc lsm
     &        ,SMFR3D=smfr3d,KEEPFR3DFLAG=keepfr3dflag,                    & ! ruc lsm

     &ts_ref=ts_ref,ps_ref=ts_ref,dtamplitude_ref=dtamplitude_ref,&
     &mavail_ref=mavail_ref,time=time,&
     &bucket_model=bucket_model,maxm=maxm,minm=minm,erate=erate,&
     &calc_ust=calc_ust,z_o=z_o,z_t=z_t,z_q=z_q&
     &                                                              )

      IF ( debug2 ) then
         write(75,*)'AFTER SFC'
         write(75,*)'rad ',gsw,glw
         write(75,*)'u,v ',u_g(1),v_g(1)
         write(75,*)'p,p8 ',p_phy(1,1,1),p8w(1,1,1)
         WRITE(75,*)'l1 ',sf_sfclay_physics,sf_surface_physics,bl_pbl_physics
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

      IF (debug3 ) THEN
         write(66,*)'BEFORE PBL'
         write(66,*)'l1', bl_pbl_physics,sf_surface_physics
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
         WRITE(66,*)'l16',P_QI,P_QV,P_QC,PARAM_FIRST_SCALAR
         WRITE(66,*)'l17',uflux(1),vflux(1),hflux(1),qflux(1), &
              &k_t(1),k_m(1)
      ENDIF


      IF (sf_surface_physics=='FLUXFORCE') THEN
         qvg=qsfc
      ELSEIF (sf_surface_physics=='SIMPLESCHEME') THEN
         qvg=qsfc
      ELSEIF (sf_surface_physics=='SLABSCHEME') THEN
         qvg=qsfc/(1-qsfc)
      ELSEIF (sf_surface_physics=='FRSCHEME') THEN
         qvg=qsfc/(1-qsfc)
      ELSEIF (sf_surface_physics=='LSMSCHEME') THEN
         qvg=qsfc/(1-qsfc)
      ENDIF

!      debug4=.TRUE.


      CALL pbl_driver(                                                    &
     &         AKHS=akhs          ,AKMS=akms                              &
     &        ,BL_PBL_PHYSICS=bl_pbl_physics                 &
     &,sf_surface_physics=sf_surface_physics &
     &        ,BR=br              ,CHKLOWQ=chklowq    ,CT=ct              &
     &        ,DT=dt              ,DX=dx              ,DZ8W=dz8w          &
     &        ,EL_MYJ=el_myj      ,EXCH_H=exch_h      ,GRDFLX=grdflx      &
     &        ,GZ1OZ0=gz1oz0      ,HFX=hfx            ,HT=ht              &
     &        ,ITIMESTEP=itimestep                    ,KPBL=kpbl          &
     &        ,LH=lh              ,LOWLYR=lowlyr      ,P8W=p8w            &
     &        ,PBLH=pblh          ,PI_PHY=pi_phy      ,PSIH=psih          &
     &        ,PSIM=psim          ,P_PHY=p_phy        ,QFX=qfx            &
     &        ,QSFC=qsfc          ,QZ0=qz0                                &
     &        ,RA_LW_PHYSICS=ra_lw_physics                   &
     &        ,RHO=rho            ,RQCBLTEN=rqcblten  ,RQIBLTEN=rqiblten  &
     &        ,RQVBLTEN=rqvblten  ,RTHBLTEN=rthblten  ,RUBLTEN=rublten    &
     &        ,RVBLTEN=rvblten    ,SNOW=snow          ,STEPBL=stepbl      &
     &        ,THZ0=thz0          ,TH_PHY=th_phy      ,TKE_MYJ=tke_myj    &
     &        ,TSK=tsk            ,T_PHY=t_phy        ,UST=ust            &
     &        ,UZ0=uz0            ,U_FRAME=u_frame    ,U_PHY=u_phy        &
     &        ,VZ0=vz0            ,V_FRAME=v_frame    ,V_PHY=v_phy        &
     &        ,WARM_RAIN=warm_rain                    ,WSPD=wspd          &
     &        ,XICE=xice          ,XLAND=xland        ,Z=z                &
     &        ,ZNT=znt                                                    &
     &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde          &
     &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme          &
     &        ,I_START=i_start,I_END=min(i_end, ide)          &
     &        ,J_START=j_start,J_END=min(j_end, jde)          &
     &        ,KTS=k_start(1), KTE=kde-1                          &
     &        ,NUM_TILES=num_tiles                                   &
          ! optional
     &        ,QV_CURR=moist(ims,kms,jms,P_QV), F_QV=F_QV                 &
     &        ,QC_CURR=moist(ims,kms,jms,P_QC), F_QC=F_QC                 &
     &        ,QR_CURR=moist(ims,kms,jms,P_QR), F_QR=F_QR                 &
     &        ,QI_CURR=moist(ims,kms,jms,P_QI), F_QI=F_QI                 &
     &        ,QS_CURR=moist(ims,kms,jms,P_QS), F_QS=F_QS                 &
     &        ,QG_CURR=moist(ims,kms,jms,P_QG), F_QG=F_QG                 &
     &        ,HOL=HOL, MOL=MOL, REGIME=REGIME                            &
           &,uflux=uflux,vflux=vflux,hflux=hflux,qflux=qflux, &
           &k_t=k_t,k_m=k_m&
           ,calc_ust=calc_ust&
     &                                                          )


      IF (debug4 ) then
         write(76,*)'AFTER PBL'
         write(76,*)'l1', bl_pbl_physics,sf_surface_physics
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

! dynamics time step here.
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
               rho(i,k,j)=p_phy(i,k,j)/(r_d*&
                    &t_phy(i,k,j)*(1.+SVP1*moist(i,k,j,P_QV)))
            ENDDO
         ENDDO
      ENDDO
      
! nudge advective tendencies if requested (follow Ghan et al 1999)
      IF ( t_advection ) THEN
      DO i=ims,ime
         DO j=jms,jme
            DO k=kts,kte
               IF ( tau_u(k) > -9998.0 .and. tau_v(k) > -9998.0 .and. &
                  th_upstream_x(k) > -9998.0 .and. th_upstream_y(k) > -9998.0)&
                  THEN
                  th_phy(i,k,j) = th_phy(i,k,j) &
                          +  dt*(th_upstream_x(k) - th_phy(i,k,j))/tau_u(k) &
                          +  dt*(th_upstream_y(k) - th_phy(i,k,j))/tau_v(k) 
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ENDIF

      IF ( qv_advection ) THEN
      DO i=ims,ime
         DO j=jms,jme
            DO k=kts,kte
               IF ( tau_u(k) > -9998.0 .and. tau_v(k) > -9998.0 .and. &
                  qv_upstream_x(k) > -9998.0 .and. qv_upstream_y(k) > -9998.0)&
                  THEN
                  moist(i,k,j,P_QV) = moist(i,k,j,P_QV) &
                         +  dt*(qv_upstream_x(k) - moist(i,k,j,P_QV))/tau_u(k)&
                         +  dt*(qv_upstream_y(k) - moist(i,k,j,P_QV))/tau_v(k)
                         
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ENDIF

! update t
      t_phy=th_phy*pi_phy

      IF (debug5 ) then
         write(85,*)'AFTER PBL and TENDS'
         write(85,*)'COEFFS ',a,b,c,d
         write(85,*)'l1', bl_pbl_physics,sf_surface_physics
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

!  IF (MOD(NINT(time-timeo),outinterval)==0) THEN
!     WRITE(isfcunit,'(30f16.7)')time/3600.,tsk(1,1),tslb(1,:,1),&
!          &tmn(1,1),qsfc(1,1),SMOIS(1,:,1),gsw(1,1),ust,hfx,xlv*qfx
!
!      rhosfc=psfc(1,1)/(r_d*tsk(1,1)*(1.+svp1*qsfc(1,1)))
!
!     CLOSE(iprofunit)
!
!     DO k=kts,kte-1
!        WRITE(iprofunit,'(20f16.7)')z(1,k,1),z8w(1,k,1),u_phy(1,k,1),&
!             &v_phy(1,k,1),SQRT(u_phy(1,k,1)**2+v_phy(1,k,1)**2),&
!             &th_phy(1,k,1),moist(1,k,1,P_QV),tke_myj(1,k,1),&
!             &k_t(k),k_m(k),z8w(1,k,1)*cor/ust(1,1),&
!             &-uflux(k)/ust**2,vflux(k)/ust**2,hflux(k)/hfx*cp*rhosfc,&
!             &SQRT(uflux(k)**2+vflux(k)**2)/ust**2,qflux(k)/qfx*rhosfc
!     ENDDO
!
!
!  ENDIF


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
     WRITE(isfcunit,'(30f16.7)')time/3600.,tsk(1,1),tmn(1,1),qsfc(1,1),&
          &ust,hfx,xlv*qfx,gsw(1,1),glw(1,1),tslb(1,:,1),SMOIS(1,:,1)
  ENDIF


  IF (itimestep==ntime+1) THEN

     CLOSE(isfcunit)
     CLOSE(ncunit)

     rhosfc=psfc(1,1)/(r_d*tsk(1,1)*(1.+svp1*qsfc(1,1)))

     IF (qfx(1,1) < epsilon) qfx=epsilon

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
        
        u_simil(k)=ust(1,1)/karman*(LOG(z_level/z_o)-psimzl)*&
             &u_phy(1,1,1)/&
             &SQRT(u_phy(1,1,1)**2+v_phy(1,1,1)**2)
        v_simil(k)=ust(1,1)/karman*(LOG(z_level/z_o)-psimzl)*&
             &v_phy(1,1,1)/&
             &SQRT(u_phy(1,1,1)**2+v_phy(1,1,1)**2)
        q_simil(k)=qstar/karman*(LOG(z_level/z_q)-psihzl)+qsfc(1,1)
        th_simil(k)=(thvstar/karman*(LOG(z_level/z_t)-psihzl)+&
             &tsk(1,1)*(1.+SVP1*qsfc(1,1)))/(1.+SVP1*moist(1,k,1,P_QV))
        
        WRITE(isimil,'(20f16.7)')&
             
             &LOG(z_level/z_o),&
             
             &LOG(z_level/z_t),&
             
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
!          

             &(q_simil(k)-qsfc(1,1))/qstar,&
             &(moist(1,k,1,P_QV)-qsfc(1,1))/qstar !,&
        

!           &(q_simil(k)-qsfc(1,1))/(moist(1,k,1,P_QV)-qsfc(1,1))
        
        
     ENDDO

     CLOSE(isimil)
     
  ENDIF
  
END SUBROUTINE output_wrf_profiles


END MODULE module_wrf

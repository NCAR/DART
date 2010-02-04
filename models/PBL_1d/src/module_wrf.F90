MODULE module_wrf
!
! DART $Id$
!
! Primary WRF column model driver.  All memory allocations are done in here.

  USE module_model_constants
  USE module_namelist
  USE module_ideal
  USE module_initialize
  USE module_init_soil_ideal
  USE module_init_soil_real
  USE module_init_soil_real_fluxforce
  USE module_soil_pre
  USE module_stochastic_cloud ! added...RT
  USE module_radiation_driver
  USE module_ra_sw
  USE module_ra_rrtm
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
  USE module_microphysics_driver
  USE module_mp_etanew
  USE module_mp_kessler
  USE module_mp_lin
  USE module_mp_ncloud3
  USE module_mp_ncloud5
  USE module_mp_thompson
  USE module_mp_wsm3
  USE module_mp_wsm5
  USE module_mp_wsm6
  USE module_wrf_init_and_bc
  USE module_snd_init_and_bc
  USE module_sfc_init_and_bc
  USE module_uvg_force
  USE module_luse_init
  USE module_getsm
  USE module_interpolations ! added...RT
  USE map_utils, only : proj_info, map_init, map_set, latlon_to_ij, &
                        PROJ_LATLON, PROJ_MERC, PROJ_LC, PROJ_PS,   &
                        gridwind_to_truewind

  USE module_nr_procedures    ! added ...RT
  
  IMPLICIT NONE
 
   INTEGER                             :: isfcunit=55,&
                                          iprofunit=56,&
                                          isimil=57,&
                                          ncunit=58,&
                                          solattunit=59,& ! added...RT
                                          stochaunit=60  ! added...RT

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
                                          nsplinetimes_soil, &
                                          n1dsplines, &
                                          iitime,&
                                          imin,imax, &
                                          imin_adv,imax_adv,&
                                          imin_flux,imax_flux,&
                                          imin_smos,imax_smos,&
                                          imin_sfc,imax_sfc, &
                                          imin_soil,imax_soil

   INTEGER, DIMENSION(num_tiles)       :: i_start=1,i_end=1, &
                                          j_start=1 ,j_end=1,&
                                          k_start=1
 
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
                                          nt_f_sfc 
! DO. I added here the following
   INTEGER                             :: nz_stag_uvg,nrecords_uvg

   INTEGER                             :: itimestep,STEPBL,itime_f=1

   REAL                                :: terrain_f

   INTEGER, DIMENSION( ims:ime , jms:jme )  :: LOWLYR

   LOGICAL                             ::   warm_rain = .false.!(T for Kessler)

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
                                           th_phy, &
                                           t8w !added for call to radiation...RT

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
                                          EXCH_H,&
                                          EL_MYJ
  
   REAL                                ::  u_frame=0.,v_frame=0.
   
   INTEGER, DIMENSION( ims:ime , jms:jme ) :: KPBL
 
   REAL, DIMENSION( ims:ime , jms:jme ) :: landmask_input,sst_input,&
                                           cpm,chs

   LOGICAL                              :: flag_sst=.TRUE.

   REAL, DIMENSION( ims:ime , jms:jme ) :: XICE, SEAMASK, CT, SNOW, LH
   

   REAL                                 :: DTMIN,DTBL

   INTEGER                              :: i,J,K,NK,jj,ij

   REAL                                 :: dx,time,cor,&
                    fract,fract_flux,fract_sfc, fract_adv, fract_soil

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
                                           RAINC, RAINNC,        &
                                           SFCEVP,SFCEXC,THC,cs,TMN,&
                                           SFCRUNOFF,SNOWH,POTEVP,&
                                           SNOPCX,SOILTB,SOILT1,&
                                           TSNAV,QSG,QVG,QCG,FLHC,&
                                           FLQC,SNOALB,SMSTAV,&
                                           SMSTOT,TSHLTR,UDRUNOFF,&
                                           VEGFRA,ALBBCK,SHDMAX,SHDMIN,&
                                           Maxm,Minm,Erate

   REAL, DIMENSION( ims:ime, jms:jme )  :: TADVECTION_SCALE
   REAL, DIMENSION( ims:ime, jms:jme )  :: QADVECTION_SCALE
   REAL, DIMENSION( ims:ime, jms:jme )  :: UADVECTION_SCALE
   REAL, DIMENSION( ims:ime, jms:jme )  :: VADVECTION_SCALE
   
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

   REAL, ALLOCATABLE, DIMENSION(:)     :: uflux,vflux,hflux,qflux,&
                                          k_t,k_m

   REAL, ALLOCATABLE, DIMENSION(:,:)   :: u_f,v_f,p_f,t_f,q_f,p8w_f
   REAL, ALLOCATABLE, DIMENSION(:)     :: u_fi,v_fi,p_fi,t_fi,q_fi
   REAL, ALLOCATABLE, DIMENSION(:)     :: nudge_func
   REAL, ALLOCATABLE, DIMENSION(:)     :: u_g,v_g, &
                                          tau_u,tau_v, &
                                          th_upstream_x, th_upstream_y, &
                                          qv_upstream_x, qv_upstream_y, &
                                          qc_upstream_x, qc_upstream_y, &
                                          qr_upstream_x, qr_upstream_y, &
                                          qi_upstream_x, qi_upstream_y, &
                                          qg_upstream_x, qg_upstream_y, &
                                          u_upstream_x, u_upstream_y,   &
                                          v_upstream_x, v_upstream_y
  
   REAL                                :: a,b,c,d,e
   REAL, PARAMETER                     :: deltat = 2.*3600

   LOGICAL, PARAMETER                  :: FNDSOILW=.TRUE.,&
                                          FNDSNOWH=.TRUE.,&
                                          restart=.FALSE.

   LOGICAL ::           F_QV=.FALSE.,     & !default to F
                        F_QC=.FALSE.,     & 
                        F_QR=.FALSE.,     &
                        F_QI=.FALSE.,     &
                        F_QS=.FALSE.,     &
                        F_QG=.FALSE.


! additional variables for microphysics interface
   INTEGER, parameter :: spec_zone = 0
   INTEGER, parameter :: P_QT=1, P_QNI=1
   LOGICAL, parameter :: specified_bdy = .FALSE.
   LOGICAL, parameter :: channel_bdy   = .FALSE.
   LOGICAL ::           F_QT=.FALSE.,     &
                        F_QNI=.FALSE.
   REAL, ALLOCATABLE, DIMENSION(:)     :: MP_RESTART_STATE, &
                                         &TBPVS_STATE,      &
                                         &TBPVS0_STATE
   REAL, ALLOCATABLE, DIMENSION(:,:)   :: SR
   REAL, ALLOCATABLE, DIMENSION(:,:)   :: SNOWNC, &
                                         &SNOWNCV, &
                                         &GRAUPELNC, &
                                         &GRAUPELNCV
   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: QT_CURR, &
                                         &QNI_CURR
   REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SCALAR

   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: w_phy ! not used in these dynamics

! variables controlling the random selection of profiles, etc
   INTEGER                             :: idum, start_seconds 

! number of available forcing times from initialization
   REAL, DIMENSION(:), ALLOCATABLE     :: times_f, times_f_flux, &
                                          times_f_soil, times_f_smos,&
                                          times_f_sfc

! time series of profiles for forcing
   REAL, DIMENSION(:,:), ALLOCATABLE   :: u_init_f,v_init_f,&
                                          t_init_f,th_init_f,&
                                          exn_init_f, p_init_f,&
                                          qv_init_f,qc_init_f, &
                                          qr_init_f,qi_init_f, &
                                          qg_init_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: tau_init_u, tau_init_v
   REAL, DIMENSION(:,:), ALLOCATABLE :: t_init_upstream_x, t_init_upstream_y
   REAL, DIMENSION(:,:), ALLOCATABLE :: qv_init_upstream_x, qv_init_upstream_y
   REAL, DIMENSION(:,:), ALLOCATABLE :: qc_init_upstream_x, qc_init_upstream_y
   REAL, DIMENSION(:,:), ALLOCATABLE :: qr_init_upstream_x, qr_init_upstream_y
   REAL, DIMENSION(:,:), ALLOCATABLE :: qi_init_upstream_x, qi_init_upstream_y
   REAL, DIMENSION(:,:), ALLOCATABLE :: qg_init_upstream_x, qg_init_upstream_y
   REAL, DIMENSION(:,:), ALLOCATABLE :: u_init_upstream_x, u_init_upstream_y
   REAL, DIMENSION(:,:), ALLOCATABLE :: v_init_upstream_x, v_init_upstream_y
   REAL, DIMENSION(:),   ALLOCATABLE :: t2_init_upstream_x,t2_init_upstream_y, &
                                      q2_init_upstream_x,q2_init_upstream_y,   &
                                      u10_init_upstream_x, u10_init_upstream_y,&
                                      v10_init_upstream_x, v10_init_upstream_y
   REAL, DIMENSION(:),   ALLOCATABLE :: tau_init_u10, tau_init_v10
 
! results of temporal and spatial interpolations
   REAL, ALLOCATABLE, DIMENSION(:)     :: splinetimes, splinetimes_flux,&
        splinetimes_smos, splinetimes_sfc, splinetimes_advection, &
        splinetimes_soil
   REAL, ALLOCATABLE, DIMENSION(:,:)   :: u_g_f,v_g_f,&
                                          t_f_uadv,t_f_vadv,t_f_wadv,&
                                          q_f_uadv,q_f_vadv,q_f_wadv,&
                                          u_f_uadv,u_f_vadv,u_f_wadv,&
                                          v_f_uadv,v_f_vadv,v_f_wadv
   REAL, DIMENSION(:,:), ALLOCATABLE :: tau_u_f, tau_v_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: th_upstream_x_f, th_upstream_y_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: qv_upstream_x_f, qv_upstream_y_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: qc_upstream_x_f, qc_upstream_y_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: qr_upstream_x_f, qr_upstream_y_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: qi_upstream_x_f, qi_upstream_y_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: qg_upstream_x_f, qg_upstream_y_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: u_upstream_x_f, u_upstream_y_f
   REAL, DIMENSION(:,:), ALLOCATABLE :: v_upstream_x_f, v_upstream_y_f
    
   REAL :: h_gamma
   REAL, ALLOCATABLE, DIMENSION(:)     :: glw_f,gsw_f,precip_f
   REAL, ALLOCATABLE, DIMENSION(:,:)   :: smois_f,tslb_f

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

   INTEGER :: ra_lw_physics_i=1,ra_sw_physics_i=1,sst_update=0

! DO eofs related stuff *** Maybe this is not needed at all
  INTEGER  :: ncid_eofs_init, nz_eo, ns_eo, nt_eo

! DO I add the following to pass information from init to uvg
  INTEGER                         :: rnd_sign
  INTEGER, DIMENSION(2)           :: control_index
  REAL, DIMENSION(2)              :: control_w
  REAL, DIMENSION(:), ALLOCATABLE :: gasdo

! Map projection information
  type(proj_info)                :: my_projection
  INTEGER                        :: projcode
  REAL                           :: cent_lat,cent_lon,knowni,knownj, &
                                    truelat1,truelat2,stdlon,        &
                                    sw_corner_lon,sw_corner_lat,     &
                                    lat,lon



! Begin - Additional declarations for interactive radiation ...RT

  LOGICAL :: allowed_to_read

  LOGICAL :: cu_rad_feedback

  LOGICAL :: first_time

  INTEGER :: NPHS, ra_call_offset, julyr

! dummy values
  INTEGER, PARAMETER                  :: n_ozmixm=1,     &
                                         levsiz=1,       &
                                         n_aerosolc=1,   &
                                         paerlev=1,      & 
                                         cam_abs_dim1=1, &
                                         cam_abs_dim2=1

  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: RTHRATENLW, &
                                         RTHRATENSW, &
                                         RTHRATEN, &
                                         CLDFRA

  REAL,  ALLOCATABLE, DIMENSION(:)    :: z_var, &
                                         var

  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: TAUCLDI, &
                                         TAUCLDC, &
                                         PB

  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: F_ICE_PHY, &
                                         F_RAIN_PHY, &
                                         F_RIMEF_PHY
! Optional (only used by CAM lw scheme)
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: abstot, absnxt

  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: emstot
                                         
  REAL, DIMENSION( ims:ime, jms:jme ) :: SWDOWN, SWDOWNC

  REAL, DIMENSION( ims:ime, jms:jme ) :: HTOP, &
                                         HBOT, &
                                         HTOPR, &
                                         HBOTR, &
                                         CUPPT
  REAL, DIMENSION( ims:ime, jms:jme ) ::             &
                                         RLWTOA,     & 
                                         RSWTOA,     & 
                                         ACFRST,     & 
                                         ACFRCV        
  REAL, DIMENSION( ims:ime, jms:jme ) ::             &
                                         CFRACH,     & 
                                         CFRACL,     & 
                                         CFRACM,     & 
                                         CZMEAN        

  INTEGER,DIMENSION( ims:ime, jms:jme ) ::  &
                                    NCFRST, &  
                                    NCFRCV     
! stuff
  REAL, DIMENSION( ims:ime, jms:jme ) ::                         &
                      ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC,          &
                      ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC,          &
                      ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC,          &
                      ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC
  REAL, DIMENSION( ims:ime, jms:jme ) ::                         &
                        SWUPT,  SWUPTC,  SWDNT,  SWDNTC,          &
                        SWUPB,  SWUPBC,  SWDNB,  SWDNBC,          &
                        LWUPT,  LWUPTC,  LWDNT,  LWDNTC,          &
                        LWUPB,  LWUPBC,  LWDNB,  LWDNBC

  REAL, DIMENSION( ims:ime, jms:jme ) ::                          &
                                          SWCF,                   &
                                          LWCF,                   &
                                          OLR

  REAL,  DIMENSION( ims:ime, levsiz, jms:jme, n_ozmixm ) ::       &
                                          OZMIXM

  REAL,  DIMENSION(levsiz) :: PIN

  REAL,  DIMENSION(ims:ime,jms:jme) :: m_ps_1,m_ps_2

  REAL,  DIMENSION( ims:ime, paerlev, jms:jme, n_aerosolc ) ::    &
      aerosolc_1, aerosolc_2

  REAL,  DIMENSION(paerlev) :: m_hybi0

  REAL :: JULIAN,SOLCON, GMT, XTIME, swrad_scat
  REAL :: w1, w2

  INTEGER :: STEPRA, ICLOUD
  INTEGER :: itf, jtf, ktf, m

  REAL, DIMENSION(13) :: solattjdayref
  REAL, DIMENSION(32) :: solattzref
  REAL :: solatt(1:32,1:14)
  REAL :: t, u, solatttop

  REAL ::   cam_abs_freq_s

! End - Additional declarations for interactive radiation ...RT

! Begin declarations for stochastic clouds ...RT
  INTEGER, ALLOCATABLE, DIMENSION(:) :: stoch_run_in, stoch_yn_in
  REAL, ALLOCATABLE, DIMENSION(:) :: stoch_time_in, stoch_cbh_in, stoch_lwp_in

  INTEGER, ALLOCATABLE, DIMENSION(:) :: stoch_yn
  REAL, ALLOCATABLE, DIMENSION(:) :: stoch_time, stoch_cbh, stoch_lwp


  INTEGER :: nbstoch, indexrealization, indexrun, data_dt_sec, nb_data_intervals
  REAL :: randraw

! END declarations for stochastic clouds ...RT


contains

SUBROUTINE STATIC_INIT_WRF(allocate_wrf)

  LOGICAL,INTENT(INOUT)           :: allocate_wrf

  kde = nz
  kme = nz
  kte = nz - 1

  IF (.NOT.init_f) allocate_wrf=.TRUE.

  IF ( allocate_wrf ) THEN
    ALLOCATE(moist( ims:ime, kms:kme, jms:jme, n_moist))
    ALLOCATE(scalar( ims:ime, kms:kme, jms:jme, 1)) ! don't need qni and qt

    ALLOCATE(p_phy(ims:ime, kms:kme, jms:jme),            &
             pi_phy(ims:ime, kms:kme, jms:jme), &
             rho(ims:ime, kms:kme, jms:jme), &
             p8w(ims:ime, kms:kme, jms:jme), & 
             t_phy(ims:ime, kms:kme, jms:jme), &
             u_phy(ims:ime, kms:kme, jms:jme), &
             v_phy(ims:ime, kms:kme, jms:jme), &
             dz8w(ims:ime, kms:kme, jms:jme), &
             z(ims:ime, kms:kme, jms:jme), &
             th_phy(ims:ime, kms:kme, jms:jme) )
    ALLOCATE(z8w(ims:ime, kms:kme+1, jms:jme) )
    ALLOCATE(RUBLTEN(ims:ime, kms:kme, jms:jme), &
             RVBLTEN(ims:ime, kms:kme, jms:jme), &
             RTHBLTEN(ims:ime, kms:kme, jms:jme), &
             RQVBLTEN(ims:ime, kms:kme, jms:jme), &
             RQCBLTEN(ims:ime, kms:kme, jms:jme), &
             RQIBLTEN(ims:ime, kms:kme, jms:jme), &
             TKE_MYJ(ims:ime, kms:kme, jms:jme), &
             &EXCH_H(ims:ime, kms:kme, jms:jme), &
             &EL_MYJ(ims:ime, kms:kme, jms:jme) )

! begin - radiation ...RT
    ALLOCATE(RTHRATENLW(ims:ime, kms:kme, jms:jme),&
             RTHRATENSW(ims:ime, kms:kme, jms:jme),&
             RTHRATEN(ims:ime, kms:kme, jms:jme), & 
             t8w(ims:ime, kms:kme, jms:jme), & 
             CLDFRA(ims:ime, kms:kme, jms:jme) )
    ALLOCATE(TAUCLDC(ims:ime, kms:kme, jms:jme),&
             TAUCLDI(ims:ime, kms:kme, jms:jme),&
             PB(ims:ime, kms:kme, jms:jme) )
    ALLOCATE(F_ICE_PHY(ims:ime, kms:kme, jms:jme), &
             F_RAIN_PHY(ims:ime, kms:kme, jms:jme) )
    ALLOCATE(z_var(1:nz+1), var(1:nz+1))

! Optional (only used by CAM lw scheme)
    ALLOCATE(abstot(ims:ime, kms:kme, cam_abs_dim2, jms:jme), &
             absnxt(ims:ime, kms:kme, cam_abs_dim1, jms:jme), &
             emstot(ims:ime, kms:kme, jms:jme))

! end - radiation ...RT

    ALLOCATE(z_init(1:nz), u_init(1:nz), v_init(1:nz), t_init(1:nz), &
             th_init(1:nz), exn_init(1:nz), q_init(1:nz), p_init(1:nz), &
             rho_init(1:nz))
    ALLOCATE(z8w_init(1:nz+1), p8w_init(1:nz+1))
    ALLOCATE(z_grid_stag(1:nz+1), z_grid(1:nz))
    ALLOCATE(uflux(1:nz), vflux(1:nz), hflux(1:nz), &
             qflux(1:nz), k_t(1:nz),k_m(1:nz))
    ALLOCATE(u_g(1:nz), v_g(1:nz))
    
    ALLOCATE(gasdo(n_eo))

!   go ahead and read the grid here - redundant
    OPEN(55,file=gridfile) ! grid now specified in namelist ...RT

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
    t8w = 0.0 ! added for radiation ...RT
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
     ! JPH kludge soil depth for DART
     zs = (/.05,.25,.70,1.50/)
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
!  ZS = 0.0

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
  nsplinetimes_soil=INT((REAL(forecast_length)/REAL(splineinterval_soil))+&
       &.999)+1

! this is needed only for the DART interface
  n1dsplines = 2*nsplinetimes_flux + nsplinetimes_smos

   IF (init_f) THEN
      ntime=NINT(REAL(forecast_length)/dt)+1
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
         CALL map_init(my_projection)
         CALL wrf_f_dims(ncid_f, nz_f,ns_f, nt_f,            &
                        dx,cent_lat,cent_lon,stdlon,         &
                        truelat1,truelat2,mminlu,julday,     &
                        sw_corner_lon,sw_corner_lat,projcode,&
                        lat, lon, cor)
         nt_f_flux = nt_f
         nt_f_soil = nt_f
         nt_f_smos = nt_f
         nt_uvg    = nt_f
         nz_uvg    = nz_f
          
         ! need map info for rotation
         knowni = -9999.0
         knownj = -9999.0  ! don't actually need these
         CALL map_set(projcode, sw_corner_lat, sw_corner_lon,  &
                      knowni, knownj,                          &
                      dx, stdlon, truelat1, truelat2, my_projection)

      CASE('OBS' , 'SFC')
         print*,'Map information set to trivial for init type ',init_f_type
         my_projection%code = PROJ_LATLON

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
           splinetimes_advection(nsplinetimes_advection),&
           splinetimes_flux(nsplinetimes_flux),&
           splinetimes_smos(nsplinetimes_smos),&
           splinetimes_soil(nsplinetimes_soil),&
           u_g_f(nz,nsplinetimes),v_g_f(nz,nsplinetimes),&
           glw_f(nsplinetimes_flux),gsw_f(nsplinetimes_flux),&
           precip_f(nsplinetimes_smos)) 
        IF ( force_soil ) THEN
           ALLOCATE (smois_f(num_soil_layers,nsplinetimes_soil), &
           tslb_f(num_soil_layers,nsplinetimes_soil))
           tslb_f = 0.
           smois_f = 0.
        ENDIF
        ALLOCATE(u_f(nz,nsplinetimes),    &
                  v_f(nz,nsplinetimes),   &
                  t_f(nz,nsplinetimes),   &
                  p_f(nz,nsplinetimes),   &
                  p8w_f(nz,nsplinetimes), &
                  q_f(nz,nsplinetimes))
        ALLOCATE(u_fi(nz),  &
                  v_fi(nz), &
                  t_fi(nz), &
                  p_fi(nz), &
                  q_fi(nz))
        IF ( trim(nudge_f_type) == 'WRF' ) ALLOCATE(nudge_func(nz))
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
              exn_init_f(nz_f,nt_f),qv_init_f(nz_f,nt_f), &
              qc_init_f(nz_f,nt_f),qr_init_f(nz_f,nt_f), &
              qi_init_f(nz_f,nt_f),qg_init_f(nz_f,nt_f), &
              p_init_f(nz_f,nt_f), z_f(nz_f,nt_f)) 
        ALLOCATE(t_init_upstream_x(nz_f,nt_f), &
                     t_init_upstream_y(nz_f,nt_f), &
                     qv_init_upstream_x(nz_f,nt_f), &
                     qv_init_upstream_y(nz_f,nt_f), &
                     qc_init_upstream_x(nz_f,nt_f), &
                     qc_init_upstream_y(nz_f,nt_f), &
                     qr_init_upstream_x(nz_f,nt_f), &
                     qr_init_upstream_y(nz_f,nt_f), &
                     qi_init_upstream_x(nz_f,nt_f), &
                     qi_init_upstream_y(nz_f,nt_f), &
                     qg_init_upstream_x(nz_f,nt_f), &
                     qg_init_upstream_y(nz_f,nt_f), &
                     u_init_upstream_x(nz_f,nt_f), &
                     u_init_upstream_y(nz_f,nt_f), &
                     v_init_upstream_x(nz_f,nt_f), &
                     v_init_upstream_y(nz_f,nt_f), &
                     t2_init_upstream_x(nt_f), &
                     t2_init_upstream_y(nt_f), &
                     q2_init_upstream_x(nt_f), &
                     q2_init_upstream_y(nt_f), &
                     u10_init_upstream_x(nt_f), &
                     u10_init_upstream_y(nt_f), &
                     v10_init_upstream_x(nt_f), &
                     v10_init_upstream_y(nt_f), &
                     tau_init_u(nz_f, nt_f), &
                     tau_init_v(nz_f, nt_f), &
                     tau_init_u10(nt_f), &
                     tau_init_v10(nt_f), &
                     th_upstream_x_f(nz,nsplinetimes_advection), &
                     th_upstream_y_f(nz,nsplinetimes_advection), &
                     qv_upstream_x_f(nz,nsplinetimes_advection), &
                     qv_upstream_y_f(nz,nsplinetimes_advection), &
                     u_upstream_x_f(nz,nsplinetimes_advection), &
                     u_upstream_y_f(nz,nsplinetimes_advection), &
                     v_upstream_x_f(nz,nsplinetimes_advection), &
                     v_upstream_y_f(nz,nsplinetimes_advection), &
                     tau_u_f(nz, nsplinetimes_advection), &
                     tau_v_f(nz, nsplinetimes_advection)) 
        ALLOCATE(tau_u(1:nz), tau_v(1:nz), &
                     th_upstream_x(1:nz),th_upstream_y(1:nz), &
                     qv_upstream_x(1:nz),qv_upstream_y(1:nz), &
                     u_upstream_x(1:nz),u_upstream_y(1:nz),   &
                     v_upstream_x(1:nz),v_upstream_y(1:nz))
        IF ( P_QC > 1 ) THEN
           ALLOCATE(qc_init_upstream_x(nz_f,nt_f), &
                     qc_init_upstream_y(nz_f,nt_f), &
                     qc_upstream_x_f(nz,nsplinetimes_advection), &
                     qc_upstream_y_f(nz,nsplinetimes_advection), &
                     qc_upstream_x(1:nz),qc_upstream_y(1:nz))
           qc_init_upstream_x = 0.0
           qc_init_upstream_y = 0.0
           qc_upstream_x_f = 0.0
           qc_upstream_y_f = 0.0
           qc_upstream_x = 0.0
           qc_upstream_y = 0.0
        ENDIF
        IF ( P_QR > 1 ) THEN
           ALLOCATE(qr_init_upstream_x(nz_f,nt_f), &
                    qr_init_upstream_y(nz_f,nt_f), &
                    qr_upstream_x_f(nz,nsplinetimes_advection), &
                    qr_upstream_y_f(nz,nsplinetimes_advection), &
                    qr_upstream_x(1:nz),qr_upstream_y(1:nz))
           qr_init_upstream_x = 0.0
           qr_init_upstream_y = 0.0
           qr_upstream_x_f = 0.0
           qr_upstream_y_f = 0.0
           qr_upstream_x = 0.0
           qr_upstream_y = 0.0
        ENDIF
        IF ( P_QI > 1 ) THEN
           ALLOCATE(qi_init_upstream_x(nz_f,nt_f), &
                     qi_init_upstream_y(nz_f,nt_f), &
                     qi_upstream_x_f(nz,nsplinetimes_advection), &
                     qi_upstream_y_f(nz,nsplinetimes_advection), &
                     qi_upstream_x(1:nz),qi_upstream_y(1:nz))
           qi_init_upstream_x = 0.0
           qi_init_upstream_y = 0.0
           qi_upstream_x_f = 0.0
           qi_upstream_y_f = 0.0
           qi_upstream_x = 0.0
           qi_upstream_y = 0.0
        ENDIF
        IF ( P_QG > 1 ) THEN
           ALLOCATE(qg_init_upstream_x(nz_f,nt_f), &
                     qg_init_upstream_y(nz_f,nt_f), &
                     qg_upstream_x_f(nz,nsplinetimes_advection), &
                     qg_upstream_y_f(nz,nsplinetimes_advection), &
                     qg_upstream_x(1:nz),qg_upstream_y(1:nz))
           qg_init_upstream_x = 0.0
           qg_init_upstream_y = 0.0
           qg_upstream_x_f = 0.0
           qg_upstream_y_f = 0.0
           qg_upstream_x = 0.0
           qg_upstream_y = 0.0
        ENDIF
        tau_init_u = 0.0
        tau_init_v = 0.0
        tau_init_u10 = 0.0
        tau_init_v10 = 0.0
        tau_u_f = 0.0
        tau_v_f = 0.0
        tau_u = 0.0
        tau_v = 0.0
        t_init_upstream_x = 0.0
        t_init_upstream_y = 0.0
        qv_init_upstream_x = 0.0
        qv_init_upstream_y = 0.0
        u_init_upstream_x = 0.0
        u_init_upstream_y = 0.0
        v_init_upstream_x = 0.0
        v_init_upstream_y = 0.0
        t2_init_upstream_x = 0.0
        t2_init_upstream_y = 0.0
        q2_init_upstream_x = 0.0
        q2_init_upstream_y = 0.0
        u10_init_upstream_x = 0.0
        u10_init_upstream_y = 0.0
        v10_init_upstream_x = 0.0
        v10_init_upstream_y = 0.0
        th_upstream_x_f = 0.0
        th_upstream_y_f = 0.0
        th_upstream_x = 0.0
        th_upstream_y = 0.0
        qv_upstream_x_f = 0.0
        qv_upstream_y_f = 0.0
        qv_upstream_x = 0.0
        qv_upstream_y = 0.0
        u_upstream_x_f = 0.0
        u_upstream_y_f = 0.0
        v_upstream_x_f = 0.0
        v_upstream_y_f = 0.0
        u_upstream_x = 0.0
        u_upstream_y = 0.0
        v_upstream_x = 0.0
        v_upstream_y = 0.0
        u_f = 0.0
        v_f = 0.0
        p_f = 0.0
        p8w_f = 0.0
        t_f = 0.0
        q_f = 0.0
        IF ( trim(nudge_f_type) == 'WRF' ) nudge_func = 0.0

        ALLOCATE(z_f_stag(nz_f+1,nt_f))
        ALLOCATE(zs_f(ns_f), dzs_f(ns_f))
        ALLOCATE(tslb_init_f(ns_f,nt_f_soil), &
              smois_init_f(ns_f,nt_f_soil))
         
        IF (force_flux) THEN
            ALLOCATE(&
                 ts_init_f(nsplinetimes_sfc),&
                 qvs_init_f(nsplinetimes_sfc),&
                 ustar_init_f(nsplinetimes_sfc),&
                 hflux_init_f(nsplinetimes_sfc),&
                 qvflux_init_f(nsplinetimes_sfc))
            
           ALLOCATE(times_f_sfc(nt_f_sfc),&
                 ts_f(nt_f_sfc),qvs_f(nt_f_sfc),&
                 ustar_f(nt_f_sfc),hflux_f(nt_f_sfc),qvflux_f(nt_f_sfc))

           ALLOCATE(splinetimes_sfc(nsplinetimes_sfc))

        ELSE

! otherwise would need optional parameters in surface_driver

           ALLOCATE(&
                 ts_init_f(1),&
                 qvs_init_f(1),&
                 ustar_init_f(1),&
                 hflux_init_f(1),&
                 qvflux_init_f(1))
            
           ALLOCATE(times_f_sfc(1),&
                 ts_f(1),qvs_f(1),&
                 ustar_f(1),hflux_f(1),qvflux_f(1))

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
! --------------------------
! for stochastic cloud ...RT
! --------------------------
      IF (forc_stochastic_cloud) then
! JPH key length of cloud records off namelist
         nb_data_intervals = INT(cloud_realization_length/dt)
         print*, 'nbrealizations= ', nbrealizations
         nbstoch = (nb_data_intervals)*nbrealizations
         print*, 'nbstoch= ', nbstoch

         print*, 'Allocating arrays for stochastic cloud with length: ', nbstoch
         ALLOCATE( stoch_yn_in(nbstoch),  &
                   stoch_cbh_in(nbstoch), &
                   stoch_lwp_in(nbstoch), &
                   stoch_time_in(nbstoch), &
                   stoch_run_in(nbstoch) )

         ALLOCATE( stoch_yn(ntime),  &
                   stoch_cbh(ntime), &
                   stoch_lwp(ntime), &
                   stoch_time(ntime) )

! read lookup table...for coeffs in mean/ARMA model...(not yet)
! or
! read data file (obs or simul)
         print*, 'Reading file of stochastic cloud data... '
         OPEN(stochaunit,file=cld_file)
         DO k=1, nbstoch
            read(stochaunit,*) stoch_run_in(k),stoch_time_in(k), &
             stoch_yn_in(k),stoch_cbh_in(k),stoch_lwp_in(k)
         ENDDO
         CLOSE(stochaunit)         

         indexrealization = 0

! checking data time intervals and model time step

         data_dt_sec = NINT((stoch_time_in(2) - stoch_time_in(1))*3600.0) ! in sec
         IF (data_dt_sec.NE.INT(dt)) then
            print*, 'Cloud data time interval do not match model time step.'
            print*, 'Stopping model'
            stop
         endif

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

            CALL wrf_init_and_bc(my_projection,lon,ncid_f,nz_f, ns_f,     &
                    nt_f,n_moist,                                         &
                    z_f,z_f_stag,t_init_f,u_init_f,v_init_f,              &
                    qv_init_f,qc_init_f,qr_init_f,qi_init_f,qg_init_f,    &
                    p_init_f,                                    &
                    t_init_upstream_x, t_init_upstream_y,                 &
                    qv_init_upstream_x, qv_init_upstream_y,               &
                    qc_init_upstream_x, qc_init_upstream_y,               &
                    qr_init_upstream_x, qr_init_upstream_y,               &
                    qi_init_upstream_x, qi_init_upstream_y,               &
                    qg_init_upstream_x, qg_init_upstream_y,               &
                    u_init_upstream_x, u_init_upstream_y,                 &
                    v_init_upstream_x, v_init_upstream_y,                 &
                    t2_init_upstream_x, t2_init_upstream_y,               &
                    q2_init_upstream_x, q2_init_upstream_y,               &
                    u10_init_upstream_x, u10_init_upstream_y,             &
                    v10_init_upstream_x, v10_init_upstream_y,             &
                    tau_init_u, tau_init_v, tau_init_u10, tau_init_v10,   &
                    th2_init_f,t2_init_f,tsk_init_f,                      &
                    u10_init_f,v10_init_f,                                &
                    q2_init_f,glw_init_f,gsw_init_f,qsfc_init_f,          &
                    tslb_init_f,smois_init_f,tmn_init_f,                  &
                    precip_init_f,                                        &
                    vegfra_f,isltyp_f,lu_index_f,ivgtyp_f, terrain_f, dx, &
                    times_f,times_f_flux,times_f_soil,times_f_smos,idum,  &
                    control_index,rnd_sign,control_w,gasdo)


         CASE('OBS' , 'SFC')

            CALL snd_init_and_bc(ncid_f, ncid_flux, ncid_soil, ncid_smos, &
                    nz_f, ns_f, nt_f, nt_f_flux, nt_f_soil, nt_f_smos,&
                    z_f,t_init_f,u_init_f,v_init_f,&
                    qv_init_f,p_init_f,&
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

         CASE('OBS' , 'SFC')

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

     IF (force_soil) THEN
       DO i=1,nsplinetimes_soil
           splinetimes_soil(i)=times_f(1)+(i-1)*splineinterval_soil
       ENDDO
     ENDIF

     IF (force_flux) THEN
        DO i=1,nsplinetimes_smos
           splinetimes_sfc(i)=times_f(1)+(i-1)*splineinterval_sfc
        ENDDO
     ENDIF

     CALL initf(itime_f,u_init_f(:,1:nt_f),v_init_f(:,1:nt_f),&
           t_init_f(:,1:nt_f),qv_init_f(:,1:nt_f),qc_init_f(:,1:nt_f),&
           qr_init_f(:,1:nt_f),qi_init_f(:,1:nt_f),qg_init_f(:,1:nt_f), &
           p_init_f(:,1:nt_f), &
           t_init_upstream_x(:,1:nt_f), t_init_upstream_y(:,1:nt_f), &
           qv_init_upstream_x(:,1:nt_f), qv_init_upstream_y(:,1:nt_f), &
           qc_init_upstream_x(:,1:nt_f), qc_init_upstream_y(:,1:nt_f), &
           qr_init_upstream_x(:,1:nt_f), qr_init_upstream_y(:,1:nt_f), &
           qi_init_upstream_x(:,1:nt_f), qi_init_upstream_y(:,1:nt_f), &
           qg_init_upstream_x(:,1:nt_f), qg_init_upstream_y(:,1:nt_f), &
           u_init_upstream_x(:,1:nt_f), u_init_upstream_y(:,1:nt_f), &
           v_init_upstream_x(:,1:nt_f), v_init_upstream_y(:,1:nt_f), &
           t2_init_upstream_x(1:nt_f), t2_init_upstream_y(1:nt_f), &
           q2_init_upstream_x(1:nt_f), q2_init_upstream_y(1:nt_f), &
           u10_init_upstream_x(1:nt_f), u10_init_upstream_y(1:nt_f), &
           v10_init_upstream_x(1:nt_f), v10_init_upstream_y(1:nt_f), &
           tau_init_u(:,1:nt_f), tau_init_v(:,1:nt_f), &
           tau_init_u10(1:nt_f), tau_init_v10(1:nt_f), &
           glw_init_f(1:nt_f_flux),gsw_init_f(1:nt_f_flux),&
           precip_init_f(1:nt_f_smos),&
           tslb_init_f(:,1:nt_f_soil), smois_init_f(:,1:nt_f_soil), &
           ts_init_f(1:nt_f_sfc),&
           qvs_init_f(1:nt_f_sfc),&
           ustar_init_f(1:nt_f_sfc),&
           hflux_init_f(1:nt_f_sfc),&
           qvflux_init_f(1:nt_f_sfc),&
           th2_init_f,q2_init_f,&
           u10_init_f,v10_init_f,&
           tsk_init_f,qsfc_init_f,&
           z_f(:,1:nt_f),nz_f,z_g,nt_f,nt_f_flux,nt_f_smos,nt_f_sfc,&
           zs_f, ns_f, nt_f_soil, &
           times_f(1:nt_f),times_f_flux(1:nt_f_flux),&
           times_f_smos(1:nt_f_smos),&
           times_f_sfc(1:nt_f_sfc), times_f_soil(1:nt_f_soil),&
           nsplinetimes,splinetimes,&
           nsplinetimes_advection,splinetimes_advection,&
           nsplinetimes_flux,splinetimes_flux,&
           nsplinetimes_smos,splinetimes_smos,&
           nsplinetimes_sfc,splinetimes_sfc,&
           nsplinetimes_soil,splinetimes_soil,&
           pblh_ref,stepbl,lowlyr,ht,pblh,&
           th_phy,t_phy,moist,u_phy,v_phy,p_phy, p8w, tke_myj,&
           u_g_f,v_g_f,glw_f,gsw_f,precip_f, &
           u_f, v_f, t_f, p_f,p8w_f,q_f,&
           tslb_f, smois_f, &
           th_upstream_x_f, th_upstream_y_f, &
           qv_upstream_x_f, qv_upstream_y_f, &
           qc_upstream_x_f, qc_upstream_y_f, &
           qr_upstream_x_f, qr_upstream_y_f, &
           qi_upstream_x_f, qi_upstream_y_f, &
           qg_upstream_x_f, qg_upstream_y_f, &
           u_upstream_x_f, u_upstream_y_f, &
           v_upstream_x_f, v_upstream_y_f, &
           tau_u_f, tau_v_f, &
           ts_f,qvs_f,ustar_f,hflux_f,qvflux_f,&
           z,z8w,dz8w,zs,num_soil_layers, &
           ims,ime,jms,jme,kms,kme)

     zl2=z(:,2,:)

! initialize the diagnostic variables in case you want a prior
     t2 = t2_init_f(1)
     q2 = q2_init_f(1)
     u10 = u10_init_f(1)
     v10 = v10_init_f(1)

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
                      control_index,rnd_sign,control_w,gasdo,nz_stag_uvg, &
                      nrecords_uvg)

       CASE DEFAULT
        print*, 'Do not know how to force from type ',force_f_type
        stop 'module_wrf'
       
       END SELECT

        wrf_rnd_seed = idum

! initialize wind profile to geostrophic ...RT
!            DO i=ims,ime
!               DO j=jms,jme
!                  DO k=kme, kms, -1
!                     u_phy(i,k,j) = u_g_f(k,1)
!                     v_phy(i,k,j) = v_g_f(k,1)
!                  ENDDO
!               ENDDO
!            ENDDO

     ENDIF


  ENDIF ! end if init ideal or real

   k_t=0.
   k_m=0.
   uflux=0.
   vflux=0.
   hflux=0.
   qflux=0.

!! -----------------------------------------
!! initialize stochastic cloud process ...RT
!! -----------------------------------------
      IF (forc_stochastic_cloud) then

         print*,'cld_sequence= ', cld_sequence 
         IF (cld_sequence.EQ.'rnd') then ! random selection
            randraw = ran1(idum)*(nbrealizations-1) + 1
            indexrealization = AINT(randraw)                  
!            IF (indexrealization.eq.0) indexrealization = 1
         ELSE 
            IF (cld_sequence.EQ.'seq') then ! sequential selection
            indexrealization = indexrealization + 1
            ELSE !unknown option...stop
               print*, 'Option not recognized for cld_sequence ',cld_sequence
               stop 'module_wrf'
            ENDIF

         ENDIF

         print*, 'nb_data_intervals= ',nb_data_intervals
         print*, 'indexrealization= ', indexrealization
         indexrun = (indexrealization*nb_data_intervals)-nb_data_intervals+1

! JPH gets incremented further to start of simulation time
         IF ( cloud_seq_start > start_forecast ) THEN
           print*,'You are trying to start the model before cloud sequences are available, make start_forecast >= cloud_seq_start'
           stop 'module_wrf'
         ENDIF
         indexrun = indexrun + AINT((start_forecast-cloud_seq_start)/dt)
         print*,'indexrun = ',indexrun

! assign cloud parameters to vector used in the simulation
         DO k=1, ntime
            stoch_time(k) = stoch_time_in(indexrun+k-1)
            stoch_yn(k) = stoch_yn_in(indexrun+k-1)
            stoch_cbh(k) = stoch_cbh_in(indexrun+k-1)
            stoch_lwp(k) = stoch_lwp_in(indexrun+k-1)
         ENDDO

         wrf_rnd_seed = idum

      ENDIF

! Advection scaling is one at start if nothing else done in DART
   tadvection_scale = 1.0
   qadvection_scale = 1.0
   uadvection_scale = 1.0
   vadvection_scale = 1.0

!! -------------------------------
!! initialize radiation begin...RT
!! -------------------------------
 IF (ra_type == 'INT') THEN

   allowed_to_read = .TRUE.
   cu_rad_feedback = .FALSE.
   first_time = .TRUE.
   swrad_scat = 1. ! scattering tuning parameter...from WRF namelist

   GMT = start_hour_f + start_minute_f/60.0 + start_forecast / 3600.0 

!---------------------
   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)
!---------------------

!-- calculate radiation time step

   STEPRA = nint(RADT*60./DT)
   STEPRA = max(STEPRA,1)

!-- initialization

   DO j=jts,jtf
   DO k=kts,ktf
   DO i=its,itf
      RTHRATEN(i,k,j)=0.
      RTHRATENLW(i,k,j)=0.
      RTHRATENSW(i,k,j)=0.
      CLDFRA(i,k,j)=0.
   ENDDO
   ENDDO
   ENDDO

!-- chose long wave radiation scheme
 
   lwrad_select: SELECT CASE(ra_lw_physics)

        CASE ('RRTMSCHEME')
             CALL rrtminit(                                 &
                           allowed_to_read ,                &
                           ids, ide, jds, jde, kds, kde,    &
                           ims, ime, jms, jme, kms, kme,    &
                           its, ite, jts, jte, kts, kte     )

        CASE ('CAMLWSCHEME')
           WRITE( wrf_err_message , * )'This LW option is not supported' !RT
           PRINT*,' This LW option is not supported'
           PRINT*,'Stopping'
           STOP
!RT#ifdef MAC_KLUDGE
!RT             CALL wrf_error_fatal ( 'CAM radiation scheme not supported under the chosen build configuration' )
!RT#endif
!RT             IF ( PRESENT( OZMIXM ) .AND. PRESENT( PIN ) .AND. &
!RT                  PRESENT(M_PS_1) .AND. PRESENT(M_PS_2) .AND.  &
!RT                  PRESENT(M_HYBI) .AND. PRESENT(AEROSOLC_1)    &
!RT                  .AND. PRESENT(AEROSOLC_2)) THEN
!RT             CALL camradinit(                                  &
!RT                         R_D,R_V,CP,G,STBOLT,EP_2,shalf,pptop, &
!RT                         ozmixm,pin,levsiz,XLAT,n_ozmixm,      &
!RT                         m_ps_1,m_ps_2,m_hybi,aerosolc_1,aerosolc_2,&
!RT                         paerlev, n_aerosolc,              &
!RT                         ids, ide, jds, jde, kds, kde,     &
!RT                         ims, ime, jms, jme, kms, kme,     &
!RT                         its, ite, jts, jte, kts, kte      )
!RT             ELSE
!RT                CALL wrf_error_fatal ( 'arguments not present for calling cam radiation' )
!RT             ENDIF
!RT
!RT            camlw = .true.

        CASE ('GFDLLWSCHEME')
           WRITE( wrf_err_message , * )'This LW option is not supported' !RT
           PRINT*,' This LW option is not supported'
           PRINT*,'Stopping'
           STOP
!RT             CALL nl_get_start_month(id,month)
!RT             CALL nl_get_start_day(id,iday)
!RT             CALL gfdletainit(emiss,sfull,shalf,pptop,      &
!RT                              julyr,month,iday,gmt,         &
!RT                              config_flags,allowed_to_read, &
!RT                              ids, ide, jds, jde, kds, kde, &
!RT                              ims, ime, jms, jme, kms, kme, &
!RT                              its, ite, jts, jte, kts, kte  )
!RT             etalw = .true.
        CASE DEFAULT

   END SELECT lwrad_select
!-- initialize short wave radiation scheme
 
   swrad_select: SELECT CASE(ra_sw_physics)

        CASE ('SWRADSCHEME')
!RT             CALL swinit(                                  &
!RT                         swrad_scat,                       &
!RT                         allowed_to_read ,                 &
!RT                         ids, ide, jds, jde, kds, kde,     &
!RT                         ims, ime, jms, jme, kms, kme,     &
!RT                         its, ite, jts, jte, kts, kte      )

           cssca = swrad_scat * 1.e-5

        CASE ('CAMSWSCHEME')
           WRITE( wrf_err_message , * )'This SW option is not supported' !RT
           PRINT*,' This SW option is not supported'
           PRINT*,'Stopping'
           STOP
!RT#ifdef MAC_KLUDGE
!RT             CALL wrf_error_fatal ( 'CAM radiation scheme not supported under the chosen build configuration' )
!RT#endif
!RT             IF(.not.camlw)THEN
!RT             CALL camradinit(                              &
!RT                         R_D,R_V,CP,G,STBOLT,EP_2,shalf,pptop,               &
!RT                         ozmixm,pin,levsiz,XLAT,n_ozmixm,     &
!RT                         m_ps_1,m_ps_2,m_hybi,aerosolc_1,aerosolc_2,&
!RT                         paerlev, n_aerosolc,              &
!RT                         ids, ide, jds, jde, kds, kde,     &
!RT                         ims, ime, jms, jme, kms, kme,     &
!RT                         its, ite, jts, jte, kts, kte      )
!RT             ENDIF

        CASE ('GSFCSWSCHEME')
           WRITE( wrf_err_message , * )'This SW option is not supported' !RT
           PRINT*,' This SW option is not supported'
           PRINT*,'Stopping'
           STOP
!RT             CALL gsfc_swinit(cen_lat, allowed_to_read )

        CASE ('GFDLSWSCHEME')
           WRITE( wrf_err_message , * )'This SW option is not supported' !RT
           PRINT*,' This SW option is not supported'
           PRINT*,'Stopping'
           STOP
!RT             IF(.not.etalw)THEN
!RT             CALL nl_get_start_month(id,month)
!RT             CALL nl_get_start_day(id,iday)
!RT             CALL gfdletainit(emiss,sfull,shalf,pptop,      &
!RT                              julyr,month,iday,gmt,         &
!RT                              config_flags,allowed_to_read, &
!RT                              ids, ide, jds, jde, kds, kde, &
!RT                              ims, ime, jms, jme, kms, kme, &
!RT                              its, ite, jts, jte, kts, kte  )
!RT             ENDIF

        CASE DEFAULT

   END SELECT swrad_select


 ENDIF
!! -----------------------------
!! initialize radiation end...RT
!! -----------------------------

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
      PRINT *,'Unknown boundary layer option ',blscheme
      PRINT *,'Stopping'
      STOP
   ENDIF
   
   OPEN(isfcunit ,file=TRIM(outdir)//'/'//&
        &TRIM(blscheme)//TRIM(sfscheme)//'_sfc.txt')
   OPEN(iprofunit,file=TRIM(outdir)//'/'//&
        &TRIM(blscheme)//TRIM(sfscheme)//'_prof.txt')
   OPEN(isimil,  file=TRIM(outdir)//'/'//&
        &TRIM(blscheme)//TRIM(sfscheme)//'_simil.txt')

   bl_select: SELECT CASE(trim(bl_pbl_physics))
      CASE('MYJPBLSCHEME') 


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
     
      CASE('YSUPBLSCHEME') 

      CALL sfclayinit(.TRUE.)
      CALL ysuinit(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,    &
           RQCBLTEN,RQIBLTEN,P_QI,               &
           PARAM_FIRST_SCALAR,                   &
           .FALSE.,.TRUE.,                              &
           ids, ide, jds, jde, kds, kde,         &
           ims, ime, jms, jme, kms, kme,         &
           its, ite, jts, jte, kts, kte          )

      CASE('MRFPBLSCHEME') 

      CALL sfclayinit(.TRUE.)
      CALL mrfinit(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,    &
           RQCBLTEN,RQIBLTEN,P_QI,               &
           PARAM_FIRST_SCALAR,                   &
           .FALSE.,.TRUE.,                       &
           ids, ide, jds, jde, kds, kde,         &
           ims, ime, jms, jme, kms, kme,         &
           its, ite, jts, jte, kts, kte          )

      CASE('GFSPBLSCHEME')

      CALL gfsinit(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,    &
           RQCBLTEN,RQIBLTEN,P_QI,               &
           PARAM_FIRST_SCALAR,                   &
           .FALSE.,.TRUE.,                       &
           ids, ide, jds, jde, kds, kde,         &
           ims, ime, jms, jme, kms, kme,         &
           its, ite, jts, jte, kts, kte          )
      
      CASE DEFAULT
    
      print*,bl_pbl_physics,' not a valid choice for PBL scheme, exiting!'
      stop 'module_wrf'

   END SELECT bl_select

 
! initialise microphysics
   IF ( n_moist > 1 ) then
     IF ( P_QV > 1 ) f_qv=.TRUE.
     IF ( P_QC > 1 ) f_qc=.TRUE.
     IF ( P_QR > 1 ) f_qr=.TRUE.
     IF ( P_QI > 1 ) f_qi=.TRUE.
     IF ( P_QS > 1 ) f_qs=.TRUE.
     IF ( P_QG > 1 ) f_qg=.TRUE.
   ENDIF

   micro_select: SELECT CASE (trim(mp_physics))
     CASE ('KESSLERSCHEME') 
       print*,mp_physics, ' not yet a valid choice, exiting'
       stop 'module_wrf'
     CASE ('THOMPSON') 
       print*,mp_physics, ' not yet a valid choice, exiting'
       stop 'module_wrf'
     CASE ('LINSCHEME') 
       ! how nice, no init necessary
     CASE ('WSM3SCHEME') 
       print*,mp_physics, ' not yet a valid choice, exiting'
       stop 'module_wrf'
     CASE ('WSM5SCHEME') 
       print*,mp_physics, ' not yet a valid choice, exiting'
       stop 'module_wrf'
     CASE ('WSM6SCHEME') 
       print*,mp_physics, ' not yet a valid choice, exiting'
       stop 'module_wrf'
     CASE ('ETAMPNEW') 
       print*,mp_physics, ' not yet a valid choice, exiting'
       stop 'module_wrf'
     CASE ('NCEPCLOUD3') 
       print*,mp_physics, ' not yet a valid choice, exiting'
       stop 'module_wrf'
     CASE ('NCEPCLOUD5') 
       print*,mp_physics, ' not yet a valid choice, exiting'
       stop 'module_wrf'
     CASE DEFAULT
        print*,'No recognized microphysics selected, will not use'  
     END SELECT micro_select

! initialize a nudging (to control) function if asked
   IF ( trim(nudge_f_type) == 'WRF' ) then
     DO k = 1,nz
       nudge_func(k) = max(  0.0,                               &
       (nudge_f_hwpress**2-(p_phy(1,k,1)-p_phy(1,nz,1))**2)/    &
       (nudge_f_hwpress**2+(p_phy(1,k,1)-p_phy(1,nz,1))**2) )
     ENDDO
   ENDIF

END SUBROUTINE INIT_WRF

!------------------------------------------------------------------
SUBROUTINE WRF(seconds,days,jday,seconds_in_day)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: days, seconds
   INTEGER, INTENT(IN) :: jday, seconds_in_day
   LOGICAL             :: debug1, debug2, debug3, debug4, debug5


   debug1 = .false.  ! before sfc
   debug2 = .false.  ! before pbl
   debug3 = .false.  ! after sfc
   debug4 = .false.  ! after pbl
   debug5 = .false.  ! after pbl and tendencies applied

! A few kludges to get it working in DART restart mode
   LOWLYR = 1 ! required here only for SCM restart

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
     
!----------------------------------------------------------
! Set up all the external forcing if using 
!----------------------------------------------------------

      IF (init_f) THEN

!   error checking to make sure DART does not over shoot the forcing
         IF ( time > start_forecast + forecast_length ) THEN
            PRINT*,"You do not have a long enough forcing series to "
            PRINT*,"reach the end of the observation sequence."
            PRINT*,time,start_forecast + forecast_length
            STOP "module_wrf"
         ENDIF

         DO i=ims,ime
            DO j=jms,jme
               PSFC(I,J)=p8w(I,kts,J)
               pi_phy(i,:,j)=(p_phy(i,:,j)/p1000mb)**rcp
            ENDDO
         ENDDO
         
!--------------------------------------------------------------------
! temporal interpolation weights
!--------------------------------------------------------------------
         imin = 1+INT(REAL(time-start_forecast)/REAL(splineinterval))
         imin = MIN(imin,nsplinetimes)
         imax = MIN(imin+1,nsplinetimes)
         fract = (time-splinetimes(imin))/REAL(splineinterval)

!-------------------------------------------------------------------
! geostrophic wind
!-------------------------------------------------------------------
         u_g(:)=u_g_f(:,imin)+&
              &(u_g_f(:,imax)-u_g_f(:,imin))*fract
         v_g(:)=v_g_f(:,imin)+&
              &(v_g_f(:,imax)-v_g_f(:,imin))*fract

!-------------------------------------------------------------------
! precip
!-------------------------------------------------------------------

         imin_smos=1+INT(REAL(time-start_forecast)/REAL(splineinterval_smos))
         imin_smos = min(imin_smos,nsplinetimes_smos)
         imax_smos = min(imin_smos+1,nsplinetimes_smos)


         ! put in forcing precip only if not running microphysics
         if ( trim(mp_physics) == 'NONE' ) then
            RAINCV=0.
            RAINNCV=precip_f(imax_smos)*dt/REAL(splineinterval_smos)
         endif

!-------------------------------------------------------------------
! soil state
!-------------------------------------------------------------------

         if ( force_soil ) then
           imin_soil=1+INT(REAL(time-start_forecast)/REAL(splineinterval_soil))
           imin_soil = min(imin_soil,nsplinetimes_soil)
           imax_soil = min(imin_soil+1,nsplinetimes_soil)
           fract_soil = (time-splinetimes_soil(imin_soil))/REAL(splineinterval_soil)
           tslb(1,:,1)=tslb_f(:,imin_soil)+&
              (tslb_f(:,imax_soil)-tslb_f(:,imin_soil))*fract_soil
           smois(1,:,1)=smois_f(:,imin_soil)+&
              (smois_f(:,imax_soil)-smois_f(:,imin_soil))*fract_soil
           sh2o = smois ! all liquid
         endif

!--------------------------------------------------------------------
! nudge advective tendencies if requested (follow Ghan et al 1999)
!--------------------------------------------------------------------

         IF ( t_advection .or. qv_advection .or. u_advection ) THEN

            imin_adv=1+INT(REAL(time-start_forecast)/REAL(splineinterval_advection))
            imin_adv = MIN(imin_adv,nsplinetimes_advection)
            imax_adv = MIN(imin_adv+1,nsplinetimes_advection)
            fract_adv = (time-splinetimes_advection(imin_adv))/REAL(splineinterval_advection)

               u_upstream_x(:)=u_upstream_x_f(:,imin_adv)+&
                 (u_upstream_x_f(:,imax_adv)-u_upstream_x_f(:,imin_adv))*fract_adv
               u_upstream_y(:)=u_upstream_y_f(:,imin_adv)+&
                 (u_upstream_y_f(:,imax_adv)-u_upstream_y_f(:,imin_adv))*fract_adv
               v_upstream_x(:)=v_upstream_x_f(:,imin_adv)+&
                 (v_upstream_x_f(:,imax_adv)-v_upstream_x_f(:,imin_adv))*fract_adv
               v_upstream_y(:)=v_upstream_y_f(:,imin_adv)+&
                 (v_upstream_y_f(:,imax_adv)-v_upstream_y_f(:,imin_adv))*fract_adv


            tau_u(:) = dx/abs(u_upstream_x(:))
            tau_v(:) = dx/abs(v_upstream_y(:))


            IF ( t_advection ) THEN

               th_upstream_x(:)=th_upstream_x_f(:,imin_adv)+&
                 (th_upstream_x_f(:,imax_adv)-th_upstream_x_f(:,imin_adv)) &
                   *fract_adv
               th_upstream_y(:)=th_upstream_y_f(:,imin_adv)+&
                 (th_upstream_y_f(:,imax_adv)-th_upstream_y_f(:,imin_adv)) &
                   *fract_adv

               DO i=ims,ime
                  DO j=jms,jme
                     DO k=kts,kte+1
                        IF ( tau_u(k) /= -9999.0 .and. tau_v(k) /= -9999.0 &
                            .and. th_upstream_x(k) /= -9999.0             &
                            .and. th_upstream_y(k) /= -9999.0) THEN
                     
                            th_phy(i,k,j) = th_phy(i,k,j) &
                            +  (dt*(th_upstream_x(k) - th_phy(i,k,j))/tau_u(k) &
                            +  dt*(th_upstream_y(k) - th_phy(i,k,j))/tau_v(k)) &
                            *  tadvection_scale(i,j) 

                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

            IF ( qv_advection ) THEN

               qv_upstream_x(:)=qv_upstream_x_f(:,imin_adv)+&
                 (qv_upstream_x_f(:,imax_adv)-qv_upstream_x_f(:,imin_adv)) &
                  *fract_adv
               qv_upstream_y(:)=qv_upstream_y_f(:,imin_adv)+&
                 (qv_upstream_y_f(:,imax_adv)-qv_upstream_y_f(:,imin_adv)) &
                  *fract_adv

               IF ( P_QC > 1 .and. qc_advection ) THEN
                  qc_upstream_x(:)=qc_upstream_x_f(:,imin_adv)+&
                    (qc_upstream_x_f(:,imax_adv)-qc_upstream_x_f(:,imin_adv)) &
                     *fract_adv
                  qc_upstream_y(:)=qc_upstream_y_f(:,imin_adv)+&
                    (qc_upstream_y_f(:,imax_adv)-qc_upstream_y_f(:,imin_adv)) &
                     *fract_adv
               ENDIF

               IF ( P_QR > 1 .and. qr_advection ) THEN
                  qr_upstream_x(:)=qr_upstream_x_f(:,imin_adv)+&
                    (qr_upstream_x_f(:,imax_adv)-qr_upstream_x_f(:,imin_adv)) &
                     *fract_adv
                  qr_upstream_y(:)=qr_upstream_y_f(:,imin_adv)+&
                    (qr_upstream_y_f(:,imax_adv)-qr_upstream_y_f(:,imin_adv)) &
                     *fract_adv
               ENDIF

               IF ( P_QI > 1 .and. qi_advection ) THEN
                  qi_upstream_x(:)=qi_upstream_x_f(:,imin_adv)+&
                    (qi_upstream_x_f(:,imax_adv)-qi_upstream_x_f(:,imin_adv)) &
                     *fract_adv
                  qi_upstream_y(:)=qi_upstream_y_f(:,imin_adv)+&
                    (qi_upstream_y_f(:,imax_adv)-qi_upstream_y_f(:,imin_adv)) &
                     *fract_adv
               ENDIF

               IF ( P_QG > 1 .and. qg_advection ) THEN
                  qg_upstream_x(:)=qg_upstream_x_f(:,imin_adv)+&
                    (qg_upstream_x_f(:,imax_adv)-qg_upstream_x_f(:,imin_adv)) &
                     *fract_adv
                  qg_upstream_y(:)=qg_upstream_y_f(:,imin_adv)+&
                    (qg_upstream_y_f(:,imax_adv)-qg_upstream_y_f(:,imin_adv)) &
                     *fract_adv
               ENDIF

               DO i=ims,ime
                  DO j=jms,jme
                     DO k=kts,kte+1
                        IF ( tau_u(k) /= -9999.0 .and. tau_v(k) /= -9999.0 &
                             .and. qv_upstream_x(k) /= -9999.0  &
                             .and. qv_upstream_y(k) /= -9999.0) THEN

                             moist(i,k,j,P_QV) = moist(i,k,j,P_QV) &
                       +  (dt*(qv_upstream_x(k) - moist(i,k,j,P_QV))/tau_u(k)&
                       +  dt*(qv_upstream_y(k) - moist(i,k,j,P_QV))/tau_v(k))&
                       *  qadvection_scale(i,j)
                        ENDIF

                        IF ( tau_u(k) /= -9999.0 .and. tau_v(k) /= -9999.0 &
                             .and. P_QC > 1 .and. qc_advection ) THEN

                             moist(i,k,j,P_QC) = moist(i,k,j,P_QC) &
                       +  (dt*(qc_upstream_x(k) - moist(i,k,j,P_QC))/tau_u(k)&
                       +  dt*(qc_upstream_y(k) - moist(i,k,j,P_QC))/tau_v(k))&
                       *  qadvection_scale(i,j)
                          
                        ENDIF

                        IF ( tau_u(k) /= -9999.0 .and. tau_v(k) /= -9999.0 &
                             .and. P_QR > 1  .and. qr_advection) THEN

                             moist(i,k,j,P_QR) = moist(i,k,j,P_QR) &
                       +  (dt*(qr_upstream_x(k) - moist(i,k,j,P_QR))/tau_u(k)&
                       +  dt*(qr_upstream_y(k) - moist(i,k,j,P_QR))/tau_v(k))&
                       *  qadvection_scale(i,j)
                          
                        ENDIF

                        IF ( tau_u(k) /= -9999.0 .and. tau_v(k) /= -9999.0 &
                             .and. P_QI > 1 .and. qi_advection ) THEN

                             moist(i,k,j,P_QI) = moist(i,k,j,P_QI) &
                       +  (dt*(qi_upstream_x(k) - moist(i,k,j,P_QI))/tau_u(k)&
                       +  dt*(qi_upstream_y(k) - moist(i,k,j,P_QI))/tau_v(k))&
                       *  qadvection_scale(i,j)
                          
                        ENDIF

                        IF ( tau_u(k) /= -9999.0 .and. tau_v(k) /= -9999.0 &
                             .and. P_QG > 1 .and. qg_advection ) THEN

                             moist(i,k,j,P_QG) = moist(i,k,j,P_QG) &
                       +  (dt*(qg_upstream_x(k) - moist(i,k,j,P_QG))/tau_u(k)&
                       +  dt*(qg_upstream_y(k) - moist(i,k,j,P_QG))/tau_v(k))&
                       *  qadvection_scale(i,j)
                          
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF ! q advection

            IF ( u_advection ) THEN

               u_upstream_x(:)=u_upstream_x_f(:,imin_adv)+&
                 (u_upstream_x_f(:,imax_adv)-u_upstream_x_f(:,imin_adv))*fract_adv
               u_upstream_y(:)=u_upstream_y_f(:,imin_adv)+&
                 (u_upstream_y_f(:,imax_adv)-u_upstream_y_f(:,imin_adv))*fract_adv
               v_upstream_x(:)=v_upstream_x_f(:,imin_adv)+&
                 (v_upstream_x_f(:,imax_adv)-v_upstream_x_f(:,imin_adv))*fract_adv
               v_upstream_y(:)=v_upstream_y_f(:,imin_adv)+&
                 (v_upstream_y_f(:,imax_adv)-v_upstream_y_f(:,imin_adv))*fract_adv

               DO i=ims,ime
                  DO j=jms,jme
                     DO k=kts,kte+1
                        IF ( tau_u(k) /= -9999.0 .and. tau_v(k) /= -9999.0 )THEN

                           u_phy(i,k,j) = u_phy(i,k,j) &
                          +  (dt*(u_upstream_x(k) - u_phy(i,k,j))/tau_u(k) &
                          +  dt*(u_upstream_y(k) - u_phy(i,k,j))/tau_v(k)) &
                          *  uadvection_scale(i,j)

                           v_phy(i,k,j) = v_phy(i,k,j) &
                          +  (dt*(v_upstream_x(k) - v_phy(i,k,j))/tau_u(k) &
                          +  dt*(v_upstream_y(k) - v_phy(i,k,j))/tau_v(k)) &
                          *  vadvection_scale(i,j)

                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF ! u advection
         ENDIF ! any advection

!--------------------------------------------------------
! This is done here to accommodate the advection and to correct
! anything funny that may have happened during DA
!--------------------------------------------------------
         ! update moist, t and rho
         where (moist(:,:,:,P_QV) < 1e-20) moist(:,:,:,P_QV) = 1e-20
         t_phy=th_phy*pi_phy
         rho=p_phy/(r_d*t_phy*(1.+SVP1*moist(:,:,:,P_QV)))

      ENDIF ! end of main forcing block

!--------------------------------------------------------
! Main radiation block, can be run in either forced or ideal mode
!--------------------------------------------------------
 
      SELECT CASE (trim(ra_type))
! Options for radiation ...RT
      CASE ('INT') ! interactive full column radiative transfer

! safeguard for radiation calculations...top required to be "high enough"...
         IF (z_grid(nz).LT.7000.0) THEN
            print*, '********************************************* '
            print*, '********************************************* '
            print*, 'Top of grid too low for interactive radiation'
            print*, 'Stopping simulation'
            print*, '********************************************* '
            print*, '********************************************* '
            STOP
         ENDIF

         IF (first_time) then

! read solar attenuation file
! and find appropriate value / to detmine the "solar constant" at the top of the model...RT

            OPEN(solattunit,file='ARMSGP_SolarClearAirAtt.dat')
            read(solattunit, *) solattjdayref
            DO k=1,32
               read(solattunit, '(F4.1, 13F12.9)') solatt(k,1), &
  &                 (solatt(k,m), m=2,14)
            ENDDO
            CLOSE(solattunit)

            first_time=.FALSE.
               
! bi-linear interpolation (julian day and height of top level)

            k = 1
            m = 1
            do while (solattjdayref(m).LT.jday)
               m = m + 1
            enddo
            do while (solatt(k,1).LT.(z_grid(nz)/1000.0))
               k = k + 1
            enddo
               
            t = (REAL(jday) - solattjdayref(m-1))/(solattjdayref(m) - solattjdayref(m-1))
            u = (z_grid(nz)/1000.0) - solatt(k-1,1)/(solatt(k,1) - solatt(k-1,1))
               
            solatttop = (1-t)*(1-u)*solatt(k,m) + t*(1-u)*solatt(k,m+1) &
            & + t*u*solatt(k-1,m+1) + (1-t)*u*solatt(k-1,m)
               
!               solatttop = 0.0
               
         ENDIF
            
!   ========  need to calculate (interpolate) t8w from updated t_phy... ============ 
!   ========  what about humidity and water variables ????              ============ 

         var(1) = tsk(1,1)
         z_var(1) = 0.
         DO k=2,nz+1
            z_var(k) = z(1,k-1,1)
            var(k) = t_phy(1,k-1,1)
         ENDDO

         DO i=ims,ime
            DO j=jms,jme
               DO k=kms,kme
                  t8w(i,k,j)=linear(z_var,var,nz+1,z8w(i,k,j))
                  rho(i,k,j)=p_phy(i,k,j)/(r_d* &
                  & t_phy(i,k,j)*(1.+SVP1*moist(i,k,j,P_QV)))
               ENDDO
            ENDDO
         ENDDO

!-------------------------------------------------------------------
! either interpolate radiation in time or compute it with parameterization
!-------------------------------------------------------------------

         julian = REAL(jday)

         if ( .not. fixed_timeofday ) then
           xtime = (time - start_forecast)/60. ! time since start of simul. (in min)
         else
           xtime = timeofday_ref*60.0 - start_forecast/60.
         endif

! STOCHASTIC CLOUD ...RT

         IF ((forc_stochastic_cloud).and.((xtime/60.).gt.spin_up_period)) then !no clouds during a spin-up period

            CALL cloud_stochastic(MOIST,p_qc,icloud,z,PBLH,rho,t_phy, &
            & p_phy,stoch_yn,stoch_cbh,stoch_lwp,itimestep,ntime,   &
            & ims,ime,jms,jme,kms,kme,n_moist,indexrealization)


         ENDIF
            
         CALL radiation_driver( ITIMESTEP=itimestep                           &
     &           ,JULDAY=julday       ,JULIAN=julian      ,JULYR=julyr        &
     &           ,GMT=gmt             ,DT=dt              ,XTIME=xtime        &
     &           ,XLAT=(/lat_ref/)    ,XLONG=(/lon_ref/)                      &
     &           ,RADT=radt           ,RA_CALL_OFFSET=ra_call_offset          &
     &           ,STEPRA=stepra       ,NPHS=1     ,soltopfrac=solatttop       &
     &           ,RA_LW_PHYSICS=ra_lw_physics, RA_SW_PHYSICS=ra_sw_physics    &
     &           ,ALBEDO=albedo       ,EMISS=emiss                            &
     &           ,XICE=xice           ,XLAND=xland        ,VEGFRA=vegfra      &
     &           ,SNOW=snow           ,WARM_RAIN=warm_rain                    &
     &           ,TSK=tsk             ,T8W=t8w            ,T=t_phy            &
     &           ,P8W=p8w             ,P=p_phy            ,PI=pi_phy          &
     &           ,RHO=rho                                                     &
     &           ,Z=z                 ,DZ8W=dz8w                              &
     &           ,HBOT=hbot           ,HTOP=htop          ,HBOTR=hbotr        &
     &           ,HTOPR=htopr         ,ICLOUD=icloud                          &
     &           ,GLW=glw             ,GSW=gsw                                &
     &           ,RSWTOA=rswtoa       ,RLWTOA=rlwtoa                          &  
     &           ,RTHRATEN=rthraten   ,RTHRATENLW=rthratenlw                  &
     &           ,RTHRATENSW=rthratensw                                       &
     &           ,SWDOWN=swdown       ,SWDOWNC=swdownc                        &
     &           ,TAUCLDC=taucldc     ,TAUCLDI=taucldi                        &
     &           ,ACFRCV=acfrcv       ,ACFRST=acfrst                          &
     &           ,CFRACH=cfrach       ,CFRACL=cfracl      ,CFRACM=cfracm      &
     &           ,CUPPT=cuppt         ,CZMEAN=czmean                          &
     &           ,NCFRCV=ncfrcv       ,NCFRST=ncfrst                          &
!     Optional urban
!     &           ,DECLIN_URB=0.   ,COSZ_URB2D=(/0.,0./) ,OMG_URB2D=(/0.,0./)  &
!     Optional CAM
     &           ,LEVSIZ=levsiz, N_OZMIXM=n_ozmixm                            &
     &           ,N_AEROSOLC=n_aerosolc                                       &
     &           ,PAERLEV=paerlev                                             &
     &           ,CAM_ABS_DIM1=cam_abs_dim1, CAM_ABS_DIM2=cam_abs_dim2        &
     &           ,CAM_ABS_FREQ_S=cam_abs_freq_s                               &
!     indexes
     &           ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde           &
     &           ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme           &
     &           ,I_START=i_start,I_END=min(i_end, ide)                       &
     &           ,J_START=j_start,J_END=min(j_end, jde)                       &
     &           ,KTS=k_start(1), KTE=kde-1                                   &
     &           ,num_tiles=num_tiles                                         &
!     Optional                          
     &           , CLDFRA=cldfra                                              &
     &           , pb=pb                                                      &
     &           , F_ICE_PHY=f_ice_phy,F_RAIN_PHY=f_rain_phy                  &
     &           , QV=moist(ims,kms,jms,P_QV), F_QV=F_QV                      &
     &           , QC=moist(ims,kms,jms,P_QC), F_QC=F_QC                      &
     &           , QR=moist(ims,kms,jms,P_QR), F_QR=F_QR                      &
     &           , QI=moist(ims,kms,jms,P_QI), F_QI=F_QI                      &
     &           , QS=moist(ims,kms,jms,P_QS), F_QS=F_QS                      &
     &           , QG=moist(ims,kms,jms,P_QG), F_QG=F_QG                      &
#ifdef ACFLUX
     &           ,ACSWUPT=acswupt    ,ACSWUPTC=acswuptc                       &
     &           ,ACSWDNT=acswdnt    ,ACSWDNTC=acswdntc                       &
     &           ,ACSWUPB=acswupb    ,ACSWUPBC=acswupbc                       &
     &           ,ACSWDNB=acswdnb    ,ACSWDNBC=acswdnbc                       &
     &           ,ACLWUPT=aclwupt    ,ACLWUPTC=aclwuptc                       &
     &           ,ACLWDNT=aclwdnt    ,ACLWDNTC=aclwdntc                       &
     &           ,ACLWUPB=aclwupb    ,ACLWUPBC=aclwupbc                       &
     &           ,ACLWDNB=aclwdnb    ,ACLWDNBC=aclwdnbc                       &
     &           ,SWUPT=swupt        ,SWUPTC=swuptc                           &
     &           ,SWDNT=swdnt        ,SWDNTC=swdntc                           &
     &           ,SWUPB=swupb        ,SWUPBC=swupbc                           &
     &           ,SWDNB=swdnb        ,SWDNBC=swdnbc                           &
     &           ,LWUPT=lwupt        ,LWUPTC=lwuptc                           &
     &           ,LWDNT=lwdnt        ,LWDNTC=lwdntc                           &
     &           ,LWUPB=lwupb        ,LWUPBC=lwupbc                           &
     &           ,LWDNB=lwdnb        ,LWDNBC=lwdnbc                           &
#endif
     &           ,LWCF=lwcf                                                   &
     &           ,SWCF=swcf                                                   &
     &           ,OLR=olr                                                     &
     &           ,OZMIXM=ozmixm, PIN=pin                                      &
     &           ,M_PS_1=m_ps_1, M_PS_2=m_ps_2, AEROSOLC_1=aerosolc_1         &
     &           ,AEROSOLC_2=aerosolc_2, M_HYBI0=m_hybi0                      &
     &           ,ABSTOT=abstot, ABSNXT=absnxt, EMSTOT=emstot                 &
#ifdef WRF_CHEM
     &           ,QC_ADJUST=GD_CLOUD , QI_ADJUST=GD_CLOUD2                    &
     &           ,PM2_5_DRY=pm2_5_dry, PM2_5_WATER=pm2_5_water                &
     &           ,PM2_5_DRY_EC=pm2_5_dry_ec                                   &
     &           ,TAUAER300=tauaer1, TAUAER400=tauaer2                        & 
     &           ,TAUAER600=tauaer3, TAUAER999=tauaer4                        &
     &           ,GAER300=gaer1, GAER400=gaer2, GAER600=gaer3, GAER999=gaer4  &  
     &           ,WAER300=waer1, WAER400=waer2, WAER600=waer3, WAER999=waer4  & 
#endif
                 ,cu_rad_feedback=cu_rad_feedback                             &
     &           )

      CASE ('WRF') ! sfc rad. fluxes from WRF output 

               imin_flux=1+INT(REAL(time-start_forecast)/REAL(splineinterval_flux))
               imin_flux = MIN(imin_flux,nsplinetimes_flux)
               imax_flux = MIN(imin_flux+1,nsplinetimes_flux)
               fract_flux = (time-splinetimes_flux(imin_flux))/REAL(splineinterval_flux)
               
               glw=glw_f(imin_flux)+(glw_f(imax_flux)-glw_f(imin_flux))*fract_flux
               gsw=gsw_f(imin_flux)+(gsw_f(imax_flux)-gsw_f(imin_flux))*fract_flux

      CASE ('IDL') ! specified sfc fluxes
         if ( fixed_timeofday ) then
           ! account for GMT offset
           glw=glwfunc((timeofday_ref-6)*3600.0)
           gsw=gswfunc((timeofday_ref-6)*3600.0,albedo(1,1))
         else
           glw=glwfunc(time)
           gsw=gswfunc(time,albedo(1,1))
         endif
         ! put in forcing precip only if not running microphysics
         if ( trim(mp_physics) == 'NONE' ) then
            RAINCV=0.
            RAINNCV=prate_ref*dt
         endif

      CASE DEFAULT 
         print*,ra_type,' not a valid radiation option, exiting'
         stop
      END SELECT ! parameterized radiation, specified from WRF, or ideal


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


!print*,tslb(1,1,1),tslb_f(1,imin_soil:imax_soil)
!print*,smois(1,1,1),smois_f(1,imin_soil:imax_soil)
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
     &        ,RA_LW_PHYSICS_I=ra_lw_physics_i        ,RHO=rho            &
     &        ,RMOL=rmol          ,SFCEVP=sfcevp      ,SFCEXC=sfcexc      &
     &        ,SFCRUNOFF=sfcrunoff                                        &
     &         ,BL_PBL_PHYSICS=bl_pbl_physics&
     &        ,SF_SFCLAY_PHYSICS=sf_sfclay_physics                        &
     &        ,SF_SURFACE_PHYSICS=sf_surface_physics  ,SH2O=sh2o          &
     &        ,SHDMAX=shdmax      ,SHDMIN=shdmin      ,SMOIS=smois        &
     &        ,SMSTAV=smstav      ,SMSTOT=smstot      ,SNOALB=snoalb      &
     &        ,SNOW=snow          ,SNOWC=snowc        ,SNOWH=snowh        &
     &        ,SST=sst_input         ,SST_UPDATE=sst_update               &
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
     &        ,DECLIN_URB=0.  ,COSZ_URB2D=(/0.,0./)                       & !I urban
     &        ,OMG_URB2D=(/0.,0./)    ,xlat_urb2d=(/0.,0./)               & !I urban
     &        ,NUM_ROOF_LAYERS=0                                          & !I urban
     &        ,NUM_WALL_LAYERS=0                                          & !I urban
     &        ,NUM_ROAD_LAYERS=0                                          &
     &        ,DZR=DZR ,DZB=DZB ,DZG=DZG                                  & !I urban
     &        ,TR_URB2D=(/0.,0./) ,TB_URB2D=(/0.,0./)                     &
     &        ,TG_URB2D=(/0.,0./)                                         & !H urban
     &        ,TC_URB2D=(/0.,0./) ,QC_URB2D=(/0.,0./)                     & !H urban
     &        ,UC_URB2D=(/0.,0./)                                         & !H urban
     &        ,XXXR_URB2D=(/0.,0./)                                       &
     &        ,XXXB_URB2D=(/0.,0./)                                       & !H urban
     &        ,XXXG_URB2D=(/0.,0./)                                       &
     &        ,XXXC_URB2D=(/0.,0./)                                       & !H urban
     &        ,TRL_URB3D=TRL_URB3D   ,TBL_URB3D=TBL_URB3D                 & !H urban
     &        ,TGL_URB3D=TBL_URB3D                                        & !H urban
     &        ,SH_URB2D=(/0.,0./)     ,LH_URB2D=(/0.,0./)                 &
     &        ,G_URB2D=(/0.,0./)                                          & !H urban
     &        ,RN_URB2D=(/0.,0./)     , TS_URB2D=(/0.,0./)                & !H urban
     &        ,FRC_URB2D=(/0.,0./)                                        & !H urban
     &        ,UTYPE_URB2D=(/0,0/)                                        & !H urban
     &        ,ucmcall=ucmcall                                            & !H urban

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
     &        ,SMFR3D=smfr3d,KEEPFR3DFLAG=keepfr3dflag,                   & ! ruc lsm

     &ts_ref=ts_ref,ps_ref=ts_ref,dtamplitude_ref=dtamplitude_ref,&
     &mavail_ref=mavail_ref,time=time,&
     &bucket_model=bucket_model,maxm=maxm,minm=minm,erate=erate,&
     &calc_ust=calc_ust,z_o=z_o,z_t=z_t,z_q=z_q&
     &                                                              )

! compute energy balance (positive down)
!     print*,'albedo, emiss, stbolt ',albedo,emiss,stbolt
!     print*,'tsk ',tsk
!     print*,'components ',gsw,-albedo*gsw,emiss*glw,-emiss*stbolt*tsk**4,-lh,-hfx,grdflx
!     print*,'residual ',(1-albedo)*gsw + emiss*glw - emiss*stbolt*tsk**4 - lh - hfx - grdflx
!
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
     &        ,BL_PBL_PHYSICS=bl_pbl_physics                              &
     &        ,sf_surface_physics=sf_surface_physics                      &
     &        ,BR=br              ,CHKLOWQ=chklowq    ,CT=ct              &
     &        ,DT=dt              ,DX=dx              ,DZ8W=dz8w          &
     &        ,EL_MYJ=el_myj      ,EXCH_H=exch_h      ,GRDFLX=grdflx      &
     &        ,GZ1OZ0=gz1oz0      ,HFX=hfx            ,HT=ht              &
     &        ,ITIMESTEP=itimestep                    ,KPBL=kpbl          &
     &        ,LH=lh              ,LOWLYR=lowlyr      ,P8W=p8w            &
     &        ,PBLH=pblh          ,PI_PHY=pi_phy      ,PSIH=psih          &
     &        ,PSIM=psim          ,P_PHY=p_phy        ,QFX=qfx            &
     &        ,QSFC=qsfc          ,QZ0=qz0                                &
     &        ,RA_LW_PHYSICS_I=ra_lw_physics_i                            &
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
           &  ,uflux=uflux,vflux=vflux,hflux=hflux,qflux=qflux, &
           &   k_t=k_t,k_m=k_m&
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
            DO k=kts,kte+1
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
               th_phy(i,k,j)=th_phy(i,k,j)+(RTHBLTEN(i,k,j)+RTHRATEN(i,k,j))*dt !added tendendcy from radiation...RT
               moist(i,k,j,P_QV)=moist(i,k,j,P_QV)+RQVBLTEN(i,k,j)*dt
               rho(i,k,j)=p_phy(i,k,j)/(r_d*&
                    &t_phy(i,k,j)*(1.+SVP1*moist(i,k,j,P_QV)))
            ENDDO
         ENDDO
      ENDDO

      ! update t
      t_phy=th_phy*pi_phy
      
! now nudge state to forcing data set
      IF ( trim(nudge_f_type) == 'WRF' ) then
         u_fi(:)=u_f(:,imin)+(u_f(:,imax)-u_f(:,imin))*fract
         u_phy(1,:,1) = u_phy(1,:,1) +           &
                   nudge_f_coeff*nudge_func*(u_fi-u_phy(1,:,1))

         v_fi(:)=v_f(:,imin)+(v_f(:,imax)-v_f(:,imin))*fract
         v_phy(1,:,1) = v_phy(1,:,1) +           &
                   nudge_f_coeff*nudge_func*(v_fi-v_phy(1,:,1))

         t_fi(:)=t_f(:,imin)+(t_f(:,imax)-t_f(:,imin))*fract
         t_phy(1,:,1) = t_phy(1,:,1) +           &
                   nudge_f_coeff*nudge_func*(t_fi-t_phy(1,:,1))

         p_fi(:)=p_f(:,imin)+(p_f(:,imax)-p_f(:,imin))*fract
         p_phy(1,:,1) = p_phy(1,:,1) +           &
                   nudge_f_coeff*nudge_func*(p_fi-p_phy(1,:,1))

         q_fi(:)=q_f(:,imin)+(q_f(:,imax)-q_f(:,imin))*fract
         moist(1,:,1,P_QV) = moist(1,:,1,P_QV) + &
                   nudge_f_coeff*nudge_func*(q_fi-moist(1,:,1,P_QV))
! update theta, etc
         pi_phy(1,:,1)=(p_phy(1,:,1)/p1000mb)**rcp
         th_phy=t_phy/pi_phy
      ELSEIF ( trim(nudge_f_type) == 'OBS' ) then
         print*,'Not ready to nudge to obs yet'
         stop 'module_wrf'
      ELSE
      ENDIF

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

! microphysics is time split in the WRF.  The tendencies are
! not explicitly computed here as in WRF.  Rather, we just let th 
! vary here.  In these dynamics is it identical.

  IF ( mp_physics /= 'NONE' ) then

             CALL microphysics_driver(                                    &
     &         DT=dt              ,DX=dx              ,DY=dx              &
     &        ,DZ8W=dz8w          ,F_ICE_PHY=f_ice_phy                    &
     &        ,ITIMESTEP=itimestep                    ,LOWLYR=lowlyr      &
     &        ,P8W=p8w            ,P=p_phy            ,PI_PHY=pi_phy      &
     &        ,RHO=rho            ,SPEC_ZONE=spec_zone                    &
     &        ,SR=sr              ,TH=th_phy                              &
     &        ,WARM_RAIN=warm_rain                    ,XLAND=xland        &
     &        ,SPECIFIED=specified_bdy, CHANNEL_SWITCH=channel_bdy        &
     &        ,F_RAIN_PHY=f_rain_phy                                      &
     &        ,F_RIMEF_PHY=f_rimef_phy                                    &
     &        ,MP_PHYSICS=mp_physics                                      &
     &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde          &
     &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme          &
     &        ,I_START=i_start,I_END=min(i_end, ide)                      &
     &        ,J_START=j_start,J_END=min(j_end, jde)                      &
     &        ,KTS=k_start(1), KTE=kde-1                                     &
     &        ,NUM_TILES=num_tiles                                        &
                 ! Optional
     &        , RAINNC=rainnc, RAINNCV=rainncv                            &
     &        , W=w_phy, Z=z, HT=ht                                       &
     &        , MP_RESTART_STATE=mp_restart_state                         &
     &        , TBPVS_STATE=tbpvs_state                                   & !  etampnew
     &        , TBPVS0_STATE=tbpvs0_state                                 & !  etampnew
     &        , QV_CURR=moist(ims,kms,jms,P_QV), F_QV=F_QV              &
     &        , QC_CURR=moist(ims,kms,jms,P_QC), F_QC=F_QC              &
     &        , QR_CURR=moist(ims,kms,jms,P_QR), F_QR=F_QR              &
     &        , QI_CURR=moist(ims,kms,jms,P_QI), F_QI=F_QI              &
     &        , QS_CURR=moist(ims,kms,jms,P_QS), F_QS=F_QS              &
     &        , QG_CURR=moist(ims,kms,jms,P_QG), F_QG=F_QG              &
     &        , QNI_CURR=scalar(ims,kms,jms,P_QNI), F_QNI=F_QNI         &
     &        , QT_CURR=scalar(ims,kms,jms,P_QT), F_QT=F_QT             &
                                                                          )

     ! update t
     t_phy=t_phy*pi_phy
  ENDIF ! call microphysics?

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
     IF ( outfinterval == 0.0 .or. MOD(time-start_forecast,outfinterval) < epsilon) THEN
        WRITE(ncunit,'(20f16.7)')time/3600.
        print*,'Writing output ',time
        where(z > 20000.0)      z = -9999.0
        where(u_phy > 20000.0)  u_phy = -9999.0
        where(v_phy > 20000.0)  v_phy = -9999.0
        where(th_phy > 20000.0) th_phy = -9999.0
        where(moist > 20000.0)  moist = -9999.0
        WRITE(ncunit,'(20f16.7)'),0.0,u10(1,1),v10(1,1),t2(1,1),q2(1,1)
        DO k=kts,kte-1
           WRITE(ncunit,'(20f16.7)')z(1,k,1),u_phy(1,k,1),&
                &v_phy(1,k,1),t_phy(1,k,1),moist(1,k,1,P_QV)
        ENDDO
     ENDIF
  ENDIF


  IF (MOD(NINT(time-timeo),outinterval)==0) THEN 
     WRITE(isfcunit,'(30f16.7)')time/3600.,tsk(1,1),tmn(1,1),qsfc(1,1),&
          &ust,hfx,xlv*qfx,gsw(1,1),glw(1,1),tslb(1,:,1),SMOIS(1,:,1)
  ENDIF

  IF (itimestep==ntime) THEN

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


! FESOM 2 (Finite-volumE Sea ice-Ocean Model)
! multi-resolution ocean general circulation model
! FESOM/fesom2 is licensed under the GNU General Public License v2.0
! Copyright (C) 2018  FESOM team
!
! This module was constructed from pieces of the FESOM V1.4 modules and
! extended to work with DART.

module fesom_modules

  ! Bits from DART
  use utilities_mod, only: get_unit, do_output

  implicit none
  save

  ! *** Modelname ***
  character(5)                 :: runid                ! a model/setup name
  logical                      :: ensemble_run=.true.
  character(5)                 :: ensid='ENS01'


  namelist /modelname/ runid,ensemble_run,ensid


  ! *** time step ***
  integer                      :: step_per_day=12               !number of steps per day
  integer                      :: run_length=1                  !run length
  character                    :: run_length_unit='y'           !unit: y, d, s
  integer                      :: day2ext
  character(4)                 :: runyear='2008'

  namelist /timestep/ step_per_day, run_length, run_length_unit, runyear

  ! *** time series *
  integer                       :: iniday=1
  integer                       :: endday=1

  namelist /timeseries/ IniDay, EndDay

  ! *** Paths for all in and out ***
  character(100)                :: MeshPath='/users/home/ans051/FEOM_PREPROC/mesh-T2G1.5L110b/'
  character(100)                :: OpbndPath='./opbnd/'
  character(100)                :: ClimateDataPath='./hydrography/'
  character(100)                :: ForcingDataPath='./forcing/'
  character(100)                :: TideForcingPath='./tide_forcing/'
  character(100)                :: ResultPath='./result/'

  namelist /paths/  MeshPath, OpbndPath, ClimateDataPath, ForcingDataPath, &
       TideForcingPath, ResultPath

  ! *** ocean climatology data name ***
  character(100)                :: OceClimaDataName='annual_woa01_ts.out'
  logical                       :: use_prepared_init_ice=.false.     !how to initial. ice at the beginning 

  namelist /initialization/ OceClimaDataName, use_prepared_init_ice

  ! *** in out ***
  character*4                   :: restartflag='last'                !restart from which saved record,'#','last'
  integer                       :: output_length=1                   !valid for d,h,s
  character                     :: output_length_unit='m'            !output period: y, m, d, h, s 
  integer                       :: logfile_outfreq=1                 !in logfile info. output frequency, # steps

  namelist /inout/ restartflag, output_length, output_length_unit, logfile_outfreq

  ! *** mesh ***
  integer                       :: grid_type=1                       ! z-level, 2 sigma, 3 sigma + z-level

  namelist /mesh_def/ grid_type

  ! *** model geometry
  logical                      :: cartesian=.false.
  logical                      :: fplane=.false.
  logical                      :: betaplane=.false.
  real(kind=8)                 :: f_fplane=-1.4e-4            ![1/s]
  real(kind=8)                 :: beta_betaplane=2.0e-11      ![1/s/m]
  real(kind=8)                 :: domain_length=360.        ![degree]
  !
  logical                      :: rotated_grid=.true.        !option only valid for coupled model case now
  real(kind=8)                 :: alphaEuler=-30.            ![degree] Euler angles, convention:
  real(kind=8)                 :: betaEuler=-90.             ![degree] first around z, then around new x,
  real(kind=8)                 :: gammaEuler=-90.            ![degree] then around new z.

  namelist /geometry/  cartesian, fplane, betaplane, f_fplane, beta_betaplane, &
       domain_length, rotated_grid, alphaEuler, betaEuler, gammaEuler

  ! *** fleap_year ***
  logical                       :: include_fleapyear=.false.
  
!  namelist /calendar/ include_fleapyear
  
   ! *** machine ***
  integer                       :: system=1                    ! XD1 2(byte), HLRN 1(word)
  
  namelist /machine/ system


  ! *** others ***
  real(kind=8)                 :: dt, dt_inv
  integer                      :: istep, nsteps
  integer                       :: n, m, n1, ind, fileID
  integer                       :: save_count
  logical                       :: r_restart

  real(kind=8)             :: timeold, timenew     !time in a day, unit: sec
  integer                  :: dayold, daynew       !day in a year
  integer                  :: yearold, yearnew     !year before and after time step
  integer                  :: month, day_in_month  !month and day in a month
  integer                  :: fleapyear            !1 fleapyear, 0 not 
  integer                  :: ndpyr                !number of days in yearnew 
  integer                  :: num_day_in_month(0:1,12)
  character(4)             :: cyearold, cyearnew   !year as character string      
! character(2)             :: cmonthnew            ! month as character string      
  character(2)             :: cmonth, cday_in_month!month and day in a month as a character string
  data num_day_in_month(0,:) /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  data num_day_in_month(1,:) /31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/


  real(kind=8), parameter      :: pif=3.141592653589793, rad=pif/180.0  
  real(kind=8), parameter      :: omega=2.0*pif/(24.0*60.0*60.0)
  real(kind=8), parameter      :: g=9.81                       ![m/s^2]
  real(kind=8), parameter      :: r_earth=6.3675e6             ![m]
  real(kind=8), parameter      :: rho0=1028.5                  ![kg/m^3]marmara sea deep water density
  real(kind=8), parameter      :: rho0r=1.0/rho0 
  real(kind=8), parameter      :: vcpw=4.2e6                   ![J/m^3/K]volum. heat cap. of water
  real(kind=8), parameter      :: small=1.0e-8                 !small value

  ! *** mixing and friction setting ***
  real                         :: Ah0=6000.                      !Lapl. hori. visc, [m^2/s]
  real                         :: Ahb0=2.7e13                    !Bihar. hori. visc,[m^4/s] <8.0e12/1degree
  real                         :: Kh0=600.                       !lateral diff 
  !
  logical                      :: biharmonic_visc=.false.        !bihar. true will turn off Lapl.
  logical                      :: smagorinsky_visc=.false. 
  !
  logical                      :: increase_equ_zonal_visc=.true. !increase zonal viscosity at equator
  real                         :: fac_visc_increase=3.0          !increase factor
  !
  logical                      :: scale_mixing_h=.true.          !scale hor. diff/visc. coef. by resol.
  integer                      :: scale_mixing_type=2            !1:scale by \delt x^2; 2:scale by \delta x
  real(kind=8)                 :: scalevol=5.e9                  !associated with 100km resolution
  !
  logical                      :: Redi_GM=.true.                 !flag for Redi/GM scheme
  logical                      :: ODM95=.true.                   !taper fcn ODM95
  logical                      :: LDD97=.true.                   !taper fcn LDD97
  real(kind=8)                 :: ratio_K_GM=1.0                 !ratio of K_GM to Kh
  real(kind=8)                 :: Slope_critical=4.0e-3          !upper slope limit for applying Redi/GM
  integer                      :: nslope_version=1               !1-neutral slope over prism,2-over tetrahedra
  !
  real                         :: Av0=1.e-4                  !background (or internal wave) vert. mixing
  real                         :: Kv0=1.e-5                    !m^2/s
  real                         :: visc_conv_limit=0.1          !visc due to convective instability
  real                         :: diff_conv_limit=0.1          !diff due to convective instability
  !
  character(5)                 :: mix_scheme='KPP'        !'KPP','PP', 'MY2p5', 'no'
  !
  real                         :: visc_sh_limit=5.0e-3         !for kpp,max visc due to shear instability
  real                         :: diff_sh_limit=5.0e-3         !for kpp,max diff due to shear instability
  logical                      :: double_diffusion=.true.      !for KPP,dd switch
  logical                      :: smooth_blmc=.true.           !for KPP,hori. smooth of b.l. mixing coeff.
  !
  real                         :: PP_max_mix_coeff=5.0e-3         !for PP, max Kv/Av 
  real                         :: wndmix=1.0e-3                !for PP, to simulate missing high frequency wind
  logical                      :: allow_convect_global=.true.  !for PP, convection for global or only NH
  !
  logical                      :: add_TB04_to_PP=.false.       !for PP, TB04 switch
  real(kind=8)                 :: modiff=0.01                  !for PP, vert. mix. coeff. for TB04
  !
  logical            :: tidal_mixing=.false.        !switch for tidal mixing
  logical            :: use_drag_dissipation=.true.     !barotropic
  logical            :: use_wave_dissipation=.false.    !baroclinic
  logical            :: read_tide_speed=.true.    !read tide speed or use default
  real(kind=8)                 :: default_tide_speed=0.01      !(m/s)
  real(kind=8)                 :: max_drag_diffusivity=5.e-3   !m2/s
  real(kind=8)                 :: max_wave_diffusivity=5.e-3   !m2/s
  character(2)                 :: Tmix_tidalconstituent='M2'   !which tidal constituent 
  character(15)                :: Tmix_tidalmodelname='tpxo71' !source model name
  !
  real(kind=8)                 :: C_d=0.0025                   !Bottom fri. coeff.
  real(kind=8)                 :: aegflx_lat=39.0              !aegean sea flux latitude
  real(kind=8)                 :: blkflx_lat=42.5              !black sea flux latitude
  real(kind=8)                 :: km3yr2m3sec=1e+09/3.1536e+07

  namelist /viscdiff/ Ah0, Ahb0, Kh0, Av0, Kv0, &
       biharmonic_visc, smagorinsky_visc, &
       increase_equ_zonal_visc, fac_visc_increase, &
       scale_mixing_h, scale_mixing_type, scalevol, &
       Redi_GM, ODM95, LDD97, ratio_K_GM, Slope_critical, nslope_version, &
       mix_scheme, visc_sh_limit, diff_sh_limit, visc_conv_limit, diff_conv_limit, &
       double_diffusion, smooth_blmc, PP_max_mix_coeff, wndmix, &
       allow_convect_global, add_TB04_to_PP, modiff, &
       tidal_mixing, use_drag_dissipation, &
       use_wave_dissipation, read_tide_speed, &
       default_tide_speed, max_drag_diffusivity, max_wave_diffusivity, &
       Tmix_tidalconstituent, Tmix_tidalmodelname, &
       C_d

  ! *** surface and open boundary setting ***
  logical                      :: ts_surfbd=.true.     
  !
  logical                      :: ref_sss_local=.false.    !virtual salt flux using local SSS or ref_sss
  real(kind=8)                 :: ref_sss=34.7            !ref. value for cal. virtual salt flux
  !
  real(kind=8)                 :: restore_s_surf=10./(180.*86400.)    ! m/s
  real(kind=8)                 :: restore_t_surf=0.0
  !
  logical                      :: balance_salt_water=.true.    !balance virtual salt or water flux or not
  !
  logical                      :: buffer_zone=.false.
  real(kind=8)                 :: restore_ts_buff= 1./(86400.*5.)         ! timescale for buffer zone [1/s]

  namelist /boundary/ ts_surfbd, ref_sss_local, ref_sss, restore_s_surf, &
       restore_t_surf, balance_salt_water, buffer_zone, restore_ts_buff

  ! *** numerical schemes
  real(KIND=8)                 :: gamma_stab=0.99              !stabilization for ssh
  real(KIND=8)                 :: gamma_stab_nh=0.5            !stabilization for nh-pressure
  real(kind=8)                 :: gamma_fct=0.4                !param for tracer FCT scheme
  real(kind=8)                 :: alpha_AB=1.55                !when Adams-Bashforth Coriolis
  real(kind=8)                 :: alpha_trapez=0.55            !when semi-implicit Coriolis
  real(kind=8)                 :: theta_ssh=0.5        !semi-implicit ssh when semiimpl sheme    
  real(kind=8)                 :: theta_vel=0.5        !semi-implicit baro. vel
  !
  logical                      :: use_vertvisc_impl=.true.     !if implicit vertical viscosity,keep true
  logical                      :: use_vertdiff_impl=.true.     !if implicit vertical diff., keep ture
  logical                      :: use_cori_semi=.false.        !if semiimplicit coriolis force
  !
  logical                      :: lump_uv_matrix=.true.       !for mass vel. matrix case, keep true!
  logical                      :: lump_ts_matrix=.true.       !for mass T/S matrix case, keep true!
  integer                      :: num_iter_solve=3            !iteration # for mass matrix case

  namelist /oce_scheme/ gamma_stab, gamma_stab_nh, gamma_fct, alpha_AB, alpha_trapez, &
       theta_ssh, theta_vel, use_cori_semi

  ! *** density and pressure force ***
  logical                      :: density_linear=.false.
  logical                      :: use_ref_density=.true.     

  namelist /denspress/ density_linear, use_ref_density

  ! *** parameters for nonlinear free surface cases ***
  real(kind=8)                 :: max_ice_loading=5.0        !m, maximal pressure from ice felt by the ocean

  namelist /param_freesurf/ max_ice_loading

  ! *** tide configuration ***
  integer                       :: nmbr_tidal_cons=4
  character(20)                 :: tidal_constituent='M2S2K1O1' !M2 S2 N2 K2 K1 O1 P1 Q1
  character(15)                 :: tidemodelname='tpxo71'
  character(10)                 :: tide_opbnd_type='Flather'     !ssh, vel, Flather 
  real(kind=8)                  :: tide_amplify_coeff=1.0       !amplify tidal amplitude by this factor  

  namelist /tide_obc/ nmbr_tidal_cons, tidal_constituent, tidemodelname, tide_opbnd_type, tide_amplify_coeff
  
  ! *** passive tracer ***
  logical                       :: use_passive_tracer=.false.    !use passive tracer
  integer                       :: num_passive_tracer=1          !only 1 before update
  integer                       :: ptr_start_year=1948           !when to start having ptr
  logical                       :: passive_tracer_flux=.false.   !ptr enters from sfc flux
  logical                       :: passive_tracer_restore=.true. !add ptr from restoring 
  logical                       :: ptr_restore_in_volume=.true.  !restoring in the 3D region
  real(kind=8)                  :: ptr_background_value=0.0      !ptr init. background value
  real(kind=8)                  :: ptr_restore_value=1.0         !restore value for ptr
 
  namelist /passive_tracer/ use_passive_tracer, num_passive_tracer, &
       ptr_start_year, passive_tracer_flux, passive_tracer_restore, &
       ptr_restore_in_volume, ptr_background_value, ptr_restore_value
  
  ! *** age tracer ***
  logical                       :: use_age_tracer=.false.        !use age tracer
  integer                       :: num_age_tracer=1              !only 1 before update
  integer                       :: age_tracer_start_year=1948    !when to start having age tr.
  logical                       :: zero_age_at_release=.true.    !keep zero age in rel. zone
  logical                       :: age_release_in_volume=.false. !age tr. release in 3D regions
  real(kind=8)                  :: age_tracer_restore_time=864000.!restore time scale in rel. zone

  namelist /age_tracer/ use_age_tracer, age_release_in_volume, zero_age_at_release, &
       age_tracer_restore_time, num_age_tracer, age_tracer_start_year       

  ! *** others ***
  integer                       :: num_tracer

  !solver index
  integer                      :: solve_u=1
  integer                      :: solve_v=2
  integer                      :: solve_w=6
  integer                      :: solve_ssh=7
  integer                      :: solve_nhp=8
  integer                      :: solve_tra=10
  logical                      :: iter_first=.true.
  logical                      :: iteruv_first=.true.
  !
!  ##ifdef PETSC
!#include "finclude/petsc.h"
!#else
!  include 'mpif.h'
!#endif

  ! communication part
  type com_struct
     integer    :: rPEnum
     integer, dimension(:), pointer :: rPE
     integer, dimension(:), pointer :: rptr
     integer, dimension(:), pointer :: rlist
     integer    :: sPEnum
     integer, dimension(:), pointer :: sPE
     integer, dimension(:), pointer :: sptr
     integer, dimension(:), pointer :: slist
  end type com_struct
  type(com_struct)                          :: com_nod2D, com_nod3D
  type com_array
     real(kind=8), dimension(:), pointer :: array
  end  type com_array
  type(com_array), allocatable              :: s_buff_2d(:), r_buff_2d(:)
  type(com_array), allocatable              :: s_buff_3d(:), r_buff_3d(:)

  ! general MPI part
  integer                                   :: MPIERR
  integer                                   :: npes
  integer                                   :: mype
  integer, allocatable, dimension(:)        :: part2D, part3D       

  ! Mesh partition
  integer                                   :: myDim_nod2D, eDim_nod2D, ToDim_nod2D
  integer, allocatable, dimension(:)        :: myList_nod2D
  integer                                   :: myDim_nod3D, eDim_nod3D, ToDim_nod3D
  integer, allocatable, dimension(:)        :: myList_nod3D
  integer                                   :: myDim_elem2D
  integer, allocatable, dimension(:)        :: myList_elem2D
  integer                                   :: myDim_elem3D  
  integer, allocatable, dimension(:)        :: myList_elem3D

  type sparse_matrix 
     integer :: nza
     integer :: dim
     real(kind=8),    pointer,   dimension(:) :: values
     integer(KIND=4), pointer,   dimension(:) :: colind
     integer(KIND=4), pointer,   dimension(:) :: rowptr
  end type sparse_matrix
  !
  type addresstype
     integer                                :: nmb
     integer(KIND=4), dimension(:), pointer :: addresses
  end type addresstype
  !
!
!----------------------------------------------------------------------
!
  ! 
  type(sparse_matrix)                             :: uvstiff, sshstiff
#ifndef use_non_hydrostatic  
  real(kind=8), allocatable, dimension(:,:,:)     :: wpot_matrix
#else
  type(sparse_matrix)                             :: nhpstiff
#endif
  type(sparse_matrix)                             :: tsstiff
  !
  real(kind=8), allocatable, dimension(:)         :: uv_lump 
  real(kind=8), allocatable, dimension(:)         :: uv_lump_prev 
  real(kind=8), allocatable, dimension(:)         :: ts_lump
!
!--------------------------------------------------------------------
!
  !
  integer, allocatable, dimension(:)           :: mapping
  integer, allocatable, dimension(:)           :: col_pos
  !
  integer                                      :: nod2D        
  real(kind=8), allocatable, dimension(:,:)    :: coord_nod2D  
  integer, allocatable, dimension(:)           :: index_nod2D  
  integer                                      :: nod3D        
  real(kind=8), allocatable, dimension(:,:)    :: coord_nod3D  
  integer(KIND=4), allocatable, dimension(:)   :: index_nod3D  
  real(kind=8), allocatable, dimension(:)      :: cos_elem2D 
  real(kind=8), allocatable, dimension(:)      :: geolat
  !
  type(addresstype), allocatable, dimension(:) :: nod_in_elem3D     
  type(addresstype), allocatable, dimension(:) :: nod_in_elem2D     
  type(addresstype), allocatable, dimension(:) :: nod_in_opbnd_tri
  type(addresstype), allocatable, dimension(:) :: nghbr_nod3D
  type(addresstype), allocatable, dimension(:) :: nghbr_nod2D
  type(addresstype), allocatable, dimension(:) :: nghbr_elem3D
  type(addresstype), allocatable, dimension(:) :: nghbr_elem2D
  !
  integer                                      :: max_num_layers
  integer, allocatable, dimension(:)           :: num_layers_below_nod2D
  integer, allocatable, dimension(:)           :: layerdepth
  integer(KIND=4), allocatable, dimension(:,:) :: nod3D_below_nod2D  
  integer(KIND=4), allocatable, dimension(:)   :: nod2D_corresp_to_nod3D 
  !
  integer(KIND=4), allocatable, dimension(:)   :: bt_nds
  !
  integer                                      :: nmbr_opbnd_n2D, nmbr_opbnd_t2D
  integer                                      :: nmbr_opbnd_n3D, nmbr_opbnd_tri
  integer                                      :: nmbr_opbnd_edg
  integer, allocatable, dimension(:)           :: opbnd_n2D, opbnd_n3D
  integer, allocatable, dimension(:)           :: mapping_opbnd_n2d
  integer, allocatable, dimension(:,:)         :: opbnd_tri, opbnd_edg
  real(kind=8), allocatable, dimension(:,:)    :: opbnd_nv, opbnd_edg_nv 

  ! for cases using cavity
  integer, allocatable, dimension(:)           :: cavity_flag_nod2d

  ! for sigma or hybrid grids
  integer, allocatable, dimension(:)           :: grid_type_elem2d
  integer, allocatable, dimension(:,:,:)       :: dens_interp_nodes
  integer, allocatable, dimension(:)           :: elem3d_layer
  real(kind=8), allocatable, dimension(:,:,:)  :: grid_slope
!
!----------------------------------------------------------------------------
!
  integer                                      :: elem2D
  integer(KIND=4), allocatable, dimension(:,:) :: elem2D_nodes 
  integer(KIND=4), allocatable, dimension(:,:) :: elem2D_nghbrs 
  integer                                      :: elem3D
  integer(KIND=4), allocatable, dimension(:,:) :: elem3D_nodes 
  integer(KIND=4), allocatable, dimension(:,:) :: elem3D_nghbrs  
  integer(KIND=4), allocatable, dimension(:)   :: elem2D_corresp_to_elem3D 
  !
  real(kind=8)                                 :: Vol2D, Vol3D   
  real(kind=8)                                 :: sProd_2Di, sProd_3Di 
  real(kind=8), allocatable, dimension(:,:)    :: sProd_2Dij, sProd_3Dij 
  real(kind=8), allocatable, dimension(:,:)    :: derivative_stdbafu_x_2D 
  real(kind=8), allocatable, dimension(:,:)    :: derivative_stdbafu_x_3D  
  real(kind=8), allocatable, dimension(:,:)    :: bafux_3d, bafuy_3d, bafuz_3d
  real(kind=8), allocatable, dimension(:)      :: voltetra
  real(kind=8), allocatable, dimension(:,:)    :: bafux_2d, bafuy_2d
  real(kind=8), allocatable, dimension(:)      :: voltriangle
  !
#ifdef use_fullfreesurf
  integer(kind=4), allocatable, dimension(:)   :: map_elem
  real(kind=8), allocatable, dimension(:)      :: voltetra_new
  real(kind=8), allocatable, dimension(:,:)    :: bafux_3d_new, bafuy_3d_new, bafuz_3d_new
#endif
  !
  real(kind=8)                                 :: ocean_area
  real(kind=8)                                 :: blkflx_area, aegflx_area
  real(kind=8),allocatable, dimension(:)       :: cluster_area_2D
  
!
  integer                  :: vert_nodes(120)
  character(4), allocatable, dimension(:)         :: prog_tracer_name

  real(kind=8), allocatable, dimension(:)         :: coriolis_param_elem2D
  real(kind=8), allocatable, dimension(:)         :: coriolis_param_nod2D  
  real(kind=8), allocatable, dimension(:)         :: stress_x
  real(kind=8), allocatable, dimension(:)         :: stress_y
  real(kind=8), allocatable, dimension(:,:)       :: stress_x_t
  real(kind=8), allocatable, dimension(:,:)       :: stress_y_t
  real(kind=8), allocatable, dimension(:)         :: Tsurf, Ssurf
  real(kind=8), allocatable, dimension(:)        :: heat_flux
  real(kind=8), allocatable, dimension(:,:)       :: heat_flux_t  
  real(kind=8), allocatable, dimension(:)         :: water_flux
  real(kind=8), allocatable, dimension(:,:)       :: water_flux_t         

  real(kind=8), allocatable, target, dimension(:) :: uv_rhs, uf, uf0, duf 
  real(kind=8), allocatable, target, dimension(:) :: ssh, ssh0, dssh, ssh_rhs
#ifndef use_non_hydrostatic
  real(kind=8), allocatable, target, dimension(:) :: wrhs, w  
#else
  real(kind=8), allocatable, target, dimension(:) :: nhp_rhs, nhp, nhp0 
#endif
  real(kind=8), allocatable, target, dimension(:,:) :: tracer_rhs, tracer, dtracer
  
  real(kind=8), allocatable, target, dimension(:,:) :: ts_sfc_force
  real(kind=8), allocatable, target, dimension(:,:) :: uv_sfc_force, uv_bott_force

  real(kind=8), allocatable, dimension(:)         :: virtual_salt, relax_salt
  real(kind=8), allocatable, dimension(:)         :: real_salt_flux

  real(kind=8), allocatable, target, dimension(:) :: ucori, vcori, ucori_back, vcori_back
 
  real(kind=8), allocatable, dimension(:)         :: density_ref, density_insitu  
  real(kind=8), allocatable, dimension(:)         :: hpressure
  real(kind=8), allocatable, dimension(:,:,:)     :: PGF

  ! buffer zone restoring 
  real(kind=8), allocatable, target, dimension(:,:) :: tracer0
  real(kind=8), allocatable, dimension(:)         :: tracer_restore_coeff

  ! restoring open boundary
  real(kind=8), allocatable, dimension(:)         :: opbnd_ssh_rhs

  ! tidal open boundary
  real(kind=8), allocatable, dimension(:)         :: opbnd_u_tide, opbnd_v_tide
  real(kind=8), allocatable, dimension(:)         :: opbnd_z_tide, opbnd_z0_tide
  real(kind=8), allocatable, dimension(:)          :: tide_period_coeff
  real(kind=8), allocatable, dimension(:,:)      :: tide_u_amp, tide_u_pha
  real(kind=8), allocatable, dimension(:,:)      :: tide_v_amp, tide_v_pha
  real(kind=8), allocatable, dimension(:,:)      :: tide_z_amp, tide_z_pha
  real(kind=8), allocatable, dimension(:)         :: opbnd_dep

  ! vertical mixing 
  real, allocatable, dimension(:)       :: Av
  real, allocatable, dimension(:,:)  :: Kv 

  !FCT advection scheme
  real(kind=8), allocatable, dimension(:,:)       :: tral   
  real(kind=8), allocatable, dimension(:,:)       :: trafluxes 
  real(kind=8), allocatable, dimension(:)         :: pplus, pminus

  !Redi/GM flag
  real(kind=8), allocatable, dimension(:,:,:)        :: neutral_slope 
  real(kind=8), allocatable, dimension(:,:)        :: neutral_slope_elem

  contains

subroutine read_namelist

  ! Routine reads namelist files to overwrite default parameters.

  character(len=100)   :: nmlfile
  namelist /clockinit/ timenew, daynew, yearnew

  nmlfile ='namelist.config'    ! name of general configuration namelist file

  fileID=get_unit()
  open (fileID,file=nmlfile)
  read (fileID,NML=modelname)
  read (fileID,NML=timestep)
  read (fileID,NML=timeseries)
  read (fileID,NML=clockinit) 
  read (fileID,NML=paths)
  read (fileID,NML=initialization)  
  read (fileID,NML=inout)
  read (fileID,NML=mesh_def)
  read (fileID,NML=geometry)
!  read (fileID,NML=calendar)
  close (fileID)
  if ( do_output() ) print*, "namelist is read"
end subroutine read_namelist
subroutine read_elem
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !!!!!! READ elem3d.out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!   mype=0
!   fileID=mype+10
!   open(fileID, file=(trim(MeshPath)//'elem3d.out'))
!     read(fileID,*) myDim_elem3D   
!     allocate(elem3D_nodes(4, myDim_elem3D))
!       read(fileID,*) elem3d_nodes
!   close(fileID)
! 
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !!!!!! READ elem2d.out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!   fileID=mype+10
!   open(fileID, file=(trim(MeshPath)//'elem2d.out'))
!     read(fileID,*) myDim_elem2D   
!     allocate(elem2D_nodes(3, myDim_elem2D))
!       read(fileID,*) elem2d_nodes
!   close(fileID)
!   if ( do_output() ) print*, 'Element tables are read: '

end subroutine read_elem
subroutine read_node
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! READ nod3d.out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real*8     :: x,y,z

  mype=0
  fileID=get_unit()
  open(fileID, file=(trim(MeshPath)//'nod3d.out'))
    read(fileID,*) myDim_nod3D
    allocate(coord_nod3D(3,myDim_nod3D))
    allocate(index_nod3D(myDim_nod3D))
    allocate(myList_nod3D(myDim_nod3D))
      do n=1,myDim_nod3D
        read(fileID,*) m, x, y, z, ind
        coord_nod3D(1,n)=x
        coord_nod3D(2,n)=y
        coord_nod3D(3,n)=z
        index_nod3D(n)=ind
        myList_nod3D(n)=m
     end do
  close(fileID)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! READ nod2d.out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fileID=get_unit()
  open(fileID, file=(trim(MeshPath)//'nod2d.out'))
    read(fileID,*) myDim_nod2D
    allocate(coord_nod2D(2,myDim_nod2D))
    allocate(index_nod2D(myDim_nod2D))
    allocate(myList_nod2D(myDim_nod3D))
      do n=1,myDim_nod2D
        read(fileID,*) m, x, y, ind
        coord_nod2D(1,n)=x
        coord_nod2D(2,n)=y
        index_nod2D(n)=ind
        myList_nod2D(n)=m
      end do
  close(fileID)
  if ( do_output() ) print*, 'Node tables are read: '

end subroutine read_node

subroutine read_aux3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! READ aux3d.out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fileID=get_unit()
  open(fileID, file=trim(meshpath)//'aux3d.out')
  read(fileID, *) max_num_layers

!!!!===========================================================!!!!
!!!!=========nod3D_below_nod2D=================================!!!!
!!!!===========================================================!!!!

  allocate(nod3D_below_nod2D(max_num_layers,myDim_nod2D))       
  do n=1, myDim_nod2D
     read(fileID, *) vert_nodes(1:max_num_layers)
        nod3D_below_nod2D(:,n)=vert_nodes(1:max_num_layers)
  end do

!!!!===========================================================!!!!
!!!!=========nod2D_corresp_to_nod3D============================!!!!
!!!!===========================================================!!!!

  allocate(nod2D_corresp_to_nod3D(myDim_nod3D)) 

  do n=1, myDim_nod3D
     read(fileID, *) m
        nod2D_corresp_to_nod3D(n)=m
  end do

  do n=1,myDim_nod3D
     m=nod2D_corresp_to_nod3D(n)
     nod2D_corresp_to_nod3D(n)=m
  end do

!!!!===========================================================!!!!
!!!!========elem2D_corresp_to_elem3D===========================!!!!
!!!!===========================================================!!!!

  allocate(elem2D_corresp_to_elem3D(myDim_elem3D)) 
  do n=1, myDim_elem3D
     read(fileID,*) m
        elem2D_corresp_to_elem3D(n)=m
  end do

  close(fileID)
  if ( do_output() ) print*, 'Aux3D tables are read: '

end subroutine read_aux3

subroutine read_depth
  real*8 :: depth

  fileID=get_unit()
  allocate(layerdepth(max_num_layers)) 
  open(fileID, file=(trim(MeshPath)//'m3d.ini'))
  do n=1,max_num_layers-1
     read(fileID,*) depth
       layerdepth(n)= depth
  end do 
  print*, layerdepth
  close(fileID)
end subroutine read_depth

subroutine oce_input
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j, dimid_rec, nrec
  integer                   :: ssh_varid, tra_varid(2)
  integer                   :: u_varid, v_varid, w_varid, wpot_varid
  integer                   :: istart(2), icount(2), n3
  character(100)            :: filename
  real(kind=8), allocatable :: aux2(:), aux3(:) 

  allocate(aux2(myDim_nod2D), aux3(myDim_nod3D)) 
  !n3=ToDim_nod3D           
  n3=myDim_nod3D           
  nrec=day2ext

  ! open files
  filename=trim(ResultPath)//runid//'.'//runyear//'.oce.nc'
  if ( do_output() ) print*,'Input filename: ',filename
  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) call handle_err(status)
  ! inquire variable id
  status=nf_inq_varid(ncid, 'ssh', ssh_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'u', u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'v', v_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#ifdef use_non_hydrostatic
  status=nf_inq_varid(ncid, 'w', w_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#else
  status=nf_inq_varid(ncid, 'wpot', wpot_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#endif
  status=nf_inq_varid(ncid, 'temp', tra_varid(1))
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'salt', tra_varid(2))
  if (status .ne. nf_noerr) call handle_err(status)

  ! read variables

  ! which record to read
!  if(restartflag=='last') then
!     status = nf_inq_dimid(ncid, 'T', dimid_rec)
!     if(status .ne. nf_noerr) call handle_err(status)
!     status = nf_inq_dimlen(ncid, dimid_rec, nrec)
!     if(status .ne. nf_noerr) call handle_err(status)
!  else
!     read(restartflag,'(i4)') nrec
!  end if

  ! 2d fields
  istart=(/1,nrec/)
  icount=(/myDim_nod2D, 1/)
  status=nf_get_vara_double(ncid, ssh_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  ssh=aux2(1:myDim_nod2D)         

  ! 3d fields
  istart=(/1,nrec/)
  icount=(/myDim_nod3D, 1/)

  do j=1,2
     status=nf_get_vara_double(ncid, tra_varid(j), istart, icount, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
     tracer(:,j)=aux3(1:myDim_nod3D)  
  end do

  status=nf_get_vara_double(ncid, u_varid, istart, icount, aux3)
  if (status .ne. nf_noerr) call handle_err(status)
  uf(1:n3)=aux3(1:myDim_nod3D)     

  status=nf_get_vara_double(ncid, v_varid, istart, icount, aux3) 
  if (status .ne. nf_noerr) call handle_err(status)
  uf(1+n3:2*n3)=aux3(1:myDim_nod3D)  

#ifdef use_non_hydrostatic
  status=nf_get_vara_double(ncid, w_varid, istart, icount, aux3) 
  uf(1+2*n3:3*n3)=aux3(myDim_nod3D)
#else
  status=nf_get_vara_double(ncid, wpot_varid, istart, icount, aux3)
  w=aux3(myDim_nod3D)             
#endif
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! the next record to be saved
  save_count=nrec+1

  deallocate(aux3, aux2)   
!  print*, 'Inputs are allocated: '

end subroutine oce_input

subroutine handle_err(errcode) 
implicit none 
       
#include "netcdf.inc"  

integer errcode 
         
write(*,*) 'Error: ', nf_strerror(errcode) 
!  call par_ex 
stop
end subroutine handle_err
end module fesom_modules

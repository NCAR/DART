
! FESOM 2 (Finite-volumE Sea ice-Ocean Model)
! multi-resolution ocean general circulation model
! FESOM/fesom2 is licensed under the GNU General Public License v2.0
! Copyright (C) 2018  FESOM team
!
! This source file was taken from  the FESOM V1.4 modules

module o_param
  implicit none
  save
  !
  integer :: fileID
  ! *** Fixed parameters ***
  real(kind=8), parameter  	:: pi=3.141592653589793, rad=pi/180.0
  real(kind=8), parameter  	:: omega=2.0*pi/(24.0*60.0*60.0)
  real(kind=8), parameter  	:: g=9.81                       ![m/s^2]
  real(kind=8), parameter  	:: r_earth=6.3675e6             ![m]
  real(kind=8), parameter  	:: rho0=1028.5                  ![kg/m^3]marmara sea deep water density
  real(kind=8), parameter 	:: rho0r=1.0/rho0
  real(kind=8), parameter  	:: vcpw=4.2e6                   ![J/m^3/K]volum. heat cap. of water
  real(kind=8), parameter	:: small=1.0e-8                 !small value

  ! *** mixing and friction setting ***
  real             	        :: Ah0=6000.            	!Lapl. hori. visc, [m^2/s]
  real             	        :: Ahb0=2.7e13             	!Bihar. hori. visc,[m^4/s] <8.0e12/1degree
  real             	        :: Kh0=600.                     !lateral diff
  !
  logical                  	:: biharmonic_visc=.false.      !bihar. true will turn off Lapl.
  logical                  	:: smagorinsky_visc=.false.
  !
  logical                       :: increase_equ_zonal_visc=.true.!increase zonal viscosity at equator
  real                          :: fac_visc_increase=3.0         !increase factor
  !
  logical                  	:: scale_mixing_h=.true.   	!scale hor. diff/visc. coef. by resol.
  integer			:: scale_mixing_type=2		!1:scale by \delt x^2; 2:scale by \delta x
  real(kind=8)             	:: scalevol=5.e9          	!associated with 100km resolution
  !
  logical	  	   	:: Redi_GM=.true.         	!flag for Redi/GM scheme
  logical                       :: ODM95=.true.                 !taper fcn ODM95
  logical                       :: LDD97=.true.                 !taper fcn LDD97
  real(kind=8)		   	:: ratio_K_GM=1.0          	!ratio of K_GM to Kh
  real(kind=8)                  :: Slope_critical=4.0e-3        !upper slope limit for applying Redi/GM
  integer                       :: nslope_version=1             !1-neutral slope over prism,2-over tetrahedra
  !
  real             	        :: Av0=1.e-4      	        !background (or internal wave) vert. mixing
  real              	        :: Kv0=1.e-5                    !m^2/s
  real                          :: visc_conv_limit=0.1          !visc due to convective instability
  real                          :: diff_conv_limit=0.1          !diff due to convective instability
  !
  character(5)                 	:: mix_scheme='KPP'		!'KPP','PP', 'MY2p5', 'no'
  !
  real                          :: visc_sh_limit=5.0e-3         !for kpp,max visc due to shear instability
  real                          :: diff_sh_limit=5.0e-3         !for kpp,max diff due to shear instability
  logical                       :: double_diffusion=.true.      !for KPP,dd switch
  logical                       :: smooth_blmc=.true.           !for KPP,hori. smooth of b.l. mixing coeff.
  !
  real             	        :: PP_max_mix_coeff=5.0e-3     	!for PP, max Kv/Av
  real                          :: wndmix=1.0e-3                !for PP, to simulate missing high frequency wind
  logical                       :: allow_convect_global=.true.  !for PP, convection for global or only NH
  !
  logical                  	:: add_TB04_to_PP=.false.   	!for PP, TB04 switch
  real(kind=8)             	:: modiff=0.01                  !for PP, vert. mix. coeff. for TB04
  !
  logical			:: tidal_mixing=.false.		!switch for tidal mixing
  logical			:: use_drag_dissipation=.true. 	!barotropic
  logical			:: use_wave_dissipation=.false.	!baroclinic
  logical			:: read_tide_speed=.true.	!read tide speed or use default
  real(kind=8)                  :: default_tide_speed=0.01      !(m/s)
  real(kind=8)                  :: max_drag_diffusivity=5.e-3   !m2/s
  real(kind=8)                  :: max_wave_diffusivity=5.e-3   !m2/s
  character(2) 	                :: Tmix_tidalconstituent='M2'   !which tidal constituent
  character(15)	                :: Tmix_tidalmodelname='tpxo71' !source model name
  !
  real(kind=8)             	:: C_d=0.0025               !Bottom fri. coeff.
  real(kind=8)             	:: aegflx_lat=39.0          !aegean sea flux latitude
  real(kind=8)             	:: blkflx_lat=42.5          !black sea flux latitude
  real(kind=8)             	:: marmin_lat=40.0          !marmarasea min latitude
  real(kind=8)             	:: marmax_lat=41.2          !marmarasea max latitude
  real(kind=8)             	:: marmin_lon=27.0          !marmarasea min longitude
  real(kind=8)             	:: marmax_lon=30.0          !marmarasea max longitude
  real(kind=8)             	:: km3yr2m3sec=1e+09/3.1536e+07

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
  logical                  	:: ts_surfbd=.true.
  !
  logical                       :: ref_sss_local=.false.	!virtual salt flux using local SSS or ref_sss
  real(kind=8)             	:: ref_sss=34.7			!ref. value for cal. virtual salt flux
  !
  real(kind=8)             	:: restore_s_surf=10./(180.*86400.)	! m/s
  real(kind=8)			:: restore_t_surf=0.0
  !
  logical                       :: balance_salt_water=.true.    !balance virtual salt or water flux or not
  !
  logical                       :: buffer_zone=.false.
  real(kind=8)             	:: restore_ts_buff= 1./(86400.*5.)     	! timescale for buffer zone [1/s]

  namelist /boundary/ ts_surfbd, ref_sss_local, ref_sss, restore_s_surf, &
       restore_t_surf, balance_salt_water, buffer_zone, restore_ts_buff

  ! *** numerical schemes
  real(KIND=8)             	:: gamma_stab=0.99          	!stabilization for ssh
  real(KIND=8)             	:: gamma_stab_nh=0.5       	!stabilization for nh-pressure
  real(kind=8)             	:: gamma_fct=0.4           	!param for tracer FCT scheme
  real(kind=8)             	:: alpha_AB=1.55            	!when Adams-Bashforth Coriolis
  real(kind=8)             	:: alpha_trapez=0.55        	!when semi-implicit Coriolis
  real(kind=8)			:: theta_ssh=0.5		!semi-implicit ssh when semiimpl sheme
  real(kind=8)			:: theta_vel=0.5		!semi-implicit baro. vel
  !
  logical                       :: use_vertvisc_impl=.true.     !if implicit vertical viscosity,keep true
  logical                       :: use_vertdiff_impl=.true.     !if implicit vertical diff., keep ture
  logical                       :: use_cori_semi=.false.        !if semiimplicit coriolis force
  !
  logical                  	:: lump_uv_matrix=.true.   	!for mass vel. matrix case, keep true!
  logical                  	:: lump_ts_matrix=.true.   	!for mass T/S matrix case, keep true!
  integer                  	:: num_iter_solve=3        	!iteration # for mass matrix case

  namelist /oce_scheme/ gamma_stab, gamma_stab_nh, gamma_fct, alpha_AB, alpha_trapez, &
       theta_ssh, theta_vel, use_cori_semi

  ! *** density and pressure force ***
  logical                  	:: density_linear=.false.
  logical                 	:: use_ref_density=.true.

  namelist /denspress/ density_linear, use_ref_density

  ! *** parameters for nonlinear free surface cases ***
  real(kind=8)             	:: max_ice_loading=5.0		!m, maximal pressure from ice felt by the ocean

  namelist /param_freesurf/ max_ice_loading

  ! *** tide configuration ***
  integer     	   	        :: nmbr_tidal_cons=4
  character(20)	                :: tidal_constituent='M2S2K1O1' !M2 S2 N2 K2 K1 O1 P1 Q1
  character(15)			:: tidemodelname='tpxo71'
  character(10)                 :: tide_opbnd_type='Flather' 	!ssh, vel, Flather
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

end module o_param
  !
  !no user specification beyond this line!
  !---------------------------------------------------------------------------------------
  !
module o_array
  implicit none
  save

  character(4), allocatable, dimension(:)         :: prog_tracer_name

  real(kind=8), allocatable, dimension(:)         :: coriolis_param_elem2D
  real(kind=8), allocatable, dimension(:)         :: coriolis_param_nod2D
  real(kind=8), allocatable, dimension(:)    	  :: stress_x
  real(kind=8), allocatable, dimension(:)    	  :: stress_y
  real(kind=8), allocatable, dimension(:,:)  	  :: stress_x_t
  real(kind=8), allocatable, dimension(:,:)  	  :: stress_y_t
  real(kind=8), allocatable, dimension(:)    	  :: Tsurf, Ssurf
  real(kind=8), allocatable, dimension(:)  	  :: heat_flux
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
  real(kind=8), allocatable, dimension(:)    	  :: tide_period_coeff
  real(kind=8), allocatable, dimension(:,:)	  :: tide_u_amp, tide_u_pha
  real(kind=8), allocatable, dimension(:,:)	  :: tide_v_amp, tide_v_pha
  real(kind=8), allocatable, dimension(:,:)	  :: tide_z_amp, tide_z_pha
  real(kind=8), allocatable, dimension(:)         :: opbnd_dep

  ! vertical mixing
  real, allocatable, dimension(:) 	  :: Av
  real, allocatable, dimension(:,:)  :: Kv

  !FCT advection scheme
  real(kind=8), allocatable, dimension(:,:)       :: tral
  real(kind=8), allocatable, dimension(:,:)       :: trafluxes
  real(kind=8), allocatable, dimension(:)         :: pplus, pminus

  !Redi/GM flag
  real(kind=8), allocatable, dimension(:,:,:)  	  :: neutral_slope
  real(kind=8), allocatable, dimension(:,:)  	  :: neutral_slope_elem

end module o_array
!
!----------------------------------------------------------------------
!
module o_solver
  implicit none
  save
  !solver index
  integer                  	:: solve_u=1
  integer       	 	:: solve_v=2
  integer                  	:: solve_w=6
  integer                  	:: solve_ssh=7
  integer                  	:: solve_nhp=8
  integer                  	:: solve_tra=10
  logical                  	:: iter_first=.true.
  logical                  	:: iteruv_first=.true.
end module o_solver
!
!----------------------------------------------------------------------
!
module g_parfe
  implicit none
  save
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
  type(com_struct)                      	:: com_nod2D, com_nod3D
  type com_array
     real(kind=8), dimension(:), pointer :: array
  end  type com_array
  type(com_array), allocatable          	:: s_buff_2d(:), r_buff_2d(:)
  type(com_array), allocatable          	:: s_buff_3d(:), r_buff_3d(:)

  ! general MPI part
  integer                               	:: MPIERR
  integer           				:: npes
  integer        				:: mype
  integer, allocatable, dimension(:)  		:: part2D, part3D

  ! Mesh partition
  integer                             		:: myDim_nod2D, eDim_nod2D, ToDim_nod2D
  integer, allocatable, dimension(:)  		:: myList_nod2D
  integer                             		:: myDim_nod3D, eDim_nod3D, ToDim_nod3D
  integer, allocatable, dimension(:)  		:: myList_nod3D
  integer                             		:: myDim_elem2D
  integer, allocatable, dimension(:)  		:: myList_elem2D
  integer                             		:: myDim_elem3D
  integer, allocatable, dimension(:)  		:: myList_elem3D
end module g_PARFE
!
!----------------------------------------------------------------------
!
module o_DATA_TYPES
  implicit none
  save
  !
  type sparse_matrix
     integer :: nza
     integer :: dim
     real(kind=8), pointer, dimension(:)      :: values
     integer(KIND=4), pointer,   dimension(:) :: colind
     integer(KIND=4), pointer,   dimension(:) :: rowptr
  end type sparse_matrix
  !
  type addresstype
     integer                                :: nmb
     integer(KIND=4), dimension(:), pointer :: addresses
  end type addresstype
  !
end module o_DATA_TYPES
!
!----------------------------------------------------------------------
!
module o_MATRICES
  !
  use o_DATA_TYPES
  implicit none
  save
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
end module o_MATRICES
!
!--------------------------------------------------------------------
!
module o_mesh
  !
  use o_DATA_TYPES
  implicit none
  save
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
  real(kind=8), allocatable, dimension(:)           :: layerdepth
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
end module o_mesh
!
!----------------------------------------------------------------------------
!
module o_elements
  implicit none
  save
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
  real(kind=8)                                 :: ocean_area, domain_area
  real(kind=8)                                 :: blkflx_area, aegflx_area
  real(kind=8)                                 :: marm_area, bosp_area, dard_area
  real(kind=8),allocatable, dimension(:)       :: cluster_area_2D

end module o_elements
!
module o_read
  character(100),parameter :: MeshPath='/users/home/ans051/FEOM_PREPROC/mesh-T2G1.5L110b/'
  integer                  :: vert_nodes(120)
end module o_read


! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8, missing_r8
use time_manager_mod, only : time_type, set_time, get_time, &
                             increment_time, print_time, set_date, &
                             set_calendar_type, GREGORIAN, &
                             operator(==), operator(<=), &
                             operator(-), operator(+)
use     location_mod, only : location_type, get_dist, set_location, &
                             get_location, query_location, &
                             LocationDims, LocationName, LocationLName, &
                             vert_is_surface, vert_is_pressure, &
                             vert_is_level, vert_is_height, &
                             get_close_maxdist_init, get_close_obs_init, get_close_obs

use    utilities_mod, only : file_exist, open_file, close_file, &
                             find_namelist_in_file, check_namelist_read, &
                             register_module, error_handler, E_ERR, E_MSG, logfileunit
use         sort_mod, only : sort   
use   random_seq_mod, only : random_seq_type, random_gaussian, &
                             init_random_seq, random_uniform
use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, &
                             KIND_V_WIND_COMPONENT, &
                             KIND_SURFACE_PRESSURE, &
                             KIND_TEMPERATURE, &
                             KIND_SPECIFIC_HUMIDITY 
use        map_utils, only : proj_info, map_init, map_set, latlon_to_ij, &
                             PROJ_LATLON, PROJ_MERC, PROJ_LC, PROJ_PS, &
                             gridwind_to_truewind

!---- WRF modules
use module_model_constants
use module_initialize
use module_namelist
use module_wrf

!use            netcdf

implicit none
private

public :: get_model_size, &
          adv_1step, &
          get_state_meta_data, &
          model_interpolate, &
          get_model_time_step, &
          end_model, &
          static_init_model, &
          init_time, &
          init_conditions, &
          model_get_close_states, &
          nc_write_model_atts, &
          nc_write_model_vars, &
          nc_read_model_vars, &
          pert_model_state, &
          get_pblh, &
          pert_params_time, &
          adjust_param_spread, &
          real_obs_period, &
          start_real_obs, &
          synchronize_mavail, &
          ok_to_nudge, &
          get_close_maxdist_init, get_close_obs_init, get_close_obs, ens_mean_for_model


! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! Basic model parameters controlled by namelist

integer                :: num_est_params  = 0
integer, dimension(10) :: est_param_types = 600
real   , dimension(10) :: pert_init_sd    = 0.0
real   , dimension(10) :: pert_param_sd   = 0.0
integer, dimension(10) :: pert_init_beta_1 = -1.0
integer, dimension(10) :: pert_init_beta_2 = -1.0
real   , dimension(10) :: pert_param_min  = 0.00001
real   , dimension(10) :: pert_param_max  = 0.99999
logical                :: maintain_initial_spread = .false.
character(len=4), dimension(10) :: dist_shape = 'logn'
integer                :: real_obs_period = 1800
integer                :: start_real_obs  = 1800

namelist /model_nml/ num_est_params, est_param_types, pert_param_sd, &
         pert_init_sd, pert_init_beta_1, pert_init_beta_2, &
         maintain_initial_spread, dist_shape, pert_param_min, &
         pert_param_max, real_obs_period, start_real_obs


! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step, initialization_time

! Private definition of model variable types

integer, parameter :: TYPE_U   = 1,   TYPE_V   = 2,  TYPE_W  = 3,  &
                      TYPE_GZ  = 4,   TYPE_T   = 5,  TYPE_MU = 6,  &
                      TYPE_QV  = 7,   TYPE_QC  = 8,  TYPE_QR = 9,  &
                      TYPE_QI  = 10,  TYPE_QS  = 11, TYPE_QG = 12, &
                      TYPE_U10 = 13,  TYPE_V10 = 14, TYPE_T2 = 15, &
                      TYPE_Q2  = 16,  TYPE_PSFC= 17, TYPE_TSLB = 18, &
                      TYPE_TSK = 19,  TYPE_TH = 20,  TYPE_TKE = 21 , &
                      TYPE_P   = 22

! additional "far away" types that will not be part of the assimilation
! GROUP THESE BY SIZE OF ARRAYS!  this will make life easier below
integer, parameter :: TYPE_GSWF  = 100,    TYPE_GLWF = 101, & !1D       
                      TYPE_PRECIPF = 102

integer, parameter :: TYPE_UF = 103,       TYPE_VF = 104, &   !2D...
                      TYPE_PF = 105,       TYPE_P8WF = 106

integer, parameter :: TYPE_TH_UPSTREAM_X = 107, TYPE_TH_UPSTREAM_Y = 108, &
                      TYPE_QV_UPSTREAM_X = 109, TYPE_QV_UPSTREAM_Y = 110, &
                      TYPE_TAU_U = 111,        TYPE_TAU_V = 112

! far away types that we may want to add to the state vector at some
! point
integer, parameter :: TYPE_SMOIS = 200,      TYPE_KEEPFR3DFLAG = 201, & ! nsoil
                      TYPE_SMFR3D = 202,     TYPE_SH2O = 203

integer, parameter :: TYPE_UZ0 = 300,        TYPE_VZ0 = 301, &   !scalars
                      TYPE_THZ0 = 302,       TYPE_QZ0 = 303, &
                      TYPE_QVG  = 304,       TYPE_QSG = 305, &
                      TYPE_QCG  = 306,       TYPE_QSFC = 307, &
                      TYPE_AKMS = 308,       TYPE_AKHS = 309, &
                      TYPE_HOL = 310,        TYPE_MOL = 311, &
                      TYPE_GRDFLX = 312,     TYPE_HFX = 313, &
                      TYPE_QFX = 314,        TYPE_PSHLTR = 315, &
                      TYPE_QSHLTR = 316,     TYPE_FLHC = 317, &
                      TYPE_FLQC = 318,       TYPE_Q10 = 319, &
                      TYPE_UDRUNOFF = 320,   TYPE_ACSNOM = 321, &
                      TYPE_ACSNOW = 322,     TYPE_UST = 323, &
                      TYPE_PBLH = 324,       TYPE_TMN = 325, &
                      TYPE_CS = 326,         TYPE_ERATE = 327, &
                      TYPE_MAXM = 328,       TYPE_MINM = 329
integer, parameter :: TYPE_RMOL = 330,       TYPE_ZNT = 331, &
                      TYPE_GAMMA = 332

! parameter section - first some strange ones
integer, parameter :: TYPE_VEGFRA = 400

! far away profiles may or may not be part of the external forcing - perhaps
! diagnostic (U_G and V_G are only for output - redundant)
integer, parameter :: TYPE_RHO = 500, TYPE_U_G  = 501, &
                      TYPE_V_G  = 502, TYPE_EL = 503

! parameter section - some out of state (constant) and are far away
! others are in-state and will be assigned a true location
integer, parameter ::  TYPE_EMISS = 600, TYPE_ALBEDO = 601, &
                       TYPE_Z0    = 602, TYPE_THC    = 603, &
                       TYPE_MAVAIL= 604

integer, parameter :: calendar_type = GREGORIAN
integer            :: internal_ensemble_counter = 0
integer            :: number_ensemble_members = 0

! these are parameters that will be stuffed in the the meta data type
! soil variables in state vector to dart
integer, parameter  ::   number_of_soil_variables = 1 
                         ! (TSLB)
! the usual profile variables in state to dart
integer, parameter  ::   number_of_profile_variables =     8 
                         ! (U, V, Z, T, QV, TH, TKE,P)
! screen height vars in state vector to dart
integer, parameter  ::   number_of_scalars =               6      
                         ! (U10, V10, T2, Q2, TSK, PSFC) 
! external forcing
integer, parameter  ::   number_of_1d_f =                3
                         ! (GSWF, GLWF)
integer, parameter  ::   number_of_2d_f =               4 
                         !(UF, VF, PF, P8WF
integer, parameter  ::   number_of_2d_advection =        6
                         ! T_UPSTREAM_X, T_UPSTREAM_Y
                         ! TAU_U, TAU_V
! profile variables NOT in state vector to dart
integer, parameter  ::   number_noassim_profiles =         4
                         ! RHO, U_G, V_G
! soil variables NOT in state vector to dart
integer, parameter  ::   number_noassim_soil_vars =        4 
                         ! (SMOIS, KEEPFR3DFLAG, 
                         !  SMFR3D, SH2O)
! scalars NOT in state vector to dart
integer, parameter  ::   number_noassim_scalars = 33
                         ! (UZ0, VZ0, THZ0, QZ0, 
                         !  QVG, QSG, QCG, QSFC, 
                         ! AKMS, AKHS, HOL, MOL, 
                         ! GRDFLX, HFX, QFX, 
                         ! PSHLTR, QSHLTR, 
                         ! FLQC, FLHC, Q10, UDRUNOFF, 
                         ! ACSNOM, ACSNOW, UST, PBLH, TMN,
                         ! RMOL, ZNT, GAMMA)

! parameters that depend on initialization date
integer, parameter  ::   number_dependent_params = 1
                         ! (VEGFRA)
! parameters that are independent, and MAY be in the state
integer, parameter  ::   number_total_params = 5
                         ! (EMISS, ALBEDO, Z0, THC, MAVAIL)

TYPE state_vector_meta_data
! additional necessary vars
   integer             ::   model_size
   integer             ::   total_number_of_vars
! soil variables in state vector to dart
   integer             ::   number_of_soil_variables   
! the usual profile variables in state to dart
   integer             ::   number_of_profile_variables
! screen height vars in state vector to dart
   integer             ::   number_of_scalars         
! external forcing
   integer             ::   number_of_1d_f 
   integer             ::   number_of_2d_f 
   integer             ::   number_of_2d_advection 
! profile variables NOT in state vector to dart
   integer             ::   number_noassim_profiles
! soil variables NOT in state vector to dart
   integer             ::   number_noassim_soil_vars   
! scalars NOT in state vector to dart
   integer             ::   number_noassim_scalars
! parameters that depend on initialization date
   integer             ::   number_dependent_params
! parameters that are independent, and are not in the state
   integer             ::   number_independent_params
! parameters in the state, will be estimated
   integer             ::   number_state_params

! list of parameter types and properties to estimate
   integer,    dimension(:), pointer :: est_param_types
   real(r8),   dimension(:), pointer :: pert_init_sd
   real(r8),   dimension(:), pointer :: pert_param_sd
   real(r8),   dimension(:), pointer :: pert_param_min
   integer,    dimension(:), pointer :: pert_init_beta_1
   integer,    dimension(:), pointer :: pert_init_beta_2
   real(r8),   dimension(:), pointer :: pert_param_max
   character(len=4),dimension(:), pointer :: dist_shape

   integer, dimension(:), pointer :: var_type
   integer, dimension(:), pointer :: var_size
   integer, dimension(:), pointer :: var_index
end TYPE state_vector_meta_data

TYPE domain_static_data
   real(r8) :: latitude, longitude
end TYPE domain_static_data

type(state_vector_meta_data)    :: wrf_meta

type(proj_info)                 :: my_projection

type(domain_static_data)        :: column

! A flag for allocation, starts as true then goes to false with first
! call to init_wrf, which occurs in init_conditions.

integer                         :: wrf_rnd_seed
logical                         :: allocate_wrf = .true.

! random sequence for the parameter estimates
type (random_seq_type)          :: param_ran_seq
type (random_seq_type)          :: pert_ran_seq

contains

!==================================================================


subroutine static_init_model()
!------------------------------------------------------------------
! Initializes class data for this model. For now, simply outputs the
! identity info, sets the location of the state variables, and initializes
! the time type for the time stepping (is this general enough for time???)
implicit none

integer  :: io, iunit
real(r8) :: x_loc
integer :: dart_index, var_cnt, k, bot, top, unit_nml, vcnt, projcode
LOGICAL :: is_it_there = .FALSE.
REAL :: timeo,timetot
    
real(r8) :: center_i, center_j, truelat1, truelat2
real(r8) :: sw_corner_lat, sw_corner_lon, dx, stdlon
real(r8) :: long_far

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
write(logfileunit, nml=model_nml)
write(     *     , nml=model_nml)

! Begin by reading the namelist input
if(file_exist('wrf1d_namelist.input')) then
   unit_nml = open_file(fname = 'wrf1d_namelist.input', action = 'read')
   call do_namelist_wrf1d(unit_nml,logfileunit)
   close(unit_nml)
endif

! initialize the random sequences
wrf_rnd_seed = rnd_seed_val
call init_random_seq(param_ran_seq)
call init_random_seq(pert_ran_seq)

! need some dimension information
call static_init_wrf(allocate_wrf)
allocate_wrf = .false.

! initialize some timing stuff
time_step = set_time(int(dt), 0)
call set_calendar_type(calendar_type)
call init_time(initialization_time)

! numbers
wrf_meta%number_of_scalars = number_of_scalars
wrf_meta%number_of_profile_variables = number_of_profile_variables
wrf_meta%number_of_soil_variables = number_of_soil_variables
wrf_meta%number_of_1d_f = number_of_1d_f
wrf_meta%number_of_2d_f = number_of_2d_f
wrf_meta%number_of_2d_advection = number_of_2d_advection
wrf_meta%number_noassim_profiles = number_noassim_profiles
wrf_meta%number_noassim_soil_vars = number_noassim_soil_vars
wrf_meta%number_noassim_scalars = number_noassim_scalars
wrf_meta%number_dependent_params = number_dependent_params
wrf_meta%number_state_params = num_est_params !from namelist

! these are not appended or retained in memory because they are constant
wrf_meta%number_independent_params = number_total_params - num_est_params 
                                         

! compute the model size
wrf_meta%model_size = num_soil_layers * wrf_meta%number_of_soil_variables  &
                    + nz              * wrf_meta%number_of_profile_variables &
                    + 1               * wrf_meta%number_of_scalars &
                    + nz*nsplinetimes * wrf_meta%number_of_2d_f &
                    + nz*nsplinetimes_advection * &
                         wrf_meta%number_of_2d_advection &
                    + n1dsplines                                & 
                    + nz              * wrf_meta%number_noassim_profiles &
                    + num_soil_layers * wrf_meta%number_noassim_soil_vars &
                    + 1               * wrf_meta%number_noassim_scalars &
                    + 1               * wrf_meta%number_dependent_params &
                    + 1               * wrf_meta%number_state_params

! state vector locations and wrf grid meta data
wrf_meta%total_number_of_vars = wrf_meta%number_of_scalars &
                              + wrf_meta%number_of_profile_variables &
                              + wrf_meta%number_of_soil_variables &
                              + wrf_meta%number_of_1d_f &
                              + wrf_meta%number_of_2d_f &
                              + wrf_meta%number_of_2d_advection &
                              + wrf_meta%number_noassim_profiles &
                              + wrf_meta%number_noassim_soil_vars &
                              + wrf_meta%number_noassim_scalars &
                              + wrf_meta%number_dependent_params &
                              + wrf_meta%number_state_params

allocate(wrf_meta%est_param_types(wrf_meta%number_state_params))
allocate(wrf_meta%pert_init_sd(wrf_meta%number_state_params))
allocate(wrf_meta%pert_param_sd(wrf_meta%number_state_params))
allocate(wrf_meta%pert_init_beta_1(wrf_meta%number_state_params))
allocate(wrf_meta%pert_init_beta_2(wrf_meta%number_state_params))
allocate(wrf_meta%pert_param_min(wrf_meta%number_state_params))
allocate(wrf_meta%pert_param_max(wrf_meta%number_state_params))
allocate(wrf_meta%dist_shape(wrf_meta%number_state_params))

allocate(wrf_meta%var_type(wrf_meta%total_number_of_vars))
allocate(wrf_meta%var_size(wrf_meta%total_number_of_vars))
allocate(wrf_meta%var_index(wrf_meta%total_number_of_vars))
allocate(state_loc(wrf_meta%model_size))

! fill the parameter estimation types (no error checking)
do i = 1, wrf_meta%number_state_params
  wrf_meta%est_param_types(i) = est_param_types(i)
  wrf_meta%pert_init_sd(i) = pert_init_sd(i)
  wrf_meta%pert_init_beta_1(i) = pert_init_beta_1(i)
  wrf_meta%pert_init_beta_2(i) = pert_init_beta_2(i)
  wrf_meta%pert_param_sd(i) = pert_param_sd(i)
  wrf_meta%pert_param_min(i) = pert_param_min(i)
  wrf_meta%pert_param_max(i) = pert_param_max(i)
  wrf_meta%dist_shape(i) = dist_shape(i)
enddo

! domain info - change lon to [0,360] for DART compliance
if ( lon_ref < 0.0_r8 ) lon_ref = lon_ref + 360.0_r8
column%longitude = lon_ref
column%latitude  = lat_ref
long_far = lon_ref+180.0_r8
if ( long_far > 360.0_r8 ) long_far = long_far - 360.0_r8

! fill locations for state variables
dart_index = 1
var_cnt    = 1

wrf_meta%var_type(var_cnt) = TYPE_U
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location(column%longitude,column%latitude,z_grid(k - bot + 1),3)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_V
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location(column%longitude,column%latitude,z_grid(k - bot + 1),3)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_GZ
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location(column%longitude,column%latitude,z_grid(k - bot + 1),3)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_T
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location(column%longitude,column%latitude,z_grid(k - bot + 1),3)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TH
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location(column%longitude,column%latitude,z_grid(k - bot + 1),3)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_QV
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location(column%longitude,column%latitude,z_grid(k - bot + 1),3)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TKE
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location(column%longitude,column%latitude,z_grid(k - bot + 1),3)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_P
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location(column%longitude,column%latitude,z_grid(k - bot + 1),3)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_U10
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
state_loc(dart_index) = set_location(column%longitude,column%latitude,10.0_r8,3)
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_V10
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
state_loc(dart_index) = set_location(column%longitude,column%latitude,10.0_r8,3)
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_T2
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
state_loc(dart_index) = set_location(column%longitude,column%latitude,2.0_r8,3)
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_Q2
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
state_loc(dart_index) = set_location(column%longitude,column%latitude,2.0_r8,3)
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_PSFC
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
state_loc(dart_index) = set_location(column%longitude,column%latitude,0.0_r8,3)
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TSLB
wrf_meta%var_size(var_cnt) = num_soil_layers
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location(column%longitude,column%latitude,0.0_r8,3)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TSK
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
state_loc(dart_index) = set_location(column%longitude,column%latitude,0.0_r8,3)
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

! here are the parameters appended to the state variable
do i = 1, wrf_meta%number_state_params
  wrf_meta%var_type(var_cnt) = wrf_meta%est_param_types(i)
  wrf_meta%var_size(var_cnt) = 1
  wrf_meta%var_index(var_cnt) = dart_index
  state_loc(dart_index) = set_location(column%longitude,column%latitude,0.0_r8,3)  !only surface here
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! all the stuff that we need to carry but don't want to affect
! THE ORDER IS IMPORTANT - make sure loop variables are correct
! 1D F vars
do vcnt = 1, wrf_meta%number_of_1d_f
  wrf_meta%var_type(var_cnt) = 100 + vcnt - 1 !TYPE_GSWF ... TYPE_PRECIPF
  select case ( wrf_meta%var_type(var_cnt) )
    case ( TYPE_GSWF )
      wrf_meta%var_size(var_cnt) = nsplinetimes_flux
    case ( TYPE_GLWF )
      wrf_meta%var_size(var_cnt) = nsplinetimes_flux
    case ( TYPE_PRECIPF )
      wrf_meta%var_size(var_cnt) = nsplinetimes_smos
    case default
      call error_handler(E_ERR, 'static_init_model', &
       'cannot size 1-d forcing type', source, revision, revdate)
  end select
  wrf_meta%var_index(var_cnt) = dart_index
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    state_loc(k) = set_location(long_far,column%latitude,(10000.0_r8+k),-1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! 2d F vars (nz,nsplinetimes)
do vcnt = 1, wrf_meta%number_of_2d_f
  wrf_meta%var_type(var_cnt) = 100 + wrf_meta%number_of_1d_f + vcnt - 1 !TYPE_UF ... TYPE_P8WF
  wrf_meta%var_size(var_cnt) = nsplinetimes*nz
  wrf_meta%var_index(var_cnt) = dart_index
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    state_loc(k) = set_location(long_far,column%latitude,(10000.0_r8+k),-1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! 2d advection-related vars (nz,nsplinetimes_advection)
do vcnt = 1, wrf_meta%number_of_2d_advection
  wrf_meta%var_type(var_cnt) = 100 + wrf_meta%number_of_1d_f &
          + wrf_meta%number_of_2d_f + vcnt - 1 !TYPE_TH_UPSTREAM_X...TYPE_TAU_V
  wrf_meta%var_size(var_cnt) = nsplinetimes_advection*nz
  wrf_meta%var_index(var_cnt) = dart_index
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    state_loc(k) = set_location(long_far,column%latitude,(10000.0_r8+k),-1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! noassim profile vars (nz)
do vcnt = 1, wrf_meta%number_noassim_profiles
  wrf_meta%var_type(var_cnt) = 500 + vcnt - 1 !RHO...EL
  wrf_meta%var_size(var_cnt) = nz
  wrf_meta%var_index(var_cnt) = dart_index
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    state_loc(k) = set_location(long_far,column%latitude,(10000.0_r8+k),1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! noassim soil vars (num_soil_layers)
do vcnt = 1, wrf_meta%number_noassim_soil_vars
  wrf_meta%var_type(var_cnt) = 200 + vcnt - 1 !TYPE_SMOIS ... TYPE_SH2O
  wrf_meta%var_size(var_cnt) = num_soil_layers
  wrf_meta%var_index(var_cnt) = dart_index
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    state_loc(k) = set_location(long_far,column%latitude,-(k - bot + 1.0_r8),1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! scalars (might want these in the assimilation vector)
do vcnt = 1, wrf_meta%number_noassim_scalars
  wrf_meta%var_type(var_cnt) = 300 + vcnt - 1 !TYPE_UZ0 ... TYPE_GAMMA
  wrf_meta%var_size(var_cnt) = 1
  wrf_meta%var_index(var_cnt) = dart_index
  state_loc(dart_index) = set_location(long_far,column%latitude,10000.0_r8,-1)
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! parameters that depend on something (strange ones)
do vcnt = 1, wrf_meta%number_dependent_params
  wrf_meta%var_type(var_cnt) = 400 + vcnt - 1 !TYPE_VEGFRA...?
  wrf_meta%var_size(var_cnt) = 1
  wrf_meta%var_index(var_cnt) = dart_index
  state_loc(dart_index) = set_location(long_far,column%latitude,10000.0_r8,-1)
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! map information 

call map_init(my_projection)

select case (init_f_type)

   case ('WRF')
      call get_projection(projcode,sw_corner_lat, sw_corner_lon,&
           center_i, center_j, &
           dx, stdlon, truelat1, truelat2)

      call map_set(projcode, sw_corner_lat, sw_corner_lon, center_i, center_j, &
                   dx, stdlon, truelat1, truelat2, my_projection)

   case ('OBS')
      print*,'Map information set to trivial for init type ',init_f_type
      my_projection%code = PROJ_LATLON
      
   case default
end select

END SUBROUTINE static_init_model

subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Initial conditions for PBL 1D model is achived via module_wrf
  implicit none

  real(r8), intent(out)  :: x(:)
  integer :: i

  call init_wrf(wrf_rnd_seed)

  call wrf_to_vector(x)

! add a tiny bit of noise to the initial conditions
!  do i = 1, wrf_meta%model_size
!    x(i) = x(i)+random_gaussian(pert_ran_seq,0.0,0.01*x(i))
!  enddo

! initialize any parameter distributions
  if ( wrf_meta%number_state_params > 0 ) then
     call init_params(x)
  endif

! count the ensemble members
  internal_ensemble_counter = internal_ensemble_counter + 1
  number_ensemble_members = max(number_ensemble_members,internal_ensemble_counter)

end subroutine init_conditions

!--------------------------------------------------------

subroutine wrf_to_vector(x)
!---------------------------------------------------------
! unrolls wrf state into dart vector
!---------------------------------------------------------
  implicit none

  real(r8), intent(out)  :: x(:)

  integer               :: dart_index, bot, top, var_cnt
  integer               :: ispline, iz, vcnt
  character(len=129)    :: errstring
  real(r8)              :: dartparam

  dart_index = 1
  var_cnt = 1

! u
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    x(k) = u_phy(1,k-bot+1,1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! v
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    x(k) = v_phy(1,k-bot+1,1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! z
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    x(k) = z(1,k-bot+1,1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! t
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    x(k) = t_phy(1,k-bot+1,1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! th
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    x(k) = th_phy(1,k-bot+1,1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! qv
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    x(k) = moist(1,k-bot+1,1,P_QV)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! tke
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    x(k) = tke_myj(1,k-bot+1,1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! p
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    x(k) = p_phy(1,k-bot+1,1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! u10
  x(wrf_meta%var_index(var_cnt)) = u10(1,1)
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! v10
  x(wrf_meta%var_index(var_cnt)) = v10(1,1)
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! t2
  x(wrf_meta%var_index(var_cnt)) = t2(1,1)
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! q2
  x(wrf_meta%var_index(var_cnt)) = q2(1,1)
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! psfc
  x(wrf_meta%var_index(var_cnt)) = psfc(1,1)
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! tslb
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    x(k) = tslb(1,k-bot+1,1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! tsk
  x(wrf_meta%var_index(var_cnt)) = tsk(1,1)
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! parameters
  do vcnt = 1, wrf_meta%number_state_params

    select case (wrf_meta%est_param_types(vcnt))
      case (TYPE_EMISS)
        dartparam = emiss(1,1)
      case (TYPE_ALBEDO)
        dartparam = albedo(1,1)
      case (TYPE_Z0)
        dartparam = z0(1,1)
      case (TYPE_THC)
        dartparam = thc(1,1)
      case (TYPE_MAVAIL)
        dartparam = mavail(1,1)
      case default
        call error_handler(E_ERR, 'wrf_to_vector', &
         'problem with est_param unroll', source, revision, revdate)
    end select

    x(wrf_meta%var_index(var_cnt)) = dartparam

    if ( wrf_meta%dist_shape(vcnt) == 'logn' ) &
         x(wrf_meta%var_index(var_cnt)) = dlog(dartparam)
    if ( wrf_meta%dist_shape(vcnt) == 'beta' ) &
         x(wrf_meta%var_index(var_cnt)) = dlog(dartparam/(1.0_r8-dartparam))

    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! all the far away ones
! gsw_f...glw_f
  do vcnt = 1, wrf_meta%number_of_1d_f
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_GSWF) 
          x(k) = gsw_f(k-bot+1) 
        case (TYPE_GLWF)
          x(k) = glw_f(k-bot+1) 
        case (TYPE_PRECIPF)
          x(k) = precip_f(k-bot+1) 
        case default
          call error_handler(E_ERR, 'wrf_to_vector', &
           'problem with 1d gen unroll', source, revision, revdate)
      end select
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! u_g_f...p8w_f
  do vcnt = 1, wrf_meta%number_of_2d_f
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    ispline = 1
    iz      = 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_UF)
          x(k) = u_g_f(iz,ispline) 
        case (TYPE_VF)
          x(k) = v_g_f(iz,ispline) 
        case (TYPE_PF)
          x(k) = p_f(iz,ispline) 
        case (TYPE_P8WF)
          x(k) = p8w_f(iz,ispline) 
        case default
          call error_handler(E_ERR, 'wrf_to_vector', &
           'problem with 2d gen unroll', source, revision, revdate)
      end select
      iz = iz + 1
      if ( iz > nz ) then
        ispline = ispline + 1
        iz = 1
      endif
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! th_upstream_x...tau_v
  do vcnt = 1, wrf_meta%number_of_2d_advection
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    ispline = 1
    iz      = 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_TH_UPSTREAM_X)
          x(k) = th_upstream_x_f(iz,ispline) 
        case (TYPE_TH_UPSTREAM_Y)
          x(k) = th_upstream_y_f(iz,ispline) 
        case (TYPE_QV_UPSTREAM_X)
          x(k) = qv_upstream_x_f(iz,ispline) 
        case (TYPE_QV_UPSTREAM_Y)
          x(k) = qv_upstream_y_f(iz,ispline) 
        case (TYPE_TAU_U)
          x(k) = tau_u_f(iz,ispline) 
        case (TYPE_TAU_V)
          x(k) = tau_v_f(iz,ispline) 
        case default
          call error_handler(E_ERR, 'wrf_to_vector', &
           'problem with 2d gen unroll', source, revision, revdate)
      end select
      iz = iz + 1
      if ( iz > nz ) then
        ispline = ispline + 1
        iz = 1
      endif
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo


! rho...el
  do vcnt = 1,wrf_meta%number_noassim_profiles
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_RHO)
          x(k) = rho(1,k-bot+1,1)
        case (TYPE_U_G)
          x(k) = u_g(k-bot+1)
        case (TYPE_V_G)
          x(k) = v_g(k-bot+1)
        case (TYPE_EL)
          x(k) = el_myj(1,k-bot+1,1)
        case default
          call error_handler(E_ERR, 'wrf_to_vector', &
           'problem with nossim_profile unroll', source, revision, revdate)
      end select
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! smois...sh2o
  do vcnt = 1,wrf_meta%number_noassim_soil_vars
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_SMOIS)
          x(k) = smois(1,k-bot+1,1)
        case (TYPE_KEEPFR3DFLAG)
          x(k) = keepfr3dflag(1,k-bot+1,1)
        case (TYPE_SMFR3D)
          x(k) = smfr3d(1,k-bot+1,1)
        case (TYPE_SH2O)
          x(k) = sh2o(1,k-bot+1,1)
        case default
          call error_handler(E_ERR, 'wrf_to_vector', &
           'problem with nossim_soil unroll', source, revision, revdate)
      end select
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! uz0...minm
  do vcnt = 1,wrf_meta%number_noassim_scalars
    select case (wrf_meta%var_type(var_cnt))
      case (TYPE_UZ0)
        x(wrf_meta%var_index(var_cnt)) = uz0(1,1)
      case (TYPE_VZ0)
        x(wrf_meta%var_index(var_cnt)) = vz0(1,1)
      case (TYPE_THZ0)
        x(wrf_meta%var_index(var_cnt)) = thz0(1,1)
      case (TYPE_QZ0)
        x(wrf_meta%var_index(var_cnt)) = qz0(1,1)
      case (TYPE_QVG)
        x(wrf_meta%var_index(var_cnt)) = qvg(1,1)
      case (TYPE_QSG)
        x(wrf_meta%var_index(var_cnt)) = qsg(1,1)
      case (TYPE_QCG)
        x(wrf_meta%var_index(var_cnt)) = qcg(1,1)
      case (TYPE_QSFC)
        x(wrf_meta%var_index(var_cnt)) = qsfc(1,1)
      case (TYPE_AKMS)
        x(wrf_meta%var_index(var_cnt)) = akms(1,1)
      case (TYPE_AKHS)
        x(wrf_meta%var_index(var_cnt)) = akhs(1,1)
      case (TYPE_HOL)
        x(wrf_meta%var_index(var_cnt)) = hol(1,1)
      case (TYPE_MOL)
        x(wrf_meta%var_index(var_cnt)) = mol(1,1)
      case (TYPE_RMOL)
        x(wrf_meta%var_index(var_cnt)) = rmol(1,1)
      case (TYPE_GRDFLX)
        x(wrf_meta%var_index(var_cnt)) = grdflx(1,1)
      case (TYPE_HFX)
        x(wrf_meta%var_index(var_cnt)) = hfx(1,1)
      case (TYPE_QFX)
        x(wrf_meta%var_index(var_cnt)) = qfx(1,1)
      case (TYPE_PSHLTR)
        x(wrf_meta%var_index(var_cnt)) = pshltr(1,1)
      case (TYPE_QSHLTR)
        x(wrf_meta%var_index(var_cnt)) = qshltr(1,1)
      case (TYPE_FLHC)
        x(wrf_meta%var_index(var_cnt)) = flhc(1,1)
      case (TYPE_FLQC)
        x(wrf_meta%var_index(var_cnt)) = flqc(1,1)
      case (TYPE_Q10)
        x(wrf_meta%var_index(var_cnt)) = q10(1,1)
      case (TYPE_UDRUNOFF)
        x(wrf_meta%var_index(var_cnt)) = udrunoff(1,1)
      case (TYPE_ACSNOM)
        x(wrf_meta%var_index(var_cnt)) = acsnom(1,1)
      case (TYPE_ACSNOW)
        x(wrf_meta%var_index(var_cnt)) = acsnow(1,1)
      case (TYPE_UST)
        x(wrf_meta%var_index(var_cnt)) = ust(1,1)
      case (TYPE_PBLH)
        x(wrf_meta%var_index(var_cnt)) = pblh(1,1)
      case (TYPE_TMN)
        x(wrf_meta%var_index(var_cnt)) = tmn(1,1)
      case (TYPE_CS)
        x(wrf_meta%var_index(var_cnt)) = cs(1,1)
      case (TYPE_ERATE)
        x(wrf_meta%var_index(var_cnt)) = erate(1,1)
      case (TYPE_MAXM)
        x(wrf_meta%var_index(var_cnt)) = maxm(1,1)
      case (TYPE_MINM)
        x(wrf_meta%var_index(var_cnt)) = minm(1,1)
      case (TYPE_ZNT)
        x(wrf_meta%var_index(var_cnt)) = znt(1,1)
      case (TYPE_GAMMA)
        x(wrf_meta%var_index(var_cnt)) = h_gamma
      case default
        call error_handler(E_ERR, 'wrf_to_vector', &
         'problem with noassim_scalar unroll', source, revision, revdate)
    end select
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! vegfra...
  do vcnt = 1,wrf_meta%number_dependent_params
    select case (wrf_meta%var_type(var_cnt))
      case (TYPE_VEGFRA)
        x(wrf_meta%var_index(var_cnt)) = vegfra(1,1)
      case default
        call error_handler(E_ERR, 'wrf_to_vector', &
         'problem with dependent_params unroll', source, revision, revdate)
    end select
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

  if ( dart_index-1 /= wrf_meta%model_size ) then
    errstring = "did not successfully unroll the dart vector"
    call error_handler(E_ERR,'wrf_to_vector', errstring, &
         source, revision, revdate)
  endif

end subroutine wrf_to_vector

subroutine adv_1step(x, dart_time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does single time step advance for PBL model - this is a call to wrf.F
! The Time argument is needed for compatibility with more complex models
! that need to know the time to compute their time tendency and is not
! used in L04. Is there a better way to do this in F90 than to just hang
! this argument out everywhere?
   implicit none

   real(r8), intent(inout) :: x(:)
   type(time_type), intent(inout) :: dart_time

   type(time_type)          :: time_into_forecast
   integer                  :: dart_seconds,dart_days

! perturb the parameters at every time step
   call pert_params_time(x)

   call impose_param_limits(x)

   call vector_to_wrf(x)

   time_into_forecast = dart_time - initialization_time
   call get_time(time_into_forecast,dart_seconds,dart_days)

! figure out what time step we are at
   call wrf(dart_seconds,dart_days)

!   call output_wrf_profiles()
!   print*,t2,t_phy(1,1,1),t_phy(1,10,1),t_phy(1,20,1)

   call wrf_to_vector(x)


end subroutine adv_1step

subroutine vector_to_wrf(x)
!---------------------------------------------------------
! rolls state vector into wrf state 
!---------------------------------------------------------
  implicit none

  real(r8), intent(in)  :: x(:)

  integer               :: dart_index, var_cnt, bot, top
  integer               :: ispline, iz, vcnt
  character(len=129)    :: errstring
  real(r8)              :: wrfparam

  dart_index = 1
  var_cnt = 1

! u
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    u_phy(1,k-bot+1,1) = x(k)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! v
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    v_phy(1,k-bot+1,1) = x(k)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! z
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    z(1,k-bot+1,1) = x(k)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! t
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    t_phy(1,k-bot+1,1) = x(k)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! th
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    th_phy(1,k-bot+1,1) = x(k)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! qv
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    moist(1,k-bot+1,1,P_QV) = x(k)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
  where (moist < 0.0_r8) moist = 0.0_r8

! tke
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    tke_myj(1,k-bot+1,1) = x(k)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
  where (tke_myj < 0.01_r8) tke_myj = 0.01_r8

! p
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    p_phy(1,k-bot+1,1) = x(k)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! u10
  u10(1,1) = x(wrf_meta%var_index(var_cnt))
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! v10
  v10(1,1) = x(wrf_meta%var_index(var_cnt))
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! t2
  t2(1,1) = x(wrf_meta%var_index(var_cnt))
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! q2
  q2(1,1) = x(wrf_meta%var_index(var_cnt))
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
  where (q2 < 0.0_r8) q2 = 0.0_r8

! psfc
  psfc(1,1) = x(wrf_meta%var_index(var_cnt))
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
  where (q2 < 0.0_r8) q2 = 0.0_r8

! tslb
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    tslb(1,k-bot+1,1) = x(k)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! tsk
  tsk(1,1) = x(wrf_meta%var_index(var_cnt))
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1

! parameters
  do vcnt = 1, wrf_meta%number_state_params

    wrfparam = x(wrf_meta%var_index(var_cnt))

    if ( wrf_meta%dist_shape(vcnt) == 'logn' ) &
       wrfparam = dexp(x(wrf_meta%var_index(var_cnt)))
    if ( wrf_meta%dist_shape(vcnt) == 'beta' ) &
       wrfparam = dexp(x(wrf_meta%var_index(var_cnt)))/ &
          (1.0_r8+dexp(x(wrf_meta%var_index(var_cnt))))

    select case (wrf_meta%est_param_types(vcnt))
      case (TYPE_EMISS)
        emiss(1,1)  = wrfparam
      case (TYPE_ALBEDO)
        albedo(1,1) = wrfparam
      case (TYPE_Z0)
        z0(1,1)     = wrfparam
      case (TYPE_THC)
        thc(1,1)    = wrfparam
      case (TYPE_MAVAIL)
        mavail(1,1) = wrfparam
      case default
        call error_handler(E_ERR, 'vector_to_wrf', &
         'problem with est_param roll', source, revision, revdate)
    end select

    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! all the far away ones
! gsw_f...precip_f
  do vcnt = 1, wrf_meta%number_of_1d_f
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_GSWF) 
          gsw_f(k-bot+1) = x(k)
        case (TYPE_GLWF)
          glw_f(k-bot+1) = x(k)
        case (TYPE_PRECIPF)
          precip_f(k-bot+1) = x(k)
        case default
          call error_handler(E_ERR, 'vector_to_wrf', &
           'problem with 1d gen roll', source, revision, revdate)
      end select
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! u_g_f...p8w_f
  do vcnt = 1, wrf_meta%number_of_2d_f
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    ispline = 1
    iz      = 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_UF)
          u_g_f(iz,ispline) = x(k)
        case (TYPE_VF)
          v_g_f(iz,ispline) = x(k)
        case (TYPE_PF)
          p_f(iz,ispline) = x(k)
        case (TYPE_P8WF)
          p8w_f(iz,ispline) = x(k)
        case default
          call error_handler(E_ERR, 'vector_to_wrf', &
           'problem with 2d f roll', source, revision, revdate)
      end select
      iz = iz + 1
      if ( iz > nz ) then
        ispline = ispline + 1
        iz = 1
      endif
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! th_upstream_x...tau_u
  do vcnt = 1, wrf_meta%number_of_2d_advection
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    ispline = 1
    iz      = 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_TH_UPSTREAM_X)
          th_upstream_x_f(iz,ispline) = x(k)
        case (TYPE_TH_UPSTREAM_Y)
          th_upstream_y_f(iz,ispline) = x(k)
        case (TYPE_QV_UPSTREAM_X)
          qv_upstream_x_f(iz,ispline) = x(k)
        case (TYPE_QV_UPSTREAM_Y)
          qv_upstream_y_f(iz,ispline) = x(k)
        case (TYPE_TAU_U)
          tau_u_f(iz,ispline) = x(k)
        case (TYPE_TAU_V)
          tau_v_f(iz,ispline) = x(k)
        case default
          call error_handler(E_ERR, 'vector_to_wrf', &
           'problem with 2d f roll', source, revision, revdate)
      end select
      iz = iz + 1
      if ( iz > nz ) then
        ispline = ispline + 1
        iz = 1
      endif
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! rho...el
  do vcnt = 1,wrf_meta%number_noassim_profiles
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_RHO)
          rho(1,k-bot+1,1) = x(k)
        case (TYPE_U_G)
          u_g(k-bot+1) = x(k)
        case (TYPE_V_G)
          v_g(k-bot+1) = x(k)
        case (TYPE_EL)
          el_myj(1,k-bot+1,1) = x(k)
        case default
          call error_handler(E_ERR, 'vector_to_wrf', &
           'problem with nossim_profile roll', source, revision, revdate)
      end select
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! smois...sh2o
  do vcnt = 1,wrf_meta%number_noassim_soil_vars
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_SMOIS)
          smois(1,k-bot+1,1) = x(k)
        case (TYPE_KEEPFR3DFLAG)
          keepfr3dflag(1,k-bot+1,1) = x(k)
        case (TYPE_SMFR3D)
          smfr3d(1,k-bot+1,1) = x(k)
        case (TYPE_SH2O)
          sh2o(1,k-bot+1,1) = x(k)
        case default
          call error_handler(E_ERR, 'vector_to_wrf', &
           'problem with nossim_soil roll', source, revision, revdate)
      end select
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! uz0...tmn
  do vcnt = 1,wrf_meta%number_noassim_scalars
    select case (wrf_meta%var_type(var_cnt))
      case (TYPE_UZ0)
        uz0(1,1) = x(wrf_meta%var_index(var_cnt)) 
      case (TYPE_VZ0)
        vz0(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_THZ0)
        thz0(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_QZ0)
        qz0(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_QVG)
        qvg(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_QSG)
        qsg(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_QCG)
        qcg(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_QSFC)
        qsfc(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_AKMS)
        akms(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_AKHS)
        akhs(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_HOL)
        hol(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_MOL)
        mol(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_RMOL)
        rmol(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_GRDFLX)
        grdflx(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_HFX)
        hfx(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_QFX)
        qfx(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_PSHLTR)
        pshltr(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_QSHLTR)
        qshltr(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_FLHC)
        flhc(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_FLQC)
        flqc(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_Q10)
        q10(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_UDRUNOFF)
        udrunoff(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_ACSNOM)
        acsnom(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_ACSNOW)
        acsnow(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_UST)
        ust(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_PBLH)
        pblh(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_TMN)
        tmn(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_CS)
        cs(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_ERATE)
        erate(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_MAXM)
        maxm(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_MINM)
        minm(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_ZNT)
        znt(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_GAMMA)
        h_gamma = x(wrf_meta%var_index(var_cnt))
      case default
        call error_handler(E_ERR, 'vector_to_wrf', &
         'problem with noassim_scalar roll', source, revision, revdate)
    end select
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! vegfra...
  do vcnt = 1,wrf_meta%number_dependent_params
    select case (wrf_meta%var_type(var_cnt))
      case (TYPE_VEGFRA)
        vegfra(1,1) = x(wrf_meta%var_index(var_cnt)) 
      case default
        call error_handler(E_ERR, 'vector_to_wrf', &
         'problem with dependent_params roll', source, revision, revdate)
    end select
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

  if ( dart_index-1 /= wrf_meta%model_size ) then
    errstring = "did not successfully roll from dart vector"
    call error_handler(E_ERR,'vector_to_wrf', errstring, &
         source, revision, revdate)
  endif

end subroutine vector_to_wrf


function get_model_size()
!------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model
implicit none

integer :: get_model_size

get_model_size = wrf_meta%model_size

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
!
! Gets the initial time for a state from the model. Where should this info
! come from in the most general case?

type(time_type), intent(out) :: time

integer                      :: start_seconds
character(len=129)  :: errstring

! based on namelist
select case ( init_f_type ) 
   case ('WRF')
      time = set_date(start_year_f, start_month_f,&
                      start_day_f, start_hour_f+int(start_forecast/3600), &
                      start_minute_f,0)
   case ('OBS')
      time = set_date(start_year_f, start_month_f,&
                      start_day_f, start_hour_f, &
                      start_minute_f,0)
   case default
      write(errstring,*) "Do not know how to initialize ",init_f_type
      call error_handler(E_ERR,'init_time', errstring, &
           source, revision, revdate)
end select

end subroutine init_time



subroutine model_interpolate(x, location, obs_kind, obs_val, istatus)
!------------------------------------------------------------------
!
! Interpolates from state vector x to the location. It's not particularly
! happy dumping all of this straight into the model. Eventually some
! concept of a grid underlying models but above locations is going to
! be more general. May want to wait on external infrastructure projects
! for this?

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_kind
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus
 
logical, parameter  :: debug = .false.
real(r8)            :: zloc_ind, zloc_val, xyz_loc(3)
integer             :: k, k1, k2
real(r8)            :: dz,dzm,lon_ref_deg
real(r8)            :: a1,utrue,vtrue,ugrid,vgrid
integer             :: in, ii, my_type
character(len=129)  :: errstring

real(r8), dimension(2) :: fld

! radians to degrees
lon_ref_deg = lon_ref / DEGRAD

! All forward operators supported   
istatus = 0

! Convert location to real
xyz_loc(:) = get_location(location)
zloc_val = xyz_loc(3)

! Don't care about surface pressure explicitly, so set location to 0.0
if ( vert_is_surface(location) ) zloc_val = 0.0_r8

if ( vert_is_level(location) .or. vert_is_surface(location) ) then      ! model level
  zloc_ind = zloc_val
elseif ( vert_is_pressure(location) ) then                      ! pressure
  istatus = 1
  errstring = "obs in pressure coords not yet implemented"
  call error_handler(E_ERR,'model_interpolate', errstring, &
       source, revision, revdate)
elseif ( vert_is_height(location) ) then                      ! height
  zloc_ind = 1
else
  istatus = 1
  errstring = "undefined obs vertical coords not yet implemented"
  call error_handler(E_ERR,'model_interpolate', errstring, &
       source, revision, revdate)
endif

if(zloc_ind /= missing_r8) then

k = max(1,int(zloc_ind))

! TJH attempt at logic to use generic obs kinds 
! Since there are no 'generic' kinds for most surface variables,
! we simply look at the generic kind and query the location to
! determine if the obs is at the surface ...

if      ( (obs_kind == KIND_U_WIND_COMPONENT) .and. (vert_is_surface(location) .or. zloc_val == 10.0) ) then
    my_type = TYPE_U10

else if ( (obs_kind == KIND_V_WIND_COMPONENT) .and. (vert_is_surface(location) .or. zloc_val == 10.0) ) then
    my_type = TYPE_V10

else if ( (obs_kind == KIND_TEMPERATURE)      .and. (vert_is_surface(location) .or. zloc_val == 2.0) ) then
    my_type = TYPE_T2

else if ( obs_kind == KIND_SPECIFIC_HUMIDITY  .and. (vert_is_surface(location) .or. zloc_val == 2.0) ) then 
    my_type = TYPE_Q2
else if ( obs_kind == KIND_TEMPERATURE ) then 
    my_type = TYPE_T
else if ( obs_kind == KIND_U_WIND_COMPONENT ) then 
    my_type = TYPE_U
else if ( obs_kind == KIND_V_WIND_COMPONENT ) then 
    my_type = TYPE_V
else if ( obs_kind == KIND_SPECIFIC_HUMIDITY ) then 
    my_type = TYPE_QV
else
    write(errstring,*) "No such obs kind in the state vector: ",obs_kind
    call error_handler(E_ERR,'model_interpolate', errstring, &
       source, revision, revdate)
endif

! If it is NOT a surface variable .... do something
! This could be replaced by the 'vert_is_surface' function.

if ( my_type .ne. TYPE_U10 .and. my_type .ne. TYPE_V10 .and. &
     my_type .ne. TYPE_T2  .and. my_type .ne. TYPE_Q2 ) then

  if ( my_type == TYPE_U .or. my_type == TYPE_V ) then
     do k2 = 1,2

       ugrid = x(get_wrf_index(k+k2-1,TYPE_U))
       vgrid = x(get_wrf_index(k+k2-1,TYPE_V))

       call gridwind_to_truewind(lon_ref_deg,my_projection,ugrid,vgrid,utrue,vtrue)

       if ( my_type == TYPE_U ) then
         fld(k2) = utrue
       else
         fld(k2) = vtrue
       endif
     enddo
  else ! must not be U or V

     k1 = get_wrf_index(k,my_type) 
     k2 = get_wrf_index(k+1,my_type) 
     fld(1) = x(k1)
     fld(2) = x(k2)
     ! convert from mixing ratio to specific humidity in here
     if ( my_type == TYPE_QV .and. obs_kind == KIND_SPECIFIC_HUMIDITY ) then
       fld(1) = fld(1) / ( 1.0_r8 + fld(1) )
       fld(2) = fld(2) / ( 1.0_r8 + fld(2) )
     endif

  endif !end if U or V

  call toGrid(zloc_ind, k, dz, dzm)

  if( k >= 1 ) then
     obs_val = dzm*fld(1) + dz*fld(2)
  elseif( k == 0 ) then
     obs_val = fld(1) - (fld(2)-fld(1))*dzm
  else
     obs_val = missing_r8
  endif

elseif ( my_type == TYPE_U10 .or. my_type == TYPE_V10 ) then !screen U,V

 ugrid = x(get_wrf_index(k,TYPE_U10))
 vgrid = x(get_wrf_index(k,TYPE_V10))

 call gridwind_to_truewind(lon_ref_deg,my_projection,ugrid,vgrid,utrue,vtrue)

 if ( my_type == TYPE_U10 ) then
   obs_val = utrue
 else
   obs_val = vtrue
 endif

elseif ( my_type == TYPE_T2 .or. my_type == TYPE_Q2 ) then ! screen T, Q

 k1 = get_wrf_index(k,my_type)
 obs_val = x(k1)
 ! convert from mixing ratio to specific humidity here
 if ( my_type == TYPE_Q2 .and. obs_kind == KIND_SPECIFIC_HUMIDITY ) then
   obs_val = obs_val / ( 1.0_r8 + obs_val )
 endif

else  !can't find it

   obs_val = missing_r8
   fld(:) = missing_r8

endif !end if type of obs

else  !don't have a valid vertical location

   obs_val = missing_r8
   fld(:) = missing_r8

endif !end if have a valid vertical location

end subroutine model_interpolate

!#######################################################################

subroutine toGrid (x, j, dx, dxm)

!  Transfer obs. x to grid j and calculate its
!  distance to grid j and j+1

  real(r8), intent(in)  :: x
  real(r8), intent(out) :: dx, dxm
  integer,  intent(out) :: j

  j = int (x)

  dx = x - real (j)

  dxm= 1.0_r8 - dx

end subroutine toGrid

!-----------------------------------------------------------------

function get_wrf_index( k,var_type )

integer, intent(in) :: k,var_type

integer :: get_wrf_index
integer :: in, ii

character(len=129) :: errstring

in = 0
do ii = 1, wrf_meta%total_number_of_vars
   if(var_type == wrf_meta%var_type(ii) ) in = ii
enddo

if ( k < 1 .or. k > wrf_meta%var_size(in) ) then
  write(errstring,*)'Index ',k,' exceeds grid dimension: ', &
       wrf_meta%var_size(in)
  call error_handler(E_ERR,'get_wrf_index', errstring, &
       source, revision, revdate)
endif

get_wrf_index = wrf_meta%var_index(in)-1 + k 

end function get_wrf_index

function get_model_time_step()
!------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?


integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type                                      

integer              :: i, bot, top
character(len=129)   :: errstring


location = state_loc(index_in)

if (present(var_type)) then
  var_type = -1
  do i = 1, wrf_meta%total_number_of_vars 
    bot = wrf_meta%var_index(i)
    top = bot + wrf_meta%var_size(i) - 1
    if ( index_in >= bot .and. index_in <=top ) then
      var_type = wrf_meta%var_type(i)
    endif
  enddo

  if ( var_type < 0 ) then
     errstring = "could not find variable type"
     call error_handler(E_ERR,'get_state_meta_data', errstring, &
                source, revision, revdate)
  endif
endif


end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Nothing for pbl_1d for now.


end subroutine end_model


subroutine model_get_close_states(o_loc, radius, inum, indices, dist, x)
!------------------------------------------------------------------
!
! Gets all states within radius of o_loc
   
implicit none

type(location_type), intent(in) :: o_loc
real(r8), intent(in) :: radius  
integer, intent(out) :: inum, indices(:)
real(r8), intent(out) :: dist(:)
real(r8), intent(in) :: x(:)

integer               :: i
integer               :: which_vert
real(r8)              :: obs_location(3)
real(r8)              :: state_location(3)
real(r8)              :: dist_tmp


which_vert = query_location(o_loc,'which_vert')
if ( which_vert > 1 ) then
     call error_handler(E_ERR, 'model_get_close_states; column', &
         'which_vert is invalid', source, revision, revdate)
endif

obs_location(:) = get_location(o_loc)

! Only two possibilities here: in the column or out of it.  Distance
! is entirely dependent on vertical location if it is in the column.
! Horizontal location is allowed a small tolerance to account for 
! precision.

inum = 0
do i = 1, wrf_meta%model_size
  state_location(:) = get_location(state_loc(i))
   if ( abs(obs_location(1) - state_location(1)) < 0.001 .and. (obs_location(2) - state_location(2)) < 0.001 ) then
      if ( which_vert == -1 ) then ! surface
        dist_tmp = state_location(3)
      elseif ( which_vert == 1 ) then ! model level
        dist_tmp = abs(obs_location(3)-state_location(3))
      endif
      if ( dist_tmp <= radius ) then
        inum = inum + 1
        indices(inum) = i
        dist(inum) = dist_tmp
      endif
   endif

enddo

end subroutine model_get_close_states


function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH Jan 24 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! JPH 22 Nov 04 -- for the PBL model, this is a whole bunch of info
! about specific model configs
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!--------------------------------------------------------------------
! General netCDF variables
!--------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!--------------------------------------------------------------------
! netCDF variables for Location
!--------------------------------------------------------------------

integer :: LocationVarID
integer :: StateVarDimID, StateVarVarID
integer :: StateVarID, MemberDimID, TimeDimID

!------------------------------------------
! same for physical space
!------------------------------------------
integer :: ZVarVarID, SLVarVarID, PVarVarID
integer :: ZVarDimID, SLVarDimID, PVarDimID
integer :: ZVarID(10), &
           SLVarID(2), ScrVarID(5), SkinVarID(5), SclrVarID(4)
! U, V, Z, T, RHO, QV, TKE, P, U_G, V_G
integer :: PVarID(wrf_meta%number_state_params) 

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1,str2

! strings for parameters
character(len=10)     :: pName, pUnits      
character(len=30)     :: pLongName      
character(len=129)   :: errstring

integer             :: i
type(location_type) :: lctn
ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_sync(ncFileID)) ! Ensure netCDF file is current
call check(nf90_Redef(ncFileID))

!--------------------------------------------------------------------
! Determine ID's from stuff already in the netCDF file
!--------------------------------------------------------------------

! make sure time is unlimited dimid

call check(nf90_inq_dimid(ncFileID,"copy",dimid=MemberDimID))
call check(nf90_inq_dimid(ncFileID,"time",dimid=TimeDimID))

!--------------------------------------------------------------------
! Write Global Attributes
!--------------------------------------------------------------------
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source", source ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision", revision ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate", revdate ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "PBL_1D"))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "PBL_type", bl_pbl_physics ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "sfclay_physics", sf_sfclay_physics ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "surface_physics", sf_surface_physics ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "dt", dt ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "dx", dx ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "deep_soil_moisture", deep_soil_moisture ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "mavail_ref", mavail_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "ivgtyp_ref", ivgtyp_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "isltyp_ref", isltyp_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "zo_ref", zo_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "emiss_ref", emiss_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "thc_ref", thc_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "albedo_ref", albedo_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "ts_ref", ts_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "tmn_ref", tmn_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "ps_ref", ps_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "lat_ref", lat_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "lon_ref", lon_ref ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "julday", julday ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "parameters_estimated", &
           wrf_meta%number_state_params ))

!--------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
! StateVariable is needed by some of our post-processing routines.
!--------------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                      len=wrf_meta%model_size, dimid = StateVarDimID))
call check(nf90_def_dim(ncid=ncFileID, name="z_level", &
                      len=nz, dimid = ZVarDimID))
call check(nf90_def_dim(ncid=ncFileID, name="sl_level", &
                      len=num_soil_layers, dimid = SLVarDimID))

!--------------------------------------------------------------------
! Define the Location Variable and add Attributes
! Some of the atts come from location_mod (via the USE: stmnt)
! CF standards for Locations:
! http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-working.html#ctype
!--------------------------------------------------------------------

if ( output_state_vector ) then
  call check(NF90_def_var(ncFileID, name=trim(adjustl(LocationName)), xtype=nf90_double, &
                dimids = StateVarDimID, varid=LocationVarID) )
  call check(nf90_put_att(ncFileID, LocationVarID, "long_name", trim(adjustl(LocationLName))))
  call check(nf90_put_att(ncFileID, LocationVarID, "dimension", LocationDims ))
  call check(nf90_put_att(ncFileID, LocationVarID, "units", "nondimensional"))
  call check(nf90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)))
else
endif

!--------------------------------------------------------------------
! Define either the "state vector" variables -OR- the "prognostic" variables.
!--------------------------------------------------------------------

if ( output_state_vector ) then
! Define the state vector coordinate variable
  call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
           dimids=StateVarDimID, varid=StateVarVarID))
  call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"))
  call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical") )
  call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, wrf_meta%model_size /)))

! Define the actual state vector
  call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_double, &
           dimids = (/ StateVarDimID, MemberDimID, TimeDimID /), varid=StateVarID))
  call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"))

else

  call check(NF90_def_var(ncFileID, name="z_level", xtype=nf90_int, &
                dimids = ZVarDimID, varid=ZVarVarID) )
  call check(NF90_def_var(ncFileID, name="sl_level", xtype=nf90_int, &
                dimids = SLVarDimID, varid=SLVarVarID) )
  call check(nf90_def_var(ncid=ncFileID, name="U", xtype=nf90_double, &
           dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(1)))
  call check(nf90_put_att(ncFileID, ZVarID(1), "long_name", "Zonal Wind"))
  call check(nf90_put_att(ncFileID, ZVarID(1), "units", "m/s"))
  call check(nf90_def_var(ncid=ncFileID, name="V", xtype=nf90_double, &
           dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(2)))
  call check(nf90_put_att(ncFileID, ZVarID(2), "long_name", "Meridional Wind"))
  call check(nf90_put_att(ncFileID, ZVarID(2), "units", "m/s"))
  call check(nf90_def_var(ncid=ncFileID, name="Z", xtype=nf90_double, &
           dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(3)))
  call check(nf90_put_att(ncFileID, ZVarID(3), "long_name", "Height AMSL"))
  call check(nf90_put_att(ncFileID, ZVarID(3), "units", "m"))
  call check(nf90_def_var(ncid=ncFileID, name="T", xtype=nf90_double, &
           dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(4)))
  call check(nf90_put_att(ncFileID, ZVarID(4), "long_name", "Temperature"))
  call check(nf90_put_att(ncFileID, ZVarID(4), "units", "K"))
  call check(nf90_def_var(ncid=ncFileID, name="RHO", xtype=nf90_double, &
           dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(5)))
  call check(nf90_put_att(ncFileID, ZVarID(5), "long_name", "Density"))
  call check(nf90_put_att(ncFileID, ZVarID(5), "units", "kg/m^3"))
  call check(nf90_def_var(ncid=ncFileID, name="QV", xtype=nf90_double, &
           dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(6)))
  call check(nf90_put_att(ncFileID, ZVarID(6), "long_name", "Mixing Ratio"))
  call check(nf90_put_att(ncFileID, ZVarID(6), "units", "kg/kg"))
  call check(nf90_def_var(ncid=ncFileID, name="TKE", xtype=nf90_double, &
           dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(7)))
  call check(nf90_put_att(ncFileID, ZVarID(7), "long_name", "Turb. Kinetic Energy"))
  call check(nf90_put_att(ncFileID, ZVarID(7), "units", "m2/s2"))
  call check(nf90_def_var(ncid=ncFileID, name="P", xtype=nf90_double, &
           dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(8)))
  call check(nf90_put_att(ncFileID, ZVarID(8), "long_name", "Pressure"))
  call check(nf90_put_att(ncFileID, ZVarID(8), "units", "Pa"))
  call check(nf90_def_var(ncid=ncFileID, name="U_G", xtype=nf90_double, &
           dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(9)))
  call check(nf90_put_att(ncFileID, ZVarID(9), "long_name", "U Forcing"))
  call check(nf90_put_att(ncFileID, ZVarID(9), "units", "m/s"))
  call check(nf90_def_var(ncid=ncFileID, name="V_G", xtype=nf90_double, &
           dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(10)))
  call check(nf90_put_att(ncFileID, ZVarID(10), "long_name", "V Forcing"))
  call check(nf90_put_att(ncFileID, ZVarID(10), "units", "m/s"))

  call check(nf90_def_var(ncid=ncFileID, name="TSLB", xtype=nf90_double, &
           dimids = (/ SLVarDimID, MemberDimID, TimeDimID /), varid=SLVarID(1)))
  call check(nf90_put_att(ncFileID, SLVarID(1), "long_name", "Soil Temperature"))
  call check(nf90_put_att(ncFileID, SLVarID(1), "units", "K"))
  call check(nf90_def_var(ncid=ncFileID, name="SMOIS", xtype=nf90_double, &
           dimids = (/ SLVarDimID, MemberDimID, TimeDimID /), varid=SLVarID(2)))
  call check(nf90_put_att(ncFileID, SLVarID(2), "long_name", "Soil Moisture"))
  call check(nf90_put_att(ncFileID, SLVarID(2), "units", "Volume Fraction"))

  call check(nf90_def_var(ncid=ncFileID, name="U10", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=ScrVarID(1)))
  call check(nf90_put_att(ncFileID, ScrVarID(1), "long_name", "10-m U-wind"))
  call check(nf90_put_att(ncFileID, ScrVarID(1), "units", "m/s"))
  call check(nf90_def_var(ncid=ncFileID, name="V10", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=ScrVarID(2)))
  call check(nf90_put_att(ncFileID, ScrVarID(2), "long_name", "10-m V-wind"))
  call check(nf90_put_att(ncFileID, ScrVarID(2), "units", "m/s"))
  call check(nf90_def_var(ncid=ncFileID, name="T2", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=ScrVarID(3)))
  call check(nf90_put_att(ncFileID, ScrVarID(3), "long_name", "2-m T"))
  call check(nf90_put_att(ncFileID, ScrVarID(3), "units", "m/s"))
  call check(nf90_def_var(ncid=ncFileID, name="Q2", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=ScrVarID(4)))
  call check(nf90_put_att(ncFileID, ScrVarID(4), "long_name", "2-m Mixing Ratio"))
  call check(nf90_put_att(ncFileID, ScrVarID(4), "units", "kg/kg"))
  call check(nf90_def_var(ncid=ncFileID, name="UST", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=ScrVarID(5)))
  call check(nf90_put_att(ncFileID, ScrVarID(5), "long_name", "U-star"))
  call check(nf90_put_att(ncFileID, ScrVarID(5), "units", "m/s"))

  call check(nf90_def_var(ncid=ncFileID, name="HFX", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=SkinVarID(1)))
  call check(nf90_put_att(ncFileID, SkinVarID(1), "long_name", "surface sensible heat flux"))
  call check(nf90_put_att(ncFileID, SkinVarID(1), "units", "W/m^2"))
  call check(nf90_def_var(ncid=ncFileID, name="QFX", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=SkinVarID(2)))
  call check(nf90_put_att(ncFileID, SkinVarID(2), "long_name", "surface moisture flux"))
  call check(nf90_put_att(ncFileID, SkinVarID(2), "units", "kg/m^2/s"))
  call check(nf90_def_var(ncid=ncFileID, name="TSK", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=SkinVarID(3)))
  call check(nf90_put_att(ncFileID, SkinVarID(3), "long_name", "skin temperature"))
  call check(nf90_put_att(ncFileID, SkinVarID(3), "units", "K"))
  call check(nf90_def_var(ncid=ncFileID, name="QVG", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=SkinVarID(4)))
  call check(nf90_put_att(ncFileID, SkinVarID(4), "long_name", "surface mixing ratio"))
  call check(nf90_put_att(ncFileID, SkinVarID(4), "units", "kg/kg"))
  call check(nf90_def_var(ncid=ncFileID, name="PRECIP", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=SkinVarID(5)))
  call check(nf90_put_att(ncFileID, SkinVarID(5), "long_name", "Accumulated precip"))
  call check(nf90_put_att(ncFileID, SkinVarID(5), "units", "m"))

  call check(nf90_def_var(ncid=ncFileID, name="PBLH", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=SclrVarID(1)))
  call check(nf90_put_att(ncFileID, SclrVarID(1), "long_name", "PBL height"))
  call check(nf90_put_att(ncFileID, SclrVarID(1), "units", "m"))
  call check(nf90_def_var(ncid=ncFileID, name="z_o", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=SclrVarID(2)))
  call check(nf90_put_att(ncFileID, SclrVarID(2), "long_name", "Surface roughness for momentum"))
  call check(nf90_put_att(ncFileID, SclrVarID(2), "units", "m"))
  call check(nf90_def_var(ncid=ncFileID, name="z_t", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=SclrVarID(3)))
  call check(nf90_put_att(ncFileID, SclrVarID(3), "long_name", "Surface roughness for heat"))
  call check(nf90_put_att(ncFileID, SclrVarID(3), "units", "m"))

  call check(nf90_def_var(ncid=ncFileID, name="z_q", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=SclrVarID(4)))
  call check(nf90_put_att(ncFileID, SclrVarID(4), "long_name", "Surface roughness for moisture"))
  call check(nf90_put_att(ncFileID, SclrVarID(4), "units", "m"))

  ! and state parameters
  do i = 1, wrf_meta%number_state_params
    select case (wrf_meta%est_param_types(i))
      case (TYPE_EMISS)
        pName = "EMISS"
        pLongName = "Surface emissivity"
        pUnits = "fraction"
      case (TYPE_ALBEDO)
        pName = "ALBEDO"
        pLongName = "Surface albedo"
      case (TYPE_Z0)
        pName = "Z0"
        pLongName = "Surface roughness"
        pUnits = "cm"
      case (TYPE_THC)
        pName = "THC"
        pLongName = "Surface thermal heat capacity"
        pUnits = ""
      case (TYPE_MAVAIL)
        pName = "MAVAIL"
        pLongName = "Surface moisture availability"
        pUnits = ""
      case default
        write(errstring,*) "No such parameter in the state vector"
        call error_handler(E_ERR,'nc_write_model_atts', errstring, &
           source, revision, revdate)
    end select
    call check(nf90_def_var(ncid=ncFileID, name=trim(pName), xtype=nf90_double,&
             dimids = (/ MemberDimID, TimeDimID /), varid=PVarID(i)))
    call check(nf90_put_att(ncFileID, PVarID(i), "long_name", trim(pLongName)))
    call check(nf90_put_att(ncFileID, PVarID(i), "units", trim(pUnits)))
  enddo

endif

! Leave define mode so we can fill
call check(nf90_enddef(ncfileID))

! Fill the state variable coordinate variable
if ( output_state_vector ) then
  call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,wrf_meta%model_size) /) ))
else
  call check(nf90_put_var(ncFileID, ZVarVarID, (/ (i,i=1,nz) /) ))
  call check(nf90_put_var(ncFileID, SLVarVarID, (/ (i,i=1,num_soil_layers) /) ))

endif

!--------------------------------------------------------------------
! Fill the location variable
!--------------------------------------------------------------------

if ( output_state_vector ) then
  do i = 1,wrf_meta%model_size
     call get_state_meta_data(i,lctn)
     call check(nf90_put_var(ncFileID, LocationVarID, get_location(lctn), (/ i /) ))
  enddo
else
endif

!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------
call check(nf90_sync(ncFileID))

write (*,*)'Model attributes written, netCDF file synched ...'

contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
      integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_atts',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check
end function nc_write_model_atts




function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH 23 May 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_04 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!-------------------------------------------------------------------------------
! General netCDF variables
!-------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, ZVarID(10), &
           ScrVarID(5), SLVarID(2), SkinVarID(5), SclrVarID(4)
integer :: PVarID(wrf_meta%number_state_params)

!-------------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------------

character(len=129)   :: errstring
character(len=10)  :: pName
ierr = 0                      ! assume normal termination

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))


if ( output_state_vector ) then
  call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
  call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
               start=(/ 1, copyindex, timeindex /)))
else
  call vector_to_wrf(statevec)

  call check(NF90_inq_varid(ncFileID, "U", ZVarID(1)) )
  call check(NF90_put_var(ncFileID, ZVarID(1), u_phy(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "V", ZVarID(2)) )
  call check(NF90_put_var(ncFileID, ZVarID(2), v_phy(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "Z", ZVarID(3)) )
  call check(NF90_put_var(ncFileID, ZVarID(3), z(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "T", ZVarID(4)) )
  call check(NF90_put_var(ncFileID, ZVarID(4), t_phy(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "RHO", ZVarID(5)) )
  call check(NF90_put_var(ncFileID, ZVarID(5), rho(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "QV", ZVarID(6)) )
  call check(NF90_put_var(ncFileID, ZVarID(6), moist(1,:,1,P_QV),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "TKE", ZVarID(7)) )
  call check(NF90_put_var(ncFileID, ZVarID(7), tke_myj(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "P", ZVarID(8)) )
  call check(NF90_put_var(ncFileID, ZVarID(8), p_phy(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "U_G", ZVarID(9)) )
  call check(NF90_put_var(ncFileID, ZVarID(9), u_g(:),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "V_G", ZVarID(10)) )
  call check(NF90_put_var(ncFileID, ZVarID(10), v_g(:),  &
               start=(/ 1, copyindex, timeindex /)))

  call check(NF90_inq_varid(ncFileID, "TSLB", SLVarID(1)) )
  call check(NF90_put_var(ncFileID, SLVarID(1), tslb(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "SMOIS", SLVarID(2)) )
  call check(NF90_put_var(ncFileID, SLVarID(2), smois(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))

  call check(NF90_inq_varid(ncFileID, "U10", ScrVarID(1)) )
  call check(NF90_put_var(ncFileID, ScrVarID(1), u10(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "V10", ScrVarID(2)) )
  call check(NF90_put_var(ncFileID, ScrVarID(2), v10(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "T2", ScrVarID(3)) )
  call check(NF90_put_var(ncFileID, ScrVarID(3), t2(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "Q2", ScrVarID(4)) )
  call check(NF90_put_var(ncFileID, ScrVarID(4), q2(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "UST", ScrVarID(5)) )
  call check(NF90_put_var(ncFileID, ScrVarID(5), ust(1,1),  &
               start=(/ copyindex, timeindex /)))

  call check(NF90_inq_varid(ncFileID, "HFX", SkinVarID(1)) )
  call check(NF90_put_var(ncFileID, SkinVarID(1), hfx(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "QFX", SkinVarID(2)) )
  call check(NF90_put_var(ncFileID, SkinVarID(2), qfx(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "TSK", SkinVarID(3)) )
  call check(NF90_put_var(ncFileID, SkinVarID(3), tsk(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "QVG", SkinVarID(4)) )
  call check(NF90_put_var(ncFileID, SkinVarID(4), qvg(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "PRECIP", SkinVarID(5)) )
  call check(NF90_put_var(ncFileID, SkinVarID(5), rainncv,  &
               start=(/ copyindex, timeindex /)))

  call check(NF90_inq_varid(ncFileID, "PBLH", SclrVarID(1)) )
  call check(NF90_put_var(ncFileID, SclrVarID(1), pblh(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "z_o", SclrVarID(2)) )
  call check(NF90_put_var(ncFileID, SclrVarID(2), z_o,  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "z_t", SclrVarID(3)) )
  call check(NF90_put_var(ncFileID, SclrVarID(3), z_t,  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "z_q", SclrVarID(4)) )
  call check(NF90_put_var(ncFileID, SclrVarID(4), z_q,  &
               start=(/ copyindex, timeindex /)))

  ! and state parameters
  do i = 1, wrf_meta%number_state_params
    select case (wrf_meta%est_param_types(i))
      case (TYPE_EMISS)
        pName = "EMISS"
        call check(NF90_inq_varid(ncFileID, trim(pName), PVarID(i)) )
        call check(NF90_put_var(ncFileID, PVarID(i), emiss(1,1),  &
                     start=(/ copyindex, timeindex /)))
      case (TYPE_ALBEDO)
        pName = "ALBEDO"
        call check(NF90_inq_varid(ncFileID, trim(pName), PVarID(i)) )
        call check(NF90_put_var(ncFileID, PVarID(i), albedo(1,1),  &
                     start=(/ copyindex, timeindex /)))
      case (TYPE_Z0)
        pName = "Z0"
        call check(NF90_inq_varid(ncFileID, trim(pName), PVarID(i)) )
        call check(NF90_put_var(ncFileID, PVarID(i), z0(1,1),  &
                     start=(/ copyindex, timeindex /)))
      case (TYPE_THC)
        pName = "THC"
        call check(NF90_inq_varid(ncFileID, trim(pName), PVarID(i)) )
        call check(NF90_put_var(ncFileID, PVarID(i), thc(1,1),  &
                     start=(/ copyindex, timeindex /)))
      case (TYPE_MAVAIL)
        pName = "MAVAIL"
        call check(NF90_inq_varid(ncFileID, trim(pName), PVarID(i)) )
        call check(NF90_put_var(ncFileID, PVarID(i), mavail(1,1),  &
                     start=(/ copyindex, timeindex /)))
      case default
        write(errstring,*) "No such parameter in the state vector"
        call error_handler(E_ERR,'nc_write_model_vars', errstring, &
           source, revision, revdate)
    end select
 enddo

endif

! write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
! write (*,*)'netCDF file is synched ...'

contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
   integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_vars',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check
end function nc_write_model_vars

!----------------------------------------------------------------------

function nc_read_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!------------------------------------------------------------------
! Reads the state vector 
! JPH
!

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(out):: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!-------------------------------------------------------------------------------
! General netCDF variables
!-------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, ZVarID(10), &
           ScrVarID(5),SLVarID(2),SkinVarID(4),SclrVarID(4)

!-------------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------------

ierr = 0                      ! assume normal termination

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

if ( output_state_vector ) then
  call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
  call check(NF90_get_var(ncFileID, StateVarID, statevec,  &
               start=(/ 1, copyindex, timeindex /)))
else

  call check(NF90_inq_varid(ncFileID, "U", ZVarID(1)) )
  call check(NF90_get_var(ncFileID, ZVarID(1), u_phy(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "V", ZVarID(2)) )
  call check(NF90_get_var(ncFileID, ZVarID(2), v_phy(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "Z", ZVarID(3)) )
  call check(NF90_get_var(ncFileID, ZVarID(3), z(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "T", ZVarID(4)) )
  call check(NF90_get_var(ncFileID, ZVarID(4), t_phy(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "RHO", ZVarID(5)) )
  call check(NF90_get_var(ncFileID, ZVarID(5), rho(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "QV", ZVarID(6)) )
  call check(NF90_get_var(ncFileID, ZVarID(6), moist(1,:,1,P_QV),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "TKE", ZVarID(7)) )
  call check(NF90_get_var(ncFileID, ZVarID(7), tke_myj(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "P", ZVarID(8)) )
  call check(NF90_get_var(ncFileID, ZVarID(8), p_phy(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "U_G", ZVarID(9)) )
  call check(NF90_get_var(ncFileID, ZVarID(9), u_g(:),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "V_G", ZVarID(10)) )
  call check(NF90_get_var(ncFileID, ZVarID(10), v_g(:),  &
               start=(/ 1, copyindex, timeindex /)))

  call check(NF90_inq_varid(ncFileID, "TSLB", SLVarID(1)) )
  call check(NF90_get_var(ncFileID, SLVarID(1), tslb(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "SMOIS", SLVarID(2)) )
  call check(NF90_get_var(ncFileID, SLVarID(2), smois(1,:,1),  &
               start=(/ 1, copyindex, timeindex /)))

  call check(NF90_inq_varid(ncFileID, "U10", ScrVarID(1)) )
  call check(NF90_get_var(ncFileID, ScrVarID(1), u10(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "V10", ScrVarID(2)) )
  call check(NF90_get_var(ncFileID, ScrVarID(2), v10(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "T2", ScrVarID(3)) )
  call check(NF90_get_var(ncFileID, ScrVarID(3), t2(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "Q2", ScrVarID(4)) )
  call check(NF90_get_var(ncFileID, ScrVarID(4), q2(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "UST", ScrVarID(5)) )
  call check(NF90_get_var(ncFileID, ScrVarID(5), ust(1,1),  &
               start=(/ copyindex, timeindex /)))

  call check(NF90_inq_varid(ncFileID, "HFX", SkinVarID(1)) )
  call check(NF90_get_var(ncFileID, SkinVarID(1), hfx(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "QFX", SkinVarID(2)) )
  call check(NF90_get_var(ncFileID, SkinVarID(2), qfx(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "TSK", SkinVarID(3)) )
  call check(NF90_get_var(ncFileID, SkinVarID(3), tsk(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "QVG", SkinVarID(4)) )
  call check(NF90_get_var(ncFileID, SkinVarID(4), qvg(1,1),  &
               start=(/ copyindex, timeindex /)))

  call check(NF90_inq_varid(ncFileID, "PBLH", SclrVarID(1)) )
  call check(NF90_get_var(ncFileID, SclrVarID(1), pblh(1,1),  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "z_o", SclrVarID(2)) )
  call check(NF90_get_var(ncFileID, SclrVarID(2), z_o,  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "z_t", SclrVarID(3)) )
  call check(NF90_get_var(ncFileID, SclrVarID(3), z_t,  &
               start=(/ copyindex, timeindex /)))
  call check(NF90_inq_varid(ncFileID, "z_q", SclrVarID(4)) )
  call check(NF90_get_var(ncFileID, SclrVarID(4), z_q,  &
               start=(/ copyindex, timeindex /)))


  call wrf_to_vector(statevec)
endif
! write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
! write (*,*)'netCDF file is synched ...'

contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
   integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_read_model_vars',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end function nc_read_model_vars

!----------------------------------------------------------------------

subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
! subroutine pert_model_state(state, pert_state, interf_provided)
!
! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

interf_provided = .true.

call init_conditions(pert_state)

end subroutine pert_model_state

SUBROUTINE get_projection(projcode,sw_corner_lat, sw_corner_lon, &
                          center_i, center_j, &
                          dx, stdlon, truelat1, truelat2)

   use netcdf
! Gets grid-related global attributes from input file

   integer, intent(out)  :: projcode
   real(r8), intent(out) :: center_i, center_j, truelat1, truelat2
   real(r8), intent(out) :: sw_corner_lat, sw_corner_lon, dx, stdlon

   integer  :: istatus, ncID

   ! these are unknown and unnecessary
   center_i = -999.
   center_j = -999.

   ! get the rest from the file
   call check(nf90_open(trim(indir)//'/'//init_f_file,NF90_NOWRITE,ncID))
 
   call check(nf90_get_att(ncid,NF90_GLOBAL,'MAP_PROJ',projcode))
   call check(nf90_get_att(ncid,NF90_GLOBAL,'SW_LAT',sw_corner_lat))
   call check(nf90_get_att(ncid,NF90_GLOBAL,'SW_LON',sw_corner_lon))
   call check(nf90_get_att(ncid,NF90_GLOBAL,'DX',dx))
   call check(nf90_get_att(ncid,NF90_GLOBAL,'CEN_LON',stdlon))
   call check(nf90_get_att(ncid,NF90_GLOBAL,'TRUELAT1',truelat1))
   call check(nf90_get_att(ncid,NF90_GLOBAL,'TRUELAT2',truelat2))
   call check(nf90_close(ncid))

   RETURN
   contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
      integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'get_projection',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check
END SUBROUTINE get_projection

!----------------------------------------------------------------------
!************* Miscellaneous functions/subroutines to do specific jobs
!----------------------------------------------------------------------

subroutine init_params(x)
! Initializes parameters with a normal, lognormal, or beta distribution

   real(r8), intent(inout)  :: x(:)
 
   integer :: iP, xloc, iR, i, j
   real(r8) :: a,b
   real(r8)              :: sd, m, s
   character(len=129)    :: errstring
   real(r8)              :: logP, logS, mP
   real(r8), dimension(:), allocatable :: udist
   integer, parameter    :: nseek=100
   integer, dimension(2) :: mloc
   real(r8), dimension(nseek,nseek) :: dist 

!  if it is in the state vector and is a logn distribution, x should
!  already be log(x)).  Likewise for beta.
   do iP = 1, wrf_meta%number_state_params

      xloc = get_wrf_index(1,wrf_meta%est_param_types(iP))

      select case ( wrf_meta%dist_shape(iP))

         case ( 'norm' )
            x(xloc) = x(xloc)+random_gaussian(param_ran_seq,0.0,wrf_meta%pert_init_sd(iP))
         case ( 'logn' )

            x(xloc) = dexp(x(xloc))

            x(xloc) = x(xloc)+random_gaussian(param_ran_seq,0.0,wrf_meta%pert_init_sd(iP))
            x(xloc) = dlog(x(xloc))

         case ( 'beta' ) 

            ! don't do anything if specified sd is 0.0 (better way to do this?)
            
            if ( wrf_meta%pert_init_sd(iP) > 0.0_r8 ) then

            select case (wrf_meta%est_param_types(iP))
              case (TYPE_EMISS)
                mP = emiss(1,1)  
              case (TYPE_ALBEDO)
                mP = albedo(1,1) 
              case (TYPE_Z0)
                mP = z0(1,1)    
              case (TYPE_THC)
                mP = thc(1,1)    
              case (TYPE_MAVAIL)
                mP = mavail(1,1) 
              case default
                call error_handler(E_ERR, 'init_params', &
                 'problem initializing parameter', source, revision, revdate)
            end select

!           If not specified, find something according to the mean and std
            if ( wrf_meta%pert_init_beta_1(iP) > 0.0_r8 ) then
               a = wrf_meta%pert_init_beta_1(iP)
               b = wrf_meta%pert_init_beta_2(iP)
            else
               b = 1.0_r8/wrf_meta%pert_init_sd(iP)
               a = mP*b / (1.0_r8-mP)
!               do i = 1, nseek
!                  do j = 1, nseek
!                    m = float(i)/float(i+j)
!                    s = float(i*j)/float(((i+j)**2)*(i+j+1))
!                    dist(i,j) = (m-mP)**2+(s-wrf_meta%pert_init_sd(iP)**2)
!                  enddo
!               enddo
!               mloc = minloc(dist)
!               a = mloc(1)
!               b = mloc(2)
            endif
            print*,'****USING BETA PARAMETERS: ',a,b
            print*,'****FOR mean, variance ',mP,wrf_meta%pert_init_sd(iP)
            allocate(udist(nint(a+b-1+0.5))) 
            do iR = 1, nint(a+b-1+0.5)
               udist(iR) = random_uniform(param_ran_seq)
            enddo
            udist = sort(udist)
            x(xloc) = dlog(udist(nint(a+0.5))/(1.0_r8-udist(nint(a+0.5))))
            deallocate(udist)

            endif ! don't do anything if sd = 0.0

         case default

            errstring = "Do not know how to initialize a "//wrf_meta%dist_shape(iP)//" distribution"
            call error_handler(E_ERR,'init_params', errstring, &
                    source, revision, revdate)

      end select

   enddo   

   call impose_param_limits(x)

end subroutine init_params

!----------------------------------------------------------------------

subroutine pert_params_time(x)
! Adds a little noise to any stoch parameters
! This is separate in case we want a different distribution from the
! initial perturbations.  

   real(r8), intent(inout)  :: x(:)
 
   integer :: iP, xloc
   real(r8)              :: sd, logP, logS
   character(len=129)    :: errstring

   do iP = 1, wrf_meta%number_state_params

      xloc = get_wrf_index(1,wrf_meta%est_param_types(iP))

         x(xloc) = x(xloc)+random_gaussian(param_ran_seq,0.0,wrf_meta%pert_param_sd(iP))

   enddo

end subroutine pert_params_time

!---------------------------------------------------------------------
subroutine adjust_param_spread(ens, nens)
! given the ensemble, find the parameters and adjust the spread to the initial
! values IF IT IS LESS.  In other words, this is only a lower bound
   integer, intent(in) :: nens
   real(r8), intent(inout) :: ens(:,:)
   integer :: iP, xloc, iEns
   real(r8)              :: logP(nens)
   real(r8)              :: mean
   real(r8)              :: std, std_new, stmp
   character(len=129)    :: errstring

   
   if ( maintain_initial_spread ) then ! should already be normal dist

      do iP = 1, wrf_meta%number_state_params

         xloc = get_wrf_index(1,wrf_meta%est_param_types(iP))
         std_new = wrf_meta%pert_init_sd(iP)

            mean = sum(ens(:,xloc))/nens
            std = sqrt(sum((ens(:,xloc)-mean)**2) / (nens-1))
            if ( std <= std_new ) then
               print*,'ADJUSTING SPREAD TO:::',std_new
               print*,'*****MIN/MAX before adjust*****',minval(ens(:,xloc)),maxval(ens(:,xloc))
               ens(:,xloc) = (ens(:,xloc)-mean) * (std_new/std) + mean
               print*,'*****MIN/MAX after adjust*****',minval(ens(:,xloc)),maxval(ens(:,xloc))
            endif
         
         ! constrain parameters 
         do iEns = 1, nens
            call impose_param_limits(ens(iEns,:)) 
         enddo
      enddo

    else
      return
    endif

end subroutine adjust_param_spread

!---------------------------------------------------------------------
subroutine impose_param_limits(x)
! Imposes max/min specified parameter limits.
! 
   real(r8), intent(inout) :: x(:)
   integer :: iP, xloc
   real(r8) :: dartparam_min, dartparam_max

   
   do iP = 1, wrf_meta%number_state_params

      xloc = get_wrf_index(1,wrf_meta%est_param_types(iP))

      ! constrain some parameters to [0,1]
      select case (wrf_meta%dist_shape(iP) )
         case ('norm')
           dartparam_min = wrf_meta%pert_param_min(iP)
           dartparam_max = wrf_meta%pert_param_max(iP)
           x(xloc) = min(x(xloc),dartparam_max)
           x(xloc) = max(x(xloc),dartparam_min)
         case ('logn')
           dartparam_min = dlog(wrf_meta%pert_param_min(iP))
           dartparam_max = dlog(wrf_meta%pert_param_max(iP))
           x(xloc) = min(x(xloc),dartparam_max)
           x(xloc) = max(x(xloc),dartparam_min)
         case ('beta')
           dartparam_min = dlog(wrf_meta%pert_param_min(iP)/ &
                        (1.0_r8-wrf_meta%pert_param_min(iP)))
           dartparam_max = dlog(wrf_meta%pert_param_max(iP)/ &
                        (1.0_r8-wrf_meta%pert_param_max(iP)))
           x(xloc) = min(x(xloc),dartparam_max)
           x(xloc) = max(x(xloc),dartparam_min)

         case default
            call error_handler (E_ERR,'impose_param_limits',&
                 'do not know distribution '//wrf_meta%dist_shape(iP), &
                 source, revision, revdate)
      end select

   enddo

end subroutine impose_param_limits

!-----------------------------------------------------------------

function get_pblh(x)
! Get diagnosed PBL height.

real(r8), intent(in)  :: x(:)

real(r8) :: get_pblh, pblh_m
real(r8), dimension(nz) :: z_tmp
integer  :: pblh_ind

pblh_m = x(get_wrf_index(1,TYPE_PBLH))

do i = 1, nz
  z_tmp(i) = x(get_wrf_index(i,TYPE_GZ))
enddo

get_pblh = -999.0_r8
do i = 1, nz-1
  if ( z_tmp(i) <= pblh_m .and. z_tmp(i+1) > pblh_m) then
    get_pblh = (pblh_m-z_tmp(i))/(z_tmp(i+1)-z_tmp(i)) + dble(i)
  endif
enddo 

if ( get_pblh < 0.0_r8 ) then
  if ( pblh_m < z_tmp(1) ) then
    get_pblh = pblh_m/z_tmp(1)
  else
    get_pblh = dble(nz)
  endif
endif

end function get_pblh

!-----------------------------------------------------------------

subroutine synchronize_mavail(obs,ens,ens_size,domain_size)

! argument ens is dimension(n_ensemble,n_state).  want to either update
! them all independently (they are all being nudged) or pick the first one
   integer,                  intent(in)    :: ens_size, domain_size
   real(r8),                 intent(in)    :: obs
   real(r8), dimension(ens_size,domain_size), intent(inout) :: ens
  
   integer                                 :: i
   real(r8), dimension(domain_size)         :: state_temp

! error check
   if ( domain_size /= wrf_meta%model_size ) stop 'error_model_size'

! loop is over ensemble members - could be removed to deal with first only
   do i = 1, ens_size
      state_temp = ens(i,:)
!  unwrap the state vector into model variables
      call vector_to_wrf(state_temp)

! here you have access to the model variables
!  solve the equation and change mavail
! need: pblh = h, h_gamma = h * gamma_q
! mavail, pblh are 2D horizontal arrays (1x1)
!
!      print*,obs,pblh(1,1),h_gamma,mavail(1,1)

!  wrap the model variables back into dart state vector
      call wrf_to_vector(state_temp)
      ens(i,:) = state_temp
    end do

end subroutine synchronize_mavail

!-----------------------------------------------------------------
 
function ok_to_nudge(obs_kind,state_index)
! checks for physical quantity consistency between obs and state
!  -- lots of hard coding
   integer, intent(in)  :: obs_kind, state_index
   logical         :: ok_to_nudge

   integer             :: var_type
   type(location_type) :: location

   ok_to_nudge = .false.
   call get_state_meta_data(state_index, location, var_type)

!  simple tests for pairs
   if ( obs_kind == 23 .and. var_type == TYPE_U   ) ok_to_nudge = .true.
   if ( obs_kind == 24 .and. var_type == TYPE_V   ) ok_to_nudge = .true.
   if ( obs_kind == 25 .and. var_type == TYPE_T   ) ok_to_nudge = .true.
   if ( obs_kind == 26 .and. var_type == TYPE_QV  ) ok_to_nudge = .true.
   
end function ok_to_nudge






subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

end subroutine ens_mean_for_model

!===================================================================
! End of model_mod
!===================================================================

end module model_mod

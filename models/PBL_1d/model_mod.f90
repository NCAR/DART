! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

use        types_mod, only : r8, missing_r8
use time_manager_mod, only : time_type, set_time, get_time, &
                             increment_time, print_time, &
                             set_calendar_type, NO_CALENDAR, &
                             operator(==), operator(<=)
use     location_mod, only : location_type, get_dist, set_location, &
                             get_location, query_location, &
                             LocationDims, LocationName, LocationLName, &
                             vert_is_surface, vert_is_pressure, &
                             vert_is_level, vert_is_height
use    utilities_mod, only : file_exist, open_file, close_file, &
                             register_module, error_handler, E_ERR, E_MSG, logfileunit
!use     obs_kind_mod, only : KIND_U, KIND_V, KIND_PS, KIND_T, KIND_QV, &
!                             KIND_P, KIND_W, KIND_QR, KIND_TD, KIND_RHO, &
!                             KIND_U10, KIND_V10, KIND_T2, &
!                             KIND_Q2, KIND_TD2
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
          get_pblh

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

! Basic model parameters controlled by namelist; no defaults

!-----------------------------------------
! namelists are declared in module_namelist
!-----------------------------------------

! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step
type(time_type) :: last_advance_time

! Private definition of model variable types

integer, parameter :: TYPE_U   = 1,   TYPE_V   = 2,  TYPE_W  = 3,  &
                      TYPE_GZ  = 4,   TYPE_T   = 5,  TYPE_MU = 6,  &
                      TYPE_QV  = 7,   TYPE_QC  = 8,  TYPE_QR = 9,  &
                      TYPE_QI  = 10,  TYPE_QS  = 11, TYPE_QG = 12, &
                      TYPE_U10 = 13,  TYPE_V10 = 14, TYPE_T2 = 15, &
                      TYPE_Q2  = 16,  TYPE_PS  = 17, TYPE_TSLB = 18, &
                      TYPE_TSK = 19,  TYPE_TH = 20,  TYPE_TKE = 21 , &
                      TYPE_P   = 22

! additional "far away" types that will not be part of the assimilation
! GROUP THESE BY SIZE OF ARRAYS!  this will make life easier below
integer, parameter :: TYPE_GSWGEN  = 100,    TYPE_GLWGEN = 101, & !1D       
                      TYPE_UGEN = 102,       TYPE_VGEN = 103, &   !2D...
                      TYPE_PGEN = 104,       TYPE_P8WGEN = 105, & 
                      TYPE_TGEN_UA = 106,    TYPE_TGEN_VA = 107, & 
                      TYPE_QGEN_UA = 108,    TYPE_QGEN_VA = 109, & 
                      TYPE_UGEN_UA = 110,    TYPE_UGEN_VA = 111, & 
                      TYPE_VGEN_UA = 112,    TYPE_VGEN_VA = 113, & 
                      TYPE_UGEN_WA = 114,    TYPE_VGEN_WA = 115, & 
                      TYPE_TGEN_WA = 116,    TYPE_QGEN_WA = 117

! far away profiles may or may not be part of the external forcing - perhaps
! diagnostic
integer, parameter :: TYPE_RHO = 500, TYPE_U_G  = 501, &
                      TYPE_V_G  = 502

! yet more far away types that we may want to add to the state vector at some
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
                      TYPE_PBLH = 324,       TYPE_TMN = 325

integer, parameter :: TYPE_VEGFRA = 400

integer, parameter :: calendar_type = NO_CALENDAR
integer            :: internal_ensemble_counter = 0
integer            :: number_ensemble_members = 0

TYPE state_vector_meta_data ! additional necessary vars
   integer             ::   model_size
   integer             ::   total_number_of_vars
! soil variables in state vector to dart
   integer             ::   number_of_soil_variables   
! the usual profile variables in state to dart
   integer             ::   number_of_profile_variables
! screen height vars in state vector to dart
   integer             ::   number_of_scalars         
! external forcing
   integer             ::   number_of_2d_gen 
   integer             ::   number_of_1d_gen 
! profile variables NOT in state vector to dart
   integer             ::   number_noassim_profiles
! soil variables NOT in state vector to dart
   integer             ::   number_noassim_soil_vars   
! scalars NOT in state vector to dart
   integer             ::   number_noassim_scalars
! parameters that depend on initialization date
   integer             ::   number_dependent_params

   integer, dimension(:), allocatable :: var_type
   integer, dimension(:), allocatable :: var_size
   integer, dimension(:), allocatable :: var_index
end TYPE state_vector_meta_data

type(state_vector_meta_data)    :: wrf_meta

type(proj_info)                 :: my_projection
!******my_projection is hard-coded below - put in wrf profile file?

! A flag for allocation, starts as true then goes to false with first
! call to init_wrf, which occurs in init_conditions.

integer                         :: wrf_rnd_seed
logical                         :: allocate_wrf = .true.

contains

!==================================================================


subroutine static_init_model()
!------------------------------------------------------------------
! Initializes class data for this model. For now, simply outputs the
! identity info, sets the location of the state variables, and initializes
! the time type for the time stepping (is this general enough for time???)
implicit none

real(r8) :: x_loc
integer :: dart_index, var_cnt, k, bot, top, unit_nml, vcnt, projcode
LOGICAL :: is_it_there = .FALSE.
REAL :: timeo,timetot
    
real(r8) :: center_i, center_j, truelat1, truelat2
real(r8) :: sw_corner_lat, sw_corner_lon, dx, stdlon

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Begin by reading the namelist input
if(file_exist('wrf1d_namelist.input')) then
   unit_nml = open_file(fname = 'wrf1d_namelist.input', action = 'read')
   call do_namelist_wrf1d(unit_nml,logfileunit)
   close(unit_nml)
endif

! initialize the random sequence
wrf_rnd_seed = rnd_seed_val

! need some dimension information
call static_init_wrf()

! initialize some timing stuff
time_step = set_time(int(dt), 0)
call set_calendar_type(calendar_type)

! initialize the last time
last_advance_time = set_time(start_forecast*3600,0)

! numbers
wrf_meta%number_of_scalars = 5          ! (U10, V10, T2, Q2, TSK)
wrf_meta%number_of_profile_variables = 8! (U, V, Z, T, QV, TH, TKE,P)
wrf_meta%number_of_soil_variables = 1   ! (TSLB)
wrf_meta%number_of_1d_gen = 2           ! (GSWGEN, GLWGEN)
wrf_meta%number_of_2d_gen = 16          !(UGEN, VGEN, PGEN, P8WGEN,
                                        ! TGEN_UA, TGEN_VA,
                                        ! QGEN_UA, QGEN_VA,
                                        ! UGEN_UA, UGEN_VA,
                                        ! VGEN_UA, VGEN_VA,
                                        ! UGEN_WA, VGEN_WA,
                                        ! TGEN_WA, QGEN_WA)
wrf_meta%number_noassim_profiles = 3    ! RHO, U_G, V_G
wrf_meta%number_noassim_soil_vars = 4   ! (SMOIS, KEEPFR3DFLAG, &
                                        !  SMFR3D, SH2O)
wrf_meta%number_noassim_scalars = 26    ! (UZ0, VZ0, THZ0, QZ0 &
                                        !  QVG, QSG, QCG, QSFC, &
                                        ! AKMS, AKHS, HOL, MOL, &
                                        ! GRDFLX, HFX, QFX, &
                                        ! PSHLTR, QSHLTR, &
                                        ! FLQC, FLHC, Q10, UDRUNOFF, &
                                        ! ACSNOM, ACSNOW, UST, PBLH, TMN)
wrf_meta%number_dependent_params = 1    ! (VEGFRA)
                                         

! compute the model size
wrf_meta%model_size = num_soil_layers * wrf_meta%number_of_soil_variables  &
                    + nz              * wrf_meta%number_of_profile_variables &
                    + 1               * wrf_meta%number_of_scalars &
                    + nz*nsplinetimes * wrf_meta%number_of_2d_gen &
                    + nsplinetimes    * wrf_meta%number_of_1d_gen &
                    + nz              * wrf_meta%number_noassim_profiles &
                    + num_soil_layers * wrf_meta%number_noassim_soil_vars &
                    + 1               * wrf_meta%number_noassim_scalars &
                    + 1               * wrf_meta%number_dependent_params

! state vector locations and wrf grid meta data
wrf_meta%total_number_of_vars = wrf_meta%number_of_scalars &
                              + wrf_meta%number_of_profile_variables &
                              + wrf_meta%number_of_soil_variables &
                              + wrf_meta%number_of_1d_gen &
                              + wrf_meta%number_of_2d_gen &
                              + wrf_meta%number_noassim_profiles &
                              + wrf_meta%number_noassim_soil_vars &
                              + wrf_meta%number_noassim_scalars &
                              + wrf_meta%number_dependent_params

allocate(wrf_meta%var_type(wrf_meta%total_number_of_vars))
allocate(wrf_meta%var_size(wrf_meta%total_number_of_vars))
allocate(wrf_meta%var_index(wrf_meta%total_number_of_vars))
allocate(state_loc(wrf_meta%model_size))

dart_index = 1
var_cnt    = 1

wrf_meta%var_type(var_cnt) = TYPE_U
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location((k - bot + 1.0_r8),1)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_V
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location((k - bot + 1.0_r8),1)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_GZ
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location((k - bot + 1.0_r8),1)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_T
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location((k - bot + 1.0_r8),1)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TH
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location((k - bot + 1.0_r8),1)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_QV
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location((k - bot + 1.0_r8),1)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TKE
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location((k - bot + 1.0_r8),1)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_P
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location((k - bot + 1.0_r8),1)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_U10
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
state_loc(dart_index) = set_location(0.0_r8,1)
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_V10
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
state_loc(dart_index) = set_location(0.0_r8,1)
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_T2
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
state_loc(dart_index) = set_location(0.0_r8,1)
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_Q2
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
state_loc(dart_index) = set_location(0.0_r8,1)
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TSLB
wrf_meta%var_size(var_cnt) = num_soil_layers
wrf_meta%var_index(var_cnt) = dart_index
bot = wrf_meta%var_index(var_cnt)
top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
do k = bot,top
  state_loc(k) = set_location(-(k - bot + 1.0_r8),1)
enddo
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TSK
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
state_loc(dart_index) = set_location(0.0_r8,1)
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

! all the stuff that we need to carry but don't want to affect
! THE ORDER IS IMPORTANT - make sure loop variables are correct
! 1D GEN vars
do vcnt = 1, wrf_meta%number_of_1d_gen
  wrf_meta%var_type(var_cnt) = 100 + vcnt - 1 !TYPE_GSWGEN ... TYPE_GLWGEN
  wrf_meta%var_size(var_cnt) = nsplinetimes
  wrf_meta%var_index(var_cnt) = dart_index
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    state_loc(k) = set_location((10000.0_r8+k),1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! 2d GEN vars (nz,nsplinetimes)
do vcnt = 1, wrf_meta%number_of_2d_gen
  wrf_meta%var_type(var_cnt) = 102 + vcnt - 1 !TYPE_UGEN ... TYPE_QGEN_WA
  wrf_meta%var_size(var_cnt) = nsplinetimes*nz
  wrf_meta%var_index(var_cnt) = dart_index
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    state_loc(k) = set_location((10000.0_r8+k),1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! noassim profile vars (nz)
do vcnt = 1, wrf_meta%number_noassim_profiles
  wrf_meta%var_type(var_cnt) = 500 + vcnt - 1 !RHO...V_G
  wrf_meta%var_size(var_cnt) = nz
  wrf_meta%var_index(var_cnt) = dart_index
  bot = wrf_meta%var_index(var_cnt)
  top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
  do k = bot,top
    state_loc(k) = set_location((10000.0_r8+k),1)
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
    state_loc(k) = set_location((10000.0_r8+k),1)
  enddo
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! scalars (might want these in the assimilation vector)
do vcnt = 1, wrf_meta%number_noassim_scalars
  wrf_meta%var_type(var_cnt) = 300 + vcnt - 1 !TYPE_UZ0 ... TYPE_TMN
  wrf_meta%var_size(var_cnt) = 1
  wrf_meta%var_index(var_cnt) = dart_index
  state_loc(dart_index) = set_location(10000.0_r8,1)
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! parameters that depend on something (strange ones)
do vcnt = 1, wrf_meta%number_dependent_params
  wrf_meta%var_type(var_cnt) = 400 + vcnt - 1 !TYPE_VEGFRA...?
  wrf_meta%var_size(var_cnt) = 1
  wrf_meta%var_index(var_cnt) = dart_index
  state_loc(dart_index) = set_location(10000.0_r8,1)
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! map information 

call map_init(my_projection)

call get_projection(projcode,sw_corner_lat, sw_corner_lon, center_i, center_j, &
     dx, stdlon, truelat1, truelat2)

!sw_corner_lat = 19.88053
!sw_corner_lon = -124.0885
!center_i      = 129.0
!center_j      = 82.0
!dx            = 22.e3
!stdlon        = -98.0
!truelat1      = 30.0
!truelat2      = 60.0

CALL map_set(projcode, sw_corner_lat, sw_corner_lon, center_i, center_j, &
     dx, stdlon, truelat1, truelat2, my_projection)


!center_i      = 129.0 
!center_j      = 82.0
!stdlon        = cent_lon
!
!IF (projcode==0) THEN
!   projcode=PROJ_LATLON
!ELSEIF (projcode==1) THEN 
!   projcode=PROJ_LC
!ELSEIF (projcode==2) THEN
!   projcode=PROJ_PS
!ELSEIF (projcode==3) THEN
!   projcode=PROJ_MERC
!ENDIF
!
!CALL map_set(projcode, sw_corner_lat, sw_corner_lon, &
!     &center_i, center_j, dx, &
!     &stdlon, truelat1, truelat2, my_projection)


END SUBROUTINE static_init_model

subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Initial conditions for PBL 1D model is achived via module_wrf
  implicit none

  real(r8), intent(out)  :: x(:)
  integer :: i

  call init_wrf(allocate_wrf,wrf_rnd_seed)

  allocate_wrf = .false.

  call wrf_to_vector(x)

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

! all the far away ones
! gsw_gen...glw_gen
  do vcnt = 1, wrf_meta%number_of_1d_gen
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_GSWGEN) 
          x(k) = gsw_gen(k-bot+1) 
        case (TYPE_GLWGEN)
          x(k) = glw_gen(k-bot+1) 
        case default
          call error_handler(E_ERR, 'wrf_to_vector', &
           'problem with 1d gen unroll', source, revision, revdate)
      end select
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! u_g_gen...qgen_wa
  do vcnt = 1, wrf_meta%number_of_2d_gen
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    ispline = 1
    iz      = 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_UGEN)
          x(k) = u_g_gen(iz,ispline) 
        case (TYPE_VGEN)
          x(k) = v_g_gen(iz,ispline) 
        case (TYPE_PGEN)
          x(k) = p_gen(iz,ispline) 
        case (TYPE_P8WGEN)
          x(k) = p8w_gen(iz,ispline) 
        case (TYPE_TGEN_UA)
          x(k) = t_gen_uadv(iz,ispline) 
        case (TYPE_TGEN_VA)
          x(k) = t_gen_vadv(iz,ispline) 
        case (TYPE_QGEN_UA)
          x(k) = q_gen_uadv(iz,ispline) 
        case (TYPE_QGEN_VA)
          x(k) = q_gen_vadv(iz,ispline) 
        case (TYPE_UGEN_UA)
          x(k) = u_gen_uadv(iz,ispline) 
        case (TYPE_UGEN_VA)
          x(k) = u_gen_vadv(iz,ispline) 
        case (TYPE_VGEN_UA)
          x(k) = v_gen_uadv(iz,ispline) 
        case (TYPE_VGEN_VA)
          x(k) = v_gen_vadv(iz,ispline) 
        case (TYPE_UGEN_WA)
          x(k) = u_gen_wadv(iz,ispline) 
        case (TYPE_VGEN_WA)
          x(k) = v_gen_wadv(iz,ispline) 
        case (TYPE_TGEN_WA)
          x(k) = t_gen_wadv(iz,ispline) 
        case (TYPE_QGEN_WA)
          x(k) = q_gen_wadv(iz,ispline) 
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

! rho...v_g
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

! uz0...tmn
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
      case (TYPE_TMN)
        x(wrf_meta%var_index(var_cnt)) = tmn(1,1)
      case (TYPE_PBLH)
        x(wrf_meta%var_index(var_cnt)) = pblh(1,1)
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

   integer                  :: dart_seconds,dart_days

! count ensemble member if needed (not now)
   if ( dart_time <= last_advance_time ) then
     internal_ensemble_counter = internal_ensemble_counter + 1
   endif
   if ( internal_ensemble_counter > number_ensemble_members ) then
     internal_ensemble_counter = 1
   endif

   call vector_to_wrf(x)

   call get_time(dart_time,dart_seconds,dart_days)

   call wrf(dart_seconds,dart_days)

!   print*,t2,t_phy(1,1,1),t_phy(1,10,1),t_phy(1,20,1)

   call wrf_to_vector(x)

   last_advance_time = dart_time

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

! all the far away ones
! gsw_gen...glw_gen
  do vcnt = 1, wrf_meta%number_of_1d_gen
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_GSWGEN) 
          gsw_gen(k-bot+1) = x(k)
        case (TYPE_GLWGEN)
          glw_gen(k-bot+1) = x(k)
        case default
          call error_handler(E_ERR, 'wrf_to_vector', &
           'problem with 1d gen roll', source, revision, revdate)
      end select
    enddo
    dart_index = dart_index + wrf_meta%var_size(var_cnt)
    var_cnt = var_cnt + 1
  enddo

! u_g_gen...qgen_wa
  do vcnt = 1, wrf_meta%number_of_2d_gen
    bot = wrf_meta%var_index(var_cnt)
    top = wrf_meta%var_index(var_cnt) + wrf_meta%var_size(var_cnt) - 1
    ispline = 1
    iz      = 1
    do k = bot,top
      select case (wrf_meta%var_type(var_cnt))
        case (TYPE_UGEN)
          u_g_gen(iz,ispline) = x(k)
        case (TYPE_VGEN)
          v_g_gen(iz,ispline) = x(k)
        case (TYPE_PGEN)
          p_gen(iz,ispline) = x(k)
        case (TYPE_P8WGEN)
          p8w_gen(iz,ispline) = x(k)
        case (TYPE_TGEN_UA)
          t_gen_uadv(iz,ispline) = x(k)
        case (TYPE_TGEN_VA)
          t_gen_vadv(iz,ispline) = x(k)
        case (TYPE_QGEN_UA)
          q_gen_uadv(iz,ispline) = x(k)
        case (TYPE_QGEN_VA)
          q_gen_vadv(iz,ispline) = x(k)
        case (TYPE_UGEN_UA)
          u_gen_uadv(iz,ispline) = x(k)
        case (TYPE_UGEN_VA)
          u_gen_vadv(iz,ispline) = x(k)
        case (TYPE_VGEN_UA)
          v_gen_uadv(iz,ispline) = x(k)
        case (TYPE_VGEN_VA)
          v_gen_vadv(iz,ispline) = x(k)
        case (TYPE_UGEN_WA)
          u_gen_wadv(iz,ispline) = x(k)
        case (TYPE_VGEN_WA)
          v_gen_wadv(iz,ispline) = x(k)
        case (TYPE_TGEN_WA)
          t_gen_wadv(iz,ispline) = x(k)
        case (TYPE_QGEN_WA)
          q_gen_wadv(iz,ispline) = x(k)
        case default
          call error_handler(E_ERR, 'wrf_to_vector', &
           'problem with 2d gen roll', source, revision, revdate)
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

! rho...v_g
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
        case default
          call error_handler(E_ERR, 'wrf_to_vector', &
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
          call error_handler(E_ERR, 'wrf_to_vector', &
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
      case (TYPE_TMN)
        tmn(1,1) = x(wrf_meta%var_index(var_cnt))
      case (TYPE_PBLH)
        pblh(1,1) = x(wrf_meta%var_index(var_cnt))
      case default
        call error_handler(E_ERR, 'wrf_to_vector', &
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
        call error_handler(E_ERR, 'wrf_to_vector', &
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

start_seconds = start_forecast * 3600
! based on namelist
time = set_time(start_seconds, 0)

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
real(r8)            :: zloc_ind, zloc_val
integer             :: k, k1, k2
real(r8)            :: dz,dzm,cent_lon_deg
real(r8)            :: a1,utrue,vtrue,ugrid,vgrid
integer             :: in, ii, my_type
character(len=129)  :: errstring

real(r8), dimension(2) :: fld

! radians to degrees
cent_lon_deg = cent_lon / DEGRAD

! All forward operators supported   
istatus = 0

! Convert location to real
zloc_val = get_location(location)

if ( vert_is_level(location) .or. vert_is_surface(location) ) then      ! model level
  zloc_ind = zloc_val
elseif ( vert_is_pressure(location) ) then                      ! pressure
  istatus = 1
  errstring = "obs in pressure coords not yet implemented"
  call error_handler(E_ERR,'model_interpolate', errstring, &
       source, revision, revdate)
elseif ( vert_is_height(location) ) then                      ! height
  istatus = 1
  errstring = "obs in height coords not yet implemented"
  call error_handler(E_ERR,'model_interpolate', errstring, &
       source, revision, revdate)
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

if      ( (obs_kind == KIND_U_WIND_COMPONENT) .and. vert_is_surface(location)) then
    my_type = TYPE_U10

else if ( (obs_kind == KIND_V_WIND_COMPONENT) .and. vert_is_surface(location)) then
    my_type = TYPE_V10

else if ( (obs_kind == KIND_TEMPERATURE)      .and. vert_is_surface(location)) then
    my_type = TYPE_T2

else if ( obs_kind == KIND_SPECIFIC_HUMIDITY  .and. vert_is_surface(location)) then 
    ! convert water vapor mixing ratio to specific humidity
    my_type = TYPE_Q2


else if ( obs_kind == KIND_U_WIND_COMPONENT ) then 
    my_type = TYPE_U
else if ( obs_kind == KIND_V_WIND_COMPONENT ) then 
    my_type = TYPE_V
else if ( obs_kind == KIND_SPECIFIC_HUMIDITY ) then 
    ! convert water vapor mixing ratio to specific humidity
    my_type = TYPE_QV


else
    write(errstring,*) "No such obs kind in the state vector: ",obs_kind
    call error_handler(E_ERR,'model_interpolate', errstring, &
       source, revision, revdate)
endif


! find the WRF type that corresponds to the obs_kind
! The new generic obs kind supercedes this block.  
!select case(obs_kind)
!  case(KIND_U)
!    my_type = TYPE_U
!  case(KIND_V)
!    my_type = TYPE_V
!  case(KIND_T)
!    my_type = TYPE_T
!  case(KIND_QV)
!    my_type = TYPE_QV
!  case(KIND_U10)
!    my_type = TYPE_U10
!  case(KIND_V10)
!    my_type = TYPE_V10
!  case(KIND_T2)
!    my_type = TYPE_T2
!  case(KIND_Q2)
!    my_type = TYPE_Q2
!  case default
!    write(errstring,*) "No such obs kind in the state vector: ",obs_kind
!    call error_handler(E_ERR,'model_interpolate', errstring, &
!       source, revision, revdate)
!end select

! If it is NOT a surface variable .... do something
! This could be replaced by the 'vert_is_surface' function.

!if(obs_kind /= KIND_PS .and. obs_kind /= KIND_U10 .and. obs_kind /= KIND_V10 &
!     .and. obs_kind /= KIND_T2 .and. obs_kind /= KIND_Q2) then

if( .not. vert_is_surface(location) ) then

  if ( my_type == TYPE_U .or. my_type == TYPE_V ) then
     do k2 = 1,2

       ugrid = x(get_wrf_index(k+k2-1,TYPE_U))
       vgrid = x(get_wrf_index(k+k2-1,TYPE_V))

       call gridwind_to_truewind(cent_lon_deg,my_projection,ugrid,vgrid,utrue,vtrue)

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

 call gridwind_to_truewind(cent_lon_deg,my_projection,ugrid,vgrid,utrue,vtrue)

 if ( my_type == TYPE_U10 ) then
   obs_val = utrue
 else
   obs_val = vtrue
 endif

elseif ( my_type == TYPE_T2 .or. my_type == TYPE_Q2 ) then ! screen T, Q

 k1 = get_wrf_index(k,my_type)
 obs_val = x(k1)

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
! Does any shutdown and clean-up needed for model. Nothing for L04 for now.


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
real(r8)              :: obs_location
real(r8)              :: state_location
real(r8)              :: dist_tmp


which_vert = query_location(o_loc,'which_vert')
if ( which_vert > 1 ) then
     call error_handler(E_ERR, 'model_get_close_states; column', &
         'which_vert is invalid', source, revision, revdate)
endif

obs_location = get_location(o_loc)

inum = 0
do i = 1, wrf_meta%model_size
  state_location = get_location(state_loc(i))
  dist_tmp = abs(obs_location-state_location)
  if ( dist_tmp <= radius ) then
    inum = inum + 1
    indices(inum) = i
    dist(inum) = dist_tmp
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
integer :: ZVarVarID, SLVarVarID
integer :: ZVarDimID, SLVarDimID
integer :: ZVarID(10), &
           SLVarID(2), ScrVarID(5), SkinVarID(4), SclrVarID(1)
! U, V, Z, T, RHO, QV, TKE, P, U_G, V_G

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1,str2

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
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "PBL_type", pbltype ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "sfc_type", sfctype ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "surface_physics", surface_physics ))
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


!--------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!--------------------------------------------------------------------
if ( output_state_vector ) then
  call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=wrf_meta%model_size, dimid = StateVarDimID))
else
  call check(nf90_def_dim(ncid=ncFileID, name="z_level", &
                        len=nz, dimid = ZVarDimID))
  call check(nf90_def_dim(ncid=ncFileID, name="sl_level", &
                        len=num_soil_layers, dimid = SLVarDimID))
endif


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

  call check(nf90_def_var(ncid=ncFileID, name="PBLH", xtype=nf90_double, &
           dimids = (/ MemberDimID, TimeDimID /), varid=SclrVarID(1)))
  call check(nf90_put_att(ncFileID, SclrVarID(1), "long_name", "PBL height"))
  call check(nf90_put_att(ncFileID, SclrVarID(1), "units", "m"))
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
           ScrVarID(5), SLVarID(2), SkinVarID(4), SclrVarID(1)

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

  call check(NF90_inq_varid(ncFileID, "PBLH", SclrVarID(1)) )
  call check(NF90_put_var(ncFileID, SclrVarID(1), pblh(1,1),  &
               start=(/ copyindex, timeindex /)))

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
           ScrVarID(5),SLVarID(2),SkinVarID(4),SclrVarID(1)

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
   call check(nf90_open(trim(indir)//'/'//init_gen_file,NF90_NOWRITE,ncID))
 
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
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_atts',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check
END SUBROUTINE get_projection

!----------------------------------------------------------------------
!************* Miscellaneous functions/subroutines to do specific jobs
!----------------------------------------------------------------------

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

!===================================================================
! End of model_mod
!===================================================================

end module model_mod

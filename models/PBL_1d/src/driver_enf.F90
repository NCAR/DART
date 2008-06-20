! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

PROGRAM driver_enf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! 1d WRF incorporates original WRF PBL routines for idealized runs.
! It can be initialized with idealized profiles of wind, temperature,
! and mixing ratio or with profiles derived from 3d WRF runs.
! All input is driven by parameters specified in wrf1d_namelist.input
! file. This code is compatible with DART software used for Data
! Assimilatiuon studies.

! The code currently has a problem when compiled with ifc
! as it crashes when calculating MINVAL of an array for unknown reason.
! It runs OK when compiled with pgf90


! The code is available as is but had been checked quite thoroughlly.
! It is of course free to use and modify but I would appreciate
! if you let me know your interest, acknowledge the authorship,
! and, even better if you contribute
! to the development and make the modifications available.

! File README gives more details on input parameters in the namelist.
! File PBLcolumn.pdf gives some general information on the model
! and its coupling with DART.

! 1d WRF was originally developed by Mariusz Pagowski 
! NOAA/GSD/CIRA, Boulder, CO, Mariusz.Pagowski@noaa.gov
! Josh Hacker,NCAR/RAL, hacker@ucar.edu
! contributed reading routines to real data and 3dWRF output and 
! interface to DART and made the code more flexible

! Based on WRFV2.1.2

  USE netcdf
  USE module_wrf
  USE types_mod
  USE time_manager_mod, only : time_type, set_time, get_time, &
                             increment_time, print_time, set_date, &
                             set_calendar_type, GREGORIAN, get_date, &
                             operator(==), operator(<=), &
                             operator(-), operator(+), julian_day 
  USE  utilities_mod, only : file_exist, open_file, close_file, &
                             find_namelist_in_file, check_namelist_read, &
                             register_module, error_handler, E_ERR, E_MSG, &
                             logfileunit

  IMPLICIT NONE

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! These are the hooks for parameter estimation - they are used in DART

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
                      TYPE_U_UPSTREAM_X = 111, TYPE_U_UPSTREAM_Y = 112, &
                      TYPE_V_UPSTREAM_X = 113, TYPE_V_UPSTREAM_Y = 114, &
                      TYPE_TAU_U = 115,        TYPE_TAU_V = 116

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
integer, parameter  ::   number_of_2d_advection =        10
                         ! T_UPSTREAM_X, T_UPSTREAM_Y
                         ! QV_UPSTREAM_X, QV_UPSTREAM_Y
                         ! U_UPSTREAM_X, U_UPSTREAM_Y
                         ! V_UPSTREAM_X, V_UPSTREAM_Y
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

type(state_vector_meta_data)    :: wrf_meta

real(r8), dimension(:), allocatable :: wrf_state

  CHARACTER(len=120) :: &
       &namelistfile='wrf1d_namelist.input',&
       &outlogfile='wrf1d_log.out'
  INTEGER  :: unit_nml=151,unit_log=152
  LOGICAL :: is_it_there = .FALSE.

  INTEGER :: wrf_rnd_seed,itime ! equivalent to itimestep in module_wrf.F

 ! some stuff for netCDF IO
  TYPE netcdf_file_type
     integer :: ncid                       ! the "unit" -- sorta
     character(len=80)        :: fname     ! filename ...
  end TYPE netcdf_file_type

  type(netcdf_file_type)       :: ncFileID

  integer ::   MemberDimID,   MemberVarID     ! for each "copy" or ensemble member
  integer :: nc_time_index
  
  LOGICAL :: allocate_wrf = .TRUE.
  INTEGER :: seconds, days
  INTEGER                  :: sim_seconds,sim_days
  INTEGER :: iyr,imo,idy,ihr,imm,iss,seconds_in_day
  REAL(digits12) :: realtime
  TYPE(time_type) :: initialization_time, time_step, &
                     sim_time, current_time
  integer, parameter :: calendar_type = GREGORIAN

  INQUIRE ( FILE = namelistfile , EXIST = is_it_there )

  IF ( is_it_there ) THEN
     OPEN ( FILE   = namelistfile     , UNIT   =  unit_nml        ,&
          & STATUS = 'OLD'            , FORM   = 'FORMATTED'      ,&
          & ACTION = 'READ'           , ACCESS = 'SEQUENTIAL'     )
  ELSE
     PRINT '(A)','Could not find the namelist: ',namelistfile
     STOP 'No_namelist_found'
  ENDIF

  OPEN ( FILE   = outlogfile        , UNIT   =  unit_log        ,&
       & STATUS = 'REPLACE'         , FORM   = 'FORMATTED'      ,&
       & ACTION = 'WRITE'           , ACCESS = 'SEQUENTIAL'     )

  CALL do_namelist_wrf1d(unit_nml,unit_log)

  CLOSE(unit_log)
  
  CALL static_init_wrf(allocate_wrf)
  CALL state_vector_init()
! state vector allocation
  allocate(wrf_state(wrf_meta%model_size))
  
  wrf_rnd_seed = rnd_seed_val
  allocate_wrf=.FALSE.

  CALL init_wrf(wrf_rnd_seed)
  
! initialize some timing stuff
  time_step = set_time(int(dt), 0)
  call set_calendar_type(calendar_type)
  call init_time(initialization_time)
  current_time = initialization_time

! initialize some of the netCDF file info
  ncFileID%fname = "SCM_output.nc"

! open output netCDF file and define a few things
  call check(nf90_create(path = trim(ncFileID%fname), cmode = nf90_share, ncid = ncFileID%ncid))

  call check(nf90_enddef(ncFileID%ncid))

  call check(nc_write_model_atts(ncFileID%ncid))

  call check(nf90_sync(ncFileID%ncid))

! main time loop in model
  nc_time_index = 1
  current_time = initialization_time
  DO itime=1,ntime+1

     call print_time(current_time,"Taking time step at: ")
! output before taking a time step - currently outputting every time step
     call get_time(current_time,seconds,days)
     realtime = days + seconds/86400.0_digits12
     call check(nc_write_model_vars(ncFileID%ncid,wrf_state,1,nc_time_index))
     nc_time_index = nc_time_index + 1

     IF (init_f) THEN
        sim_time = current_time - initialization_time
        call get_time(sim_time,sim_seconds,sim_days)
        call get_date(current_time,iyr,imo,idy,ihr,imm,iss)
        current_time = increment_time(current_time,int(dt),0)

     ELSE
        sim_days=0
        sim_seconds=(itime-1)*dt
     ENDIF
     IF ( forecast_length == 0 ) stop '0 time steps'
     IF ( init_f ) THEN
       seconds_in_day = iss + 60*imm + 3600 * ihr
       CALL wrf(sim_seconds,sim_days,julian_day(iyr,imo,idy),seconds_in_day)
     ELSE
       iyr = 1970
       imo = 1
       idy = 0
       seconds_in_day = 0
       CALL wrf(sim_seconds,sim_days,julian_day(iyr,imo,idy),seconds_in_day)
!       CALL wrf(dart_seconds,dart_days)
     ENDIF
     if ( forecast_length == 0 ) stop '0 time steps'
     CALL wrf(sim_seconds,sim_days,julian_day(iyr,imo,idy),seconds_in_day)
!     CALL output_wrf_profiles()
    current_time = current_time + time_step
  ENDDO

  call check(nf90_close(ncFileID%ncid))
  
!*******************************************************************
  CONTAINS

  subroutine state_vector_init()

  implicit none


  integer :: dart_index, var_cnt, ivar, vcnt

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

allocate(wrf_meta%var_type(wrf_meta%total_number_of_vars))
allocate(wrf_meta%var_size(wrf_meta%total_number_of_vars))
allocate(wrf_meta%var_index(wrf_meta%total_number_of_vars))

! fill locations for state variables
dart_index = 1
var_cnt    = 1

wrf_meta%var_type(var_cnt) = TYPE_U
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_V
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_GZ
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_T
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TH
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_QV
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TKE
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_P
wrf_meta%var_size(var_cnt) = nz
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_U10
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_V10
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_T2
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_Q2
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_PSFC
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TSLB
wrf_meta%var_size(var_cnt) = num_soil_layers
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

wrf_meta%var_type(var_cnt) = TYPE_TSK
wrf_meta%var_size(var_cnt) = 1
wrf_meta%var_index(var_cnt) = dart_index
dart_index = dart_index + wrf_meta%var_size(var_cnt)
var_cnt = var_cnt + 1

! here are the parameters appended to the state variable
do i = 1, wrf_meta%number_state_params
  wrf_meta%var_type(var_cnt) = wrf_meta%est_param_types(i)
  wrf_meta%var_size(var_cnt) = 1
  wrf_meta%var_index(var_cnt) = dart_index
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
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! 2d F vars (nz,nsplinetimes)
do vcnt = 1, wrf_meta%number_of_2d_f
  wrf_meta%var_type(var_cnt) = 100 + wrf_meta%number_of_1d_f + vcnt - 1 !TYPE_UF ... TYPE_P8WF
  wrf_meta%var_size(var_cnt) = nsplinetimes*nz
  wrf_meta%var_index(var_cnt) = dart_index
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! 2d advection-related vars (nz,nsplinetimes_advection)
do vcnt = 1, wrf_meta%number_of_2d_advection
  wrf_meta%var_type(var_cnt) = 100 + wrf_meta%number_of_1d_f &
          + wrf_meta%number_of_2d_f + vcnt - 1 !TYPE_TH_UPSTREAM_X...TYPE_TAU_V
  wrf_meta%var_size(var_cnt) = nsplinetimes_advection*nz
  wrf_meta%var_index(var_cnt) = dart_index
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! noassim profile vars (nz)
do vcnt = 1, wrf_meta%number_noassim_profiles
  wrf_meta%var_type(var_cnt) = 500 + vcnt - 1 !RHO...EL
  wrf_meta%var_size(var_cnt) = nz
  wrf_meta%var_index(var_cnt) = dart_index
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! noassim soil vars (num_soil_layers)
do vcnt = 1, wrf_meta%number_noassim_soil_vars
  wrf_meta%var_type(var_cnt) = 200 + vcnt - 1 !TYPE_SMOIS ... TYPE_SH2O
  wrf_meta%var_size(var_cnt) = num_soil_layers
  wrf_meta%var_index(var_cnt) = dart_index
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! scalars (might want these in the assimilation vector)
do vcnt = 1, wrf_meta%number_noassim_scalars
  wrf_meta%var_type(var_cnt) = 300 + vcnt - 1 !TYPE_UZ0 ... TYPE_GAMMA
  wrf_meta%var_size(var_cnt) = 1
  wrf_meta%var_index(var_cnt) = dart_index
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

! parameters that depend on something (strange ones)
do vcnt = 1, wrf_meta%number_dependent_params
  wrf_meta%var_type(var_cnt) = 400 + vcnt - 1 !TYPE_VEGFRA...?
  wrf_meta%var_size(var_cnt) = 1
  wrf_meta%var_index(var_cnt) = dart_index
  dart_index = dart_index + wrf_meta%var_size(var_cnt)
  var_cnt = var_cnt + 1
enddo

END SUBROUTINE state_vector_init
!*******************************************************************


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
        case (TYPE_U_UPSTREAM_X)
          x(k) = u_upstream_x_f(iz,ispline) 
        case (TYPE_U_UPSTREAM_Y)
          x(k) = u_upstream_y_f(iz,ispline) 
        case (TYPE_V_UPSTREAM_X)
          x(k) = v_upstream_x_f(iz,ispline) 
        case (TYPE_V_UPSTREAM_Y)
          x(k) = v_upstream_y_f(iz,ispline) 
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

!------------------------------------------------------------------

subroutine vector_to_wrf(x)
!---------------------------------------------------------
! rolls state vector into wrf state 
!---------------------------------------------------------
  implicit none

  real(r8), intent(inout)  :: x(:)

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
        case (TYPE_U_UPSTREAM_X)
          u_upstream_x_f(iz,ispline) = x(k)
        case (TYPE_U_UPSTREAM_Y)
          u_upstream_y_f(iz,ispline) = x(k)
        case (TYPE_V_UPSTREAM_X)
          v_upstream_x_f(iz,ispline) = x(k)
        case (TYPE_V_UPSTREAM_Y)
          v_upstream_y_f(iz,ispline) = x(k)
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
integer :: StateVarID, MemberDimID

!------------------------------------------
! same for physical space
!------------------------------------------
integer :: ZVarVarID, SLVarVarID, PVarVarID
integer :: ZVarDimID, SLVarDimID, PVarDimID
integer :: ZVarID(10), &
           SLVarID(2), ScrVarID(5), SkinVarID(5), SclrVarID(4)
! U, V, Z, T, RHO, QV, TKE, P, U_G, V_G
integer :: PVarID(wrf_meta%number_state_params) 
integer ::     TimeDimID,     TimeVarID

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
ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_sync(ncFileID)) ! Ensure netCDF file is current
call check(nf90_Redef(ncFileID))

call check(nf90_def_dim(ncid=ncFileID, &
             name="copy",           len=1, dimid = MemberDimID))
call check(nf90_def_var(ncid=ncFileID, name="copy", xtype=nf90_int, dimids=MemberDimID, &
             varid=MemberVarID)) 
call check(nf90_put_att(ncFileID, MemberVarID, "long_name", "ensemble member or copy"))
call check(nf90_put_att(ncFileID, MemberVarID, "units", "nondimensional") )
call check(nf90_put_att(ncFileID, MemberVarID, "valid_range", (/ 1, 1 /)))

call check(nf90_def_dim(ncid=ncFileID,  &
        name="time",           len = nf90_unlimited,         dimid = TimeDimID))
call check(nf90_def_var(ncFileID, name="time", xtype=nf90_double, dimids=TimeDimID, &
       varid =TimeVarID) ) 
  


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

! Leave define mode so we can fill
call check(nf90_enddef(ncfileID))

call check(nf90_put_var(ncFileID, MemberVarID,   (/ (i,i=1,1) /) ))

!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------
call check(nf90_sync(ncFileID))

write (*,*)'Model attributes written, netCDF file synched ...'

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
integer :: PVarID(wrf_meta%number_state_params), TimeVarID

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

! right now no need to get state from vector - alrady have it
!  call vector_to_wrf(statevec)

  call check(NF90_inq_varid(ncFileID, "time", TimeVarID) )
  call check(NF90_put_var(ncFileID, TimeVarID, realtime, start = (/timeindex/)))
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

! write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
! write (*,*)'netCDF file is synched ...'

end function nc_write_model_vars

!----------------------------------------------------------------------

function nc_read_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!------------------------------------------------------------------
! Reads the state vector 
! JPH
!

use typeSizes

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

! write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
! write (*,*)'netCDF file is synched ...'

end function nc_read_model_vars

!----------------------------------------------------------------------

   subroutine check(istatus)
      integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_atts',&
      trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

END PROGRAM driver_enf

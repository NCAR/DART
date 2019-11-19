
! FESOM 2 (Finite-volumE Sea ice-Ocean Model)
! multi-resolution ocean general circulation model
! FESOM/fesom2 is licensed under the GNU General Public License v2.0
! Copyright (C) 2018  FESOM team
!
! This module was constructed from pieces of the FESOM V1.4 modules and
! extended to work with DART.

module g_config
  implicit none
  save

  ! *** Modelname ***
  character(5)                  :: runid                ! a model/setup name
  character(5)                  :: ensid='ENS01'        ! a model/setup name
  integer                       :: ensmem=1

  namelist /modelname/ runid, ensid, ensmem

  integer                       :: dart_days, dart_secs

  namelist /dart_nml/ dart_days, dart_secs
  ! *** time step ***
  integer                       :: step_per_day=12           !number of steps per day
  integer                       :: run_length=1                !run length
  character                     :: run_length_unit='y'          !unit: y, d, s
  character(4)                  :: runyear='2008'

  namelist /timestep/ step_per_day, run_length, run_length_unit, runyear

  ! *** time series *
  integer                       :: iniday=1
  integer                       :: endday=1

  namelist /timeseries/ IniDay, EndDay

  ! *** Paths for all in and out ***
  character(100)                :: MeshPath='./mesh/'
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
  character*4                   :: restartflag='last'               !restart from which saved record,'#','last'
  integer                       :: output_length=1                   !valid for d,h,s
  character                     :: output_length_unit='m'           !output period: y, m, d, h, s 
  integer                       :: logfile_outfreq=1                 !in logfile info. output frequency, # steps

  namelist /inout/ restartflag, output_length, output_length_unit, logfile_outfreq

  ! *** mesh ***
  integer                       :: grid_type=1              ! z-level, 2 sigma, 3 sigma + z-level

  namelist /mesh_def/ grid_type

  ! *** model geometry
  logical                  :: cartesian=.false.
  logical                  :: fplane=.false.
  logical                  :: betaplane=.false.
  real(kind=8)             :: f_fplane=-1.4e-4        ![1/s]
  real(kind=8)             :: beta_betaplane=2.0e-11  ![1/s/m]
  real(kind=8)             :: domain_length=360.    ![degree]
  !
  logical                  :: rotated_grid=.true.    !option only valid for coupled model case now
  real(kind=8)             :: alphaEuler=-30.        ![degree] Euler angles, convention:
  real(kind=8)             :: betaEuler=-90.         ![degree] first around z, then around new x,
  real(kind=8)             :: gammaEuler=-90.        ![degree] then around new z.

  namelist /geometry/  cartesian, fplane, betaplane, f_fplane, beta_betaplane, &
       domain_length, rotated_grid, alphaEuler, betaEuler, gammaEuler

  ! *** fleap_year ***
  logical                       :: include_fleapyear=.false.
  
  namelist /calendar/ include_fleapyear
  
   ! *** machine ***
  integer                       :: system=1                     ! XD1 2(byte), HLRN 1(word)
  
  namelist /machine/ system

  integer                       :: tool=1              
  integer                       :: level_number
  character*2                   :: flux_section
  character*100                 :: thalweg_directory

  namelist /postproc/ tool, level_number, flux_section, thalweg_directory

  ! *** others ***
  real(kind=8)                  :: dt, dt_inv
  integer                       :: istep, nsteps
  integer                       :: save_count
  integer                       :: day2ext
  logical                       :: r_restart

end module g_config

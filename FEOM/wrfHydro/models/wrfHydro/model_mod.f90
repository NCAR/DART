! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! Interface for wrfHydro
! This is modeled on an earlier noah interface, so the term noah maybe used when 
! wrfHydro is more appropriate.
! We are attempting to accomodate use of either Noah or NoahMP as LSMs.
! We are also attempting to accomodate runing with/without the hydro sstate_component, so assim would
! just be for the LSM without the hydro component.
! See input.nml and it's model_nml to specfiy LSM and if hydro compnent is active.
! The state variables of interest can also be specified here. (Debatable if desired to show
! state variables for inactive components (potential confusion), perhaps can be filtered out for log?)

module model_mod

use        types_mod, only : r4, r8, MISSING_R8, obstypelength
use time_manager_mod, only : time_type, set_time, set_date, get_time,          &
     print_time, print_date, set_calendar_type,        &
     operator(*),  operator(+), operator(-),           &
     operator(>),  operator(<), operator(/),           &
     operator(/=), operator(<=)

! TJH FIXME ... to accomodate the fact that watersheds are 'close' but not physically
! 'close' - we need to have a function that returns the watershed ID for any arbitrary
! location. Our local get_close_obs() can then get everything within the (large) 
! localization radius and for all locations in different watersheds - make them very 
! far away. get_close() is used not only state vector elements, but also for 
! observations close to a location.

use     location_mod, only : location_type, get_dist, query_location,          &
     get_close_maxdist_init, get_close_type,           &
     set_location, get_location, horiz_dist_only,      &
     vert_is_surface,  VERTISSURFACE,                  &
     vert_is_height,   VERTISHEIGHT,                   &
     vert_is_level,    VERTISLEVEL,                    &
     get_close_obs_init, get_close_obs,                &
     set_location_missing, write_location

use    utilities_mod, only : register_module, error_handler, nc_check,         &
     E_ERR, E_MSG, logfileunit, get_unit,              &
     nmlfileunit, do_output, do_nml_file, do_nml_term, &
     find_namelist_in_file, check_namelist_read,       &
     file_exist, find_textfile_dims, file_to_text

use    obs_kind_mod, only :                   &
     KIND_SOIL_MOISTURE,                      &
     KIND_SOIL_LIQUID_WATER,                  &
     KIND_SURFACE_HEAD,                       &
     KIND_STREAM_FLOW,                        &
     KIND_STREAM_HEIGHT,                      &
     KIND_DEEP_GROUNDWATER_LEVEL,             &
     KIND_SOIL_TEMPERATURE,                   &
     KIND_GROUND_SURF_TEMPERATURE,            &
     KIND_CANOPY_TEMPERATURE,                 &
     KIND_LEAF_AREA_INDEX,                    &
     KIND_SNOW_THICKNESS,                     &
     KIND_SNOW_WATER,                         &
     KIND_SNOWCOVER_FRAC,                     &
     KIND_2D_PARAMETER,        &     !! eg precip multiplier, soil texture parameter
     KIND_GEOPOTENTIAL_HEIGHT, &     !! maybe be used for model_interpolate
     paramname_length,         &
     get_raw_obs_kind_index

use mpi_utilities_mod, only: my_task_id
use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use typesizes
use netcdf

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
public :: &
     get_model_size,         &
     adv_1step,              &
     get_state_meta_data,    &
     model_interpolate,      &
     get_model_time_step,    &
     end_model,              &
     static_init_model,      &
     init_time,              &
     init_conditions,        &
     nc_write_model_atts,    &
     nc_write_model_vars,    &
     pert_model_state,       &
     get_close_maxdist_init, &
     get_close_obs_init,     &
     get_close_obs,          &
     ens_mean_for_model

! not required by DART but for larger models can be useful for
! utility programs that are tightly tied to the other parts of
! the model_mod code.

public :: &
     model_to_dart_vector,            & 
     dart_vector_to_model_files,      & 
     get_lsm_restart_filename,        &     
     get_hydro_restart_filename,      &   
     get_assimOnly_restart_filename,  &
     get_model_timestepping,          & 
     get_debug_level

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
     "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer :: nfields
! The NSOLDX (number of soil layers) parameter comes from the NOAH source code. We need it
! because we have to read the NOAH namelist for timestep information.
integer, parameter :: NSOLDX = 100  !!jlm typically set to 4, this is an upper bound
integer, parameter :: MAX_STATE_VARIABLES = 40
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 2

!------------------------------------------------------------------
! Things which can/should be in the DART model_nml
! The variables in the noah restart file that are used to create the
! DART state vector are specified in the input.nml:model_nml namelist.
! jlm - are these dummy values? they should be read from file, this is somewhat confusing, should be
! stated that these are dummies.
character(len=128)    :: lsm_model_choice             = 'noahMP'
logical               :: hydro_model_active           = .true.
character(len=256)    :: lsm_netcdf_filename          = 'restart.nc'
character(len=256)    :: hydro_netcdf_filename        = 'restart.hydro.nc'
character(len=256)    :: assimOnly_netcdf_filename    = ''  !! intentionally blank
integer               :: assimilation_period_days     = 0
integer               :: assimilation_period_seconds  = 3600
real(r8)              :: model_perturbation_amplitude = 0.2
logical               :: output_state_vector          = .true.
integer               :: debug    = 0  ! turn up for more and more debug messages

! These are the state variables to be considered in DART.
! jlm - How does one specify if a certain parameter set should be considered as uncertain? skip in model_nml?
character(len=obstypelength), allocatable, dimension(:,:) :: all_state_variables
character(len=obstypelength), dimension(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) :: lsm_state_variables = ' '
character(len=obstypelength), dimension(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) :: hydro_state_variables = ' '
character(len=obstypelength), dimension(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) :: assimOnly_state_variables = ' '


namelist /model_nml/ lsm_model_choice, hydro_model_active,         &
     lsm_netcdf_filename, hydro_netcdf_filename, assimOnly_netcdf_filename,  &
     assimilation_period_days, assimilation_period_seconds,   &
     model_perturbation_amplitude, output_state_vector,       &
     debug, lsm_state_variables, hydro_state_variables,       &
     assimOnly_state_variables

logical               :: assimOnly_active           = .false.

!------------------------------------------------------------------
! From the models' namelists we get everything needed to recreate
! specify these (for restarts after DART filters).
! For both noah and noahMP the namelist is called namelist.hrldas
character(len=256) :: lsm_namelist_filename   = 'namelist.hrldas' ! mandate
character(len=256) :: hydro_namelist_filename = 'hydro.namelist'  ! mandate

! jlm
! Small conundrum of how to handle the noah and noahMp namelists.
! They are both called NOAHLSM_OFFLINE but they have different variables.
! Cannot rename them in namelist.hrldas b/c that would require rewriting HRLDAS
! or noah/noahMP.
! Seems like the best way is to specify agnostic defaults, then attempt to
! read the union of their variables from the file. This should not change variables not
! found in the file (which would remain to their default, agnostic values).

!Variables in both noah and noahMP (spaces reflect the general grouping in our nml files)
character(len=256) :: hrldas_constants_file = " "
character(len=256) :: indir = "."

character(len=256) :: outdir = "."

integer            :: start_year, start_month, start_day
integer            :: start_hour, start_min

character(len=256) :: restart_filename_requested = " "

integer            :: kday  = 0
integer            :: khour = 0

integer            :: forcing_timestep = -999
integer            :: noah_timestep = -999
integer            :: output_timestep  = -999

integer            :: restart_frequency_hours = -999
integer            :: split_output_count = 1

integer            :: nsoil ! number of soil layers in use.

real(r8)           :: zlvl

integer            :: iz0tlnd = 0
integer            :: sfcdif_option = 0
logical            :: update_snow_from_forcing = .true.

!! wrfHydro specific ??
integer            :: FORC_TYP = 3
character(len=256) :: GEO_STATIC_FLNM = "DOMAIN/geo_em.d03.nc"
integer            :: HRLDAS_ini_typ = 0
integer            :: SNOW_assim = 0

!! only Noah
real(r8), dimension(NSOLDX) :: zsoil
integer            :: subwindow_xstart = 1
integer            :: subwindow_ystart = 1
integer            :: subwindow_xend = 0
integer            :: subwindow_yend = 0
real(r8)           :: zlvl_wind

!! only NoahMP
character(len=256) :: MMF_RUNOFF_FILE = ""
integer            :: DYNAMIC_VEG_OPTION                = 4
integer            :: CANOPY_STOMATAL_RESISTANCE_OPTION = 1
integer            :: BTR_OPTION                        = 4
integer            :: RUNOFF_OPTION                     = 3
integer            :: SURFACE_DRAG_OPTION               = 1
integer            :: FROZEN_SOIL_OPTION                = 1
integer            :: SUPERCOOLED_WATER_OPTION          = 1
integer            :: RADIATIVE_TRANSFER_OPTION         = 3
integer            :: SNOW_ALBEDO_OPTION                = 2
integer            :: PCP_PARTITION_OPTION              = 1
integer            :: TBOT_OPTION                       = 1
integer            :: TEMP_TIME_SCHEME_OPTION           = 1
real(r8), dimension(NSOLDX) :: soil_thick_input         = MISSING_R8

!! Not in either of our noahlsm_offline nmls but in the earlier noah-DART model_mod
!character(len=256) :: external_fpar_filename_template = " "
!character(len=256) :: external_lai_filename_template = " "

namelist /NOAHLSM_OFFLINE/ hrldas_constants_file, indir, outdir, &
     start_year, start_month, start_day, start_hour, start_min, &
     restart_filename_requested, kday, khour, forcing_timestep, &
     noah_timestep, output_timestep, restart_frequency_hours, split_output_count, &
     nsoil, zlvl, iz0tlnd, sfcdif_option, update_snow_from_forcing, &
     FORC_TYP, GEO_STATIC_FLNM, HRLDAS_ini_typ, SNOW_assim, zsoil, &
     subwindow_xstart, subwindow_ystart, subwindow_xend, subwindow_yend, zlvl_wind, &
     MMF_RUNOFF_FILE, DYNAMIC_VEG_OPTION, CANOPY_STOMATAL_RESISTANCE_OPTION, BTR_OPTION, &
     RUNOFF_OPTION, SURFACE_DRAG_OPTION, FROZEN_SOIL_OPTION, SUPERCOOLED_WATER_OPTION, &
     RADIATIVE_TRANSFER_OPTION, SNOW_ALBEDO_OPTION, PCP_PARTITION_OPTION, TBOT_OPTION, &
     TEMP_TIME_SCHEME_OPTION, soil_thick_input

!&URBAN_OFFLINE
! This is in namelist.hrldas.
integer  :: UCMCALL = 0
real(r8) :: ZLVL_URBAN = 15.0
namelist /URBAN_OFFLINE/ UCMCALL,  ZLVL_URBAN

!! &HYDRO_nlist
!! The noah and noahMP models have some same/repeated variables in their respective namelists.
!! I note repeated varaibles and any related issues here.

integer            :: sys_cpl = 1  !! this is hrldas, should be enforced
!! character(len=256) :: GEO_STATIC_FLNM = "" !! repeated in the two namelists, but equal
character(len=256) :: GEO_FINEGRID_FLNM = ""
character(len=256) :: RESTART_FILE  = ''
integer            :: IGRID = 3
integer            :: rst_dt = 1440
integer            :: out_dt = 1440
logical            :: HISTORY_OUTPUT = .true.
!! integer            :: SPLIT_OUTPUT_COUNT = 1  !! repeated but equal
integer            :: rst_typ = 1
integer            :: RSTRT_SWC = 0
integer            :: HIRES_OUT = 2
integer            :: order_to_write = 1
integer            :: TERADJ_SOLAR = 0
!! integer            :: NSOIL=4  !! repeated but equal
real(r8), dimension(NSOLDX) :: zsoil8  !! this is for the hydro component (bad name)
real(r8)           :: DXRT = -999.0_r8
integer            :: AGGFACTRT = -999
integer            :: DTRT = 2
integer            :: SUBRTSWCRT = 1
integer            :: OVRTSWCRT = 1
integer            :: rt_option    = 1
integer            :: CHANRTSWCRT = 1
integer            :: channel_option =3
character(len=256) :: route_link_f = ""
integer            :: GWBASESWCRT = 2
integer            :: GW_RESTART = 1
character(len=256) :: gwbasmskfil = "DOMAIN/basn_msk1k_frng_ohd.txt"

namelist /HYDRO_nlist/ sys_cpl, GEO_STATIC_FLNM, GEO_FINEGRID_FLNM, RESTART_FILE, &
     IGRID, rst_dt, out_dt, HISTORY_OUTPUT, SPLIT_OUTPUT_COUNT, rst_typ, RSTRT_SWC, HIRES_OUT, &
     order_to_write, TERADJ_SOLAR, NSOIL, zsoil8, DXRT, AGGFACTRT, DTRT, SUBRTSWCRT, &
     OVRTSWCRT, rt_option, CHANRTSWCRT, channel_option, route_link_f, GWBASESWCRT, &
     GW_RESTART, gwbasmskfil

!------------------------------------------------------------------
! Everything needed to describe a variable and
! its relationship to the DART vector.

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME) :: component  !! lsm vs hydro
   character(len=NF90_MAX_NAME) :: grid  
   character(len=NF90_MAX_NAME) :: MemoryOrder
   character(len=NF90_MAX_NAME) :: gridMemOrd  ! combo of the previous 2
   character(len=NF90_MAX_NAME) :: description
   character(len=NF90_MAX_NAME) :: stagger
   character(len=obstypelength), dimension(NF90_MAX_VAR_DIMS) :: dimnames
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer  :: numdims
   integer  :: maxlevels
   integer  :: xtype       ! the ncdf variable type
   integer  :: varsize     ! prod(dimlens(1:numdims))
   integer  :: index1      ! location in dart state vector of first occurrence
   integer  :: indexN      ! location in dart state vector of last  occurrence
   integer  :: dart_kind
   integer  :: rangeRestricted
   real(r8) :: maxvalue
   real(r8) :: minvalue
   character(len=paramname_length) :: kind_string
end type progvartype

!! jlm fixme? can't this be allocated to the appropriate size after reading input.nml?
type(progvartype), dimension(MAX_STATE_VARIABLES) :: progvar

! State/'progvar' location information is stored here!
type(location_type),allocatable, dimension(:) :: state_loc

!--------------------------------------------------------------------------
! model parameters
type(time_type)     :: time_step

!------------------------------------------------------------------------------
! These are the metadata arrays that are the same size as the state vector.
real(r8), allocatable, dimension(:) :: ens_mean ! may be needed for forward ops

!------------------------------------------------------------------
! module storage
!------------------------------------------------------------------

integer            :: model_size       ! the state vector length
type(time_type)    :: model_time       ! valid time of the model state
type(time_type)    :: model_time_step  ! smallest time to adv model
character(len=512) :: string1, string2, string3
logical, save      :: module_initialized = .false.
character(len=32)  :: calendar = 'Gregorian'

real(r8), allocatable, dimension(:,:) :: xlong, xlat, hlong, hlat
real(r8), allocatable, dimension(:) :: linkLat, linkLong, channelIndsX, channelIndsY
real(r8), allocatable, dimension(:) :: basnMask, basnLon, basnLat
integer :: south_north, west_east, n_hlong, n_hlat, n_link, n_basn
integer, dimension(2)               :: fine2dShape, coarse2dShape
integer, dimension(3)               :: fine3dShape, coarse3dShape
!! Following are global because they are used in multiple subroutines.
real(R8), dimension(:,:,:), allocatable :: smc, sice, sh2oMaxRt, sh2oWltRt
real(R8), dimension(:,:),   allocatable ::              smcMax1, smcWlt1
!! Used to hold the fine res variables to be adjusted. Not allocated if not used.
real(R8), dimension(:,:,:), allocatable :: sh2oDisag
real(R8), dimension(:,:),   allocatable :: sfcHeadDisag
real(r8) :: fineGridArea, coarseGridArea
integer  :: hydroSmcPresent, hydroSfcHeadPresent

interface vector_to_prog_var
   module procedure vector_to_1d_prog_var
   module procedure vector_to_2d_prog_var
   module procedure vector_to_3d_prog_var
end interface vector_to_prog_var

!==================================================================
contains
!==================================================================


!===============================================================================
subroutine static_init_model()
! one time initialization of the model

! Local variables - all the important ones have module scope

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
character(len=NF90_MAX_NAME)          :: dimname
character(len=paramname_length)       :: kind_string
character(len=NF90_MAX_NAME)          :: component

! these are not needed globally, state_loc should handle that.
real(r8), allocatable, dimension(:) :: state_lon, state_lat, state_level
real(r8), allocatable, dimension(:) :: dumGridLon, dumGridLat, dumGridLevel
real(r8), allocatable, dimension(:) :: zsoilComp  
integer  :: VarID, dimlen, varsize
integer  :: iunit, io, ivar, iunit_lsm, iunit_hydro, iunit_assimOnly, igrid, iState
integer  :: ilat, ilon, ilev, myindex,  i, index1
integer  :: nLayers, n_lsm_fields, n_hydro_layers, n_hydro_fields, n_assimOnly_fields
integer  :: dumNLon, dumNLat, dumSize, wp, dumNumDims, lsmSmcPresent
integer, allocatable, dimension(:)  :: whichVars, keepLsmVars0, keepLsmVars
integer  :: whVar1

character(len=32), allocatable, dimension(:) :: uniqueGridMemOrd
character(len=32) :: dumGridMemOrd
integer  :: nUniqueGrids

if ( module_initialized ) return ! only need to do this once.

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the DART namelist
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the DART namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Check to make sure the required LSM input files exist
if ( .not. file_exist(lsm_netcdf_filename) ) then
   write(string1,*) 'cannot open LSM restart file [', trim(lsm_netcdf_filename),'] for reading.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif
if ( .not. file_exist(lsm_namelist_filename) ) then
   write(string1,*) 'cannot open LSM namelist file [', trim(lsm_namelist_filename),'] for reading.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif

! Check to make sure the required HYDRO input files exist
if (hydro_model_active) then
   if ( .not. file_exist(hydro_netcdf_filename) ) then
      write(string1,*) 'cannot open HYDRO restart file [', trim(hydro_netcdf_filename),'] for reading.'
      call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
   endif
   if ( .not. file_exist(hydro_namelist_filename) ) then
      write(string1,*) 'cannot open HYDRO namelist file [', trim(hydro_namelist_filename),'] for reading.'
      call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
   endif
endif

! Determine if the "assimOnly" or parameter component of the assimilation is active based on 
! two namelist variables. 
assimOnly_active = ( (trim(assimOnly_state_variables(1,1)) .ne. '') .and. &
     (trim(assimOnly_netcdf_filename) .ne. '') )
if (assimOnly_active) then 
   ! check to make sure that the required assimOnly input/restart file exists. 
   if ( .not. file_exist(assimOnly_netcdf_filename) ) then
      write(string1,*) 'cannot open assimOnly restart file ', & 
           trim(assimOnly_netcdf_filename),' for reading.'
      call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
   endif
endif

! Read the LSM namelist
call find_namelist_in_file(lsm_namelist_filename, 'NOAHLSM_OFFLINE', iunit)
read(iunit, nml = NOAHLSM_OFFLINE, iostat = io)
call check_namelist_read(iunit, io, 'NOAHLSM_OFFLINE')

! Read the hydro namelist
if (hydro_model_active) then
   call find_namelist_in_file(hydro_namelist_filename, 'HYDRO_nlist', iunit)
   read(iunit, nml = HYDRO_nlist, iostat = io)
   call check_namelist_read(iunit, io, 'HYDRO_nlist')
endif

! Dxrt is necessary for aggregation. it may not be needed but... if it dosent need to 
! exist there are other potential problems to ponder. 
if (dxrt .lt. 0.) then
   write(string1,*) 'dxrt was not specified in the hydro.namelist. If you are not running an equal area projected grid, please proceede with extreme caution.'
   call error_handler(E_ERR, 'static_init_model', string1, source, revision, revdate)
endif
!! this is in square meters, dxrt in meters
fineGridArea = dxrt**2.
coarseGridArea = (dxrt*AGGFACTRT)**2.

! If the zsoils dont match between the models, throw an error
if (trim(lsm_model_choice) .eq. 'noah') then
   if( .not.  all( zsoil(1:nsoil) == zsoil8(1:nsoil) )) then
      write(string1,*) 'soils layer specifications in the two namelist files are not identical '
      write(string2,*) 'zsoil   from noah  [',trim(lsm_namelist_filename),']'
      write(string3,*) 'zsoil8  from hydro [',trim(hydro_namelist_filename),']'
      call error_handler(E_ERR, 'static_init_model', string1, &
                 source, revision, revdate, text2=string2, text3=string3)
   endif
soil_thick_input = zsoil
do i=nsoil,2,-1 
   soil_thick_input(i)=zsoil(i)-zsoil(i-1)
enddo
endif

if (trim(lsm_model_choice) .eq. 'noahMP') then
   allocate(zsoilComp(nsoil))
   zsoilComp = zsoil8*0-999   ! TJH FIXME i.e. -999 ?
   do i = 1,nsoil
      zsoilComp(i) = -1.0_r8 * sum(soil_thick_input(1:i))
   enddo
   if( .not.  all( zsoil8(1:nsoil) == zsoilComp(1:nsoil) )) then
      write(string1,*) 'soil layer specifications in the two namelist files are not identical '
      write(string2,*) 'zsoil   from noahMP [',trim(lsm_namelist_filename),']'
      write(string3,*) 'zsoil8  from hydro  [',trim(hydro_namelist_filename),']'
      call error_handler(E_ERR, 'static_init_model', string1, &
                 source, revision, revdate, text2=string2, text3=string3)
   endif
   deallocate(zsoilComp)
endif 



! Record the NOAH namelist
if (do_nml_file()) write(nmlfileunit, nml=NOAHLSM_OFFLINE)
if (do_nml_term()) write(     *     , nml=NOAHLSM_OFFLINE)

! Record the hydro namelist
if (hydro_model_active) then
   if (do_nml_file()) write(nmlfileunit, nml=HYDRO_nlist)
   if (do_nml_term()) write(     *     , nml=HYDRO_nlist)
end if

! Check to make sure the hrldasconstants file exists
if ( .not. file_exist(hrldas_constants_file) ) then
   write(string1,*) &
      'cannot open NOAH constants file [',trim(hrldas_constants_file),'] for reading.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif

! Check to make sure the required NOAH namelist items are set:
if ( (kday             < 0    ) .or. &
     (khour            < 0    ) .or. &
     (forcing_timestep /= 3600) .or. &
     (noah_timestep    /= 3600) .or. &
     (output_timestep  /= 3600) .or. &
     (restart_frequency_hours /= 1) ) then
   write(string3,*)'the only configuration supported is for hourly timesteps &
               &(kday, khour, forcing_timestep==3600, noah_timestep=3600, &
               &output_timestep=3600, restart_frequency_hours=1)'
   write(string2,*)'restart_frequency_hours must be equal to the noah_timestep'
   write(string1,*)'unsupported noah namelist settings'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate,&
        text2=string2,text3=string3)
endif

! This gets the LSM geospatial information for the module:
!   south_north, west_east, xlong, xlat
call get_hrldas_constants(hrldas_constants_file)

if (hydro_model_active) then
   !  though all non-soil variables are "surface" it may be advisable to extract 
   ! elevation at this point?  for localization routines?
   ! **** NOTE that all variables from this file (Fulldom) must  ****
   ! ****      be FLIPPED in y to match the noah/wrf model.      ****
   call get_hydro_constants(GEO_FINEGRID_FLNM) !!
   !! global variables defining the aggregation/disaggregation dimension
   !! these are always xyz since they are for dis/agg wrfHydro variables 
   fine2dShape   = (/ n_hlong, n_hlat /)  !! for disaggregating
   fine3dShape   = (/ n_hlong, n_hlat, nsoil /)  !! for disaggregating
   coarse2dShape = (/ west_east, south_north /) !! for reaggregating
   coarse3dShape = (/ west_east, south_north, nsoil /) !! for reaggregating
else
   ! TJH FIXME fine3dShape coarse3dShape are used later ... even when (logically) hydro_model_active may be false.
   ! JLM fixme - that shouldnt happen since there's no "fine" grid for disag, will try to fix... 
   fine3dShape   = (/ n_hlong, n_hlat, nsoil /)  !! for disaggregating
   coarse3dShape = (/ west_east, south_north, nsoil /) !! for reaggregating
endif

! The time_step in terms of a time type must also be initialized.
call set_calendar_type( calendar )
call nc_check(nf90_open(adjustl(lsm_netcdf_filename), NF90_NOWRITE, iunit_lsm), &
     'static_init_model', 'open '//trim(lsm_netcdf_filename))

model_time = get_state_time(iunit_lsm, trim(lsm_netcdf_filename))

! FIXME ... make sure model_time_step is attainable given OUTPUT_TIMESTEP
model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)
if (do_output() .and. (debug > 0)) then
   call print_date(model_time     ,' static_init_model:model date')
   call print_time(model_time     ,' static_init_model:model time')
   call print_time(model_time_step,' static_init_model:model timestep')
endif

! Make sure the number of soil layers is as we expect
! Check that LSM and hydro have the same sub-surface dimensions
call nc_check(nf90_inq_dimid(iunit_lsm, 'soil_layers_stag', dimIDs(1)), &
     'static_init_model','inq_dimid soil_layers_stag '//trim(lsm_netcdf_filename))
call nc_check(nf90_inquire_dimension(iunit_lsm, dimIDs(1), len=nLayers), &
     'static_init_model','inquire_dimension soil_layers_stag '//trim(lsm_netcdf_filename))

if (hydro_model_active) then
   call nc_check(nf90_open(adjustl(hydro_netcdf_filename), NF90_NOWRITE, iunit_hydro), &
        'static_init_model', 'open '//trim(hydro_netcdf_filename))
   call nc_check(nf90_inq_dimid(iunit_hydro, 'depth', dimIDs(1)), &
        'static_init_model','inq_dimid soil_layers_stag '//trim(hydro_netcdf_filename))
   call nc_check(nf90_inquire_dimension(iunit_hydro, dimIDs(1), len=n_hydro_layers), &
        'static_init_model','inquire_dimension soil_layers_stag '//trim(hydro_netcdf_filename))
else
   n_hydro_layers = nLayers
endif

if (nsoil /= nLayers .or. nsoil /= n_hydro_layers) then
   if (nsoil /= nLayers) write(string1,*) 'Expected ',nsoil,' soil layers [', &
        trim(lsm_netcdf_filename),'] has ',nLayers
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
   if (n_hydro_layers /= nLayers) &
      write(string1,*) 'LSM and HYDRO components do not have the same number of soil layers.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif

!---------------------------------------------------------------
! Compile the list of NOAH variables to use in the creation
! of the DART state vector. Required to determine model_size.
!
! Verify all variables are in the NOAH netcdf file.
! Compute the offsets into the state vector for each variable type.
! Record the extent of the variable type in the state vector.

!! nfields = number of state variables (NOT number of total locations * state variables)
call verify_state_variables(iunit_lsm, lsm_netcdf_filename, lsm_state_variables, n_lsm_fields)

n_hydro_fields = 0
if (hydro_model_active) then
   call verify_state_variables(iunit_hydro, hydro_netcdf_filename, &
        hydro_state_variables, n_hydro_fields)
endif

n_assimOnly_fields = 0
if (assimOnly_active) then 
   call nc_check(nf90_open(adjustl(assimOnly_netcdf_filename), NF90_NOWRITE, iunit_assimOnly), &
        'static_init_model', 'open '//trim(assimOnly_netcdf_filename))
   call verify_state_variables(iunit_assimOnly, assimOnly_netcdf_filename, &
        assimOnly_state_variables, n_assimOnly_fields)
endif

!! Dealing with soil mositure is complicated by
!!  1) duplication between LSM and HYDRO restarts
!!  2) change of scale from LSM to HYDRO involves several additional restart
!!     variables which are then required to be brought in (but kept from DART).
!! If rst_type=1 then LSM restart soil mositure is ignored.
!! However we want to be able to work with either the LSM or the LSM+hydro.

!! 1) deal with duplication.
!! If SMC is in *both* LSM and HYDRO, we want to only keep the HYDRO values for DART. (Though
!! we will write these back in to both restart files, it's not necessary if rst_typ=1 in the
!! hydro.namelist). Getting rid of one copy essentially means we are assuming rst_typ=1 in
!! hydro.namelist, so we'll enforce this when removing the LSM copy.

allocate(keepLsmVars0(n_lsm_fields))
keepLsmVars0 = (/ (i, i=1,n_lsm_fields) /)   ! TJH do you mean 1, i=1,n_lsm_fields ... summing ...
! use any( )  jlm fixme
lsmSmcPresent =  sum( keepLsmVars0 , mask = (lsm_state_variables(1,:) .eq. 'SOIL_W') .or. &
     (lsm_state_variables(1,:) .eq. 'SH2O') )

hydroSmcPresent = 0
hydroSfcHeadPresent = 0
if (hydro_model_active) then
   hydroSmcPresent     = any( hydro_state_variables(1,:) .eq. 'sh2ox' )
   hydroSfcHeadPresent = any( hydro_state_variables(1,:) .eq. 'sfcheadrt' )
   if (hydroSmcPresent)     allocate(sh2oDisag(fine3dShape(1),fine3dShape(2),fine3dShape(3)))
   if (hydroSfcHeadPresent) allocate(sfcHeadDisag(fine2dShape(1),fine2dShape(2)))
   if (hydroSmcPresent .OR. hydroSfcHeadPresent) call disagHydro()
endif


if (lsmSmcPresent > 0 .and. hydroSmcPresent > 0) then
   if (rst_typ /= 1) then
      write(string1,*) 'Seems BAD: Using hydro SMC but rst_type != 1!'
      call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
   endif
   allocate(keepLsmVars(n_lsm_fields-1))
   keepLsmVars = pack( keepLsmVars0, mask = (keepLsmVars0 /= lsmSmcPresent) )
   if (do_output() .and. (debug > 0)) then
      write(logfileunit,*)
      write(     *     ,*)
      write(logfileunit,*)"Removing LSM SMC from combined state vector."
      write(     *     ,*)"Removing LSM SMC from combined state vector."
   endif
else
   allocate(keepLsmVars(n_lsm_fields))
   keepLsmVars = keepLsmVars0
end if
n_lsm_fields = size(keepLsmVars)

! 2) deal with change of spatial resolution (disag) for liquid soil moisture in the
!    hydro restart file. Soil ice is not treated by the hydro model, so it's not allowed
!    as an uncertain state. (Would have to adjust in the LSM instead).
!    Right now this only deals with soil moisture, it may be desirable to add other variables
!    within this section(?).

if (hydro_model_active) then
   if (any(hydro_state_variables .eq. 'smc')) then
      write(string1,*) 'cannot adjust total soil moisture (smc) in the hydro model'
      call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
   end if
endif


!! make a combined list of: state variables, types, component
nfields = n_lsm_fields + n_hydro_fields + n_assimOnly_fields
allocate(all_state_variables(3,nfields))  !! additional column to identify the model component
all_state_variables(1:2,1:n_lsm_fields) = lsm_state_variables(1:2,keepLsmVars)
all_state_variables(3,1:n_lsm_fields) = 'LSM'
if (hydro_model_active) then
   all_state_variables(1:2,(n_lsm_fields+1):(n_lsm_fields+n_hydro_fields)) = hydro_state_variables
   all_state_variables(  3,(n_lsm_fields+1):(n_lsm_fields+n_hydro_fields)) = 'HYDRO'
endif
if (assimOnly_active) then 
   all_state_variables(1:2,(n_lsm_fields+n_hydro_fields+1): &
        (n_lsm_fields+n_hydro_fields+n_assimOnly_fields)) = &
        assimOnly_state_variables
   all_state_variables(  3, (n_lsm_fields+n_hydro_fields+1): &
        (n_lsm_fields+n_hydro_fields+n_assimOnly_fields)) = 'assimOnly'
endif


if (do_output() .and. (debug > 0)) then
   write(logfileunit,*)
   write(     *     ,*)
   write(logfileunit,*)"Combined state variables:"
   write(     *     ,*)"Combined variables:"
   do i=1,nfields
      write(logfileunit,*)i,': ',trim(all_state_variables(1,i)), '  ', &
           trim(all_state_variables(2,i)), '  ',trim(all_state_variables(3,i))
      write(     *     ,*)i,': ',trim(all_state_variables(1,i)), '  ', &
           trim(all_state_variables(2,i)), '  ',trim(all_state_variables(3,i))
   end do
end if

index1  = 1
FILL_PROGVAR : do ivar = 1, nfields

   varname                   = trim(all_state_variables(1,ivar))
   kind_string               = trim(all_state_variables(2,ivar))
   component                 = trim(all_state_variables(3,ivar))
   progvar(ivar)%varname     = varname
   progvar(ivar)%kind_string = kind_string
   progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
   progvar(ivar)%component   = component
   progvar(ivar)%grid   = ' '
   progvar(ivar)%dimlens     = 0
   progvar(ivar)%dimnames    = ' '
   progvar(ivar)%maxlevels   = 0  !jlm fix this

   if (component .eq. 'HYDRO') then
      iunit = iunit_hydro
      string2 = 'Hydro Restart File - '//trim(varname)
   elseif (component .eq. 'assimOnly') then 
      iunit = iunit_assimOnly
      string2 = 'assimOnly Restart File - '//trim(varname)
   else
      iunit = iunit_lsm
      string2 = 'LSM restart file - '//varname
   end if

   call nc_check(nf90_inq_varid(iunit, trim(varname), VarID), &
        'static_init_model', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(iunit, VarID, dimids=dimIDs, &
        ndims=progvar(ivar)%numdims, xtype=progvar(ivar)%xtype), &
        'static_init_model', 'inquire '//trim(string2))

   ! If the long_name and/or units attributes are set, get them.
   ! They are not REQUIRED to exist but are nice to use if they are present.

   if( nf90_inquire_attribute(    iunit, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(iunit, VarID, 'long_name' , progvar(ivar)%long_name), &
           'static_init_model', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = varname
   endif

   if( nf90_inquire_attribute(    iunit, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(iunit, VarID, 'units' , progvar(ivar)%units), &
           'static_init_model', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = '-'
   endif

   if( nf90_inquire_attribute(    iunit, VarID, 'stagger') == NF90_NOERR )  then
      call nc_check( nf90_get_att(iunit, VarID, 'stagger' , progvar(ivar)%stagger), &
           'static_init_model', 'get_att stagger '//trim(string2))
   else
      progvar(ivar)%stagger = '-'
   endif

   ! This is fundamentally hardcoded clamping. See WRF/model_mod.f90 for a namelist-driven
   ! example.  It would be nice to use the netCDF file valid_range attribute ...
   ! JLM - It would be nice if this information were somehow built into the OBS_KIND, since that's
   ! required anyway. I'd prefer not to have another namelist.
   !
   ! if the variable is bounded, then we need to know how to restrict it.
   ! rangeRestricted == 0 is unbounded
   ! rangeRestricted == 1 is bounded below
   ! rangeRestricted == 2 is bounded above           
   ! rangeRestricted == 3 is bounded above and below 
   ! rangeRestricted == 4 means apply nonscalar restrictions to this field as 
   !                    coded in vector_to_nd_progvar for the appropriate dimension. 

   progvar(ivar)%rangeRestricted = 0
   progvar(ivar)%minvalue        = -1.0_r8*huge(0.0_r8)
   progvar(ivar)%maxvalue        = huge(0.0_r8)

   ! set minvalue for certain variables which deserve it.
   if ( varname == 'SNODEP'        .or. &    ! noah
        varname == 'WEASD'         .or. &
        varname == 'SNOWH'         .or. &    ! noahMP
        varname == 'SNEQV'         .or. &
        varname == 'SNOWC'         .or. &
        varname == 'SOIL_T'        .or. &
        varname == 'CANLIQ'        .or. &
        varname == 'LAI'           .or. &
        varname == 'WA'            .or. &
        varname == 'qlink1'        .or. &    ! hydro
        varname == 'hlink'         .or. &
        varname == 'cvol'          .or. &
        varname == 'z_gwsubbas'    .or. &
        varname == 'sfcheadrt'     .or. &
        varname == 'infxswgt'      .or. &
        varname == 'precipMult'    .or. &    ! assimOnly
        varname == 'OVROUGHRTFAC'  .or. &
        varname == 'RETDEPRTFAC'   .or. &
        varname == 'gwCoeff'       .or. &
        varname == 'gwExpon'       .or. &
        varname == 'ksatMult'      .or. &
        varname == 'maxSmcMult'    .or. &
        varname == 'bbMult'        .or. &
        varname == 'satPsiMult'    .or. &
        varname == 'slope'         .or. &
        varname == 'refkdt'         .or. &
        varname == 'rsMult'        .or. &
        varname == 'ch2opMult'     .or. &
        varname == 'czil'               &
        ) then
      progvar(ivar)%rangeRestricted = 1
      ! Most vars have minvalue of 0.
      progvar(ivar)%minvalue     = 0
      ! Write exceptions here, this is a dummy example.
      if ( varname == 'SOIL_T' )        progvar(ivar)%minvalue     = 100  !! Kelvin. 
      if ( varname == 'precipMult' )    progvar(ivar)%minvalue     = .01
      if ( varname == 'OVROUGHRTFAC' )  progvar(ivar)%minvalue     = .01
      if ( varname == 'RETDEPRTFAC' )   progvar(ivar)%minvalue     = .01
      if ( varname == 'gwCoeff' )       progvar(ivar)%minvalue     = .01
      if ( varname == 'gwExpon' )       progvar(ivar)%minvalue     = .01
      if ( varname == 'ksatMult' )      progvar(ivar)%minvalue     = .01
      if ( varname == 'maxSmcMult' )    progvar(ivar)%minvalue     = .01
      if ( varname == 'bbMult' )        progvar(ivar)%minvalue     = .1
      if ( varname == 'satPsiMult' )    progvar(ivar)%minvalue     = .001
      if ( varname == 'slope' )         progvar(ivar)%minvalue     = .01
      if ( varname == 'refkdt' )         progvar(ivar)%minvalue    = .5
      if ( varname == 'rsMult' )        progvar(ivar)%minvalue     = .3  
      if ( varname == 'ch2opMult' )     progvar(ivar)%minvalue     = 0.001  
      if ( varname == 'czil' )          progvar(ivar)%minvalue     = .01 
   end if

   ! set maxvalue for those deserving variables. 
   if ( varname == 'SOIL_T'       .or. &  ! assimOnly
        varname == 'precipMult'   .or. &  ! assimOnly
        varname == 'OVROUGHRTFAC' .or. &   
        varname == 'RETDEPRTFAC'  .or. &
        varname == 'gwCoeff'      .or. &
        varname == 'gwExpon'      .or. &
        varname == 'ksatMult'     .or. &
        varname == 'maxSmcMult'   .or. &
        varname == 'bbMult'       .or. &
        varname == 'satPsiMult'   .or. &
        varname == 'slope'        .or. &
        varname == 'refkdt'        .or. &
        varname == 'rsMult'       .or. &
        varname == 'ch2opMult'    .or. &
        varname == 'czil'              &
        ) then
      progvar(ivar)%rangeRestricted = progvar(ivar)%rangeRestricted + 2
      progvar(ivar)%maxvalue        = 14.2  !! total dummy. 
      if ( varname == 'SOIL_T' )        progvar(ivar)%maxvalue = 273.15+75.  ! Kelvin. 
      if ( varname == 'precipMult' )    progvar(ivar)%maxvalue = 10.
      if ( varname == 'OVROUGHRTFAC' )  progvar(ivar)%maxvalue = 20.  ! fairly arbitrary
      if ( varname == 'RETDEPRTFAC' )   progvar(ivar)%maxvalue = 20.
      if ( varname == 'gwCoeff' )       progvar(ivar)%maxvalue = 100. ! fairly arbitrary
      if ( varname == 'gwExpon' )       progvar(ivar)%maxvalue = 100.
      if ( varname == 'ksatMult' )      progvar(ivar)%maxvalue = 10.
      if ( varname == 'maxSmcMult' )    progvar(ivar)%maxvalue = 1.25
      if ( varname == 'bbMult' )        progvar(ivar)%maxvalue = 3.
      if ( varname == 'satPsiMult' )    progvar(ivar)%maxvalue = 7.
      if ( varname == 'slope' )         progvar(ivar)%maxvalue = 1.
      if ( varname == 'refkdt' )        progvar(ivar)%maxvalue = 5.
      if ( varname == 'rsMult' )        progvar(ivar)%maxvalue = 1.2
      if ( varname == 'ch2opMult' )     progvar(ivar)%maxvalue = 5.
      if ( varname == 'czil' )          progvar(ivar)%maxvalue = 1
   end if

   ! specify if array bounds are to be used. 
   if ( varname == 'SOIL_W' .or. &   ! noah
        varname == 'SH2O'   .or. &   ! noahMP
        varname == 'sh2ox'       &   ! hydro
        ) then
      progvar(ivar)%rangeRestricted = 4
   end if
   
   ! These variables have a Time dimension. We only want the most recent time.
   varsize = 1
   dimlen  = 1
   DimensionLoop : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)

      call nc_check(nf90_inquire_dimension(iunit, dimIDs(i), name=dimname, len=dimlen), &
           'static_init_model', string1)

      ! These variables have a Time dimension. We only want the most recent time.
      if ((trim(dimname) == 'Time') .or. (trim(dimname) == 'time')) dimlen = 1

      ! the hack on the following line allows us to change dimensions from file storage to 
      ! disaggregated dimensions for hydro sh2ox and use the disaggregated in the assim.
      ! use dimnames to identify changing the size instead... jlm fixme

      if (trim(component) .eq. 'HYDRO') then 
         if (trim(varname) .eq. 'sh2ox')       dimlen = fine3dShape(i)
         if (trim(varname) .eq. 'sfcheadrt')   dimlen = fine2dShape(i)
      endif

      progvar(ivar)%dimlens(i) = dimlen
      progvar(ivar)%dimnames(i) = trim(dimname)
      varsize = varsize * dimlen
   enddo DimensionLoop

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1

   ! Determine the grid type, but fill the matching state vectors outside this loop.
   ! FIXME jlm - what if the different grids are of equal size? 
   ! Use dimension names instead of size.
   ! FIXME TJH - the dimnames are supposed to cover that.

   dumNumDims = progvar(ivar)%numdims
   if (progvar(ivar)%component == 'LSM') dumNumDims = dumNumDims-1  !LSM has time dim of length 1
   if (dumNumDims .eq. 1) then
      if(progvar(ivar)%component == 'HYDRO') then 
         if (progvar(ivar)%varsize .eq. n_link) progvar(ivar)%grid = 'link1d'
         if (progvar(ivar)%varsize .eq. n_basn) progvar(ivar)%grid = 'basn1d'
      end if
      if(progvar(ivar)%component == 'assimOnly') then 
         if (progvar(ivar)%varsize .eq. 1) progvar(ivar)%grid = 'scalar'
      end if
   endif
   if (dumNumDims .eq. 2) then
      if (progvar(ivar)%varsize .eq. n_hlong*n_hlat) progvar(ivar)%grid = 'fine2d'
      if (progvar(ivar)%varsize .eq. south_north*west_east) progvar(ivar)%grid = 'coarse2d'
   endif
   if (dumNumDims .eq. 3) then
      if (progvar(ivar)%varsize .eq. n_hlong*n_hlat*nsoil) progvar(ivar)%grid = 'fine3d'
      if (progvar(ivar)%varsize .eq. south_north*west_east*nsoil) progvar(ivar)%grid = 'coarse3d'
   endif

   ! Memory ORDER
   ! TJH doesn't this change if the variable comes from noah vs. noahMP

   if( nf90_inquire_attribute(    iunit, VarID, 'MemoryOrder') == NF90_NOERR )  then
      call nc_check( nf90_get_att(iunit, VarID, 'MemoryOrder' , progvar(ivar)%MemoryOrder), &
           'static_init_model', 'get_att MemoryOrder '//trim(string2))
   else
      progvar(ivar)%MemoryOrder = ''   ! lazy else statment
      ! TJH could be a parameter ... or something ... true or ...
      if(progvar(ivar)%numdims == 1) progvar(ivar)%MemoryOrder = 'X'
      if(progvar(ivar)%numdims == 2) progvar(ivar)%MemoryOrder = 'XY'
      if(progvar(ivar)%numdims == 3) then 
         ! jlm fixme - use dimname instead of length - or just use lsm_model_choice
         ! This is Noah
         if (progvar(ivar)%dimlens(3) .eq. nsoil) progvar(ivar)%MemoryOrder = 'XYZ'
         ! This is NoahMP
         if (progvar(ivar)%dimlens(2) .eq. nsoil) progvar(ivar)%MemoryOrder = 'XZY'
      endif
   endif

   progvar(ivar)%gridMemOrd = trim(progvar(ivar)%grid) // trim(progvar(ivar)%MemoryOrder)

   if ( (do_output()) .and. debug > 10 ) then
      write(logfileunit,*)
      write(logfileunit,*) 'variable number ', ivar, ' is ',trim(progvar(ivar)%varname)
      write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) '  component   ',trim(progvar(ivar)%component)
      write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
      write(logfileunit,*) '  MemoryOrder ',trim(progvar(ivar)%MemoryOrder)
      write(logfileunit,*) '  grid        ',trim(progvar(ivar)%grid)
      write(logfileunit,*) '  GridMemOrd  ',trim(progvar(ivar)%gridMemOrd)
      write(logfileunit,*) '  stagger     ',trim(progvar(ivar)%stagger)
      write(logfileunit,*) '  dimnames    ',progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
      write(logfileunit,*) '  dimlens     ',progvar(ivar)%dimlens( 1:progvar(ivar)%numdims)
      write(logfileunit,*) '  numdims     ',progvar(ivar)%numdims
      write(logfileunit,*) '  maxlevels   ',progvar(ivar)%maxlevels
      write(logfileunit,*) '  xtype       ',progvar(ivar)%xtype
      write(logfileunit,*) '  varsize     ',progvar(ivar)%varsize
      write(logfileunit,*) '  index1      ',progvar(ivar)%index1
      write(logfileunit,*) '  indexN      ',progvar(ivar)%indexN
      write(logfileunit,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(logfileunit,*) '  restriction ',progvar(ivar)%rangeRestricted
      write(logfileunit,*) '  minvalue    ',progvar(ivar)%minvalue
      write(logfileunit,*) '  maxvalue    ',progvar(ivar)%maxvalue
      write(logfileunit,*) '  kind_string ',progvar(ivar)%kind_string

      write(     *     ,*)
      write(     *     ,*) 'variable number ', ivar, ' is ',trim(progvar(ivar)%varname)
      write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(     *     ,*) '  component   ',trim(progvar(ivar)%component)
      write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
      write(     *     ,*) '  MemoryOrder ',trim(progvar(ivar)%MemoryOrder)
      write(     *     ,*) '  grid        ',trim(progvar(ivar)%grid)
      write(     *     ,*) '  GridMemOrd  ',trim(progvar(ivar)%gridMemOrd)
      write(     *     ,*) '  stagger     ',trim(progvar(ivar)%stagger)
      write(     *     ,*) '  dimnames    ',progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
      write(     *     ,*) '  dimlens     ',progvar(ivar)%dimlens( 1:progvar(ivar)%numdims)
      write(     *     ,*) '  numdims     ',progvar(ivar)%numdims
      write(     *     ,*) '  maxlevels   ',progvar(ivar)%maxlevels
      write(     *     ,*) '  xtype       ',progvar(ivar)%xtype
      write(     *     ,*) '  varsize     ',progvar(ivar)%varsize
      write(     *     ,*) '  index1      ',progvar(ivar)%index1
      write(     *     ,*) '  indexN      ',progvar(ivar)%indexN
      write(     *     ,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(     *     ,*) '  restriction ',progvar(ivar)%rangeRestricted
      write(     *     ,*) '  minvalue    ',progvar(ivar)%minvalue
      write(     *     ,*) '  maxvalue    ',progvar(ivar)%maxvalue
      write(     *     ,*) '  kind_string ',progvar(ivar)%kind_string
   endif

   ! sets up for next variable
   index1 = index1 + varsize

enddo FILL_PROGVAR

call nc_check(nf90_close(iunit_lsm), 'static_init_model', &
                  'close '//trim(lsm_netcdf_filename))
if (hydro_model_active) then
   call nc_check(nf90_close(iunit_hydro), 'static_init_model', &
                      'close '//trim(hydro_netcdf_filename))
end if
if (assimOnly_active) then
   call nc_check(nf90_close(iunit_assimOnly), 'static_init_model', &
                     'close '//trim(assimOnly_netcdf_filename))
end if

model_size = progvar(nfields)%indexN

!-------------------------------------------------------------------------------
! State location and grid information
! According to DART/lanai/models/WRF_Hydro/work/path_names_model_mod_check
! we are using a 3dsphere location module.
! Convert soil thicknesses (from namelist.hrldas) to "heights" (VERTISHEIGHT)..
! Closer to the center of the earth is an increasingly large negative number
! Do the appropriate summing of layer thicknesses and multiplying by -1 for noahMP.

allocate(state_lon(model_size),state_lat(model_size),state_level(model_size))
!-------------------------------------------------------------------------------
! Fill the dart/state coordinate vectors: state_lon, state_lat, and state_level (also elevation?)
! TJH FIXME - I cannot envision a scenario where the elevation relative to the geoid would be useful.

! TJH left off here ... 
! The filling of the coordinate vectors is highly sensitive to the manner in which the 
! variables are packed into the state vector.

! Make a unique list of progvar%gridMemOrd
allocate(uniqueGridMemOrd(nfields))
nUniqueGrids = 0
do ivar = 1,nfields
   dumGridMemOrd = progvar(ivar)%gridMemOrd
   if (.not. any(uniqueGridMemOrd == dumGridMemOrd)) then 
      nUniqueGrids = nUniqueGrids + 1      
      uniqueGridMemOrd(nUniqueGrids) = dumGridMemOrd
   endif
enddo

!print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!print*,nUniqueGrids
!do i=1,nUniqueGrids 
!   print*,trim(uniqueGridMemOrd(i))
!end do

! Break this into sub routines named after grids, clarify the overall strategy, 
! jlm fixme - ditch identifying unique grids, and just loop through progvar
do igrid = 1, nUniqueGrids
   
   ! identify the progvars with that grid
   allocate(whichVars( sum( (/ (1, i=1,nfields) /), &
                           mask=progvar%gridMemOrd .eq. trim(uniqueGridMemOrd(igrid)) )) )
   whichVars = pack( (/ (i, i=1,size(progvar%gridMemOrd)) /), &
                     mask = progvar%gridMemOrd .eq. trim(uniqueGridMemOrd(igrid)) )  !! shouldnt need a trim on progvar%grid

      ! create a dummy copy of those coords in 1D for the current grid
   if (size(whichVars) .gt. 0) then

      whVar1 = whichVars(1)  !! since all variables identified have same grid and memoryOrder, reference against the first.
      ! 1D
      if( trim(progvar(whVar1)%grid) .eq. 'link1d' ) then  
         dumSize=n_link
         allocate(dumGridLon(dumSize),dumGridLat(dumSize),dumGridLevel(dumSize))
         dumGridLon = linkLong
         dumGridLat = linkLat
         dumGridLevel = linkLong*0.0_R8 ! surface
      endif

      ! these "grids" dont have coordinates.
      if( (trim(progvar(whVar1)%grid) .eq. 'basn1d') .or. (trim(progvar(whVar1)%grid) .eq. 'scalar') ) then
         if( trim(progvar(whVar1)%grid) .eq. 'basn1d' ) dumSize=n_basn
         if( trim(progvar(whVar1)%grid) .eq. 'scalar' ) dumSize=1
         allocate(dumGridLon(dumSize),dumGridLat(dumSize),dumGridLevel(dumSize))
         dumGridLon(:) = basnLon(1:dumSize)
         dumGridLat(:) = basnLat(1:dumSize)
         dumGridLevel = 0.0_R8 ! surface
      endif

      !! 2D grids
      if (trim(progvar(whVar1)%grid) .eq. 'fine2d' .or. trim(progvar(whVar1)%grid) .eq. 'coarse2d') then

         if (trim(progvar(whVar1)%grid) .eq. 'coarse2d') then
            dumNLon=west_east
            dumNLat=south_north
         endif
         if (trim(progvar(whVar1)%grid) .eq. 'fine2d') then
            dumNLon = n_hlong !fine grid is default
            dumNLat = n_hlat
         end if

         dumSize=dumNLon*dumNLat
         allocate(dumGridLon(dumSize),dumGridLat(dumSize),dumGridLevel(dumSize))
         myindex=1
         do ilat=1,dumNLat    !Y
         do ilon=1,dumNLon    !X
            if (trim(progvar(whVar1)%grid) .eq. 'fine2d') then
               dumGridLon(myindex) = hlong(ilon, ilat)
               dumGridLat(myindex) =  hlat(ilon, ilat)
            endif
            if (trim(progvar(whVar1)%grid) .eq. 'coarse2d') then
               dumGridLon(myindex) = xlong(ilon, ilat)
               dumGridLat(myindex) =  xlat(ilon, ilat)
            endif
            myindex = myindex + 1
         end do              !X
         end do              !Y
         dumGridLevel = dumGridLat*0.0_R8 ! surface
      end if

      ! 3D
      if (trim(progvar(whVar1)%grid) .eq. 'fine3d'  .or.   &
          trim(progvar(whVar1)%grid) .eq. 'coarse3d') then

         if (trim(progvar(whVar1)%grid) .eq. 'coarse3d') then
            dumNLon=west_east
            dumNLat=south_north
         endif
         if (trim(progvar(whVar1)%grid) .eq. 'fine3d') then
            dumNLon = n_hlong !fine grid is default
            dumNLat = n_hlat
         end if

         dumSize=dumNLon*dumNLat*nsoil
         allocate(dumGridLon(dumSize),dumGridLat(dumSize),dumGridLevel(dumSize))
         myindex=1
         
         ! XYZ order - Noah or HYDRO
         if (trim(progvar(whVar1)%MemoryOrder) .eq. 'XYZ') then           
            do ilev=1,nsoil     !Z
            do ilat=1,dumNLat   !Y
            do ilon=1,dumNLon   !X
               if (trim(progvar(whVar1)%grid) .eq. 'fine3d') then
                  dumGridLon(myindex) = hlong(ilon, ilat)
                  dumGridLat(myindex) =  hlat(ilon, ilat)            
               endif
               if (trim(progvar(whVar1)%grid) .eq. 'coarse3d') then
                  dumGridLon(myindex) = xlong(ilon, ilat)
                  dumGridLat(myindex) =  xlat(ilon, ilat)
               endif
               dumGridLevel(myindex) = zsoil8(ilev)
               myindex = myindex + 1
            end do              !X
            end do              !Y
            end do              !Z
         end if !XYZ

         ! XZY order - NoahMP
         if (trim(progvar(whVar1)%MemoryOrder) .eq. 'XZY') then
            do ilat=1,dumNLat   !Y
            do ilev=1,nsoil     !Z
            do ilon=1,dumNLon   !X
               if (trim(progvar(whVar1)%grid) .eq. 'fine3d') then
                  dumGridLon(myindex) = hlong(ilon, ilat)
                  dumGridLat(myindex) =  hlat(ilon, ilat)
               endif
               if (trim(progvar(whVar1)%grid) .eq. 'coarse3d') then
                  dumGridLon(myindex) = xlong(ilon, ilat)
                  dumGridLat(myindex) =  xlat(ilon, ilat)
               endif
               dumGridLevel(myindex) = zsoil8(ilev)
               !dumGridLevel(myindex) = -1*sum(soil_thick_input(1:ilev))
               myindex = myindex + 1
            end do             !X
            end do             !Z
            end do             !Y
         end if !XZY

      end if !3D

      ! should we enforce that the grids are ([levels,] lat, lon) for each (non 1-d) progvar?

      ! Fill the state vector using the dumGrid variables
      do ivar=1,size(whichVars)
         wp=whichVars(ivar)
         state_lon(progvar(wp)%index1:progvar(wp)%indexN) = dumGridLon
         state_lat(progvar(wp)%index1:progvar(wp)%indexN) = dumGridLat
         state_level(progvar(wp)%index1:progvar(wp)%indexN) = dumGridLevel
      end do

      deallocate(dumGridLon,dumGridLat,dumGridLevel)
   end if
enddo

! Make sure we are [0,360] and [-90,90]
! jlm fixme I thin it's better to require these to be done for the 
! individdual coordinate systems upon read than now, just in the state
! vector because the coords read in get pushed to the Post/Prior diagnostic files. 
! I think it's better to throw an error here than to correct one, since it may leave
! inconsistencies in the output.
if (any(state_lon < 0.0_r8)) then
   write(string1,*)'longitudes in "state_lon" still negative.' // &
        ' Please fix these immediately upon reading in the coordinates or prior.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif
if (any(state_lat < -90.0_r8) .or. any(state_lat > 90.0_r8) ) then
   write(string1,*)'Latitudes in "state_lat" outside [-90,90].' // &
        ' Please fix these immediatley upon reading in the coordinates or prior.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif

allocate(state_loc(model_size))
do iState=1,model_size
   state_loc(iState) = set_location(   state_lon(iState), state_lat(iState), &
                                     state_level(iState), VERTISHEIGHT)
enddo


deallocate(state_lon, state_lat, state_level, &
           keepLsmVars0, keepLsmVars, uniqueGridMemOrd, whichVars)

end subroutine static_init_model


!===============================================================================
subroutine init_conditions(x)
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.
real(r8), intent(out) :: x(:)
if ( .not. module_initialized ) call static_init_model
write(string1,*) 'PROBLEM: no known way to set arbitrary initial conditions.'
write(string2,*) 'start_from_restart must be .true. in this model.'
call error_handler(E_ERR,'init_conditions',string1,source,revision,revdate, &
     text2=string2)
x = MISSING_R8

end subroutine init_conditions



!===============================================================================
subroutine adv_1step(x, time)
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to
! compute a timestep, for instance for radiation computations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a
! NULL INTERFACE.
real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time
if ( .not. module_initialized ) call static_init_model  
write(string1,*) 'PROBLEM: cannot advance model with async == 0.'
write(string2,*) 'async == 2 is a good choice.'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate, &
     text2=string2)
end subroutine adv_1step


!===============================================================================
function get_model_size()
! Returns the size of the model as an integer. Required for all
! applications.
integer :: get_model_size
if ( .not. module_initialized ) call static_init_model
get_model_size = model_size
end function get_model_size


!===============================================================================
subroutine init_time(time)
! Companion interface to init_conditions. Returns a time that is somehow
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.
type(time_type), intent(out) :: time
if ( .not. module_initialized ) call static_init_model
time = model_time
end subroutine init_time


!===============================================================================
subroutine model_interpolate(x, location, itype, obs_val, istatus)
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field to that location and returns
! the value in obs_val. The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the type of field (for
! instance temperature, zonal wind component, etc.).
real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

real(r8)               :: loc_lon, loc_lat, loc_depth
real(r8),  allocatable, dimension(:) :: stateVarLon, stateVarLat, stateVarDepth, dists
real(r8), dimension(3) :: loc
integer,  dimension(1) :: closestIndArr
integer                :: closestInd
integer                :: zlev, n, ivar, indx, nPts, index1, indexN, iPt

if ( .not. module_initialized ) call static_init_model

! FIXME - for the single column case - there is no obvious way to
! determine the extent of the domain ... EVERYTHING matches.

! can likely take this out, though I'm not sure what the above comment about 
! everything matching means - a check of some sort?
!if( west_east*south_north /= 1 ) then
!     write(string1,*) 'PROBLEM: not set up for a case with multiple locations yet.'
!     call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
!endif

! The return code for successful return should be 0.
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.

istatus = 1
obs_val = MISSING_R8

! if observation is outside region encompassed in the history file - fail
loc       = get_location(location) ! loc is in DEGREES
loc_lon   = loc(1)
loc_lat   = loc(2)
loc_depth = loc(3)

! need to know what variable we are interpolating.
ivar = 0
FindVariable : do n = 1,nfields
   if( progvar(n)%dart_kind == itype ) then
      ivar = n
      exit FindVariable
   endif
enddo FindVariable

if (ivar == 0) then
   write(string1,*) 'unable to find state vector component matching type ',itype
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
endif

index1 = progvar(ivar)%index1 ! in the DART state vector, start looking here
indexN = progvar(ivar)%indexN ! in the DART state vector, stop  looking here
nPts = indexN-index1+1

!jlm not currently setting maxlevels but I should.
! BOMBPROOFING - check for a vertical dimension for this variable
!if (progvar(ivar)%maxlevels < 2) then
!   write(string1, *)'Variable '//trim(varstring)//' should not use this routine.'
!   write(string2, *)'use compute_gridcell_value() instead.'
!   call error_handler(E_ERR,'get_grid_vertval', string1, &
!                  source, revision, revdate, text2=string2)
!endif

allocate(stateVarLon(nPts), stateVarLat(nPts), &
     stateVarDepth(nPts), dists(nPts))
do iPt=1,nPts
   loc = get_location(state_loc(index1+iPt-1))
   stateVarLon(iPt)   = loc(1)
   stateVarLat(iPt)   = loc(2)
   stateVarDepth(iPt) = loc(3)
enddo

dists = sqrt( (stateVarLat  -loc_lat  )**(2.) + &
              (stateVarLon  -loc_lon  )**(2.) )
!             (stateVarDepth-loc_depth)**(2.) )  !! 3d
closestIndArr = minloc( dists )
closestInd = closestIndArr(1) + index1 - 1

obs_val = x( closestInd )

if (obs_val /= MISSING_R8) istatus = 0

if ( (do_output()) .and. debug > 20 ) then
   write(*,*)'model_interpolate : progvar%kind_string is ',trim(progvar(ivar)%kind_string)
   write(*,*)'model_interpolate : state index         is ',closestInd
   write(*,*)'model_interpolate : variable index         is ',closestInd-progvar(ivar)%index1+1
   write(*,*)'model_interpolate : value               is ',obs_val     
   call write_location(n,state_loc(closestInd),charstring=string1)
   write(*,*)'state location ',trim(string1)
   call write_location(n,location,charstring=string1)
   write(*,*)'observation location ',trim(string1)
endif

deallocate(stateVarLat, stateVarLon, stateVarDepth, dists)

end subroutine model_interpolate


!===============================================================================
function get_model_time_step()
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.
type(time_type) :: get_model_time_step
if ( .not. module_initialized ) call static_init_model
! The NOAH model can only be advanced in multiples of the restart frequency.
get_model_time_step = set_time(khour*3600,kday)
end function get_model_time_step


!===============================================================================
subroutine get_state_meta_data(index_in, location, var_type)
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.
integer,             intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: var_type

integer :: n, subindex, layer, ivar

if ( .not. module_initialized ) call static_init_model

FindIndex : do n = 1,nfields
   if( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) then
      ivar     = n
      exit FindIndex
   endif
enddo FindIndex

if (present(var_type)) var_type = progvar(ivar)%dart_kind


!if ( (do_output()) .and. debug > 30 ) then
!   write(*,*)'get_state_meta_data: index_in is ',index_in
!   write(*,*)'get_state_meta_data: ivar     is ',ivar
!   !write(*,*)'get_state_meta_data: layer    is ',layer
!   write(*,*)'get_state_meta_data: type    is ',var_type
!   write(*,*)
!endif

!if( layer == -1 ) then
!     write(string1,*) 'Problem, cannot find base_offset, index_in is: ', index_in
!     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
!endif

!if (progvar(ivar)%varsize == 1) layer = 0
location = state_loc(index_in)

end subroutine get_state_meta_data


!===============================================================================
subroutine end_model()
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.
! good style ... perhaps you could deallocate stuff (from static_init_model?).
if ( .not. module_initialized ) call static_init_model

deallocate(sh2oDisag, sfcHeadDisag,   &
           all_state_variables,       &
           channelIndsX, channelIndsY, linkLong, linkLat, &
           state_loc, &
           smcWlt1, smcMax1, sh2oWltRt, sh2oMaxRt, smc, sice, &
           hlong, hlat, basnMask, basnLon, basnLat)

end subroutine end_model


!===============================================================================
function nc_write_model_atts( ncFileID ) result (ierr)
! This routine writes all the netCDF 'infrastructure' and sets up the
! global attributes, dimensions, coordinate variables, and output variables.
! The actuall filling of the output variables is done by
! nc_write_model_vars() which can be called repeatedly for each
! assimilation cycle.
!
!!jlm this routine initis the ncdf file to hold the apriori and aposteriori
!! state vector either as a dum vector (output_state_vector) or as a smarter version.
!! but I'm confused, I' not seeing where this is ever used.
!! DART 'restart(?)' files are being written by assim_model_mod:awrite_restart in
!! wrfHydro_to_dart
!
! All errors are fatal, so the return code is always '0 == normal'
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN           ! open existing netCDF dataset
!    NF90_redef       ! put into define mode
!    NF90_def_dim     ! define additional dimensions (if any)
!    NF90_def_var     ! define variables: from name, type, and dims
!    NF90_put_att     ! assign attribute values
! NF90_ENDDEF         ! end definitions: leave define mode
!    NF90_put_var     ! provide values for variable
! NF90_CLOSE          ! close: save updated netCDF dataset

!! jlm fix : this routine is fairly rushed. there are many places where hydro_model_active
!! and assimOnly_active should be used. further more, hydro model could be active but 
!! a given variable/dimension might not be. Similarly for assimOnly_active. 
!! perhaps this is best handled by querying the progvar grid.

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!----------------------------------------------------
! variables if we just blast out one long state vector
!----------------------------------------------------

integer :: StateVarDimID    ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID      ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID        ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

!----------------------------------------------------
! variables if we parse the state vector into individual variables
!----------------------------------------------------
integer :: ixDimID, iyDimID, ixrtDimID, iyrtDimID, depthDimID
integer :: linksDimID, basnsDimID, scalarDimID
integer :: myndims
integer :: ivar, varID
integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
character(len=NF90_MAX_NAME) :: varname

!----------------------------------------------------
! variables for the namelist output
!----------------------------------------------------
! TJH FIXME ... the 129 has to match the declaration of the input.nml length,
! TJH FIXME .... which comes from someplace else ... 
character(len=129), allocatable, dimension(:) :: lsmTextblock, hydroTextblock
integer :: LineLenDimID, nlinesDimID, lsmNmlVarID, hydroNmlVarID
integer :: nlines, linelen
logical :: has_noah_namelist, has_hydro_namelist

!----------------------------------------------------
! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.
!----------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer :: i

character(len=128) :: filename

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
! and then put into define mode.
!-------------------------------------------------------------

ierr = -1 ! assume things go poorly

!--------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------
write(filename,*) 'ncFileID', ncFileID

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
     'nc_write_model_atts', 'inquire '//trim(filename))
call nc_check(nf90_redef(ncFileID), 'nc_write_model_atts', 'redef '//trim(filename))

!-------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension.
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='NMLlinelen', dimid = linelenDimID), &
     'nc_write_model_atts', 'inq_dimid NMLlinelen '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &  
     'nc_write_model_atts', 'inq_dimid copy '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid= TimeDimID), &
     'nc_write_model_atts', 'inq_dimid time '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
        ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!-------------------------------------------------------------
! Define the model size / state variable dimension 
!-------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable',  &
     len=model_size, dimid=StateVarDimID), &
     'nc_write_model_atts', 'def_dim state '//trim(filename))


!-------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------
call date_and_time(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
     values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,str1), &
     'nc_write_model_atts', 'put_att creation_date '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source), &
     'nc_write_model_atts', 'put_att model_source '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
     'nc_write_model_atts', 'put_att model_revision '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate), &
     'nc_write_model_atts', 'put_att model_revdate '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model','wrfHydro'), &  
     'nc_write_model_atts', 'put_att model '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'calendar',trim(calendar)), &
     'nc_write_model_atts', 'put_att calendar '//trim(filename))

!! LSM namelist -  check noah vs noahMP differences
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'HRLDAS_CONSTANTS_FILE',trim(hrldas_constants_file)), &
     'nc_write_model_atts', 'put_att HRLDAS_CONSTANTS_FILE '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'HRLDAS_INDIR',trim(INDIR)), &
     'nc_write_model_atts', 'put_att HRLDAS_INDIR '//trim(filename))

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'RESTART_FILENAME_REQUESTED ',trim(restart_filename_requested)), &
     'nc_write_model_atts', 'put_att RESTART_FILENAME_REQUESTED '//trim(filename))

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'KDAY',KDAY), &
     'nc_write_model_atts', 'put_att KDAY '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'KHOUR',KHOUR), &
     'nc_write_model_atts', 'put_att KHOUR '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'FORCING_TIMESTEP',forcing_timestep), &
     'nc_write_model_atts', 'put_att FORCING_TIMESTEP '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'NOAH_TIMESTEP',noah_timestep), &
     'nc_write_model_atts', 'put_att NOAH_TIMESTEP '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'OUTPUT_TIMESTEP',output_timestep), &
     'nc_write_model_atts', 'put_att OUTPUT_TIMESTEP '//trim(filename))

!! for noahMP will want to output options. 
!! jlm - more attributes for hydro namelist

!-------------------------------------------------------------
! Determine shape of namelist.
! long lines are truncated when read into textblock
!-------------------------------------------------------------
!! LSM
call find_textfile_dims(lsm_namelist_filename, nlines, linelen)
if (nlines > 0) then
   has_noah_namelist = .true.
else
   has_noah_namelist = .false.
endif
if (has_noah_namelist) then
   allocate(lsmTextblock(nlines))
   lsmTextblock = ''
   call nc_check(nf90_def_dim(ncid=ncFileID, name='noahNMLnlines', &
        len = nlines, dimid = nlinesDimID), &
        'nc_write_model_atts', 'def_dim noahNMLnlines '//trim(filename))
   call nc_check(nf90_def_var(ncFileID,name=trim(lsm_namelist_filename), xtype=nf90_char,    &
        dimids = (/ linelenDimID, nlinesDimID /),  varid=lsmNmlVarID), &
        'nc_write_model_atts', 'def_var noah_namelist '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, lsmNmlVarID, &
        'long_name', 'contents of '//trim(lsm_namelist_filename)), &
        'nc_write_model_atts', 'put_att noah_namelist '//trim(filename))
endif

!! Hydro
call find_textfile_dims(hydro_namelist_filename, nlines, linelen)
if (nlines > 0) then
   has_hydro_namelist = .true.
else
   has_hydro_namelist = .false.
endif
if (has_hydro_namelist) then
   allocate(hydroTextblock(nlines))
   hydroTextblock = ''
   call nc_check(nf90_def_dim(ncid=ncFileID, name='hydroNMLnlines', &
        len = nlines, dimid = nlinesDimID), &
        'nc_write_model_atts', 'def_dim hydroNMLnlines '//trim(filename))
   call nc_check(nf90_def_var(ncFileID,name=trim(hydro_namelist_filename), xtype=nf90_char,    &
        dimids = (/ linelenDimID, nlinesDimID /),  varid=hydroNmlVarID), &
        'nc_write_model_atts', 'def_var hydro_namelist '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, hydroNmlVarID, 'long_name', &
        'contents of '//trim(hydro_namelist_filename)), &
        'nc_write_model_atts', 'put_att hydro_namelist '//trim(filename))
endif

!-------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.
!-------------------------------------------------------------

if ( output_state_vector ) then

   !----------------------------------------------------------
   ! Create a variable into which the state vector is blasted
   !----------------------------------------------------------
   ! Define the state vector coordinate variable and some attributes.
   call nc_check(nf90_def_var(ncid=ncFileID,name='StateVariable', xtype=NF90_INT, &
        dimids=StateVarDimID, varid=StateVarVarID), &
        'nc_write_model_atts', 'def_var StateVariable')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID,'long_name','State Variable ID'), &
        'nc_write_model_atts', 'put_att StateVariable long_name')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'units',     'indexical'), &
        'nc_write_model_atts', 'put_att StateVariable units')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'valid_range', (/ 1, model_size /)), &
        'nc_write_model_atts', 'put_att StateVariable valid_range')

   ! Define the actual (3D) state vector, which gets filled as time goes on ...
   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=NF90_REAL, &
        dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
        varid=StateVarID), 'nc_write_model_atts', 'def_var state')
   call nc_check(nf90_put_att(ncFileID, StateVarID, 'long_name', 'model state or fcopy'), &
        'nc_write_model_atts', 'put_att state long_name')

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts', 'state_vector enddef')

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /)), &
        'nc_write_model_atts', 'put_var state')

else

   !----------------------------------------------------------
   ! We need to output the state vector as individual variables.
   !----------------------------------------------------------
   ! Define the additional dimensions IDs
   !----------------------------------------------------------
   !! I dont particularly like these variable names but I'm using the conventions in 
   !! the wrfHydro HYDRO_RST files

   !! Check which grids are actually being used and only output those. 
   !print*,'-------------------------------------------------------------------'
   !print*,progvar%grid     use this... TBD soon, fixme jlm

   !! coarse grid
   call nc_check(nf90_def_dim(ncid=ncFileID, name='ix', len=west_east, dimid=ixDimID),&
        'nc_write_model_atts','ix def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='iy', len=south_north, dimid=iyDimID), &
        'nc_write_model_atts', 'iy def_dim '//trim(filename))

   !! fine grid
   call nc_check(nf90_def_dim(ncid=ncFileID, name='ixrt', len=n_hlong, dimid=ixrtDimID), &
        'nc_write_model_atts','ixrt def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='iyrt', len=n_hlat, dimid=iyrtDimID), &
        'nc_write_model_atts','iyrt def_dim '//trim(filename))

   !! vertical
   call nc_check(nf90_def_dim(ncid=ncFileID, name='depth', len=nsoil, dimid=depthDimID), &
        'nc_write_model_atts', 'depth def_dim'//trim(filename))

   !! channel grid
   call nc_check(nf90_def_dim(ncid=ncFileID, name='links', len=n_link, dimid=linksDimID), &
        'nc_write_model_atts', 'links def_dim'//trim(filename))

   !! basin grid
   call nc_check(nf90_def_dim(ncid=ncFileID, name='basns', len=n_basn, dimid=basnsDimID), &
        'nc_write_model_atts', 'basns def_dim'//trim(filename))

   !! scalar "grid"
   call nc_check(nf90_def_dim(ncid=ncFileID, name='scalar', len=1, dimid=scalarDimID), &
        'nc_write_model_atts', 'scalar def_dim'//trim(filename))

   !----------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------
   ! Coarse Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lon', xtype=nf90_real, &
        dimids=(/ ixDimID, iyDimID /), varid=VarID),&
        'nc_write_model_atts', 'lon def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate longitude'), &
        'nc_write_model_atts', 'lon long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'coordinates', 'lon lat'),  &
        'nc_write_model_atts', 'lon coordinates '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'FieldType', 104),  &
        'nc_write_mode0l_atts', 'lon  FieldType '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'MemoryOrder', 'XY'),  &
        'nc_write_model_atts', 'lon MemoryOrder '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, '_FillValue', -999.0 ), &
        'nc_write_model_atts', 'lon _FillValue '//trim(filename))
   ! Coarse Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='lat', xtype=nf90_real, &
        dimids=(/ ixDimID, iyDimID /), varid=VarID),&
        'nc_write_model_atts', 'lat def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate latitude'), &
        'nc_write_model_atts', 'lat long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'coordinates', 'lon lat'),  &
        'nc_write_model_atts', 'lat coordinates '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'FieldType', 104),  &
        'nc_write_model_atts', 'lat FieldType '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'MemoryOrder', 'XY'),  &
        'nc_write_model_atts', 'lat MemoryOrder '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, '_FillValue', -999.0 ), &
        'nc_write_model_atts', 'lat _FillValue '//trim(filename))

   ! Fine/routing Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lonRt', xtype=nf90_real, &
        dimids=(/ ixrtDimID, iyrtDimID /), varid=VarID),&
        'nc_write_model_atts', 'lonRt def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'routing grid coordinate longitude'), &
        'nc_write_model_atts', 'lonRt long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'coordinates', 'lonRt latRt'),  &
        'nc_write_model_atts', 'lonRt coordinates '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'FieldType', 104),  &
        'nc_write_mode0l_atts', 'lonRt  FieldType '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'MemoryOrder', 'XY'),  &
        'nc_write_model_atts', 'lonRt MemoryOrder '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, '_FillValue', -999.0 ), &
        'nc_write_model_atts', 'lonRt _FillValue '//trim(filename))
   ! Fine/routing Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='latRt', xtype=nf90_real, &
        dimids=(/ ixrtDimID, iyrtDimID /), varid=VarID),&
        'nc_write_model_atts', 'latRt def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'routing grid coordinate latitude'), &
        'nc_write_model_atts', 'latRt long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'coordinates', 'lonRt latRt'),  &
        'nc_write_model_atts', 'latRt coordinates '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'FieldType', 104),  &
        'nc_write_model_atts', 'latRt FieldType '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'MemoryOrder', 'XY'),  &
        'nc_write_model_atts', 'latRt MemoryOrder '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, '_FillValue', -999.0 ), &
        'nc_write_model_atts', 'latRt _FillValue '//trim(filename))

   ! subsurface levels common to both grids 
   call nc_check(nf90_def_var(ncFileID,name='depth', xtype=nf90_real, &
        dimids=(/ depthDimID /), varid=VarID),&
        'nc_write_model_atts', 'depth def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate soil levels'), &
        'nc_write_model_atts', 'depth long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'm'),  &
        'nc_write_model_atts', 'depth units '//trim(filename))

   ! channel routing grid "links" - only output the indices relative to the
   ! full routing grid. 
   ! x-index
   call nc_check(nf90_def_var(ncFileID,name='linkIndX', xtype=nf90_int, &
        dimids=(/ linksDimID /), varid=VarID),&
        'nc_write_model_atts', 'links X def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'links x-index on routing grid '), &
        'nc_write_model_atts', 'links long_name '//trim(filename))
   ! y-index
   call nc_check(nf90_def_var(ncFileID,name='linkIndY', xtype=nf90_int, &
        dimids=(/ linksDimID /), varid=VarID),&
        'nc_write_model_atts', 'links Y def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'links y-index on routing grid '), &
        'nc_write_model_atts', 'links long_name '//trim(filename))

   ! basn mask index
   ! fix - this is conditional on gw being active and or in the DART state vector.
   call nc_check(nf90_def_var(ncFileID,name='basns', xtype=nf90_int, &
        dimids=(/ basnsDimID /), varid=VarID),&
        'nc_write_model_atts', 'basns def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'basin mask index '), &
        'nc_write_model_atts', 'basin mask index long_name '//trim(filename))

   ! scalarInd - seems unnecessary... but at least it's consistent and low-cost.
   if(assimOnly_active) then 
      call nc_check(nf90_def_var(ncFileID,name='scalarInd', xtype=nf90_int, &
           dimids=(/ scalarDimID /), varid=VarID),&
           'nc_write_model_atts', 'scalarInd def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'scalar index'), &
           'nc_write_model_atts', 'scalar index long_name '//trim(filename))
   endif

   !----------------------------------------------------------
   ! Create the (empty) state variables and their attributes
   !----------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(ivar, ncFileID, MemberDimID, unlimitedDimID, myndims, mydimids)

      ! define the variable and set the attributes

!print*,trim(varname)
!print*,mydimids(1:myndims)

      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), xtype=progvar(ivar)%xtype, &
           dimids = mydimids(1:myndims), varid=VarID),&
           'nc_write_model_atts', trim(string1)//' def_var' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', trim(progvar(ivar)%long_name)), &
           'nc_write_model_atts', trim(string1)//' put_att long_name' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'DART_kind', trim(progvar(ivar)%kind_string)), &
           'nc_write_model_atts', trim(string1)//' put_att dart_kind' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'units', trim(progvar(ivar)%units)), &
           'nc_write_model_atts', trim(string1)//' put_att units' )

      ! Preserve the original missing_value/_FillValue code.
      !     if (  progvar(ivar)%xtype == NF90_INT ) then
      !        call nc_check(nf90_put_att(ncFileID, VarID, 'missing_value', progvar(ivar)%spvalINT), &
      !             'nc_write_model_atts', trim(string1)//' put_att missing_value' )
      !        call nc_check(nf90_put_att(ncFileID, VarID, '_FillValue',  progvar(ivar)%spvalINT), &
      !             'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

      !     elseif (  progvar(ivar)%xtype == NF90_FLOAT ) then
      !        call nc_check(nf90_put_att(ncFileID, VarID, 'missing_value', progvar(ivar)%spvalR4), &
      !             'nc_write_model_atts', trim(string1)//' put_att missing_value' )
      !        call nc_check(nf90_put_att(ncFileID, VarID, '_FillValue',  progvar(ivar)%spvalR4), &
      !             'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

      !     elseif (  progvar(ivar)%xtype == NF90_DOUBLE ) then
      !        call nc_check(nf90_put_att(ncFileID, VarID, 'missing_value', progvar(ivar)%spvalR8), &
      !             'nc_write_model_atts', trim(string1)//' put_att missing_value' )
      !        call nc_check(nf90_put_att(ncFileID, VarID, '_FillValue',  progvar(ivar)%spvalR8), &
      !             'nc_write_model_atts', trim(string1)//' put_att _FillValue' )
      !     endif

   enddo

   !----------------------------------------------------------
   ! Finished with dimension/variable definitions, must end 'define' mode to fill.
   !----------------------------------------------------------

   call nc_check(nf90_enddef(ncfileID), 'nc_write_model_atts', 'prognostic enddef')

   !----------------------------------------------------------
   ! Fill the coordinate variables
   !----------------------------------------------------------
   ! coarse grid lon
   call nc_check(nf90_inq_varid(ncFileID, 'lon', VarID), &
        'nc_write_model_atts', 'inq_varid lon '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, xlong ), &
        'nc_write_model_atts', 'put_var lon '//trim(filename))
   ! coarse grid lat
   call nc_check(nf90_inq_varid(ncFileID, 'lat', VarID), &
        'nc_write_model_atts', 'inq_varid lat '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, xlat ), &
        'nc_write_model_atts', 'put_var lat '//trim(filename))

   ! fine grid lon
   call nc_check(nf90_inq_varid(ncFileID, 'lonRt', VarID), &
        'nc_write_model_atts', 'inq_varid lonRt '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, hlong ), &
        'nc_write_model_atts', 'put_var lon '//trim(filename))
   ! fine grid lat
   call nc_check(nf90_inq_varid(ncFileID, 'latRt', VarID), &
        'nc_write_model_atts', 'inq_varid latRt '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, hlat ), &
        'nc_write_model_atts', 'put_var latRt '//trim(filename))

   ! subsurface
   call nc_check(nf90_inq_varid(ncFileID, 'depth', VarID), &
        'nc_write_model_atts', 'inq_varid depth '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, zsoil8(1:nsoil)), &
        'nc_write_model_atts', 'put_var depth '//trim(filename))

   ! link indices X
   call nc_check(nf90_inq_varid(ncFileID, 'linkIndX', VarID), &
        'nc_write_model_atts', 'inq_varid linkIndX '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, channelIndsX ), &
        'nc_write_model_atts', 'put_var channelIndsX '//trim(filename))
   ! link indices Y
   call nc_check(nf90_inq_varid(ncFileID, 'linkIndY', VarID), &
        'nc_write_model_atts', 'inq_varid linkIndY '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, channelIndsY ), &
        'nc_write_model_atts', 'put_var channelIndsY '//trim(filename))

   ! basin mask index
   call nc_check(nf90_inq_varid(ncFileID, 'basns', VarID), &
        'nc_write_model_atts', 'inq_varid basns '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, basnMask ), &
        'nc_write_model_atts', 'put_var basns '//trim(filename))

   ! scalar index
   if(assimOnly_active) then 
      call nc_check(nf90_inq_varid(ncFileID, 'scalarInd', VarID), &
           'nc_write_model_atts', 'inq_varid scalarInd '//trim(filename))
      call nc_check(nf90_put_var(ncFileID, VarID, 1 ), &
           'nc_write_model_atts', 'put_var scalarInd '//trim(filename))
   endif
endif

!-------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------
if (has_noah_namelist) then
   call file_to_text(lsm_namelist_filename, lsmTextblock)
   call nc_check(nf90_put_var(ncFileID, lsmNmlVarID, lsmTextblock ), &
        'nc_write_model_atts', 'put_var lsmNmlVarID')
   deallocate(lsmTextblock)
endif

if (has_hydro_namelist) then
   call file_to_text(hydro_namelist_filename, hydroTextblock)
   call nc_check(nf90_put_var(ncFileID, hydroNmlVarID, hydroTextblock ), &
        'nc_write_model_atts', 'put_var hydroNmlVarID')
   deallocate(hydroTextblock)
endif

!-------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------
call nc_check(nf90_sync(ncFileID),'nc_write_model_atts', 'sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts
!===============================================================================


!===============================================================================
function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)
! JLM - This writes the Prior_Diag.nc and Posterior_Diag.nc files?? That should 
!       be stated. Maybe this does other things?
!
! All errors are fatal, so the return code is always '0 == normal'
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN           ! open existing netCDF dataset
!    NF90_redef       ! put into define mode
!    NF90_def_dim     ! define additional dimensions (if any)
!    NF90_def_var     ! define variables: from name, type, and dims
!    NF90_put_att     ! assign attribute values
! NF90_ENDDEF         ! end definitions: leave define mode
!    NF90_put_var     ! provide values for variable
! NF90_CLOSE          ! close: save updated netCDF dataset

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: state_vec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname
character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames
integer :: i, ivar, VarID, ncNdims, dimlen, numdims
integer :: TimeDimID, CopyDimID

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array
real(r8), allocatable, dimension(:,:,:)   :: data_3d_array

character(len=128) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!---------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!---------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!--------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
!--------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
     'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
     'nc_write_model_vars', 'inq_dimid time '//trim(filename))

if ( output_state_vector ) then

   !-----------------------------------------------------------
   ! simply blast out the vector
   !-----------------------------------------------------------

   call nc_check(nf90_inq_varid(ncFileID, 'state', VarID), &
        'nc_write_model_vars', 'inq_varid state' )
   call nc_check(nf90_put_var(ncFileID, VarID, state_vec,  &
        start=(/ 1, copyindex, timeindex /)), &
        'nc_write_model_vars', 'put_var state')

else

   !-----------------------------------------------------------
   ! We need to process the individual variables in the state vector.
   ! We already put their coordinates in above. 
   !-----------------------------------------------------------
   do ivar = 1,nfields

      varname = trim(progvar(ivar)%varname)
      string2 = trim(filename)//' '//trim(varname)

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried.
      ! This requires that Time is the unlimited dimension (the last one in Fortran),
      ! and that 'copy' is the second-to-last. The variables declared in the DART
      ! diagnostic files are required to have the same shape as in the source
      ! restart file. If Time is present there, it must also be the 'last' one.

      ! FIXME ... somewhere I should ensure that IF time is present in the original
      ! prognostic variable from the model, it is the last/unlimited dimension.

      call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
           'nc_write_model_vars', 'inq_varid '//trim(string2))

      call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
           'nc_write_model_vars', 'inquire '//trim(string2))


      mystart(:) = 1
      mycount(:) = 1
      DimCheck : do i = 1,ncNdims

         write(string1,'(A,i2,A)') 'inquire dimension ',i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), name=dimnames(i), len=dimlen), &
              'nc_write_model_vars', trim(string1))

         if (dimIDs(i) == CopyDimID) cycle DimCheck
         if (dimIDs(i) == TimeDimID) cycle DimCheck

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
            write(string2,*)' but it should be.'
            call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
                 source, revision, revdate, text2=trim(string2))
         endif

         mycount(i) = dimlen

      enddo DimCheck

      where(dimIDs == CopyDimID) mystart = copyindex
      where(dimIDs == CopyDimID) mycount = 1
      where(dimIDs == TimeDimID) mystart = timeindex
      where(dimIDs == TimeDimID) mycount = 1

      if ( (do_output()) .and. debug > 10 ) then
         write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',mystart(1:ncNdims)
         write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',mycount(1:ncNdims)
         write(*,'(A20)')'nc_write_model_vars ',dimnames(1:progvar(ivar)%numdims)
      endif


      ! this dimension count is relative to what's store in progvar
      ! the LSM variables have a extra time dimension, which we'll take off in the count.
      numdims = progvar(ivar)%numdims
      if (progvar(ivar)%component == 'LSM') numDims = numDims-1

      if ( numdims == 1 ) then

         if ( ncNdims /= 3 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 3 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                 source, revision, revdate, text2=string2)
         endif

         allocate(data_1d_array( progvar(ivar)%dimlens(1) ))
         call vector_to_prog_var(state_vec, ivar, data_1d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
              start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
              'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_1d_array)

      elseif (  numdims == 2 ) then

         if ( ncNdims /= 4 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 4 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                 source, revision, revdate, text2=string2)
         endif

         allocate(data_2d_array( progvar(ivar)%dimlens(1),  &
                                 progvar(ivar)%dimlens(2) ))
         call vector_to_prog_var(state_vec, ivar, data_2d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
              start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
              'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_2d_array)

      elseif ( numdims == 3 ) then

         if ( ncNdims /= 5 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 5 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                 source, revision, revdate, text2=string2)
         endif

         allocate(data_3d_array( progvar(ivar)%dimlens(1),  &
                                 progvar(ivar)%dimlens(2),  &
                                 progvar(ivar)%dimlens(3) ))
         call vector_to_prog_var(state_vec, ivar, data_3d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
              start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
              'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_3d_array)

      else

         write(string1,*)'do not know how to handle NOAH variables with more than 3 non-time dimensions'
         write(string2,*)trim(progvar(ivar)%varname),' has dimensions ', progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
         call error_handler(E_ERR,'nc_write_model_vars',string1, &
              source,revision,revdate,text2=string2)

      endif

   enddo

endif

!--------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars



!===============================================================================
subroutine pert_model_state(state, pert_state, interf_provided)
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.  The returned pert_state should in any
! case be valid, since it will be read by filter even if
! interf_provided is .false.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'pert_model_state', &
     'NOAH cannot be started from a single vector', &
     source, revision, revdate, &
     text2='see comments in noah/model_mod.f90::pert_model_state()',&
     text3='or noah/model_mod.html#pert_model_state')

interf_provided = .false.

end subroutine pert_model_state



!===============================================================================
subroutine ens_mean_for_model(ens_mean)
real(r8), intent(in) :: ens_mean(:)  !! this variable is module scope
if ( .not. module_initialized ) call static_init_model
end subroutine ens_mean_for_model


!==================================================================
! PUBLIC interfaces that aren't required by the DART code but are
! generally useful for other related utility programs.
! (less necessary for small models; generally used for larger models
! with predefined file formats and control structures.)
!==================================================================

function get_debug_level()
integer :: get_debug_level

get_debug_level = debug

end function get_debug_level


!===============================================================================
subroutine verify_state_variables( ncid, filename, stateVarList, ngood )
! for an ncfile handle and name, return the number of good variables.
integer,                      intent(in)  :: ncid
character(len=*),             intent(in)  :: filename
character(len=obstypelength), intent(in)  :: stateVarList(:,:)
integer,                      intent(out) :: ngood

integer :: i, VarID
character(len=NF90_MAX_NAME) :: varname, dartstr

ngood = 0
MyLoop : do i = 1, size(stateVarList,2)
   varname = trim(stateVarList(1,i))
   dartstr = trim(stateVarList(2,i))

   if ( varname == ' ' .and. dartstr == ' ' ) exit MyLoop ! Found end of list.

   if ( varname == ' ' .or.  dartstr == ' ' ) then
      string1 = 'model_nml:stateVarList not fully specified'
      string2 = 'reading from ['//trim(filename)//']'
      call error_handler(E_ERR,'verify_state_variables', string1, &
              source, revision, revdate, text2=string2)
   endif

   ! Make sure variable exists in netCDF file
   write(string1,'(''there is no variable '',a,'' in '',a)') trim(varname), trim(filename)
   call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
        'verify_state_variables', trim(string1))

   ! Make sure DART kind is valid
   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector
   if (do_output() .and. (debug > 0)) then
      if(i==1) then
         write(logfileunit,*)''
         write(     *     ,*)''
         write(logfileunit,*)'Variables in: [',trim(filename),']'
         write(     *     ,*)'Variables in: [',trim(filename),']'
      end if
      write(logfileunit,*)'variable ',i,' is ',trim(varname), '   ', trim(dartstr)
      write(     *     ,*)'variable ',i,' is ',trim(varname), '   ', trim(dartstr)
   endif
   ngood = ngood + 1
enddo MyLoop

if (ngood == MAX_STATE_VARIABLES) then
   string1 = 'WARNING: There is a possibility you need to increase ''MAX_STATE_VARIABLES'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate,text2=string2)
endif

end subroutine verify_state_variables

!===============================================================================
! TJH FIXME ... how are these used - do they need to be public?
!! why not just rename the variables *_restart_filename instead of *_netcdf_filename?
subroutine get_lsm_restart_filename( lsm_restart_filename )
! why not just rename the variables *_restart_filename instead of *_netcdf_filename?
character(len=*), intent(out) :: lsm_restart_filename
if ( .not. module_initialized ) call static_init_model
lsm_restart_filename = lsm_netcdf_filename
end subroutine get_lsm_restart_filename

!===============================================================================
! TJH FIXME ... how are these used - do they need to be public?
subroutine get_hydro_restart_filename( hydro_restart_filename )
character(len=*), intent(out) :: hydro_restart_filename
if ( .not. module_initialized ) call static_init_model
hydro_restart_filename = hydro_netcdf_filename
end subroutine get_hydro_restart_filename

!===============================================================================
subroutine get_assimOnly_restart_filename( assimOnly_restart_filename )
character(len=*), intent(out) :: assimOnly_restart_filename
if ( .not. module_initialized ) call static_init_model
assimOnly_restart_filename = assimOnly_netcdf_filename
end subroutine get_assimOnly_restart_filename


!===============================================================================
subroutine model_to_dart_vector(filenameLsm, filenameHydro, filenameAssimOnly, &
     state_vector, restart_time)
! Reads the current time and state variables from a model data
! file and packs them into a dart state vector.

character(len=*), intent(in)  :: filenameLsm, filenameHydro, filenameAssimOnly
real(r8),         intent(out) :: state_vector(:)
type(time_type),  intent(out) :: restart_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, ncstart, nccount
character(len=NF90_MAX_NAME) :: dimname, varname, filename

integer  :: ncid, ncNdims, dimlen, VarID
integer  :: ncidLsm, ncidHydro, ncidAssimOnly, ncomponents
integer  :: i, indx1, indx2, indx3, indx4, indx, ivar, ntimes, ifile

! jlm fixme 2 variables for extreme debugging
!real(r8), dimension(3)  :: dumLoc
!integer  :: ii

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array
real(r8), allocatable, dimension(:,:,:)   :: data_3d_array
real(r8), allocatable, dimension(:,:,:,:) :: data_4d_array

if ( .not. module_initialized ) call static_init_model

state_vector(:) = MISSING_R8

! Check that the input file exists ...

ncomponents=1
if (hydro_model_active) ncomponents=2
if (assimOnly_active)   ncomponents=3
if ( (.not. hydro_model_active) .and. assimOnly_active ) then
   write(string1,*) 'Not configured to have assimOnly_active while not hydro_model_active'
   call error_handler(E_ERR,'model_to_dart_vector',string1,source,revision,revdate)
endif

do ifile=1,ncomponents

   filename = filenameLsm
   if (ifile==2) filename = filenameHydro
   if (ifile==3) filename = filenameAssimOnly

   if ( .not. file_exist(filename) ) then
      write(string1,*) 'file <', trim(filename),'> does not exist.'
      call error_handler(E_ERR,'model_to_dart_vector',string1,source,revision,revdate)
   endif

   call nc_check(nf90_open(adjustl(filename), NF90_NOWRITE, ncid), &
        'model_to_dart_vector', 'open '//trim(filename))

   if (ifile==1) restart_time = get_state_time(ncid, trim(filename))

   if ( (do_output()) .and. debug > 99 ) then
      call print_date(restart_time,'model_to_dart_vector:date of restart file '//trim(filename))
      call print_time(restart_time,'model_to_dart_vector:time of restart file '//trim(filename))
   endif

   if (ifile==1) ncidLsm = ncid
   if (ifile==2) ncidHydro = ncid
   if (ifile==3) ncidAssimOnly = ncid

end do

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.

do ivar=1, nfields  !! jlm - going to need nfieldsLsm and nfieldsHydro

   ntimes     = -1
   ncstart(:) = -1
   nccount(:) = -1
   varname    = trim(progvar(ivar)%varname)
   string3    = trim(filename)//' '//trim(varname)

   !! if statement for LSM file vs hydro file which assigns
   ncid = ncidLsm
   if (progvar(ivar)%component == 'HYDRO') ncid = ncidHydro
   if (progvar(ivar)%component == 'assimOnly') ncid = ncidAssimOnly

   call nc_check(nf90_inq_varid(ncid, varname, VarID), &
        'model_to_dart_vector', 'inq_varid '//trim(string3))
   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
        'model_to_dart_vector', 'inquire '//trim(string3))

   ! Check the shape of the variable
   if ( ncNdims /= progvar(ivar)%numdims ) then
      write(string1, *) 'netCDF rank of '//trim(varname)//' does not agree with internal rank.'
      write(string2, *) 'netCDF rank is ',ncNdims,' expected ',progvar(ivar)%numdims
      write(string3, *) 'should not happen'
      call error_handler(E_ERR,'model_to_dart_vector', string1, &
           source,revision,revdate,text2=string2,text3=string3)
   endif

   ! Check the memory order of the variable
   ! making sure we only compare the last timestep ...
   do i = 1,ncNdims
      write(string1,'(''inquire dimension'',i2,A)') i,trim(string3)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimlen), &
           'model_to_dart_vector',string1)

      if (trim(dimname) /= trim(progvar(ivar)%dimnames(i))) then
         write(string1, *) 'netCDF dimnames of '//trim(varname)//' does not match expected dimname'
         write(string2, *) 'netCDF dimname is ',trim(dimname),' expected ',trim(progvar(ivar)%dimnames(i))
         write(string3, *) 'should not happen'
         call error_handler(E_ERR,'model_to_dart_vector', string1, &
              source,revision,revdate,text2=string2,text3=string3)
      endif

      ncstart(i) = 1
      if (trim(varname) .eq. 'sh2ox')     dimlen = fine3dShape(i)
      if (trim(varname) .eq. 'sfcheadrt') dimlen = fine2dShape(i)
      nccount(i) = dimlen

      if ( trim(dimname) == 'Time' ) then
         ntimes     = dimlen
         ncstart(i) = dimlen
         nccount(i) = 1
      elseif ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string3),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'model_to_dart_vector',string1,source,revision,revdate)
      endif
   enddo

   if ( (do_output()) .and. debug > 99 ) then
      write(*,*)
      write(*,*)'model_to_dart_vector: component ',trim(progvar(ivar)%component)
      write(*,*)'model_to_dart_vector: variable ' ,trim(varname)
      write(*,*)'model_to_dart_vector: ncstart '  ,ncstart(1:ncNdims)
      write(*,*)'model_to_dart_vector: nccount '  ,nccount(1:ncNdims)
      write(*,*)'model_to_dart_vector: state_vector start:end ',&
           progvar(ivar)%index1, progvar(ivar)%indexN
      write(*,*)
   endif

   ! FIXME - this is probably the place to ensure that the Time dimension is the last
   ! dimension if it is present. unlimited dimension

   !-----------------------------------
   ! Pack the variable into the DART state vector
   indx = progvar(ivar)%index1

   if (ncNdims == 1) then

      allocate(data_1d_array(nccount(1)))

      call nc_check(nf90_get_var(ncid, VarID, data_1d_array,  &
           start=ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
           'model_to_dart_vector', 'get_var '//trim(string3))

      do indx1 = 1, nccount(1)
         state_vector(indx) = data_1d_array(indx1)
         indx = indx + 1
      enddo
      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      allocate(data_2d_array(nccount(1), nccount(2)))

      if (trim(varname) .eq. 'sfcheadrt') then 
         !! disags sfcheadrt  this is in the wrong dimension
         write(logfileunit,*)"Disaggregating HYDRO surface head (sfcheadrt)."
         write(     *     ,*)"Disaggregating HYDRO surface head (sfcheadrt)."
         data_2d_array = sfcHeadDisag
      else
         call nc_check(nf90_get_var(ncid, VarID, data_2d_array,  &
              start=ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
              'model_to_dart_vector', 'get_var '//trim(string3))
      end if

      do indx2 = 1, nccount(2)
      do indx1 = 1, nccount(1)
         state_vector(indx) = data_2d_array(indx1,indx2)
         indx = indx + 1
      enddo
      enddo
      deallocate(data_2d_array)

   elseif (ncNdims == 3) then

      allocate(data_3d_array(nccount(1), nccount(2), nccount(3)))

      if (trim(varname) .eq. 'sh2ox') then 
         !! disags liquid water content (sh2oRt) into module memory for use elsewhere.
         write(logfileunit,*)"Disaggregating HYDRO liquid soil moisture."
         write(     *     ,*)"Disaggregating HYDRO liquid soil moisture."
         data_3d_array = sh2oDisag
      else 
         call nc_check(nf90_get_var(ncid, VarID, data_3d_array,  &
              start=ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
              'model_to_dart_vector', 'get_var '//trim(string3))
      end if

      do indx3 = 1, nccount(3)
      do indx2 = 1, nccount(2)
      do indx1 = 1, nccount(1)
         state_vector(indx) = data_3d_array(indx1,indx2,indx3)
         indx = indx + 1
      enddo
      enddo
      enddo
      deallocate(data_3d_array)

   elseif (ncNdims == 4) then

      allocate(data_4d_array(nccount(1), &
           nccount(2), &
           nccount(3), &
           nccount(4)))

      call nc_check(nf90_get_var(ncid, VarID, data_4d_array,  &
           start=ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
           'model_to_dart_vector', 'get_var '//trim(string3))

      ! TJH RAFAEL transform goes here.

      do indx4 = 1, nccount(4)
      do indx3 = 1, nccount(3)
      do indx2 = 1, nccount(2)
      do indx1 = 1, nccount(1)
         state_vector(indx) = data_4d_array(indx1,indx2,indx3,indx4)
         indx = indx + 1
      enddo
      enddo
      enddo
      enddo
      deallocate(data_4d_array)

   else

      write(string1, *)'Variable '//trim(varname)//' has ',ncNdims,' dimensions.'
      write(string2, *)'cannot handle that.'
      call error_handler(E_ERR,'model_to_dart_vector', string1, &
           source,revision,revdate,text2=string2)

   endif

! jlm fixme - extreme debugging
!do ii = progvar(ivar)%index1, progvar(ivar)%indexN 
!   dumLoc = get_location( state_loc(ii) )
!   state_vector(ii) = dumLoc(2)  ! lon, lat, ele
!end do
! end extreme debugging

   indx = indx - 1
   if ( indx /= progvar(ivar)%indexN ) then
      write(string1, *)'Variable '//trim(varname)//' filled wrong.'
      write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',indx
      call error_handler(E_ERR,'model_to_dart_vector', string1, &
           source,revision,revdate,text2=string2)
   endif

enddo

if ( (do_output()) .and. debug > 100 ) then
   write(*,*)'newest time is ',ntimes
   do i = 1,size(state_vector)
      write(*,*)'state vector(',i,') is',state_vector(i)
   enddo
endif

!!! jlm fix ? where is/are the file/s opened above closed??

end subroutine model_to_dart_vector


!===============================================================================-------------
subroutine disagHydro()
! There are two hydro quantities which are disaggregated. Soil moisture content (SMC) and 
! surface head (sfchead). sfchead depends on smc in the disaggregation, so they are 
! done together if sfchead is present. If sfchead is not present, then only sh2ort is disaggregated. 
! This is possible since their re-aggregation does not have dependence. 
!
! The hydro code used here starts at (roughly) line 2327 in Routing/Noah_distr_routing.F
! Taken from HYDRO_drv/module_HYDRO_drv.F near line 536
! Several changes in naming conventions are mentioned below. 
!
! --- SMC ---
! Disaggregating of sh2ox can be done simply using weights.
! HOWEVER, we would like to physically constrain perturbations to sh2ox using
! the wilt and max water holding of the soil and the amount of ice in a given layer (sice=smc-sh2ox). 
! total soil moisture coarse and fine res:    smc, smcrt
! liquid soil moisture coarse and fine res: sh2ox, sh2ort
! low res ice: sice=smc-sh2ox
! coarse res, 2D, total soil moisture wilt and max soil water: smcWlt1, smcMax1
! fine res, 3D, *liquid* soil moisture wilt & max soil water: sh2oWltRt, sh2oxMaxRt (adj for ice)
! fine res weights for disag: sh2owgts
! NOTE I renamed: smcrt->sh2oRt, smcWltRt->sh2oWltRt, and smcMaxRt->sh2oMaxRt compared to hydro code.

! --- SFCHEAD ---
! In this routine 
! sfcHeadSubRt - fine resolution surface head
! sfcHeadRt    - coarse resolution surface head
! infxswgts    - fine resolution weights to go between the above surface head. 
! In the hydro code infxsrt is disaggregated to infxsubrt. 
! This is because 
! sfcheadrt -> LSM -> infxsrt -> HYDRO -> infxsubrt -> sfcheadsubrt -> infxswgt -> sfcheadrt
! So that infiltration excess is surface head after it's been allowed to infiltrate.
! sfcheadrt and infxswgt are the prognostic variables saved to restart file at the end of 
! the advance. We do the following disaggregation here, outside the model advance:
! sfcheadrt, infxwgt -> sfcheadsubrt
! Then sfcheadsubrt is adjusted in the assimilation. 
! (The aggregation routine reverses this: sfcheadsubrt -> sfcheadrt, infxwgt)

! IN/OUT
! Done by setting global variables.
! These are set by the disag
!real(R8), dimension(:,:,:), allocatable :: sh2oDisag
!real(R8), dimension(:,:),   allocatable :: sfcHeadDisag

! global 
! real(R8), dimension(:,:,:), allocatable :: sice  !! global
! real(R8), dimension(:,:,:), ALLOCATABLE :: sh2oMaxRt, sh2oWltRt !! global
! real(R8), dimension(west_east,south_north)   :: smcWlt1, smcMax1

! read from hydro restart file.
real(R8), dimension(west_east,south_north,nsoil) :: sh2ox
real(R8), dimension(n_hlong,n_hlat,nsoil)        :: sh2oWgt,   sh2ort
real(R8), dimension(:,:), allocatable            :: sfcHeadRt, infxsWgt, sfcHeadSubRt

! local
real(R8)           :: LSMVOL, SMCEXCS, WATHOLDCAP, area_lsm
character(len=250) :: errString
integer            :: ncid, varId, i, j, ix, iy, jx, ixrt, iyrt
integer            :: krt, ixxrt, jyyrt, ncNdims, kf
integer            :: AGGFACXRT, AGGFACYRT, AGGFACTRT
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
real(r8), dimension(nsoldx) :: SLDPTH

!! this code starts at (roughly) line 2327 in Routing/Noah_distr_routing.F
iy = south_north
jx = iy
ix = west_east
ixrt = n_hlong
iyrt = n_hlat
!nsoil=nsoil
aggfactrt = n_hlong/west_east
sldpth = soil_thick_input
area_lsm = coarseGridArea

!! global scope variables, so be more careful
!! allocated here but deallocated in aggSh2ox. These should not be allocated when arriving here...
if (.not. allocated(smcWlt1))    allocate(  smcWlt1(west_east,south_north))
if (.not. allocated(smcMax1))    allocate(  smcMax1(west_east,south_north))
if (.not. allocated(sh2oWltRt))  allocate(sh2oWltRt(ixrt,iyrt,      nsoil))
if (.not. allocated(sh2oMaxRt))  allocate(sh2oMaxRt(ixrt,iyrt,      nsoil))
if (.not. allocated(smc)    )    allocate(      smc(ix,   iy,       nsoil))
if (.not. allocated(sice)    )   allocate(     sice(ix,   iy,       nsoil))

!! open the hydro restart file
call nc_check(nf90_open(adjustl(hydro_netcdf_filename), NF90_NOWRITE, ncid), &
     'disagHydro', 'open '//trim(hydro_netcdf_filename))

!! SMC
errString    = 'disagHydro: smc'
call nc_check(nf90_inq_varid(ncid, 'smc', varId), 'varId: '//trim(errString))
call nc_check(nf90_get_var(ncid, VarID, smc), 'disagHydro', 'get_var: '//trim(errString))

!! sh2ox
errString    = 'disagHydro: sh2ox'
call nc_check(nf90_inq_varid(ncid, 'sh2ox', varId), 'varId: '//trim(errString))
call nc_check(nf90_get_var(ncid, VarID, sh2ox),'disagHydro', 'get_var: '//trim(errString))

!! smcWlt1
errString    = 'disagHydro: smcWlt1'
call nc_check(nf90_inq_varid(ncid, 'smcwlt1', varId), 'varId: '//trim(errString))
call nc_check(nf90_get_var(ncid, VarID, smcWlt1),'disagHydro', 'get_var: '//trim(errString))

!! smcMax1
errString    = 'disagHydro: smcmax1'
call nc_check(nf90_inq_varid(ncid, 'smcmax1', varId), 'varId: '//trim(errString))
call nc_check(nf90_get_var(ncid, VarID, smcMax1),'disagHydro', 'get_var: '//trim(errString))

!! sh2owgt
errString    = 'disagHydro: sh2owgt'
call nc_check(nf90_inq_varid(ncid, 'sh2owgt', varId), 'varId: '//trim(errString))
call nc_check(nf90_get_var(ncid, VarID, sh2owgt),'disagHydro', 'get_var: '//trim(errString))

if (hydroSfcHeadPresent) then 
   allocate(sfcHeadRt(west_east,south_north), &
            infxsWgt(n_hlong,n_hlat), sfcHeadSubRt(n_hlong,n_hlat) )

   !! sfcHeadRt
   errString    = 'disagHydro: sfcHeadRt'
   call nc_check(nf90_inq_varid(ncid, 'sfcheadrt', varId), 'varId: '//trim(errString))
   call nc_check(nf90_get_var(ncid, VarID, sfcHeadRt), 'disagHydro', 'get_var: '//trim(errString))

   !! infxswgt
   errString    = 'disagHydro: infxsWgt'
   call nc_check(nf90_inq_varid(ncid, 'infxswgt', varId), 'varId: '//trim(errString))
   call nc_check(nf90_get_var(ncid, VarID, infxsWgt),'disagHydro', 'get_var: '//trim(errString))
end if

!! close the hydro restart file
call nc_check(nf90_close(ncid), 'disagHydro', 'close '//trim(hydro_netcdf_filename))

!! calculate sice
sice=smc-sh2ox
!! scmrefrt not needed unless we want to adjust the disaggregated lksat.

!! disag
do J=1,JX  !! also know as y
   do I=1,IX

      if (hydroSfcHeadPresent) &
           !LSMVOL=INFXSRT(I,J)*area_lsm(I,J)
           LSMVOL=sfcHeadRt(I,J)*area_lsm

      do AGGFACYRT=AGGFACTRT-1,0,-1
         do AGGFACXRT=AGGFACTRT-1,0,-1

            IXXRT=I*AGGFACTRT-AGGFACXRT
            JYYRT=J*AGGFACTRT-AGGFACYRT

            if (hydroSfcHeadPresent) &
                 !INFXSUBRT(IXXRT,JYYRT)=LSMVOL*INFXSWGT(IXXRT,JYYRT)/fineGridArea
                 sfcHeadSubRt(IXXRT,JYYRT)=LSMVOL*INFXSWGT(IXXRT,JYYRT)/fineGridArea

            !! note that this block could be moved out of the agg do loop as nothing
            !! depends on high-res indices and the assingments to high-res can be vectorized.
            !! In other words, ice is uniform at coarse resolution.
            do KRT=1,NSOIL  !Do for soil profile loop
               if(SICE(I,J,KRT).gt.0) then  !...adjust for soil ice
                  !DJG Adjust SMCMAX for SICE when subsfc routing...make 3d variable
                  sh2oMaxRt(IXXRT,JYYRT,KRT)=SMCMAX1(I,J)-SICE(I,J,KRT)
                  WATHOLDCAP = SMCMAX1(I,J) - SMCWLT1(I,J)
                  if (SICE(I,J,KRT) .le. WATHOLDCAP)    then
                     sh2oWltRt(IXXRT,JYYRT,KRT) = SMCWLT1(I,J)
                  else
                     if(SICE(I,J,KRT).lt.SMCMAX1(I,J)) &
                          sh2oWltRt(IXXRT,JYYRT,KRT) = SMCWLT1(I,J) - &
                          (SICE(I,J,KRT)-WATHOLDCAP)
                     if(SICE(I,J,KRT).ge.SMCMAX1(I,J)) sh2oWltRt(IXXRT,JYYRT,KRT) = 0.
                  end if
               else
                  sh2oMaxRt(IXXRT,JYYRT,KRT)= SMCMAX1(I,J)
                  sh2oWltRt(IXXRT,JYYRT,KRT) = SMCWLT1(I,J)
                  ! watholdcap not used after here... 
               end if   !endif adjust for soil ice...

               !Now Adjust soil moisture
               !This block does depend on high-res grid.
               if(sh2oMaxRt(IXXRT,JYYRT,KRT).gt.0) then !Check for smcmax data (=0 over water)
                  SH2ORT(IXXRT,JYYRT,KRT)=SH2OX(I,J,KRT)*SH2OWGT(IXXRT,JYYRT,KRT)
               else
                  SH2ORT(IXXRT,JYYRT,KRT) = 0.001  !will be skipped w/ landmask
                  sh2oMaxRt(IXXRT,JYYRT,KRT) = 0.001
               end if

               if (hydroSfcHeadPresent) then
                  !DJG Check/Adjust so that subgrid cells do not exceed saturation...
                  if (sh2oRt(IXXRT,JYYRT,KRT).gt.sh2oMaxRt(IXXRT,JYYRT,KRT)) then
                     SMCEXCS = (sh2oRt(IXXRT,JYYRT,KRT) - sh2oMaxRt(IXXRT,JYYRT,KRT)) &
                          * SLDPTH(KRT)*1000.  !Excess soil water in units of (mm)
                     sh2oRt(IXXRT,JYYRT,KRT) = sh2oMaxRt(IXXRT,JYYRT,KRT)
                     do KF = KRT-1,1, -1  !loop back upward to redistribute excess water from disagg.
                        sh2oRt(IXXRT,JYYRT,KF) = sh2oRt(IXXRT,JYYRT,KF) + SMCEXCS/(SLDPTH(KF)*1000.) 
                        if (sh2oRt(IXXRT,JYYRT,KF).gt.sh2oMaxRt(IXXRT,JYYRT,KF)) then  !Recheck new lyr sat.
                           SMCEXCS = (sh2oRt(IXXRT,JYYRT,KF) - sh2oMaxRt(IXXRT,JYYRT,KF)) &
                                * SLDPTH(KF)*1000.  !Excess soil water in units of (mm)
                           sh2oRt(IXXRT,JYYRT,KF) = sh2oMaxRt(IXXRT,JYYRT,KF)
                        else  ! Excess soil water expired
                           SMCEXCS = 0.
                           exit
                        end if
                     end do
                     if (SMCEXCS.gt.0) then  !If not expired by sfc then add to Infil. Excess
                        !INFXSUBRT(IXXRT,JYYRT) = INFXSUBRT(IXXRT,JYYRT) + SMCEXCS
                        sfcHeadSubRt(IXXRT,JYYRT) = sfcHeadSubRt(IXXRT,JYYRT) + SMCEXCS
                        SMCEXCS = 0.
                     end if
                  end if  !End if for soil moisture saturation excess
               end if
            end do !KRT

            do KRT=1,NSOIL  !debug loop
               if (sh2oRt(IXXRT,JYYRT,KRT).gt.sh2oMaxRt(IXXRT,JYYRT,KRT)) then
                  string1 = &
                       "SMCMAX exceeded upon disagg. Inds & values: ixxrt,jyyrt,krt,sh2oRt,sh2oMaxRt"
                  write(string2,*) ixxrt,jyyrt,krt, sh2oRt(IXXRT,JYYRT,KRT),sh2oMaxRt(IXXRT,JYYRT,KRT)
                  call error_handler(E_ERR, 'disaggregateHydro', string1, &
                                     source, revision, revdate, text2=string2)
               else if (sh2oRt(IXXRT,JYYRT,KRT).le.0.) then
                  string1 = &
                       "SMCRT depleted on disag. Inds &values:ixxrt,jyyrt,krt,sh2oRt,sh2oMaxRt,sh2ox"
                  write(string2,*) ixxrt,jyyrt,krt, &
                       sh2oRt(IXXRT,JYYRT,KRT),sh2oMaxRt(IXXRT,JYYRT,KRT),sh2ox(I,J,KRT)
                  write(string3,*) "i,j,krt, nsoil",i,j,krt,nsoil
                  call error_handler(E_ERR, 'disaggregateHydro', string1, &
                                     source, revision, revdate, text2=string2, text3=string3)
               end if
            end do !debug loop

         end do !AGGFACXRT
      end do !AGGFACYRT
   end do !IX
end do !JX

!! this is not vegas? (what happens here stays here or not).
if (hydroSmcPresent)     sh2oDisag    = sh2oRt
if (hydroSfcHeadPresent) sfcHeadDisag = sfcHeadSubRt

if (hydroSfcHeadPresent) then 
   deallocate(sfcHeadRt, infxsWgt, sfcHeadSubRt)
endif

end subroutine disagHydro


!===============================================================================
subroutine dart_vector_to_model_files(state_vector, &
     filenameLsm, filenameHydro, filenameAssimOnly, &
     dart_time, skip_variables)
! Writes the current time and state variables from a dart state
! vector (1d array) into a noah netcdf restart file.
!
! This is VERY similar to nc_write_model_vars() for this model.
! If it were not for the 'copy' dimension, it would be identical, I think.

! this should print message about variables being skipped.

real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filenameLsm, filenameHydro, filenameAssimOnly
type(time_type),  intent(in) :: dart_time
character(len=*), intent(in) :: skip_variables(:)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname, component, filename
character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: ncFileID, ncidLsm, ncidHydro, ncidAssimOnly
integer :: VarID, ncNdims, TimeDimID
integer :: timeindex, dimlen, numdims, timedimcounter, idum

type(time_type) :: lsm_time, hydro_time

! temp space to hold data while we are writing it
integer :: i, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array
real(r8), allocatable, dimension(:,:,:)     :: smcAgg, smcAggWeights
real(r8), allocatable, dimension(:,:)       :: sfcHeadAgg, infxsWeights

! variables related to screwing around with soil moisture states
integer :: varidOrig, ncNdimsOrig
integer, dimension(NF90_MAX_VAR_DIMS) :: mycountSh2ox, mystartSh2ox, mycountSfcHead, mystartSfcHead

if ( .not. module_initialized ) call static_init_model

! Check that the output file exists ...
! LSM
if ( .not. file_exist(filenameLsm) ) then
   write(string1,*) 'cannot open file ', trim(filenameLsm),' for writing.'
   call error_handler(E_ERR,'dart_vector_to_model_files',string1,source,revision,revdate)
endif
call nc_check(nf90_open(trim(filenameLsm), NF90_WRITE, ncidLsm), &
     'dart_vector_to_model_files','open '//trim(filenameLsm))

! Hydro
if ( .not. file_exist(filenameHydro) ) then
   write(string1,*) 'cannot open file ', trim(filenameHydro),' for writing.'
   call error_handler(E_ERR,'dart_vector_to_model_files',string1,source,revision,revdate)
endif
call nc_check(nf90_open(trim(filenameHydro), NF90_WRITE, ncidHydro), &
     'dart_vector_to_model_files','open '//trim(filenameHydro))

! AssimOnly
if(assimOnly_active) then
   if ( .not. file_exist(filenameAssimOnly) ) then
      write(string1,*) 'cannot open file ', trim(filenameAssimOnly),' for writing.'
      call error_handler(E_ERR,'dart_vector_to_model_files',string1,source,revision,revdate)
   endif
   call nc_check(nf90_open(trim(filenameAssimOnly), NF90_WRITE, ncidAssimOnly), &
        'dart_vector_to_model_files','open '//trim(filenameAssimOnly))
endif

! none of the returned quantities are used in this routine.
! comment out for lsm and see if it breaks anything...
!call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
!                  'dart_vector_to_model_files', 'inquire '//trim(filename))

! There is no time dimension in the hydro restarts. so I skip this for hydro
call nc_check(nf90_inq_dimid(ncidLsm, 'Time', TimeDimID), &
     'dart_vector_to_model_files','inq_dimid Time '//trim(filenameLsm))

! make sure the time in the file is the same as the time on the data
! we are trying to insert.  we are only updating part of the contents
! of the NOAH restart file, and state vector contents from a different
! time won't be consistent with the rest of the file.

lsm_time = get_state_time(ncidLsm, trim(filenameLsm), timeindex)
! needs implimented for the new hydro global att restart time.
hydro_time = lsm_time  !! dummy for now
!hydro_time = get_state_time(ncidHydro, trim(filenameHydro), timeindex)
!! jlm fix - can now get the time for hydro and should also have the 
!! time in the restart.assimOnly.nc as well. 

if ( lsm_time /= dart_time ) then 
   if ( lsm_time /= dart_time ) then
      call print_time(dart_time,'DART current time',logfileunit)
      call print_time(lsm_time,'LSM current time',logfileunit)
      call print_time(dart_time,'DART current time')
      call print_time(lsm_time,'LSM current time')
   end if
   if ( hydro_time /= dart_time ) then
      call print_time(dart_time,'DART current time',logfileunit)
      call print_time(hydro_time,'HYDRO current time',logfileunit)
      call print_time(dart_time,'DART current time')
      call print_time(hydro_time,'HYDRO current time')
   end if
   write(string1,*),' current time /= model time. FATAL error.'
   call error_handler(E_ERR,'dart_vector_to_model_files',string1,source,revision,revdate)
endif

if (do_output() .and. (debug > 0)) then
   call print_date(lsm_time,'dart_vector_to_model_files: date of restart files '&
        //trim(filenameLsm)//' & '//trim(filenameHydro))
   call print_time(lsm_time,'dart_vector_to_model_files: time of restart files '&
        //trim(filenameLsm)//' & '//trim(filenameHydro))
endif

! The DART prognostic variables are only defined for a single time.
! IF the netCDF variable has a TIME dimension, it must be the last dimension.

UPDATE : do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   component = trim(progvar(ivar)%component)
   filename = filenameLsm
   if (component == 'HYDRO') filename = filenameHydro
   if (component == 'assimOnly') filename = filenameAssimOnly
   !! set the ncid based on the component
   ncFileId = ncidLsm
   if (trim(component) == 'HYDRO') ncFileID = ncidHydro
   if (trim(component) == 'assimOnly') ncFileID = ncidAssimOnly
   string2 = trim(filename)//' '//trim(varname)  ! for diagnostics

   ! If this variable is on the skip list ... skip it.
   SKIPME : do i = 1,size(skip_variables)
      if (len_trim(skip_variables(i)) < 1) cycle SKIPME
      if (skip_variables(i) == varname) cycle UPDATE
   enddo SKIPME

   ! Ensure netCDF variable is conformable with DART progvar quantity.
   call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
        'dart_vector_to_model_files', 'inq_varid '//trim(string2))
   call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
        'dart_vector_to_model_files', 'inquire '//trim(string2))

   timedimcounter = 0
   mystart(:) = 1
   mycount(:) = 1
   DimCheck : do i = 1,progvar(ivar)%numdims
      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), name=dimnames(i), len=dimlen), &
           'dart_vector_to_model_files', string1)

      !! These vars wont check because they are different in file than used in DART.
      if (trim(varname) .eq. 'sh2ox') dimlen = fine3dShape(i) !! completely cheating
      if (trim(varname) .eq. 'sfcheadrt') dimlen = fine2dShape(i) !! completely cheating

      if (component == 'LSM' .and. dimIDs(i) == TimeDimID) timedimcounter = 1
      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         write(string2,*)' but it should be.'
         call error_handler(E_ERR, 'dart_vector_to_model_files', string1, &
              source, revision, revdate, text2=string2)
      endif
      mycount(i) = dimlen
   enddo DimCheck

   if (component == 'LSM' .and. dimIDs(ncNdims) /= TimeDimID) then
      write(string1,*) trim(string2),' required to have "Time" as the last/unlimited dimension'
      write(string2,*)' last dimension is ',trim(dimnames(ncNdims))
      call error_handler(E_ERR, 'dart_vector_to_model_files', string1, &
                         source, revision, revdate, text2=string2)
      where(dimIDs == TimeDimID) mystart = timeindex
      where(dimIDs == TimeDimID) mycount = 1
   endif

   if ( (do_output()) .and. debug > 99 ) then
      write(*,*)
      write(*,*)trim(varname)
      write(*,*)'dart_vector_to_model_files '//trim(varname)//' start is ',mystart(1:ncNdims)
      write(*,*)'dart_vector_to_model_files '//trim(varname)//' count is ',mycount(1:ncNdims)
      do idum=1,progvar(ivar)%numdims
         write(*,*)'dart_vector_to_model_files ',trim(dimnames(idum))
      end do
   endif

   numdims = progvar(ivar)%numdims - timedimcounter

   if ( numdims == 1 ) then

      allocate(data_1d_array(progvar(ivar)%dimlens(1)))
      call vector_to_prog_var(state_vector, ivar, data_1d_array,limit=.true.)
      call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
           start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
           'dart_vector_to_model_files', 'put_var '//trim(string2))
      deallocate(data_1d_array)

   elseif ( numdims == 2 ) then

      allocate(data_2d_array(progvar(ivar)%dimlens(1), &
                             progvar(ivar)%dimlens(2)) )
      call vector_to_prog_var(state_vector, ivar, data_2d_array,limit=.true.)

      !! I really need to do these in a better way, just dont have time right now.
      if (trim(varname) .eq. 'sfcheadrt') then
         allocate(sfcHeadAgg(coarse2dShape(1),coarse2dShape(2)), &
                  infxsWeights(fine2dShape(1),fine2dShape(2))     )
         sfcHeadDisag = data_2d_array ! global, but declared yet??
         deallocate(data_2d_array)
         call aggSfcHead(sfcHeadAgg, infxsWeights)
         allocate(data_2d_array(coarse2dShape(1),coarse2dShape(2)))
         !! sfcHeadRt onfile = sfcHeadAgg
         data_2d_array = sfcHeadAgg
         deallocate(sfcHeadAgg)
         varidOrig = varid !! this will be written to file after this if statment
         ncNdimsOrig = ncNdims
         
         !! output sh2owgts to the ncid file
         !! i'm kinda cheating by not doing all the dimension checks...
         varname = 'infxswgt'
         mycount(1:2)=fine2dShape(1:2)
         string2 = trim(filename)//' '//trim(varname)  ! for diagnostics
         call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
              'dart_vector_to_model_files', 'inq_varid '//trim(string2))
         call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
              'dart_vector_to_model_files', 'inquire '//trim(string2))
         call nc_check(nf90_put_var(ncFileID, VarID, infxsWeights, &
              start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
              'dart_vector_to_model_files', 'put_var '//trim(string2))
         deallocate(infxsWeights)
         ! Make note that the variable has been updated by DART
         call nc_check(nf90_Redef(ncFileID),'dart_vector_to_model_files', 'redef '//trim(filename))
         call nc_check(nf90_put_att(ncFileID, VarID,'DART','variable modified by DART'),&
              'dart_vector_to_model_files', 'modified '//trim(varname))
         call nc_check(nf90_enddef(ncfileID),'dart_vector_to_model_files','state enddef '//trim(filename))

         !! reset the current variable to sfcheadrt
         varname = 'sfcheadrt'
         string2 = trim(filename)//' '//trim(varname)  ! for diagnostics
         varid = varidOrig
         ncNdims = ncNdimsOrig
         mycount(1:2) = coarse2dShape(1:2)
         print*,shape(data_2d_array)
      endif  !! if sfcHeadRt

      call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
           start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
           'dart_vector_to_model_files', 'put_var '//trim(string2))
      deallocate(data_2d_array)

   elseif ( numdims == 3) then

      allocate(data_3d_array(progvar(ivar)%dimlens(1), &
                             progvar(ivar)%dimlens(2), &
                             progvar(ivar)%dimlens(3)) )

      ! This applies physical limits 
      call vector_to_prog_var(state_vector, ivar, data_3d_array,limit=.true.)
      ! In the case sh2ox, range restriction applies only to the high res field only. 
      ! However, the above range restriction is not applied to low resolution  
      ! "diagnostic" variables, such as smc in smc = sh2ox + sice. 
      ! How to handle that? Doing this currently in aggSh2ox. 
      ! Generally, I think range restriction needs to be a separate subroutine which 
      ! is not embedded in vector_to_prog_var. Vector_to_prog_var should do only that, 
      ! then (somthing named like) restrict_var_ranges should be called, then 
      ! the variable(s) should be written to file. 

      ! Ad hoc handling of "dual resolution" variable
      ! Aggregate then put weights and coarseRes liquid soil moisture in the restart. 
      ! Should this be moved to vector_to_3d_prog_var? also have to deal with change of 
      ! resolution above and somehow with writing the weights here. 
      if (trim(varname) .eq. 'sh2ox') then

         ! sh2ox is properly constrained in vector_to_prog_var (...,limit=.true.).
         ! However it lives as coarse res soil moisture and weights in the restart files. 
         ! Spatially aggregate and solve the weights and the total soil moisture. 
         allocate(smcAgg(coarse3dShape(1),coarse3dShape(2),coarse3dShape(3)), &
                  smcAggWeights(fine3dShape(1),fine3dShape(2),fine3dShape(3)) )
         sh2oDisag = data_3d_array
         deallocate(data_3d_array)
         call aggSh2ox(smcAgg, smcAggWeights)  !! also sets the global smc variable
         allocate(data_3d_array(coarse3dShape(1),coarse3dShape(2),coarse3dShape(3)))
         !! sh2ox onfile = smcAgg
         data_3d_array = smcAgg  !! this gets set just outside this if block
         deallocate(smcAgg)

         varidOrig = varid !! this will be written to file after this if statment
         ncNdimsOrig = ncNdims
         
         !! output sh2owgts to the ncid file
         !! i'm kinda cheating by not doing all the dimension checks...
         varname = 'sh2owgt'
         mycount(1:3)=fine3dShape(1:3)
         string2 = trim(filename)//' '//trim(varname)  ! for diagnostics
         call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
              'dart_vector_to_model_files', 'inq_varid '//trim(string2))
         call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
              'dart_vector_to_model_files', 'inquire '//trim(string2))
         call nc_check(nf90_put_var(ncFileID, VarID, smcAggWeights, &
              start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
              'dart_vector_to_model_files', 'put_var '//trim(string2))
         deallocate(smcAggWeights)
         ! Make note that the variable has been updated by DART
         call nc_check(nf90_Redef(ncFileID),'dart_vector_to_model_files', 'redef '//trim(filename))
         call nc_check(nf90_put_att(ncFileID, VarID,'DART','variable modified by DART'),&
              'dart_vector_to_model_files', 'modified '//trim(varname))
         call nc_check(nf90_enddef(ncfileID),'dart_vector_to_model_files','state enddef '//trim(filename))

         !! output smc to the ncid file. still not doing the dim checks. Note this is coarse res.
         varname = 'smc'
         mycount(1:3)=coarse3dShape(1:3)
         string2 = trim(filename)//' '//trim(varname)  ! for diagnostics
         call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
              'dart_vector_to_model_files', 'inq_varid '//trim(string2))
         call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
              'dart_vector_to_model_files', 'inquire '//trim(string2))
         call nc_check(nf90_put_var(ncFileID, VarID, smc, &
              start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
              'dart_vector_to_model_files', 'put_var '//trim(string2))
         ! Make note that the variable has been updated by DART
         call nc_check(nf90_Redef(ncFileID),'dart_vector_to_model_files', 'redef '//trim(filename))
         call nc_check(nf90_put_att(ncFileID, VarID,'DART','variable modified by DART'),&
              'dart_vector_to_model_files', 'modified '//trim(varname))
         call nc_check(nf90_enddef(ncfileID),'dart_vector_to_model_files','state enddef '//trim(filename))

         !! reset the current variable to sh2ox
         varname = 'sh2ox'
         string2 = trim(filename)//' '//trim(varname)  ! for diagnostics
         varid = varidOrig
         ncNdims = ncNdimsOrig
         mycount(1:3) = coarse3dShape(1:3)

      endif

      call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
           start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
           'dart_vector_to_model_files', 'put_var '//trim(string2))
      deallocate(data_3d_array)

   else

      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'dart_vector_to_model_files', string1, &
           source,revision,revdate)
   endif

   ! Make note that the variable has been updated by DART
   call nc_check(nf90_Redef(ncFileID),'dart_vector_to_model_files', 'redef '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VarID,'DART','variable modified by DART'),&
        'dart_vector_to_model_files', 'modified '//trim(varname))
   call nc_check(nf90_enddef(ncfileID),'dart_vector_to_model_files','state enddef '//trim(filename))

enddo UPDATE

call nc_check(nf90_close(ncFileID),'dart_vector_to_model_files','close '//trim(filename))
ncFileID = 0


end subroutine dart_vector_to_model_files

!===============================================================================
subroutine aggSh2ox(sh2ox, sh2owgt)
! Taken from HYDRO_drv/module_HYDRO_drv.F near line 536
! NOTE I renamed: smcrt->sh2oRt, smcWltRt->sh2oWltRt, and 
!                 smcMaxRt->sh2oMaxRt compared to that code

! Here we are updating (lsm grid:) sh2ox & smc and (rt grid:) sh2owgt
! requires (rt grid:) smcrt, sh2oRt, sh2oWltRt, sh2oMaxRt and (lsm grid:) sice
! Global/module variables: sh2oWltRt, sh2oMaxRt, sice

! Because this is independent of sfchead, I use separate subroutines for aggregation.

! Global
! real(R8), dimension(:,:,:)  :: sh2oDisag, smc

! liquid water content on the low res grid
real(R8), dimension(west_east,south_north,nsoil), intent(out) :: sh2ox
! weights to disaggregate 
real(R8), dimension(n_hlong,n_hlat,nsoil),        intent(out) :: sh2owgt

! local
real(r8), dimension(nsoil)                     :: sh2oaggrt
real(r8), dimension(n_hlong,n_hlat,nsoil)      :: sh2oRt
integer            :: i, j, ix, jx, krt, ixxrt, jyyrt, ss
integer            :: AGGFACXRT, AGGFACYRT, AGGFACTRT

jx = south_north
ix = west_east
aggfactrt = n_hlong/west_east

sh2oRt = sh2oDisag  ! global

!---------------------------------------------------------------------
! The following block could be done here, but it is done in 
! vector_to_3d_prog_var instead because that's the place constraints 
! are placed on assim state vectors prior to writing to file. 
!---------------------------------------------------------------------
! This happens because disaggregation happens in wrfHydro_to_dart
! but                     aggregation happens in dart_to_wrfHydro
! so these variables are lost from memory in the meanwhile b/c 
! these are different programs. Here, disaggregation is happening on 
! the prior (just as in wrfHydro_to_dart) because the fields are all read 
! from file. The soil moisture is different but not used, while all the other
! fields are the same (sice, sh2oWltRt, and sh2oMaxRt).
!if ((.not. allocated(sh2oWltRt)) .or. (.not. allocated(sh2oMaxRt))) then
!   call disagHydro(sh2oRt)  !! This allocates and calculates both of these
!   sh2oRt = sh2oDisag
!end if
! These constraints use the sice to constrain the liquid fraction.
!sh2oRt = max( min( sh2oRt,  sh2oMaxRt ), sh2oWltRt )
!---------------------------------------------------------------------

do J=1,JX
   do I=1,IX

      !! sh2oaggrt- layer totals per lsm grid cell
      do KRT=1,NSOIL
         SH2OAGGRT(KRT) = 0.
      end do !nsoil
      do AGGFACYRT=AGGFACTRT-1,0,-1
         do AGGFACXRT=AGGFACTRT-1,0,-1
            IXXRT=I*AGGFACTRT-AGGFACXRT
            JYYRT=J*AGGFACTRT-AGGFACYRT
            do KRT=1,NSOIL
               SH2OAGGRT(KRT) = SH2OAGGRT(KRT) &
                    + SH2ORT(IXXRT,JYYRT,KRT)
            end do !nsoil
         end do !aggfacxrt
      end do !aggfacyrt

      ! because sice is uniform & the liquid water (sh2oRt) was 
      ! conatrained above by a constant amt of ice
      ! aggregation will not violate smcWilt or smcMax at the coarse scale. 

      !! sh2ox - lsm grid layer is layer average of rt grid
      do KRT=1,NSOIL
         SH2OX(I,J,KRT) = SH2OAGGRT(KRT) &
              / (AGGFACTRT**2)
      end do !nsoil

      !! sh2owgt - rt grid. Weights are independent of total
      do AGGFACYRT=AGGFACTRT-1,0,-1
         do AGGFACXRT=AGGFACTRT-1,0,-1
            IXXRT=I*AGGFACTRT-AGGFACXRT
            JYYRT=J*AGGFACTRT-AGGFACYRT
            do KRT=1,NSOIL
               SH2OWGT(IXXRT,JYYRT,KRT) &
                    = SH2ORT(IXXRT,JYYRT,KRT) &
                    / SH2OX(I,J,KRT)
            end do !nsoil
         end do !aggfacxrt
      end do !aggfacyrt

   end do
end do

!! update smc = sice + sh2ox at the lsm/coarse grid resolution.
!! calculate sice
!! smc is global
smc=sh2ox+sice
do ss=1,nsoil
   smc(:,:,ss) = max( min(smc(:,:,ss), smcMax1), smcWlt1)
end do

end subroutine aggSh2ox


!===============================================================================
subroutine aggSfcHead(sfcHeadRt, infxsWgt)

! (The aggregation: sfcheadsubrt -> sfcheadrt, infxwgt)
! sfcHeadDisag == global sfcHeadDisag - fine resolution surface head
! sfcHeadRt    - coarse resolution surface head
! infxswgts    - fine resolution weights to go between the above surface head. 

! Taken from HYDRO_drv/module_HYDRO_drv.F near line 536

! Because this is independent of soil moisture, I use separate subroutines for aggregation.

! surface head (aka infiltration excess after infiltration) on the rt grid
! real(R8), dimension(:,:)  :: sfcHeadDisag

! sfc head on the low res grid
real(R8), dimension(west_east,south_north), intent(out) :: sfcHeadRt
! weights to disaggregate 
real(R8), dimension(n_hlong,n_hlat),        intent(out) :: infxsWgt

! local
real(r8)                     :: sfcHeadAggRt, lsmVol
integer                      :: i, j, ix, jx, krt, ixxrt, jyyrt, ss
integer                      :: aggFacXRt, aggFacYRt, aggFactRt

jx = south_north
ix = west_east
aggfactrt = n_hlong/west_east

do J = 1,JX
   do I = 1,IX

      sfcHeadAggRt = 0.
      lsmVol=0.

      do aggFacYRt = aggFactRt-1,0,-1
         do aggFacXRt = aggFactRt-1,0,-1
            IXXRT=I*aggFactRt-aggFacXRt
            JYYRT=J*aggFactRt-aggFacYRt

            sfcHeadAggRt = sfcHeadAggRt + sfcHeadDisag(IXXRT,JYYRT)
            lsmVol = lsmVol + sfcHeadDisag(IXXRT,JYYRT)*fineGridArea

         end do  !! aggFacYRt
      end do  !! aggFacXRt

      sfcHeadRt(I,J) = sfcHeadAggRt / (aggFactRt**2)

      do aggFacYRt = aggFactRt-1,0,-1
         do aggFacXRt = aggFactRt-1,0,-1
            IXXRT=I*aggFactRt-aggFacXRt
            JYYRT=J*aggFactRt-aggFacYRt

            if (lsmVol .gt. 0.) then
               infxsWgt(IXXRT,JYYRT) = sfcHeadDisag(IXXRT,JYYRT) * fineGridArea/lsmVol
            else
               infxsWgt(IXXRT,JYYRT) = 1./FLOAT(aggFactRt**2)
            end if
      
         end do  !! aggFacYRt
      end do  !! aggFacXRt
      
   end do
end do


end subroutine aggSfcHead


!===============================================================================
function get_state_time(ncid, filename, timeindex)
! The LSM restart files have "time".
! We are always using the 'most recent' which is, by defn, the last one.
! The time in the restart file is NOT the time at which the state is valid.
! It is one noah_timestep AHEAD of the valid time.
! Example:
! the noah_timestep = 3600 seconds &
! restart_frequency_hours = 1
! LSM restart filename is RESTART.2004010102_DOMAIN1 has
!     Times = '2004-01-01_02:00:00' ;
! The data is valid @ 2004-01-01_01:00:00, the previous noah_timestep
! (The time in the restart file name is the time of the restart, not the valid state time.
!  The valid state time was one noah_timestep previously.)

type(time_type) :: get_state_time
integer,           intent(in)  :: ncid
character(len=*),  intent(in)  :: filename
integer, optional, intent(out) :: timeindex

character(len=19), allocatable, dimension(:) :: datestring
integer               :: year, month, day, hour, minute, second
integer               :: DimID, VarID, strlen, ntimes
type(time_type)       :: filetime, timestep

! Get the dimensions for the strings of times

call nc_check(nf90_inq_dimid(ncid, 'Time', DimID), &
     'get_state_time','inq_dimid Time '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=ntimes), &
     'get_state_time','inquire_dimension Time '//trim(filename))

call nc_check(nf90_inq_dimid(ncid, 'DateStrLen', DimID), &
     'get_state_time','inq_dimid DateStrLen '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=strlen), &
     'get_state_time','inquire_dimension DateStrLen '//trim(filename))

if (strlen /= len(datestring)) then
   write(string1,*)'DatStrLen string length ',strlen,' /= ',len(datestring)
   call error_handler(E_ERR,'get_state_time', string1, source, revision, revdate)
endif

! Get all the Time strings, use the last one.
call nc_check(nf90_inq_varid(ncid, 'Times', VarID), &
     'get_state_time', 'inq_varid Times '//trim(filename))

allocate(datestring(ntimes))

call nc_check(nf90_get_var(ncid, VarID, datestring), &
     'get_state_time', 'get_var Times '//trim(filename))

read(datestring(ntimes),'(i4,5(1x,i2))') year, month, day, hour, minute, second

timestep       = set_time(noah_timestep, 0)
filetime       = set_date(year, month, day, hours=hour, minutes=minute, seconds=second)
get_state_time = filetime - timestep

if (present(timeindex)) timeindex = ntimes

if ( (do_output()) .and. debug > 99 ) write(*,*)'get_state_time: Last time string is '//trim(datestring(ntimes))
if ( (do_output()) .and. debug > 99 ) call print_date(get_state_time,' get_state_time: means valid time is ')

deallocate(datestring)

end function get_state_time



!===============================================================================
subroutine get_hrldas_constants(filename)
! Read the 'wrfinput' netCDF file for grid information, etc.
! This is all time-invariant, so we can mostly ignore the Time coordinate.
!
! MODULE variables set by this routine:
!    south_north
!    west_east
!    xlong
!    xlat

character(len=*), intent(in) :: filename

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, ncstart, nccount
character(len=NF90_MAX_NAME)          :: dimname

integer :: i, iunit, DimID, VarID, numdims, dimlen, xtype

call nc_check(nf90_open(adjustl(filename), NF90_NOWRITE, iunit), &
     'get_hrldas_constants', 'open '//trim(filename))

call nc_check(nf90_inq_dimid(iunit, 'south_north', DimID), &
     'get_hrldas_constants','inq_dimid south_north '//trim(filename))
call nc_check(nf90_inquire_dimension(iunit, DimID, len=south_north), &
     'get_hrldas_constants','inquire_dimension south_north '//trim(filename))

call nc_check(nf90_inq_dimid(iunit, 'west_east', DimID), &
     'get_hrldas_constants','inq_dimid west_east '//trim(filename))
call nc_check(nf90_inquire_dimension(iunit, DimID, len=west_east), &
     'get_hrldas_constants','inquire_dimension west_east '//trim(filename))

! Require that the xlong and xlat are the same shape.
allocate(xlong(west_east,south_north), xlat(west_east,south_north))

call nc_check(nf90_inq_varid(iunit, 'XLONG', VarID), &
     'get_hrldas_constants','inq_varid XLONG '//trim(filename))
call nc_check(nf90_inquire_variable(iunit, VarID, dimids=dimIDs, &
     ndims=numdims, xtype=xtype), &
     'get_hrldas_constants', 'inquire_variable XLONG '//trim(filename))

! Form the start/count such that we always get the 'latest' time.
ncstart(:) = 0
nccount(:) = 0
do i = 1,numdims
   write(string1,'(''inquire dimension'',i2,A)') i,trim(filename)
   call nc_check(nf90_inquire_dimension(iunit, dimIDs(i), name=dimname, len=dimlen), &
        'get_hrldas_constants', string1)
   ncstart(i) = 1
   nccount(i) = dimlen
   if ((trim(dimname) == 'Time') .or. (trim(dimname) == 'time')) then
      ncstart(i) = dimlen
      nccount(i) = 1
   endif
enddo

if ( (do_output()) .and. debug > 99 ) write(*,*)'DEBUG get_hrldas_constants ncstart is',ncstart(1:numdims)
if ( (do_output()) .and. debug > 99 ) write(*,*)'DEBUG get_hrldas_constants nccount is',nccount(1:numdims)

!get the longitudes
call nc_check(nf90_get_var(iunit, VarID, xlong, &
     start=ncstart(1:numdims), count=nccount(1:numdims)), &
     'get_hrldas_constants', 'get_var XLONG '//trim(filename))
where(xlong < 0.0_r8) xlong = xlong + 360.0_r8
where(xlong == 360.0_r8) xlong = 0.0_r8

! finally get the latitudes
call nc_check(nf90_inq_varid(iunit, 'XLAT', VarID), &
     'get_hrldas_constants','inq_varid XLAT '//trim(filename))
call nc_check(nf90_get_var(iunit, VarID, xlat, &
     start=ncstart(1:numdims), count=nccount(1:numdims)), &
     'get_hrldas_constants', 'get_var XLAT '//trim(filename))
where (xlat < -90.0_r8) xlat = -90.0_r8
where (xlat >  90.0_r8) xlat =  90.0_r8


end subroutine get_hrldas_constants

!===============================================================================
subroutine get_hydro_constants(filename)
! Read the 'GEO_FINEGRID_FLNM' netCDF file for grid information, etc.
! This is all time-invariant, so we can mostly ignore the Time coordinate.
!
! MODULE variables set by this routine:
!    n_hlat, n_hlong, n_link, n_basn, basnMask
!    hlong, hlat, linkLat, linkLong

! **** NOTE that all variables from this file (Fulldom) must  ****
! ****      be FLIPPED in y to match the noah/wrf model.      ****

character(len=*), intent(in) :: filename

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, ncstart, nccount
character(len=NF90_MAX_NAME)          :: dimname
integer, allocatable, dimension(:,:) ::  basnGrid
real(r8), allocatable, dimension(:,:) ::  hlongFlip, hlatFlip ! local dummies
integer, allocatable, dimension(:) ::  basnMaskTmp
real(r8), allocatable, dimension(:) ::  channelLon1D, channelLat1D
integer :: i, ii, jj, iunit, DimID, VarID, numdims, dimlen, xtype
integer :: indx, indx1, indx2, indx3, indx4, dumSum

call nc_check(nf90_open(adjustl(filename), NF90_NOWRITE, iunit), &
     'get_hydro_constants', 'open '//trim(filename))

call nc_check(nf90_inq_dimid(iunit, 'y', DimID), &
     'get_hydro_constants','inq_dimid y '//trim(filename))
call nc_check(nf90_inquire_dimension(iunit, DimID, len=n_hlat), &
     'get_hydro_constants','inquire_dimension y '//trim(filename))

call nc_check(nf90_inq_dimid(iunit, 'x', DimID), &
     'get_hydro_constants','inq_dimid x '//trim(filename))
call nc_check(nf90_inquire_dimension(iunit, DimID, len=n_hlong), &
     'get_hydro_constants','inquire_dimension x '//trim(filename))

!! module allocation
allocate(hlong(n_hlong, n_hlat), &
         hlat(n_hlong, n_hlat)) 
!! local allocation
allocate( basnGrid(n_hlong, n_hlat), &
         hlongFlip(n_hlong, n_hlat), &
          hlatFlip(n_hlong, n_hlat))

! Require that the xlong and xlat are the same shape.??
call nc_check(nf90_inq_varid(iunit, 'LONGITUDE', VarID), &
     'get_hydro_constants','inq_varid LONGITUDE '//trim(filename))
call nc_check(nf90_inquire_variable(iunit, VarID, dimids=dimIDs, &
     ndims=numdims, xtype=xtype), &
     'get_hydro_constants', 'inquire_variable LONGITUDE '//trim(filename))

! numdims = 2, these are all 2D fields
! Form the start/count such that we always get the 'latest' time.
ncstart(:) = 0
nccount(:) = 0
do i = 1,numdims
   write(string1,'(''inquire dimension'',i2,A)') i,trim(filename)
   call nc_check(nf90_inquire_dimension(iunit, dimIDs(i), name=dimname, len=dimlen), &
        'get_hydro_constants', string1)
   ncstart(i) = 1
   nccount(i) = dimlen
   if ((trim(dimname) == 'Time') .or. (trim(dimname) == 'time')) then
      ncstart(i) = dimlen
      nccount(i) = 1
   endif
enddo !i

if ( (do_output()) .and. debug > 99 ) write(*,*)'DEBUG get_hydro_constants ncstart is',ncstart(1:numdims)
if ( (do_output()) .and. debug > 99 ) write(*,*)'DEBUG get_hydro_constants nccount is',nccount(1:numdims)

!get the longitudes
call nc_check(nf90_get_var(iunit, VarID, hlong, &
     start=ncstart(1:numdims), count=nccount(1:numdims)), &
     'get_hydro_constants', 'get_var LONGITUDE '//trim(filename))
where(hlong < 0.0_r8) hlong = hlong + 360.0_r8
where(hlong == 360.0_r8) hlong = 0.0_r8

!get the latitudes
call nc_check(nf90_inq_varid(iunit, 'LATITUDE', VarID), &
     'get_hydro_constants','inq_varid LATITUDE '//trim(filename))
call nc_check(nf90_get_var(iunit, VarID, hlat, &
     start=ncstart(1:numdims), count=nccount(1:numdims)), &
     'get_hydro_constants', 'get_var LATITUDE '//trim(filename))
where (hlat < -90.0_r8) hlat = -90.0_r8
where (hlat >  90.0_r8) hlat =  90.0_r8

! Flip the longitues and latitudes
do ii=1,n_hlong
   do jj=1,n_hlat
      hlongFlip(ii,jj) = hlong(ii,n_hlat-jj+1)
       hlatFlip(ii,jj) =  hlat(ii,n_hlat-jj+1)
   end do
end do
hlong = hlongFlip
hlat  = hlatFlip
deallocate(hlongFlip, hlatFlip)

! get the basin grid - this wont need flipped as unique values are  
! packed in to a 1D array without geolocaiton information.
call nc_check(nf90_inq_varid(iunit, 'basn_msk', VarID), &
     'get_hydro_constants','inq_varid basn_msk '//trim(filename))
call nc_check(nf90_get_var(iunit, VarID, basnGrid, &
     start=ncstart(1:numdims), count=nccount(1:numdims)), &
     'get_hydro_constants', 'get_var basn_msk '//trim(filename))

! make it a 1D array of single values
! qeustion is how to localize this, since it has no coordinates.
! for now going to use the average lat and lon of each basin
allocate(basnMaskTmp(maxval(basnGrid)))  !local

basnMaskTmp(:) = -9999
indx2=0
do indx = 1,maxval(basnGrid)
   if(any(basnGrid == indx)) then
      indx2=indx2+1
      basnMaskTmp(indx2) = indx
   end if
enddo
n_basn = indx2
allocate(basnMask(n_basn),basnLon(n_basn),basnLat(n_basn))
where(basnMaskTmp .ne. -9999)
   basnMask(:) = basnMaskTmp(1:n_basn)
end where

! geolocate the basins
do indx = 1, n_basn
      basnLon = sum(hlong, basnGrid .eq. indx) / count(basnGrid .eq. indx)
      basnLat = sum( hlat, basnGrid .eq. indx) / count(basnGrid .eq. indx)
enddo 



!get the channelgrid
! i'm doing this exactly to match how it's done in the wrfHydro code 
! (line 1044 of /home/jamesmcc/WRF_Hydro/ndhms_wrf_hydro/trunk/NDHMS/Routing/module_HYDRO_io.F)
! so that the output set of inidces correspond to the grid in the Fulldom file 
! and therefore these can be used to grab other channel attributes in that file. 
! but the code is really long so I've put it in a module subroutine. 
! Dont need to flip lat and lon in this (already done) but will flip other vars from Fulldom file.
call getChannelCoords(filename, iunit, numdims, ncstart, nccount)


call nc_check(nf90_close(iunit), 'get_hydro_constants '//trim(filename))

deallocate(basnMaskTmp, basnGrid)

end subroutine get_hydro_constants


!===============================================================================
subroutine define_var_dims(ivar, ncid, memberdimid, unlimiteddimid, ndims, dimids)
! I am trying to preserve the shape of the NOAH variable as much as possible.
!
! noah_variable(Time,       soil_layers_stag, south_north, west_east) becomes
! noah_variable(time, copy, soil_layers_stag, south_north, west_east)
!
! Since 'time' or 'Time' is the unlimited dimension in both ... I can skip it
! in the DEFDIM loop.

integer,               intent(in)  :: ivar, ncid, memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

integer :: i, mydimid
character(len=obstypelength) :: theDimname

ndims = 0

DEFDIM : do i = 1,progvar(ivar)%numdims

   if ((trim(progvar(ivar)%dimnames(i)) == 'Time') .or. &
        (trim(progvar(ivar)%dimnames(i)) == 'time')) cycle DEFDIM


   theDimname = progvar(ivar)%dimnames(i)
   !! upshot of having two restart files with different dimnames. 
   !! I adopt the wrfHydro restart dimnames with goal of preserving the 
   !! the Noah vs NoahMP dimension order.
   !! Only have to chance the dimnames of LSM variable components. 
   if (progvar(ivar)%component == 'LSM' ) then 
      if (trim(theDimname) == 'south_north')      theDimname = 'iy'
      if (trim(theDimname) == 'west_east')        theDimname = 'ix'
      if (trim(theDimname) == 'soil_layers_stag') theDimname = 'depth'
   endif


   !! other exceptions....
   !! though hydro sh2ox is output on the coarse grid with the weights sh2owgts
   !! i'd rather have the fine grid sh2ox in the diagnostic files. 
   if (progvar(ivar)%component == 'HYDRO' ) then 
      if (progvar(ivar)%varname == 'sh2ox' .OR. &
          progvar(ivar)%varname == 'sfcheadrt' ) then 
         if (trim(theDimname) == 'ix')      theDimname = 'ixrt'
         if (trim(theDimname) == 'iy')      theDimname = 'iyrt'
      endif
   endif

   if ( (do_output()) .and. debug > 0) then
      write(*,*)trim(progvar(ivar)%component), &
           trim(progvar(ivar)%varname), &
           progvar(ivar)%numdims, &
           trim(theDimname)
   endif

   call nc_check(nf90_inq_dimid(ncid=ncid, name=theDimname, dimid=mydimid), &
        'define_var_dims','inq_dimid '//trim(progvar(ivar)%dimnames(i)))

   ndims         = ndims + 1
   dimids(ndims) = mydimid

enddo DEFDIM

ndims = ndims + 1
dimids(ndims) = memberdimid
ndims = ndims + 1
dimids(ndims) = unlimitedDimid

if ( (do_output()) .and. debug > 99 ) then
   write(logfileunit,*)
   write(logfileunit,*)'define_var_dims knowledge'
   write(logfileunit,*)trim(progvar(ivar)%varname),' has dimnames ', &
        progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
   write(logfileunit,*)' thus dimids ',dimids(1:ndims)
   write(     *     ,*)
   write(     *     ,*)'define_var_dims knowledge'
   write(     *     ,*)trim(progvar(ivar)%varname),' has dimnames ', &
        progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
   write(     *     ,*)' thus dimids ',dimids(1:ndims)

endif

return
end subroutine define_var_dims


!===============================================================================
subroutine vector_to_1d_prog_var(x, ivar, data_1d_array, limit)
! convert the values from a 1d array, starting at an offset, into a 1d array.
!
! If the optional argument (ncid) is specified, some additional
! processing takes place. The variable in the netcdf is read.
! This must be the same shape as the intended output array.
! Anywhere the DART MISSING code is encountered in the input array,
! the corresponding (i.e. original) value from the netCDF file is
! used.

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:),   intent(out) :: data_1d_array
logical,  optional,       intent(in)  :: limit

integer :: i,ii

! unpack the right part of the DART state vector into a 1D array.

ii = progvar(ivar)%index1

do i = 1, progvar(ivar)%dimlens(1)
   data_1d_array(i) = x(ii)
   ii = ii + 1
enddo

if (present(limit)) then
   if ( limit ) then
      if (    progvar(ivar)%rangeRestricted == 1) then
         where(data_1d_array < progvar(ivar)%minvalue) data_1d_array = progvar(ivar)%minvalue
      elseif (progvar(ivar)%rangeRestricted == 2) then
         where(data_1d_array > progvar(ivar)%maxvalue) data_1d_array = progvar(ivar)%maxvalue
      elseif (progvar(ivar)%rangeRestricted == 3) then
         where(data_1d_array < progvar(ivar)%minvalue) data_1d_array = progvar(ivar)%minvalue
         where(data_1d_array > progvar(ivar)%maxvalue) data_1d_array = progvar(ivar)%maxvalue
      elseif (progvar(ivar)%rangeRestricted == 4) then
         write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' has rangeRestricted==4.'
         write(string2, *)'No code written to restrict its range, however.'
         call error_handler(E_ERR,'vector_to_1d_prog_var', string1, &
              source, revision, revdate, text2=string2)
      endif
   endif
endif

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' packed wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_1d_prog_var', string1, &
        source, revision, revdate, text2=string2)
endif

end subroutine vector_to_1d_prog_var


!===============================================================================
subroutine vector_to_2d_prog_var(x, ivar, data_2d_array, limit)
! convert the values from a 1d array, starting at an offset,
! into a 2d array.
!
real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:,:), intent(out) :: data_2d_array
logical,  optional,       intent(in)  :: limit

integer :: i,j,ii

! unpack the right part of the DART state vector into a 2D array.

ii = progvar(ivar)%index1

do j = 1,progvar(ivar)%dimlens(2)
do i = 1,progvar(ivar)%dimlens(1)
   data_2d_array(i,j) = x(ii)
   ii = ii + 1
enddo
enddo

if (present(limit)) then
   if ( limit ) then
      if (    progvar(ivar)%rangeRestricted == 1) then
         where(data_2d_array < progvar(ivar)%minvalue) data_2d_array = progvar(ivar)%minvalue
      elseif (progvar(ivar)%rangeRestricted == 2) then
         where(data_2d_array > progvar(ivar)%maxvalue) data_2d_array = progvar(ivar)%maxvalue
      elseif (progvar(ivar)%rangeRestricted == 3) then
         where(data_2d_array < progvar(ivar)%minvalue) data_2d_array = progvar(ivar)%minvalue
         where(data_2d_array > progvar(ivar)%maxvalue) data_2d_array = progvar(ivar)%maxvalue
      elseif (progvar(ivar)%rangeRestricted == 4) then
         write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' has rangeRestricted==4.'
         write(string2, *)'No code written to restrict its range, however.'
         call error_handler(E_ERR,'vector_to_2d_prog_var', string1, &
              source, revision, revdate, text2=string2)
      endif
   endif
endif

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_2d_prog_var', string1, &
        source, revision, revdate, text2=string2)
endif

end subroutine vector_to_2d_prog_var


!===============================================================================
subroutine vector_to_3d_prog_var(x, ivar, data_3d_array, limit)
! convert the values from a 1d array, starting at an offset,
! into a 3d array.
!
real(r8), dimension(:),     intent(in)  :: x
integer,                    intent(in)  :: ivar
real(r8), dimension(:,:,:), intent(out) :: data_3d_array
logical,  optional,         intent(in)  :: limit

real(R8),dimension(:,:,:),allocatable :: dum !! dum to disagHydro to init global sh2oWltRt and sh2oMaxRt
integer :: i,j,k,ii

! unpack the right part of the DART state vector into a 3D array.

ii = progvar(ivar)%index1

do k = 1,progvar(ivar)%dimlens(3)
do j = 1,progvar(ivar)%dimlens(2)
do i = 1,progvar(ivar)%dimlens(1)
   data_3d_array(i,j,k) = x(ii)
   ii = ii + 1
enddo
enddo
enddo

if (present(limit)) then
   if ( limit ) then
      if (    progvar(ivar)%rangeRestricted == 1) then
         where(data_3d_array < progvar(ivar)%minvalue) data_3d_array = progvar(ivar)%minvalue
      elseif (progvar(ivar)%rangeRestricted == 2) then
         where(data_3d_array > progvar(ivar)%maxvalue) data_3d_array = progvar(ivar)%maxvalue
      elseif (progvar(ivar)%rangeRestricted == 3) then
         where(data_3d_array < progvar(ivar)%minvalue) data_3d_array = progvar(ivar)%minvalue
         where(data_3d_array > progvar(ivar)%maxvalue) data_3d_array = progvar(ivar)%maxvalue
      elseif (progvar(ivar)%rangeRestricted == 4) then
         if (trim(progvar(ivar)%varname) == 'sh2ox') then 
            ! Constrain the liquid soil moisture on high-res grid using
            ! the wilt and  max which assume ice is NOT adjusted.       
            ! This "if block" sets sh2oWltRt and sh2oMaxRt if they happen to be unallocated
            ! which happens every time because disaggregation happens in wrfHydro_to_dart
            ! but                                 aggregation happens in dart_to_wrfHydro
            ! so these variables are lost from memory in the meanwhile b/c 
            ! these are different programs. Here, disaggregation is happening on 
            ! the prior (just as in wrfHydro_to_dart) because the fields are all read 
            ! from file. The soil moisture is different but not used, while all the other
            ! fields are the same (sice, sh2oWltRt, and sh2oMaxRt).
            if ( (.not. allocated(sh2oWltRt)) .or. (.not. allocated(sh2oMaxRt)) ) then
               call disagHydro()  
            endif
            data_3d_array = max( min( data_3d_array, sh2oMaxRt ), sh2oWltRt )
         else 
            write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' has rangeRestricted==4.'
            write(string2, *)'No code written to restrict its range, however.'
            call error_handler(E_ERR,'vector_to_2d_prog_var', string1, &
                 source, revision, revdate, text2=string2)
         endif
      endif
   endif
endif

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_3d_prog_var', string1, &
        source, revision, revdate, text2=string2)
endif

end subroutine vector_to_3d_prog_var


!===============================================================================
subroutine get_model_timestepping(day,hour,dynamical,output,forcing,restart)
integer,          intent(out) :: day,hour,dynamical,output,forcing,restart
day       = kday
hour      = khour
dynamical = noah_timestep
output    = output_timestep
forcing   = forcing_timestep
restart   = restart_frequency_hours*3600
end subroutine get_model_timestepping

!===============================================================================
! Painful amount of code for getting the channell lat/lon/ele which matches
! the wrfHydro state variable
subroutine getChannelCoords(filename, iunit, numdims, ncstart, nccount)

integer,               intent(in) :: iunit
integer,               intent(in) :: numdims
character(len=*),      intent(in) :: filename
integer, dimension(:), intent(in) :: ncstart
integer, dimension(:), intent(in) :: nccount

integer                               :: IXRT,JXRT
real(r8), allocatable, dimension(:,:) :: ELRT, ELRT_in
integer,  allocatable, dimension(:,:) :: DIRECTION, LAKE_MSKRT, channelGrid
integer,  allocatable, dimension(:,:) :: DIRECTION_in, LAKE_MSKRT_in, channelGrid_in

integer :: VarID, tmp, cnt, i, j, jj
character(len=155) :: header

!integer                                  :: NLAKES
!real(r4), allocatable, dimension(NLAKES) :: HRZAREA, LAKEMAXH, WEIRC, WEIRL
!real(r4), allocatable, dimension(NLAKES) :: ORIFICEC, ORIFICEA, ORIFICEE
!real(r4), allocatable, dimension(NLAKES) :: LATLAKE,LONLAKE,ELEVLAKE

! allocate the local variables
! these grid ones have to be flipped on y.
allocate(channelGrid_in(n_hlong,n_hlat), LAKE_MSKRT_in(n_hlong,n_hlat), &
           DIRECTION_in(n_hlong,n_hlat),       ELRT_in(n_hlong,n_hlat) )
allocate(   channelGrid(n_hlong,n_hlat),    LAKE_MSKRT(n_hlong,n_hlat), &
              DIRECTION(n_hlong,n_hlat),          ELRT(n_hlong,n_hlat) )

call nc_check(nf90_inq_varid(iunit, 'CHANNELGRID', VarID), &
     'getChannelCoords','inq_varid CHANNELGRID '//trim(filename))
call nc_check(nf90_get_var(iunit, VarID, channelGrid_in, &
     start=ncstart(1:numdims), count=nccount(1:numdims)), &
     'getChannelCoords', 'get_var CHANNELGRID '//trim(filename))

call nc_check(nf90_inq_varid(iunit, 'LAKEGRID', VarID), &
     'getChannelCoords','inq_varid LAKEGRID '//trim(filename))
call nc_check(nf90_get_var(iunit, VarID, LAKE_MSKRT_in, &
     start=ncstart(1:numdims), count=nccount(1:numdims)), &
     'getChannelCoords', 'get_var LAKEGRID '//trim(filename))

call nc_check(nf90_inq_varid(iunit, 'FLOWDIRECTION', VarID), &
     'getChannelCoords','inq_varid FLOWDIRECTION '//trim(filename))
call nc_check(nf90_get_var(iunit, VarID, DIRECTION_in, &
     start=ncstart(1:numdims), count=nccount(1:numdims)), &
     'getChannelCoords', 'get_var FLOWDIRECTION '//trim(filename))

call nc_check(nf90_inq_varid(iunit, 'TOPOGRAPHY', VarID), &
     'getChannelCoords','inq_varid TOPOGRAPHY '//trim(filename))
call nc_check(nf90_get_var(iunit, VarID, ELRT_in, &
     start=ncstart(1:numdims), count=nccount(1:numdims)), &
     'getChannelCoords', 'get_var TOPOGRAPHY '//trim(filename))

ixrt = n_hlong
jxrt = n_hlat

! wrfHydro flips the y dimension of the variables from the Fulldom file
do i=1,ixrt
   do j=1,jxrt
      channelGrid(i,j) = channelGrid_in(i,jxrt-j+1)
      LAKE_MSKRT(i,j)  =  LAKE_MSKRT_in(i,jxrt-j+1)
      DIRECTION(i,j)   =   DIRECTION_in(i,jxrt-j+1)
      ELRT(i,j)        =        ELRT_in(i,jxrt-j+1)
   end do
end do
deallocate(channelGrid_in, LAKE_MSKRT_in, DIRECTION_in, ELRT_in)

! subset to the 1D channel network as presented in the hydro restart file.
n_link = sum(channelGrid_in*0+1, mask = channelGrid .ge. 0)

! allocate the necessary wrfHydro variables with module scope
allocate(channelIndsX(n_link), channelIndsY(n_link), linkLong(n_link), linkLat(n_link))

! temp fix for buggy Arc export...
do j=1,jxrt
   do i=1,ixrt
      if(DIRECTION(i,j).eq.-128) DIRECTION(i,j)=128
   end do
end do

if ((CHANRTSWCRT.eq.1.or.CHANRTSWCRT.eq.2).and.channel_option .ne. 3) then
   ! not routing on grid, read from file
   write(string1, *)'Reach based routing not currently enabled in DART.'
   call error_handler(E_ERR,'getChannelCoords', string1, &
        source, revision, revdate, text2=string1)

elseif ((CHANRTSWCRT.eq.1.or.CHANRTSWCRT.eq.2).and.channel_option.eq.3) then

   cnt = 0

   do j = 1, JXRT  !rows
      do i = 1 ,IXRT   !colsumns
         if (CHANNELGRID(i, j) .ge. 0) then !get its direction and assign its elevation and order
            if ((DIRECTION(i, j) .eq. 64) .and. (j + 1 .le. JXRT) .and. &
                 (CHANNELGRID(i,j+1).ge.0) ) then !North
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j  !! again have to match the flip
            else if ((DIRECTION(i, j) .eq. 128) .and. (i + 1 .le. IXRT) &
                 .and. (j + 1 .le. JXRT) .and. (CHANNELGRID(i+1,j+1).ge.0) ) then !North East
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
            else if ((DIRECTION(i, j) .eq. 1) .and. (i + 1 .le. IXRT) &
                 .and. (CHANNELGRID(i+1,j).ge.0) ) then !East
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
            else if ((DIRECTION(i, j) .eq. 2) .and. (i + 1 .le. IXRT) &
                 .and. (j - 1 .ne. 0) .and. (CHANNELGRID(i+1,j-1).ge.0) ) then !south east
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
            else if ((DIRECTION(i, j) .eq. 4) .and. (j - 1 .ne. 0).and.(CHANNELGRID(i,j-1).ge.0) ) then !due south
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
            else if ((DIRECTION(i, j) .eq. 8) .and. (i - 1 .gt. 0) &
                 .and. (j - 1 .ne. 0) .and. (CHANNELGRID(i-1,j-1).ge.0)) then !south west
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
            else if ((DIRECTION(i, j) .eq. 16) .and. (i - 1 .gt. 0).and.(CHANNELGRID(i-1,j).ge.0) ) then !West
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
            else if ((DIRECTION(i, j) .eq. 32) .and. (i - 1 .gt. 0) &
                 .and. (j + 1 .le. JXRT) .and. (CHANNELGRID(i-1,j+1).ge.0) ) then !North West
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
            else
               write(string1,*)"NO MATCH", i,j
               call error_handler(E_MSG,'getChannelCoords',string1)
            end if

         end if !CHANNELGRID check for this node

      end do
   end do

!   print *, "after exiting the channel, this many nodes", cnt
!   write(*,*) " "

   !Find out if the boundaries are on an edge
   do j = 1,JXRT
      do i = 1 ,IXRT
         if (CHANNELGRID(i, j) .ge. 0) then !get its direction
            !-- 64's can only flow north
            if (((DIRECTION(i, j).eq. 64) .and. (j + 1 .gt. JXRT)) .or. &
                 ((DIRECTION(i, j).eq. 64) .and. (j < jxrt) .and.  &
                 (CHANNELGRID(i,j+1) .lt. 0))) then !North
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
!               print *, "Pour Point N"
            else if ( ((DIRECTION(i, j) .eq. 128) .and. (i + 1 .gt. IXRT))  &
                                !-- 128's can flow out of the North or East edge
                 .or.  ((DIRECTION(i, j) .eq. 128) .and. (j + 1 .gt. JXRT))  &
                                !   this is due north edge
                 .or.  ((DIRECTION(i, j) .eq. 128) .and. (i<ixrt .and. j<jxrt) .and. &
                 (CHANNELGRID(i + 1, j + 1).lt.0))) then !North East
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
!               print *, "Pour Point NE"
            else if (((DIRECTION(i, j) .eq. 1) .and. (i + 1 .gt. IXRT)) .or. &    !-- 1's can only flow due east
                 ((DIRECTION(i, j) .eq. 1) .and. (i<ixrt) .and. (CHANNELGRID(i + 1, j) .lt. 0))) then !East
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
!               print *, "Pour Point E"
            else if ( ((DIRECTION(i, j) .eq. 2) .and. (i + 1 .gt. IXRT))    &      !-- 2's can flow out of east or south edge
                 .or. ((DIRECTION(i, j) .eq. 2) .and. (j - 1 .eq. 0))       &      !-- this is the south edge
                 .or. ((DIRECTION(i, j) .eq. 2) .and. (i<ixrt .and. j>1) .and.(CHANNELGRID(i + 1, j - 1) .lt.0))) then !south east
               cnt = cnt + 1
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
!               print *, "Pour Point SE"
            else if (((DIRECTION(i, j) .eq. 4) .and. (j - 1 .eq. 0)) .or. &       !-- 4's can only flow due south
                 ((DIRECTION(i, j) .eq. 4) .and. (j>1) .and.(CHANNELGRID(i, j - 1) .lt. 0))) then !due south
               cnt = cnt + 1
               !ZELEV(cnt) = ELRT(i,j)
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
!               print *, "Pour Point S"
            else if ( ((DIRECTION(i, j) .eq. 8) .and. (i - 1 .le. 0))      &      !-- 8's can flow south or west
                 .or.  ((DIRECTION(i, j) .eq. 8) .and. (j - 1 .eq. 0))      &      !-- this is the south edge
                 .or.  ((DIRECTION(i, j).eq.8).and. (i>1 .and. j>1) .and.(CHANNELGRID(i - 1, j - 1).lt.0))) then !south west
               cnt = cnt + 1
               !ZELEV(cnt) = ELRT(i,j)
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
!               print *, "Pour Point SW"
            else if (((DIRECTION(i, j) .eq. 16) .and. (i - 1 .le.0) ) &                  !16's can only flow due west
                 .or.((DIRECTION(i, j).eq.16) .and. (i>1) .and.(CHANNELGRID(i - 1, j).lt.0))) then !West
               cnt = cnt + 1
               !ZELEV(cnt) = ELRT(i,j)
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
!               print *, "Pour Point W"
            else if ( ((DIRECTION(i, j) .eq. 32) .and. (i - 1 .le. 0))      &      !-- 32's can flow either west or north
                 .or.  ((DIRECTION(i, j) .eq. 32) .and. (j + 1 .gt. JXRT))   &      !-- this is the north edge
                 .or.  ((DIRECTION(i, j).eq.32) .and. (i>1 .and. j<jxrt) .and.(CHANNELGRID(i - 1, j + 1).lt.0))) then !North West
               cnt = cnt + 1
               !ZELEV(cnt) = ELRT(i,j)
               linkLat(cnt) = hlat(i,j)
               linkLong(cnt) = hlong(i,j)
               channelIndsX(cnt) = i
               channelIndsY(cnt) = j
!               print *, "Pour Point NW"
            endif
         endif !CHANNELGRID check for this node
      end do
   end do
endif

!close(79)

deallocate(channelGrid, LAKE_MSKRT, DIRECTION, ELRT)

if (cnt .ne. n_link) then
   write(string1,*) 'Error with number of links in the channel grid.'
   call error_handler(E_ERR,'getChannelCoords',string1,source,revision,revdate)
endif

end subroutine getChannelCoords

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

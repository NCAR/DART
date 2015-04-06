! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! $Id$


!>  This is the interface module between DART and the atmospheric components of CESM; 
!>  CAM, WACCM, CAM-Chem (, ...?).  
!>  It contains the required 16 interface procedures, as specified by DART.  
!>  It also contains several utility routines which help translate between CAM and DART 
!>  formats, and deal with time.
!>  It is used by filter, perfect_model_obs, $dart/models/cam/{cam_to_dart,dart_to_cam}.
!>
!>  This module handles all of the eulerian, finite volume, and spectral element (HOMME)
!>  dynamical core versions of CAM.  The first two have logically rectangular grids,
!>  while CAM-SE uses the cubed sphere (non-rectangular) horizontal grid.
!>
!>  It contains a perturburbation routine for generating initial ensembles,
!>  but does not provide adv_1step or init_conditions because CAM is a separate executable
!>  and cannot be called as a subroutine.
!>
!>  This module intercepts the get_close_obs() calls and can alter the distances
!>  for obs near the top of the model to reduce the impact on the state near the
!>  top.
!>
!>  The coordinate orders of fields stored in various forms have also been simplified.
!>  For example; various vintages of CAM 3D fields may be read in with (lon, lat, lev) or
!>  (lon, lev, lat).  These are uniformly converted to (lev, lon, lat) for use in model_mod.
!>  This latter form is different than pre MPI model_mods.  Then such fields are stored in
!>  the state vector with the same coordinate order.  They are converted to the *modern*
!>  CAM coordinate order when written to caminput.nc files.
!>
!>  If a user wants to add new TYPE_s to the state vector,
!>  then more KIND_s may be needed in the obs_kind_mod and the 'use obs_kind_mod' statement.
!> 
!>  Observations below the lowest model level (including surface observations) and above
!>  the highest model level cannot be assimilated (yet).  The spatial extent of observations
!>  can be further restricted using model_nml namelist variables.
!> 
!>  MODULE ORGANIZATION (search for the following strings to find the corresponding section)
!>
!>      'use' statements
!>       Global storage for describing cam model class
!>       Namelist variables with default values
!>       Derived parameters
!>       static_init_model section
!>       Module I/O to/from DART and files
!>       model_interpolate section
!>       Vector-field translations
!>       get_close_obs section
!>       Utility routines; called by several main subroutines
!>       Stubs not used by cam/model_mod (this is not all of them)
!>
!>  See the subversion code logs for history of this module.
!>  There is an html description of this module in ./model_mod.html.
!>

module model_mod

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  This module keeps a copy of the ensemble mean in module global storage and
!  uses it for computing the pressure-to-height conversions.
!
!  During the assimilation stage, only a piece of the state vector is available to each
!  process, and each process calls parts of model_mod.  In order to handle the conversion
!  of vertical coordinates of obs and/or state variables into a consistent coordinate,
!  an entire state vector is needed, so the ensemble mean is passed to model_mod before
!  the assimilation starts.  This is NOT done for model_interpolate; the whole vector is
!  available, and should be used.  All locations are now converted to a standard coordinate
!  (pressure or log(P0/pressure), aka scale height), instead of always converting the state 
!  vertical location to that of the ob.  The highest_obs_level and ..._height_m variables 
!  are derived from highest_obs_pressure_Pa namelist variable.
!
!  The coordinate orders of fields stored in various forms have also been simplified.
!  For example; various vintages of CAM 3D fields may be read in with (lon, lat, lev) or
!  (lon, lev, lat).  These are uniformly converted to (lev, lon, lat) for use in model_mod.
!  This latter form is different than pre MPI model_mods.  Then such fields are stored in
!  the state vector with the same coordinate order.  They are converted back to the modern
!  CAM coordinate order when written to caminput.nc files.
!  Surface pressure may be needed on the A-grid (thermodynamic variables) and grids staggered
!  relative to the A-grid.   Currently, PS for the A-grid and for the 2 staggered grids is
!  stored for global access for the (re)calculation of pressures and heights on model levels
!  as needed.  The 3d pressure on the A-grid or cubed sphere grid is calculated and stored, 
!  but 3d pressure on staggered grids is calculated as needed. 

!  If a user wants to add new TYPE_s to the state vector,
!  then more KIND_s may be needed in the obs_kind_mod and the 'use obs_kind_mod' statement.

!  The coordinates of CAM (lats, lons, ncol, etc.) and their dimensions  and attributes are
!  read into globally accessible data structures (see grid_1d_type).
!
!     MODULE ORGANIZATION (search for the following strings to find the corresponding section)
!
!         'use' statements
!          Global storage for describing cam model class
!          Namelist variables with default values
!          Derived parameters
!          static_init_model section
!          Module I/O to/from DART and files
!          model_interpolate section
!          Vector-field translations
!          get_close_obs section
!          Utility routines; called by several main subroutines
!          Stubs not used by cam/model_mod (this is not all of them)

!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

!  USE statements

use netcdf
use typeSizes

use types_mod,         only : r8, MISSING_I, MISSING_R8, gravity_const => gravity, &
                              PI, DEG2RAD, RAD2DEG, obstypelength, earth_radius
! FIXME; these constants should be consistent with CESM, not necessarily with DART.
!          add after verification against Hui's tests;  gas_constant_v,gas_constant,ps0,PI,DEG2RAD

use time_manager_mod,  only : time_type, set_time, set_date, print_time, print_date,    &
                              set_calendar_type, get_calendar_type, get_time, get_date, &
                              operator(-), operator(==)

use utilities_mod,     only : open_file, close_file, find_namelist_in_file, check_namelist_read, &
                              register_module, error_handler, file_exist, E_ERR, E_WARN, E_MSG,  &
                              logfileunit, nmlfileunit, do_output, nc_check, get_unit, do_nml_file, &
                              do_nml_term

use mpi_utilities_mod, only : my_task_id, task_count

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
use location_mod,      only : location_type, get_location, set_location, query_location,         &
                              LocationName, LocationLName, horiz_dist_only,      &
                              vert_is_level, vert_is_pressure, vert_is_height, vert_is_surface,  &
                              vert_is_undef, vert_is_scale_height,                               &
                              VERTISUNDEF, VERTISSURFACE, VERTISLEVEL,                           &
                              VERTISPRESSURE, VERTISHEIGHT, VERTISSCALEHEIGHT, write_location,   &
                              get_close_type, get_close_maxdist_init, get_close_obs_init,        &
                              get_close_obs_destroy,get_dist,loc_get_close_obs => get_close_obs

use xyz_location_mod, only : xyz_location_type, xyz_get_close_maxdist_init,          &
                             xyz_get_close_type, xyz_set_location, xyz_get_location, &
                             xyz_get_close_obs_init, xyz_get_close_obs_destroy,      &
                             xyz_find_nearest

! get_close_maxdist_init, get_close_obs_init, can be modified here (i.e. to add vertical information
! to the initial distance calcs), but will need subroutine pointers like get_close_obs.
! READ THIS SYNTAX as:
!   There's a subroutine in location_mod named 'get_close_obs'.
!   If I want to use that one in this module then refer to it as 'loc_get_close_obs'.
!   If I call 'get_close_obs', then I'll get the one in this module,
!   which does some stuff I need, AND ALSO CALLS 'loc_get_close_obs'

! FIXME
! I've put a copy of solve_quadratic in this model_mod.
! Eventually it should go into a untilities module.
! use utilities_YYY, only : solve_quadratic

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, KIND_PRESSURE,       &
                             KIND_SURFACE_PRESSURE, KIND_TEMPERATURE, KIND_SPECIFIC_HUMIDITY,   &
                             KIND_CLOUD_LIQUID_WATER, KIND_CLOUD_ICE, KIND_CLOUD_FRACTION,      &
                             KIND_GRAV_WAVE_DRAG_EFFIC, KIND_GRAV_WAVE_STRESS_FRACTION,         &
                             KIND_SURFACE_ELEVATION,                                            &
                             KIND_CO, KIND_CO2, KIND_NO, KIND_NO2, KIND_CH4, KIND_NH3, KIND_O3, &
                             get_raw_obs_kind_index, get_raw_obs_kind_name, get_obs_kind_var_type


! Other possibilities (names have changed with various CAM versions):
! Atmos
!    CLOUD:       "Cloud fraction" ;
!    QCWAT:       "q associated with cloud water" ;
!    TCWAT:       "T associated with cloud water" ;
!    CWAT:        "Total Grid box averaged Condensate Amount (liquid + ice)" ;
!   also? LCWAT

! pbl
!    PBLH:        "PBL height" ;
!    QPERT:       "Perturbation specific humidity (eddies in PBL)" ;
!    TPERT:       "Perturbation temperature (eddies in PBL)" ;

! Surface
!    LANDFRAC:    "Fraction of sfc area covered by land" ;
!    LANDM:       "Land ocean transition mask: ocean (0), continent (1), transition (0-1)" ;
!      also LANDM_COSLAT
!    ICEFRAC:     "Fraction of sfc area covered by sea-ice" ;
!    SGH:         "Standard deviation of orography" ;
!    Z0FAC:       "factor relating z0 to sdv of orography" ;
!    TS:          "Surface temperature (radiative)" ;
!    TSOCN:       "Ocean tempertare" ;
!    TSICE:       "Ice temperature" ;
!    TSICERAD:    "Radiatively equivalent ice temperature" ;

! Land/under surface
!    SICTHK:      "Sea ice thickness" ;
!    SNOWHICE:    "Water equivalent snow depth" ;
!    TS1:         "subsoil temperature" ;
!    TS2:         "subsoil temperature" ;
!    TS3:         "subsoil temperature" ;
!    TS4:         "subsoil temperature" ;

! Other fields are not included because they look more CLM oriented.

! Other fields which users may add to the CAM initial files are not listed here.
! Examples are EFGWORO, FRACLDV from the gravity wave drag parameterization study
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use data_structure_mod, only : ensemble_type

use distributed_state_mod

! end of use statements
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

! CAM global/module declarations

implicit none
private

! The first block are the 16 required interfaces.  The following block
! are additional useful interfaces that utility programs can call.
public ::                                                            &
   static_init_model, get_model_size, get_model_time_step,           &
   pert_model_state, get_state_meta_data, model_interpolate_distrib,         &
   nc_write_model_atts, nc_write_model_vars,                         &
   init_conditions, init_time, adv_1step, end_model,                 &
   get_close_maxdist_init, get_close_obs_init, get_close_obs_distrib!, convert_base_obs_location

! Why were these in public?   get_close_maxdist_init, get_close_obs_init, &
! Because assim_model needs them to be there.

public ::                                                   &
   model_type, prog_var_to_vector, vector_to_prog_var,      &
   read_cam_init,                                           &
   init_model_instance, end_model_instance, write_cam_init, &
   write_cam_times

interface get_surface_pressure
   module procedure get_surface_pressure_state, get_surface_pressure_mean
end interface

interface model_heights
   module procedure model_heights_distrib_fwd, model_heights_distrib_mean
end interface

interface interp_lonlat_distrib
   module procedure interp_lonlat_distrib_fwd, interp_lonlat_distrib_mean
end interface

!-----------------------------------------------------------------------
! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
!-----------------------------------------------------------------------

! DART form of ensemble mean, global storage for use by get_close_obs:convert_vert
! Ensemble mean is used so that the same "state" will be used for the height calculations
! on all processors, for all ensemble members.
! This is allocated in static_init_model().
! Distrbuted version does not need the mean.
!real(r8), allocatable :: ens_mean(:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Global storage for describing cam model class
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!-----------------------------------------------------------------------
! Definition of variable types
! Values will be defined in order_state_fields
! All fields will have entries in the TYPE_xD corresponding to their orders
!   in state_names_Xd.  The explicitly named TYPE_s are for convenience
integer :: TYPE_PS = MISSING_I,        &
           TYPE_T = MISSING_I,         &
           TYPE_U = MISSING_I,         &
           TYPE_V = MISSING_I,         &
           TYPE_Q = MISSING_I,         &
           TYPE_CLDICE = MISSING_I,    &
           TYPE_CLDLIQ = MISSING_I,    &
           TYPE_PHIS = MISSING_I,      &
           TYPE_SGH = MISSING_I,       &
           TYPE_PBLH = MISSING_I,      &
           TYPE_TBOT = MISSING_I,      &
           TYPE_TS = MISSING_I,        &
           TYPE_TSOCN = MISSING_I,     &
           TYPE_LCWAT = MISSING_I,     &
           TYPE_QCWAT = MISSING_I,     &
           TYPE_EFGWORO  = MISSING_I,  &
           TYPE_FRACLDV = MISSING_I,   &
           TYPE_CO  = MISSING_I,   &
           TYPE_CO2 = MISSING_I,   &
           TYPE_NO  = MISSING_I,   &
           TYPE_NO2 = MISSING_I,   &
           TYPE_CH4 = MISSING_I,   &
           TYPE_NH3 = MISSING_I,   &
           TYPE_O   = MISSING_I,   &
           TYPE_O3  = MISSING_I

integer, allocatable :: TYPE_0D(:), TYPE_1D(:), TYPE_2D(:), TYPE_3D(:)

!-----------------------------------------------------------------------

! A type for cam model.
! Each variable will be allowed to have different dimensions, even different from
! others of the same rank (i.e. 2d).
! The maximum size for each dimension (for a given rank) will be used to allocate space
! when a model_type variable is initialized.
type model_type
   private
   real(r8), allocatable :: vars_0d(:)
   real(r8), allocatable :: vars_1d(:, :)
   real(r8), allocatable :: vars_2d(:, :, :)
   real(r8), allocatable :: vars_3d(:, :, :, :)
end type model_type

integer :: model_size
! This list of dimensions used to define fields will be ordered as they are on the caminput.nc file.
integer                                   :: num_dims
integer,                      allocatable :: dim_sizes(:)
character(len=NF90_MAX_NAME), allocatable :: dim_names(:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Grid fields
! These structures are used by nc_write_model_atts.
! They are dimensioned in init_grid_1D_instance and filled in read_cam_coord.

! Should this whole type be allocatable, since scalars can be allocated?
! No, not needed.  Deallocating just the allocatable components is enough cleaning.

type grid_1d_type
   private
   character(len=8)       :: label      =  ''
   integer                :: dim_id     =  MISSING_I
   integer                :: length     =  MISSING_I
   real(r8)               :: resolution =  MISSING_R8
   real(r8), allocatable  :: vals(:)    
   integer                :: num_atts   =  MISSING_I
   character(len=NF90_MAX_NAME), allocatable :: atts_names(:) 
   character(len=NF90_MAX_NAME), allocatable :: atts_vals(:)  
end type grid_1d_type

integer :: iii
! integer :: grid_num_0d = 0              ! # of grid scalars to read from file
! P0 now a "coordinate",  and may be removed entirely
! character(len=8),dimension(100) :: grid_names_0d = (/'P0',(' ',iii=1,100)/)

integer          :: grid_num_1d = 12   ! # of 1d grid fields to read from file
character(len=8) :: grid_names_1d(100) = &
          (/ 'lon     ','lat     ','lev     ','gw      ', &
             'hyam    ','hybm    ','hyai    ','hybi    ', &
             'slon    ','slat    ','ilev    ','P0      ', &
            ('        ',iii=1,88 ) /)
! These names should match the grid_names_1d to keep things clear.
! All the possible coordinates (not dimensions) on the caminput.nc file.
type(grid_1d_type), target ::  lon ,lat ,lev ,gw ,hyam ,hybm ,hyai ,hybi, slon ,slat ,ilev, P0
! "Any non-pointer sub-object of an object with the TARGET attribute also has the TARGET attribute."
! So I can point to, e.g., lat%vals.

! Other useful 1D grid arrays (for cubed sphere)
real(r8), allocatable :: lon_rad(:), lat_rad(:)   ! longitude and latitude in radians, used by bearings()

! grid_2d_type ?
! integer :: grid_num_2d = 0              ! # of 2d grid fields to read from file
! ? should phis be in grid_names_2d?
! character (len=8),dimension(100) :: grid_names_2d = (/(' ',iii=1,100)/)


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Namelist variables with default values follow

! output_state_vector = .true.     results in a "state-vector" netCDF file
! output_state_vector = .false.    results in a "prognostic-var" netCDF file
logical :: output_state_vector = .false.

! Files where basic info about model configuration can be found
character(len=128) :: &
   model_config_file = 'caminput.nc',             & ! An example cam initial file.
   cam_phis          = 'cam_phis.nc',             & ! Separate source of PHIS/topography.
   homme_map_file    = 'HommeMapping.nc',         & ! Corners of each cubed sphere cell.
   cs_grid_file      = 'HommeMapping_cs_grid.nc', & ! Relationships among corners/nodes.
   model_version     = '5.0'


! Define location restrictions on which observations are assimilated
! (values are calculated anyway, but istatus is set to 2)
character(len=8) :: vert_coord = 'pressure'            ! or 'log_invP'
real(r8) :: max_obs_lat_degree        = 90.0_r8
real(r8) :: highest_obs_pressure_Pa   = 1000.0_r8
real(r8) :: highest_state_pressure_Pa = 9400.0_r8

! Namelist variables and default values for defining state vector.

integer :: state_num_0d = 0              ! # of scalars fields to read from file
integer :: state_num_1d = 0              ! # of 1d fields to read from file
integer :: state_num_2d = 1              ! # of 2d fields to read from file
integer :: state_num_3d = 4              ! # of 3d fields to read from file

! These can't be allocatable since they are namelist items.
! They have to have a fixed size at compile time.
! Large, arbitrary dimension could be avoided by reading in sizes from a first namelist,
! allocating, setting default values, then get values from second namelist.
! Or, allocate with defaults values, read in namelist, deallocate and reallocate.
integer, parameter :: MAX_STATE_NAMES = 100
character(len=8) :: state_names_0d(MAX_STATE_NAMES)  = ' '
character(len=8) :: state_names_1d(MAX_STATE_NAMES)  = ' '
character(len=8) :: state_names_2d(MAX_STATE_NAMES)  = ' '
character(len=8) :: state_names_3d(MAX_STATE_NAMES)  = ' '

! NOVERT
!         There's a danger of having a mismatch of which_vert_Xd with the state_names_Xd.
!         Should this definition be part of a new structure state_names_Xd, which is parsed
!         into a name and which_vert after being read?  Not for now.

integer :: which_vert_1d(MAX_STATE_NAMES) = MISSING_I
integer :: which_vert_2d(MAX_STATE_NAMES) = MISSING_I
integer :: which_vert_3d(MAX_STATE_NAMES) = MISSING_I


! Is there a way to exclude state_nums from namelist and have those filled in
! the  subroutine which sorts state_names?
! Yes, use two namelists model_nml_1 and model_nml_2 at future date.

! List of fields which this code needs to perturb because they're
! constant valued model parameters and show no spread when start_from_restart = .true.
character(len=8) :: pert_names    (MAX_STATE_NAMES) = '        '
real(r8)         :: pert_sd       (MAX_STATE_NAMES) = MISSING_R8
real(r8)         :: pert_base_vals(MAX_STATE_NAMES) = MISSING_R8

! Special for an experiment.  Specify one string kind e.g KIND_CLOUD_LIQUID and
! observations of that kind will only impact other obs and state vars of that
! same kind.  All other kinds of obs and state vars will not be impacted
! by obs of this kind.  A null string means behave as normal.
character(len=obstypelength) :: impact_only_same_kind = ' '


! Specify shortest time step that the model will support
! This is limited below by CAMs fixed time step but is also impacted
! by numerical stability concerns for repeated restarting in leapfrog.
integer :: Time_step_seconds = 21600, Time_step_days = 0

! set to .true. to get more details about the state vector and the
! CAM fields and sizes in the init code.
logical :: print_details = .false.


namelist /model_nml/ vert_coord, output_state_vector, model_version, cam_phis,        &
                       state_num_0d,   state_num_1d,   state_num_2d,   state_num_3d,  &
                     state_names_0d, state_names_1d, state_names_2d, state_names_3d,  &
                                      which_vert_1d,  which_vert_2d,  which_vert_3d,  &
                     pert_names, pert_sd, pert_base_vals,                             &
                     highest_obs_pressure_Pa, highest_state_pressure_Pa,              &
                     max_obs_lat_degree, Time_step_seconds, Time_step_days,           &
                     impact_only_same_kind, print_details,                            &
                     model_config_file, cs_grid_file, homme_map_file


!---- end of namelist (found in file input.nml) ----
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Derived parameters

! make sure static init code only called once
logical :: module_initialized = .false.

! Variable to keep track of restricting chemistry observations.
integer                      :: impact_kind_index = -1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(time_type) :: Time_step_atmos

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Random sequence and init for pert_model_state
logical                 :: first_pert_call = .true.
type(random_seq_type)   :: random_seq
integer                 :: ens_member = 0
logical                 :: output_task0

! common message string used by many subroutines
character(len=512) :: string1, string2, string3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
integer :: nflds         ! # fields to read

! f_dim_#d are the sizes of the coordinates of each variable as found on caminput File.
! s_dim_#d are the sizes of the coordinates IN THE ORDER ESTABLISHED BY trans_coord for the State.
!     1d fields are easy; whatever the 1 dimension is from f_dim_1d is what s_dim_1d will be.
!     3d fields are easy; always put the dimensions in order lev, lon, lat  (or staggered versions).
!        so s_dim_3d will be filled with 3 values chosen from the 11 dimensions/coords
!        on the caminput.nc
!     2d fields are trickier; if lev is one dimension, it will be in s_dim_2d(1,i)
!                             if lat is one dimension, it will be in s_dim_2d(2,i)
!                                lon can be in either s_dim_2d(1,i) or s_dim_2d(2,i)
!                                depending on what the other dimension is.
! Similarly for the X_dimid_#d arrays, but for dimension id numbers.
! X_dimid_#d first dimension is 1 larger than # spatial dimensions to accomodate time dimension
! on caminit.nc files.
! These are filled in trans_coord
integer              :: coord_order
integer, allocatable :: s_dim_3d(:,:), s_dim_2d(:,:), s_dim_1d(  :),  &
                        f_dim_3d(:,:), f_dim_2d(:,:), f_dim_1d(:,:),  &
                        f_dimid_3d(:,:), f_dimid_2d(:,:), f_dimid_1d(:,:),  &
                        s_dimid_3d(:,:), s_dimid_2d(:,:), s_dimid_1d(  :)
integer :: s_dim_max(3,3)
integer :: f_dim_max(4,3)


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Surface pressures, used by vertical interpolation routines.
!
! I assume that staggered grids (US and VS) are staggered only 1 direction (each),
! so that surface pressure interpolations to get staggered ps use only 2 A-grid ps values.
! The interpolations for columns of heights are more general, but will do a 2 point interp
!     if the staggering is only in one direction.
!
! Need more arrays if any future fields are doubly staggered.

! Should this be a grid_2d_type, specified above?
logical               :: allocate_ps=.true. ! Flag whether to alloc space for ps[_stagr]
real(r8), allocatable :: ps(:, :)           ! surface pressure used to calc P and height profiles.
real(r8), allocatable :: ps_stagr_lon(:, :) ! ps used to calc P profiles & heights on grid staggered
                                            !    East-West (i.e. for VS) relative to ps
real(r8), allocatable :: ps_stagr_lat(:, :) ! ps used to calc P profiles & heights on grid staggered
                                            !    North-South (i.e. for US) relative to ps
real(r8), allocatable :: p(:,:,:)           ! 3D pressure
! height
! Surface potential; used for calculation of geometric heights.
logical               :: alloc_phis=.true.    ! Flag whether to allocate space for phis[_stagr]
real(r8), allocatable :: phis(:, :)           ! surface geopotential
real(r8), allocatable :: phis_stagr_lon(:, :) ! surface geopotential staggered as for ps
real(r8), allocatable :: phis_stagr_lat(:, :) ! surface geopotential staggered as for ps

! I'd need 3 of these; 1 for A-grid and 2 for C-grids, so don't keep model_h laying around.
! real(r8), allocatable :: model_h(:, :, :) ! cartesian heights of model levels

! Columns of pressure and model level heights for use in convert_vert
! HK why are these global?
real(r8), allocatable :: p_col(:), model_h(:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CS Variables holding relationships among cubed sphere nodes.
! Read in and/or set in static_init_model (and routines it calls) at the beginning of an assimilation.
! Used by model_interpolate.
logical :: l_rectang = .true.        ! Flag to tell whether grid is logically rectangular
                                     ! (Eul, FV) or not (cubed-sphere, ...?)
                                     ! Will be set .false. if model_config_file has dimension 'ncol'.
logical :: l_refined = .false.       ! Flag to tell whether grid is a refined mesh or not.

! ne = global metadata from model_config_file, giving the number of 'elements' per face edge of the 'cube'.
! np = # of nodes/edge of each element.  Edges are shared with adjacent elements.
! FIXME? put these in a derived type to prevent accidentally using them as local variables.
integer :: ne, np

! Number of columns, or nodes, in the cubed-sphere grid.
integer :: ncol

! The nominal resolution is (30 degrees/ne), assuming np = 4 (3x3 cells per element).
real(r8) :: coarse_grid

! Dimensions of array 'corners', from homme_map_file.
integer :: ncorners, ncenters

! Maximum number of neighbors a node can have (6 in refined, 4 otherwise)
! Get from namelist, to reduce file & array sizes?
! Or derive from l_refined, after ne is read from caminput.nc?
integer, parameter :: max_neighbors = 6

! Array from homme_map_file.
integer,  allocatable :: corners(:,:)    ! The 4 corners (nodes) of each cell, from HommeMapping.nc

! 5 arrays from cs_grid_file
integer,  allocatable :: num_nghbrs(:)   ! Number of neighbors of each node/column in the cubed sphere grid.
integer,  allocatable :: centers(:,:)    ! The names of the cells that use each node as a corner.
real(r8), allocatable :: a(:,:,:)        ! Coefficients of mapping from planar to unit square space for 'x'
real(r8), allocatable :: b(:,:,:)        ! Coefficients of mapping from planar to unit square space for 'y'
real(r8), allocatable :: x_ax_bearings(:,:)  ! The directions from each node to its neighbors,
                                         ! measured from the vector pointing north.  (-PI <= bearing <= PI)

! Locations of cubed sphere nodes, in DART's location_type format.
type(location_type), allocatable :: cs_locs(:)

! Location of cubed sphere nodes, in cartesian coordinates
type(xyz_location_type), allocatable :: cs_locs_xyz(:)
type(xyz_get_close_type)             :: cs_gc_xyz

! Structure containing grid point locations, etc.,
! defined in static_init_mod after reading in CS lons, lats, and levels,
! needed in model_interpolate:interp_cubed_sphere.
type(get_close_type) :: cs_gc

! Array of KINDs of cubed sphere grid points.
! As of 2014-3-28 this is only used by location_mod, which doesn't actually use it.
integer, allocatable :: cs_kinds(:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Array 'cflds' is filled with simple loops over state_names_xxx.
! I could replace that with code which orders namelist input field names
! into cflds, regardless of their original order, and tallies how many of each.
! Is there a way to exclude state_nums from namelist and have those filled in
! the same subroutine?

character(len=8), allocatable :: cflds(:)

! Attribute values for the fields which comprise the state vector.
! These are filled by nc_read_model_atts.
character(len=nf90_max_name), allocatable :: state_long_names(:)
character(len=nf90_max_name), allocatable :: state_units(:)

! array for the linking of obs_kinds(KIND_) to model field TYPE_s
! It's filled in map_kinds
! The max size of KIND_ should come from obs_kind_mod
! These should be dimensioned the same size as the total of state_names_Nd.
integer :: dart_to_cam_types(300) = MISSING_I
integer :: cam_to_dart_kinds(300) = MISSING_I
!
!-----------------------------------------------------------------------
! These are calculated from highest_obs_pressure_Pa
integer             :: highest_obs_level    = MISSING_I
real(r8)            :: highest_obs_height_m = MISSING_R8
! Better damping
! Variables to deal with CAM's damping of the top levels of the model.
! These are set in static_init_model and used in get_close_obs.
real(r8)            :: highest_state_scale_h  = MISSING_R8
real(r8)            :: model_top              = MISSING_R8 
real(r8)            :: damp_wght              = MISSING_R8
type(location_type) :: highest_state_loc, model_top_loc

!-----------------------------------------------------------------------

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

contains

!#######################################################################

! static_init_model section

!-----------------------------------------------------------------------
!>
!> Static_init_model does many things which must be done once at the beginning 
!> of the use of model_mod:
!>   + set the calendar and time variables,
!>   + read, check and archive the model_nml namelist,
!>   + set some output level variables,
!>   + set the state vector size,
!>   + read coordinate variables from the CAM initial file,
!>   + read the model topography
!>   + read and/or generate cubed sphere grid arrays if CAM-SE is being used,
!>   + make the connection between DART KINDs and local model TYPEs

subroutine static_init_model()

! Initializes class data for CAM model (all the stuff that needs to be done once).
! For now, does this by reading info from a fixed name netcdf file.

integer :: iunit, io, i, nc_file_ID
integer :: max_levs, ierr

! only execute this code once
if (module_initialized) return

! Make sure we only come through here once.
module_initialized = .true.

! Register the module
call register_module(source, revision, revdate)

! setting calendar type
! calendar types listed in time_manager_mod.f90
! this information is NOT passed to CAM; it must be set in the CAM namelist
call set_calendar_type('GREGORIAN')

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")
call verify_namelist()

! Stand-alone CAM (not CESM) uses the first block of the if statement.
! Set the printed output logical variable to reduce printed output;

if (file_exist('element')) then
   iunit = get_unit()
   open(unit=iunit, file='element', form='formatted')
   read(iunit,*) ens_member
   close(iunit)
   output_task0 = .false.
   if (ens_member == 1) output_task0 = .true.
else
   output_task0 = do_output()
endif

! Record the namelist values
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(    *      , nml=model_nml)

! Set the model minimum time step from the namelist seconds and days input
Time_step_atmos = set_time(Time_step_seconds, Time_step_days)
if (print_details .and. output_task0) call print_time(Time_step_atmos)

! Open CAM 'initial' file to read dimensions and coordinates of fields.
call nc_check(nf90_open(path=trim(model_config_file), mode=nf90_nowrite, ncid=nc_file_ID), &
             'static_init_model', 'opening '//trim(model_config_file))

! Get sizes of dimensions/coordinates from netcdf and put in global storage.
! Also may change l_rectang to .false.
call read_cam_init_size(nc_file_ID)

! Compute overall model size and put in global storage
! s_dim_#d come from read_cam_init_size/trans_coord, and are in global storage
model_size = state_num_0d
do i=1,state_num_1d
   model_size = model_size + s_dim_1d(i)
enddo
do i=1,state_num_2d
   model_size = model_size + s_dim_2d(1,i) * s_dim_2d(2,i)
enddo
do i=1,state_num_3d
   model_size = model_size + s_dim_3d(1,i) * s_dim_3d(2,i) * s_dim_3d(3,i)
enddo
if (output_task0) then
   write(string1, '(A,I9)') 'CAM state vector size: ', model_size
   call error_handler(E_MSG, 'static_init_model', string1)
endif

!HK not in the RMA branch
!allocate(ens_mean(model_size))

! Allocate space for global coordinate arrays and read them in.
! There's a query of caminput.nc within read_cam_coord for the existence of the field.
! The second argument is a grid_1d_type structure
! CS; ncol is a dimension, but there's no coordinate variable of the same name.
!     There are lat and lon arrays for the ncol grid points.
call read_cam_coord(nc_file_ID, 'lon', lon)
call read_cam_coord(nc_file_ID, 'lat', lat)
call read_cam_coord(nc_file_ID, 'lev', lev)
call read_cam_coord(nc_file_ID, 'ilev', ilev)
call read_cam_coord(nc_file_ID, 'gw', gw)
call read_cam_coord(nc_file_ID, 'slon', slon)
call read_cam_coord(nc_file_ID, 'slat', slat)

! read hybrid vert coord coefs
call read_cam_coord(nc_file_ID, 'hyai', hyai)
call read_cam_coord(nc_file_ID, 'hybi', hybi)
call read_cam_coord(nc_file_ID, 'hyam', hyam)
call read_cam_coord(nc_file_ID, 'hybm', hybm)

! It's a scalar, but I can put it into the same coord structure as previous fields.
! It's length will be 1
call read_cam_coord(ncfileid, 'P0', P0)    ! thats a p-zero

!------------------------------------------------------------------------
! Better damping algorithm for state variables near/in the CAM damped levels
! at the top of the model.  
! See get_close_obs and models/cam/doc/highest_state_p_Pa.pptx for details.
! This section must come after the definition of P0 and hyai.
if (vert_coord == 'pressure') then
   ! CAM's model_top is 1/2 level above the highest state variable level, so
   ! hyai instead of hyam.
   ! P0 is in Pa.
   model_top = hyai%vals(1)*P0%vals(1)
   ! The (lon,lat) here must match the ones in the definition of vert_only_loc in get_close_obs.
   ! FIXME; is this hard-coding OK?
   highest_state_loc = set_location(1.0_r8,1.0_r8,highest_state_pressure_Pa,VERTISPRESSURE)
   model_top_loc     = set_location(1.0_r8,1.0_r8,model_top,                VERTISPRESSURE)
   ! damp_wght must be in the same units (dist = radians) as the distances in get_close_obs.
   if (highest_state_pressure_Pa /= model_top) then
      damp_wght = 1.0_r8/get_dist(highest_state_loc,model_top_loc,no_vert=.false.)
   end if
elseif (vert_coord == 'log_invP') then
   highest_state_scale_h = scale_height(p_surface=P0%vals(1), p_above=highest_state_pressure_Pa)
   model_top             = scale_height(p_surface=P0%vals(1), p_above=(hyai%vals(1)*P0%vals(1)) )
   highest_state_loc = set_location(1.0_r8,1.0_r8,highest_state_scale_h,VERTISSCALEHEIGHT)
   model_top_loc     = set_location(1.0_r8,1.0_r8,model_top,            VERTISSCALEHEIGHT)
   if (highest_state_scale_h /= model_top) then
      damp_wght = 1.0_r8/get_dist(highest_state_loc,model_top_loc,no_vert=.false.)
   end if
else
   write(string1, '(A,A)') 'Somehow vert_coord /= {pressure,log_invP}: ', vert_coord
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! # fields to read
nflds = state_num_0d + state_num_1d + state_num_2d + state_num_3d
if (print_details) then
   write(string1, '(A,I3,A,4I3)') '# of fields in state vector =  ', nflds, &
        ' = sum of ', state_num_0d ,state_num_1d ,state_num_2d ,state_num_3d
   call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)
endif

! Order the state vector parts into cflds.
allocate(cflds(nflds))
call order_state_fields()

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Get field attributes needed by nc_write_model_atts from caminput.nc.
allocate(state_long_names(nflds), state_units(nflds))    
call nc_read_model_atts(nc_file_ID, 'long_name', state_long_names)
call nc_read_model_atts(nc_file_ID, 'units', state_units)

if (.not. l_rectang) then
   ! Read some attributes from the cubed sphere model_config_file.
   ! ne is the number of elements/cube edge.  Usually 0 for refined grids.
   ! np is the number of nodes/element edge (shared with adjacent element.
   call nc_read_global_int_att(nc_file_ID, 'ne', ne)
   call nc_read_global_int_att(nc_file_ID, 'np', np)

   ! Calculate the nominal resolution of the (coarse) grid,
   ! for use by model_interpolate's call to get_close_obs.
   if (ne == 0) then
      ! Refined mesh; assume the coarsest grid is the default '1-degree'.
      ! Need factor of 1.5 to make sure that there are at least 2 nodes 'close' to any location.
      ! There seems to be a tricky interplay between the lon-lat boxes used in the quick search
      ! for potentially close nodes, and the cubed sphere grid, so that a coarse_grid of only
      ! slightly more than 1.0 degrees can yield 0 close nodes.
      coarse_grid = 1.2_r8 * DEG2RAD
      l_refined = .true.
   else
      ! Standard cubed sphere; there are 3x num_elements/face_edge x 4 nodes
      ! around the equator.  ne = 30 -> 3x4x30 = 360 nodes -> '1-degree'
      ! Yielded a location with only 1 close ob, but need 2.
      ! coarse_grid = (30.01_r8/ne) * DEG2RAD
      coarse_grid = 1.2_r8*(30.0_r8/ne) * DEG2RAD
   endif
   if (print_details) then
      write(string1,*),'Cubed sphere coarse_grid resolution (rad) used in cs_gc definition = ',&
                      coarse_grid,' because ne = ',ne
      call error_handler(E_MSG, 'static_init_model', string1,source,revision,revdate)
   endif

   ! Fill cs_gc for use by model_mod.  Inputs and outputs are in global storage.
   ! In particular, ncol must be defined before this call.
   ncol = dim_sizes(find_name('ncol',dim_names))
   call fill_gc()

   ! Fill arrays that are useful for bearings and distances.
   allocate(lon_rad(ncol), lat_rad(ncol))
   do i=1,ncol
      lon_rad(i) = lon%vals(i)*DEG2RAD
      lat_rad(i) = lat%vals(i)*DEG2RAD
   enddo

endif

call nc_check(nf90_close(nc_file_ID), &
              'static_init_model', 'closing '//trim(model_config_file))

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! height
! Get dimensions and surface geopotential from a new netcdf file and test for consistency.
! Open file and read PHIS from it.
! Allocate global variables which will be used in vertical interpolations
! Check for pressures on vertically staggered grid, as well as standard grid.

call read_cam_2Dreal(cam_phis, 'PHIS')

max_levs = lev%length
if (ilev%label /= '') max_levs = max(ilev%length, lev%length)
allocate(p_col(max_levs), model_h(max_levs))

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CS Cubed sphere grid data.
! Read in or create a file containing the relationships among cubed sphere nodes,
! such as neighbors, centers, and bearings, which will be used to identify the cell
! which contains an observation.
! Fields will be stored in global storage.
! Write the cubed sphere grid arrays to a new NetCDF file.

if (.not. l_rectang) then
   if (file_exist(cs_grid_file)) then
      call nc_read_cs_grid_file()
   elseif (file_exist(homme_map_file)) then
      call create_cs_grid_arrays()
      if (my_task_id() == 0) call nc_write_cs_grid_file( cs_grid_file, homme_map_file )
   else
      write(string1, *)'No cs_grid_file "',trim(cs_grid_file), &
                    '" nor homme_map_file "',trim(homme_map_file),'"'
      call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
   endif

endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Fills arrays for the linking of obs_kinds (KIND_) to model field TYPE_s
call map_kinds()

! If restricting impact of a particular kind to only obs and state vars
! of the same kind, look up and set the kind index.
if (len_trim(impact_only_same_kind) > 0) then
   impact_kind_index = get_raw_obs_kind_index(impact_only_same_kind)
endif

! This validates the namelist value and sets the module global value highest_obs_level.
call set_highest_obs_limit()

end subroutine static_init_model

!-----------------------------------------------------------------------

subroutine verify_namelist()

!  FIXME; PS must always be in the state vector;
!         always add PS in to state vector (if missing)
!         In the future we may want to let people not update PS in filter, or ...?

integer :: i
logical :: ps_present = .false.
logical :: mismatch_which   = .false.
logical :: mismatch_size    = .false.

if (state_num_0d > 0) then
   if (state_names_0d(state_num_0d)   == ' ' .or. &
       state_names_0d(state_num_0d+1) /= ' ') mismatch_size = .true.
endif

if (state_num_1d > 0) then
   if (state_names_1d(state_num_1d)   == ' ' .or. &
       state_names_1d(state_num_1d+1) /= ' ') mismatch_size = .true.
   if ( which_vert_1d(state_num_1d)   == MISSING_I .or. &
        which_vert_1d(state_num_1d+1) /= MISSING_I) mismatch_which = .true.
endif

if (state_num_2d > 0) then
   if (state_names_2d(state_num_2d) == ' ' .or. &
       state_names_2d(state_num_2d+1) /= ' ') mismatch_size = .true.
   if ( which_vert_2d(state_num_2d) == MISSING_I .or. &
        which_vert_2d(state_num_2d+1) /= MISSING_I) mismatch_which = .true.
endif

if (state_num_3d > 0) then
   if (state_names_3d(state_num_3d) == ' ' .or. &
       state_names_3d(state_num_3d+1) /= ' ') mismatch_size = .true.
   if ( which_vert_3d(state_num_3d) == MISSING_I .or. &
        which_vert_3d(state_num_3d+1) /= MISSING_I) mismatch_which = .true.
endif

if (mismatch_size) then
   write(string1,*) 'Mismatch between state_num_#d and state_names_#d in model_nml'
   call error_handler(E_ERR,'verify_namelist',string1,source,revision,revdate)
endif

if (mismatch_which) then
   write(string1,*) 'Mismatch between state_num_#d and which_vert_#d in model_nml'
   call error_handler(E_ERR,'verify_namelist',string1,source,revision,revdate)
endif

mismatch_which = .false.
do i=1,max(state_num_1d,state_num_2d,state_num_3d)
   if (which_vert_1d(i) > 1 ) mismatch_which = .true.
   if (which_vert_2d(i) > 1 ) mismatch_which = .true.
   if (which_vert_3d(i) > 1 ) mismatch_which = .true.

   ! PS can't be 0d or 3d.
   if (state_names_1d(i) == 'PS') ps_present = .true.
   if (state_names_2d(i) == 'PS') ps_present = .true.
enddo

if (mismatch_which) then
   write(string1,*) 'The CAM model state is defined on levels and the surface. ', &
                    '   which_vert_#d must be -2, -1, or 1 for each state variable.'
   call error_handler(E_ERR,'verify_namelist',string1,source,revision,revdate)
endif


if (.not. ps_present) then
   write(string1,*) '"PS" (surface pressure) must be one of the state variables, but was not found'
   call error_handler(E_ERR,'verify_namelist',string1,source,revision,revdate)
endif

if (vert_coord /= 'pressure' .and. vert_coord /= 'log_invP') then
   write(string1,*) 'vert_coord must be "pressure" or "log_invP"'
   call error_handler(E_ERR,'verify_namelist',string1,source,revision,revdate)
endif

end subroutine verify_namelist

!-----------------------------------------------------------------------

subroutine read_cam_init_size(nc_file_ID)

! Gets the number, names, and sizes of field dimensions from a CAM init netcdf file
! in file_name (regardless of dynamical core).
! Called by static_init_model (only).

integer,  intent(in)  :: nc_file_ID

integer :: i,j

if (.not. module_initialized) call static_init_model()

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! learn how many dimensions are defined in this file.
call nc_check(nf90_inquire(nc_file_ID, num_dims), 'read_cam_init_size', 'inquire num_dims')

allocate(dim_names(num_dims), dim_sizes(num_dims))

! Cycle through dimension ids until there aren't any more.
! dimension ids are sequential integers on the NetCDF file.
do i = 1,num_dims
   call nc_check(nf90_inquire_dimension(nc_file_ID, i, dim_names(i), dim_sizes(i)), &
                 'read_cam_init_size', 'inquire for '//trim(dim_names(i)))
   if (print_details .and. output_task0) then
      write(string1,*) 'Dims info = ',i, trim(dim_names(i)), dim_sizes(i)
      call error_handler(E_MSG, 'read_cam_init_size', string1,source,revision,revdate)
   endif

   if (dim_names(i) == 'ncol') l_rectang = .false.
enddo

! Find and store shapes of all the state vector fields.  Grouped by rank of fields into
! separate s_dim_RANKd arrays.
! Fields with the same rank can have different shapes and still be handled; efgworo(lat,lon) and
! frac(lat,lev) will both have their shapes stored in X_dim_2d.
!   CS; remove the following complication?
! Also keep track of whether the init file is old (lon,lev,lat) or new (lon,lat,lev).

call trans_coord(nc_file_ID)

! The arrays into which CAM fields are put are dimensioned by the largest values of
! the sizes of the dimensions listed in Y_dim_RANKd, Y=[sf], RANK=[1-3] .
! The second dimension denotes the rank of the array for which the first dim
! gives the max size(s).
if (state_num_1d > 0) then
   f_dim_max(1:2, 1) = maxval(f_dim_1d, dim=2)   ! gets the max value of f_dim_1d (1:2, :)
   s_dim_max(1  , 1) = maxval(s_dim_1d, dim=1)   ! gets the max value of s_dim_1d (:)
else
   f_dim_max(1:2, 1) = 0
   s_dim_max(1  , 1) = 0
endif

if (state_num_2d > 0) then
   f_dim_max(1:3, 2) = maxval(f_dim_2d, dim=2)   ! gets the max values of f_dim_2d (1:3, :)
   s_dim_max(1:2, 2) = maxval(s_dim_2d, dim=2)   ! gets the max values of s_dim_2d (1:2, :)
else
   f_dim_max(1:3, 2) = 0
   s_dim_max(1:2, 2) = 0
endif

if (state_num_3d > 0) then
   f_dim_max(1:4, 3) = maxval(f_dim_3d, dim=2)   ! gets the max values of f_dim_3d (1:4, :)
   s_dim_max(1:3, 3) = maxval(s_dim_3d, dim=2)   ! gets the max values of s_dim_3d (1:3, :)
else
   f_dim_max(1:4, 3) = 0
   s_dim_max(1:3, 3) = 0
endif

if (print_details .and. output_task0 ) then
   if (state_num_1d > 0) then
      write(string1,*) 's_dim_1d = ',s_dim_1d
      write(string2,*) (s_dim_max(i,1),i=1,3)
      call error_handler(E_MSG, 'read_cam_init_size', string1,source,revision,revdate, text2=string2)
   endif

   do i=1,2
      write(string1,*) 's_dim_2d = ',(s_dim_2d(i,j),j=1,state_num_2d),'s_dim_max = ',s_dim_max(i,2)
      call error_handler(E_MSG, 'read_cam_init_size', string1,source,revision,revdate)
   enddo

   do i=1,3
      write(string1,'(/A,(10I4))') 's_dim_3d = ',(s_dim_3d(i,j),j=1,state_num_3d)
      write(string2,'(A,(10I4))') 's_dim_max = ',s_dim_max(i,3)
      call error_handler(E_MSG, 'read_cam_init', string1,source,revision,revdate, text2=string2)
      write(string1,'(A,(10I4))') 'f_dim_3d = ',(f_dim_3d(i,j),j=1,state_num_3d)
      write(string2,'(A,(10I4))') 'f_dim_max = ',f_dim_max(i,3)
      call error_handler(E_MSG, 'read_cam_init', string1,source,revision,revdate, text2=string2)
   enddo
endif

end subroutine read_cam_init_size

!-----------------------------------------------------------------------

subroutine trans_coord(nc_file_ID)

! CS; simplify this by not considering old CAM coord orders.
! Figure out which coordinates are lon, lat, lev, based on CAM version
! from the namelist, which has form #.#[.#[.#]].

integer,            intent(in) :: nc_file_ID

! local workspace
character(len=4) :: form_version = '(I0)'
character(len=4) :: char_version
integer          :: part, nchars, tot_chars, i, j, k, varid, next
integer          :: int_version(4)

int_version = (/ (0,i=1,4) /)

! Choose order of coordinates based on CAM version
part = 1
nchars=0
char_version = ' '
tot_chars = len_trim(model_version)
do i=1,tot_chars+1
   if ( i == tot_chars+1 .or.  model_version(i:i) == '.' ) then
      write(form_version(3:3),'(I1)') nchars
      read(char_version,form_version) int_version(part)
      part = part + 1
      nchars = 0
      char_version = ' '
   else
      nchars = nchars + 1
      char_version(nchars:nchars) = model_version(i:i)
   endif
enddo
if (output_task0) then
   if (print_details) then
      write(string1,'(A,A10,4(I3,2X))') 'model_version, version(1:4) = ' &
                                  ,model_version,(int_version(i),i=1,4)
      call error_handler(E_MSG, 'trans_coord', string1,source,revision,revdate)
   else
      call error_handler(E_MSG, 'trans_coord', 'CAM model version: '//trim(model_version))
   endif
endif

! assume cam3.0.7 (modern) format to start
coord_order = 2
if (int_version(1) < 3) then
   coord_order = 1
elseif (int_version(1) == 3 .and. int_version(2) == 0 .and. int_version(3) < 3) then
   coord_order = 1
endif

! Cycle through each field's dimension IDs.
! Pick the dimensions needed out of dim_sizes, using the dimension names in dim_names.
! Fill the state dimids according to the order model_mod wants to see.  (lev, lon, lat).

! 3D is easy; lev, lon, lat are always the coordinates,
! and model_mod always wants them in that order.

if (state_num_3d > 0) then
   allocate(s_dim_3d(3,state_num_3d), s_dimid_3d(3,state_num_3d), &
            f_dim_3d(4,state_num_3d), f_dimid_3d(4,state_num_3d))
   s_dim_3d   = 0
   s_dimid_3d = 0
   f_dim_3d   = 0
   f_dimid_3d = 0
endif

do i = 1,state_num_3d
   ! Get variable id for a  3d field
   call nc_check(nf90_inq_varid(nc_file_ID, state_names_3d(i), varid), &
                 'trans_coord', 'inq_varid '//trim(state_names_3d(i)))
   ! Get dimension ids for the dimensions of the field
   call nc_check(nf90_inquire_variable(nc_file_ID, varid, dimids=f_dimid_3d(1:4,i)), &
                 'trans_coord', 'inquire_variable'//trim(state_names_3d(i)))

   Alldim3: do j = 1,4                          ! time and 3 space
      k = f_dimid_3d(j,i)                       ! shorthand; the dimid of this fields current dim
      f_dim_3d(j,i) = dim_sizes(k)
      ! Put the dimensions we want in the state field positions we want.
      if (dim_names(k) == 'lev' .or. dim_names(k) == 'ilev') then
         s_dim_3d  (1,i) = dim_sizes(k)
         s_dimid_3d(1,i) = k
!         s_dimid_3d(1,i) = f_dimid_3d(j,i)
      elseif (dim_names(k) == 'lon' .or. dim_names(k) == 'slon') then
         s_dim_3d  (2,i) = dim_sizes(k)
         s_dimid_3d(2,i) = k
      elseif (dim_names(k) == 'lat' .or. dim_names(k) == 'slat') then
         s_dim_3d  (3,i) = dim_sizes(k)
         s_dimid_3d(3,i) = k
      endif
!            cycle Alldim3
   enddo Alldim3
   if (   s_dim_3d(1,i) == 0 .or.  s_dim_3d(2,i) == 0 .or.  s_dim_3d(3,i) == 0 ) then
      call error_handler(E_ERR, 'trans_coord', &
          'num_[lons,lats,levs] was not assigned and = 0' , source, revision, revdate)
   endif
enddo

! Fill dimids according to the order model_mod wants to see.
! 2d is trickier (except for CS); 2 of (lev, lon, lat) in that order.

if (state_num_2d > 0) then
   allocate(s_dim_2d(2,state_num_2d), s_dimid_2d(2,state_num_2d), &
            f_dim_2d(3,state_num_2d), f_dimid_2d(3,state_num_2d))

   s_dim_2d = 0;  s_dimid_2d  = 0;
   f_dim_2d = 0;  f_dimid_2d  = 0;
endif

do i = 1,state_num_2d
   call nc_check(nf90_inq_varid(nc_file_ID, state_names_2d(i), varid), &
              'trans_coord', 'inq_varid '//trim(state_names_2d(i)))
   call nc_check(nf90_inquire_variable(nc_file_ID, varid, dimids=f_dimid_2d(1:3,i)), &
              'trans_coord', 'inquire_variable '//trim(state_names_2d(i)))

   ! Extract spatial dimids from the fields dimids
   next = 1
   Alldim2: do j = 1,3      ! time and 2 space
      k = f_dimid_2d(j,i)
      f_dim_2d(j,i) = dim_sizes(k)
      if (dim_names(k) == 'lev' .or. dim_names(k) == 'ilev' ) then
         ! Move whatever may have come in first to the second place.
         ! s_dim_2d was initialized to 0 at the start of the module.
         s_dim_2d  (2,i) = s_dim_2d  (1,i)
         s_dimid_2d(2,i) = s_dimid_2d(1,i)
         ! Put lev in the first dimension
         s_dim_2d  (1,i) = dim_sizes(k)
         s_dimid_2d(1,i) = k
!         s_dimid_2d(1,i) = f_dimid_2d(j,i)
         next = 2
      elseif (dim_names(k) == 'lon' .or. dim_names(k) == 'slon') then
         ! longitude always comes first on CAM initial files.
         ! Otherwise, I'll need a test like for levs, but more complicated.
         s_dim_2d  (next,i) = dim_sizes(k)
         s_dimid_2d(next,i) = k
         next = 2
      elseif (dim_names(k) == 'ncol' ) then
         s_dim_2d  (next,i) = dim_sizes(k)
         s_dimid_2d(next,i) = k
         next = 2
      elseif (dim_names(k) == 'lat' .or. dim_names(k) == 'slat' ) then
         s_dim_2d  (next,i) = dim_sizes(k)
         s_dimid_2d(next,i) = k
      endif
   enddo Alldim2
   if (   s_dim_2d(1,i) == 0 .or.  s_dim_2d(2,i) == 0 ) then
      call error_handler(E_ERR, 'trans_coord', &
          'num_[lons,lats,levs,ncol] was not assigned and = 0' , source, revision, revdate)
   endif
enddo

! 1d fields are easy.

if (state_num_1d > 0) then
   allocate(s_dim_1d(  state_num_1d), s_dimid_1d(  state_num_1d))
   allocate(f_dim_1d(2,state_num_1d), f_dimid_1d(2,state_num_1d))
   s_dim_1d = 0;   s_dimid_1d = 0;
   f_dim_1d = 0;   f_dimid_1d = 0;
endif

do i = 1,state_num_1d
   call nc_check(nf90_inq_varid       (nc_file_ID, state_names_1d(i), varid), &
              'trans_coord', 'inq_varid '//trim(state_names_1d(i)))
   call nc_check(nf90_inquire_variable(nc_file_ID, varid, dimids=f_dimid_1d(1:2,i)), &
              'trans_coord', 'inq_varid '//trim(state_names_1d(i)))

   Alldim1: do j = 1,2       ! time and 1 space
      k = f_dimid_1d(j,i)
      f_dim_1d(j,i) = dim_sizes(k)
      if (dim_names(k) == 'lon' .or. dim_names(k) == 'slon' .or. &
          dim_names(k) == 'ncol'.or. &
          dim_names(k) == 'lat' .or. dim_names(k) == 'slat' .or. &
          dim_names(k) == 'lev' .or. dim_names(k) == 'ilev' ) then
         s_dim_1d(i) = dim_sizes(k)
         s_dimid_1d(i) = k
!         s_dimid_1d(i) = f_dimid_1d(j,i)
      endif
   enddo Alldim1

   if ( s_dim_1d(i) == 0 ) then
      write(string1, '(A,I3,A)') ' state 1d dimension(',i,') was not assigned and = 0'
      call error_handler(E_ERR, 'trans_coord',trim(string1), source, revision, revdate)
   endif
enddo

end subroutine trans_coord

!-----------------------------------------------------------------------

subroutine read_cam_2Dreal(file_name, cfield)

! Subroutine to read in a 2D/horizontal CAM field, such as PHIS.
! Handles both logically rectangular arrays (FV and Eul) and irregular
! (SE-CAM/cubed-sphere).

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

character(len=*), intent(in)  :: file_name
character(len=*), intent(in)  :: cfield

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer :: nc_file_ID, nc_var_ID  ! NetCDF variables
integer :: field_dim_IDs(3)       ! Array of dimension IDs for cfield
                                  ! (2 space (FV) and time dimension (CAM .h0. files).
integer :: i_dim1, i_dim2         ! Variables to reference the dimension(s) of cfield
integer :: num_dim1, num_dim2     ! NetCDF file variable dimension sizes, for comparison to file_name's
integer :: slon_index, slat_index, lon_index, lat_index !indices of [s]lon and [s]lat
                                                        ! within the list of dimensions
integer :: n,m
character(len=NF90_MAX_NAME) :: name_dim1,name_dim2    ! Names of dimensions of cfield
real(r8), allocatable         :: var(:,:)               ! Temp array used by nc_get_var

field_dim_IDs = MISSING_I    ! Array of dimension IDs for cfield

if (file_name == cam_phis .and. .not.file_exist(trim(file_name))) then
   write(string1,'(2A)') trim(file_name),  &
        ' is missing; trying to find a CAM history file (h0) to provide '//cfield
   call error_handler(E_WARN, 'read_cam_2Dreal', trim(string1), source, revision, revdate)
endif

! Open the file and get dimension information.
if (file_exist(trim(file_name))) then
   call nc_check(nf90_open(path=trim(file_name), mode=nf90_nowrite, ncid=nc_file_ID), &
              'static_init_model:read_cam_2Dreal', 'opening '//trim(file_name))
   if (print_details .and. output_task0) then
      write(string1, *) 'file_name for ',cfield,' is ', trim(file_name)
      call error_handler(E_MSG, 'read_cam_2Dreal', string1,source,revision,revdate)
   endif

   ! get field id
   call nc_check(nf90_inq_varid(nc_file_ID, trim(cfield), nc_var_ID), &
              'read_cam_2Dreal', 'inq_varid: '//cfield)

   ! get dimension 'id's
   call nc_check(nf90_inquire_variable(nc_file_ID, nc_var_ID, dimids = field_dim_IDs), &
              'read_cam_2Dreal', 'inquire_variable: '//cfield)

   ! get dimension sizes
   ! The first spatial dimension is always present.
   call nc_check(nf90_inquire_dimension(nc_file_ID, field_dim_IDs(1), name_dim1, num_dim1 ), &
                 'read_cam_2Dreal', 'inquire_dimension: '//name_dim1)
   if (field_dim_IDs(2) == MISSING_I)  then
      num_dim2 = 1
      name_dim2 = 'no2nd_dim_'
   else
      call nc_check(nf90_inquire_dimension(nc_file_ID, field_dim_IDs(2), name_dim2, num_dim2 ), &
                    'read_cam_2Dreal', 'inquire_dimension: '//name_dim2)
   endif

   ! Check for consistent dimensions between initial file and cam_phis file.
   if (file_name == cam_phis) then
      i_dim1 = dim_sizes(find_name(name_dim1,dim_names))
      if (num_dim1 /= i_dim1) then
         write(string1,'(A,2I8,A)') 'i_dim1, num_dim1, name_dim1 =' ,&
                                       i_dim1, num_dim1, trim(name_dim1)
         call error_handler(E_MSG, 'read_cam_2Dreal', trim(string1), source, revision, revdate)
         write(string1,'(A,4I12)') 'horizontal dimensions mismatch of initial files and topog ' &
               ,i_dim1, num_dim1
         call error_handler(E_ERR, 'read_cam_2Dreal', trim(string1), source, revision, revdate)
      endif

      if (field_dim_IDs(2) /= MISSING_I) then
         i_dim2 = dim_sizes(find_name(name_dim2,dim_names))
         if ( num_dim2 /= i_dim2 ) then
            write(string1,'(A,2I8,A)') 'i_dim2, num_dim2, name_dim2 =', &
                                          i_dim2, num_dim2, trim(name_dim2)
            call error_handler(E_MSG, 'read_cam_2Dreal', trim(string1), source, revision, revdate)
            write(string1,'(A,4I12)') 'horizontal dimensions mismatch of initial files and topog ', &
                  i_dim2, num_dim2
            call error_handler(E_ERR, 'read_cam_2Dreal', trim(string1), source, revision, revdate)
         endif
      endif
   endif
else
   write(string1,'(2A)') trim(file_name),  &
        ' is missing; I do not know how to find it.'
   call error_handler(E_ERR, 'read_cam_2Dreal', trim(string1), source, revision, revdate)
endif

! Allocate local arrays, based on size of this variable on the file.
allocate(var(num_dim1, num_dim2))

! Read surface geopotential from cam_phis for use in vertical interpolation in height.
! Coordinate order not affected by CAM model version.
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, var, start=(/ 1, 1 /), &
              count=(/ num_dim1, num_dim2 /)), 'read_cam_2Dreal', trim(cfield))

! assign values to phis grids for use by the rest of the module.
if (cfield == 'PHIS') then

   if (alloc_phis) allocate(phis(num_dim1, num_dim2))
   ! Don't want to set alloc_phis = false yet; there may be staggered phis to set.
   phis(1:num_dim1,1:num_dim2) = var

   ! If needed, generate phis on the staggered grids.
   slon_index = find_name('slon',dim_names)
   slat_index = find_name('slat',dim_names)
   lat_index  = find_name('lat',dim_names)
   lon_index  = find_name('lon',dim_names)

   ! CS these sections will be skipped for grids with no staggered grids, e.g. cubed sphere
   if (slon_index /= 0) then
      if (alloc_phis) allocate(phis_stagr_lon(dim_sizes(slon_index), dim_sizes( lat_index)))
      do n=1,dim_sizes( lat_index)
         ! CS Would a better interpolation (e.g. splines) help anything in a meaningful way?
         ! ? Better interpolation using interp_cubed_sphere?  But it's 'bilinear', so no differenc?
         phis_stagr_lon(1,n) = 0.5_r8 * (phis(1,n) + phis(dim_sizes(lon_index),n))
         do m=2,dim_sizes(slon_index)
            phis_stagr_lon(m,n) = 0.5_r8 * (phis(m-1,n) + phis(m,n))
         enddo
      enddo
   endif

   if (slat_index /= 0) then
      if (alloc_phis) allocate(phis_stagr_lat(dim_sizes( lon_index), dim_sizes(slat_index)))
      do n=1,dim_sizes(slat_index)
         do m=1,dim_sizes( lon_index)
            phis_stagr_lat(m,n) = 0.5_r8 * (phis(m,n) + phis(m,n+1))
         enddo
      enddo
   endif
   alloc_phis = .false.

endif

call nc_check(nf90_close(nc_file_ID), 'read_cam_2Dreal', 'closing '//trim(file_name))

deallocate(var)

end subroutine read_cam_2Dreal

!-----------------------------------------------------------------------

subroutine read_cam_2Dint(file_name, cfield, field, num_dim1, num_dim2)

! Read 2d integer field from, e.g., HommeMapping.nc

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

character(len=*),     intent(in)  :: file_name
character(len=*),     intent(in)  :: cfield
integer, allocatable, intent(out) :: field(:,:)
integer,              intent(out) :: num_dim1     !The dimension(s) of cfield
integer,              intent(out) :: num_dim2

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer :: nc_file_ID, nc_var_ID                     !NetCDF variables
integer :: field_dim_IDs(2)                          !Array of dimension IDs for cfield
character(len=NF90_MAX_NAME) :: name_dim1,name_dim2  !Names of dimensions of cfield

field_dim_IDs = MISSING_I                  !Array of dimension IDs for cfield

if (file_exist(file_name)) then
   call nc_check(nf90_open(path=trim(file_name), mode=nf90_nowrite, ncid=nc_file_ID), &
              'read_cam_2Dint', 'opening '//trim(file_name))
   if (print_details .and. output_task0) then
      write(string1,*) 'file_name for ',cfield,' is ', trim(file_name)
      call error_handler(E_MSG, 'read_cam_2Dint', string1,source,revision,revdate)
   endif

   ! get field id
   call nc_check(nf90_inq_varid(nc_file_ID, trim(cfield), nc_var_ID), &
              'read_cam_2Dint', 'inq_varid: '//cfield)

   ! get dimension 'id's
   call nc_check(nf90_inquire_variable(nc_file_ID, nc_var_ID, dimids=field_dim_IDs), &
              'read_cam_2Dint', 'inquire_variable: '//cfield)

   ! get dimension sizes
   ! The first spatial dimension is always present.
   call nc_check(nf90_inquire_dimension(nc_file_ID, field_dim_IDs(1), name_dim1, num_dim1 ), &
                 'read_cam_2Dint', 'inquire_dimension: '//name_dim1)
   if (field_dim_IDs(2) /= MISSING_I)  then
      call nc_check(nf90_inquire_dimension(nc_file_ID, field_dim_IDs(2), name_dim2, num_dim2 ), &
                    'read_cam_2Dint', 'inquire_dimension: '//name_dim2)
   else
      num_dim2 = 1
      name_dim2 = 'no2nd_dim_'
   endif

   if (print_details .and. output_task0) then
      write(string1,*) cfield,' dimensions num_dim1, num_dim2 = ',num_dim1, num_dim2
      call error_handler(E_MSG, 'read_cam_2Dint', string1,source,revision,revdate)
   endif
else
   write(string1,'(3A)') 'Required file "',trim(file_name),'" is missing.'
   call error_handler(E_ERR, 'read_cam_2Dint', string1, source, revision, revdate)
endif

! Allocate array, based on size of this variable on the file.
allocate(field(num_dim1,num_dim2))

if (field_dim_IDs(2) /= MISSING_I)  then
   call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, field, start=(/ 1, 1 /), &
        count=(/num_dim1, num_dim2 /)), 'read_cam_2Dint', trim(cfield))
else
   call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, field),  &
                  'read_cam_2Dint', trim(cfield))
endif

call nc_check(nf90_close(nc_file_ID), 'read_cam_2Dint', 'closing '//trim(file_name))

end subroutine read_cam_2Dint

!-----------------------------------------------------------------------

subroutine create_cs_grid_arrays()

! Subroutine to create arrays of relationships between cubed sphere nodes (corners)
! and cell centers, including bearings between nodes.
! These will be used to identify the cell containing an observation.
! The relationships read from HommeMapping.nc will be augmented.
! All will be stored in global storage, and written to a new file for
! subsequent use.

! Local variables
integer  :: sh_corn(4), n(4)     ! Shifted corners to put closest at the origin.
integer  :: col, nbr, c, cent            ! Indices for loops.
integer  :: num_n, min_ind(1)
real(r8) :: dist, angle
real(r8) :: bearings(3), x_planar(3), y_planar(3)

! ncol = number of nodes/corners/grid points.  Global storage.
! corners = the names of the corners associated with each cell center
! neighbors = the nodes around each node which partner to make the sides of the cells
!             which may contain an observation.

! Get array of corner nodes/columns which define the cells (identified by 'center').
if (file_exist(homme_map_file)) then
   call read_cam_2Dint(homme_map_file, 'element_corners', corners,ncenters,ncorners)

   if (ncenters /= (ncol -2) ) then
      write(string1, *) trim(homme_map_file),' ncenters inconsistent with ncol-2 ', ncenters, ncol
      call error_handler(E_ERR,'create_cs_grid_arrays',string1,source,revision,revdate)
   endif

      allocate(num_nghbrs           (ncol), &
               centers(max_neighbors,ncol), &
               a          (3,ncorners,ncenters),   &
               b          (2,ncorners,ncenters),   &
               x_ax_bearings(ncorners,ncenters))

   num_nghbrs    = 0
   centers       = MISSING_I
   a             = MISSING_R8
   b             = MISSING_R8
   x_planar      = MISSING_R8
   y_planar      = MISSING_R8
   x_ax_bearings = MISSING_R8
else
   write(string1, *) 'CAM-SE grid file "',trim(homme_map_file),'" can not be found '
   call error_handler(E_ERR,'create_cs_grid_arrays',string1,source,revision,revdate)
endif

! Invert the element_corners array to compile all of the neighbors of each node (corner).
! Loop over HommeMapping cell centers.
Quads: do cent = 1,ncenters
   Corns: do c = 1,4
      ! Get the node numbers that define this cell
      ! and shift (rotate) them to create a separate mapping for each corner/node.
      ! Shift the section of corners 1 place to the 'left'/lower for the first corner,
      ! 2 for the 2nd, etc.  This will put the node closest to the ob in position 4
      ! (of the shifted corners). Then the (x,y) origin will the the closest node,
      ! and the indexing of the a,b,x_ax_bearing arrays will be easy.
      ! Shifting preserves the order of the corners as we go around the cell (clockwise).
      sh_corn = cshift(corners(cent,:), c)

      ! Check a few cells for corner consistency.
      if (print_details .and. sh_corn(4) < 10) then
         write(string1,'(A,5I7,/31X,4I7)') 'c, 4 corners, shifted = ', &
              c,(corners(cent,nbr),nbr=1,4),(sh_corn(nbr),nbr=1,4)
         call error_handler(E_MSG,'create_cs_grid_arrays',string1,source,revision,revdate)
      endif

      ! Increment the number of neighbors (and centers ) this corner(node) has.
      ! sh_corn(4) is used for all cases because the corner we're working on always
      ! ends up in that position, when c is incremented, then the corners are shifted.
      n(c) = num_nghbrs(sh_corn(4)) + 1

      ! Update the number of neighbors of each corner of the cell,
      num_nghbrs(sh_corn(4)) = n(c)

      ! Store the info that this center is associated with 4 node->neighbor pairs.
      centers(n(c),sh_corn(4)) = cent

      ! Define the planar coordinates for this center/cell and this corner/node.
      ! The 4th corner is the origin, and the cell side from the 4th to the 3rd is
      ! the x-axis of this cell's coordinate system for this corner.
      ! This is established in the definition of bearings().
      ! This choice makes mapping coefficients a(0) and b(0) = 0 (see below).
      ! It also helps make the indexing of bearings easy to use and store.

      ! Check a few cells for corner consistency
      if (print_details .and. sh_corn(4) < 10) then
         write(string1,'(A,3F10.6)') 'lon1, lat1 = ', lon_rad(sh_corn(4)), lat_rad(sh_corn(4))
         call error_handler(E_MSG, 'create_cs_grid_arrays', string1, source, revision, revdate)
      endif

      ! Descend through neighbors so that bearings(3) is already defined when needed at loop end.
      do nbr = 3,1,-1
         ! Bearings from the current origin node of the cell to the other 3.
         bearings(nbr) = bearing(lon_rad(sh_corn(4)),   lat_rad(sh_corn(4)),  &
                                 lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)) )

         dist = get_dist(cs_locs(sh_corn(4)), cs_locs(sh_corn(nbr)), 0, 0, .true.)

         if (sh_corn(4) < 10) then
            write(string1,'(A,3F10.6)') 'create_cs_grid:    lon2, lat2, bearing = ', &
                 lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)), bearings(nbr)
            call error_handler(E_MSG, 'create_cs_grid_arrays', string1, source, revision, revdate)
         endif

         ! This difference order looks wrong, but we need to change the sign of angles from the
         ! clockwise direction used by bearings to the counterclockwise direction used by
         ! trig functions.
         angle = bearings(3) - bearings(nbr)

         ! Normalize to -PI < angle <= PI.
         angle = mod(angle,PI) - PI*int(angle/PI)

         ! Set the planar location of this corner/node.
         x_planar(nbr) = dist * cos(angle)
         y_planar(nbr) = dist * sin(angle)

      enddo

      ! Store the baseline for use when interpolating to an ob location.
      x_ax_bearings(c,cent) = bearings(3)

      ! Define another bearings array to allow coord_ind_cs to find the right
      ! cell around the closest node by using a search through bearings,
      ! rather than a call to unit_square_location.
      !   Propagate sort_bearings to writing and reading of HommeMapping_cs_grid.nc file.
      !   real(r8), allocatable, :: sort_bearings(max_neighbors,ncol)
      !   sort_bearings(n(c),sh_corn(4)) = bearings(3)
      ! But see ordering of bearings, commented out below.

      ! Define quantities used to map the planar coordinate system of each cell
      ! to the unit square coordinate system.

      ! I'll use the mapping from planar space (x,y) to unit square space (l,m):
      ! x = a0 + a1*l*m + a2*m + a3*l
      ! y = b0 + b1*l*m + b2*m + b3*l
      ! The 4 corners (x,y) can be mapped to the four corners (l,m) to yield 4 equations.
      ! This can be written as vec_x = mat_A * vec_a^T.
      ! Then AI is the inverse of the mapping from physical space to the unit square space,
      ! and the coefficients of the mapping, aN and bN, can be calculated from:
      !    a = matmul(AI,x_planar(0:3))
      !    b = matmul(AI,y_planar(0:3))
      ! But the mapping from (lon,lat) to (x,y) space put corner "4" of the cell at (x4,y4) = (0,0)
      ! and corner "3" at (x3,y3) = (d3,0)
      ! This ends up making a0 = b0 = b3 = 0, and the equations simplify to the point
      ! that it doesn't make sense to encode this tranformation in a matrix.
      ! Replace matrix and matmul with simpler direct equations:
      a(3,c,cent) = x_planar(3)
      a(2,c,cent) = x_planar(1)
      a(1,c,cent) = x_planar(2) - x_planar(1) - x_planar(3)
      b(2,c,cent) = y_planar(1)
      b(1,c,cent) = y_planar(2) - y_planar(1)

      if (cent < 10) then
         write(string1,'(A,1p4E12.4)') 'create_cs_grid_arrays: a = ',(a(nbr,c,cent),nbr=1,3)
         write(string2,'(A,1p4E12.4)') 'create_cs_grid_arrays: b = ',(b(nbr,c,cent),nbr=1,2)
         call error_handler(E_MSG, 'create_cs_grid_arrays', string1, source, revision, revdate,text2=string2)
      endif

      if (a(3,c,cent)* a(2,c,cent) *a(1,c,cent) == 0.0_r8) then
         write(string1,'(A,2I8,A,1p3E12.4)') 'a(:,',c,cent,') = ',(a(nbr,c,cent),nbr=1,3)
         write(string2,'(A,(6X,2F10.6))') 'create_cs_grid:    lon, lat for nghr=1-3  = ', &
              (lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)), nbr=1,3)
         write(string3,'(A,5I7,/31X,4I7)') 'c, 4 corners, shifted = ', &
               c,(corners(cent,nbr),nbr=1,4),(sh_corn(nbr),nbr=1,4)
         call error_handler(E_MSG, 'create_cs_grid_arrays', string1, source, revision, revdate,text2=string2,text3=string3)
      endif

      if (b(2,c,cent) *b(1,c,cent) == 0.0_r8) then
         write(string1,'(A,2I8,A,1p3E12.4)') 'b(:,',c,cent,') = ',(b(nbr,c,cent),nbr=1,2)
         write(string2,'(A,(6X,2F10.6))') 'create_cs_grid:    lon, lat for nghr=1-3  = ', &
              (lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)), nbr=1,3)
         write(string3,'(A,5I7,/31X,4I7)') 'c, 4 corners, shifted = ', &
               c,(corners(cent,nbr),nbr=1,4),(sh_corn(nbr),nbr=1,4)
         call error_handler(E_MSG, 'create_cs_grid_arrays', string1, source, revision, revdate,text2=string2,text3=string3)
      endif

   enddo Corns
enddo Quads


! Check that all nodes have at least 3 neighbors and no more than 6.
do col = 1,ncol
   if (num_nghbrs(col) < 3 .or. num_nghbrs(col) > max_neighbors) then
      write(string1,'(A,I6,A,6I8)') 'num_nghbrs(',col,') <3 or >6: ', num_nghbrs(col)
      call error_handler(E_ERR,'create_cs_grid_arrays',string1,source,revision,revdate)
   endif
enddo

! There's code in earlier versions of model_mod to
! reorder the neighbors so that they are sequential around each node
! to make the search for the cell containing an ob faster.

return

end subroutine create_cs_grid_arrays

!-----------------------------------------------------------------------

subroutine nc_write_cs_grid_file(cs_grid_file, homme_map_file)

! Write out the number of neighbors, the neighbors, corners, centers, and bearings
! to a netCDF file once for this grid at the beginning of the assimilation.
! Called by static_init_model.

character(len=*), intent(in) :: cs_grid_file
character(len=*), intent(in) :: homme_map_file

integer :: nc_file_ID
integer ::                               &
        ncenters_ID,       centers_var_ID,  &
        ncorners_ID,       corners_var_ID,  &
               a_ID,             a_var_ID,  &
               b_ID,             b_var_ID,  &
            ncol_ID, x_ax_bearings_var_ID,  &
   max_neighbors_ID,    num_nghbrs_var_ID


! Create the file
call nc_check(nf90_create(path=trim(cs_grid_file), cmode=NF90_SHARE, ncid=nc_file_ID), &
              'nc_write_cs_grid_file', 'create '//trim(cs_grid_file))

write(string1,*) trim(cs_grid_file),' is nc_file_ID ',nc_file_ID
call error_handler(E_MSG,'nc_write_cs_grid_file',string1,source,revision,revdate)

! Define the dimensions
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                          &
              name="ncenters",      len = ncenters,      dimid = ncenters_ID), &
              'nc_write_cs_grid_file', 'def_dim ncenters '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                          &
              name="ncorners",      len = ncorners,      dimid = ncorners_ID), &
              'nc_write_cs_grid_file', 'def_dim ncorners '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                  &
              name="max_neighbors", len = max_neighbors, dimid = max_neighbors_ID), &
              'nc_write_cs_grid_file', 'def_dim max_neighbors'//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                  &
              name="ncol",          len = ncol,          dimid = ncol_ID), &
              'nc_write_cs_grid_file', 'def_dim ncol '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                  &
              name="ncoef_a",          len = 3,          dimid = a_ID), &
              'nc_write_cs_grid_file', 'def_dim a '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=nc_file_ID,                                  &
              name="ncoef_b",          len = 2,          dimid = b_ID), &
              'nc_write_cs_grid_file', 'def_dim b '//trim(cs_grid_file))

! Write Global Attributes
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "title", trim(cs_grid_file)), &
              'nc_write_cs_grid_file',   'put_att title '//trim(cs_grid_file))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "model_mod_source", source ), &
              'nc_write_cs_grid_file',   'put_att model_mod_source '//trim(cs_grid_file))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "model_mod_revision", revision ), &
              'nc_write_cs_grid_file',   'put_att model_mod_revision '//trim(cs_grid_file))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "model_mod_revdate", revdate ), &
              'nc_write_cs_grid_file',   'put_att model_mod_revdate '//trim(cs_grid_file))

call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "elements_per_cube_edge", ne ), &
              'nc_write_cs_grid_file',   'put_att elements_per_cube_edge '//trim(cs_grid_file))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "nodes_per_element_edge", np ), &
              'nc_write_cs_grid_file',   'put_att nodes_per_elements_edge '//trim(cs_grid_file))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "HommeMapping_file", homme_map_file ), &
              'nc_write_cs_grid_file',   'put_att HommeMapping_file '//trim(cs_grid_file))

! Create variables and attributes.
call nc_check(nf90_def_var(ncid=nc_file_ID, name="num_nghbrs", xtype=nf90_int, &
              dimids=(/ ncol_ID /), varid=num_nghbrs_var_ID),  &
              'nc_write_cs_grid_file', 'def_var num_nghbrs')
call nc_check(nf90_put_att(nc_file_ID, num_nghbrs_var_ID, "long_name", &
              "number of neighbors of each node/column"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, num_nghbrs_var_ID, "units",     "nondimensional"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(nc_file_ID, num_nghbrs_var_ID, "valid_range", &
              (/ 1, max_neighbors /)), 'nc_write_cs_grid_file', 'put_att valid_range')

call nc_check(nf90_def_var(ncid=nc_file_ID, name="centers", xtype=nf90_int, &
              dimids=(/ max_neighbors_ID, ncol_ID /), varid=centers_var_ID),  &
              'nc_write_cs_grid_file', 'def_var centers')
call nc_check(nf90_put_att(nc_file_ID, centers_var_ID, "long_name", &
              "cells which use node/column as a corner"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, centers_var_ID, "units",     "nondimensional"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(nc_file_ID, centers_var_ID, "valid_range", &
              (/ 1, ncenters /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(nc_file_ID, centers_var_ID, "missing_value", &
              (/ MISSING_I /)), 'nc_write_cs_grid_file', 'put_att missing_value')

call nc_check(nf90_def_var(ncid=nc_file_ID, name="corners", xtype=nf90_int, &
              dimids=(/ ncenters_ID, ncorners_ID /), varid=corners_var_ID),  &
              'nc_write_cs_grid_file', 'def_var corners')
call nc_check(nf90_put_att(nc_file_ID, corners_var_ID, "long_name", &
              "corners/nodes of each cell "), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, corners_var_ID, "units",     "nondimensional"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(nc_file_ID, corners_var_ID, "valid_range", &
              (/ 1, ncol /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(nc_file_ID, corners_var_ID, "missing_value", &
              (/ MISSING_I /)), 'nc_write_cs_grid_file', 'put_att missing_value')

call nc_check(nf90_def_var(ncid=nc_file_ID, name="a", xtype=nf90_double, &
              dimids=(/ a_ID, ncorners_ID, ncenters_ID /), varid=a_var_ID),  &
              'nc_write_cs_grid_file', 'def_var a')
call nc_check(nf90_put_att(nc_file_ID, a_var_ID, "long_name",  &
              "Coefficients of mapping from planar x coord to unit square"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, a_var_ID, "units",     "nondimensional"), &
              'nc_write_cs_grid_file', 'units')
!call nc_check(nf90_put_att(nc_file_ID, a_var_ID, "valid_range", &
!              (/ 1, ncol /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(nc_file_ID, a_var_ID, "missing_value", &
              (/ MISSING_R8 /)), 'nc_write_cs_grid_file', 'put_att missing_value')


call nc_check(nf90_def_var(ncid=nc_file_ID, name="b", xtype=nf90_double, &
              dimids=(/ b_ID, ncorners_ID, ncenters_ID /), varid=b_var_ID),  &
              'nc_write_cs_grid_file', 'def_var b')
call nc_check(nf90_put_att(nc_file_ID, b_var_ID, "long_name", &
              "Coefficients of mapping from planar y coord to unit square"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, b_var_ID, "units",     "nondimensional"), &
              'nc_write_cs_grid_file', 'units')
!call nc_check(nf90_put_att(nc_file_ID, b_var_ID, "valid_range", &
!              (/ 1, ncol /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(nc_file_ID, b_var_ID, "missing_value", &
              (/ MISSING_R8 /)), 'nc_write_cs_grid_file', 'put_att missing_value')

call nc_check(nf90_def_var(ncid=nc_file_ID, name="x_ax_bearings", xtype=nf90_double, &
              dimids=(/ ncorners_ID, ncenters_ID /), varid=x_ax_bearings_var_ID),  &
              'nc_write_cs_grid_file', 'def_var x_ax_bearings')
call nc_check(nf90_put_att(nc_file_ID, x_ax_bearings_var_ID, "long_name", &
              "bearing (clockwise from North) from origin node(corner 4) of each mapping to corner 3"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(nc_file_ID, x_ax_bearings_var_ID, "units",     "radians"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(nc_file_ID, x_ax_bearings_var_ID, "valid_range", &
              (/ -PI, PI /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(nc_file_ID, x_ax_bearings_var_ID, "missing_value", &
              (/ MISSING_R8 /)), 'nc_write_cs_grid_file', 'put_att missing_value')

! Leave define mode so we can fill
call nc_check(nf90_enddef(nc_file_ID), 'nc_write_cs_grid_file', 'enddef '//trim(cs_grid_file))

! sync to disk, but leave open
call nc_check(nf90_sync(nc_file_ID), 'nc_write_cs_grid_file', 'sync '//trim(cs_grid_file))

! Fill the variables
call nc_check(nf90_put_var(nc_file_ID, num_nghbrs_var_ID, num_nghbrs),  &
              'nc_write_cs_grid_file ','put_var num_nghbrs ')
call nc_check(nf90_put_var(nc_file_ID, centers_var_ID, centers),        &
              'nc_write_cs_grid_file ','put_var centers ')
call nc_check(nf90_put_var(nc_file_ID, corners_var_ID, corners),        &
              'nc_write_cs_grid_file ','put_var centers ')
call nc_check(nf90_put_var(nc_file_ID, a_var_ID, a),    &
              'nc_write_cs_grid_file ','put_var a ')
call nc_check(nf90_put_var(nc_file_ID, b_var_ID, b),    &
              'nc_write_cs_grid_file ','put_var b ')
call nc_check(nf90_put_var(nc_file_ID, x_ax_bearings_var_ID, x_ax_bearings),      &
              'nc_write_cs_grid_file ','put_var x_ax_bearings ')

call nc_check(nf90_close(nc_file_ID), 'nc_write_cs_grid_file', 'closing '//trim(cs_grid_file))

end subroutine nc_write_cs_grid_file

!-----------------------------------------------------------------------

subroutine nc_read_cs_grid_file()

! Read the number of neighbors, corners, centers, a and b coefficients, and x_ax_bearings
! from a netCDF file once for this grid at the beginning of the assimilation.

integer :: nc_file_ID, nc_var_ID, nc_size, n_dims, max_nghbrs, shp(2)
character(len=NF90_MAX_NAME) :: nc_name

! Open the cubed sphere grid relationships file
call nc_check(nf90_open(path=trim(cs_grid_file), mode=nf90_nowrite, ncid=nc_file_ID), &
      'nc_read_cs_grid_file', 'opening '//trim(cs_grid_file))

! learn how many dimensions are defined in this file.
call nc_check(nf90_inquire(nc_file_ID, n_dims), 'nc_read_cs_grid_file', 'inquire n_dims')

! Dimensions written out:
!              name="ncenters",      len = ncenters,      dimid = ncenters_ID), &
!              name="ncorners",      len = ncorners,      dimid = ncorners_ID), &
!              name="max_neighbors", len = max_neighbors, dimid = max_neighbors_ID), &
!              name="ncol",          len = ncol,          dimid = ncol_ID), &
!              name="ncoef_a",       len = 3,             dimid = a_ID), &
!              name="ncoef_b",       len = 2,             dimid = b_ID), &
call nc_check(nf90_inquire_dimension(nc_file_ID, 1, nc_name, ncenters), &
              'nc_read_cs_grid_file', 'inquire for '//trim(nc_name))

call nc_check(nf90_inquire_dimension(nc_file_ID, 2, nc_name, ncorners), &
              'nc_read_cs_grid_file', 'inquire for '//trim(nc_name))

call nc_check(nf90_inquire_dimension(nc_file_ID, 3, nc_name, max_nghbrs), &
              'nc_read_cs_grid_file', 'inquire for '//trim(nc_name))
if (trim(nc_name) /= 'max_neighbors') then
   write(string1, *) trim(cs_grid_file),' max_nghbrs does not match ', trim(model_config_file)
   call error_handler(E_ERR,'nc_read_cs_grid_file',string1,source,revision,revdate)
endif
! Check value against the namelist/parameter value.
if (max_nghbrs /= max_neighbors) then
   write(string1, *) trim(cs_grid_file),' max_nghbrs does not match max_neighbors', &
         max_nghbrs,max_neighbors
   call error_handler(E_ERR,'nc_read_cs_grid_file',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(nc_file_ID, 4, nc_name, nc_size), &
              'nc_read_cs_grid_file', 'inquire for '//trim(nc_name))
if (nc_size == ncol .and. trim(nc_name) == 'ncol') then
   allocate (corners(ncenters,ncorners),  &
             num_nghbrs           (ncol), &
             centers(max_neighbors,ncol), &
             x_ax_bearings  (ncorners,ncenters), &
             a            (3,ncorners,ncenters), &
             b            (2,ncorners,ncenters)  )
   ! Initialize the grid variables
   num_nghbrs    = MISSING_I
   centers       = MISSING_I
   a             = MISSING_R8
   b             = MISSING_R8
   x_ax_bearings = MISSING_R8

   if (allocated(centers) .and. output_task0 .and. print_details) then
      shp = shape(centers)
      write(string1,*) 'Shape of centers = ',shp
      call error_handler(E_MSG,'nc_read_cs_grid_file',string1,source,revision,revdate)
   endif
   if (allocated(corners) .and. output_task0 .and. print_details) then
      shp = shape(corners)
      write(string1,*) 'Shape of corners = ',shp
      call error_handler(E_MSG,'nc_read_cs_grid_file',string1,source,revision,revdate)
   endif
else
   write(string1,*) trim(cs_grid_file),' ncol does not match ', trim(model_config_file)
   call error_handler(E_ERR,'nc_read_cs_grid_file',string1,source,revision,revdate)
endif

call nc_check(nf90_inq_varid(nc_file_ID, 'num_nghbrs', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid num_nghbrs')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, num_nghbrs ), &
                'nc_read_cs_grid_file', 'get_var num_nghbrs')

call nc_check(nf90_inq_varid(nc_file_ID, 'centers', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid centers')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, centers ), &
                'nc_read_cs_grid_file', 'get_var centers')

call nc_check(nf90_inq_varid(nc_file_ID, 'corners', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid corners')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, corners, &
                           start=(/ 1, 1 /),count=(/ ncenters, ncorners /) ), &
                'nc_read_cs_grid_file', 'get_var corners')

call nc_check(nf90_inq_varid(nc_file_ID, 'a', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid a')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, a ), &
                'nc_read_cs_grid_file', 'get_var a')

call nc_check(nf90_inq_varid(nc_file_ID, 'b', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid b')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, b ), &
                'nc_read_cs_grid_file', 'get_var b')

call nc_check(nf90_inq_varid(nc_file_ID, 'x_ax_bearings', nc_var_ID), &
                'nc_read_cs_grid_file', 'inq_varid x_ax_bearings')
call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, x_ax_bearings ), &
                'nc_read_cs_grid_file', 'get_var x_ax_bearings')

call nc_check(nf90_close(nc_file_ID), 'nc_read_cs_grid_file', 'closing '//trim(cs_grid_file))

end subroutine nc_read_cs_grid_file

!-----------------------------------------------------------------------

real function bearing(lon1,lat1,lon2,lat2)

! Calculate the direction along the great circle from point 1 on a sphere
! to point 2, relative to north.
! All inputs should have units of radians.
! Output is radians.
! From http://www.movable-type.co.uk/scripts/latlong.html

real(r8), intent(in)    :: lon1,lat1, lon2,lat2

real(r8) :: lon1c,lon2c, cos_lat2, del_lon

real(r8), parameter :: half_PI = PI*0.5_r8

! Make sure the poles are handled consistently:
! If the pole point is the origin point, and the longitude of the pole point is
! defined as 0.0, then the bearing to a nearby point will = the longitude of the point.
! This is consistent/continuous with the bearing from points extremely near
! the pole.
if (half_PI - abs(lat1) < epsilon(lat1)) then
   lon1c = 0.0_r8
else
   lon1c = lon1
endif
if (half_PI - abs(lat2) < epsilon(lat2)) then
   lon2c = 0.0_r8
else
   lon2c = lon2
endif

cos_lat2 = cos(lat2)
del_lon  = lon2c - lon1c

! Normalize del_lon to -pi<=angle<=pi.
del_lon = mod(del_lon,PI) - PI*int(del_lon/PI)
bearing = atan2(cos_lat2*sin(del_lon),  &
                cos(lat1)*sin(lat2) - sin(lat1)*cos_lat2*cos(del_lon) )

end function bearing

!-----------------------------------------------------------------------

subroutine nc_read_model_atts(nc_file_ID, att, att_vals)

! reads the value of an attribute for each of the fields in cflds.
!
! should be called with att = one of the attributes from the program variable
! input file, which will be written to the Posterior and Prior.nc files

integer,                      intent(in)  :: nc_file_ID
character(len=*),             intent(in)  :: att
character(len=nf90_max_name), intent(out) :: att_vals(nflds) 

integer :: i, ierr
integer :: nc_var_ID, att_type

if (print_details .and. output_task0) then
   write(string1,*) 'nc_read_model_atts: reading ',trim(att)
      call error_handler(E_MSG, 'nc_read_model_atts', string1,source,revision,revdate)
endif

do i = 1,nflds
   att_vals(i) = ' '
   call nc_check(nf90_inq_varid(nc_file_ID, cflds(i), nc_var_ID), 'nc_read_model_atts', &
                 'inq_varid '//trim(cflds(i)))

   ierr = nf90_inquire_attribute(nc_file_ID, nc_var_ID, att)

   if (ierr == nf90_noerr) then
      call nc_check(nf90_get_att(nc_file_ID, nc_var_ID, att, att_vals(i)), &
                    'nc_read_model_atts', 'get_att '//trim(att))
      if (print_details .and. output_task0) then
         write(string1,'(A,1X,I6,1X,A,1X,A)') att, nc_var_ID, cflds(i), trim(att_vals(i))
         call error_handler(E_MSG, 'nc_read_model_atts', string1,source,revision,revdate)
      endif
   endif
enddo

end subroutine nc_read_model_atts

!-----------------------------------------------------------------------

subroutine nc_read_global_int_att(nc_file_ID, att, att_val)

! Reads the value of a global attribute.

integer,          intent(in)  :: nc_file_ID
character(len=*), intent(in)  :: att
integer,          intent(out) :: att_val

integer :: ierr

! NF90_GLOBAL is the psuedo-variable name used for global attributes.
ierr = nf90_inquire_attribute(nc_file_ID, NF90_GLOBAL, att)

if (ierr == nf90_noerr) then
   call nc_check(nf90_get_att(nc_file_ID, NF90_GLOBAL, att, att_val), &
                 'nc_read_global_int_att', 'get_att '//trim(att))
   if (print_details .and. output_task0) then
      write(string1,'(A,I5,2A, I6)') 'nc_read_global_int_att for file ',nc_file_ID, &
                                    ' attribute and value = ',trim(att), att_val
      call error_handler(E_MSG, 'nc_read_global_int_att', string1,source,revision,revdate)
   endif
endif

end subroutine nc_read_global_int_att

!-----------------------------------------------------------------------

subroutine read_cam_coord(nc_file_ID, cfield, var)

! read CAM 'initial' file coordinate, i.e. 'lat', 'lon', 'gw', 'hyai',...

integer,            intent(in)    :: nc_file_ID
character(len=*),   intent(in)    :: cfield
type(grid_1d_type), intent(inout) :: var

integer :: i, coord_size   ! grid/array indices
integer :: nc_var_ID         ! file and field IDs
integer :: fld_exist       ! grid field may not exist CAM initial file (e.g. slat)
integer :: ncerr           ! other nc errors; don't abort
integer :: coord_dimid(1)  ! Coordinates can have only 1 dimension,
                           ! but this must be a vector.

! Some attributes are _Fillvalue (real) which I'll ignore for now.
! The following are used to repack the attributes I want into a compact form
integer :: num_atts, keep_atts
integer :: att_type
character(len=nf90_max_name)               :: att_name
character(len=nf90_max_name), allocatable  :: att_names(:)
character(len=nf90_max_name), allocatable  :: att_vals(:)
real(r8)                                   :: resol, resol_1, resol_n

coord_dimid = MISSING_I      

fld_exist = nf90_inq_varid(nc_file_ID, cfield, nc_var_ID)
if (fld_exist /= nf90_noerr ) then
   var%label = ' '
   return
endif

ncerr = nf90_inquire_variable(nc_file_ID, nc_var_ID, dimids=coord_dimid, nAtts=num_atts)
if (ncerr /= nf90_noerr ) then
   write(string1,*) 'Variable ',cfield,' dimids = ',coord_dimid(1)
   write(string2,*) 'NetCDF error code = ',nf90_strerror(ncerr)
   call error_handler(E_MSG, 'read_cam_coord', string1,source,revision,revdate, text2=string2)
   var%label = ' '
   var%dim_id = 0
   return
endif

if (coord_dimid(1) == 0) then
   coord_size = 1                 ! to handle P0
else
   coord_size = dim_sizes(coord_dimid(1))
endif

allocate(att_names(num_atts), att_vals(num_atts))

keep_atts = 0
do i=1,num_atts
   call nc_check(nf90_inq_attname(nc_file_ID, nc_var_ID, i, att_name), &
                 'read_cam_coord', 'inq_attname '//trim(att_name))

! CAM FV initial files have coordinates with attributes that are numerical, not character
! (_FillValue).  These are not used because the coordinates are dimensioned exactly
! the right size.  I'll test for the type of att, and if it's not char, I'll ignore it.

! Otherwise I need a var%atts_type and separate var%atts_vals_YYY for each NetCDF
! external type (6 of them) I might run into.

   call nc_check(nf90_inquire_attribute(nc_file_ID, nc_var_ID, att_name, xtype=att_type), &
                 'read_cam_coord', 'inquire_attribute '//trim(att_name))

   if (att_type == nf90_char) then
      keep_atts = keep_atts + 1
      att_vals(keep_atts) = ' '
      att_names(keep_atts) = att_name
      call nc_check(nf90_get_att(nc_file_ID, nc_var_ID, att_name, att_vals(keep_atts)), &
                    'read_cam_coord', 'get_att '//trim(att_name) )

   else
      if (output_task0) then
         write(string1,*) '                ignoring attribute ',trim(att_name),    &
                    ' because it is not a character type'
         call error_handler(E_MSG, 'read_cam_coord', string1,source,revision,revdate)
      endif
   endif
enddo

call create_grid_1d_instance(coord_size, keep_atts, var)

! The rest of this routine populates 'var' with values.

var%label = cfield
var%dim_id = coord_dimid(1)

do i = 1,keep_atts
   var%atts_names(i) = att_names(i)
   var%atts_vals(i)  = att_vals(i)
enddo

call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, var%vals, start=(/ 1 /), count=(/ coord_size /)), &
              'read_cam_coord', 'get_var '//cfield)

! Determine whether coordinate is regularly spaced,
! If so, store the coordinate resolution in the grid_1d_type.
if (cfield(1:2) == 'hy' .or. cfield(1:2) == 'P0') then
   var%resolution = MISSING_R8
else
   resol_1 = var%vals(2) - var%vals(1)
   if (resol_1 /= 0.0_r8) then
      var%resolution = resol_1

      ! Test all other coordinate spacings.  If any of them differ from the first
      ! by more than epsilon (smallest meaningful number relative to the coordinate spacings)
      ! then spacing is irregular.
      resol = 1.0_r8/resol_1
      Res: do i = 3,coord_size
         resol_n = var%vals(i) - var%vals(i-1)
         if (((resol_n - resol_1) *resol) > epsilon(resol_n)) then
            var%resolution = MISSING_R8
            exit Res
         endif
      enddo Res
   else
      var%resolution = MISSING_R8
   endif
endif

if (print_details .and. output_task0) then
   write(string1,'(3A,I6,A,I8,A,1pE12.4)')  'reading ',cfield,' using id ',nc_var_ID,  &
          ' size ',coord_size,' resolution ', var%resolution
   write(string2,*) 'first, last val: ', var%vals(1),var%vals(coord_size)
   call error_handler(E_MSG, 'read_cam_init', string1,source,revision,revdate, text2=string2)
endif

deallocate(att_names, att_vals)

end subroutine read_cam_coord

!-----------------------------------------------------------------------

subroutine create_grid_1d_instance(length, num_atts, var)

! Initializes an instance of a cam grid variable

integer,            intent(in )   :: length
integer,            intent(in )   :: num_atts
type(grid_1d_type), intent(inout) :: var
! Does 'var' need to have the TARGET attribute here?
! Metcalf p 50 says 'yes'.
! But Intel says that allocating an object gives it the target attribute:
! "If an object does not have the TARGET attribute or has not been allocated 
! (using an ALLOCATE statement), no part of it can be accessed by a pointer."
! And this has worked without specifying the 'target' attribute.

! Initialize the storage space and return
allocate(var%vals      (length))
allocate(var%atts_names(num_atts))
allocate(var%atts_vals (num_atts))

var%length = length
var%num_atts = num_atts

end subroutine create_grid_1d_instance

!-----------------------------------------------------------------------

subroutine end_grid_1d_instance(var)

! Ends an instance of a cam grid_1d variable

type(grid_1d_type), intent(inout) :: var

if (var%label == ' ') return

if (.not. allocated(var%vals)) then
   write(string1,*) 'Calling end_grid_1d_instance on an uninitialized grid_1d_type'
   call error_handler(E_ERR,'end_grid_1d_instance',string1, source, revision, revdate)
endif

deallocate(var%vals, var%atts_names, var%atts_vals)

end subroutine end_grid_1d_instance

!-----------------------------------------------------------------------

subroutine order_state_fields()

! Fills cflds with state_names for use in I/O of caminput.nc.
! Also assigns TYPE_s for use various routines.

integer :: i, i1, nfld

allocate(TYPE_1D(state_num_1d),TYPE_2D(state_num_2d),TYPE_3D(state_num_3d))

nfld = 0

! 0D fields
do i=1,state_num_0d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_0d(i)
enddo

! 1D fields (1 spatial *coordinate* on the CAM initial file.
! The field may have 2 *physical* spatial dimensions.
do i=1,state_num_1d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_1d(i)
   TYPE_1D(i) = nfld
!  Cubed sphere
   if (state_names_1d(i) == 'PS')      TYPE_PS      = nfld
   if (state_names_1d(i) == 'EFGWORO') TYPE_EFGWORO = nfld
   if (state_names_1d(i) == 'FRACLDV') TYPE_FRACLDV = nfld
   if (state_names_1d(i) == 'TBOT')    TYPE_TBOT    = nfld
   if (state_names_1d(i) == 'TS')      TYPE_TS      = nfld
   if (state_names_1d(i) == 'TSOCN')   TYPE_TSOCN   = nfld
enddo

! 2D fields
do i=1,state_num_2d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_2d(i)
   TYPE_2D(i) = nfld
   if (state_names_2d(i) == 'PS')      TYPE_PS      = nfld
   if (state_names_2d(i) == 'EFGWORO') TYPE_EFGWORO = nfld
   if (state_names_2d(i) == 'FRACLDV') TYPE_FRACLDV = nfld
   if (state_names_2d(i) == 'TBOT')    TYPE_TBOT    = nfld
   if (state_names_2d(i) == 'TS')      TYPE_TS      = nfld
   if (state_names_2d(i) == 'TSOCN')   TYPE_TSOCN   = nfld
!  Cubed sphere
   if (state_names_2d(i) == 'T')       TYPE_T = nfld
   if (state_names_2d(i) == 'U')       TYPE_U = nfld
   if (state_names_2d(i) == 'V')       TYPE_V = nfld
   if (state_names_2d(i) == 'Q')       TYPE_Q = nfld
   if (state_names_2d(i) == 'CLDICE')  TYPE_CLDICE = nfld
   if (state_names_2d(i) == 'CLDLIQ')  TYPE_CLDLIQ = nfld
   if (state_names_2d(i) == 'LCWAT')   TYPE_LCWAT  = nfld
   if (state_names_2d(i) == 'QCWAT')   TYPE_QCWAT  = nfld
   ! chem
   if (state_names_2d(i) == 'CO')      TYPE_CO  = nfld
   if (state_names_2d(i) == 'CO2')     TYPE_CO2 = nfld
   if (state_names_2d(i) == 'NO')      TYPE_NO  = nfld
   if (state_names_2d(i) == 'NO2')     TYPE_NO2 = nfld
   if (state_names_2d(i) == 'CH4')     TYPE_CH4 = nfld
   if (state_names_2d(i) == 'NH3')     TYPE_NH3 = nfld
   if (state_names_2d(i) == 'O')       TYPE_O   = nfld
   if (state_names_2d(i) == 'O3')      TYPE_O3  = nfld
enddo

! 3D fields (including q)
do i=1,state_num_3d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_3d(i)
   TYPE_3D(i) = nfld
   if (state_names_3d(i) == 'T')      TYPE_T = nfld
   if (state_names_3d(i) == 'U')      TYPE_U = nfld
   if (state_names_3d(i) == 'V')      TYPE_V = nfld
   if (state_names_3d(i) == 'US')     TYPE_U = nfld
   if (state_names_3d(i) == 'VS')     TYPE_V = nfld
   if (state_names_3d(i) == 'Q')      TYPE_Q = nfld
   if (state_names_3d(i) == 'CLDICE') TYPE_CLDICE = nfld
   if (state_names_3d(i) == 'CLDLIQ') TYPE_CLDLIQ = nfld
   if (state_names_3d(i) == 'LCWAT')  TYPE_LCWAT  = nfld
   if (state_names_3d(i) == 'QCWAT')  TYPE_QCWAT  = nfld
   ! chem
   if (state_names_3d(i) == 'CO')     TYPE_CO  = nfld
   if (state_names_3d(i) == 'CO2')    TYPE_CO2 = nfld
   if (state_names_3d(i) == 'NO')     TYPE_NO  = nfld
   if (state_names_3d(i) == 'NO2')    TYPE_NO2 = nfld
   if (state_names_3d(i) == 'CH4')    TYPE_CH4 = nfld
   if (state_names_3d(i) == 'NH3')    TYPE_NH3 = nfld
   if (state_names_3d(i) == 'O')      TYPE_O   = nfld
   if (state_names_3d(i) == 'O3')     TYPE_O3  = nfld

enddo

if (nfld /= nflds) then
   write(string1, *) 'nfld = ',nfld,', nflds = ',nflds,' must be equal '
   call error_handler(E_ERR, 'order_state_fields', string1, source, revision, revdate)
endif

if (output_task0) then
   if (print_details) then
      write(string1,'(A)') 'State vector is composed of these fields: '
      call error_handler(E_MSG, 'order_state_fields', string1, source, revision, revdate)
   !   write(string1,'((8(A8,1X)))') (cflds(i),i=1,nflds)
      do i=1,state_num_0d
         write(string1,'(A,I4)') cflds(i), TYPE_1D(i)
         call error_handler(E_MSG, 'order_state_fields', string1, source, revision, revdate)
      enddo
      i1 = state_num_0d
      do i=1,state_num_1d
         write(string1,'(A,I4)') cflds(i1+i), TYPE_1D(i)
         call error_handler(E_MSG, 'order_state_fields', string1, source, revision, revdate)
      enddo
      i1 = i1 + state_num_1d
      do i=1,state_num_2d
         write(string1,'(A,I4)') cflds(i1+i), TYPE_2D(i)
         call error_handler(E_MSG, 'order_state_fields', string1, source, revision, revdate)
      enddo
      i1 = i1 + state_num_2d
      do i=1,state_num_3d
         write(string1,'(A,I4)') cflds(i1+i), TYPE_3D(i)
         call error_handler(E_MSG, 'order_state_fields', string1, source, revision, revdate)
      enddo
      write(string1,'(A)')        'TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q, TYPE_CLDLIQ, TYPE_CLDICE = '
      write(string2,'((8(I8,1X)))') TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q, TYPE_CLDLIQ, TYPE_CLDICE
      call error_handler(E_MSG, 'order_state_fields', string1, source, revision, revdate, &
                         text2=string2)
      write(string1,'(A)')        'TYPE_CO, TYPE_CO2, TYPE_NO, TYPE_NO2, TYPE_CH4, TYPE_NH3, TYPE_O, TYPE_O3 = '
      write(string2,'((8(I8,1X)))') TYPE_CO, TYPE_CO2, TYPE_NO, TYPE_NO2, TYPE_CH4, TYPE_NH3, TYPE_O, TYPE_O3
      call error_handler(E_MSG, 'order_state_fields', string1, source, revision, revdate, &
                         text2=string2)
   else
      call error_handler(E_MSG, 'order_state_fields', 'State vector is composed of these fields: ')
      do i = 1,nflds
         call error_handler(E_MSG, 'order_state_fields', trim(cflds(i)))
      enddo
   endif
endif

end subroutine order_state_fields

!-----------------------------------------------------------------------

subroutine map_kinds()

! ? Should this be a function instead; removes need to dimension obs_loc_in arbitrarily
!   and wastefully.  But then it's called millions of times, instead of accessing a small
!   array that's defined once.

! Makes an array of 'locations within the state vector' of the obs kinds
! that come from obs_kind_mod, which we anticipate CAM's model_mod will need.
! The obs kind that's needed will be the index into this array,
! the corresponding value will be the position of that field (not individual variable)
! within the state vector according to state_name_Xd.
! This subroutine will be called from static_init_model, so it will not have to be
! recomputed for every ob.
! Also maps the local model_mod TYPE_s onto the DART KIND_s by the same mechanism.

! other KIND_ possibilities are listed after the 'use obs_kind_mod' statement

integer :: i

! Physically 2D fields
dart_to_cam_types(KIND_SURFACE_PRESSURE) = TYPE_PS
if (TYPE_PS /= MISSING_I) cam_to_dart_kinds(TYPE_PS) = KIND_SURFACE_PRESSURE

dart_to_cam_types(KIND_GRAV_WAVE_DRAG_EFFIC) = TYPE_EFGWORO
if (TYPE_EFGWORO /= MISSING_I) &
   cam_to_dart_kinds(TYPE_EFGWORO) = KIND_GRAV_WAVE_DRAG_EFFIC

dart_to_cam_types(KIND_GRAV_WAVE_STRESS_FRACTION) = TYPE_FRACLDV
if (TYPE_FRACLDV /= MISSING_I) &
   cam_to_dart_kinds(TYPE_FRACLDV) = KIND_GRAV_WAVE_STRESS_FRACTION

! dart_to_cam_types(KIND_SURFACE_TEMPERATURE  ?  ) = TYPE_TS
! dart_to_cam_types(KIND_SEA_SURFACE_TEMPERATURE  ?  ) = TYPE_TSOCN

! Physically 3D fields
dart_to_cam_types(KIND_TEMPERATURE)        = TYPE_T
dart_to_cam_types(KIND_U_WIND_COMPONENT)   = TYPE_U
dart_to_cam_types(KIND_V_WIND_COMPONENT)   = TYPE_V
dart_to_cam_types(KIND_SPECIFIC_HUMIDITY)  = TYPE_Q
dart_to_cam_types(KIND_CLOUD_LIQUID_WATER) = TYPE_CLDLIQ
dart_to_cam_types(KIND_CLOUD_ICE)          = TYPE_CLDICE
! dart_to_cam_types(KIND_CLOUD_WATER  ?  ) = TYPE_LCWAT
dart_to_cam_types(KIND_CO)  = TYPE_CO
dart_to_cam_types(KIND_CO2) = TYPE_CO2
dart_to_cam_types(KIND_NO)  = TYPE_NO
dart_to_cam_types(KIND_NO2) = TYPE_NO2
dart_to_cam_types(KIND_CH4) = TYPE_CH4
dart_to_cam_types(KIND_NH3) = TYPE_NH3
! dart_to_cam_types(KIND_O)   = TYPE_O
dart_to_cam_types(KIND_O3)  = TYPE_O3

if (TYPE_T      /= MISSING_I) cam_to_dart_kinds(TYPE_T)      = KIND_TEMPERATURE
if (TYPE_U      /= MISSING_I) cam_to_dart_kinds(TYPE_U)      = KIND_U_WIND_COMPONENT
if (TYPE_V      /= MISSING_I) cam_to_dart_kinds(TYPE_V)      = KIND_V_WIND_COMPONENT
if (TYPE_Q      /= MISSING_I) cam_to_dart_kinds(TYPE_Q)      = KIND_SPECIFIC_HUMIDITY
if (TYPE_CLDLIQ /= MISSING_I) cam_to_dart_kinds(TYPE_CLDLIQ) = KIND_CLOUD_LIQUID_WATER
if (TYPE_CLDICE /= MISSING_I) cam_to_dart_kinds(TYPE_CLDICE) = KIND_CLOUD_ICE
! cam_to_dart_kinds(TYPE_LCWAT) = KIND_CLOUD_WATER  ?
if (TYPE_CO     /= MISSING_I) cam_to_dart_kinds(TYPE_CO)  = KIND_CO
if (TYPE_CO2    /= MISSING_I) cam_to_dart_kinds(TYPE_CO2) = KIND_CO2
if (TYPE_NO     /= MISSING_I) cam_to_dart_kinds(TYPE_NO)  = KIND_NO
if (TYPE_NO2    /= MISSING_I) cam_to_dart_kinds(TYPE_NO2) = KIND_NO2
if (TYPE_CH4    /= MISSING_I) cam_to_dart_kinds(TYPE_CH4) = KIND_CH4
if (TYPE_NH3    /= MISSING_I) cam_to_dart_kinds(TYPE_NH3) = KIND_NH3
! if (TYPE_O  /= MISSING_I) cam_to_dart_kinds(TYPE_O)   = KIND_O
if (TYPE_O3     /= MISSING_I) cam_to_dart_kinds(TYPE_O3)  = KIND_O3


if (print_details .and. output_task0) then
   write(string1,*) 'OBS_KIND   FIELD_TYPE'
   call error_handler(E_MSG, 'map_kinds', string1,source,revision,revdate)
   do i=1,300
      if (dart_to_cam_types(i) /= MISSING_I) then
         write(string1,'(2I8)') i, dart_to_cam_types(i)
         call error_handler(E_MSG, 'map_kinds', string1,source,revision,revdate)
      endif
   enddo
endif

end subroutine map_kinds

!-----------------------------------------------------------------------

subroutine fill_gc()

! Subroutine to generate location_types of the cubed sphere grid
! and put them into get_close_type cs_gc, with other derived components.

integer :: c

allocate(cs_locs(ncol), cs_kinds(ncol))

! CS inputs in degrees.
do c=1,ncol
   cs_locs(c)  = set_location(lon%vals(c), lat%vals(c), MISSING_R8, VERTISUNDEF)
   cs_kinds(c) = 0
enddo

! Initialize cs_gc%maxdist using the maximum grid spacing.
! There will always be at least 2 nodes within 1 coarse_grid in all directions.
call get_close_maxdist_init(cs_gc, coarse_grid)

! Use cs_gc%maxdist and node locations to define the rest of cs_gc.
call get_close_obs_init(cs_gc, ncol, cs_locs)

end subroutine fill_gc

! End of static_init_model section
!#######################################################################

! Module I/O to/from DART and files

!-----------------------------------------------------------------------

subroutine read_cam_init(file_name, var, model_time)

! Fill the model_type 'var' using fields from a CAM initial file.
! Init_model_instance must be called before this subroutine.

! CAM initial files are used instead of restart files for (at least) 6 reasons.
! 1) The contents of the restart files vary depending on both the model release version
!    and the physics packages selected.
! 2) There is no metadata on the restart files describing the variables.
!    Some information can be tracked down in the atm.log file, but not all of it.
! 3) The restart files (for non-chemistry model versions) are much larger than the
!    initial files (and we need to deal with an ensemble of them).
! 4) The temperature on the restart files is virtual equivalent potential temperature (?),
!    which requires (at least) surface pressure, specific humidity, and sensible temperature
!    to calculate.
! 5) CAM does not call the initialization routines when restart files are used,
!    so fields which are not modified by DART may be inconsistent with fields which are.
! 6) If DART modifies the contents of the .r. restart file, it might also need to 
!    modify the contents of the .rs. restart file, which has similar characteristics 
!    (1-3 above) to the .r. file.

character(len=*), intent(in)    :: file_name
type(model_type), intent(inout) :: var
type(time_type),  intent(inout) :: model_time

integer :: i, k, n, m, ifld  
integer :: nc_file_ID, nc_var_ID, dimid, varid, dimlen
integer :: iyear, imonth, iday, ihour, imin, isec, rem
integer :: timestep
integer,  allocatable :: datetmp(:), datesec(:)
real(r8), allocatable :: temp_3d(:,:,:), temp_2d(:,:)

! read CAM 'initial' file domain info
call nc_check(nf90_open(path=file_name, mode=nf90_nowrite, ncid=nc_file_ID), &
      'read_cam_init', 'opening '//trim(file_name))

! Read the time of the current state.
! CAM initial files have two variables of length 'time' (the unlimited dimension): date, datesec
! The rest of the routine presumes there is but one time in the file -

call nc_check(nf90_inq_dimid(nc_file_ID, 'time', dimid), &
        'read_cam_init', 'inq_dimid time '//trim(file_name))
call nc_check(nf90_inquire_dimension(nc_file_ID, dimid, len=dimlen), &
        'read_cam_init', 'inquire_dimension time '//trim(file_name))

if (dimlen /= 1) then
   write(string1,*)trim(file_name),' has',dimlen,'times. Require exactly 1.'
   call error_handler(E_ERR, 'read_cam_init', string1, source, revision, revdate)
endif

allocate(datetmp(dimlen), datesec(dimlen))

call nc_check(nf90_inq_varid(nc_file_ID, 'date', varid), &
       'read_cam_init', 'inq_varid date '//trim(file_name))
call nc_check(nf90_get_var(nc_file_ID, varid, values=datetmp), &
       'read_cam_init', 'get_var date '//trim(file_name))

call nc_check(nf90_inq_varid(nc_file_ID, 'datesec', varid), &
       'read_cam_init', 'inq_varid datesec '//trim(file_name))
call nc_check(nf90_get_var(nc_file_ID, varid, values=datesec), &
       'read_cam_init', 'get_var datesec '//trim(file_name))

! for future extensibility, presume we find a 'timeindex' that we want.
! Since we only support 1 timestep in the file, this is easy.

timestep = 1

! The 'date' is YYYYMMDD ... datesec is 'current seconds of current day'
iyear  = datetmp(timestep) / 10000
rem    = datetmp(timestep) - iyear*10000
imonth = rem / 100
iday   = rem - imonth*100

ihour  = datesec(timestep) / 3600
rem    = datesec(timestep) - ihour*3600
imin   = rem / 60
isec   = rem - imin*60

deallocate(datetmp, datesec)

! some cam files are from before the start of the gregorian calendar.
! since these are 'arbitrary' years, just change the offset.

if (iyear < 1601) then
   write(string1,*)' '
   write(string2,*)'WARNING - ',trim(file_name),' changing year from ',iyear,'to',iyear+1601
   call error_handler(E_MSG, 'read_cam_init', string1, source, revision, &
                revdate, text2=string2,text3='to make it a valid Gregorian date.')
   write(string1,*)' '
   call error_handler(E_MSG, 'read_cam_init', string1, source, revision)
   iyear = iyear + 1601
endif

model_time = set_date(iyear,imonth,iday,ihour,imin,isec)

if (output_task0) then
   call print_date(model_time,' read_cam_init ... input date')
   call print_time(model_time,' read_cam_init ... input time')
   call print_date(model_time,' read_cam_init ... input date',logfileunit)
   call print_time(model_time,' read_cam_init ... input time',logfileunit)
endif

! The temporary arrays into which fields are read are dimensioned by the largest values of
! the sizes of the dimensions listed in f_dim_RANKd
! f_dim_max contents assume that time is always the last dimension on NetCDF files,
! so f_dim_max(4,3) and f_dim_max(3,2) are the non-spatial dimensions to ignore here.
allocate(temp_3d(f_dim_max(1,3),f_dim_max(2,3),f_dim_max(3,3)), &
         temp_2d(f_dim_max(1,2),f_dim_max(2,2)) )

ifld = 0
!0d fields; scalars are recognized and handled differently than vectors by NetCDF
do i= 1, state_num_0d
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(nc_file_ID, cflds(ifld), nc_var_ID), &
                'read_cam_init', 'inq_varid '//trim(cflds(ifld)))
   if (print_details .and. output_task0) then
      write(string1,*) 'reading ',cflds(ifld),' using id ',nc_var_ID
      call error_handler(E_ERR, 'read_cam_init', string1,source,revision,revdate)
   endif

   ! Fields on file are 1D; TIME(=1)
   call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, var%vars_0d(i) ), &
                'read_cam_init', 'get_var '//trim(cflds(ifld)))
enddo


!1d fields
do i= 1, state_num_1d
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(nc_file_ID, cflds(ifld), nc_var_ID), &
                'read_cam_init', 'inq_varid '//trim(cflds(ifld)))
   if (print_details .and. output_task0) then
      write(string1,*) 'reading ',cflds(ifld),' using id ',nc_var_ID
      call error_handler(E_MSG, 'read_cam_init', string1,source,revision,revdate)
   endif

   ! s_dim_1d should = f_dim_1d
   call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, var%vars_1d(1:s_dim_1d(i), i), &
                 start=(/ 1, timestep /), count=(/ f_dim_1d(1,i), 1/) ), &
                 'read_cam_init', 'get_var '//trim(cflds(ifld)))
enddo

!2d fields on file are 3D; 2 spatial dimensions, then  TIME(=1).
do i= 1, state_num_2d
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(nc_file_ID, cflds(ifld), nc_var_ID), &
                'read_cam_init', 'inq_varid '//trim(cflds(ifld)))
   if (print_details .and. output_task0) then
      write(string1,*) 'reading ',cflds(ifld),' using id ',nc_var_ID
      call error_handler(E_MSG, 'read_cam_init', string1,source,revision,revdate)
   endif

   ! Need to use temp_Nd; I am coding for not knowing what the 2 spatial dimensions of this field.
   call nc_check(nf90_get_var(nc_file_ID, nc_var_ID,  temp_2d(1:f_dim_2d(1,i), 1:f_dim_2d(2,i)),  &
                        start=(/ 1, 1, timestep/) ,count=(/ f_dim_2d(1,i),   f_dim_2d(2,i), 1/) ), &
                        'read_cam_init', 'get_var '//trim(cflds(ifld)))

   if (s_dim_2d(1,i) == f_dim_2d(1,i)) then
      var%vars_2d(1:s_dim_2d(1,i), 1:s_dim_2d(2,i),i) = &
          temp_2d(1:f_dim_2d(1,i), 1:f_dim_2d(2,i)  )

   elseif (s_dim_2d(1,i) == f_dim_2d(2,i)) then
      do k=1,s_dim_2d(1,i)
      do m=1,s_dim_2d(2,i)   ! first temp dim is inner loop for faster reads
         var%vars_2d(k,m,i) = temp_2d(m,k)
      enddo
      enddo
   else
      write(string1, *)'Dimension size for ',cflds(i),' in model_mod, ',s_dim_2d(1,i), &
         ', does not match sizes on file ',f_dim_2d(1,i),f_dim_2d(1,i)
      call error_handler(E_ERR, 'read_cam_init', string1, source, revision, revdate)
   endif
enddo

! Spatially 3d fields on file are 4D; lon, lev, lat, TIME(=1)
!                                 or; lon, lat, lev, TIME
do i=1, state_num_3d
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(nc_file_ID, cflds(ifld), nc_var_ID), &
                'read_cam_init', 'inq_varid '//trim(cflds(ifld)))
   if (print_details .and. output_task0) then
      write(string1,*) 'reading ',cflds(ifld),' using id ',nc_var_ID
      call error_handler(E_MSG, 'read_cam_init', string1,source,revision,revdate)
   endif

   call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, &
        temp_3d(1:f_dim_3d(1,i), 1:f_dim_3d(2,i), 1:f_dim_3d(3,i)), start=(/ 1, 1, 1, timestep/),  &
        count=(/  f_dim_3d(1,i),   f_dim_3d(2,i),   f_dim_3d(3,i),  1 /) ), &
        'read_cam_init', 'get_var '//trim(cflds(ifld)))

! Repackage depending on coord_order; put it into lev,lon,lat order.
   if (coord_order == 1) then
!     lon,lev,lat as in CAM <= 3.0.3
      do n=1,s_dim_3d(3,i)  ! lats
      do k=1,s_dim_3d(1,i)  ! levs
      do m=1,s_dim_3d(2,i)  ! lons/first file dim   put in inner loop for faster reads
         var%vars_3d(k,m,n, i) = temp_3d(m,k,n)
      enddo
      enddo
      enddo
   elseif (coord_order == 2) then
!     lon,lat,lev as in  CAM > 3.0.3
      do n=1,s_dim_3d(3,i)  ! lats
      do k=1,s_dim_3d(1,i)  ! levs
      do m=1,s_dim_3d(2,i)  ! lons/first file dim   put in inner loop for faster reads
         var%vars_3d(k,m,n, i) = temp_3d(m,n,k)
      enddo
      enddo
      enddo
   endif
enddo

call nc_check(nf90_close(nc_file_ID), 'read_cam_init', 'closing '//trim(file_name))

deallocate(temp_3d,temp_2d)

end subroutine read_cam_init

!-----------------------------------------------------------------------

subroutine write_cam_coord_def(nc_file_ID, c_name, coord, dim_id, c_id)

integer,            intent(in)  :: nc_file_ID
character(len=*),   intent(in)  :: c_name
type(grid_1d_type), intent(in)  :: coord
integer,            intent(in)  :: dim_id
integer,            intent(out) :: c_id

integer  :: i

call nc_check(nf90_def_var(nc_file_ID, name=c_name, xtype=nf90_double, dimids=dim_id, &
                        varid=c_id), 'write_cam_coord_def', 'def_var '//trim(c_name))

do i=1,coord%num_atts
   call nc_check(nf90_put_att(nc_file_ID, c_id, coord%atts_names(i), coord%atts_vals(i)), &
                 'write_cam_coord_def', 'put_att '//trim(coord%atts_names(i)))
enddo

end subroutine write_cam_coord_def

!-----------------------------------------------------------------------

subroutine write_cam_init(file_name, model_time, var)

! Write CAM 'initial' file fields (from var) that have been updated
! to a CAM initial file.

character(len=*), intent(in)    :: file_name
type(time_type),  intent(in)    :: model_time
type(model_type), intent(inout) :: var

type(time_type)       :: CAM_time
integer               :: i, k, n, m, ifld
integer               :: nc_file_ID, nc_var_ID, f_dim1, f_dim2
integer               :: dimid, dimlen, varid
integer               :: iyear, imonth, iday, ihour, imin, isec, leftover
integer               :: itime, timeindex

integer,  allocatable :: datetmp(:), datesec(:)
real(r8), allocatable :: temp_3d(:,:,:), temp_2d(:,:)

if (.not. module_initialized) call static_init_model()

call nc_check(nf90_open(path=trim(file_name), mode=nf90_write, ncid=nc_file_ID), &
           'write_cam_init', 'opening '//trim(file_name))

! Need to figure out which timeslot to update in the CAM initial file.
! It is not likely, but possible, that the initial file will have multiple
! timesteps in it. We have to figure out which slot matches the DART model time.
! the 'date' and 'datesec' variables contain the CAM state time.

call nc_check(nf90_inq_dimid(nc_file_ID, 'time', dimid), &
          'write_cam_init', 'inq_dimid time '//trim(file_name))
call nc_check(nf90_inquire_dimension(nc_file_ID, dimid, len=dimlen), &
          'write_cam_init', 'inquire_dimension time '//trim(file_name))

if (dimlen /= 1) then
   write(string1,*)'UNUSUAL - ',trim(file_name),' has',dimlen,'times. Expected 1.'
   call error_handler(E_MSG, 'write_cam_init', string1, source, revision, revdate, &
            text2='Searching for a matching time ...')
endif

allocate(datetmp(dimlen), datesec(dimlen))

call nc_check(nf90_inq_varid(nc_file_ID, 'date', varid), &
          'write_cam_init', 'inq_varid date '//trim(file_name))
call nc_check(nf90_get_var(nc_file_ID, varid, values=datetmp), &
          'write_cam_init', 'get_var date '//trim(file_name))

call nc_check(nf90_inq_varid(nc_file_ID, 'datesec', varid), &
          'write_cam_init', 'inq_varid datesec '//trim(file_name))
call nc_check(nf90_get_var(nc_file_ID, varid, values=datesec), &
          'write_cam_init', 'get_var datesec '//trim(file_name))

timeindex = -1
TIMELOOP: do itime = 1,dimlen

   iyear    = datetmp(itime)/10000
   leftover = datetmp(itime) - iyear*10000
   imonth   = leftover/100
   iday     = leftover - imonth*100
   ihour    = datesec(itime)/3600
   leftover = datesec(itime) - ihour*3600
   imin     = leftover/60
   isec     = leftover - imin*60

   CAM_time = set_date(iyear, imonth, iday, ihour, imin, isec)

   if (CAM_time == model_time) then
      if (dimlen /= 1) then
         write(string1,*)'Found matching time at index ',itime
         call error_handler(E_MSG, 'write_cam_init', string1, source, revision, revdate)
      endif

      timeindex = itime
      exit TIMELOOP
   endif

enddo TIMELOOP

deallocate(datetmp, datesec)

if (timeindex < 1) then

   call get_date(model_time, iyear, imonth, iday, ihour, imin, isec)

   write(string1,*)trim(file_name),' had no times that matched the model time.'
   write(string2,*)'model_time (YYYY MM DD) is ',iyear, imonth, iday
   write(string3,*)'model_time      (SSSSS) is ',isec + imin*60 + ihour*3600
   call error_handler(E_ERR, 'write_cam_init', string1, source, revision, revdate, &
            text2=string2,text3=string3)
endif

! So now we know that the right timeslot is 'timeindex'.

! The temporary arrays into which fields are read are dimensioned by the largest values of
! the sizes of the dimensions listed in coord_RANKd
allocate(temp_3d(f_dim_max(1,3),f_dim_max(2,3),f_dim_max(3,3)))
allocate(temp_2d(f_dim_max(1,2),f_dim_max(2,2)))

if (print_details .and. output_task0) then
   write(string1,*) 'write_cam_init; f_dim_max(:2) = ',f_dim_max(1,2),f_dim_max(2,2)
   call error_handler(E_MSG, 'write_cam_init', string1,source,revision,revdate)
endif

ifld = 0

! 0d fields are first ... there is no concern about shape or dimensions
do i = 1, state_num_0d
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(nc_file_ID, cflds(ifld), nc_var_ID), &
                 'write_cam_init', 'inq_var '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(nc_file_ID, nc_var_ID, var%vars_0d(i) ), &
                 'write_cam_init', 'put_var '//trim(cflds(ifld)))
enddo

! 1d fields
do i = 1, state_num_1d
   ! CS added this from 2d loop below.
   ! special code:  check and error out if the PS field has gone negative
   if (state_names_1d(i) == 'PS') then
      if (minval(var%vars_1d(:,i)) < 0.0_r8) then
         write(string1, *)'PS has negative values; should not happen'
         call error_handler(E_ERR, 'write_cam_init', string1, source, revision, revdate)
      endif
   endif
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(nc_file_ID, cflds(ifld), nc_var_ID), &
                 'write_cam_init', 'inq_var '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(nc_file_ID, nc_var_ID, var%vars_1d(1:s_dim_1d(i),i),   &
                             start=(/ 1, timeindex /), count = (/ s_dim_1d(i), 1 /)), &
                 'write_cam_init', 'put_var '//trim(cflds(ifld)))
enddo

do i = 1, state_num_2d
   ! special code:  check and error out if the PS field has gone negative
   if (state_names_2d(i) == 'PS') then
      if (minval(var%vars_2d(:,:,i)) < 0.0_r8) then
         write(string1, *)'PS has negative values; should not happen'
         call error_handler(E_ERR, 'write_cam_init', string1, source, revision, revdate)
      endif
   ! CS: this section copied from 3d below, and converted to 2d.
   elseif (state_names_2d(i) == 'Q') then
      where (var%vars_2d(:,:,i) < 1.e-12_r8) var%vars_2d(:,:,i) = 1.e-12_r8
   elseif (state_names_2d(i) == 'CLDLIQ' .or. &
            state_names_2d(i) == 'CLDICE') then
      where (var%vars_2d(:,:,i) < 0.0_r8)     var%vars_2d(:,:,i) = 0.0_r8
   elseif (state_names_2d(i) == 'T') then
      if (minval(var%vars_2d(:,:,i)) < 0.0_r8) then
         write(string1, *)'T has negative values; should not happen'
         call error_handler(E_ERR, 'write_cam_init', string1, source, revision, revdate)
      endif
   endif

   ! 2d fields ; tricky because coordinates may have been rearranged.

   if (f_dimid_2d(1,i) == s_dimid_2d(1,i)) then

      ! model_mod and caminput store this variable with the same coordinate order.
      f_dim1 = s_dim_2d(1,i)
      f_dim2 = s_dim_2d(2,i)
      do n=1,f_dim2
      do m=1,f_dim1
         temp_2d(m,n) = var%vars_2d(m,n,i)
      enddo
      enddo

   elseif (f_dimid_2d(1,i) == s_dimid_2d(2,i)) then

      ! model_mod and caminput store this variable with transposed coordinate order.
      f_dim1 = s_dim_2d(2,i)
      f_dim2 = s_dim_2d(1,i)
      do m=1,f_dim1
      do n=1,f_dim2
         temp_2d(m,n) = var%vars_2d(n,m,i)
      enddo
      enddo

   else

      ! There aren't any more to try
      write(string1, *) 'The dimension ID for ',state_names_2d(i), &
                        ' did not match an state dimension IDs'
      call error_handler(E_ERR, 'write_cam_init', string1, source, revision, revdate)

   endif

   ifld = ifld + 1
   call nc_check(nf90_inq_varid(nc_file_ID, trim(cflds(ifld)), nc_var_ID),           &
                 'write_cam_init','inq_varid '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(nc_file_ID, nc_var_ID,    temp_2d(1:f_dim1,1:f_dim2),     &
                      start=(/ 1, 1, timeindex /), count = (/ f_dim1,  f_dim2, 1/)),    &
                 'write_cam_init','put_var '//trim(cflds(ifld)))
enddo

! 3d fields; all 3 coordinates are present, and the order for model_mod fields is always the same.
do i = 1, state_num_3d
   ! special code:  set a minimum threshold for certain variables
   if (state_names_3d(i) == 'Q') then
      where (var%vars_3d(:,:,:,i) < 1.0e-12_r8) var%vars_3d(:,:,:,i) = 1.0e-12_r8
   elseif (state_names_3d(i) == 'CLDLIQ' .or. &
            state_names_3d(i) == 'CLDICE') then
      where (var%vars_3d(:,:,:,i) < 0.0_r8)     var%vars_3d(:,:,:,i) = 0.0_r8
   elseif (state_names_3d(i) == 'T') then
      if (minval(var%vars_3d(:,:,:,i)) < 0.0_r8) then
         write(string1, *)'T has negative values; should not happen'
         call error_handler(E_ERR, 'write_cam_init', string1, source, revision, revdate)
      endif
   endif

! CS; simplify by removing old CAM coordinate order?
   !  Repackage depending on coord_order then write the dummy variable.
   if (coord_order == 1) then
      ! lon,lev,lat as in original CAM
      do n=1,s_dim_3d(3,i)      ! lats
      do m=1,s_dim_3d(2,i)      ! lons
      do k=1,s_dim_3d(1,i)      ! levs
         temp_3d(m,k,n) = var%vars_3d(k,m,n,i)
      enddo
      enddo
      enddo
   elseif (coord_order == 2) then
      ! lon,lat,lev as in new CAM
      do n=1,s_dim_3d(3,i)
      do m=1,s_dim_3d(2,i)
      do k=1,s_dim_3d(1,i)
         temp_3d(m,n,k) = var%vars_3d(k,m,n,i)
      enddo
      enddo
      enddo
   endif

   ifld = ifld + 1
   call nc_check(nf90_inq_varid(nc_file_ID, trim(cflds(ifld)), nc_var_ID),   &
                 'write_cam_init', 'inq_varid '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(nc_file_ID, nc_var_ID,                                           &
        temp_3d(1:f_dim_3d(1,i), 1:f_dim_3d(2,i), 1:f_dim_3d(3,i)),                   &
        start=(/ 1, 1, 1, timeindex /), &
         count=(/ f_dim_3d(1,i),   f_dim_3d(2,i),   f_dim_3d(3,i), 1/)), &
        'write_cam_init', 'put_var '//trim(cflds(ifld)))
enddo

call nc_check(nf90_close(nc_file_ID), 'write_cam_init', 'close cam initial file')

deallocate(temp_3d, temp_2d)

end subroutine write_cam_init

!-----------------------------------------------------------------------

subroutine write_cam_times(model_time, adv_time)
! Not needed in CESM+DART framework

! Writes model time and advance time into a file called 'times',
! which is simply numbers.  A script reads those and passes them to CAM's build-namelist.
!  -namelist "&camexp START_YMD=$times[3] START_TOD=$times[4]
!                     STOP_YMD=$times[1] STOP_TOD=$times[2] NHTFRQ=$times[5] "
! End time is first, then beginning time

type(time_type), intent(in) :: model_time
type(time_type), intent(in) :: adv_time

integer :: tfile_unit, cam_date, cam_tod, nhtfrq
integer :: year, month, day, hour, minute, second
type(time_type) :: forecast_length

if (.not. module_initialized) call static_init_model()

! calculate number of hours in forecast, and pass to history tape
! write frequency

forecast_length = adv_time - model_time

call get_time(forecast_length, second, day)

hour = second/3600
minute = mod(second,3600)
if (minute/=0) &
   call error_handler(E_ERR, 'write_cam_times', &
      ' not integer number of hours; nhtfrq error', source, revision, revdate);

! convert to hours, and negative to signal units are hours

nhtfrq = -1*(day*24 + hour)


tfile_unit = open_file("times", "formatted", "write")

call get_date(adv_time, year, month, day, hour, minute, second)

cam_date = year*10000 + month*100 + day
cam_tod  = hour*3600  + minute*60 + second

write(tfile_unit,'(I8.8,1X,I8)') cam_date, cam_tod


call get_date(model_time, year, month, day, hour, minute, second)

cam_date = year*10000 + month*100 + day
cam_tod  = hour*3600  + minute*60 + second

write(tfile_unit,'(I8.8,1X,I8)') cam_date, cam_tod

write(tfile_unit,'(I8)') nhtfrq

close(tfile_unit)


end subroutine write_cam_times

!-----------------------------------------------------------------------
!>
!> Subroutine get_state_meta_data 
!> Given an integer index into the state vector structure, 
!> returns the associated location and vertical location type 'which_vert'.
!> Optionally returns the DART KIND of the variable.
!> 
!> @param[in]    index_in
!> The 'index' of a variable in the state vector, whose physical location 
!> and possibly variable kind are needed,
!>
!> @param[inout] location
!> The DART location_type location of the variable denoted by 'index'
!> 
!> @param[out]   var_kind
!> The optional argument which can return the DART KIND of the variable.


subroutine get_state_meta_data(index_in, location, var_kind)

! Given an integer index into the state vector structure, returns the
! associated location.
! The location may have components that are MISSING_R8 values, since some fields
! don't have locations in all three dimensions, i.e. PS has no vertical level,
! and other fiendish fields to be devised by parameterization studies may not
! have a longitude, or latitude.  The which_vert should take care of the vertical
! coordinate (it will be ignored), but the others will require more interesting  fixes.
! See order_state_fields for the KIND_s (and corresponding model_mod TYPE_s).

integer,             intent(in)    :: index_in
type(location_type), intent(inout) :: location
integer, optional,   intent(out)   :: var_kind

integer  :: which_vert
integer  :: i, indx, index_1, index_2, index_3, nfld
integer  :: box, slice
logical  :: lfound

real(r8) :: lon_val, lat_val, lev_val

character(len=8)   :: dim_name

if (.not. module_initialized) call static_init_model()

lfound    = .false.

! In order to find what variable this is, and its location, I must subtract the individual
! components of the state vector, since they may have varying sizes.
! Save the original index.
! index_in will be < 0 if it's an identity obs (called from convert_vert)

indx = abs(index_in)
which_vert = MISSING_I
index_1 = 0
index_2 = 0
index_3 = 0
nfld = 0
lon_val = MISSING_R8
lat_val = MISSING_R8
lev_val = MISSING_R8

! Cycle through 0d state variables
State_0D: do i=1,state_num_0d
   nfld = nfld + 1
   if (indx == i ) then
      which_vert = VERTISUNDEF
      lfound = .true.
      exit State_0d
   else
      indx = indx - 1
   endif
enddo State_0D

! Cycle through 1d state variables
! Note that indices of fields can have varying dimensions.
! WARNING: For the FV, if there's a 1D state variable, then 2 of
!          the _vals will end up being MISSING_R8 at the set_location below.
!          The user needs to figure out what the 'missing' location values
!          should be.
if (.not.lfound) then
State_1D: do i=1,state_num_1d
   nfld = nfld + 1
   if (indx > s_dim_1d(i) ) then
      indx = indx - s_dim_1d(i)
   else
      ! We've found the desired field; now find lat, lon or lev of indx
      dim_name = dim_names(s_dimid_1d(i))
      which_vert = which_vert_1d(i)

      if (dim_name == 'lev') then
         ! FIXME; I think I figured out the get_state_meta_data question about lev_val = real(indx).

         ! This (unlikely) option of a state variable with only a vertical coordinate
         ! cannot be handled by location_mod (yet).  It would require a new which_vert option
         ! which has no horizontal location and/or changes to handle no horizontal location.
         write(string1,*) 'a state variable with only a vertical coordinate cannot be handled ', &
                          'by location_mod (yet).  Guilty field is ',cflds(nfld)
         call error_handler(E_ERR, 'get_state_meta_data', string1, source, revision, revdate);

         ! And what about ilev?
      endif

      if (which_vert == VERTISSURFACE .or. &
          which_vert == VERTISUNDEF) then
         call coord_val(dim_name, indx, lon_val, lat_val, lev_val)

      else
         ! ? Should this be able to return a vertical location for VERTIS{LEVEL,PRESSURE,HEIGHT,VERTISSCALEHEIGHT}?
         !   That is, should users be able to define a CAM(5 and earlier) state variable with those which_verts?
         !   NO
         write(string1,*) 'a 1-D state variable with only a vertical coordinate cannot be handled ', &
                          'by location_mod (yet).  Fix the which_vert_1d of ',cflds(nfld)
         call error_handler(E_ERR, 'get_state_meta_data', string1, source, revision, revdate);

      endif

      lfound = .true.
      exit State_1D
   endif
enddo State_1D
endif

! Cycle through 2d state variables.
! Note that indices of fields can have varying dimensions.
if (.not.lfound) then
State_2D: do i=1,state_num_2d
   nfld = nfld + 1
   slice = s_dim_2d(1,i) * s_dim_2d(2,i)
   if (indx > slice ) then
      indx = indx - slice
   else
      ! We've found the desired field.
      ! Now find lat and/or lon and/or lev of indx if called by assim_tools_mod:filter_assim

         ! FIXME;  Put in a check about whether namelist input (which_vert_Nd) is compatible
      which_vert = which_vert_2d(i)

      ! # second dimension rows to subtract off; temporary value for index_2
      index_2 = (indx -1) / s_dim_2d(1,i)
      index_1 = indx - (index_2 * s_dim_2d(1,i))
      dim_name = dim_names(s_dimid_2d(1,i))
      ! FIXME;  Put in a check about whether namelist input (which_vert_Nd) is compatible with dim_name.
      ! If any dimension is level, it will be the first one.
      if ((dim_name == 'lev' .or. dim_name == 'ilev')  .and. which_vert /= VERTISLEVEL) then
         write(string1,*) 'dim_name is ',dim_name,' but which_vert is not VERTISLEVEL'
         call error_handler(E_ERR, 'get_state_meta_data', string1,source,revision,revdate)
      endif

      if (print_details .and. indx == 1) then
         write(string1,'(A,3I7,A)') 'index_in, index_1, index_2, dim_name', &
                                     index_in, index_1, index_2, dim_name
         write(string2,'(A,I7,A,2I7)') '   s_dim_2d(1:2,',i,') = ',s_dim_2d(1,i), s_dim_2d(2,i)
         call error_handler(E_MSG, 'get_state_meta_data', string1,source,revision,revdate, text2=string2)
      endif

      ! Find the coordinate value (i.e. 270.5) of the first dimension index (i.e. 54)
      if (which_vert == VERTISLEVEL   .or. &
          which_vert == VERTISSURFACE .or. &
          which_vert == VERTISUNDEF) then
         call coord_val(dim_name, index_1, lon_val, lat_val, lev_val)

      else
         ! This section is redundant now that verify_namelist checks for unacceptable values.
         write(string1, *) 'which_vert_2d = ',which_vert_2d(i),', for ',cflds(nfld),  &
              ', cannot be handled in get_state_meta_data->coord_val.'
         call error_handler(E_ERR, 'get_state_meta_data', string1, source, revision, revdate)

      endif

      ! index_2 of the variable in question is 1 more than the # subtracted off to get index_1
      index_2 = index_2 + 1
      dim_name = dim_names(s_dimid_2d(2,i))
      call coord_val(dim_name, index_2, lon_val, lat_val, lev_val)

      lfound = .true.
      exit State_2D
   endif
enddo State_2D
endif

! Cycle through 3d state variables
! Note that indices of fields can have varying dimensions.
if (.not.lfound) then
State_3D: do i=1,state_num_3d
   nfld = nfld + 1
   box = s_dim_3d(1,i) * s_dim_3d(2,i) * s_dim_3d(3,i)
   if (indx > box ) then
      indx = indx - box
   else
      ! We've found the desired field.
      ! Now find lat and/or lon and/or lev of indx if called by assim_tools_mod:filter_assim

      which_vert = which_vert_3d(i)

      ! # of (first x second dimension) slices to subtract off in order to find the current slice
      slice = s_dim_3d(1,i) * s_dim_3d(2,i)
      index_3 = (indx -1) / slice         ! temporary value used to find index_2 and index_1
      index_2 = (indx -1 - (index_3 * slice)) / s_dim_3d(1,i)    ! same for index_2 to find index_1
      index_1 = (indx - (index_3 * slice) - (index_2 * s_dim_3d(1,i)))

      ! Should return dim_name = 'lev' and lev_val = REAL(index_1)
      dim_name = dim_names(s_dimid_3d(1,i))

      ! If any dimension is level, it will be the first one.
      if ((dim_name == 'lev' .or. dim_name == 'ilev')  .and. which_vert /= VERTISLEVEL) then
         write(string1,*) 'dim_name is ',dim_name,' but which_vert is not VERTISLEVEL' 
         call error_handler(E_ERR, 'get_state_meta_data', string1,source,revision,revdate, text2=string2)
      endif

      call coord_val(dim_name, index_1, lon_val, lat_val, lev_val)

      ! index_2 of the variable in question is one more than the number subtracted off to get index_1
      index_2 = index_2 + 1
      dim_name = dim_names(s_dimid_3d(2,i))
      call coord_val(dim_name, index_2, lon_val, lat_val, lev_val)

      ! index_3 of the variable in question is one more than the number subtracted off to get index_1
      index_3 = index_3 + 1
      dim_name = dim_names(s_dimid_3d(3,i))
      call coord_val(dim_name, index_3, lon_val, lat_val, lev_val)

      lfound = .true.
      exit State_3D
   endif
enddo State_3D
endif

! This will malfunction for fields that are filled with MISSING_R8 for lat_val or lon_val.
if (lon_val == MISSING_R8 .or. lat_val == MISSING_R8 ) then
   write(string1, *) 'Field ',cflds(nfld),' has no lon or lat dimension.  ', &
         'What should be specified for it in the call to location?'
   call error_handler(E_ERR, 'get_state_meta_data', string1, source, revision, revdate)
else
   location = set_location(lon_val, lat_val, lev_val, which_vert)
endif

! If the type is wanted, return it
if (present(var_kind)) then
   ! used by call from assim_tools_mod:filter_assim, which wants the DART KIND_
   var_kind = cam_to_dart_kinds(nfld)
endif

end subroutine get_state_meta_data

!-----------------------------------------------------------------------
!>
!> Subroutine ens_mean_for_model
!> makes the ensemble mean available to model_mod
!> and generates the ensemble mean surface pressure arrays needed in get_close_obs.
!> 
!> @param[in] filter_ens_mean
!> The ensemble mean state vector from filter.

!HK not called by the distributed version
subroutine ens_mean_for_model(filter_ens_mean)

real(r8), intent(in) :: filter_ens_mean(:)

if (.not. module_initialized) call static_init_model()

ens_mean = filter_ens_mean

! Fill ps, ps_stagr_lxx.
! newFIXME; it would be nice to have static_init_model call set_ps_arrays,
!           but we would need to call ens_mean_from_model
!           first (from outside of model_mod).
!           But static_init_model is always called before the ensemble mean is set,
!           so that won't work.
!           I could just test allocate_ps before each use of a ps variable
!           and tell the user that ens_mean_for_model must be called before
!           the subroutine (or module entry point) which triggers the message.

!HK Not calling set ps_arrays, waste of communication and memory. This subroutine
! is not called in the distrbuted version of filter.
call set_ps_arrays(ens_mean)

end subroutine ens_mean_for_model

!-----------------------------------------------------------------------
!>
!> Function get_model_size assigns the 'model_size' calculated in static_init_model
!> to the function result 'get_model_size'.

function get_model_size()

integer :: get_model_size

if (.not. module_initialized) call static_init_model()

get_model_size = model_size

end function get_model_size

!-----------------------------------------------------------------------
!>
!> Function get_model_time_step assigns the 'Time_step_atmos' calculated in 
!> static_init_model to the function result 'get_model_time_step'.

function get_model_time_step()

! Returns the time step of the model.

type(time_type) :: get_model_time_step

if (.not. module_initialized) call static_init_model()

get_model_time_step =  Time_step_atmos

end function get_model_time_step

!-----------------------------------------------------------------------
!>
!> Function nc_write_model_atts
!> writes the model-specific attributes to a netCDF file.
!> 
!> @param[in] nc_file_ID      
!>  netCDF file identifier

function nc_write_model_atts( nc_file_ID ) result(ierr)

! Writes the model-specific attributes to a netCDF file.
! TJH Fri Aug 29 MDT 2003

integer, intent(in)  :: nc_file_ID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer :: n_dims, n_vars, n_attribs, unlimited_dim_ID
integer :: member_dim_ID, state_var_dim_ID, time_dim_ID,scalar_dim_ID
integer :: x_var_ID,state_var_ID, state_var_var_ID
integer :: P_id(num_dims)
integer :: i, ifld, dim_id, g_id
integer :: grid_id(grid_num_1d)
character(len=8)  :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10) :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)  :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer           :: values(8)   ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

if (.not. module_initialized) call static_init_model()

! FIXME; bad strategy; start with failure.
ierr = 0     ! assume normal termination

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Make sure nc_file_ID refers to an open netCDF file,
! and then put into define mode.
! nf90_Inquire  returns all but the nc_file_ID; these were defined in the calling routine.
!    More dimensions, variables and attributes will be added in this routine.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write(string1,*) 'nc_file_ID', nc_file_ID
call nc_check(nf90_Inquire(nc_file_ID, n_dims, n_vars, n_attribs, unlimited_dim_ID), &
              'nc_write_model_atts', 'Inquire '//trim(string1))
call nc_check(nf90_Redef(nc_file_ID), 'nc_write_model_atts', 'Redef '//trim(string1))

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! We need the dimension ID for the number of copies
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call nc_check(nf90_inq_dimid(ncid=nc_file_ID, name="copy", dimid=member_dim_ID), &
              'nc_write_model_atts', 'inq_dimid copy')
call nc_check(nf90_inq_dimid(ncid=nc_file_ID, name="time", dimid=  time_dim_ID), &
              'nc_write_model_atts', 'inq_dimid time')

if ( time_dim_ID /= unlimited_dim_Id ) then
  write(string1,*)'Time dimension ID ',time_dim_ID,'must match Unlimited dimension ID ',unlimited_dim_Id
  call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Define the model size, state variable dimension ... whatever ...
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call nc_check(nf90_def_dim(ncid=nc_file_ID, name="StateVariable",  &
                        len=model_size, dimid = state_var_dim_ID),  &
              'nc_write_model_atts', 'def_dim StateVariable')

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Write Global Attributes
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call DATE_AND_TIME(crdate,crtime,crzone,values)

write(str1,'("YYYY MM DD HH MM SS = ",i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "creation_date",str1),        &
              'nc_write_model_atts', 'put_att creation_date'//trim(str1))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "model_revision",revision),   &
              'nc_write_model_atts', 'put_att model_revision'//trim(revision))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "model_revdate",revdate),     &
              'nc_write_model_atts', 'put_att model_revdate'//trim(revdate))
call nc_check(nf90_put_att(nc_file_ID, NF90_GLOBAL, "model","CAM"),               &
              'nc_write_model_atts','put_att model CAM')

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Define the new dimensions IDs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! They have different dimids for this file than they had for caminput.nc
! P_id serves as a map between the 2 sets.
if (print_details .and. output_task0) then
   write(string1,*) 'num_dims = ',num_dims
   write(string2,*) ' dimens,       name,  size, cam dim_id, P[oste]rior id'
   call error_handler(E_MSG, 'nc_write_model_atts', string1,source,revision,revdate, text2=string2)
endif

do i = 1,num_dims
   if (trim(dim_names(i)) /= 'time')  then
      call nc_check(nf90_def_dim(ncid=nc_file_ID, name=trim(dim_names(i)), len=dim_sizes(i),  &
                    dimid=P_id(i)), 'nc_write_model_atts','def_dim '//trim(dim_names(i)))
   else
     P_id(i) = 0
   endif
   if (print_details .and. output_task0) then
      write(string1,'(I5,1X,A13,1X,2(I7,2X))') i,trim(dim_names(i)),dim_sizes(i), P_id(i)
   endif
enddo

call nc_check(nf90_def_dim(ncid=nc_file_ID, name="scalar",   len = 1,   dimid = scalar_dim_ID) &
             ,'nc_write_model_atts', 'def_dim scalar')

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Create the (empty) Coordinate Variables and their attributes
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! grid longitudes, latitudes, levels, and other coordinates.
! grid_id() is filled here; it's the dimid of the desired coordinate *on this P_Diag.nc file*.
! It's used to write coordinates.  ! There's some overlap of names, unfortunately.
! The argument after the 'xxx    ' label is a structure with all the relevant info in it.
! The structures are defined in "Grid fields" and filled by calls to create_grid_1d_instance
! in read_cam_coord.
! CS; ncol doesn't belong here because it's just a dimension, not a coordinate variable.

grid_id = MISSING_I

if (lon%label /= ' ')  then
   dim_id = P_id(lon%dim_id)
   g_id   = find_name('lon',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'lon',lon , dim_id, grid_id(g_id))
endif
if (lat%label /= ' ')  then
   dim_id = P_id(lat%dim_id)
   g_id   = find_name('lat',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'lat',lat , dim_id, grid_id(g_id))
endif
if (lev%label /= ' ')  then
   dim_id = P_id(lev%dim_id)
   g_id   = find_name('lev',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'lev',lev , dim_id, grid_id(g_id))
! Gaussian weights -- because they're there.
endif
if (gw%label /= ' ')  then
   dim_id = P_id(gw%dim_id)
   g_id   = find_name('gw',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'gw',gw  , dim_id, grid_id(g_id))
! Hybrid grid level coefficients, parameters
endif
if (hyam%label /= ' ')  then
   dim_id = P_id(hyam%dim_id)
   g_id   = find_name('hyam',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'hyam',hyam, dim_id, grid_id(g_id))
endif
if (hybm%label /= ' ')  then
   dim_id = P_id(hybm%dim_id)
   g_id   = find_name('hybm',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'hybm',hybm, dim_id, grid_id(g_id))
endif
if (hyai%label /= ' ')  then
   dim_id = P_id(hyai%dim_id)
   g_id   = find_name('hyai',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'hyai',hyai, dim_id, grid_id(g_id))
endif
if (hybi%label /= ' ')  then
   dim_id = P_id(hybi%dim_id)
   g_id   = find_name('hybi',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'hybi',hybi, dim_id, grid_id(g_id))
endif
if (slon%label /= ' ')  then
   dim_id = P_id(slon%dim_id)
   g_id   = find_name('slon',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'slon',slon, dim_id, grid_id(g_id))
endif
if (slat%label /= ' ')  then
   dim_id = P_id(slat%dim_id)
   g_id   = find_name('slat',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'slat',slat, dim_id, grid_id(g_id))
endif
if (ilev%label /= ' ')  then
   dim_id = P_id(ilev%dim_id)
   g_id   = find_name('ilev',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'ilev',ilev, dim_id, grid_id(g_id))
endif
if (P0%label /= ' ')  then
   dim_id = P_id(P0%dim_id)
   ! dim_id here is 0; will that work?  It's what's read in from caminput.nc
   ! If not, then I'll need to (re)define grid_0d_kind, etc.
   g_id   = find_name('P0',grid_names_1d)
   call write_cam_coord_def(nc_file_ID,'P0',P0  , dim_id, grid_id(g_id))
endif

if (print_details .and. output_task0) then
   write(string1,*) '1d field#, grid_id, grid_names_1d'
   call error_handler(E_MSG, 'nc_write_model_atts', string1,source,revision,revdate)
   do i=1,grid_num_1d
      write(string1,*) 'grid_ = ', i, grid_id(i), trim(grid_names_1d(i))
      call error_handler(E_MSG, 'nc_write_model_atts', string1,source,revision,revdate)
   enddo
endif

if ( output_state_vector ) then

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Create attributes for the state vector
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Define the state vector coordinate variable
   call nc_check(nf90_def_var(ncid=nc_file_ID,name="StateVariable", xtype=nf90_int,           &
              dimids=state_var_dim_ID, varid=state_var_var_ID),                                   &
                 'nc_write_model_atts','def_var  state vector')
   call nc_check(nf90_put_att(nc_file_ID, state_var_var_ID, "long_name", "State Variable ID"),   &
                 'nc_write_model_atts','put_att long_name state vector ')
   call nc_check(nf90_put_att(nc_file_ID, state_var_var_ID, "units",     "indexical"),           &
                 'nc_write_model_atts','put_att units state vector ' )
   call nc_check(nf90_put_att(nc_file_ID, state_var_var_ID, "valid_range", (/ 1, model_size /)), &
                 'nc_write_model_atts','put_att valid range state vector ')
   ! Define the actual state vector
   call nc_check(nf90_def_var(ncid=nc_file_ID, name="state", xtype=nf90_real,                 &
              dimids = (/ state_var_dim_ID, member_dim_ID, unlimited_dim_ID /), varid=state_var_ID), &
                 'nc_write_model_atts','def_var state vector')
   call nc_check(nf90_put_att(nc_file_ID, state_var_ID, "long_name", "model state or fcopy"),   &
                 'nc_write_model_atts','put_att long_name model state or fcopy ')

   ! Leave define mode so we can fill
   call nc_check(nf90_enddef(nc_file_ID), 'nc_write_model_atts','enddef ')

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(nc_file_ID, state_var_var_ID, (/ (i,i=1,model_size) /) ),         &
                 'nc_write_model_atts','put_var state_var ')

else

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Create the (empty) Variables and the Attributes
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! 0-d fields
   ifld = 0
   do i = 1,state_num_0d
      ifld = ifld + 1
      call nc_check(nf90_def_var(ncid=nc_file_ID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ member_dim_ID, unlimited_dim_ID /),                             &
                 varid  = x_var_ID),                                                       &
                 'nc_write_model_atts','def_var 0d '//trim(cflds(ifld)))
      call nc_check(nf90_put_att(nc_file_ID, x_var_ID, "long_name", state_long_names(ifld)), &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(nc_file_ID, x_var_ID, "units", state_units(ifld)),          &
                 'nc_write_model_atts','put_att units ')
   enddo

   ! 1-d fields
   do i = 1,state_num_1d
      ifld = ifld + 1
      call nc_check(nf90_def_var(ncid=nc_file_ID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ P_id(s_dimid_1d(1)), member_dim_ID, unlimited_dim_ID /),        &
                 varid  = x_var_ID),                                                       &
                 'nc_write_model_atts','def_var 1d '//trim(cflds(ifld)))
      call nc_check(nf90_put_att(nc_file_ID, x_var_ID, "long_name", state_long_names(ifld)), &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(nc_file_ID, x_var_ID, "units", state_units(ifld)),          &
                 'nc_write_model_atts','put_att units ')
   enddo

   ! 2-d fields
   do i = 1,state_num_2d
      ifld = ifld + 1
      call nc_check(nf90_def_var(ncid=nc_file_ID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ P_id(s_dimid_2d(1,i)), P_id(s_dimid_2d(2,i)),               &
                             member_dim_ID, unlimited_dim_ID /),                             &
                 varid  = x_var_ID),                                                       &
                 'nc_write_model_atts','def_var 2d '//trim(cflds(ifld)))
      call nc_check(nf90_put_att(nc_file_ID, x_var_ID, "long_name", state_long_names(ifld)), &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(nc_file_ID, x_var_ID, "units", state_units(ifld)),          &
                 'nc_write_model_atts','put_att units ')
   enddo

   ! 3-d fields
   do i = 1,state_num_3d
      ifld = ifld + 1
      call nc_check(nf90_def_var                                                              &
           (ncid=nc_file_ID, name=trim(cflds(ifld)), xtype=nf90_real,                           &
            dimids = (/ P_id(s_dimid_3d(1,i)), P_id(s_dimid_3d(2,i)), P_id(s_dimid_3d(3,i)),  &
                        member_dim_ID, unlimited_dim_ID /),                                       &
            varid  = x_var_ID),                                                                 &
                 'nc_write_model_atts','def_var 3d'//trim(cflds(ifld)))
      call nc_check(nf90_put_att(nc_file_ID, x_var_ID, "long_name", state_long_names(ifld)),      &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(nc_file_ID, x_var_ID, "units", state_units(ifld)),               &
                 'nc_write_model_atts','put_att units ')
   enddo

   ! Leave define mode so we can fill variables
   call nc_check(nf90_enddef(nc_file_ID), 'nc_write_model_atts','enddef ')

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Fill the coordinate variables
   ! Each 'vals' vector has been dimensioned to the right size for its coordinate.
   ! The default values of 'start' and 'count'  write out the whole thing.
   ! CS; ncol doesn't belong here because it's just a dimension, not a coordinate variable.
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
   if (lon%label  /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('lon',grid_names_1d)),  lon%vals) &
                    ,'nc_write_model_atts', 'put_var lon')
   if (lat%label  /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('lat',grid_names_1d)),  lat%vals) &
                    ,'nc_write_model_atts', 'put_var lat')
   if (lev%label  /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('lev',grid_names_1d)),  lev%vals) &
                    ,'nc_write_model_atts', 'put_var lev')
   if (gw%label   /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('gw',grid_names_1d)),   gw%vals) &
                    ,'nc_write_model_atts', 'put_var gw')
   if (hyam%label /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('hyam',grid_names_1d)), hyam%vals) &
                    ,'nc_write_model_atts', 'put_var hyam')
   if (hybm%label /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('hybm',grid_names_1d)), hybm%vals) &
                    ,'nc_write_model_atts', 'put_var hybm')
   if (hyai%label /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('hyai',grid_names_1d)), hyai%vals) &
                    ,'nc_write_model_atts', 'put_var hyai')
   if (hybi%label /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('hybi',grid_names_1d)), hybi%vals) &
                    ,'nc_write_model_atts', 'put_var hybi')
   if (slon%label /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('slon',grid_names_1d)), slon%vals) &
                    ,'nc_write_model_atts', 'put_var slon')
   if (slat%label /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('slat',grid_names_1d)), slat%vals) &
                    ,'nc_write_model_atts', 'put_var slat')
   if (ilev%label /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('ilev',grid_names_1d)), ilev%vals) &
                    ,'nc_write_model_atts', 'put_var ilev')
   if (P0%label   /= ' ') &
       call nc_check(nf90_put_var(nc_file_ID, grid_id(find_name('P0',grid_names_1d)),   P0%vals) &
                    ,'nc_write_model_atts', 'put_var P0')

endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Flush the buffer and leave netCDF file open
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call nc_check(nf90_sync(nc_file_ID),'nc_write_model_atts', 'sync ')

end function nc_write_model_atts

!-----------------------------------------------------------------------
!>
!> Function nc_write_model_vars
!> writes the model-specific variables to a netCDF file.
!> 
!> @param[in] nc_file_ID
!> netCDF file identifier
!> 
!> @param[in] statevec(:)
!> The state vector to be written to 'nc_file_ID'
!> 
!> @param[in] copyindex
!> The 'copy' in the file into which the state vector will be written
!> 
!> @param[in] timeindex
!> The time slot in the file, into which the state vector will be written

function nc_write_model_vars( nc_file_ID, statevec, copyindex, timeindex ) result(ierr)

! Writes the model-specific variables to a netCDF file
! TJH 25 June 2003

integer,  intent(in) :: nc_file_ID
real(r8), intent(in) :: statevec(:)
integer,  intent(in) :: copyindex
integer,  intent(in) :: timeindex
integer   :: ierr

type(model_type) :: Var

integer :: n_dims, n_vars, n_attribs, unlimited_dim_ID
integer :: state_var_ID, nc_var_ID
integer :: ifld, i

character(len=8) :: cfield

if (.not. module_initialized) call static_init_model()

! FIXME; bad strategy; start with failure.
ierr = 0     ! assume normal termination

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! make sure nc_file_ID refers to an open netCDF file,
! then get all the Variable ID's we need.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call nc_check(nf90_Inquire(nc_file_ID, n_dims, n_vars, n_attribs, unlimited_dim_ID), &
              'nc_write_model_vars','Inquire ')

if ( output_state_vector ) then

   call nc_check(nf90_inq_varid(nc_file_ID, "state", state_var_ID),'nc_write_model_vars ','inq_varid state' )
   call nc_check(nf90_put_var(nc_file_ID, state_var_ID, statevec,  &
                start=(/ 1, copyindex, timeindex /)),'nc_write_model_vars ','put_var state')

else

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Fill the variables
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   call init_model_instance(Var)     ! Explicity released at end of routine.

   call vector_to_prog_var(statevec,  Var)

   ifld = 0
   ZeroDVars: do i = 1, state_num_0d
      ifld = ifld + 1
      cfield = trim(cflds(ifld))
      call nc_check(nf90_inq_varid(nc_file_ID, cfield, nc_var_ID),       &
                    'nc_write_model_vars ','inq_varid 0d '//cfield)
      call nc_check(nf90_put_var(nc_file_ID, nc_var_ID, Var%vars_0d(i),  &
                                 start=(/ copyindex, timeindex /) ),   &
                    'nc_write_model_vars ','put_var 0d '//cfield)
   enddo ZeroDVars

   ! 'start' and 'count' are needed here because Var%vars_Nd are dimensioned by the largest
   ! values for the dimensions of the rank N fields, but some of the fields are smaller than that.
   OneDVars: do i = 1, state_num_1d
      ifld = ifld + 1
      cfield = trim(cflds(ifld))
      call nc_check(nf90_inq_varid(nc_file_ID, cfield, nc_var_ID),                                    &
                    'nc_write_model_vars ','inq_varid 1d '//cfield)
      call nc_check(nf90_put_var(nc_file_ID, nc_var_ID,                                               &
                    Var%vars_1d(1:s_dim_1d(  i), i),                                              &
                    start=   (/ 1              ,copyindex, timeindex /),                          &
                    count=   (/   s_dim_1d(  i),1        , 1/) ),                                 &
                    'nc_write_model_vars ','put_var 1d '//cfield)
   enddo OneDVars

   ! Write out 2D variables as 2 of (lev,lon,lat), in that order, regardless of caminput.nc
   ! coordinate order
   ! The sizes can be taken from s_dim_2d, even though the s_dimid_2d don't pertain to this
   ! P_Diag.nc file, because the dimids were mapped correctly in nc_write_model_atts.
   TwoDVars: do i = 1, state_num_2d
      ifld = ifld + 1
      cfield = trim(cflds(ifld))
      call nc_check(nf90_inq_varid(nc_file_ID, cfield, nc_var_ID),                                    &
                    'nc_write_model_vars ','inq_varid 2d '//cfield)
      call nc_check(nf90_put_var(nc_file_ID, nc_var_ID,                                               &
                    Var%vars_2d(1:s_dim_2d(1,i),1:s_dim_2d(2,i), i),                              &
                    start=   (/ 1              ,1              , copyindex, timeindex /),         &
                    count=   (/   s_dim_2d(1,i),  s_dim_2d(2,i), 1        , 1/) ),                &
                    'nc_write_model_vars ','put_var 2d '//cfield)
   enddo TwoDVars

   ! Write out 3D variables as (lev,lon,lat) regardless of caminput.nc coordinate order
   ThreeDVars: do i = 1,state_num_3d
      ifld = ifld + 1
      cfield = trim(cflds(ifld))
      call nc_check(nf90_inq_varid(nc_file_ID, cfield, nc_var_ID),                                    &
                    'nc_write_model_vars ','inq_varid 3d '//cfield)
      call nc_check(nf90_put_var(nc_file_ID, nc_var_ID,                                               &
                 Var%vars_3d(1:s_dim_3d(1,i),1:s_dim_3d(2,i),1:s_dim_3d(3,i),i)                   &
                 ,start=   (/1              ,1              ,1              ,copyindex,timeindex/)&
                 ,count=   (/  s_dim_3d(1,i),  s_dim_3d(2,i),  s_dim_3d(3,i),1        ,1 /) ),    &
                    'nc_write_model_vars ','put_var 3d '//cfield)
   enddo ThreeDVars

endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Flush the buffer and leave netCDF file open
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call nc_check(nf90_sync(nc_file_ID),'nc_write_model_vars ','sync ')

call end_model_instance(Var)   ! should avoid any memory leaking

end function nc_write_model_vars


! End of Module I/O

!#######################################################################

! model_interpolate section

!-----------------------------------------------------------------------
!>
!> Subroutine model_interpolate
!> Interpolates the provided state vector (on model grid points) to an arbitrary
!> location in the atmosphere (e.g. where an observation is).
!> 
!> @param[in] :: st_vec(:)
!> The state vector which will be interpolated to 'location'
!> 
!> @param[in] :: location
!> The DART location_type 'location' of the desired state estimate.
!> 
!> @param[in] :: obs_kind
!> The DART KIND of the variable being estimated.
!> 
!> @param[out] :: interp_val
!> The state estimate of the 'obs_kind' at 'location'
!> 
!> @param[out] :: istatus
!> A flag to signal the success of the interpolation.

subroutine model_interpolate_distrib(state_ens_handle, location, obs_kind, istatus, interp_val)

! This subroutine is now a short routine that calls
! either a rectangular grid version for eul/FV
! or non-rectangular for cubed-sphere code.
! This does get KINDs from filter, not specific obs TYPEs.

! Model_interpolate must return a positive value for istatus for a failure.
! 0 means success, negative values are reserved for DART internal use.

real(r8),            intent(in) :: state_ens_handle
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_kind
real(r8),           intent(out) :: interp_val(:)
integer,            intent(out) :: istatus(:)

! Local variables for CAM-SE
! cell_corners(1) = MISSING_I tells interp_cubed_sphere to find the corners of the cell,
! instead of using ones passed to it.
integer  :: cell_corners(4)
real(r8) :: l, m

! Don't initialize these in the declaration statement or they'll get the save attribute
! and won't be re-initialized each time model_interpolate is entered.
cell_corners = MISSING_I
l = MISSING_R8
m = MISSING_R8

! FIXME; In future DARTs it may be useful to return the DART KIND too.
!        also convert to a field name (DART subroutine (get_raw_...?)).

if (.not. module_initialized) call static_init_model()

! FIXME; Tim  test for ob being out of bounds (horizontally, vertically?)
! and return if it is.
! But interp_yyy might need to be called anyway (in the future), to get the value,
! even if it won't be assimilated.
! Also, lat bounds could be enforced here with a small amount of code,
! but enforcing vertical bounds would require bringing lots of code from interp_yyy up here,
! and have if-tests to separate out the lonlat from the cubed_sphere.

! FIXME; Tim Also add an argument to inter_XXX to tell it what to do when the ob
! is out of bounds, but still calculatable.

if (l_rectang) then
   call interp_lonlat_distrib(state_ens_handle, location, obs_kind, istatus, interp_val)
else
   call error_handler(E_ERR, 'no cubed sphere version of RMA')
   !call interp_cubed_sphere(st_vec, location, obs_kind, interp_val, istatus, cell_corners, l, m)
endif

end subroutine model_interpolate_distrib

!-----------------------------------------------------------------------
! HK why is this recursive?  Poor abstraction?
recursive subroutine interp_cubed_sphere(st_vec, obs_loc, obs_kind, interp_val, istatus, cell_corners, l, m)

! Find the cell that encloses an ob at 'obs_loc'
! and interpolate the values of obs_kind from the cell's corners to that location.
! The obs_kinds passed in here are always explicitly KIND_xxxx parameters,
! except for the call from model_interpolate, which comes from filter.
! cell_corners, l, and m are inout because the call to interp_cubed_sphere in model_height
! can re-use those values many times.

real(r8),            intent(in)    :: st_vec(:)
type(location_type), intent(in)    :: obs_loc
integer,             intent(in)    :: obs_kind
real(r8),            intent(out)   :: interp_val
integer,             intent(out)   :: istatus
integer,             intent(inout) :: cell_corners(4)   ! node numbers of the corners of the enclosing cell
real(r8),            intent(inout) :: l
real(r8),            intent(inout) :: m

! FIXME; In future DARTs it may be useful to return the DART KIND too.
!        also convert to a field name (DART subroutine (get_raw_...?)).

integer  :: vstatus, i, closest

integer               :: s_type, s_type_1d, s_type_2d
character(len=8)      :: col_name, lev_name
real(r8)              :: lon_lat_lev(3), vals(4)

! Start with failure, then change to success as warranted.
istatus = 1
interp_val = MISSING_R8

! Get the observation (horizontal) position, in degrees.
lon_lat_lev = get_location(obs_loc)

! Get horizontal grid coordinate names.

! Set [ncol,lev] names to defaults, which may be overwritten for variables in the state vector,
! but not for other acceptable variables (3D pressure, surface elevation, ...?)
col_name = 'ncol'

! ? How to separate the 3D 'other' variables from 2D 'other' variables?
!   Can't do it automatically/generically because they're not part of state vector
!   and that info isn't coming from DART.
if (obs_kind == KIND_SURFACE_ELEVATION) then
   lev_name = 'none'
elseif (obs_kind == KIND_PRESSURE) then
   lev_name = 'lev'
endif

! check whether model_mod can interpolate the requested variable
! Pressure (3d) can't be specified as a state vector field (so s_type will = MISSING_I),
! but can be calculated for CAM, so obs_kind = KIND_PRESSURE is acceptable.
s_type = dart_to_cam_types(obs_kind)
if (s_type == MISSING_I .and. &
   (obs_kind /= KIND_PRESSURE) .and.  (obs_kind /= KIND_SURFACE_ELEVATION)) then
   ! FIXME; Should be warning about the first one only, don't die.
   !        Could that be done with an array that's 'save'd from call to call?
   write(string1,*) 'Wrong type of obs = ', obs_kind
   call write_location(0, obs_loc, charstring=string2)
   call error_handler(E_WARN, 'interp_cubed_sphere', string1,source,revision,revdate, text2=string2)
   return
endif

! There can't be any 0d ob fields, so subtract off the 0d state fields from 
! the list of state fields.
!    What about earth rotation obs?  Are they converted into 2D fields?

! Positions within the spatially rank 1 and 2 fields.
s_type_1d = s_type - state_num_0d
s_type_2d = s_type_1d - state_num_1d

if (s_type == MISSING_I .and. &
   (obs_kind == KIND_PRESSURE) .or. (obs_kind == KIND_SURFACE_ELEVATION)) then
   ! use defaults col_name set above

elseif (s_type <= state_num_0d ) then
   ! error; can't deal with observed variables that are 0D in model_mod.
   write(string1,*) 'DART cannot handle 0d observations of ', cflds(s_type), &
        ' because DART requires a (lon,lat) location for each observation '
   write(string2,*) 'Skipping this observation'
   call error_handler(E_WARN, 'interp_cubed_sphere', string1,source,revision,revdate, text2=string2)
   return

elseif (s_type_1d > 0 .and. s_type_1d <= state_num_1d) then
   col_name = dim_names(s_dimid_1d(s_type_1d))

elseif (s_type_2d > 0 .and. s_type_2d <= state_num_2d) then
   col_name = dim_names(s_dimid_2d(2,s_type_2d))
   lev_name = dim_names(s_dimid_2d(1,s_type_2d))

else
   ! There are no spatially rank 3 fields in the CAM-SE initial file.
   write(string1,*) 'Unexpected state type value, s_type = ', s_type
   call error_handler(E_MSG, 'interp_cubed_sphere', string1,source,revision,revdate)
   return
endif

! Need the node names of the corners of the cell which contains the observation.
! Also returns the location of the ob in the unit square space (l,m),
! since it will be calculated to find the cell, and is then needed for the interpolation.
! But don't bother calling this if cell_corners (and l and m) were found in
! a previous call and passed into this routine.

if (cell_corners(1) == MISSING_I) then
   call coord_ind_cs(obs_loc, obs_kind, .false., closest, cell_corners, l, m)
! else
! !    cell_corners, l, and m, are being passed in from convert_vert:model_heights
endif

if (print_details .and. output_task0) then
   write(string1,'(A,4I8)') 'cell_corners = ',(cell_corners(i), i=1,4)
   write(string2,'(A,I8,A,1p2E20.12) ')'{lon,lat}(',cell_corners(4),') = ', &
        lon%vals(cell_corners(4)), lat%vals(cell_corners(4))
   call error_handler(E_MSG, 'interp_cubed_sphere', string1,text2=string2)
endif

! The interpolation.
! First interpolate the field in the vertical at each of the 4 corners
! to the height of the ob.
! The subroutines and arrays appear to want indices for the lon and lat dimensions,
! while the cubed sphere has only the ncol horizontal dimension.
! This is handled by passing the ncol index as 'lon_index', and 1 as 'lat_index'.
! Then get_val (the bottom of the calling trees) uses these correctly for the cubed sphere.

if (obs_kind == KIND_SURFACE_ELEVATION) then
   ! Acceptable KIND that's not in the state vector
   ! convert from geopotential height to real height in meters
   vals(1) = phis(cell_corners(1),1) / gravity_const
   vals(2) = phis(cell_corners(2),1) / gravity_const
   vals(3) = phis(cell_corners(3),1) / gravity_const
   vals(4) = phis(cell_corners(4),1) / gravity_const
   vstatus = 0

elseif (vert_is_level(obs_loc)) then
   ! Case 1: model level specified in vertical
   ! Pobs

      call get_val_level(st_vec, cell_corners(1),1, nint(lon_lat_lev(3)), obs_kind, vals(1), vstatus)
   if (vstatus /= 1) &
      call get_val_level(st_vec, cell_corners(2),1, nint(lon_lat_lev(3)), obs_kind, vals(2), vstatus)
   if (vstatus /= 1) &
      call get_val_level(st_vec, cell_corners(3),1, nint(lon_lat_lev(3)), obs_kind, vals(3), vstatus)
   if (vstatus /= 1) &
      call get_val_level(st_vec, cell_corners(4),1, nint(lon_lat_lev(3)), obs_kind, vals(4), vstatus)
   ! Pobs end

elseif (vert_is_pressure(obs_loc)) then
   ! which_vert is pressure for this obs
      call get_val_pressure(st_vec,cell_corners(1),1,lon_lat_lev(3),obs_kind,vals(1),vstatus)
   if (vstatus /= 1) &
      call get_val_pressure(st_vec,cell_corners(2),1,lon_lat_lev(3),obs_kind,vals(2),vstatus)
   if (vstatus /= 1) &
      call get_val_pressure(st_vec,cell_corners(3),1,lon_lat_lev(3),obs_kind,vals(3),vstatus)
   if (vstatus /= 1) &
      call get_val_pressure(st_vec,cell_corners(4),1,lon_lat_lev(3),obs_kind,vals(4),vstatus)

elseif (vert_is_height(obs_loc)) then

   !HK model_heights is called 4 times on the obs_loc here. Also why call the interpolation
   ! routine over and over again?

   ! which_vert is height for this obs
      call get_val_height(st_vec, cell_corners(1),1, lon_lat_lev(3), obs_loc, obs_kind, vals(1), vstatus)
   if (vstatus /= 1) &
      call get_val_height(st_vec, cell_corners(2),1, lon_lat_lev(3), obs_loc, obs_kind, vals(2), vstatus)
   if (vstatus /= 1) &
      call get_val_height(st_vec, cell_corners(3),1, lon_lat_lev(3), obs_loc, obs_kind, vals(3), vstatus)
   if (vstatus /= 1) &
      call get_val_height(st_vec, cell_corners(4),1, lon_lat_lev(3), obs_loc, obs_kind, vals(4), vstatus)


elseif (vert_is_surface(obs_loc)) then
   ! location_mod:interactive_location asks for surface obs to have vertical coord = ps(hPa)
   ! The 'lev' argument is set to 1 because there is no level for these types, and 'lev' will be
   ! ignored.
                     call get_val(st_vec, cell_corners(1),1, 1, obs_kind, vals(1), vstatus)
   if (vstatus /= 1) call get_val(st_vec, cell_corners(2),1, 1, obs_kind, vals(2), vstatus)
   if (vstatus /= 1) call get_val(st_vec, cell_corners(3),1, 1, obs_kind, vals(3), vstatus)
   if (vstatus /= 1) call get_val(st_vec, cell_corners(4),1, 1, obs_kind, vals(4), vstatus)

! This needs to be at the end of the block.  Otherwise, it short circuits GPS
! which asks for pressures on heights.
! elseif (obs_kind == KIND_PRESSURE) then
!    ! Calculate pressures from surface pressures and A and B coeffs.
!    ! Pobs
!       write(string1,'(A)') 'No code available yet for obs_kind = KIND_PRESSURE '
!       call error_handler(E_ERR, 'interp_cubed_sphere', string1)

! Need option for vert_is_scale_height
elseif (vert_is_scale_height(obs_loc)) then
   write(string1,*)'Scale height is not an acceptable vert coord yet.  Skipping observation'
   call error_handler(E_WARN, 'interp_cubed_sphere', string1,source,revision,revdate)
   return

! Need option for vert_is_undefined

else
   write(string1,*) '   No vert option chosen!'
   call error_handler(E_WARN, 'interp_cubed_sphere', string1,source,revision,revdate)
   return
endif

! Then interpolate horizontally to the (lon,lat) of the ob.
! The following uses Jeff's recommended 'generalized quadrilateral interpolation', as in
! http://www.particleincell.com/2012/quad-interpolation/.
! Most of the work is done in create_cs_grid_arrays() and coord_ind_cs().

! Interpolate from the cell's corners to the ob location on the unit square.
! This is done by weighting the field at each corner by the rectangular area
! ((l,m) space) diagonally across the ob location from the corner.
! AKA 'linear area weighting'.

! The vals indices are consistent with how mapping of corners was done,
! and how cell_corners was assigned.
if (vstatus == 1) then
   if (print_details) then
      write(string1,'(A,2F10.6,1pE20.12)') 'istatus = 1, no interpolation'
      call error_handler(E_MSG, 'interp_cubed_sphere', string1)
   endif
   return
else
   if (abs(lon_lat_lev(2)) > max_obs_lat_degree) then
      ! Define istatus to be a combination of vstatus (= 0 or 2 (for higher than highest_obs...))
      ! and whether the ob is poleward of the limits set in the namelist (+ 4).
      istatus = 10*vstatus + 4
   else
      istatus = vstatus
   endif
endif

interp_val = vals(2) *           l *          m &
           + vals(1) * (1.0_r8 - l)*          m &
           + vals(4) * (1.0_r8 - l)*(1.0_r8 - m) &
           + vals(3) *           l *(1.0_r8 - m)

if (print_details ) then
   write(string1,'(A,2F10.6,1pE20.12)') ' l,m, interpolated val = ', &
         l,m,interp_val
   call error_handler(E_MSG, 'interp_cubed_sphere', string1)
endif


end subroutine interp_cubed_sphere

!-----------------------------------------------------------------------

subroutine unit_square_location(cell, closest, location, lon_o,lat_o, found_cell, origin, l,m)

! Subroutine based on http://www.particleincell.com/2012/quad-interpolation/.
! The idea is to derive a mapping from any convex quadrilateral(x,y) onto a unit square (l,m).
! Also map the location of the ob onto that square.
! This is a bilinear interpolation;
! x = a0 + a1*l*m + a2*m + a3*l
! y = b0 + b1*l*m + b2*m + b3*l
! so does not take into account the curvature of the quadrilateral on the sphere.
!
! That has been handled by the intermediate mapping from (lon,lat) to a flat planar
! coordinate system (x,y). The locations of the corners/nodes are converted to
! the distances and directions from one node to the other three.  See create_cs_grid_arrays.
! Distances and directions relative to the origin node are preserved, but distances and
! directions between 2 non-origin points are slightly distorted.
! Even these small errors are avoided by defining a planar coordinates system for each corner
! of each cell.
! Then the ob is never near the 'far edges', where distortion could be a problem.
!
! A higher order method exists (Nagata 2005: Simple Local Interpolation of Surfaces
! Using Normal Vectors) to map curved quadrilaterals onto the unit square,
! but the inverse map cannot be done analytically(?), so is not developed here.

integer,             intent(in)    :: cell
integer,             intent(in)    :: closest
type(location_type), intent(in)    :: location
real(r8),            intent(in)    :: lon_o
real(r8),            intent(in)    :: lat_o
logical,             intent(inout) :: found_cell
integer,             intent(out)   :: origin
real(r8),            intent(out)   :: l
real(r8),            intent(out)   :: m

! Observation location in the planar space.
real(r8) :: x_o, y_o

real(r8) :: angle, d, bearing_o  ! Locations in polar coordinate space (bearing,distance).
real(r8) :: aa, bb, cc           ! Coefficients of quadratic equation for m.
real(r8) :: det, m1, m2          ! Determinant and roots.
logical  :: neg_root             ! helpful logical variable to store usefulness of the -root.
real(r8) :: m_neg, l_neg         ! Potential alternate solutions to the m quadratic equation
integer  :: oc(1)

m1    = MISSING_R8  ! first  root returned by solve_quadratic
m2    = MISSING_R8  ! second root returned by solve_quadratic
l     = MISSING_R8  ! unit square abscissa ('x' coord)
m     = MISSING_R8  ! unit square ordinate ('y' coord)
l_neg = MISSING_R8  ! same but for the negative root of the m quadratic equation.
m_neg = MISSING_R8  ! same
neg_root = .false.

! Map the location of the ob into the planar space

! Figure out which corner (1,2,3 or 4) of cell is the closest to the ob,
! by comparing the names of the corners to the name of the node/corner closest
! to the ob, which was passed in.
! Used to get the correct x_ax_bearing and a and b coeffs (from the cs_grid_file).

oc = minloc(corners(cell,:), mask = (corners(cell,:) == closest))
origin = oc(1)

! The bearing of the observation relative to the origin/closest corner.
bearing_o = bearing(lon_rad(closest),lat_rad(closest),lon_o*DEG2RAD,lat_o*DEG2RAD )

! Calculate the difference of the ob bearing from x_axis of this cell.
! The order is opposite of what might be expected because bearings are measured clockwise,
! while angles are measured counterclockwise.
angle = x_ax_bearings(origin,cell) - bearing_o

! Normalize angle to -pi<angle<pi.
angle = mod(angle,PI) - PI*int(angle/PI)

! Calculate the distance from this cell's origin
! and then the (x,y) coordinates of the observation.
d = get_dist(cs_locs(closest), location, no_vert=.true.)
x_o = d * cos(angle)
y_o = d * sin(angle)

! Coefficients of the quadratic equation for m for this cell.
aa =   a(1,origin,cell)*b(2,origin,cell) &
     - a(2,origin,cell)*b(1,origin,cell)

bb =   a(3,origin,cell)*b(2,origin,cell) &
     - a(1,origin,cell)*y_o              &
     + b(1,origin,cell)*x_o

cc = - a(3,origin,cell)*y_o

! Calculate m from the binomial equation, given the quadratic equation coefficients.
call solve_quadratic(aa,bb,cc,m1,m2)
! newFIXME: Simplify this subroutine?
! There can only be one mapping from the (x,y) space to the unit square.
! One of the (potentially) 2 'm's generated here should be eliminated by the l calculation.

if (m1 == MISSING_R8 .and. m2 == MISSING_R8) then
   ! determinant was < 0.
   write(string1,'(A,I6,1X,1p4E12.4)') 'm b^2-4ac <0: cell, angle, d, x_o, y_o',cell, angle, d, x_o, y_o
   call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
   write(string1,'(A,1X,1p4E12.4)') '   a: a(1)*b(2) - a(2)*b(1) : ',  &
                            a(1,origin,cell),b(2,origin,cell), &
                            a(2,origin,cell),b(1,origin,cell)
   write(string2,'(A,1X,1p4E12.4)') '   b: ', bb
   write(string3,'(A,1X,1p4E12.4)') '   c: a(3)*y_o : ', a(3,origin,cell), y_o
   call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate, &
                      text2=string2, text3=string3)
else
   if (aa > 0.0_r8) then
      ! Only m values (roots) between 0 and 1 mean that the ob is in this cell.
      if (m1 >=0 .and. m1 <= 1) then
         m = m1
      elseif (m2 >=0 .and. m2 <= 1) then
         m = m2
      else
         ! Neither root is a map.  Leave m as MISSING_R8
      endif
   elseif (m1 /= MISSING_R8 .and. m2 == MISSING_R8 ) then
      ! Cell is square; solved the linear equation    m*bb + cc = 0
      m = m1
   else
      ! aa < 0; Either both or neither roots yield m>0.
      ! Start with the +root.
!       m = (-bb + sqrt_det)/(2.0_r8*aa)
      m = m1

      if (bb > 0.0_r8) then
         ! Both roots yield m > 0.
         if (m > 1.0_r8) then
            ! The +root didn't yield a usable m.  Try the -root.
            ! m = (-bb - sqrt_det)/(2.0_r8*aa)
            m = m2
         else
            ! It could be that both roots yield a usable m.  Keep track of both
            ! (for testing/debugging only).
            m_neg = m2
         endif

      elseif (bb < 0.0_r8) then
         ! aa < 0 and bb < 0 yields no roots with m>0.
         write(string1,'(A,I6,A)') 'aa < 0 and bb < 0: It appears that cell ',cell,           &
              ' is a highly distorted quadrilateral'
         write(string2,'(A)')                                                                 &
              'and no mapping is possible.  bb = a(3)*b(2) - a(1)*y_o + b(1)*x_o: '
         write(string3,'(1p,(3X,2E12.4))')                                                     &
              a(3,origin,cell),b(2,origin,cell),                                              &
              a(1,origin,cell),y_o,                                                           &
              b(1,origin,cell),x_o
         call error_handler(E_ERR, 'unit_square_location', string1,source,revision,revdate,  &
                            text2=string2, text3=string3)

      elseif (bb == 0.0_r8) then
         ! aa < 0 and bb = 0  should be excluded by the non-negativeness test on det, above.
         write(string1,'(A,1p,2(1x,E12.4))') &
              'aa < 0 and bb = 0 should have been excluded ',aa,bb
         call error_handler(E_ERR, 'unit_square_location', string1,source,revision,revdate)

      endif

   endif

endif

! If m (and maybe m_neg) is out of the possible range, return to calling program
! with found_cell still false.
if ( m < 0.0_r8 .or. m > 1.0_r8) then
   if (.not.found_cell) then
      if (m_neg < 0.0_r8 .or. m_neg > 1.0_r8) then
         ! This includes m_neg == MISSING_R8, due to only m being assigned above.
         return
      endif
   ! ? Can these 2 sections ever be entered?
   else
      ! Exit with error if m is outside valid range.
      write(string1, *) 'location of ob in unit square is out of bounds m = [0,1] ',m, &
                        'but status is "found"'
      call error_handler(E_ERR, 'unit_square_location', string1, source, revision, revdate)
   endif
endif

! Use m to calculate the other unit square coordinate value, 'l'.
det = a(3,origin,cell) + a(1,origin,cell) * m
if (det /= 0.0_r8) then
   l = (x_o - a(2,origin,cell)*m) / det
else
   write(string1,'(A,I6,1X,1p4E12.4)') 'l denominator = 0: cell, angle, d, x_o, y_o',cell, angle, d, x_o, y_o
   write(string2,'(A,1X,1p4E12.4)') '  a(3) + a(1)*m = 0 : ', a(3,origin,cell), a(1,origin,cell),m
   call error_handler(E_ERR, 'unit_square_location', string1,source,revision,revdate, text2=string2)
endif

! Repeat for the -root, if it is a possibility.
if (m_neg /= MISSING_R8) then
   det = (a(3,origin,cell) + a(1,origin,cell)*m_neg)
   if (det /= 0.0_r8) then
      l_neg = (x_o -a(2,origin,cell)*m_neg) / det
   else
      write(string1,'(A,I6,1X,1p4E12.4)') 'l_neg denominator = 0: cell, angle, d, x_o, y_o', &
           cell, angle, d, x_o, y_o
      write(string2,'(A,1X,1p4E12.4)') '  a(3) + a(1)*m = 0 : ', a(3,origin,cell), a(1,origin,cell),m
      call error_handler(E_ERR, 'unit_square_location', string1,source,revision,revdate, text2=string2)
   endif

   ! Informational output, if the observation is exactly on the m-axis
   if (l_neg == 0.0_r8 .and. output_task0) then
      write(string1,'(A,I6,1X,1p4E12.4)') 'l_neg cell, x_o - a(2)*m = ',cell, x_o ,a(2,origin,cell),m
      call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
   endif

endif

! Informational output, if the observation is exactly on the m-axis
if (l == 0.0_r8 .and. output_task0) then
   write(string1,'(A,I6,1X,1p4E12.4)') 'Ob is on x-axis: l-cell, x_o - a(2)*m = ',cell, x_o ,a(2,origin,cell),m
   call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
endif

! If l (and maybe l_neg) is out of the possible range, return to calling program
! with found_cell still false.
if (l < 0.0_r8 .or. l > 1.0_r8) then
   if (.not.found_cell) then
      if (l_neg < 0.0_r8 .or. l_neg > 1.0_r8) then
         ! This includes m_neg == MISSING_R8, due to only m being assigned above
         ! Return with found_cell still = failure (0) to test the next cell.
         return
      endif
   ! ? Can these 2 sections ever be entered?
      ! Exit with error if l is outside valid range.
   else
      ! Exit with error if l is outside valid range.
      write(string1, *) 'location of ob in unit square is out of bounds l = [0,1] ',l, &
                        'but status is "found"'
      call error_handler(E_ERR, 'unit_square_location', string1, source, revision, revdate)
   endif
endif

! If we get this far, then this cell contains the ob.

! But which root(s) of the m quadratic equation led to the mapping?
! Put the right values in l and m.
neg_root = m_neg >= 0.0_r8 .and. m_neg <= 1.0_r8 .and. &
           l_neg >= 0.0_r8 .and. l_neg <= 1.0_r8
if (m >= 0.0_r8 .and. m <= 1.0_r8 .and. &
    l >= 0.0_r8 .and. l <= 1.0_r8 ) then
   ! Both roots yield a good mapping.
   if (neg_root) then
      write(string1, *) 'BOTH roots of the m quadratic yield usable mappings.  The +root is being used.'
      call error_handler(E_MSG, 'unit_square_location', string1, source, revision, revdate)
   endif

elseif ( neg_root) then
   ! The -root yields a good mapping.  Pass along the -root m and l.
   m = m_neg
   l = l_neg
   write(string1, *) 'The negative root of the m quadratic yielded the only usable mapping.'
   call error_handler(E_MSG, 'unit_square_location', string1, source, revision, revdate)
endif

! Return with found_cell = true; success.
found_cell = .true.

end subroutine unit_square_location

!-----------------------------------------------------------------------

subroutine solve_quadratic(a, b, c, r1, r2)

real(r8), intent(in)  :: a
real(r8), intent(in)  :: b
real(r8), intent(in)  :: c
real(r8), intent(out) :: r1
real(r8), intent(out) :: r2

real(r8) :: scaling, as, bs, cs, disc

r1 = MISSING_R8
r2 = MISSING_R8

! Scale the coefficients to get better round-off tolerance
scaling = max(abs(a), abs(b), abs(c))
as = a / scaling
bs = b / scaling
cs = c / scaling

if (abs(as) < epsilon(as)) then
   ! Solve the linear equation bs*r + cs = 0
   r1 = -cs / bs
else
   ! Get discriminant of scaled equation
   disc = bs * bs - 4.0_r8 * as * cs
   if (disc >= 0.0_r8) then

      ! Calculate the largest root (+ or - determined by sign of bs)
      ! Handling of bs = 0 different from pre-review code
      !    if(bs > 0.0_r8) then
      if(bs >= 0.0_r8) then
         r1 = (-bs - sqrt(disc)) / (2.0_r8 * as)
      else
         r1 = (-bs + sqrt(disc)) / (2.0_r8 * as)
      endif

      ! Compute the second root given the larger (not most positive) one
      if (r1 == 0.0_r8) then
         ! The b AND c must have been 0: solved the equation a*r1^2 = 0 above
         ! and there's no 2nd root.
         r2 = 0.0_r8
      else
         ! 'as' and 'r1' have been tested for 0.
         r2 = cs / (as * r1)
      endif
   endif
endif

end subroutine solve_quadratic

!-----------------------------------------------------------------------
! HK why is this recursive (it is model_heights that calls this)?
recursive subroutine interp_lonlat_distrib(state_ens_handle, obs_loc, obs_kind, istatus, interp_val)

! Find the 4 corners of the lon-lat grid cell that encloses an ob at 'obs_loc'
! and interpolate the values of obs_kind to that location.

! istatus   meaning                  return expected obs?   assimilate?
! 0         obs and model are fine;  yes                    yes
! 1         fatal problem;           no                     no
! 2         exclude valid obs        yes                    no   (ob > ! highest_obs_X)
! 3         unfamiliar obs type      no                     no
! 4         ob excl by namelist(lat) yes                    no
! NM        2 digit number means more than one namelist reason to exclude from assim.

! Any value > 0 will not be assimilated (---> QC non-0).
! Do we want some istatus values to tell filter to evaluate (---> QC of 1)?
! That would be nice, but filter has no convention for understanding non-0
! values from model_mod (from all the available models).  So all non-0 values of
! istatus ---> QC = 4.


type(ensemble_type), intent(in) :: state_ens_handle
type(location_type), intent(in) :: obs_loc
integer,             intent(in) :: obs_kind
integer,            intent(out) :: istatus(:)
real(r8),           intent(out) :: interp_val(:)

! FIXME; In future DARTs it may be useful to return the DART KIND too.
!        also convert to a field name (DART subroutine (get_raw_...?)).

integer  :: i
real(r8) :: bot_lon, top_lon, delta_lon,                                &
            lon_below, lat_below, lat_above, lev_below,                 &
            lon_fract, lat_fract, temp_lon,           &
            lon_lat_lev(3)

real(r8), allocatable :: val_11(:), val_12(:), val_21(:), val_22(:), a(:, :)
integer,  allocatable :: vstatus(:), track_vstatus(:)

! FIXME: Positions within the rank 2 and 3 fields.   I don't remember the issue...
integer  :: s_type, s_type_01d,s_type_2d,s_type_3d,   &
            lon_ind_below, lon_ind_above, lat_ind_below, lat_ind_above, &
            num_lons
character(len=8)   :: lon_name, lat_name, lev_name
integer   :: ens_size

! FIXME; idea of right number of dimensions for each field...
! These are observations and will have 2d or 3d locations, but the
! corresponding state-vector component could be missing one of the dimensions.
! Surface pressure is the obvious example, but parameterization tuning might
! introduce others.
! Such artificial fields would not have observations associated with them.
! So assume that observed fields are not missing any dimensions.

ens_size = copies_in_window(state_ens_handle)
allocate(val_11(ens_size),val_12(ens_size), val_21(ens_size), val_22(ens_size))
allocate(a(ens_size, 2))
allocate(vstatus(ens_size))
allocate(track_vstatus(ens_size))


! Start with failure, then change as warranted.
istatus(:) = 1
vstatus(:) = MISSING_I
val_11 = MISSING_R8
val_12 = MISSING_R8
val_21 = MISSING_R8
val_22 = MISSING_R8

! Get the observation (horizontal) position, in degrees
lon_lat_lev = get_location(obs_loc)

! Check whether model_mod can interpolate the requested variable.
! Pressure (3d) can't be specified as a state vector field (so s_type will = MISSING_I),
! but can be calculated for CAM, so obs_kind = KIND_PRESSURE is acceptable.
! obs_kind truly is a DART KIND variable, generally passed from
! obs_def/obs_def_XXX.f90: call interpolate.
s_type = dart_to_cam_types(obs_kind)

if (s_type == MISSING_I .and. &
   (obs_kind /= KIND_PRESSURE) .and.  (obs_kind /= KIND_SURFACE_ELEVATION)) then
   ! HK these are not set because they are set at the beginning?
!   istatus = 1
!   interp_val = MISSING_R8
   write(string1,*) 'Wrong type of obs = ', obs_kind
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
   return
endif

! Get lon and lat dimension names.

! Set [lon,lat,lev] names to a default, which will be overwritten for variables
! in the state vector, but not for other acceptable variables (3D pressure, surface
! elevation, ...?)
lon_name = 'lon'
lat_name = 'lat'
if (obs_kind == KIND_SURFACE_ELEVATION) then
   lev_name = 'none'
elseif (obs_kind == KIND_PRESSURE) then
   lev_name = 'lev'
endif

! DART can't handle any 0d or 1d ob fields, so lump them together for elimination in this search.
s_type_01d = state_num_0d + state_num_1d
s_type_2d = s_type - s_type_01d
s_type_3d = s_type_2d - state_num_2d

if (s_type == MISSING_I .and. &
   (obs_kind == KIND_PRESSURE) .or.  (obs_kind == KIND_SURFACE_ELEVATION)) then
   ! use defaults lon_name and lat_name set above
elseif (s_type <= state_num_0d + state_num_1d) then
   ! error; can't deal with observed variables that are 0 or 1D in model_mod.
!   istatus = 1
!   interp_val = MISSING_R8
   write(string1,*) 'DART cannot handle 0d or 1d observations of ', cflds(s_type), &
        ' because DART requires a (lon,lat) location for each observation '
   write(string2,*) 'Skipping this observation'
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate, text2=string2)
   return
elseif (s_type_2d > 0 .and. s_type_2d <= state_num_2d) then
   lon_name = dim_names(s_dimid_2d(1,s_type_2d))
   lat_name = dim_names(s_dimid_2d(2,s_type_2d))
   lev_name = 'none'
elseif (s_type_3d > 0 .and. s_type_3d <= state_num_3d) then
   lon_name = dim_names(s_dimid_3d(2,s_type_3d))
   lat_name = dim_names(s_dimid_3d(3,s_type_3d))
   lev_name = dim_names(s_dimid_3d(1,s_type_3d))
else
!   istatus = 1
!   interp_val = MISSING_R8
   write(string1,*) 'Unexpected state type value, s_type = ', s_type
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
   return
endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! staggered longitudes; slon (4x5 fv grid) = [-2.5, 2.5,...,352.5]  !
!                        lon ( "         ) = [    0.,  5.,...,  355.]
! This is a complication for lon = 359, for example.  It's not in the range of slon.
!    coord_index handles it.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Compute bracketing lon indices
! Define a local longitude to deal with CAM-FV's staggered, longitude grid.
temp_lon = lon_lat_lev(1)

if (lon_name == 'lon') then
   num_lons  = lon%length
   bot_lon   = lon%vals(1)
   top_lon   = lon%vals(num_lons)
   delta_lon = lon%vals(2) - lon%vals(1)
elseif (lon_name == 'slon') then
   num_lons  = slon%length
   bot_lon   = slon%vals(1)
   top_lon   = slon%vals(num_lons)
   delta_lon = slon%vals(2) - slon%vals(1)
   ! Make certain longitudes conform to the CAM staggered grid.
   if ((lon_lat_lev(1) - top_lon) >= delta_lon) temp_lon = lon_lat_lev(1) - 360.0_r8
endif

if (temp_lon >= bot_lon .and. temp_lon   <  top_lon) then
   ! adding the 1 makes up for subtracting the bot_lon.
   lon_ind_below = int((temp_lon - bot_lon) / delta_lon) + 1
   lon_ind_above = lon_ind_below + 1
   lon_fract = (temp_lon - ((lon_ind_below - 1) * delta_lon + bot_lon)) / delta_lon
else
   ! At wraparound point
   lon_ind_above = 1
   lon_ind_below = num_lons
   lon_fract = (temp_lon - top_lon) / delta_lon
endif


! Next, compute neighboring lat rows
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE MAY BE WRONG
! Inefficient search used for latitudes in Gaussian grid. Might want to speed up.
! CAM-FV; lat = -90., ...   ,90.
!        slat =   -88.,...,88.

call coord_index(lat_name, lon_lat_lev(2), lat_ind_above, lat_ind_below)

! FIXME; maybe move this into coord_index
!        Probably not; coord_index sometimes returns the single closest index,
!                      which will always be the first index returned.
!                      I suppose there could be a flag argument telling coord_index
!                      whether to return 1 or a pair, with the 2nd index always > first
!                      (or vice versa).
!        calculate and return fraction too?
if (lat_ind_above == lat_ind_below) then
   if (lat_ind_above == 1) then
      lat_fract = 0.0_r8
   else                     !both must be equal to the max (s)lat index
      lat_fract = 1.0_r8
   endif
else
   if (lat_ind_above < lat_ind_below) then
      ! switch order
      i = lat_ind_above
      lat_ind_above = lat_ind_below
      lat_ind_below = i
   endif
   ! only lat_xxx is changed by these calls
   call coord_val(lat_name, lat_ind_below, lon_below, lat_below, lev_below)
   call coord_val(lat_name, lat_ind_above, lon_below, lat_above, lev_below)
   lat_fract = (lon_lat_lev(2) - lat_below) / (lat_above - lat_below)
endif

! Find the values for the four corners

! Determine the vertical coordinate: model level, pressure, or height
if (obs_kind == KIND_SURFACE_ELEVATION) then
   ! Acceptable field that's not in the state vector
   ! convert from geopotential height to real height in meters

   val_11 = phis(lon_ind_below, lat_ind_below) / gravity_const
   val_12 = phis(lon_ind_below, lat_ind_above) / gravity_const
   val_21 = phis(lon_ind_above, lat_ind_below) / gravity_const
   val_22 = phis(lon_ind_above, lat_ind_above) / gravity_const

elseif (vert_is_level(obs_loc)) then
   ! Pobs
   ! FIXME; I may want to change get_val_level to accept REAL level, not INT.
   !        What's the benefit?
   !        But it would be inconsistent with lon_ and lat_ indices,
   !           and I'd have to create an integer level anyway.
   !        May also want to handle staggered vertical grid (ilev).
   call get_val_level_distrib(state_ens_handle, lon_ind_below, lat_ind_below, nint(lon_lat_lev(3)), obs_kind, val_11, vstatus)
         track_vstatus = vstatus
   if (vstatus /= 1) &
      call get_val_level_distrib(state_ens_handle, lon_ind_below, lat_ind_above, nint(lon_lat_lev(3)), obs_kind, val_12, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   if (vstatus /= 1) &
      call get_val_level_distrib(state_ens_handle, lon_ind_above, lat_ind_below, nint(lon_lat_lev(3)), obs_kind, val_21, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   if (vstatus /= 1) &
      call get_val_level_distrib(state_ens_handle, lon_ind_above, lat_ind_above, nint(lon_lat_lev(3)), obs_kind, val_22, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo
   vstatus = track_vstatus
   ! Pobs end

elseif (vert_is_pressure(obs_loc)) then

   call get_val_pressure_distrib(val_11, state_ens_handle, lon_ind_below, lat_ind_below, lon_lat_lev(3), obs_type,vstatus)
   track_vstatus = vstatus

   call get_val_pressure_distrib(val_12, state_ens_handle, lon_ind_below, lat_ind_above, lon_lat_lev(3), obs_type,vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   call get_val_pressure_distrib(val_21, state_ens_handle, lon_ind_above, lat_ind_below, lon_lat_lev(3), obs_type,vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 ) track_vstatus(e) = vstatus(e)
   enddo

   call get_val_pressure_distrib(val_22, state_ens_handle, lon_ind_above, lat_ind_above, lon_lat_lev(3), obs_type,vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 ) track_vstatus(e) = vstatus(e)
   enddo
   vstatus = track_vstatus

elseif (vert_is_height(obs_loc)) then

   !HK model_heights is called 4 times on the obs_loc here. Also why call the interpolation
   ! routine over and over again?

   call get_val_height_distrib(val_11, state_ens_handle, lon_ind_below, lat_ind_below, lon_lat_lev(3), obs_type, vstatus)
   track_vstatus = vstatus

   call get_val_height_distrib(val_12, state_ens_handle, lon_ind_below, lat_ind_above, lon_lat_lev(3), obs_type, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   call get_val_height_distrib(val_21, state_ens_handle, lon_ind_above, lat_ind_below, lon_lat_lev(3), obs_type, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   call get_val_height_distrib(val_22, state_ens_handle, lon_ind_above, lat_ind_above, lon_lat_lev(3), obs_type, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo
   vstatus = track_vstatus

elseif (vert_is_surface(obs_loc)) then
   ! The 'lev' argument is set to 1 because there is no level for these types, and 'lev' will be
   ! ignored.
   call get_val_distrib(state_ens_handle, ens_size, lon_ind_below, lat_ind_below, 1, obs_kind, val_11, vstatus)
   track_vstatus = vstatus

   call get_val_distrib(state_ens_handle, ens_size, lon_ind_below, lat_ind_above, 1, obs_kind, val_12, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   call get_val_distrib(state_ens_handle, ens_size, lon_ind_above, lat_ind_below, 1, obs_kind, val_21, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   call get_val_distrib(state_ens_handle, ens_size, lon_ind_above, lat_ind_above, 1, obs_kind, val_22, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo
   vstatus = track_vstatus

! This needs to be at the end of the block.  Otherwise, it short circuits GPS
! which asks for pressures on heights.
! elseif (obs_kind == KIND_PRESSURE) then
!    ! Calculate pressures from surface pressures and A and B coeffs.
!    ! Pobs
!     write(string1,'(A)') 'No code available yet for obs_kind = KIND_PRESSURE '
!     call error_handler(E_ERR, 'interp_lon_lat', string1)

elseif (vert_is_scale_height(obs_loc)) then
   ! Need option for vert_is_scale_height
   write(string1,*)'Scale height is not an acceptable vert coord yet.  Skipping observation'
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
   return

! Need option for vert_is_undefined
else
   write(string1,*) '   No vert option chosen!'
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
   return

endif

! Conundrum (but unimportant for now): an ob could be excluded for > 1 reason.
! E.g. it's too far north and it's above the highest_obs_pressure_Pa.
! What istatus to return? a 2 (or more) digit number?  Like vstatus*10 + 4?
! HK loop around ensembles
do e = 1, ens_size
   ! lat is already converted to degrees by get_location
   if (abs(lon_lat_lev(2)) > max_obs_lat_degree .and. vstatus(e) /= 1) then
      ! Define istatus to be a combination of vstatus (= 0 or 2 (for higher than highest_obs...))
      ! and whether the ob is poleward of the limits set in the namelist (+ 4).
      istatus(e) = 10*vstatus + 4
   else
      istatus(e) = vstatus(e)
   end if

   ! indices of vals are (longitude, latitude)
   if (istatus(e) /= 1) then
      a(e, 1) = lon_fract * val_21(e) + (1.0_r8 - lon_fract) * val_11(e)
      a(e, 2) = lon_fract * val_22(e) + (1.0_r8 - lon_fract) * val_12(e)

      interp_val(e) = lat_fract * a(e, 2) + (1.0_r8 - lat_fract) * a(e, 1)

   else
      interp_val(e) = MISSING_R8
   end if

enddo

end subroutine interp_lonlat_distrib

!-----------------------------------------------------------------------
! HK single value (mean) version of model interpolate. Don't really want to do this.
! This is for model_heights called from convert_vert
! Now need mean versions of every call in this subroutine.
recursive subroutine interp_lonlat_distrib_mean(state_ens_handle, obs_loc, obs_kind, istatus, interp_val)

! Find the 4 corners of the lon-lat grid cell that encloses an ob at 'obs_loc'
! and interpolate the values of obs_kind to that location.

! istatus   meaning                  return expected obs?   assimilate?
! 0         obs and model are fine;  yes                    yes
! 1         fatal problem;           no                     no
! 2         exclude valid obs        yes                    no   (ob > ! highest_obs_X)
! 3         unfamiliar obs type      no                     no
! 4         ob excl by namelist(lat) yes                    no
! NM        2 digit number means more than one namelist reason to exclude from assim.

! Any value > 0 will not be assimilated (---> QC non-0).
! Do we want some istatus values to tell filter to evaluate (---> QC of 1)?
! That would be nice, but filter has no convention for understanding non-0
! values from model_mod (from all the available models).  So all non-0 values of
! istatus ---> QC = 4.


real(r8),            intent(in) :: st_vec(:)
type(location_type), intent(in) :: obs_loc
integer,             intent(in) :: obs_kind
real(r8),           intent(out) :: interp_val
integer,            intent(out) :: istatus

! FIXME; In future DARTs it may be useful to return the DART KIND too.
!        also convert to a field name (DART subroutine (get_raw_...?)).

integer  :: i, vstatus
real(r8) :: bot_lon, top_lon, delta_lon,                                &
            lon_below, lat_below, lat_above, lev_below,                 &
            lon_fract, lat_fract, vals(2, 2), temp_lon, a(2),           &
            lon_lat_lev(3)

! FIXME: Positions within the rank 2 and 3 fields.   I don't remember the issue...
integer  :: s_type, s_type_01d,s_type_2d,s_type_3d,   &
            lon_ind_below, lon_ind_above, lat_ind_below, lat_ind_above, &
            num_lons
character(len=8)   :: lon_name, lat_name, lev_name

! FIXME; idea of right number of dimensions for each field...
! These are observations and will have 2d or 3d locations, but the
! corresponding state-vector component could be missing one of the dimensions.
! Surface pressure is the obvious example, but parameterization tuning might
! introduce others.
! Such artificial fields would not have observations associated with them.
! So assume that observed fields are not missing any dimensions.

! Start with failure, then change as warranted.
istatus    = 1
vstatus    = MISSING_I
vals       = MISSING_R8
interp_val = MISSING_R8

! Get the observation (horizontal) position, in degrees
lon_lat_lev = get_location(obs_loc)

! Check whether model_mod can interpolate the requested variable.
! Pressure (3d) can't be specified as a state vector field (so s_type will = MISSING_I),
! but can be calculated for CAM, so obs_kind = KIND_PRESSURE is acceptable.
! obs_kind truly is a DART KIND variable, generally passed from
! obs_def/obs_def_XXX.f90: call interpolate.
s_type = dart_to_cam_types(obs_kind)

if (s_type == MISSING_I .and. &
   (obs_kind /= KIND_PRESSURE) .and.  (obs_kind /= KIND_SURFACE_ELEVATION)) then
!   istatus = 1
!   interp_val = MISSING_R8
   write(string1,*) 'Wrong type of obs = ', obs_kind
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
   return
endif

! Get lon and lat dimension names.

! Set [lon,lat,lev] names to a default, which will be overwritten for variables
! in the state vector, but not for other acceptable variables (3D pressure, surface
! elevation, ...?)
lon_name = 'lon'
lat_name = 'lat'
if (obs_kind == KIND_SURFACE_ELEVATION) then
   lev_name = 'none'
elseif (obs_kind == KIND_PRESSURE) then
   lev_name = 'lev'
endif

! DART can't handle any 0d or 1d ob fields, so lump them together for elimination in this search.
s_type_01d = state_num_0d + state_num_1d
s_type_2d = s_type - s_type_01d
s_type_3d = s_type_2d - state_num_2d

if (s_type == MISSING_I .and. &
   (obs_kind == KIND_PRESSURE) .or.  (obs_kind == KIND_SURFACE_ELEVATION)) then
   ! use defaults lon_name and lat_name set above
elseif (s_type <= state_num_0d + state_num_1d) then
   ! error; can't deal with observed variables that are 0 or 1D in model_mod.
!   istatus = 1
!   interp_val = MISSING_R8
   write(string1,*) 'DART cannot handle 0d or 1d observations of ', cflds(s_type), &
        ' because DART requires a (lon,lat) location for each observation '
   write(string2,*) 'Skipping this observation'
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate, text2=string2)
   return
elseif (s_type_2d > 0 .and. s_type_2d <= state_num_2d) then
   lon_name = dim_names(s_dimid_2d(1,s_type_2d))
   lat_name = dim_names(s_dimid_2d(2,s_type_2d))
   lev_name = 'none'
elseif (s_type_3d > 0 .and. s_type_3d <= state_num_3d) then
   lon_name = dim_names(s_dimid_3d(2,s_type_3d))
   lat_name = dim_names(s_dimid_3d(3,s_type_3d))
   lev_name = dim_names(s_dimid_3d(1,s_type_3d))
else
!   istatus = 1
!   interp_val = MISSING_R8
   write(string1,*) 'Unexpected state type value, s_type = ', s_type
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
   return
endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! staggered longitudes; slon (4x5 fv grid) = [-2.5, 2.5,...,352.5]  !
!                        lon ( "         ) = [    0.,  5.,...,  355.]
! This is a complication for lon = 359, for example.  It's not in the range of slon.
!    coord_index handles it.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Compute bracketing lon indices
! Define a local longitude to deal with CAM-FV's staggered, longitude grid.
temp_lon = lon_lat_lev(1)

if (lon_name == 'lon') then
   num_lons  = lon%length
   bot_lon   = lon%vals(1)
   top_lon   = lon%vals(num_lons)
   delta_lon = lon%vals(2) - lon%vals(1)
elseif (lon_name == 'slon') then
   num_lons  = slon%length
   bot_lon   = slon%vals(1)
   top_lon   = slon%vals(num_lons)
   delta_lon = slon%vals(2) - slon%vals(1)
   ! Make certain longitudes conform to the CAM staggered grid.
   if ((lon_lat_lev(1) - top_lon) >= delta_lon) temp_lon = lon_lat_lev(1) - 360.0_r8
endif

if (temp_lon >= bot_lon .and. temp_lon   <  top_lon) then
   ! adding the 1 makes up for subtracting the bot_lon.
   lon_ind_below = int((temp_lon - bot_lon) / delta_lon) + 1
   lon_ind_above = lon_ind_below + 1
   lon_fract = (temp_lon - ((lon_ind_below - 1) * delta_lon + bot_lon)) / delta_lon
else
   ! At wraparound point
   lon_ind_above = 1
   lon_ind_below = num_lons
   lon_fract = (temp_lon - top_lon) / delta_lon
endif


! Next, compute neighboring lat rows
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE MAY BE WRONG
! Inefficient search used for latitudes in Gaussian grid. Might want to speed up.
! CAM-FV; lat = -90., ...   ,90.
!        slat =   -88.,...,88.

call coord_index(lat_name, lon_lat_lev(2), lat_ind_above, lat_ind_below)

! FIXME; maybe move this into coord_index
!        Probably not; coord_index sometimes returns the single closest index,
!                      which will always be the first index returned.
!                      I suppose there could be a flag argument telling coord_index
!                      whether to return 1 or a pair, with the 2nd index always > first
!                      (or vice versa).
!        calculate and return fraction too?
if (lat_ind_above == lat_ind_below) then
   if (lat_ind_above == 1) then
      lat_fract = 0.0_r8
   else                     !both must be equal to the max (s)lat index
      lat_fract = 1.0_r8
   endif
else
   if (lat_ind_above < lat_ind_below) then
      ! switch order
      i = lat_ind_above
      lat_ind_above = lat_ind_below
      lat_ind_below = i
   endif
   ! only lat_xxx is changed by these calls
   call coord_val(lat_name, lat_ind_below, lon_below, lat_below, lev_below)
   call coord_val(lat_name, lat_ind_above, lon_below, lat_above, lev_below)
   lat_fract = (lon_lat_lev(2) - lat_below) / (lat_above - lat_below)
endif

! Find the values for the four corners

! Determine the vertical coordinate: model level, pressure, or height
if (obs_kind == KIND_SURFACE_ELEVATION) then
   ! Acceptable field that's not in the state vector
   ! convert from geopotential height to real height in meters
   vals(1,1) = phis(lon_ind_below, lat_ind_below) / gravity_const
   vals(1,2) = phis(lon_ind_below, lat_ind_above) / gravity_const
   vals(2,1) = phis(lon_ind_above, lat_ind_below) / gravity_const
   vals(2,2) = phis(lon_ind_above, lat_ind_above) / gravity_const

elseif (vert_is_level(obs_loc)) then
   ! Pobs
   ! FIXME; I may want to change get_val_level to accept REAL level, not INT.
   !        What's the benefit?
   !        But it would be inconsistent with lon_ and lat_ indices,
   !           and I'd have to create an integer level anyway.
   !        May also want to handle staggered vertical grid (ilev).
      call get_val_level(st_vec, lon_ind_below, lat_ind_below, nint(lon_lat_lev(3)), obs_kind, vals(1,1), vstatus)
   if (vstatus /= 1) &
      call get_val_level(st_vec, lon_ind_below, lat_ind_above, nint(lon_lat_lev(3)), obs_kind, vals(1,2), vstatus)
   if (vstatus /= 1) &
      call get_val_level(st_vec, lon_ind_above, lat_ind_below, nint(lon_lat_lev(3)), obs_kind, vals(2,1), vstatus)
   if (vstatus /= 1) &
      call get_val_level(st_vec, lon_ind_above, lat_ind_above, nint(lon_lat_lev(3)), obs_kind, vals(2,2), vstatus)
   ! Pobs end

elseif (vert_is_pressure(obs_loc)) then
      call get_val_pressure(st_vec,lon_ind_below,lat_ind_below,lon_lat_lev(3),obs_kind,vals(1,1),vstatus)
   if (vstatus /= 1) &
      call get_val_pressure(st_vec,lon_ind_below,lat_ind_above,lon_lat_lev(3),obs_kind,vals(1,2),vstatus)
   if (vstatus /= 1) &
      call get_val_pressure(st_vec,lon_ind_above,lat_ind_below,lon_lat_lev(3),obs_kind,vals(2,1),vstatus)
   if (vstatus /= 1) &
      call get_val_pressure(st_vec,lon_ind_above,lat_ind_above,lon_lat_lev(3),obs_kind,vals(2,2),vstatus)

elseif (vert_is_height(obs_loc)) then
      call get_val_height(st_vec, lon_ind_below, lat_ind_below, lon_lat_lev(3), obs_loc, obs_kind, vals(1,1), vstatus)
   if (vstatus /= 1) &
      call get_val_height(st_vec, lon_ind_below, lat_ind_above, lon_lat_lev(3), obs_loc, obs_kind, vals(1,2), vstatus)
   if (vstatus /= 1) &
      call get_val_height(st_vec, lon_ind_above, lat_ind_below, lon_lat_lev(3), obs_loc, obs_kind, vals(2,1), vstatus)
   if (vstatus /= 1) &
      call get_val_height(st_vec, lon_ind_above, lat_ind_above, lon_lat_lev(3), obs_loc, obs_kind, vals(2,2), vstatus)


elseif (vert_is_surface(obs_loc)) then
   ! The 'lev' argument is set to 1 because there is no level for these types, and 'lev' will be
   ! ignored.
                     call get_val(st_vec, lon_ind_below, lat_ind_below, 1, obs_kind, vals(1,1), vstatus)
   if (vstatus /= 1) call get_val(st_vec, lon_ind_below, lat_ind_above, 1, obs_kind, vals(1,2), vstatus)
   if (vstatus /= 1) call get_val(st_vec, lon_ind_above, lat_ind_below, 1, obs_kind, vals(2,1), vstatus)
   if (vstatus /= 1) call get_val(st_vec, lon_ind_above, lat_ind_above, 1, obs_kind, vals(2,2), vstatus)


! This needs to be at the end of the block.  Otherwise, it short circuits GPS
! which asks for pressures on heights.
! elseif (obs_kind == KIND_PRESSURE) then
!    ! Calculate pressures from surface pressures and A and B coeffs.
!    ! Pobs
!     write(string1,'(A)') 'No code available yet for obs_kind = KIND_PRESSURE '
!     call error_handler(E_ERR, 'interp_lon_lat', string1)

elseif (vert_is_scale_height(obs_loc)) then
   ! Need option for vert_is_scale_height
   write(string1,*)'Scale height is not an acceptable vert coord yet.  Skipping observation'
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
   return

! Need option for vert_is_undefined
else
   write(string1,*) '   No vert option chosen!'
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
   return

endif

! Conundrum (but unimportant for now): an ob could be excluded for > 1 reason.
! E.g. it's too far north and it's above the highest_obs_pressure_Pa.
! What istatus to return? a 2 (or more) digit number?  Like vstatus*10 + 4?
if (vstatus == 1) then
   return     ! Failed to get value for interpolation; return istatus = 1
              ! Error will be handled by calling routines?
else
   if (abs(lon_lat_lev(2)) > max_obs_lat_degree) then
      ! Define istatus to be a combination of vstatus (= 0 or 2 (for higher than highest_obs...))
      ! and whether the ob is poleward of the limits set in the namelist (+ 4).
      istatus = 10*vstatus + 4
   else
      istatus = vstatus
   endif
endif

! indices of vals are (longitude, latitude)
do i = 1, 2
   a(i) = lon_fract * vals(2, i) + (1.0_r8 - lon_fract) * vals(1, i)
enddo
interp_val = lat_fract * a(2) + (1.0_r8 - lat_fract) * a(1)

end subroutine interp_lonlat


!-----------------------------------------------------------------------

! Pobs
subroutine get_val_level_distrib(state_ens_handle, lon_index, lat_index, level, obs_kind, val, istatus)

! Gets the value on level for variable obs_kind
! at lon_index, lat_index horizontal grid point
!
! written by Kevin Raeder, based on code from Hui Liu 4/28/2006 and get_val_pressure
!         from Jeff Anderson
!
! This routine indicates things with the return code:
!   istatus 0 - success
!   istatus 1 - failure (e.g. above or below highest/lowest level, or can't
!                          interpolate the value)
!   istatus 2 - val is set successfully, but level is above highest_obs_level
!
! This routine assumes level is an integer value.  To make it work for
! fractional levels (some models do support this) the code would have to be
! altered to find the value at the level below and above, and interpolate in
! the vertical.

type(ensemble_type), intent(in) :: state_ens_handle
integer,  intent(in)  :: lon_index
integer,  intent(in)  :: lat_index
integer,  intent(in)  :: level
integer,  intent(in)  :: obs_kind
real(r8), intent(out) :: val(:)
integer,  intent(out) :: istatus(:)

integer               :: vstatus, num_levs, i, lowest_ok
real(r8)              :: p_surf, threshold

! Start with failure condition
istatus = 1
vstatus = 1
val = MISSING_R8

ens_size = copies_in_window(state_ens_handle)

! This assumes that all variables are defined on model levels, not on interface levels.
num_levs = dim_sizes(find_name('lev',dim_names))

! Exclude obs below the model's lowest level and above the highest level
if (level > num_levs .or. level < 1) return

! Interpolate in vertical to get two bounding levels, but treat pressure
! specially since it has to be computed from PS instead of interpolated.

if (obs_kind == KIND_PRESSURE) then

   ! p_surf is returned in pascals, which is the right units for plevs_cam() below.
   call get_val_distrib(state_ens_handle, ens_size, lon_index, lat_index, -1, KIND_SURFACE_PRESSURE, p_surf, vstatus)

   do e = 1, ens_size
      if (vstatus == 0) return
         ! Next, get the values on the levels for this PS.
         call plevs_cam(p_surf(e), num_levs, p_col)
         val = p_col(level)
   enddo
else
   call get_val_distrib(state_ens_handle, ens_size, lon_index, lat_index, level, obs_kind, val, vstatus)
endif

! if this routine is called with a location that has a vertical level above
! the pressure cutoff, go ahead and compute the value but return an istatus=2
! (unless some other error occurs later in this routine).  note that smaller
! level numbers are higher up in the atmosphere; level 1 is at the top.

if (level < highest_obs_level) then
   istatus = 2
else
   istatus = vstatus
endif

end subroutine get_val_level_distrib
! Pobs end

!-----------------------------------------------------------------------

subroutine get_val_pressure_distrib(state_ens_handle, lon_index, lat_index, pressure, obs_kind, val, istatus)

! Gets the vertically interpolated value on pressure for variable obs_kind
! at lon_index, lat_index horizontal grid point
!
! This routine indicates things with the return code:
!   istatus 0 - success
!   istatus 1 - failure (e.g. above or below highest/lowest level, or can't
!                          interpolate the value)
!   istatus 2 - val is set successfully, but vert is above highest_obs_pressure
!
! Excludes observations below lowest level pressure and above highest level pressure.

type(ensemble_type), intent(in) :: state_ens_handle
real(r8), intent(in)  :: pressure
integer,  intent(in)  :: lon_index
integer,  intent(in)  :: lat_index
integer,  intent(in)  :: obs_kind
real(r8), intent(out) :: val(:)
integer,  intent(out) :: istatus(:)

real(r8), allocatable :: bot_val(:), top_val(:), p_surf(:), frac(:)
integer,  allocatable :: bot_lev(:), top_lev(:)
real(r8), allocatable :: ps_local(:, :)
real(r8), allocatable :: p_col_distrib(:, :)
integer               :: i, vstatus, num_levs ! HK vstaus?
integer               :: fld_index
integer               :: ens_size, e

! No errors to start with
istatus = 1
vstatus = 1
val     = MISSING_R8

ens_size = copies_in_window(state_ens_handle)
allocate(bot_val(ens_size), top_val(ens_size), p_surf(ens_size), frac(ens_size))
allocate(ps_local(2, ens_size))
allocate(p_col_distrib(ens_size, num_levs))
allocate(bot_lev(ens_size), top_lev(ens_size)) !> @todo HK I don't know why you need two values, one is just + 1 to the other

! Need to get the surface pressure at this point.
! Find out whether the observed field is a staggered field in CAM.
! Check whether the state vector has wind components on staggered grids, i.e. whether CAM is FV.
! find_name returns 0 if the field name is not found in the cflds list.

fld_index   = find_name('PS',cflds)
i = index_from_grid(1,lon_index,lat_index,  fld_index)
!ps_local(1, :) = st_vec(i)
call get_state(ps_local(1, :), i, state_ens_handle)

if (obs_kind == KIND_U_WIND_COMPONENT .and. find_name('US', cflds) /= 0) then
   ! ps defined on lat grid (-90...90, nlat = nslat + 1),
   !    need it on staggered lat grid, which starts half a grid spacing north.

   i = index_from_grid(1,lon_index,lat_index+1,fld_index)
   !ps_local(2) = st_vec(i)
   call get_state(ps_local(2, :), i, state_ens_handle)
   p_surf(:) = (ps_local(1, :) + ps_local(2, :))* 0.5_r8

elseif (obs_kind == KIND_V_WIND_COMPONENT .and. find_name('VS', cflds) /= 0) then
   ! lon =     0...     255 (for 5 degree grid)
   !slon = -2.5 ... 252.5
   if (lon_index == slon%length) then
      i = index_from_grid(1,1,          lat_index ,fld_index)
   else
      i = index_from_grid(1,lon_index+1,lat_index ,fld_index)
   endif
   !ps_local(2) = st_vec(i)
   call get_state(ps_local(2, :), i, state_ens_handle)
   p_surf      = (ps_local(1, :) + ps_local(2, :))* 0.5_r8
else
   ! A-grid ps can be retrieved from state vector, which was used to define ps on entry to
   ! model_interpolate.
   p_surf     = ps_local(1, :)
endif

! Next, get the pressures on the levels for this ps
! Assuming we'll only need pressures on model mid-point levels, not interface levels.
! This pressure column will be for the correct grid for obs_kind, since p_surf was taken
!     from the grid-correct ps[_stagr] grid
num_levs = dim_sizes(find_name('lev',dim_names))
call plevs_cam_distrib(p_surf, num_levs, p_col_distrib, ens_size)

! Exclude obs below the model's lowest level and above the highest level
! We *could* possibly use ps and p(num_levs) to interpolate for points below the lowest level.
do e = 1, ens_size
   if (pressure <= p_col_distrib(1,e) .or. pressure >= p_col_distrib(num_levs,e)) then
      istatus(e) = 1
      val(e) = MISSING_R8
      !return
   endif
enddo

! Interpolate in vertical to get two bounding levels

! Search down through pressures
levloop: do i = 2, num_levs
   if (pressure < p_col(i, e)) then
      top_lev = i -1
      bot_lev = i
      frac = (p_col_distrib(i, e) - pressure) / &
             (p_col_distrib(i, e) - p_col_distirb(i - 1, e))
      exit levloop
   endif
enddo levloop

! Pobs
if (obs_kind == KIND_PRESSURE) then
   ! can't get pressure on levels from state vector; get it from p_col instead
   bot_val = p_col_distrib(bot_lev, :)
   top_val = p_col_distrib(top_lev, :)
else

   ! need to grab values for each bot_val
   do e = 1, ens_size ! HK you only need to do this for distinct bot_vals
      call get_val_distirb(state_ens_handle, lon_index, lat_index, bot_lev(e), obs_kind, bot_val, vstatus)
      if (vstatus /= 1) call get_val_distrib(state_ens_handle, lon_index, lat_index, top_lev(e), obs_kind, top_val, vstatus)
      ! Failed to get value for interpolation; return istatus = 1
      !if (vstatus == 1)
      istatus(e) = vstatus
   enddo
endif
! Pobs

! if this routine is called with a location that has a vertical pressure above
! the pressure cutoff, go ahead and compute the value but return an istatus=2
! (unless some other error occurs later in this routine).

if (pressure < highest_obs_pressure_Pa) then
   istatus(:) = 2
!else
   ! HK can't do this with an ensemble
!   istatus = 0
endif

do e = 1, ens_size
   if(istatus(e) == 0 .or. istatus(e) == 2) val = (1.0_r8 - frac) * bot_val + frac * top_val
enddo

deallocate(bot_val, top_val, p_surf, frac, p_col_distrib)
deallocate(bot_lev, top_lev)

end subroutine get_val_pressure_distrib

!-----------------------------------------------------------------------

subroutine get_val_height_distrib(st_vec, lon_index, lat_index, height, location, obs_kind, val, istatus)

! Gets the vertically interpolated value on height for variable obs_kind
! at lon_index, lat_index horizontal grid point
!
! written by Kevin Raeder, based on code from Hui Liu 4/28/2006 and get_val_pressure
!         from Jeff Anderson
!
! This routine indicates things with the return code:
!   istatus 0 - success
!   istatus 1 - failure (e.g. above or below highest/lowest level, or can't
!                          interpolate the value)
!   istatus other - val is set successfully, but obs is excluded according to namelist restrictions.

real(r8),            intent(in)  :: st_vec(:)
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
real(r8),            intent(in)  :: height
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
real(r8),            intent(out) :: val
integer,             intent(out) :: istatus

integer  :: top_lev, bot_lev, i, vstatus, num_levs, fld_index, ind
real(r8) :: bot_val, top_val, frac
real(r8) :: p_surf, ps_local(2)
logical  :: stagr_lon, stagr_lat

! No errors to start with
istatus   = 1
vstatus   = 1
val       = MISSING_R8
stagr_lon = .false.
stagr_lat = .false.

! Assuming we'll only need pressures on model mid-point levels, not interface levels.
num_levs = dim_sizes(find_name('lev',dim_names))

! Need to get the surface pressure at this point.
! Check whether the state vector has wind components on staggered grids, i.e. whether CAM is FV.
! See get_val_pressure for more documentation.
fld_index   = find_name('PS',cflds)
ind         = index_from_grid(1,lon_index,lat_index,  fld_index)
ps_local(1) = st_vec(ind)

! find_name returns 0 if the field name is not found in the cflds list.
if (obs_kind == KIND_U_WIND_COMPONENT .and. find_name('US', cflds) /= 0) then
   stagr_lat = .true.
   ind = index_from_grid(1,lon_index,lat_index+1,fld_index)
   ps_local(2) = st_vec(ind)
   p_surf = (ps_local(1) + ps_local(2))* 0.5_r8
elseif (obs_kind == KIND_V_WIND_COMPONENT .and. find_name('VS', cflds) /= 0) then
   stagr_lon = .true.
   if (lon_index == slon%length) then
      ind = index_from_grid(1,1,          lat_index ,fld_index)
   else
      ind = index_from_grid(1,lon_index+1,lat_index ,fld_index)
   endif
   ps_local(2) = st_vec(ind)
   p_surf = (ps_local(1) + ps_local(2))* 0.5_r8
else
   p_surf = ps_local(1)
endif

! Next, get the heights on the levels for this ps

! We want to use the new vec for each new ob on height because the state was updated
! for all previous obs, and we want to use the most up to date state to get the best location.
! HK ******** THE STATE IS NOT UPDATED YOU ARE USING THE MEAN COPY from before assimilation. *******
call model_heights(num_levs, st_vec, p_surf, location, model_h, vstatus)
if (vstatus == 1) return    ! Failed to get model heights; return istatus = 1

! Exclude obs below the model's lowest level and above the highest level
if (height >= model_h(1) .or. height <= model_h(num_levs)) return

! ? Implement 3Dp here?  or should/can it not use the ens mean PS field?
call plevs_cam(p_surf, num_levs, p_col)

! The highest_obs_pressure_Pa has already been checked to ensure it's a valid value.
! So this loop will always set the highest_obs_height_m before exiting.
! This could be refined to interpolate between the p_col to highest_obs_pressure_Pa.
! Also, if using the nearest lower model level is good enough, then it might be good
! enough to only calculate highest_obs_height_m once; put if (highest_obs_height_m == MISSING_R8)
! around the loop.
! Actually, I see in gph2gmh that the heights in model_h are relative to mean sea level,
! so they will be independent from the surface height and vertical coordinate system.
! They will vary slightly with surface pressure.
! So I think that highest_obs_height_m could be calculated once
if (highest_obs_height_m == MISSING_R8) then
levloop: do i=2,num_levs
   if (p_col(i) > highest_obs_pressure_Pa) then
      ! highest_obs_height_m = model_h(i)
      highest_obs_height_m = model_h(i) + (model_h(i)-model_h(i-1))*  &
                                          ((p_col(i)-highest_obs_pressure_Pa) / &
                                           (p_col(i)-p_col(i-1)))
      exit levloop
   endif
enddo levloop
endif



! Interpolate in vertical to get two bounding levels.
! Search down through heights and set the enclosing level numbers
! and the fraction between them.  There has already been a test to
! ensure the height is between the levels (and has discarded values
! exactly equal to the limits), so this will always succeed.
lev2loop: do i = 2, num_levs
   if (height > model_h(i)) then
      top_lev = i -1
      bot_lev = i
      frac = (model_h(i) - height      ) / &
             (model_h(i) - model_h(i-1))
      exit lev2loop
   endif
enddo lev2loop


if (obs_kind == KIND_PRESSURE) then
   bot_val = p_col(bot_lev)
   top_val = p_col(top_lev)
else
                     call get_val(st_vec, lon_index, lat_index, bot_lev, obs_kind, bot_val, vstatus)
   if (vstatus == 0) call get_val(st_vec, lon_index, lat_index, top_lev, obs_kind, top_val, vstatus)
   ! Failed to get a value to use in interpolation
   if (vstatus == 1) return   
end if

if (height > highest_obs_height_m ) then
   ! if this routine is called with a location that has a vertical height above
   ! the pressure cutoff, go ahead and compute the value but return an istatus=2
   ! (unless some other error occurs later in this routine).
   istatus = 2
else
   istatus = 0
endif

val = (1.0_r8 - frac) * bot_val + frac * top_val

end subroutine get_val_height_distrib

!-----------------------------------------------------------------------

subroutine get_val_distrib(state_ens_handle, ens_size, lon_index, lat_index, level, obs_kind, val, istatus)

integer,             intent(in)  :: ens_size !< how may pieces of state to grab
real(r8),            intent(out) :: val(ens_size)
type(ensemble_type), intent(in)  :: state_ens_handle
integer, intent(in)   :: lon_index
integer, intent(in)   :: lat_index
integer, intent(in)   :: level
integer, intent(in)   :: obs_kind
integer, intent(out)  :: istatus

integer :: indx, field_type

! Start with error condition.
istatus = 1
val = MISSING_R8

field_type = dart_to_cam_types(obs_kind)
if (field_type <= 0 .or. field_type > nflds) return

indx = index_from_grid(level, lon_index, lat_index, field_type)
if (indx > 0 .and. indx <= model_size) then
   istatus = 0
   !val = st_vec(indx)
   call get_state(val, indx, state_ens_handle)
endif

end subroutine get_val_distrib

!-----------------------------------------------------------------------

subroutine set_highest_obs_limit()

! Verify that the value for highest_obs_pressure_Pa in the namelist is ok.
!
! If this routine detects an error it calls the error handler with a
! fatal error.  If it returns, the namelist value is ok.
!
! Sets the module global variable 'highest_obs_level', and references
! the hybm array.

integer  :: num_levs, i, lowest_ok
real(r8) :: p_surf

! This assumes that all variables are defined on model levels, not on interface levels.
num_levs = dim_sizes(find_name('lev',dim_names))

! This code determines the model level that is below but closest to the
! 'highest obs pressure' limit that was set in the namelist.  It is counting
! on the limit being set high enough that the level heights are
! determined solely by the pressure values with no contribution from the terrain.
! Instead of computing a surface pressure from an ensemble member at a particular
! longitude and latitude, assume a surface pressure of 1000 mb and compute
! a pressure column based on that.   Then, verify that the 'hybm' value at
! the selected level is 0 - otherwise the levels still have a contribution
! based on terrain and cannot be solely determined based on pressure.
! Then we can use this single level value for any lat/lon/ensemble member.

! Compute a pressure column based on a 1000 mb (*100 = pascals) surface pressure
call plevs_cam(P0%vals(1), num_levs, p_col)

! Loop downwards through pressure values (1 = model top, num_levs = bottom).
! The level is set to the highest level which is below the given threshold.


! highest_obs_level = 1
! do while ((p_col(highest_obs_level)) < highest_obs_pressure_Pa .and. highest_obs_level < num_levs)
!    highest_obs_level = highest_obs_level + 1
! enddo
High: do highest_obs_level=1,num_levs
   if (p_col(highest_obs_level) > highest_obs_pressure_Pa) exit High
enddo High

! Test to be sure user hasn't set level so low that contributions from
! terrain are going to be an issue.  If they do, tell them in the error
! message what the valid limit is.

if (hybm%vals(highest_obs_level) > 0.0_r8) then
   lowest_ok = 1
   findzero: do i=2, num_levs
      if (hybm%vals(i) > 0.0_r8) then
         lowest_ok = i-1
         exit findzero
      endif
   enddo findzero
   write(string1, '(A)') 'invalid value for namelist "highest_obs_pressure_Pa"'
   write(string2, '(A)') 'value is too large (and so located too low in atmosphere)'
   write(string3, '(A,F9.3,A)') 'must specify a value located above pressure ', p_col(lowest_ok), ' Pascals'
   call error_handler(E_ERR, 'set_highest_obs_limit', string1, source, revision, revdate, &
                      text2=string2, text3=string3)
endif

end subroutine set_highest_obs_limit

! End of model_interpolate section

!#######################################################################

! Vector-field translations

!-----------------------------------------------------------------------

subroutine prog_var_to_vector(var, st_vec)

type(model_type), intent(in)  :: var
real(r8),         intent(out) :: st_vec(:)

integer :: i, j, k, nf, indx

! Load components of state vector, starting with scalars (0D) and finishing with 3D
! A whole field will be loaded (by columns for 3D) before the next field is started.
! This is completely different than the B-grid organization, which loaded all the fields
! at a point before moving on to the next point.  The motivations for this change are:
! 1) This easily allows fields with the same rank, but different sizes to be loaded into
!    the vector (i.e. U_staggered  and T in the cam-fv)
! 2) The dominant form of access into the state vector is vertical interpolations in
!    get_expected_val and computation of columns of virtual temperature from T and Q
!    in model_heights.  model_get_close_states, which searched for all variables close
!    to an obs, is not part of the MPI DART, so spatially co-located variables don't
!    need to be close to each other in memory.

if (.not. module_initialized) call static_init_model()

indx = 0

!  0d variables
do nf = 1, state_num_0d
   indx = indx + 1
   st_vec(indx) = var%vars_0d(nf)
enddo

!  1d variables
do nf = 1, state_num_1d
   do i=1,s_dim_1d(nf)
      indx = indx + 1
      st_vec(indx) = var%vars_1d(i, nf)
   enddo
enddo

!  2d variables
do nf = 1, state_num_2d
   do j=1,s_dim_2d(2,nf)
   do i=1,s_dim_2d(1,nf)
      indx = indx + 1
      st_vec(indx) = var%vars_2d(i, j, nf)
   enddo
   enddo
enddo

!  3D fields, loaded by columns (note the coordinate order).
!  This section is only entered for models with logically rectangular grids,
!  which will have dimensions level, longitude, and latitude.
!  This is also looping over latitude values in the outer loop,
!  so if there is some spatial searching; all the pole points will be less scattered
!  than longitude in memory.
do nf= 1, state_num_3d
   do i=1,s_dim_3d(3,nf)   !lats
   do j=1,s_dim_3d(2,nf)   !lons
   do k=1,s_dim_3d(1,nf)   !levs  both reads and writes will be contiguous in this case
      indx = indx + 1
      st_vec(indx) = var%vars_3d(k,j,i, nf)
   enddo
   enddo
   enddo
enddo

if (indx /= model_size) then
   write(string1, *) 'Number of elements copied = ',indx,', must equal model_size, ',model_size
   call error_handler(E_ERR, 'prog_var_to_vector', string1, source, revision, revdate)
endif

end subroutine prog_var_to_vector

!-----------------------------------------------------------------------

subroutine vector_to_prog_var(st_vec, var)

real(r8),         intent(in)    :: st_vec(:)
type(model_type), intent(inout) :: var

integer :: i, j, k, nf, indx

if (.not. module_initialized) call static_init_model()

indx = 0

! 0d arrays
do nf = 1, state_num_0d
   indx = indx + 1
   var%vars_0d(nf) = st_vec(indx)
enddo

!  1d fields
do nf = 1, state_num_1d
   do i=1,s_dim_1d(nf)
      indx = indx + 1
      var%vars_1d(i, nf) = st_vec(indx)
   enddo
enddo

!  2d fields
do nf = 1, state_num_2d
   do j = 1, s_dim_2d(2,nf)
   do i = 1, s_dim_2d(1,nf)
      indx = indx + 1
      var%vars_2d(i, j, nf) = st_vec(indx)
   enddo
   enddo
enddo

! 3D fields; see comments in prog_var_to_vect
do nf = 1, state_num_3d
   do i = 1, s_dim_3d(3,nf)
   do j = 1, s_dim_3d(2,nf)
   do k = 1, s_dim_3d(1,nf)
      indx = indx + 1
      var%vars_3d(k,j,i, nf) = st_vec(indx)
   enddo
   enddo
   enddo
enddo

if (indx /= model_size) then
   write(string1, *) 'Number of elements copied = ',indx,', must equal model_size, ',model_size
   call error_handler(E_ERR, 'vector_to_prog_var', string1, source, revision, revdate)
endif

end subroutine vector_to_prog_var

! End of Vector-field translations


!#######################################################################
! get_close_obs section

!-----------------------------------------------------------------------
!>
!> Subroutine get_close_obs
!> 
!> get_close_obs takes as input an "observation" location, a DART TYPE (not KIND),
!> and a list of all potentially close locations and KINDS on this task.
!>
!> get_close_obs
!>    *) converts vertical coordinates as needed to vert_coord,
!>    *) calls location_mod/threed_sphere:get_close_obs,
!>       to which it sends this (converted) array of locations,
!>    *) gets back the distances and indices of those locations that are
!>       "close" to the base observation.
!>    *) tests for being above the highest_obs_pressure_Pa threshold,
!>       and increases the vertical distance based on height above highest_*.
!> 
!> @param[in]    filt_gc
!> The DART get_close_type containing the state variables which are potentially close to 'location'
!> 
!> @param[in]    base_obs_loc
!> The DART location_type location of the observation, which is the target of *get_close_obs*
!> 
!> @param[in]    base_obs_type 
!> The DART TYPE (not KIND) of the observation
!> 
!> @param[inout] locs(:)
!> The DART location_type locations of the potentially close state variables
!> 
!> @param[in]    kinds(:)
!> The DART KINDs of the potentially close state variables
!> 
!> @param[out]   num_close
!> The number of state variables which are deemed to be close to the observation
!> after get_close_obs has evaluated them
!> 
!> @param[out]   close_indices(:)
!> The state vector indices of the close state variables.
!> 
!> @param[out]   distances(:)
!> The distances of the close state variables from the observation.


subroutine get_close_obs_distrib(filt_gc, base_obs_loc, base_obs_type, locs, kinds, &
                            num_close, close_indices, distances)

! get_close_obs takes as input an "observation" location, a DART TYPE (not KIND),
! and a list of all potentially close locations and KINDS on this task.
!
! get_close_obs
!    *) converts vertical coordinates as needed to vert_coord,
!    *) calls location_mod/threed_sphere:get_close_obs,
!       to which it sends this (converted) array of locations,
!    *) gets back the distances and indices of those locations that are
!       "close" to the base observation.
!    *) tests for being above the highest_obs_pressure_Pa threshold,
!       and increases the vertical distance based on height above highest_*.
!
! get_close_obs will use the ensemble average to convert the obs and/or state
!               vertical location(s) to a standard (vert_coord) vertical location

type(get_close_type), intent(in)    :: filt_gc
type(location_type),  intent(in)    :: base_obs_loc
integer,              intent(in)    :: base_obs_type
type(location_type),  intent(inout) :: locs(:)
integer,              intent(in)    :: kinds(:)
integer,              intent(out)   :: num_close
integer,              intent(out)   :: close_indices(:)
real(r8),             intent(out)   :: distances(:)

! FIXME remove some (unused) variables?
integer  :: k, t_ind
integer  :: base_which, local_base_which, obs_which, local_obs_which
integer  :: base_obs_kind
real(r8) :: base_array(3), local_base_array(3), obs_array(3), local_obs_array(3)
real(r8) :: damping_dist, threshold, thresh_wght
type(location_type) :: local_base_obs_loc, local_loc

if (.not. module_initialized) call static_init_model()

! If base_obs vert type is not pressure; convert it to pressure
base_which    = nint(query_location(base_obs_loc))
base_array    = get_location(base_obs_loc)
base_obs_kind = get_obs_kind_var_type(base_obs_type)

! Upgrading convert_vert to use field profiles at the actual ob location is
! probably not worthwhile: that approx horiz location of the obs is used only to
! convert its vert coord to pressure (if necessary),
! which, in turn, is used to modify the distance if the ob or model variable
! is higher than highest_XXX_Pa.  That modification tapers to 0,
! so any errors introduced by this approx will be continuous and random,
! introducing no bias.
if (base_which == VERTISPRESSURE .and. vert_coord == 'pressure') then
   local_base_obs_loc = base_obs_loc
   local_base_array   = get_location(base_obs_loc)  ! needed in num_close loop
   local_base_which   = base_which
else
   call convert_vert(base_array, base_which, base_obs_loc, base_obs_kind, &
                     local_base_array, local_base_which)
   local_base_obs_loc = set_location(base_array(1), base_array(2), local_base_array(3), &
                                     local_base_which)
endif

! Get all the potentially close obs but no distances. 
call loc_get_close_obs(filt_gc, local_base_obs_loc, base_obs_type, locs, kinds, &
                       num_close, close_indices)

do k = 1, num_close

   ! The indices in close_obs refer to the subset of (state) vars or obs ON 1 TASK.
   ! That subset is (re)labeled 1...num_vars_task#, where num_vars_task# ~ state_vec_size/ntasks.
   ! So those indices can't tell me which state vector element I'm dealing with.
   ! I need to use the location of each returned close_indices to learn anything about it.

   t_ind = close_indices(k)
   obs_array = get_location(locs(t_ind))
   ! query_location returns location%which_vert, if not 'attr' argument is given.
   obs_which = nint(query_location(locs(t_ind)))

   ! FIXME Nancy; what about 'ob's on scale height, but vert_coord is pressure.
   !              KDR: the base ob was converted to pressure, if necessary, in the first section, 
   !                   before the loop over num_close.
   !              And can these if blocks be collapsed by defining local_obs_array(1:2 at least)
   !              before the if tests.
   if ((obs_which == VERTISPRESSURE    .and. vert_coord == 'pressure') .or. &
       (obs_which == VERTISSCALEHEIGHT .and. vert_coord == 'log_invP')) then
      ! put the vertical (pressure) of the state/ob in local storage
      local_obs_array(3) = obs_array(3)
      local_obs_which    = obs_which

   elseif (obs_which == VERTISUNDEF) then
      ! obs_which = -2 (VERTISUNDEF) means this ob is vertically close to base_obs, no matter what.
      ! if (local_obs_array(3) == MISSING_R8) then
      local_obs_array(3) = local_base_array(3)
      local_obs_which    = local_base_which

   else
      call convert_vert(obs_array, obs_which, locs(t_ind), kinds(t_ind), &
                        local_obs_array, local_obs_which)

      ! save the converted location back into the original list.
      ! huge improvement in speed since we only do the vertical convert
      ! once per location, instead of num_close * nobs times.
      locs(t_ind) = set_location( local_obs_array(1), local_obs_array(2), &
                                  local_obs_array(3), local_obs_which)

   endif

!  FIXME: I think this line could be replaced by moving 'locs(t_ind) = ' 
!         out of the last if-block above, and referencing locs(t_ind) below.
!         This is because the lon and lat are not changing: obs_array(1) = local_obs_array(1),...
   local_loc = set_location(obs_array(1), obs_array(2), local_obs_array(3), &
                                   local_obs_which)

!  allow a namelist specified kind string to restrict the impact of those
!  obs kinds to only other obs and state vars of the same kind.
   if ((impact_kind_index >= 0)                .and. &
       (impact_kind_index == base_obs_kind)    .and. &
       (impact_kind_index /= kinds(t_ind))) then
      distances(k) = 999999.0_r8     ! arbitrary very large distance

   else
      ! Need to damp the influence of all obs (VERTISUNDEF, VERTISSURFACE too) on model state vars
      ! above highest_state_pressure_Pa.

      ! The which vert of local_base_obs_loc determines how vertical distance to local_loc is calculated.
      ! It can be VERTISSCALEHEIGHT.
      distances(k) = get_dist(local_base_obs_loc, local_loc, base_obs_type, kinds(t_ind))

      ! Damp the influence of obs, which are below the namelist variable highest_OBS_pressure_Pa,
      ! on variables above highest_STATE_pressure_Pa.
      ! This section could also change the distance based on the KIND_s of the base_obs and obs.

      ! distances = 0 for some for synthetic obs.

      ! Better damping
      ! Should be units of distance (radians), so normalize the distance added to the existing dist(k),
      ! below, by the vert_normalization_{pressure,scale_height}.
      ! Vert_norm is not public, so call get_dist with 2 locations having the same
      ! horiz location, but different verticals, and the appropriate which_vert.

      if ((vert_coord == 'pressure' .and. (local_obs_array(3) < highest_state_pressure_Pa)) .or. &
          (vert_coord == 'log_invP' .and. (local_obs_array(3) > highest_state_scale_h))    ) then
         ! The (lon,lat) here must match the definition of highest_state_loc in static_init_mod.
         ! FIXME; is this hard-coding OK?
         ! local_obs_which should be consistent with local_base_obs_which, (and vert_coord).
         vert_only_loc = set_location(1.0_r8,1.0_r8,local_obs_array(3),local_obs_which)

         ! This gets the vertical distance (> 0) only, and uses the appropriate 
         ! vert_normalization to convert from pressure or scale_height to radians.
         damping_dist = get_dist(highest_state_loc,vert_only_loc,no_vert=.false.)

         ! This (new) added distance varies smoothly from 0 at highest_state_pressure_Pa 
         ! to > 2*cutoff*vert_normalization at the levels where CAM has extra damping 
         ! (assuming that highest_state_pressure_Pa has been chosen low enough).
   
         distances(k) = distances(k) + damping_dist * damping_dist * damp_wght

      endif
   endif

enddo

end subroutine get_close_obs_distrib

!-----------------------------------------------------------------------

subroutine convert_vert_distrib(state_ens_handle, old_array, old_which, old_loc, old_kind, new_array, new_which)

! Uses model information and subroutines to convert the vertical location of an ob
! (prior, model state variable, or actual ob) into the standard vertical coordinate
! (pressure or log_invP = log(P0/ps)).
! Kevin Raeder 10/26/2006
! updated 2014 for WACCM use; log_invP vertical coordinate.

type(ensemble_type),    intent(in)    :: state_ens_handle
real(r8), intent(in)    :: old_array(3)
integer,  intent(in)    :: old_which
type(location_type),    intent(in)    :: old_loc
integer,  intent(in)    :: old_kind
real(r8), intent(inout) :: new_array(3)
integer,  intent(out)   :: new_which

integer  :: num_levs, top_lev, bot_lev
integer  :: istatus, closest
integer  :: cell_corners(4)
integer  :: lon_ind, lat_ind, cam_type
real(r8) :: p_surf, frac, l, m, lon_lat_vert(3)
type(location_type)   :: temp_loc


character(len=8) :: cam_varname

! Don't initialize these in the declaration statements, or they'll get the save attribute
! and won't be initialized each time this routine is entered.
cell_corners = MISSING_I    ! corners of the cell which contains the ob

!HK not building ps arrays.

! this code does not alter the lat/lon, only the vertical.
! but still return a full location for subsequent use.
new_array(1) = old_array(1)
new_array(2) = old_array(2)

! these should be set by the code below; it's an error if not.
new_which    = MISSING_I
new_array(3) = MISSING_R8
num_levs     = lev%length

if (.not. (old_which == VERTISPRESSURE .or. old_which == VERTISHEIGHT  .or. &
           old_which == VERTISLEVEL    .or. old_which == VERTISSURFACE .or. &
           old_which == VERTISUNDEF    .or. old_which == VERTISSCALEHEIGHT) ) then
   ! There's no procedure to translate a which_vert value into text.
   ! So I'll just point users to location_mod.
   write(string1,'(A,3(F12.5,1x),A,I2)') 'obs at (', old_array,  &
        ') has unsupported vertical type = ',old_which
   write(string2,*) 'See location_mod.f90; VERTISxxx to decode this vertical type'
   call error_handler(E_ERR, 'convert_vert', string1,source,revision,revdate,text2=string2)
endif

! Need lon and lat indices to select ps for calc of p_col for vertical conversion.

! Find the surface pressure and the column of pressures at this location.
if (old_which == VERTISLEVEL ) then
   ! This assumes that if VERTISLEVEL, then the potentially close 'ob' is a 
   ! model state variable, defined at a grid point.  So we can figure out the
   ! grid indices and grab the surface pressure from the global ps array.
   ! This may not be true; there can be observations on a model level which
   ! don't lie on a grid point.  Then this vertical coordinate conversion
   ! will be more approximate than if we interpolated the pressure to the 
   ! actual 'ob' horizontal location.
   cam_type = dart_to_cam_types(old_kind)
   if ( cam_type < 0 ) then
      write(string1,*)'old_kind  is ',old_kind,' | cam_type is ',cam_type
      write(string2,*)'get_raw_obs_kind_name of old_kind ', trim(get_raw_obs_kind_name(old_kind))
      call error_handler(E_ERR,'convert_vert',string1,source,revision,revdate,text2=string2)
   endif

   if (l_rectang) then
      ! Assumes 2D obs locations are (lon, lat) and 3D are (lev,lon,lat).

      ! Get the column of pressures at this location, from the ensemble mean.

      cam_varname = trim(cflds(cam_type))
      if (cam_varname == 'US') then
         call coord_index('lon', old_array(1), lon_ind)
         call coord_index('slat', old_array(2), lat_ind)
         !p_surf = ps_stagr_lat(lon_ind,lat_ind)
         p_surf = 0.5*(get_surface_pressure(state_ens_handle, lon_index, lat_index) + &
                       get_surface_pressure(state_ens_handle, lon_index, lat_index +1) )
         ! WHAT ABOUT FIELDS THAT MIGHT COME ON ilevS ?   have lev_which_dimid from above;
         !     test = ilev%dim_id or lev%dim_id
         call plevs_cam(p_surf, num_levs, p_col)
      elseif (cam_varname == 'VS') then
         call coord_index('slon', old_array(1), lon_ind)
         call coord_index('lat', old_array(2), lat_ind)
         !p_surf = ps_stagr_lon(lon_ind,lat_ind)
         if ( lon_index == 1 ) then
            p_surf = 0.5*(get_surface_pressure(state_ens_handle, lon_index, lat_index) + &
                  get_surface_pressure(state_ens_handle, dim_sizes(slon_index), lat_index) )
         else
            p_surf = 0.5*(get_surface_pressure(state_ens_handle, lon_index -1, lat_index) + &
                  get_surface_pressure(state_ens_handle, lon_index, lat_index) )
         endif
         call plevs_cam(p_surf, num_levs, p_col)
      else
         ! Already have the column of 3D pressure. HK - not anymore.
         call coord_index('lon', old_array(1), lon_ind)
         call coord_index('lat', old_array(2), lat_ind)
         !p_surf = ps(lon_ind,lat_ind)
         p_surf = get_surface_pressure(state_ens_handle, lon_index, lat_index)
         p_col(1:num_levs) = p(1:num_levs,lon_ind,lat_ind)
      endif
   else
      ! Cubed sphere; more complicated search for indices of this location.
      ! The 3D index into the state vector is NOT known.
      ! We have the index into the subset of state vars ON 1 TASK.
      ! That subset has its own indexing, which is what's available here.
      ! The relationship to global indices is stuck back in filter.
      ! Distributed state version: we *will* know this overall state vector index
      !    in the future.

      call coord_ind_cs(old_loc, old_kind, .true., closest, cell_corners, l, m)

      ! Use surface and 3D pressure from ens_mean.
      p_surf = ps(closest,1)
      p_col(1:num_levs) = p(1:num_levs,closest,1)
   endif
else !HK so if oldwhich==VERTISLEVEL, it is surface?

   ! Make a vertical that has a vert type of surface.
   lon_lat_vert = get_location(old_loc)
   temp_loc = set_location(lon_lat_vert(1), lon_lat_vert(2), 0.0_r8, VERTISSURFACE)
   ! Find ps at the ob point.  Need to interpolate.
   ! HK don't want to call the interpolation routine.
   if (l_rectang) then
      ! Only interested in P (columns), so don't need to worry about staggered grids here.
      !call interp_lonlat(ens_mean, temp_loc, KIND_SURFACE_PRESSURE, p_surf, istatus)
       call interp_lonlat_distrib(state_ens_handle, location, obs_kind, istatus, interp_val)


   else
      call error_handler(E_ERR, 'no rma cubed sphere')
      !call interp_cubed_sphere(ens_mean, temp_loc, KIND_SURFACE_PRESSURE, p_surf, istatus, cell_corners, l, m)
   endif

   if (istatus == 1) then
      write(string1,'(A,I8)') 'interp_X failed for KIND_SURFACE_PRESSURE.'
      call write_location(0, old_loc, charstring=string2)
      call error_handler(E_ERR, 'convert_vert', string1,source,revision,revdate, text2=string2)
   endif

   call plevs_cam(p_surf, num_levs, p_col)

endif

! Convert vertical coordinate to vert_coord (pressure or log_invP).
if (old_which == VERTISUNDEF) then
   ! Field with no vertical location; get_dist will only calculate horiz dist unless
   ! this case is handled by the calling routine.

   ! If a parameter/state variable is supposed to be close to everything,
   ! then I would need to have the/an other location to set it to,
   ! Send back new_array empty and test for that in the calling routine,
   ! where the other location exists.
   ! For model variables user specifies which_vert for each state field,
   ! so when user specifies undefined, then this should return;
   new_array(3) = MISSING_R8
   new_which    = old_which

elseif (old_which == VERTISSURFACE) then
   ! surface field; change which_vert for the distance calculation
   if (vert_coord == 'pressure') then
      new_array(3) =  p_surf
      new_which = VERTISPRESSURE
   elseif (vert_coord == 'log_invP') then
    ! Scale height at the surface is 0.0_r8 by definition [log(p_surf/p_surf)]
      new_array(3) = 0.0_r8
      new_which = VERTISSCALEHEIGHT
   endif

elseif (old_which == VERTISPRESSURE) then
   if (vert_coord == 'log_invP') then
      new_array(3) = scale_height(p_surface=p_surf, p_above=old_array(3))
      new_which = VERTISSCALEHEIGHT
   endif

elseif (old_which == VERTISSCALEHEIGHT) then
   if (vert_coord == 'pressure') then
      new_array(3) = p_surf / exp(old_array(3))
      new_which = VERTISPRESSURE
   endif

elseif (old_which == VERTISLEVEL) then
   ! FIXME
   ! WHAT ABOUT FIELDS THAT MIGHT COME ON ilevS ?   have lev_which_dimid from above;
   !     test = ilev%dim_id or lev%dim_id
   ! OR do this for all columns in static_init_model_dist, which would make PS (and P)
   ! globally available for all regions?
   if (vert_coord == 'pressure') then
      new_array(3) =            p_col(nint(old_array(3)))
      new_which = VERTISPRESSURE
   elseif (vert_coord == 'log_invP') then
      new_array(3) = scale_height(p_surface=p_surf, p_above=p_col(nint(old_array(3))))
      new_which = VERTISSCALEHEIGHT
   endif

elseif (old_which == VERTISHEIGHT) then

   ! Ens_mean is global storage that should have been filled
   ! by a call from filter_assim to ens_mean_for_model.
   ! HK model_h is a global. Why are you passing it in?
   call model_heights_distrib(num_levs, state_ens_handle, p_surf, old_loc,  model_h, istatus)
   if (istatus == 1) then
      write(string1, *) 'model_heights failed'
      call error_handler(E_ERR, 'convert_vert', string1)
      ! return
   endif

   ! Search down through heights
   ! This assumes linear relationship of pressure to height over each model layer,
   ! when really it's exponential.  How bad is that?
!   bot_lev = 2
!   do while (old_array(3) <= model_h(bot_lev) .and. bot_lev <= num_levs)
!      bot_lev = bot_lev + 1
!   enddo
   Bottom: do bot_lev = 2,num_levs+1
      if (old_array(3) <= model_h(bot_lev)) exit Bottom
   enddo Bottom
   top_lev = bot_lev - 1

   ! Write warning message if not found within model level heights.
   ! Maybe this should return failure somehow?
   if (top_lev == 1 .and. old_array(3) > model_h(1)) then
      ! above top of model
      frac = 1.0_r8
      write(string1, *) 'ob height ',old_array(3),' above CAM levels at ' &
                          ,old_array(1) ,old_array(2) ,' for kind',old_kind
      call error_handler(E_MSG, 'convert_vert', string1,source,revision,revdate)
   elseif (bot_lev <= num_levs) then
      ! within model levels
      frac = (old_array(3) - model_h(bot_lev)) / (model_h(top_lev) - model_h(bot_lev))
   else
      ! below bottom of model
      frac = 0.0_r8
      write(string1, *) 'ob height ',old_array(3),' below CAM levels at ' &
                          ,old_array(1) ,old_array(2) ,' for kind',old_kind
      call error_handler(E_MSG, 'convert_vert', string1,source,revision,revdate)
   endif

   new_array(3) = (1.0_r8 - frac) * p_col(bot_lev) + frac * p_col(top_lev)

   if (vert_coord == 'pressure') then
      new_which = VERTISPRESSURE
   elseif (vert_coord == 'log_invP') then
      new_array(3) = scale_height(p_surface=p_surf, p_above=new_array(3))
      new_which = VERTISSCALEHEIGHT
   endif

else
   write(string1, *) 'model which_vert = ',old_which,' not handled in convert_vert '
   call error_handler(E_ERR, 'convert_vert', string1,source,revision,revdate)
endif

return

end subroutine convert_vert_distrib

!------------------------------------------------------------
! Subroutines from mpas_atm/model_mod.f90, for using cartesian coordinates to
! find closest node to an ob.

subroutine init_closest_node()

! use ncol, lats and lons of nodes (corners) to initialize a get_close structure
! to be used later in find_closest_node().

integer :: i

allocate(cs_locs_xyz(ncol))

do i=1, ncol
   cs_locs_xyz(i) = xyz_set_location(lon%vals(i), lat%vals(i), 0.0_r8, earth_radius)
enddo

! the width (2nd arg of ...init) really isn't used anymore, but it's part of the
! interface so we have to pass some number in.
call xyz_get_close_maxdist_init(cs_gc_xyz, 1.0_r8)
call xyz_get_close_obs_init    (cs_gc_xyz, ncol, cs_locs_xyz)

end subroutine init_closest_node

!------------------------------------------------------------

function find_closest_node(lat, lon)

! Determine the index for the closest node to the given point
! 2D calculation only.

real(r8), intent(in)  :: lat
real(r8), intent(in)  :: lon
integer               :: find_closest_node

type(xyz_location_type) :: pointloc
integer                 :: closest_node, rc
! This 'save' is redundant with initializing the variable here in the declaration statement.
logical, save           :: search_initialized = .false.

! do this exactly once.
if (.not. search_initialized) then
   call init_closest_node()
   search_initialized = .true.
endif

pointloc = xyz_set_location(lon, lat, 0.0_r8, earth_radius)

call xyz_find_nearest(cs_gc_xyz, pointloc, cs_locs_xyz, closest_node, rc)

! decide what to do if we don't find anything.
if (rc /= 0 .or. closest_node < 0) then
   if (output_task0) then
      write(string1,*) 'cannot find a nearest node to lon, lat: ', lon, lat
      call error_handler(E_WARN, 'find_closest_node', string1,source,revision,revdate)
      ! newFIXME; should this be E_ERR instead?
   endif
   find_closest_node = -1
   return
endif

! this is the cell index for the closest center
find_closest_node = closest_node

end function find_closest_node

!------------------------------------------------------------

subroutine finalize_closest_node()

! get rid of storage associated with GC for cell centers.

call xyz_get_close_obs_destroy(cs_gc_xyz)

end subroutine finalize_closest_node


! End of get_close_obs section

!#######################################################################

! Initial conditions for DART

!-----------------------------------------------------------------------
!>
!> Subroutine pert_model_state is used for generating initial ensembles
!> by perturbing a model state.  Each ensemble member calls this
!> routine separately through the ensemble_manager.
!> Returning interf_provided means go ahead and do this with
!> small independent perturbations. 
!> It is controlled by model_nml namelist variables.
!> There are two modes of perturbation.  The most common will perturb
!> every state variable by a small random amount.  
!> See model_mod.html for details.
!> 
!> @param[in]    state(:)
!> The model state which will be perturbed
!> 
!> @param[out]   pert_state(:)
!> The perturbed model state
!> 
!> @param[out]   interf_provided
!> A flag to tell filter that this perturbation interface has been provided to it.

  subroutine pert_model_state(state, pert_state, interf_provided)

! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with
! small independent perturbations. Each ensemble member calls this
! routine separately through the ensemble_manager.
!
! If module storage variable 'pert_sd' is positive, then we will randomly
! perturb the fields (based on the values of pert_sd) for each of
! the variable names listed in pert_names.
!
! If 'pert_sd' is negative (which includes MISSING) then the field (only one)
! listed in pert_names is set to a different constant value for each
! ensemble member.  Those values come from 'pert_base_vals'.

real(r8), intent(in)    :: state(:)
real(r8), intent(out)   :: pert_state(:)
logical,  intent(out)   :: interf_provided

type(model_type)        :: var_temp
integer                 :: i, j, k, m, pert_fld, mode, field_num
integer                 :: dim1, dim2, dim3, member
integer, save           :: seed
logical                 :: perturbed

! FIX for 1D 0D  fields?

! Perturb model state vector values for the filter initial conditions.
! The input is a single model state vector that has (different) gaussian
! noise added to each member to generate an initial ensemble.

if (.not. module_initialized) call static_init_model()

interf_provided = .true.

call init_model_instance(var_temp)
call vector_to_prog_var(state,var_temp)

! If first call, then initialize a seed to use for initializing random sequences.
if (first_pert_call) then
   ! This line generates a unique base number, and subsequent calls add 1
   ! each time (which happens if there are multiple ensemble members/task).
   ! It is assuming there are no more than 1000 ensembles/task, which seems safe
   ! given the current sizes of state vecs and hardware memory.  This will make
   ! the results reproduce for runs with the same number of MPI tasks.  It will
   ! NOT give the same random sequence if you change the task count.
   seed = (my_task_id()+1) * 1000
   first_pert_call = .false.
endif

! After using the seed, increment by one so if this routine is called again
! for a different ensemble member it will generate a different set of nums.
call init_random_seq(random_seq, seed)
seed = seed + 1

pert_fld = 1

Vars2Perturb : do pert_fld=1,100
   if (pert_names(pert_fld) == ' ') exit Vars2Perturb

   ! Keep track of whether or not this field is matched and was perturbed.
   perturbed = .false.

   ExistingVars : do m=1,nflds

      if (pert_names(pert_fld) /= cflds(m)) cycle ExistingVars

      perturbed = .true.

      call error_handler(E_MSG,'pert_model_state', 'Perturbing '//trim(pert_names(pert_fld)))

      ! FIXME : below is not robust. ens_member is always 0 in CESM context.
      !         Probably should remove this option from future versions; hasn't been used for years.

      !  Choose mode of perturbations/resets;
      if (pert_sd(pert_fld) <= 0.0_r8 ) then
         ! Set each ensemble member to its own constant value,
         ! as found in pert_base_vals(ens_member).
         ! This only works when setting a single field = to a different constant value
         ! for each ensemble member.
         ! Could add more fields by overloading pert_base_vals and
         ! adding code to find those values.
         mode = ens_member
      else
         ! Set each *field* to it's own pert_base_val +/- pert_sd
         mode = pert_fld
      endif

      ! Handle the 0d fields
      if (m <= state_num_0d) then
         field_num = m

         ! reset base values to value provided in namelist.
         if ( pert_base_vals(mode) /= MISSING_R8 ) then
            if (print_details) then
               write(string1,*) 'Using a new base value ',pert_base_vals(mode), 'for ',cflds(m)
               call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
            endif
            var_temp%vars_0d(field_num) = pert_base_vals(mode)
         endif

         if (print_details) then
             write(string1,'(A,A8,A3,1x,2E24.15)') 'org first and last state for ',&
                    cflds(m),' = ', var_temp%vars_0d(field_num)
             call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
         endif

         ! randomly perturb each point around its base value.
         if (pert_sd(pert_fld) > 0.0_r8 ) then
            do i = 1, dim1
               var_temp%vars_0d(field_num) = &
                  random_gaussian(random_seq, var_temp%vars_0d(field_num),pert_sd(mode))
            enddo
         endif

         if (print_details) then
            write(string1,'(A,A8,A3,1x,2E24.15)') 'new first and last state for ',&
                    cflds(m),' = ', var_temp%vars_0d(field_num)
            call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
         endif

      ! Handle the 1d fields
      elseif (m <= state_num_1d + state_num_0d) then
         field_num = m - state_num_0d
         dim1 = dim_sizes(s_dimid_1d(field_num))

         ! reset base values to value provided in namelist.
         if ( pert_base_vals(mode) /= MISSING_R8 ) then
            if (print_details) then
               write(string1,*) 'Using a new base value ',pert_base_vals(mode), 'for ',cflds(m)
               call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
            endif
            var_temp%vars_1d(1:dim1,field_num) = pert_base_vals(mode)
         endif

         if (print_details) then
             write(string1,'(A,A8,A3,1x,2E24.15)') 'org first and last state for ',&
                    cflds(m),' = ', var_temp%vars_1d(   1,field_num), &
                                    var_temp%vars_1d(dim1,field_num)
             call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
         endif

         ! randomly perturb each point around its base value.
         if (pert_sd(pert_fld) > 0.0_r8 ) then
            do i = 1, dim1
               var_temp%vars_1d(i,field_num) = &
                   random_gaussian(random_seq,var_temp%vars_1d(i,field_num),pert_sd(mode))
            enddo
         endif

         if (print_details) then
            write(string1,'(A,A8,A3,1x,2E24.15)') 'new first and last state for ',&
                    cflds(m),' = ', var_temp%vars_1d(   1,field_num), &
                                    var_temp%vars_1d(dim1,field_num)
            call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
         endif

      elseif (m <= state_num_2d + state_num_1d + state_num_0d) then
         field_num = m - state_num_1d - state_num_0d
         dim1 = dim_sizes(s_dimid_2d(1,field_num))
         dim2 = dim_sizes(s_dimid_2d(2,field_num))

         ! reset base values to value provided in namelist.
         if ( pert_base_vals(mode) /= MISSING_R8 ) then
            if (print_details) then
               write(string1,*) 'Using a new base value ',pert_base_vals(mode), 'for ',cflds(m)
               call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
            endif
            var_temp%vars_2d(1:dim1,1:dim2,field_num) = pert_base_vals(mode)
         endif

         if (print_details) then
             write(string1,'(A,A8,A3,1x,2E24.15)') 'org first and last state for ',&
                    cflds(m),' = ', var_temp%vars_2d(   1,   1,field_num), &
                                    var_temp%vars_2d(dim1,dim2,field_num)
             call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
         endif

         ! randomly perturb each point around its base value.
         if (pert_sd(pert_fld) > 0.0_r8 ) then
            do j = 1, dim2
            do i = 1, dim1
               var_temp%vars_2d(i,j,field_num) = &
                  random_gaussian(random_seq, var_temp%vars_2d(i,j,field_num),pert_sd(mode))
            enddo
            enddo
         endif

         if (print_details) then
            ! newFIXME; Nancy wants to see min and max instead of (addition to?)first and last.
            ! 3d and 1d  and 0d too.
            write(string1,'(A,A8,A3,1x,2E24.15)') 'new first and last state for ',&
                    cflds(m),' = ', var_temp%vars_2d(   1,   1,field_num), &
                                    var_temp%vars_2d(dim1,dim2,field_num)
            call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
         endif

      else ! do the 3D fields

         field_num = m - state_num_2d - state_num_1d - state_num_0d
         dim1 = dim_sizes(s_dimid_3d(1,field_num))
         dim2 = dim_sizes(s_dimid_3d(2,field_num))
         dim3 = dim_sizes(s_dimid_3d(3,field_num))

         if (print_details) then
            write(string1,'(A,A8,A3,1x,2E24.15)') 'org first and last state for ',&
                    cflds(m), ' = ', var_temp%vars_3d(   1,   1,   1,field_num),  &
                                     var_temp%vars_3d(dim1,dim2,dim3,field_num)
            call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
         endif

         ! reset base values to value provided in namelist.
         if ( pert_base_vals(mode) /= MISSING_R8 ) then
            if (print_details) then
               write(string1,*) '  uses a new base value ',pert_base_vals(mode),' for ',cflds(m)
               call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
            endif
            var_temp%vars_3d(1:dim1,1:dim2,1:dim3,field_num) = pert_base_vals(mode)
         endif

         ! randomly perturb each point around its base value.
         if (pert_sd(pert_fld) > 0.0_r8 ) then
            if (print_details) then
               write(string1,*) 'Perturbing base value of ',cflds(m),' by st dev ',pert_sd(mode)
               call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
            endif

            do j = 1, dim3
            do i = 1, dim2
            do k = 1, dim1
               ! new val = rand#(O(0-1)) * standard dev  + mean
               var_temp%vars_3d(k,i,j,field_num) = &
                  random_gaussian(random_seq, var_temp%vars_3d(k,i,j,field_num),pert_sd(mode))
            enddo
            enddo
            enddo
         endif

         if (print_details) then
            write(string1,'(A,A8,A3,1x,2E24.15)') 'new first and last state for ',&
                    cflds(m), ' = ', var_temp%vars_3d(   1,   1,   1,field_num),  &
                                     var_temp%vars_3d(dim1,dim2,dim3,field_num)
            call error_handler(E_MSG, 'pert_model_state', string1, source, revision, revdate)
         endif

      endif

   enddo ExistingVars

   if (.not. perturbed) then
      write(string1,*)trim(pert_names(pert_fld)),' not found in list of state variables.'
      write(string2,*)'but was supposed to be used to perturb.'
      call error_handler(E_ERR,'pert_model_state', string1, source, revision, revdate, text2=string2)
   endif

enddo Vars2Perturb

call prog_var_to_vector(var_temp,pert_state)
call end_model_instance(var_temp)

end subroutine pert_model_state


!-----------------------------------------------------------------------
!>
!> Subroutine init_conditions
!> reads in restart initial conditions  -- noop for CESM atmospheric components.
!> 
!> @param[inout] st_vec(:)
!> The state vector which is NOT read from a file by this routine.

subroutine init_conditions(st_vec)

! Reads in restart initial conditions  -- noop for CAM

real(r8), intent(inout) :: st_vec(:)

if (.not. module_initialized) call static_init_model()

call error_handler(E_ERR,"init_conditions", &
                  "WARNING!!  CAM model has no built-in default state", &
                  source, revision, revdate, &
                  text2="cannot run with 'start_from_restart = .false.'", &
                  text3="use 'cam_to_dart' to create a CAM state vector file")

st_vec = MISSING_R8  ! just to silence compiler messages

end subroutine init_conditions



! End of initial model state section

!#######################################################################

! Utility routines; called by several main subroutines

!-----------------------------------------------------------------------

function index_from_grid(lev_ind, lon_ind, lat_ind, ifld)

! Calculate the index into the state vector, given the coordinate indices
! and the field number (out of nflds state vector fields).

! This could be more efficient by calculating and globally storing the beginning index
! of each ifld and/or the size of each ifld.

integer, intent(in) :: lev_ind
integer, intent(in) :: lon_ind
integer, intent(in) :: lat_ind
integer, intent(in) :: ifld
integer             :: index_from_grid

integer :: i, j, fld_ind(2), done

index_from_grid = 0
done = 0

! Cycle through 0d state variables
do i=1,state_num_0d
   index_from_grid = index_from_grid + 1
   if (ifld == i) return
enddo
done = done + state_num_0d

! Cycle through 1d state variables
! Note that indices of fields can have varying dimensions.
! FIXME; replace multiple ifs with a case structure.  Or will this be replaced?
do i=1,state_num_1d
   if (ifld - done == i) then
         ! FIXME: use select cases here, and add a failure case
      if (dim_names(s_dimid_1d(i)) == 'lon' .or. &
          dim_names(s_dimid_1d(i)) == 'slon') index_from_grid = index_from_grid + lon_ind
      if (dim_names(s_dimid_1d(i)) == 'lat' .or. &
          dim_names(s_dimid_1d(i)) == 'slat') index_from_grid = index_from_grid + lat_ind
      if (dim_names(s_dimid_1d(i)) == 'lev' .or. &
          dim_names(s_dimid_1d(i)) == 'ilev') index_from_grid = index_from_grid + lev_ind
      ! CS lon_ind has been pirated by ncol.
      ! FIXME; replace xxx_ind with non-specific names; ind1...
      !        Trace the xxx_ind names back through calls (get_val_...).
      !        Or will this be replaced by separate lonlat and cubed-sphere versions?
      if (dim_names(s_dimid_1d(i)) == 'ncol') index_from_grid = index_from_grid + lon_ind
      return
   else
      index_from_grid = index_from_grid + s_dim_1d(i)
   endif
enddo
done = done + state_num_1d

! Cycle through 2d state variables.
! Note that indices of fields can have varying dimensions.
do i=1,state_num_2d
   if (ifld - done == i) then
      ! We've found the desired field; now find index of lev and/or lon and/or lat
      do j=1,2
         ! FIXME: use select cases here, and add a failure case
         if (dim_names(s_dimid_2d(j,i)) == 'lon' .or. &
             dim_names(s_dimid_2d(j,i)) == 'slon'     ) fld_ind(j) = lon_ind
         if (dim_names(s_dimid_2d(j,i)) == 'lat' .or. &
             dim_names(s_dimid_2d(j,i)) == 'slat'     ) fld_ind(j) = lat_ind
         if (dim_names(s_dimid_2d(j,i)) == 'lev' .or. &
             dim_names(s_dimid_2d(j,i)) == 'ilev'     ) fld_ind(j) = lev_ind
         ! CS lon_ind has been pirated by ncol.
         if (dim_names(s_dimid_2d(j,i)) == 'ncol')     fld_ind(j) = lon_ind
      enddo

      index_from_grid = index_from_grid + (fld_ind(2)-1) * s_dim_2d(1,i) + fld_ind(1)
      return
   else
      index_from_grid = index_from_grid +  s_dim_2d(2,i) * s_dim_2d(1,i)
   endif
enddo
done = done + state_num_2d


! Cycle through 3d state variables
! Note that indices of fields can have varying dimensions.
! CS There won't be any 3d fields for the cubed sphere.
do i=1,state_num_3d
   if (ifld - done == i) then
      ! We've found the desired field; now find index of lat, lon, lev
      index_from_grid = index_from_grid  &
                      + (lat_ind-1) * s_dim_3d(2,i) * s_dim_3d(1,i) &
                      + (lon_ind-1) * s_dim_3d(1,i)                 &
                      + lev_ind
      return
   else
      index_from_grid = index_from_grid + s_dim_3d(3,i) * s_dim_3d(2,i) * s_dim_3d(1,i)
   endif
enddo

end function index_from_grid

!-----------------------------------------------------------------------

function find_name(nam, list)

character(len=*), intent(in) :: nam
character(len=*), intent(in) :: list(:)
integer                      :: find_name

integer :: i

find_name = 0
do i = 1,size(list,1)
if (nam == list(i)) then
   find_name = i
   return
endif
enddo

end function find_name

!-----------------------------------------------------------------------

subroutine coord_val(dim_name, indx, lon_val, lat_val, lev_val)

! Given the name of the coordinate to be searched and the index into that array,
! returns the coordinate value  in either lon_val, lat_val, or lev_val.
! All 3 _val arguments are present so that this routine can return the value
! in the coordinate that the calling routine wants it to be, and search/placement doesn't
! have to be done there.

character(len=*), intent(in)    :: dim_name
integer,          intent(in)    :: indx
real(r8),         intent(inout) :: lon_val
real(r8),         intent(inout) :: lat_val
real(r8),         intent(inout) :: lev_val

! Check for acceptable value of indx?
! FIXME; replace these ifs with select case and add a failure case.

if (dim_name == 'lon') lon_val = lon%vals(indx)
if (dim_name == 'lat') lat_val = lat%vals(indx)
if (dim_name == 'slon') then
   lon_val = slon%vals(indx)
   ! CAM staggered longitude grid -2.5, ..., 352.5 (FV4x5)
   ! but DART wants to see 0.,...,360.  only.
   if (lon_val < 0.0_r8) lon_val = lon_val + 360.0_r8
endif
if (dim_name == 'slat') lat_val = slat%vals(indx)
if (dim_name == 'ncol') then
   lon_val = lon%vals(indx)
   lat_val = lat%vals(indx)
endif

if (lat_val <= -90.0_r8) lat_val = -89.9999999_r8
if (lat_val >=  90.0_r8) lat_val =  89.9999999_r8

! FIXME this is returning the NOMINAL vertical location (to get_state_meta_data)
!       Is that good enough?  Or do I need to calculate the actual vertical location?
!       This IS good enough for the calls in interp_lonlat because only lat_val is set by those calls.
! 2FIXME: lev from the initial file is for PS = 1000 hPa (not 10^5 Pa); missing topography and weather.
! So it's useless for our purposes.  (It's not even consistent with units of PS).
! if (dim_name == 'lev') lev_val = lev%vals(indx) * 100.0_r8
! if (dim_name == 'ilev') lev_val = ilev%vals(indx) * 100.0_r8
! Any need for the lev pressure values will be calculated in get_close_obs:convert_vert.
if (dim_name == 'lev' .or. dim_name == 'ilev') then
   lev_val = real(indx)
endif
! Add more for other coords?  hyam...?  Not for now; never referenced indirectly

end subroutine coord_val

!-----------------------------------------------------------------------

subroutine coord_ind_cs(obs_loc, obs_kind, closest_only, closest , cell_corners, l, m)

! Find the node closest to a location, and the possibly the corners of the cell which contains 
! the location.

! Variables needed by loc_get_close_obs:
type(location_type),  intent(in)  :: obs_loc
integer,              intent(in)  :: obs_kind
logical,              intent(in)  :: closest_only
integer,              intent(out) :: closest
integer,              intent(out) :: cell_corners(4)
real(r8),             intent(out) :: l
real(r8),             intent(out) :: m

! Output from loc_get_close_obs
integer  :: num_close

! It would be nice if these could be smaller, but I don't know what number would work.
! It has to be large enough to accommodate all of the grid points that might lie
! within 2xcutoff; resolution and location dependent.
! The size must be specified here; (:) yields an error, and 'allocatable' doesn't help.
integer, allocatable  :: close_ind(:)
real(r8), allocatable :: dist(:)

! Local Variables
! dist_# in radians (Can't be initialized here or they will get the 'save' property,
! and will not be reset during subsequent entries to this subroutine.)
real(r8) :: dist_1, dist_2
real(r8) :: lon_lat_lev(3)
integer  :: k, k1, k2, closest2, origin
logical  :: found_cell

lon_lat_lev = get_location(obs_loc)

! See whether this obs_ is a state variable.
! This could be done by 2 calls to minloc(dist), with the 2nd call using a mask
! to prevent finding the closest, which was found in the first call.
! But would those 2 intrinsic searches through dist be faster than my 1 explicit search?

if (closest_only) then
   ! Use xyz/cartesian coordinates to quickly find the closest node.
   ! If convert_vert only needs the closest node, don't find the l,m weights.
   closest = find_closest_node(lon_lat_lev(2), lon_lat_lev(1))

   ! Can return without deallocating close_ind and dist
   ! because they haven't been allocated yet.
   return
endif

! Allocate space for the potentially close nodes.
allocate(close_ind(ncol), dist(ncol))

! Look for the 2 closest nodes, using slower way of getting all of the close obs
! and searching for the 2 closest.
! --------------
! FIXME: Nancy has a location_xyz:find_closest_???? which will return the N closest points,
! which may be significantly faster than threed_sphere/location_mod.f90:get_close_obs.
! --------------
! FIXME; can the closest node not be a corner of the containing cell in grids generated by SQuadGen?
! --------------
! For a refined grid (from 1 degree to 1/8 degree) loc_get_close_obs is going to return lists
! that are 64x larger in the refined region than in the coarse region

!   obs_'kind' is passed to location.f90:get_close_obs.
!   There it is passed to only get_dist, which only uses it if special_vert_norm is used,
!   and gc%special_maxdist.
!   Model_mod is not using either of those.
call loc_get_close_obs(cs_gc, obs_loc, obs_kind, cs_locs, cs_kinds, &
                       num_close, close_ind, dist)

dist_1 = 10.0_r8
dist_2 = 10.0_r8
closest = MISSING_I
k1 = MISSING_I

! Keep track of k1, k2, and distances in this search.
! Assign closest and closest2 afterwards.
if (num_close <= 0) then
   write(string1,*) 'Unusable num_close, obs_kind : ',num_close, obs_kind
   call write_location(0, obs_loc, charstring=string2)
   write(string3,*) 'dist(1) = ',dist(1)
   call error_handler(E_ERR, 'coord_ind_cs', string1,source,revision,revdate,text2=string2, text3=string3)
endif

do k = 1,num_close
   if (dist(k) < dist_2) then
      ! Replace 2nd with new one.
      k2 = k
      dist_2 = dist(k)
      if (dist_2 <= dist_1) then
         ! Switch 1st and new 2nd.   '<=' To make sure k1 is filled, even for the first k.
         k2 = k1
         k1 = k
         dist_2 = dist_1
         dist_1 = dist(k)
      endif
   endif
enddo
closest  = close_ind(k1)

if (k2 == MISSING_I) then
   write(string1,'(A)') 'Did not find a second closest node to ob:'
   write(string2,'(A,3F10.2,3I6,1p2E12.4)')                    &
        'lon_lat_lev, obs_kind, num_close, closest, dist_1, dist_2 = ', &
         lon_lat_lev, obs_kind, num_close, closest, dist_1, dist_2
   call write_location(0, cs_locs(closest), charstring=string3)
   string3 = 'closest node location = '//string3
   call error_handler(E_ERR, 'coord_ind_cs', string1,source,revision,revdate,text2=string2,text3=string3)
else
   closest2 = close_ind(k2)
endif

! Find the cell which contains the ob.
! First search the cells which have 'closest' as 1 corner.
! If that fails, search the cells around closest2.
! The search consists of passing the ob location to unit_square_location
! and letting it determine whether the ob location maps into the unit square.

! Initial value of success flag.
found_cell = .false.

! FIXME; debug in verify_namelist
! write(string1,*) 'STARTING Cloop num_nghbrs = ',num_nghbrs(closest)
! call error_handler(E_MSG,'coord_ind_cs',string1,source,revision,revdate)

Cloop: do k=1,num_nghbrs(closest)
   ! centers(k,closest) refers to the cell center name associated with neighboring node k
   ! of the closest node.  It is used to retrieve mapping coefficients for the cell being tested.

   call unit_square_location(centers(k,closest), closest, obs_loc,          &
                             lon_lat_lev(1),lon_lat_lev(2), found_cell, origin, l,m)
   if (found_cell) exit Cloop
enddo Cloop

! Try the 2nd closest point, if the first failed.
if ((.not.found_cell) .and. closest2 /= MISSING_I) then

   Second_closest: do k=1,num_nghbrs(closest2)
      call unit_square_location(centers(k,closest2), closest2, obs_loc,         &
                                lon_lat_lev(1),lon_lat_lev(2), found_cell, origin, l,m)
      if (found_cell) then
         ! Put '2nd closest' information into 'closest'.
         dist_1 = dist_2
         closest = closest2

         write(string1,'(A,2F10.7,2I8,1p2E12.4)') &
              'Using 2nd closest node to the ob: l, m, closest2, origin2 = ', &
              l, m, closest, origin
         call error_handler(E_MSG, 'coord_ind_cs', string1,source,revision,revdate)

         exit Second_closest
      endif
   enddo Second_closest
endif

if (found_cell) then
   ! Need to shift corners according to which was chosen as the origin corner
   ! in num_nghbrs loop, above.  The weighted interp calculation assumes, as in
   ! the create_cs_grid_arrays mapping scheme, that the origin node is corner 4.
   cell_corners(1:4) = cshift(corners(centers(k,closest),1:4), origin)

else
   ! Both closest nodes failed; abort
   write(string1, '(A,2I8,A,2F10.4)') &
         'Neither of the 2 closest nodes ',  closest,closest2, &
         ' is a corner of the cell containing ob at ', lon_lat_lev(1),lon_lat_lev(2)
   call error_handler(E_ERR, 'coord_ind_cs', string1,source,revision,revdate)
endif

deallocate(close_ind, dist)

end subroutine coord_ind_cs

!-----------------------------------------------------------------------

subroutine coord_index(dim_name, val, indx, other_indx)

! Given the name of the (Eulerian or FV) coordinate to be searched and the value,
! Returns the index of the closest coordinate value.
! Optionally returns the next closest index too, which may be < or > the closest.

character(len=*), intent(in)   :: dim_name
real(r8),         intent(in)   :: val
integer,          intent(out)  :: indx
integer, optional, intent(out) :: other_indx

real(r8), pointer :: coord(:)
real(r8)          :: diff_upper, diff_lower, val_local, resol
integer           :: coord_len, i

nullify(coord)
val_local = val

if (dim_name == 'lon') then
   coord     => lon%vals
   coord_len =  lon%length
   resol     =  lon%resolution
elseif (dim_name == 'lat') then
   coord     => lat%vals
   coord_len =  lat%length
   resol     =  lat%resolution
elseif (dim_name == 'lev') then
   coord     => lev%vals
   coord_len =  lev%length
   resol     =  lev%resolution
elseif (dim_name == 'slon') then
   coord     => slon%vals
   coord_len =  slon%length
   resol     =  slon%resolution
   ! Make sure longitudes conform to the CAM staggered grid.
   if ((val - coord(coord_len)) >= (coord(coord_len)-coord(coord_len-1))) &
      val_local = val_local - 360.0_r8
elseif (dim_name == 'slat') then
   coord     => slat%vals
   coord_len =  slat%length
   resol     =  slat%resolution
elseif (dim_name == 'ilev') then
   coord     => ilev%vals
   coord_len =  ilev%length
   resol     =  ilev%resolution
else
   ! should not happen; fatal error.
   write(string1, *) 'unexpected dim_name, ', trim(dim_name)
   call error_handler(E_ERR, 'coord_index', string1,source,revision,revdate)
endif

! Assumes that coordinates are monotonic.

if (val_local <= coord(1)) then
   indx = 1
   if (present(other_indx)) other_indx = 1
   nullify (coord)
   return
elseif (val_local >= coord(coord_len)) then
   indx = coord_len
   if (present(other_indx)) other_indx = coord_len
   nullify (coord)
   return
else
   if (resol > 0.0_r8) then
      ! regularly spaced; calculate the index
      ! NINT is used because some calls to this routine want the single closest indx,
      ! regardless of whether val_local is < or > coord(indx).
      indx = NINT((val_local - coord(1))/resol) + 1

      if (present(other_indx)) then
         if (val_local > coord(indx)) then
            other_indx = indx + 1
         else
            other_indx = indx - 1
         endif
      endif
   else
      ! IRregularly spaced (but still monotonically increasing); search for the index
      ! Replace with a binary search?
      do i=1, coord_len - 1
         diff_upper = coord(i+1) - val_local
         if (diff_upper >= 0.0_r8) then
            diff_lower = val_local - coord(i)
            ! Alway return the closer coord point in the first (non-optional) argument
            if (diff_upper > diff_lower) then
               indx = i
               if (present(other_indx)) other_indx = i + 1
            else
               indx = i + 1
               if (present(other_indx)) other_indx = i
            endif
            nullify (coord)
            return
         endif
      enddo
   endif
endif
! Try reclaiming coord memory before returning.
nullify (coord)

end subroutine coord_index

!-----------------------------------------------------------------------

subroutine set_ps_arrays(vec)

! Assign values to pressure arrays for use by the rest of the module.
! The form(s) of the arrays will depend on the CAM dynamical core being used.
! Spectral/eulerian: ps(lon,lat), p(lev,lon,lat)
! Finite Volume    : ps(lon,lat), p(lev,lon,lat), ps_stagr_lon(slon,lat), ps_stagr_lat(lon,slat)
! Spectral Element : ps(ncol), p(lev,ncol)

real(r8), intent(in) :: vec(:)

integer :: ind, m,n, slon_index, slat_index, lon_index, lat_index, lev_index, &
           fld_index_ps, fld_index, dim1, dim2, dim_lev


if (l_rectang) then
   ! Rectangular grid; 2 horizontal dimensions.
   lon_index  = find_name('lon',dim_names)
   lat_index  = find_name('lat',dim_names)
   slon_index = find_name('slon',dim_names)
   slat_index = find_name('slat',dim_names)
   lev_index  = find_name('lev',dim_names)

   fld_index_ps = find_name('PS',state_names_2d)
   dim2 = s_dim_2d(2,fld_index_ps)
   dim1 = s_dim_2d(1,fld_index_ps)
   dim_lev = dim_sizes(lev_index)

! newFIXME: change p and ps to mean_p mean_ps and ...  throughout
   if (allocate_ps) then
      allocate(ps(         dim_sizes(lon_index), dim_sizes(lat_index)))
      allocate(p (dim_lev, dim_sizes(lon_index), dim_sizes(lat_index)))
      if (slon_index /= 0) &
         allocate(ps_stagr_lon (dim_sizes(slon_index), dim_sizes( lat_index)))
      if (slat_index /= 0) &
         allocate(ps_stagr_lat (dim_sizes( lon_index), dim_sizes(slat_index)))
      allocate_ps = .false.
   endif

   fld_index = find_name('PS',cflds)
   ind = index_from_grid(1,1,1,fld_index)
   ps = reshape(vec(ind:ind+(dim1*dim2)-1),(/dim1,dim2/))

   ! Fill the p(:,:,:) array.
   do n=1,dim2
   do m=1,dim1
      call plevs_cam(ps(m,n), dim_lev, p(1:dim_lev,m,n))
   enddo
   enddo
   ! write(*,'(A,1pe12.4)') 'p(1,      1,   1) = ',p(1,      1,   1)
   ! write(*,'(A,1pe12.4)') 'p(dim_lev,1,   1) = ',p(dim_lev,1,   1)
   ! write(*,'(A,1pe12.4)') 'p(1,      dim1,1) = ',p(1,      dim1,1)
   ! write(*,'(A,1pe12.4)') 'p(dim_lev,dim1,1) = ',p(dim_lev,dim1,1)
   ! write(*,'(A,1pe12.4)') 'p(1,      1,   dim2) = ',p(1,      1,   dim2)
   ! write(*,'(A,1pe12.4)') 'p(dim_lev,1,   dim2) = ',p(dim_lev,1,   dim2)
   ! write(*,'(A,1pe12.4)') 'p(1,      dim1,dim2) = ',p(1,      dim1,dim2)
   ! write(*,'(A,1pe12.4)') 'p(dim_lev,dim1,dim2) = ',p(dim_lev,dim1,dim2)


   if (slon_index /= 0) then
      do n=1,dim_sizes( lat_index)
         ! The first element of slon is 1/2 grid obs *west* of lon(1) (which
         ! = 0. longitude)
         ! Make p[hi]s_stagr_lon line up with the slon array.
         ! The index to the second ps can be either slon or lon because num_slon = num_lon
         ps_stagr_lon(1,n) = 0.5_r8 * (ps(1,n) + ps(dim_sizes(slon_index),n))
         do m=2,dim_sizes(slon_index)
            ps_stagr_lon(m,n) = 0.5_r8 * (ps(m-1,n) + ps(m,n))
         enddo
      enddo
   endif

   if (slat_index /= 0) then
      do n=1,dim_sizes(slat_index)
         do m=1,dim_sizes(lon_index)
            ps_stagr_lat(m,n) = 0.5_r8 * (ps(m,n) + ps(m,n+1))
         enddo
      enddo
   endif
else
   ! Non-rectangular grid; 1 horizontal dimension.

   !CS PS is a 1d field (cubed sphere).
   dim2    = 1
   fld_index_ps = find_name('PS',state_names_1d)
   dim1    = s_dim_1d(fld_index_ps)
   lev_index  = find_name('lev',dim_names)
   dim_lev = dim_sizes(lev_index)

   if (allocate_ps) then
      allocate(ps        (dim1, dim2))
      allocate(p(dim_lev, dim1, dim2))
      allocate_ps = .false.
   endif

   fld_index = find_name('PS',cflds)
   ! Index_from_grid returns first (1,1,1) index of field fld_index,
   ! no matter where it is in the state vector.
   ind       = index_from_grid(1,1,1,fld_index) -1
   do m=1,dim1
      ind = ind + 1
      ps(m,1) = vec(ind)
      call plevs_cam(ps(m,1), dim_lev, p(1:dim_lev,m,1))
   enddo

   if (output_task0) then
      write(string1,*) 'Finished assignment of PS for ',dim1,' elements.  Will write ',dim_lev,' levels'
      call error_handler(E_MSG, 'set_ps_arrays', string1,source,revision,revdate)
      do n=1,dim_lev
         write(string1,'(3X,(1p8E12.4))') (p(n,m,1),m=1,dim1,dim1/7)
         ! Slow1 This is NOT written out
         call error_handler(E_MSG, 'set_ps_arrays', string1,source,revision,revdate)
      enddo
   endif
endif

end subroutine set_ps_arrays

!-----------------------------------------------------------------------

function scale_height(p_surface, p_above)

! Function to calculate scale height, given a surface pressure and a pressure.
! Using the surface pressure instead of, e.g., mean sea level as the reference pressure
! ensures that scale height is always positive.
! FIXME; But is it a distortion to have the scale heights follow the terrain?

real(r8), intent(in) :: p_surface
real(r8), intent(in) :: p_above
real(r8)             :: scale_height

scale_height = 5000.0_r8  ! arbitrary impossibly large number of scale heights.
if (p_above > 0.0_r8) scale_height = log(p_surface/p_above)

end function scale_height

!-----------------------------------------------------------------------

subroutine plevs_cam (p_surf, num_levs, pmid )

! Define the pressures of the layer midpoints from the
! coordinate definitions and the surface pressure.

real(r8), intent(in)  :: p_surf    ! Surface pressure (pascals)
integer,  intent(in)  :: num_levs
real(r8), intent(out) :: pmid(:)   ! Pressure at model levels

integer :: k

! Set midpoint pressures and layer thicknesses

do k=1,num_levs
   pmid(k) = hyam%vals(k)*P0%vals(1) + hybm%vals(k)*p_surf
enddo

end subroutine plevs_cam

!-----------------------------------------------------------------------

subroutine plevs_cam_distrib(p_surf, num_levs, pmid, ens_size)

! Define the pressures of the layer midpoints from the
! coordinate definitions and the surface pressure.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: p_surf(ens_size)    ! Surface pressure (pascals)
integer,  intent(in)  :: num_levs
real(r8), intent(out) :: pmid(ens_size, num_levs)   ! Pressure at model levels

integer :: k

! Set midpoint pressures and layer thicknesses

do k=1,num_levs
   pmid(:, k) = hyam%vals(k)*P0%vals(1) + hybm%vals(k)*p_surf
enddo

end subroutine plevs_cam_distrib

!-----------------------------------------------------------------------
!HK do you need a forward and a mean version of this?
! I think you need both.
subroutine model_heights_distrib_fwd(num_levs, state_ens_handle, p_surf, base_obs_loc, model_h, istatus)

! This routine calculates geometrical height (m) at mid-layers of the CAM model
!
! was Hui's dcz2ccm1
!    has globally defined inputs:
!          hyam(num_levs),hybm(num_levs),hyai(num_levs),hybi(num_levs) =
!          hybrid vertical coefficients, top to bottom.
!          (P = P0*hyam + ps*hybm)
!          P0 - Hybrid base pressure (pascals)
!
! Kevin Raeder converted to single column version 4/28/2006
!              removed longitude dimension entirely and extra arrays 10/2006
!   5/31/2013; Rewritten to adapt to convert_vert handling obs TYPEs,
!              not obs KINDS, and to handle lonlat and cubed sphere
!              grids/interpolations.

integer,             intent(in) :: num_levs
type(ensemble_type), intent(in) :: state_ens_handle
real(r8),            intent(in) :: p_surf(:) ! ens_size
type(location_type), intent(in) :: base_obs_loc

real(r8), intent(out) :: model_h(:, :) ! HK This is a global, why are you passing it?
integer,  intent(out) :: istatus

! local variables; pterm must be dimensioned as an array because dcz2 has it that way
real(r8), dimension(num_levs) :: phi, tv, q, t, pterm
real(r8) :: pmln(num_levs+1), hybrid_As(num_levs+1,2), hybrid_Bs(num_levs+1,2)
real(r8) :: phi_surf, ht_tmp
real(r8) :: l               ! location of ob in unit square space.
real(r8) :: m               ! location of ob in unit square space.
integer  :: cell_corners(4) ! corners of the cell which contains the ob

! CS Should these come from common_mod?
! That might be inconsistent with how levels, etc were defined in CAM originally.
! DART's values are 287.0_r8 and 461.6_r8.
real(r8), parameter :: rd = 287.05_r8
real(r8), parameter :: rv = 461.51_r8
real(r8), parameter :: rr_factor = (rv/rd) - 1.0_r8

real(r8) :: lon_lat_lev(3)
type(location_type) :: temp_obs_loc

integer :: k, i, e
integer, allocatable :: track_status(:), vstatus(:)

ens_size = copies_in_window(state_ens_handle)
allocate(track_status(ens_size, vstate(ens_size)))

istatus = 1
vstatus = 1

! Don't initialize these in the declaration statements, or they'll get the save attribute
! and won't be initialized each time this routine is entered.
l = MISSING_R8              ! location of ob in unit square space.
m = MISSING_R8              ! location of ob in unit square space.
cell_corners = MISSING_I    ! corners of the cell which contains the ob

model_h(:) = MISSING_R8
phi(:)     = MISSING_R8
pterm(:)   = MISSING_R8

! lat, lon and vertical in height
lon_lat_lev = get_location(base_obs_loc)

! copy to temporary arrays

! All arrays except hybrid_As, hybrid_Bs are oriented top to bottom.

! The 'interface' levels have an 'extra' level at model bottom, compared to the midpoint levels.
! Initialize this extra level, before filling the rest in a loop.
k = num_levs +1
hybrid_As(1,1) = hyai%vals(k)
hybrid_Bs(1,1) = hybi%vals(k)

!   hyam(num_levs) = 0 -> hybrid_As(2,2) = 0, so it
!   would be safe to set  hybrid_As(1,2) = 0.
!   It's safe because this element is used to set pmln in dcz2, but that element of pmln is never used.
hybrid_As(1,2) = 0.0_r8

! hyb[im]  have non-0 values at the bottom, 0s at the top;
!      hyb[im] coeffs multiply sigma in the calculation of pressure on levels,
!      and CAM's vertical coord is pure sigma at the bottom, so hybrid_Bs = 1.0 there.
hybrid_Bs(1,2) = 1.0_r8

! mid-points: 2nd dimension of hybrid_[AB]s = 2
! note that hyXm(num_levs + 1) is not defined (= MISSING_R8)
do k = 2,num_levs +1
   i = num_levs +2 - k
   hybrid_As(k,1) = hyai%vals(i)
   hybrid_Bs(k,1) = hybi%vals(i)
   hybrid_As(k,2) = hyam%vals(i)
   hybrid_Bs(k,2) = hybm%vals(i)
enddo

! Calculate phi_surf and tv for this column, for use by dcz2.
! HK This seems like way too much work. Same horizonal location for each level.
! HK What happened to get_interp_prof?
if (l_rectang) then

   call interp_lonlat(state_ens_handle, base_obs_loc, KIND_SURFACE_ELEVATION, vstatus,phi_surf)
   ! newFIXME; put Fail message other places like this   ! Failure; istatus = 1
   if (vstatus == 1) then
      write(string1,'(A,1p3F12.6)') 'surface elevation could not be interpolated in interp_lonlat at ', &
                           lon_lat_lev
      call error_handler(E_WARN, 'model_heights', string1)
      return
   endif

   ! loop through all levels to get the temperature and moisture.
   ! the interp routine will return a vstatus of 2 when the level is
   ! above the 'highest obs' threshold but we don't care about that here.
   ! error out for other return code values, but continue if vstatus is
   ! either 0 (all ok) or 2 (too high)
   do k = 1, num_levs
      ! construct a location with the same lat/lon but cycle though the model levels
      ! HK if they are on levels why not just get the columns and use the horizontal 
      ! interpolation?
      temp_obs_loc = set_location(lon_lat_lev(1), lon_lat_lev(2), real(k,r8), VERTISLEVEL)

      call interp_lonlat(state_ens_handle, temp_obs_loc, KIND_TEMPERATURE, vstatus, t(k))
      track_status = vstatus

      do e = 1, ens_size
         if (vstatus(e) == 1) then
            write(string1,'(A,I2,A)') 'Temperature level ',k, &
               ' could not be interpolated in interp_lonlat'
            call error_handler(E_WARN, 'model_heights', string1)
            !return
         endif
      enddo

      call interp_lonlat(state_ens_handle, temp_obs_loc, KIND_SPECIFIC_HUMIDITY, vstatus, q(k))

      do e = 1, ens_size
         if (vstatus(e) == 1 ) then
            write(string1,'(A,I2,A)') 'specific humidity level ',k, &
               ' could not be interpolated in interp_lonlat'
            call error_handler(E_WARN, 'model_heights', string1)
            !return
         endif
         if (vstatus(e) /= 0) track_status(e) = vstatus(e)
      enddo

      do e = 1, ens_size
         if (track_status(e) ==0) tv(k) = t(k)*(1.0_r8 + rr_factor*q(k))
      enddo
   enddo

else ! for cubed sphere:

   ! FIXME; It would be nice if the unit square locations (= weighting functions)
   !        found in this call could be used in the loop over levels, below.
   !        It can now, but addition of "cell_corners, l, m" to call needs to be tested.
   call error_handler(E_ERR, 'no RMA cubed sphere')

endif

do e = 1, ens_size
   call dcz2(num_levs, p_surf, phi_surf, tv, P0%vals(1) ,hybrid_As, hybrid_Bs, pmln, pterm, phi)

   ! used; hybrid_Bs, hybrid_As, hprb
   ! output from dcz2;  pmln, pterm , phi

   ! Conversion from geopotential height to geometric height depends on latitude
   ! Convert to kilometers for gph2gmh call, then back to meters for return value.
   do k = 1,num_levs
      ht_tmp = phi(k) * 0.001_r8        ! convert to km for following call only
      model_h(k, e) = gph2gmh(ht_tmp, lon_lat_lev(2)) * 1000.0_r8
   enddo
enddo

! model_heights returns only istatus 0 or 1
!istatus = 0 !HK This is annoying.
istatus = track_status

allocate(track_status, vstatus)

end subroutine  model_heights_fwd


!-----------------------------------------------------------------------
! HK Once again, model_h is global and used as an arguement.
! Why is this routine using interp_lonlat?  It did not before.
subroutine model_heights_distrib_mean(num_levs, state_ens_handle, p_surf, base_obs_loc, model_h, istatus)

! This routine calculates geometrical height (m) at mid-layers of the CAM model
!
! was Hui's dcz2ccm1
!    has globally defined inputs:
!          hyam(num_levs),hybm(num_levs),hyai(num_levs),hybi(num_levs) =
!          hybrid vertical coefficients, top to bottom.
!          (P = P0*hyam + ps*hybm)
!          P0 - Hybrid base pressure (pascals)
!
! Kevin Raeder converted to single column version 4/28/2006
!              removed longitude dimension entirely and extra arrays 10/2006
!   5/31/2013; Rewritten to adapt to convert_vert handling obs TYPEs,
!              not obs KINDS, and to handle lonlat and cubed sphere
!              grids/interpolations.

integer,             intent(in) :: num_levs
type(ensemble_type), intent(in) :: state_ens_handle
real(r8),            intent(in) :: p_surf
type(location_type), intent(in) :: base_obs_loc

real(r8), intent(out) :: model_h(:)
integer,  intent(out) :: istatus

! local variables; pterm must be dimensioned as an array because dcz2 has it that way
real(r8), dimension(num_levs) :: phi, tv, q, t, pterm
real(r8) :: pmln(num_levs+1), hybrid_As(num_levs+1,2), hybrid_Bs(num_levs+1,2)
real(r8) :: phi_surf, ht_tmp
real(r8) :: l               ! location of ob in unit square space.
real(r8) :: m               ! location of ob in unit square space.
integer  :: cell_corners(4) ! corners of the cell which contains the ob

! CS Should these come from common_mod?
! That might be inconsistent with how levels, etc were defined in CAM originally.
! DART's values are 287.0_r8 and 461.6_r8.
real(r8), parameter :: rd = 287.05_r8
real(r8), parameter :: rv = 461.51_r8
real(r8), parameter :: rr_factor = (rv/rd) - 1.0_r8

real(r8) :: lon_lat_lev(3)
type(location_type) :: temp_obs_loc

integer :: k, i, vstatus

istatus = 1
vstatus = 1

! Don't initialize these in the declaration statements, or they'll get the save attribute
! and won't be initialized each time this routine is entered.
l = MISSING_R8              ! location of ob in unit square space.
m = MISSING_R8              ! location of ob in unit square space.
cell_corners = MISSING_I    ! corners of the cell which contains the ob

model_h(:) = MISSING_R8
phi(:)     = MISSING_R8
pterm(:)   = MISSING_R8

! lat, lon and vertical in height
lon_lat_lev = get_location(base_obs_loc)

! copy to temporary arrays

! All arrays except hybrid_As, hybrid_Bs are oriented top to bottom.

! The 'interface' levels have an 'extra' level at model bottom, compared to the midpoint levels.
! Initialize this extra level, before filling the rest in a loop.
k = num_levs +1
hybrid_As(1,1) = hyai%vals(k)
hybrid_Bs(1,1) = hybi%vals(k)

!   hyam(num_levs) = 0 -> hybrid_As(2,2) = 0, so it
!   would be safe to set  hybrid_As(1,2) = 0.
!   It's safe because this element is used to set pmln in dcz2, but that element of pmln is never used.
hybrid_As(1,2) = 0.0_r8

! hyb[im]  have non-0 values at the bottom, 0s at the top;
!      hyb[im] coeffs multiply sigma in the calculation of pressure on levels,
!      and CAM's vertical coord is pure sigma at the bottom, so hybrid_Bs = 1.0 there.
hybrid_Bs(1,2) = 1.0_r8

! mid-points: 2nd dimension of hybrid_[AB]s = 2
! note that hyXm(num_levs + 1) is not defined (= MISSING_R8)
do k = 2,num_levs +1
   i = num_levs +2 - k
   hybrid_As(k,1) = hyai%vals(i)
   hybrid_Bs(k,1) = hybi%vals(i)
   hybrid_As(k,2) = hyam%vals(i)
   hybrid_Bs(k,2) = hybm%vals(i)
enddo

! Calculate phi_surf and tv for this column, for use by dcz2.
if (l_rectang) then

   ! HK why do you need interp_lonlat? You did not it previous versions of this routine.
   call interp_lonlat_distrib(state_ens_handle, base_obs_loc, KIND_SURFACE_ELEVATION, vstatus, phi_surf)
   ! newFIXME; put Fail message other places like this   ! Failure; istatus = 1
   if (vstatus == 1) then
      write(string1,'(A,1p3F12.6)') 'surface elevation could not be interpolated in interp_lonlat at ', &
                           lon_lat_lev
      call error_handler(E_WARN, 'model_heights', string1)
      return
   endif

   ! loop through all levels to get the temperature and moisture.
   ! the interp routine will return a vstatus of 2 when the level is
   ! above the 'highest obs' threshold but we don't care about that here.
   ! error out for other return code values, but continue if vstatus is
   ! either 0 (all ok) or 2 (too high)
   do k = 1, num_levs
      ! construct a location with the same lat/lon but cycle though the model levels
      temp_obs_loc = set_location(lon_lat_lev(1), lon_lat_lev(2), real(k,r8), VERTISLEVEL)

      call interp_lonlat_distrib(state_ens_handle, temp_obs_loc, KIND_TEMPERATURE, vstatus, t(k))
      if (vstatus == 1) then
         write(string1,'(A,I2,A)') 'Temperature level ',k, &
              ' could not be interpolated in interp_lonlat'
         call error_handler(E_WARN, 'model_heights', string1)
         return
      endif
      call interp_lonlat_distrib(state_ens_handle, temp_obs_loc, KIND_SPECIFIC_HUMIDITY, vstatus, q(k))
      if (vstatus == 1 ) then
         write(string1,'(A,I2,A)') 'specific humidity level ',k, &
              ' could not be interpolated in interp_lonlat'
         call error_handler(E_WARN, 'model_heights', string1)
         return
      endif

      tv(k) = t(k)*(1.0_r8 + rr_factor*q(k))
   enddo

else ! for cubed sphere:

   ! FIXME; It would be nice if the unit square locations (= weighting functions)
   !        found in this call could be used in the loop over levels, below.
   !        It can now, but addition of "cell_corners, l, m" to call needs to be tested.
   call error_handler(E_ERR, 'no RMA cubed sphere')

endif

call dcz2(num_levs, p_surf, phi_surf, tv, P0%vals(1) ,hybrid_As, hybrid_Bs, pmln, pterm, phi)

! used; hybrid_Bs, hybrid_As, hprb
! output from dcz2;  pmln, pterm , phi

! Conversion from geopotential height to geometric height depends on latitude
! Convert to kilometers for gph2gmh call, then back to meters for return value.
do k = 1,num_levs
   ht_tmp = phi(k) * 0.001_r8        ! convert to km for following call only
   model_h(k) = gph2gmh(ht_tmp, lon_lat_lev(2)) * 1000.0_r8
enddo

! model_heights returns only istatus 0 or 1
istatus = 0

end subroutine  model_heights_distrib_mean


!-----------------------------------------------------------------------

subroutine dcz2(kmax,p_surf,phis0,tv,hprb,hybrid_As,hybrid_Bs,pmln,pterm,z2)

! Compute geopotential height for a CESM hybrid coordinate column.
! All arrays except hybrid_As, hybrid_Bs are oriented top to bottom.
! hybrid_[AB]s first subscript:
!  = 1 for layer interfaces
!  = 2 for layer midpoints
! hybrid_As coord coeffs for P0 reference pressure term in plevs_cam
! hybrid_Bs coord coeffs for surf pressure term in plevs_cam (in same format as hybrid_As)

integer,  intent(in)  :: kmax                ! Number of vertical levels
real(r8), intent(in)  :: p_surf              ! Surface pressure           (pascals)
real(r8), intent(in)  :: phis0               ! Surface geopotential
real(r8), intent(in)  :: tv(kmax)            ! Virtual temperature, top to bottom
real(r8), intent(in)  :: hprb                ! Hybrid base pressure       (pascals)
real(r8), intent(in)  :: hybrid_As(kmax+1,2)
real(r8), intent(in)  :: hybrid_Bs(kmax+1,2)
real(r8), intent(out) :: pmln(kmax+1)        ! logs of midpoint pressures
real(r8), intent(out) :: pterm(kmax)         ! pressure profile
real(r8), intent(out) :: z2(kmax)            ! Geopotential height, top to bottom

! Local variables
real(r8), parameter :: r = 287.04_r8    ! Different than model_heights !
real(r8), parameter :: g0 = 9.80616_r8  ! Different than model_heights:gph2gmh:G !
real(r8), parameter :: rbyg=r/g0

integer  :: i,k,l
real(r8) :: ARG

! Compute intermediate quantities using scratch space

! Invert vertical loop
! Compute top only if top interface pressure is nonzero.
!
! newFIXME; p_col could be used here, instead of (re)calculating it in ARG
do K = kmax+1, 1, -1
   i = kmax-K+2
   ARG = hprb*hybrid_As(i,2) + p_surf *hybrid_Bs(i,2)
   if (ARG > 0.0_r8) THEN
       pmln(K) = LOG(ARG)
   else
       pmln(K) = 0.0_r8
   endif
enddo

do K = 2,kmax - 1
   pterm(k) = rbyg*tv(k)*0.5_r8* (pmln(k+1)-pmln(k-1))
enddo

! Initialize z2 to sum of ground height and thickness of top half-layer
do K = 1,kmax - 1
   z2(k) = phis0/g0 + rbyg*tv(k)*0.5_r8* (pmln(K+1)-pmln(K))
enddo
z2(kmax) = phis0/g0 + rbyg*tv(kmax)* (log(p_surf*hybrid_Bs(1,1))-pmln(kmax))

do k = 1,kmax - 1
    z2(k) = z2(k) + rbyg*tv(kmax)* (log(p_surf*hybrid_Bs(1,1))-0.5_r8* &
                                       (pmln(kmax-1)+pmln(kmax)))
enddo

! Add thickness of the remaining full layers
! (i.e., integrate from ground to highest layer interface)

do K = 1,kmax - 2
    do L = K+1, kmax-1
       z2(K) = z2(K) + pterm(L)
    enddo
enddo

end subroutine dcz2

!-----------------------------------------------------------------------

function gph2gmh(h, lat)

!  Convert a list of geopotential altitudes to mean sea level altitude.

real(r8), intent(in) :: h         ! geopotential altitude (in km)
real(r8), intent(in) :: lat       ! latitude  of profile in degrees.
real(r8)             :: gph2gmh   ! MSL altitude, in km.

real(r8), parameter ::  be = 6356.7516_r8             ! min earth radius, km
real(r8), parameter ::  ae = 6378.1363_r8             ! max earth radius, km
real(r8), parameter ::  pi = 3.14159265358979_r8
! FIXME; another definition of gravitational acceleration.  See g0 and gravity_constant elsewhere.
real(r8), parameter ::  G = 0.00980665_r8          ! WMO reference g value, km/s**2, at 45.542N(S)

real(r8) :: g0
real(r8) :: r0
real(r8) :: latr

latr = lat * (pi/180.0_r8)           ! in radians
call gravity(latr, 0.0_r8, g0)

! compute local earth's radius using ellipse equation

r0 = sqrt( ae**2 * cos(latr)**2 + be**2 * sin(latr)**2)

! Compute altitude above sea level
gph2gmh = (r0 * h) / (((g0*r0)/G) - h)

end function gph2gmh

!-----------------------------------------------------------------------

subroutine gravity(xlat,alt,galt)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! This subroutine computes the Earth's gravity at any altitude
! and latitude.  The model assumes the Earth is an oblate
! spheriod rotating at a the Earth's spin rate.  The model
! was taken from "Geophysical Geodesy, Kurt Lambeck, 1988".
!
!  input:    xlat, latitude in radians
!            alt,  altitude above the reference ellipsiod, km
!  output:   galt, gravity at the given lat and alt, km/sec**2
!
! Compute acceleration due to the Earth's gravity at any latitude/altitude
! author     Bill Schreiner   5/95
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real(r8), intent(in)  :: xlat
real(r8), intent(in)  :: alt
real(r8), intent(out) :: galt

real(r8),parameter :: xmu = 398600.4415_r8         ! km^3/s^2
real(r8),parameter :: ae  = 6378.1363_r8           ! km
real(r8),parameter :: f   = 1.0_r8/298.2564_r8
real(r8),parameter :: w   = 7.292115e-05_r8        ! rad/s
real(r8),parameter :: xm  = 0.003468_r8            !
real(r8),parameter :: f2  = 5.3481622134089e-03_r8 ! f2 = -f + 5.0* 0.50*xm - 17.0/14.0*f*xm + 15.0/4.0*xm**2
real(r8),parameter :: f4  = 2.3448248012911e-05_r8 ! f4 = -f**2* 0.50 + 5.0* 0.50*f*xm

real(r8) :: ge
real(r8) :: g


! compute gravity at the equator, km/s2
ge = xmu/ae**2/(1.0_r8 - f + 1.5_r8*xm - 15.0_r8/14.0_r8*xm*f)

! compute gravity at any latitude, km/s2
g = ge*(1.0_r8 + f2*(sin(xlat))**2 - 1.0_r8/4.0_r8*f4*(sin(2.0_r8*xlat))**2)

! compute gravity at any latitude and at any height, km/s2
galt = g - 2.0_r8*ge*alt/ae*(1.0_r8 + f + xm + (-3.0_r8*f + 5.0_r8* 0.50_r8*xm)*  &
                          (sin(xlat))**2) + 3.0_r8*ge*alt**2/ae**2

end subroutine gravity

!-----------------------------------------------------------------------

subroutine init_model_instance(var)

! Initializes an instance of a cam model state variable

type(model_type), intent(inout) :: var

if (.not. module_initialized) call static_init_model()

! Initialize the storage space and return

! The temporary arrays into which fields are read are dimensioned by the largest values of
! the sizes of the dimensions listed in s_dim_RANKd.  Those are stored in s_dim_max.

allocate(var%vars_0d(                                              state_num_0d))
allocate(var%vars_1d(s_dim_max(1,1),                               state_num_1d))
allocate(var%vars_2d(s_dim_max(1,2),s_dim_max(2,2),                state_num_2d))
allocate(var%vars_3d(s_dim_max(1,3),s_dim_max(2,3),s_dim_max(3,3), state_num_3d))

end subroutine init_model_instance

!-----------------------------------------------------------------------

subroutine end_model_instance(var)

! Ends an instance of a cam model state variable

type(model_type), intent(inout) :: var

if (.not. module_initialized) call static_init_model()

if (.not. allocated(var%vars_0d)) then
   write(string1,*) 'Calling end_model_instance on an uninitialized state structure'
   call error_handler(E_ERR,'end_model_instance',string1, source, revision, revdate)
endif

deallocate(var%vars_0d, var%vars_1d, var%vars_2d, var%vars_3d)

end subroutine end_model_instance


! End of utility routines

!#######################################################################

!-----------------------------------------------------------------------
!>
!> Subroutine adv_1step
!> advances model 1 forecast length  -- noop for CESM atmospheric components.
!> 
!> @param[inout] st_vec(:)
!> The state vector which is NOT advanced by this routine.
!> 
!> @param[in] Time
!> The DART time_type which is NOT the end of a forecast.


subroutine adv_1step(st_vec, Time)

real(r8), intent(inout) :: st_vec(:)

! Time is needed for more general models like this; need to add in to low-order models.
type(time_type), intent(in) :: Time

! This is a no-op for CAM; only asynch integration
! Can be used to test the assim capabilities with a null advance

if (.not. module_initialized) call static_init_model()

! make it an error by default; comment these calls out to actually
! test assimilations with null advance.

call error_handler(E_ERR,'adv_1step', &
                  'CAM model cannot be called as a subroutine; async cannot = 0', &
                  source, revision, revdate)
! newFIXME; add code here to silence compiler warnings about unused variables.

end subroutine adv_1step

!-----------------------------------------------------------------------
!>
!> Subroutine end_model
!> deallocates arrays that are in module global storage.

subroutine end_model()

! HK no ens_mean in distributed
deallocate(dim_names, dim_sizes, phis)
deallocate(TYPE_1D,TYPE_2D,TYPE_3D)
deallocate(state_long_names, state_units)
deallocate(cflds)

if (allocated(s_dim_3d)) then
   deallocate(s_dim_3d, s_dimid_3d, f_dim_3d, f_dimid_3d)
endif
if (allocated(s_dim_2d)) then
   deallocate(s_dim_2d, s_dimid_2d, f_dim_2d, f_dimid_2d)
endif
if (allocated(s_dim_1d)) then
   deallocate(s_dim_1d, s_dimid_1d, f_dim_1d, f_dimid_1d)
endif

if (allocated(phis_stagr_lon)) deallocate(phis_stagr_lon)
if (allocated(phis_stagr_lat)) deallocate(phis_stagr_lat)

deallocate (ps, p, p_col, model_h)
if (allocated(ps_stagr_lon)) deallocate(ps_stagr_lon)
if (allocated(ps_stagr_lat)) deallocate(ps_stagr_lat)

if (.not. l_rectang) then
   deallocate(cs_locs)
   deallocate(cs_kinds)
   deallocate(cs_locs_xyz)
   deallocate(lon_rad, lat_rad)
   deallocate(corners)
   deallocate(num_nghbrs, centers, a, b, x_ax_bearings)
endif

call end_grid_1d_instance(lon)
call end_grid_1d_instance(lat)
call end_grid_1d_instance(lev)
call end_grid_1d_instance(gw)
call end_grid_1d_instance(hyam)
call end_grid_1d_instance(hybm)
call end_grid_1d_instance(hyai)
call end_grid_1d_instance(hybi)
call end_grid_1d_instance(slon)
call end_grid_1d_instance(slat)
call end_grid_1d_instance(ilev)
call end_grid_1d_instance(P0)

! Deallocate _gc variables; cs_gc_xyz and cs_gc
call finalize_closest_node()
call get_close_obs_destroy(cs_gc)

end subroutine end_model

!-----------------------------------------------------------------------
!>
!> Subroutine init_time
!> reads in initial time  -- noop for CESM atmospheric components.
!> 
!> @param[inout] time
!> The DART time_type time which is NOT initialized here.

subroutine init_time(time)

! For now returns value of Time_init which is set in initialization routines.

type(time_type), intent(inout) :: time

if (.not. module_initialized) call static_init_model()

call error_handler(E_ERR,"init_conditions", &
                  "WARNING!!  CAM model has no built-in default time", &
                  source, revision, revdate, &
                  text2="cannot run with 'start_from_restart = .false.'", &
                  text3="use 'cam_to_dart' to create a CAM state vector file")

! To silence the compiler warnings:
time = set_time(0, 0)

end subroutine init_time

!-----------------------------------------------------------------------

subroutine set_print_details(how)

! reset the print_details module global variable to control
! how much output there is

logical, intent(in) :: how

print_details = how

end subroutine set_print_details

!#######################################################################
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$


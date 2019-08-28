! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! This was copied from model_mod_chem_xyz.f90,
! which was tested in SE30_{Og16_4mem_4,Og16_osse5,loc_test*,Hyb*}.
! Earlier versions of model_mod with the algorithms for finding the
!    cell which contains an ob were tested in check_model_mod.f90,
!    and print out from the search process confirmed that these algorithms
!    find the right cell.
!    That test hasn't been repeated for this model_mod, but will be soon (12/7/2013).
! It is believed to have all of the vert_coord changes to accomodate WACCM,
!    and chemistry species to enable CAMchem assimilations.
!    Nick Pedatella tested it for the FV core version of WACCM.
! It also has the cartesian 'xyz' search for the nearest node (only 1) in get_close_obs.
! And, of course, it handles CAM-SE state vectors.

module model_mod

!----------------------------------------------------------------------
! Interface code between CAM and DART.  Contains the required 16 interfaces
!  as specified by DART.  Also contains several utility routines which help
!  translate between CAM and DART formats, and deal with time.
!
!  Contains a perturb routine for generating initial ensembles.  Does not
!  provide adv_1step or init_conditions because CAM is a separate executable
!  and cannot be called as a subroutine.
!
!  This module intercepts the get_close_obs() calls and alters the distances
!  for obs at the top of the model so they do not impact the state.
!
!  This module keeps a copy of the ensemble mean in module global storage and
!  uses it for computing the pressure-to-height conversions.
!
!  See the subversion code logs for history of this module.
!
!  During the assimilation stage, only a piece of the state vector is available to each
!  process, and each process calls parts of model_mod.  In order to handle the conversion
!  of vertical coordinates of obs and/or state variables into a consistent coordinate,
!  an entire state vector is needed, so the ensemble mean is passed to model_mod before
!  the assimilation starts.  This is NOT done for model_interpolate; the whole vector is
!  available, and should be used.  All locations are now converted to a standard coordinate
!  (pressure or log(pressure)), instead of always converting the state vertical location to
!  that of the ob.  The highest_obs_level and ..._height_m variables are derived from
!  highest_obs_pressure_Pa namelist variable.
!
!  This module has been updated to handle both the eulerian and finite volume core versions
!  of CAM (they have different logically rectangular grids),
!  and the HOMME dynamical core of Spectral Element CAM which uses the cubed sphere
!  (non-rectangular) horizontal grid.
!
!  The coordinate orders of fields stored in various forms have also been simplified.
!  For example; various vintages of CAM 3D fields may be read in with (lon, lat, lev) or
!  (lon, lev, lat).  These are uniformly converted to (lev, lon, lat) for use in model_mod.
!  This latter form is different than pre MPI model_mods.  Then such fields are stored in
!  the state vector with the same coordinate order.  They are converted back to the modern
!  CAM coordinate order when written to caminput.nc files.
!
!     These may be needed on the regular A-grid (thermodynamic variables) and grids staggered
!     relative to the A-grid.   Currently, PS for the A-grid and for the 2 staggered grids is
!     stored for global access for the (re)calculation of pressures and heights on model levels
!     as needed.  In the future it may be deemed worthwhile to store the 3d pressures and
!     heights on the 3 grids, but for now that seemed like too much memory to be worthwhile.
!
!     If a user wants to add new TYPE_s to the state vector,
!     then more QTY_s may be needed in the obs_kind_mod and the 'use obs_kind_mod' statement.

!     The coordinates of CAM (lats, lons, ncol, etc.) and their dimensions  and attributes are
!     now read into globally accessible data structures (see grid_1d_type).
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

!==============================================================================================
! TO DO
!--------------------------------------------------------------------
! DISTRIBUTION;
! > instructions for use
! > new filter_ics for the cases I've distributed
! > init2ud script for them to create their own.
!      list of files that script/model_mod needs;
!          cam_to_dart+restart_file_tool (was trans_pv_sv_time0)  (the executable)
!          input.nml:perfect_model_nml (with the date+time to put in the filter_ics files)
!          cam_phis.nc  (for static_init_model to get PHIS)
!          caminput_#.nc  (to convert into filter_ic.####)

!--------------------------------------------------------------------
!
! Describe:
!    - The coordinate orders and translations; CAM initial file, model_mod, and DART _Diag.nc.
!      Motivations
!    - There need to be 2 sets of arrays for dimensions and dimids;
!        one describing the caminput file (f_...)
!        and one for the state (s_...) (storage in this module).
!             Call them f_dim_Nd , f_dimid_Nd
!                       s_dim_Nd , s_dimid_Nd
!

!   0   > Convert read_cam_2Dreal too?   convert phis to a grid_2d_type?
!    x  > See 'premature' below for something else to fix
!       > change (private only) subroutine argument lists; structures first, regardless of in/out
!         then output, and input variables.
!       > change declarations to have dummy argument integers used as dimensions first
!       > deallocate grid_1d_arrays using end_1d_grid_instance in end_model
!       > test with a single ob (to see how get_state_meta_data indexing works)
!    (x) > document before compiling (and forgetting!)
!          - model_interpolate (and called routines?)
!       > new qc; are my istatus = [012] still good?
!       > More error checking?
!       > debug statements?  Test with a small obs set, small resolution, (limit further?)
!       > recollect lists of subroutines & functions
!           document them in web page (and here?)

! EFFICIENCY   (see model_mod.f90 from 2/1/07 for speed-ups that were implemented already)
!       + index_from_grid (and others?) could be more efficient by calculating and
!         globally storing the beginning index of each cfld and/or the size of each cfld.
!       + global storage of height fields?  but need them on staggered grids (only sometimes)
!       ! Some compilers can't handle passing a section of an array to a subroutine/function;
!         I do this in nc_write_model_vars(?) and/or write_cam_init(?); replace with an
!         exactly sized array?
!     x + replace lon (and lat) index finders in model_interpolate with (new) coord_index.
!            But cubed sphere will need a different sort of search for location indices?
!         CS:
!            Convert lon,lat refs into dim1,dim2 in more subroutines.
!               get_val_heights is called with (column_ind,1) by CAM-SE code, and (lon_ind, lat_ind) otherwise).
!            Change units of lons, lats read from model_config_file from degrees to radians?
!               Helpful for calls to bearing, but what about elsewhere?
!               Model_heights may be expecting it in degrees.
!               > Test what I have, then maybe change later when looking for speed-up.
!            Remove ability to handle (lon,lev,lat) coordinate order of old CAMs.
!            Add powerpoint of finding cell algorithm to the models/cam/doc directory.
!               First add slides about the 'wrong cell' solution,
!               and the interpolation.
!            Replace some do loops with forall (constructs)
!          x lfound  initialized, scoped correctly?
!          x Comment/remove all these and embedded notes and questions.
!            Add code to scripts and namelists to locate a h0 file and link to it
!              in read_cam_2Dreal, in case cam_phis.nc is not found.
!            subroutine write_cam_times(model_time, adv_time)
!               Not needed in CESM+DART framework?
!            Remove the code that accomodates old CAM coordinate order (lon,lev,lat).




! See model_mod.f90 from ~2/1/07 for more ISSUES and resolutions
!
! ISSUE; model_heights was hardwired to use ens_mean, but 2 different uses of
!        model_heights passed ens_mean and st_vec to it.  model_heights now uses the
!        array passed into 'vec'.  Confirm that this is correct.

! ISSUE; 3d P calculation assumes that PS is the first (spatially) 2D field in
!        the state_names list and in the state vector.
!        It might not be in either, or be in the wrong place.

! ISSUE; In P[oste]rior_Diag.nc ensemble members are written out *between* the field mean/spread
!        pair and the inflation mean/sd pair.  Would it make more sense to put members after
!        both pairs?  Easy to do?

! ISSUE?; model_interpolate assumes that obs with a vertical location have 2 horizontal locations
!          too.  The state vector may have fields for which this isn't true, but no obs we've seen
!          so far violate this assumption.  It would have to be a synthetic/perfect_model obs,
!          like some sort of average or parameter value.

! ISSUE; In convert_vert, if a 2D field has dimensions (lev, lat) then how is p_surf defined?
!        I've set it to P0, but is this correct or meaningful?

! ISSUE; The QTY_ list from obs_def_mod must be updated when new fields are added to state vector.
!        This could be done by the preprocessor when it inserts the code bits corresponding to the
!        lists of observation types, but it currently (10/06) does not.  Document accordingly.

! ISSUE: Should get_val_xxxx return MISSING_VALUE instead of 0 for istat = 1 cases?

! ISSUE: The CCM code (and Hui's packaging) for geopotentials and heights  use different
!        values of the physical constants than DART's.  In one case Shea changed g from
!        9.81 to 9.80616, to get agreement with CCM(?...), so it may be important.
!        Also, matching with Hui's tests may require using his values;  change to DART
!        after verifying?

! ISSUE: It's possible to figure out the model_version from the NetCDF file
!        itself, rather than have that be user-provided (sometimes incorrect and hard
!        to debug) meta-data.  model_version is also misnamed; it's really the
!        caminput.nc model version.  The actual model might be a different version(?)
!        The problem with removing it from the namelist is that the scripts need it
!        too, so some rewriting there would be needed.

! "CS" marks changes to accommodate Spectral Element CAM's cubed sphere grid;
!      logically non-rectangular.
!      ==============================
!      New/Replaced/Renamed procedures:
!      read_topog_size ->                  folded into read_cam_2Dreal
!      read_cam_horiz  -> x x x subroutine read_cam_2Dreal (so far only for PHIS)
!      _______         -> x x x subroutine read_cam_2Dint  (so far only for HommeMapping.nc)
!                                 read_cam_2D might be made more general by:
!                                     splitting off the opening and dimension reading from read_cam_2Dreal
!                                     adding functionality to read_cam_2Dint?
!                                 Then they could be called by nc_write_cs_grid_file.
!      _______         -> x x x subroutine create_cs_grid_arrays (to generate and write out CS neighbor arrays)
!      _______         -> x x x function   nc_write_cs_grid_file
!      _______         -> x x x function   nc_read_cs_grid_file
!      ______          -> x x x function   bearing
!      _______         -> x x x subroutine nc_read_global_att
!      _______         -> x x x subroutine fill_gc
!      model_interpola -> x x x subroutine interp_lonlat
!      _______         -> x x x subroutine model_interpolate (just a fork to either interp_lonlat or interp_cubed_sphere)
!      _______         -> x x x subroutine interp_cubed_sphere
!      _______         -> x x x subroutine coord_ind_cs
!      _______         -> x x x subroutine unit_square_location
!      _______         -> x x x subroutine convert_vert (calls coord_ind_cs, and ...?)
!                         | | |
!                         | | careful read through and add debugging/evaluation output
!                         | compile
!                         test offline

! "Pobs" marks changes for providing expected obs of P
!        break from past philosophy; P is not a native CAM variable (but is already calced here)

! NOVERT marks modifications for fields with no vertical location,
! i.e. GWD parameters.



!==============================================================================================
!  USE statements

use netcdf
use typeSizes

use types_mod,         only : r8, MISSING_I, MISSING_R8, gravity_const => gravity, PI, DEG2RAD, RAD2DEG
!          add after verification against Hui's tests;  gas_constant_v,gas_constant,ps0,PI,DEG2RAD

use time_manager_mod,  only : time_type, set_time, set_date, print_time, print_date,  &
                              set_calendar_type, get_calendar_type, operator(-),      &
                              get_time, get_date
use utilities_mod,     only : open_file, close_file, find_namelist_in_file, check_namelist_read, &
                              register_module, error_handler, file_exist, E_ERR, E_WARN, E_MSG,  &
                              logfileunit, nmlfileunit, do_output, nc_check, get_unit, do_nml_file, &
                              do_nml_term
use mpi_utilities_mod, only : my_task_id, task_count

!-------------------------------------------------------------------------
use location_mod,      only : location_type, get_location, set_location, query_location,         &
                              LocationDims, LocationName, LocationLName, horiz_dist_only,        &
                              vert_is_level, vert_is_pressure, vert_is_height, vert_is_surface,  &
                              VERTISUNDEF, VERTISSURFACE, VERTISLEVEL,                           &
                              VERTISPRESSURE, VERTISHEIGHT,                                      &
                              get_close_type, get_close_maxdist_init, get_close_obs_init,        &
                              get_dist,loc_get_close_obs => get_close_obs

use xyz_location_mod, only : xyz_location_type, xyz_get_close_maxdist_init,    &
                             xyz_get_close_type, xyz_set_location, xyz_get_location, &
                             xyz_get_close_obs_init, xyz_get_close_obs_destroy, &
                             xyz_find_nearest

! get_close_maxdist_init, get_close_obs_init, can be modified here (i.e. to add vertical information
! to the initial distance calcs), but will need subroutine pointers like get_close_obs.
! READ THIS SYNTAX as:
!   There's a subroutine in location_mod named 'get_close_obs'.
!   If I want to use that one in this module then refer to it as 'loc_get_close_obs'.
!   If I call 'get_close_obs', then I'll get the one in this module,
!   which does some stuff I need, AND ALSO CALLS 'loc_get_close_obs'
!   ? How does filter know to call the get_close_obs in model_mod, instead of in location_mod?
!     Does the compiler figure out that any call to get_close_obs means the one in model mod?


!-----------------------------------------------------------------------------
! these PREPROCESS comment lines are not currently used, but are one
! proposed way to automatically extract the kinds needed by a model.
! the idea is that only those actually in use will be defined
! in obs_kind_mod.f90.
! BEGIN DART PREPROCESS USED KINDS
use     obs_kind_mod, only : QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, QTY_PRESSURE,     &
                             QTY_SURFACE_PRESSURE, QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY, &
                             QTY_CLOUD_LIQUID_WATER, QTY_CLOUD_ICE, QTY_CLOUD_FRACTION,    &
                             QTY_GRAV_WAVE_DRAG_EFFIC, QTY_GRAV_WAVE_STRESS_FRACTION,       &
                             QTY_SURFACE_ELEVATION,                                          &
                             QTY_CO, QTY_CO2, QTY_NO, QTY_NO2, QTY_CH4, QTY_NH3, QTY_O3, &
                             get_index_for_quantity, get_name_for_quantity, get_quantity_for_type_of_obs

! END DART PREPROCESS USED KINDS


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
!-----------------------------------------------------------------------------

use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, get_var_owner_index

use distributed_state_mod

! end of use statements
!=========================================================================================
!
! CAM global/module declarations

implicit none
private

! The first block are the 16 required interfaces.  The following block
! are additional useful interfaces that utility programs can call.
public ::                                                            &
   static_init_model, get_model_size, get_model_time_step,           &
   pert_model_state, get_state_meta_data_distrib, model_interpolate_distrib,         &
   nc_write_model_atts, nc_write_model_vars,                         &
   init_conditions, init_time, adv_1step, end_model,                 &
   get_close_maxdist_init, get_close_obs_init, get_close_obs_distrib,        &
   ens_mean_for_model, convert_base_obs_location

public ::                                                            &
   model_type, prog_var_to_vector, vector_to_prog_var,               &
   read_cam_init, read_cam_init_size,                                &
   init_model_instance, end_model_instance, write_cam_init,          &
   write_cam_times


interface get_surface_pressure
   module procedure get_surface_pressure_state, get_surface_pressure_mean
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
real(r8), allocatable :: ens_mean(:)

!----------------------------------------------------------------------
! Global storage for describing cam model class
!----------------------------------------------------------------------

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

! CAM3; additional types, or can TRACER handle the new ones?

!----------------------------------------------------------------------

! A type for cam model.
! Each variable will be allowed to have different dimensions, even different from
! others of the same rank (i.e. 2d).
! The maximum size for each dimension (for a given rank) will be used to allocate space
! when a model_type variable is initialized.
type model_type
    private
   real(r8), pointer :: vars_0d(:)
   real(r8), pointer :: vars_1d(:, :)
   real(r8), pointer :: vars_2d(:, :, :)
   real(r8), pointer :: vars_3d(:, :, :, :)
end type model_type

integer :: model_size
! This list of dimensions used to define fields will be ordered as they are on the caminput.nc file.
integer                                    :: num_dims
integer,                       allocatable :: dim_ids(:)
integer,                       allocatable :: dim_sizes(:)
character (len=NF90_MAX_NAME), allocatable :: dim_names(:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Grid fields
! These structures are used by nc_write_model_atts.
! They are dimensioned in init_grid_1D_instance and filled in read_cam_coord
type grid_1d_type
   private
   character (len=8)            :: label
   integer                      :: dim_id
   integer                      :: length
   real(r8)                     :: resolution
   real(r8), pointer            :: vals(:)
   integer                      :: num_atts
   character (len=128), pointer :: atts_names(:)
   character (len=128), pointer :: atts_vals(:)
end type grid_1d_type

! temporary output
integer :: num_calced = 0, num_searched = 0

integer :: iii
! integer :: grid_num_0d = 0              ! # of grid scalars to read from file
! P0 now a "coordinate",  and may be removed entirely
! character (len=8),dimension(100) :: grid_names_0d = (/'P0      ',('        ',iii=1,100)/)

integer                          :: grid_num_1d = 12   ! # of 1d grid fields to read from file
character (len=8),dimension(100) :: grid_names_1d = &
          (/ 'lon     ','lat     ','lev     ','gw      '        &
            ,'hyam    ','hybm    ','hyai    ','hybi    '        &
            ,'slon    ','slat    ','ilev    ','P0      '        &
            ,('        ',iii=1,88 )/)
! These names should match the grid_names_1d to keep things clear.
! All the possible coordinates (not dimensions) on the caminput.nc file.
type(grid_1d_type), target ::  lon ,lat ,lev ,gw ,hyam ,hybm ,hyai ,hybi, slon ,slat ,ilev, P0

! Other useful 1D grid arrays (for cubed sphere)
real(r8), allocatable :: lon_rad(:), lat_rad(:)   ! longitude and latitude in radians, used by bearings()

! grid_2d_type ?
! integer :: grid_num_2d = 0              ! # of 2d grid fields to read from file
! ? should phis be in grid_names_2d?
! character (len=8),dimension(100) :: grid_names_2d = (/('        ',iii=1,100)/)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!
!----------------------------------------------------------------------
! Namelist variables with default values follow

! output_state_vector = .true.     results in a "state-vector" netCDF file
! output_state_vector = .false.    results in a "prognostic-var" netCDF file
logical :: output_state_vector = .false.

! Files where basic info about model configuration can be found
character(len = 128) :: model_config_file = 'caminput.nc',             & ! An example cam initial file.
                        cam_phis          = 'cam_phis.nc',             & ! Separate source of PHIS/topography.
                        homme_map_file    = 'HommeMapping.nc',         & ! Corners of each cubed sphere cell.
                        cs_grid_file      = 'HommeMapping_cs_grid.nc', & ! Relationships among corners/nodes.
                        model_version     = '4.0'


! Define location restrictions on which observations are assimilated
! (values are calculated anyway, but istatus is set to 2)
character(len = 8) :: vert_coord = 'pressure'            ! or 'log_invP'
real(r8) :: max_obs_lat_degree        = 90.0_r8
real(r8) :: highest_obs_pressure_Pa   = 15000.0_r8
real(r8) :: highest_state_pressure_Pa = 15000.0_r8
! These are not namelist variables, but are related, and calculated from highest_obs_pressure_Pa
real(r8) :: highest_obs_level         = MISSING_R8
real(r8) :: highest_obs_height_m      = MISSING_R8

! Namelist variables and default values for defining state vector.
! Large, arbitrary dimension could be avoided by
! reading in sizes from a first namelist, allocating, setting default values,
! then get values from second namelist.
! Or, allocate with defaults values, read in namelist, deallocate and reallocate.

integer :: state_num_0d = 0              ! # of scalars fields to read from file
integer :: state_num_1d = 0              ! # of 1d fields to read from file
integer :: state_num_2d = 1              ! # of 2d fields to read from file
integer :: state_num_3d = 4              ! # of 3d fields to read from file

! These can't be allocatable since they are namelist items.  They have to
! have a fixed size at compile time.
integer, parameter :: MAX_STATE_NAMES = 100
character(len=8),dimension(MAX_STATE_NAMES) :: state_names_0d = (/('        ',iii=1,MAX_STATE_NAMES)/)
character(len=8),dimension(MAX_STATE_NAMES) :: state_names_1d = (/('        ',iii=1,MAX_STATE_NAMES)/)
character(len=8),dimension(MAX_STATE_NAMES) :: state_names_2d = (/'PS      ',('        ',iii=1,MAX_STATE_NAMES-1)/)

character(len=8),dimension(MAX_STATE_NAMES) :: state_names_3d =  &
          (/'T       ','U       ','V       ', &
            'Q       ', ('        ',iii=1,MAX_STATE_NAMES-4)/)

! NOVERT  need (a) namelist parameter(s) to define which_vert for each 2D (xD?) field
!         There's a danger of having a mismatch with the state_names_Xd; should this
!         definition be part of state_names_Xd, which is parsed into a name and which_vert
!         after being read?  Not for now.

integer , dimension(MAX_STATE_NAMES) :: which_vert_1d = (/(-2,iii=1,MAX_STATE_NAMES)/)
integer , dimension(MAX_STATE_NAMES) :: which_vert_2d = (/(-1,iii=1,MAX_STATE_NAMES)/)
integer , dimension(MAX_STATE_NAMES) :: which_vert_3d = (/( 1,iii=1,MAX_STATE_NAMES)/)


! Is there a way to exclude state_nums from namelist and have those filled in
! the  subroutine which sorts state_names?
! Yes, use two namelists model_nml_1 and model_nml_2 at future date

! list of fields which this code needs to perturb because they're
! constant valued model parameters and show no spread when start_from_restart = .true.
character (len=8),dimension(MAX_STATE_NAMES) :: pert_names     = (/('        ',iii=1,MAX_STATE_NAMES)/)
real(r8)         ,dimension(MAX_STATE_NAMES) :: pert_sd        = (/(-888888.0d0,iii=1,MAX_STATE_NAMES)/)
real(r8)         ,dimension(MAX_STATE_NAMES) :: pert_base_vals = (/(-888888.0d0,iii=1,MAX_STATE_NAMES)/)

! Special for an experiment.  Specify one string kind e.g QTY_CLOUD_LIQUID and
! observations of that kind will only impact other obs and state vars of that
! same kind.  All other kinds of obs and state vars will not be impacted
! by obs of this kind.  A null string means behave as normal.  Kind strings
! are limited by compilers to be 32 chars, since they are declared as params.
character(len=32) :: impact_only_same_kind = ' '
integer           :: impact_kind_index = -1


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
!----------------------------------------------------------------------
! Derived parameters

! make sure static init code only called once
logical :: module_initialized = .false.

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(time_type) :: Time_step_atmos

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Random sequence and init for pert_model_state
logical                 :: first_pert_call = .true.
type(random_seq_type)   :: random_seq
integer                 :: ens_member = 0
logical                 :: do_out

! common message string used by many subroutines
character(len=129) :: string1, string2, string3

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
! on caminit.nc files
! These are filled in trans_coord
integer              :: coord_order
integer, allocatable :: s_dim_3d(:,:), s_dim_2d(:,:), s_dim_1d(  :),  &
                        f_dim_3d(:,:), f_dim_2d(:,:), f_dim_1d(:,:),  &
                        f_dimid_3d(:,:), f_dimid_2d(:,:), f_dimid_1d(:,:),  &
                        s_dimid_3d(:,:), s_dimid_2d(:,:), s_dimid_1d(  :)
integer, dimension(3,3) :: s_dim_max
integer, dimension(4,3) :: f_dim_max


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Surface pressures, used by vertical interpolation routines.
!
! I assume that staggered grids (US and VS) are staggered only 1 direction (each),
! so that surface pressure interpolations to get staggered ps use only 2 A-grid ps values.
! The interpolations for columns of heights are more general, but will do a 2 point interp
!     if the staggering is only in one direction.
!
! Need more arrays if any future fields are doubly staggered.

! ? merge; should this be a grid_2d_type, specified above?
logical               :: alloc_ps=.true.    ! Flag whether to alloc space for ps[_stagr]
real(r8), allocatable :: ps(:, :)           ! surface pressure used to calc P and height profiles.
real(r8), allocatable :: ps_stagr_lon(:, :) ! ps used to calc P profiles & heights on grid staggered
                                            !    East-West (i.e. for VS) relative to ps
real(r8), allocatable :: ps_stagr_lat(:, :) ! ps used to calc P profiles & heights on grid staggered
                                            !    North-South (i.e. for US) relative to ps

! HK DISTRUBUTED
logical               :: alloc_ps_distrib=.true.    ! Flag whether to alloc space for ps[_stagr]
real(r8), allocatable :: ps_distrib(:, :, :)           ! surface pressure used to calc P and height profiles.
real(r8), allocatable :: ps_stagr_lon_distrib(:, :, :) ! ps used to calc P profiles & heights on grid staggered
                                            !    East-West (i.e. for VS) relative to ps
real(r8), allocatable :: ps_stagr_lat_distrib(:, :, :) ! ps used to calc P profiles & heights on grid staggered
                                            !    North-South (i.e. for US) relative to ps
!> @todo 3D pressure
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
real(r8), allocatable :: p_col(:), model_h(:)
real(r8), allocatable :: p_col_distrib(:, :), model_h_distrib(:,:)

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
integer :: ne, np

! Number of columns/nodes in the cubed-sphere grid
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
! Variables used for debugging/diagnostics.  They could be reduced in rank, and
! defined in the subroutines where used, if no mapping diagnostics are needed.
real(r8), allocatable :: x_planar(:,:,:), y_planar(:,:,:)

! Locations of cubed sphere nodes, in DART's location_type format.
type (location_type), allocatable :: cs_locs(:)

! Location of cubed sphere nodes, in cartesian coordinates
type(xyz_location_type), allocatable :: cs_locs_xyz(:)
type(xyz_get_close_type)             :: cs_gc_xyz

! Structure containing grid point locations, etc.,
! defined in static_init_mod after reading in CS lons, lats, and levels,
! needed in model_interpolate:interp_cubed_sphere.
type (get_close_type) :: cs_gc

! array of KINDs of cubed sphere grid points,
integer, allocatable :: cs_kinds(:)

real(r8) :: half_PI           ! Used in function bearing.

! earth radius; needed to convert lat/lon to x,y,z cartesian coords.
real(r8), parameter :: radius = 6371229.0 ! meters


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! CAM3 array 'cflds' is filled with simple loops over state_names_xxx,
! which requires that the field names be provided in the right order.
! Later I should replace that with code which orders namelist input field names
! into cflds, regardless of their original order, and tallies how many of each.
! Is there a way to exclude state_nums from namelist and have those filled in
! the same subroutine?

character (len=8), allocatable :: cflds(:)

! Attribute values for the fields which comprise the state vector
! These are filled by nc_read_model_atts
character (len=128), allocatable :: state_long_names(:)
character (len=128), allocatable :: state_units(:)
! character (len=128), allocatable :: state_units_long_names(:)

! array for the linking of obs_kinds (QTY_) to model field TYPE_s
! It's filled in map_kinds
! The max size of QTY_ should come from obs_kind_mod
! These should be dimensioned the same size as the total of state_names_Nd.
integer, dimension(300) :: dart_to_cam_kinds = (/(MISSING_I,iii=1,300)/)
integer, dimension(300) :: cam_to_dart_kinds = (/(MISSING_I,iii=1,300)/)
!
!-----------------------------------------------------------------------

!#######################################################################

contains

!#######################################################################

! static_init_model section

   subroutine static_init_model()
!=======================================================================
! subroutine static_init_model()
!
! Initializes class data for CAM model (all the stuff that needs to be done once).
! For now, does this by reading info from a fixed name netcdf file.

integer :: iunit, io, i, ncfileid
integer :: max_levs, ierr

! only execute this code once
if (module_initialized) return


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

! set the printed output logical variable to reduce printed output;
! depends on whether this is being called by dart_to_cam (read ens member # from file 'element' )
! or by filter (multiple processes, printout controlled by do_output())

if (file_exist('element')) then
   iunit = get_unit()
   open(unit = iunit, file='element', form = 'formatted')
   read(iunit,*) ens_member
   close(iunit)
   do_out = .false.
   if (ens_member == 1) do_out = .true.
else
   do_out = do_output()
   !write(*,*) 'do_out = ',do_out
   ! static_init_model is called once for each MPI task.
   ! There may be more or fewer ensemble members than tasks.
end if

! Record the namelist values
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(    *      , nml=model_nml)

! Set the model minimum time step from the namelist seconds and days input
Time_step_atmos = set_time(Time_step_seconds, Time_step_days)
! kdr debug
if (print_details .and. do_out) call print_time(Time_step_atmos)

! Open CAM 'initial' file to read dimensions and coordinates of fields.
call nc_check(nf90_open(path = trim(model_config_file), mode = nf90_nowrite, ncid = ncfileid), &
             'static_init_model', 'opening '//trim(model_config_file))

! Get sizes of dimensions/coordinates from netcdf and put in global storage.
! Also may change l_rectang to .false.
call read_cam_init_size(ncfileid)

! Compute overall model size and put in global storage
! s_dim_#d come from read_cam_init_size/trans_coord, and are in global storage
model_size = state_num_0d
do i=1,state_num_1d
   model_size = model_size + s_dim_1d(i)
end do
do i=1,state_num_2d
   model_size = model_size + s_dim_2d(1,i) * s_dim_2d(2,i)
end do
do i=1,state_num_3d
   model_size = model_size + s_dim_3d(1,i) * s_dim_3d(2,i) * s_dim_3d(3,i)
end do
if (do_out) then
   write(string1, '(A,I9)') 'CAM state vector size: ', model_size
   call error_handler(E_MSG, 'static_init_model', string1)
end if

! Allocate space for longitude and latitude global arrays
! and Allocate space for hybrid vertical coord coef arrays
! height; phis
! allocate(lons(num_lons), lats(num_lats), gw(num_lats), hyai(num_levs+1), &
!          hybi(num_levs+1), hyam(num_levs), hybm(num_levs), &
!          phis(num_lons, num_lats) )


! There's a query of caminput.nc within read_cam_coord for the existence of the field.
! The second argument is a grid_1d_type structure
! CS; ncol is a dimension, but there's no coordinate variable of the same name.
call read_cam_coord (ncfileid, lon,   'lon     ')
call read_cam_coord (ncfileid, lat,   'lat     ')
call read_cam_coord (ncfileid, lev,   'lev     ')
call read_cam_coord (ncfileid, ilev,  'ilev    ')
call read_cam_coord (ncfileid, gw,    'gw      ')
call read_cam_coord (ncfileid, slon,  'slon    ')
call read_cam_coord (ncfileid, slat,  'slat    ')

! read hybrid vert coord coefs
call read_cam_coord (ncfileid, hyai, 'hyai    ')
call read_cam_coord (ncfileid, hybi, 'hybi    ')
call read_cam_coord (ncfileid, hyam, 'hyam    ')
call read_cam_coord (ncfileid, hybm, 'hybm    ')

! It's a scalar, but I can put it into the same coord structure as previous fields.
! It's length will be 1
call read_cam_coord (ncfileid, P0, 'P0      ')    ! thats a p-zero

!------------------------------------------------------------------------
! # fields to read
nflds = state_num_0d + state_num_1d + state_num_2d + state_num_3d
if (print_details .and. do_out) &
   write(*, '(A,I3,A,4I3)') '# of fields in state vector =  ', nflds, &
        ' = sum of ', state_num_0d ,state_num_1d ,state_num_2d ,state_num_3d

! CAM3 subroutine to order the state vector parts into cflds
! cflds is needed for reading attributes
allocate (cflds(nflds))
call order_state_fields (cflds, nflds)

!------------------------------------------------------------------------
! CAM3 get field attributes needed by nc_write_model_atts from caminput.nc
allocate (state_long_names(nflds), state_units(nflds))    ! , state_units_long_names(nflds))
call nc_read_model_atts('long_name', state_long_names, nflds, ncfileid)
call nc_read_model_atts('units', state_units, nflds, ncfileid)
! call nc_read_model_atts('units_long_name', state_units_long_names, nflds)
if (.not. l_rectang) then
   ! Read some attributes from the cubed sphere model_config_file.
   ! ne is the number of elements/cube edge.  Usually 0 for refined grids.
   ! np is the number of nodes/element edge (shared with adjacent element.
   call nc_read_global_att(ncfileid, 'ne', ne)
   call nc_read_global_att(ncfileid, 'np', np)

   ! Calculate the nominal resolution of the (coarse) grid,
   ! for use by model_interpolate's call to get_close_obs.
   if (ne == 0) then
      ! Refined mesh; assume the coarsest grid is the default '1-degree'.
      coarse_grid = 1.001_r8 * DEG2RAD
      l_refined = .true.
   else
      ! Standard cubed sphere; there are 3x num_elements/face_edge x 4 nodes
      ! around the equator.  ne = 30 -> 3x4x30 = 360 nodes -> '1-degree'
      coarse_grid = (30.01_r8/ne) * DEG2RAD
   end if
   if (do_out) PRINT*,'Cubed sphere coarse_grid resolution (rad) used in cs_gc definition = ',&
                      coarse_grid

   half_PI = PI*0.5_r8   ! Used in function bearing.

   ! Fill cs_gc for use by model_mod.  Inputs and outputs are in global storage.
   ! ncol must be defined before it's called.
   ncol = dim_sizes(find_name('ncol    ',dim_names))
   call fill_gc

   ! Fill arrays that are useful for bearings and distances.
   allocate(lon_rad(ncol), lat_rad(ncol))
   do i=1,ncol
      lon_rad(i) = lon%vals(i)*DEG2RAD
      lat_rad(i) = lat%vals(i)*DEG2RAD
   end do


end if

call nc_check(nf90_close(ncfileid), &
              'static_init_model', 'closing '//trim(model_config_file))

!------------------------------------------------------------------------
! height
! Get dimensions and surface geopotential from a new netcdf file and test for consistency.

! Open file and read PHIS from it.
call read_cam_2Dreal (cam_phis, 'PHIS    ')

! call read_cam_2Dreal (ncfileid, phis , topog_lons, topog_lats, 'PHIS    ')

! Allocate global variables which will be used in vertical interpolations
max_levs = lev%length
if (ilev%label /= '        ') then
   max_levs = max(ilev%length, lev%length)
end if
allocate (p_col(max_levs), model_h(max_levs))

!------------------------------------------------------------------------
! CS Cubed sphere grid data.
! Read in or create a file containing the relationships among cubed sphere nodes,
! such as neighbors, centers, and bearings, which will be used to identify the cell
! which contains an observation.
! Fields will be stored in global storage.
if (.not. l_rectang) then
   if (file_exist(cs_grid_file)) then
      ! Read the cubed sphere grid file, whose name comes from the namelist
      ierr = nc_read_cs_grid_file()
   else if (file_exist(homme_map_file)) then
      ! Read the HommeMapping.nc file and create the cubed sphere grid arrays.
      call create_cs_grid_arrays
      ! Write the cubed sphere grid arrays to a new NetCDF file.
      ierr = nc_write_cs_grid_file( cs_grid_file, homme_map_file )
      if (ierr /= 0) then
         write(string1, *)'nc_write_cs_grid_file  bombed with error ', ierr
         call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
      end if
   else
      write(string1, *)'No cs_grid_file "',trim(cs_grid_file), &
                    '" nor homme_map_file "',trim(homme_map_file),'"'
      call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
   end if

end if

!------------------------------------------------------------------------
! Arrays for the linking of obs_kinds (QTY_) to model field TYPE_s;
!    dart_to_cam_kinds and cam_to_dart_kinds
call map_kinds()

! nsc fri, 13mar09
! if restricting impact of a particular kind to only obs and state vars
! of the same kind, look up and set the kind index.
if (len_trim(impact_only_same_kind) > 0) then
   impact_kind_index = get_index_for_quantity(impact_only_same_kind)
end if

! Make sure we only come through here once.
module_initialized = .true.

end subroutine static_init_model


   subroutine read_cam_init_size(ncfileid)
!=======================================================================
! subroutine read_cam_init_size(ncfileid)

!
! Gets the number, names, and sizes of field dimensions from a CAM init netcdf file
! in file_name (regardless of dynamical core).
! Called by static_init_model (only).

integer,  intent(in)  :: ncfileid

integer :: i,j

!------------------------------------
! learn how many dimensions are defined in this file.
call nc_check(nf90_inquire(ncfileid, num_dims), 'read_cam_init_size', 'inquire num_dims')

if (allocated(dim_ids) .or. allocated(dim_names) .or. allocated(dim_sizes)) then
   write(string1, *) 'dim_ids, dim_names, and/or dim_sizes already allocated'
   call error_handler(E_ERR, 'read_cam_init_size', string1, source, revision, revdate)
end if

! where to deallocate?
allocate (dim_ids(num_dims), dim_names(num_dims), dim_sizes(num_dims))

! Cycle through dimids until there aren't any more
! Dimension ids are sequential integers on the NetCDF file.
do i = 1,num_dims
   dim_ids(i) = i
   call nc_check(nf90_inquire_dimension(ncfileid, i, dim_names(i), dim_sizes(i)), &
                 'read_cam_init_size', 'inquire for '//trim(dim_names(i)))
   if (print_details .and. do_out) write(*,*) 'Dims info = ',i, trim(dim_names(i)), dim_sizes(i)
   if (dim_names(i) == 'ncol    ') l_rectang = .false.
   if (trim(dim_names(i)) == 'ncol' .and. l_rectang) PRINT*,'I should use trim(dim_names)'
end do

! Find and store shapes of all the state vector fields.  Grouped by rank of fields into
! separate s_dim_RANKd arrays.
! Fields with the same rank can have different shapes and still be handled; efgworo(lat,lon) and
! frac(lat,lev) will both have their shapes stored in X_dim_2d.
!   CS; remove the following complication?
! Also keep track of whether the init file is old (lon,lev,lat) or new (lon,lat,lev).

call trans_coord(ncfileid)

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
end if

if (state_num_2d > 0) then
   f_dim_max(1:3, 2) = maxval(f_dim_2d, dim=2)   ! gets the max values of f_dim_2d (1:3, :)
   s_dim_max(1:2, 2) = maxval(s_dim_2d, dim=2)   ! gets the max values of s_dim_2d (1:2, :)
else
   f_dim_max(1:3, 2) = 0
   s_dim_max(1:2, 2) = 0
end if

if (state_num_3d > 0) then
   f_dim_max(1:4, 3) = maxval(f_dim_3d, dim=2)   ! gets the max values of f_dim_3d (1:4, :)
   s_dim_max(1:3, 3) = maxval(s_dim_3d, dim=2)   ! gets the max values of s_dim_3d (1:3, :)
else
   f_dim_max(1:4, 3) = 0
   s_dim_max(1:3, 3) = 0
end if

!debug   where will this be written?
if (print_details .and. do_out .and. .false.) then
   if (state_num_1d > 0) then
      write(*,*) 's_dim_1d = ',s_dim_1d
      write(*,*) (s_dim_max(i,1),i=1,3)
   end if

   do i=1,2
      write(*,*) 's_dim_2d = ',(s_dim_2d(i,j),j=1,state_num_2d),'s_dim_max = ',s_dim_max(i,2)
   end do

   do i=1,3
      write(*,'(/A,(10I4))') 's_dim_3d = ',(s_dim_3d(i,j),j=1,state_num_3d)
      write(*,'(A,(10I4))') 's_dim_max = ',s_dim_max(i,3)
      write(*,'(A,(10I4))') 'f_dim_3d = ',(f_dim_3d(i,j),j=1,state_num_3d)
      write(*,'(A,(10I4))') 'f_dim_max = ',f_dim_max(i,3)
   end do
end if
!debug end

end subroutine read_cam_init_size

   subroutine trans_coord(ncfileid)
!=======================================================================
! subroutine trans_coord(ncfileid)
!

! CS; simplify this by not considering old CAM coord orders.
! Figure out which coordinates are lon, lat, lev, based on CAM version
! from the namelist, which has form #.#[.#[.#]]
! merge; add state_num_#d , and

integer,            intent(in) :: ncfileid

! local workspace
character (len=4)              :: form_version = '(I0)'
character (len=4)              :: char_version
integer                        :: part, nchars, tot_chars, i, j, k, varid, next
integer, dimension(4)          :: int_version

int_version = (/(0,i=1,4)/)


! Choose order of coordinates based on CAM version
part = 1
nchars=0
char_version = '    '
tot_chars = len_trim(model_version)
do i=1,tot_chars+1
   if ( i == tot_chars+1 .or.  model_version(i:i) == '.' ) then
      write(form_version(3:3),'(I1)') nchars
      read(char_version,form_version) int_version(part)
      part = part + 1
      nchars = 0
      char_version = '    '
   else
      nchars = nchars + 1
      char_version(nchars:nchars) = model_version(i:i)
   end if
end do
if (do_out) then
   if (print_details) then
      write(*,'(A,A10,4(I3,2X))') 'model_version, version(1:4) = ' &
                                  ,model_version,(int_version(i),i=1,4)
   else
      call error_handler(E_MSG, 'trans_coord', 'CAM model version: '//trim(model_version))
   end if
end if

! assume cam3.0.7 format to start
! test on version cam3.0.3
coord_order = 2
if (int_version(1) < 3) then
   coord_order = 1
else if (int_version(1) == 3 .and. int_version(2) == 0 .and. int_version(3) < 3) then
   coord_order = 1
end if

! Cycle through each field's dimension IDs.
! Pick the dimensions needed out of dim_sizes, using the dimension names in dim_names.
! Fill the state dimids according to the order model_mod wants to see.  (lev, lon, lat).

! 3D is easy; lev, lon, lat are always the coordinates,
! and model_mod always wants them in that order.

if (state_num_3d > 0) then
   allocate(s_dim_3d(3,state_num_3d), s_dimid_3d(3,state_num_3d), &
            f_dim_3d(4,state_num_3d), f_dimid_3d(4,state_num_3d))

   s_dim_3d = 0;    s_dimid_3d = 0;
   f_dim_3d = 0;    f_dimid_3d = 0;
end if

do i = 1,state_num_3d
   ! Get variable id for a  3d field
   call nc_check(nf90_inq_varid(ncfileid, state_names_3d(i), varid), &
                 'trans_coord', 'inq_varid '//trim(state_names_3d(i)))
   ! Get dimension ids for the dimensions of the field
   call nc_check(nf90_inquire_variable(ncfileid, varid, dimids=f_dimid_3d(1:4,i)), &
                 'trans_coord', 'inquire_variable'//trim(state_names_3d(i)))

   Alldim3: do j = 1,4                          ! time and 3 space
      k = f_dimid_3d(j,i)                       ! shorthand; the dimid of this fields current dim
      f_dim_3d(j,i) = dim_sizes(k)
!      do k = 1,num_dims                        ! all dimensions used in the CAM file
!         if (dim_ids(k) == f_dimid_3d(j,i)) then
      ! Put the dimensions we want in the state field positions we want.
      if (dim_names(k) == 'lev' .or. dim_names(k) == 'ilev') then
         s_dim_3d  (1,i) = dim_sizes(k)
         s_dimid_3d(1,i) = k
!         s_dimid_3d(1,i) = f_dimid_3d(j,i)
      else if (dim_names(k) == 'lon' .or. dim_names(k) == 'slon') then
         s_dim_3d  (2,i) = dim_sizes(k)
         s_dimid_3d(2,i) = k
      else if (dim_names(k) == 'lat' .or. dim_names(k) == 'slat') then
         s_dim_3d  (3,i) = dim_sizes(k)
         s_dimid_3d(3,i) = k
      end if
!            cycle Alldim3
!         end if
!      end do
   end do Alldim3
   if (   s_dim_3d(1,i) == 0 .or.  s_dim_3d(2,i) == 0 .or.  s_dim_3d(3,i) == 0 ) then
      call error_handler(E_ERR, 'trans_coord', &
          'num_[lons,lats,levs] was not assigned and = 0' , source, revision, revdate)
   end if
end do

! fill dimids according to the order model_mod wants to see.
! 2d is trickier; 2 of (lev, lon, lat) in that order

if (state_num_2d > 0) then
   allocate(s_dim_2d(2,state_num_2d), s_dimid_2d(2,state_num_2d), &
            f_dim_2d(3,state_num_2d), f_dimid_2d(3,state_num_2d))

   s_dim_2d = 0;  s_dimid_2d  = 0;
   f_dim_2d = 0;  f_dimid_2d  = 0;
end if

do i = 1,state_num_2d
   call nc_check(nf90_inq_varid(ncfileid, state_names_2d(i), varid), &
              'trans_coord', 'inq_varid '//trim(state_names_2d(i)))
   call nc_check(nf90_inquire_variable(ncfileid, varid, dimids=f_dimid_2d(1:3,i)), &
              'trans_coord', 'inquire_variable '//trim(state_names_2d(i)))

   ! extract spatial dimids from the fields dimids
   next = 1
   Alldim2: do j = 1,3      ! time and 2 space
      k = f_dimid_2d(j,i)
      f_dim_2d(j,i) = dim_sizes(k)
!      do k = 1,num_dims
!         if (dim_ids(k) == f_dimid_2d(j,i)) then
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
      else if (dim_names(k) == 'lon' .or. dim_names(k) == 'slon') then
         ! longitude always comes first on CAM initial files.
         ! Otherwise, I'll need a test like for levs, but more complicated.
         s_dim_2d  (next,i) = dim_sizes(k)
         s_dimid_2d(next,i) = k
         next = 2
      else if (dim_names(k) == 'ncol' ) then
         ! CS add ncol
         s_dim_2d  (next,i) = dim_sizes(k)
         s_dimid_2d(next,i) = k
         next = 2
      else if (dim_names(k) == 'lat' .or. dim_names(k) == 'slat' ) then
         s_dim_2d  (next,i) = dim_sizes(k)
         s_dimid_2d(next,i) = k
      end if
!      cycle Alldim2
!         end if
!      end do
   end do Alldim2
   if (   s_dim_2d(1,i) == 0 .or.  s_dim_2d(2,i) == 0 ) then
      ! CS add ncol
      call error_handler(E_ERR, 'trans_coord', &
          'num_[lons,lats,levs,ncol] was not assigned and = 0' , source, revision, revdate)
   end if
end do

! 1d fields are easy.

if (state_num_1d > 0) then
   allocate(s_dim_1d(  state_num_1d), s_dimid_1d(  state_num_1d))
   allocate(f_dim_1d(2,state_num_1d), f_dimid_1d(2,state_num_1d))
   s_dim_1d = 0;   s_dimid_1d = 0;
   f_dim_1d = 0;   f_dimid_1d = 0;
end if

do i = 1,state_num_1d
   call nc_check(nf90_inq_varid       (ncfileid, state_names_1d(i), varid), &
              'trans_coord', 'inq_varid '//trim(state_names_1d(i)))
   call nc_check(nf90_inquire_variable(ncfileid, varid, dimids=f_dimid_1d(1:2,i)), &
              'trans_coord', 'inq_varid '//trim(state_names_1d(i)))

   Alldim1: do j = 1,2       ! time and 1 space
      k = f_dimid_1d(j,i)
      f_dim_1d(j,i) = dim_sizes(k)
!      do k = 1,num_dims
!         if (dim_ids(k) == f_dimid_1d(j,i)) then
! CS add ncol
      if (dim_names(k) == 'lon' .or. dim_names(k) == 'slon' .or. &
          dim_names(k) == 'ncol'.or. &
          dim_names(k) == 'lat' .or. dim_names(k) == 'slat' .or. &
          dim_names(k) == 'lev' .or. dim_names(k) == 'ilev' ) then
         s_dim_1d(i) = dim_sizes(k)
         s_dimid_1d(i) = k
!         s_dimid_1d(i) = f_dimid_1d(j,i)
      end if
!            cycle Alldim1
!         end if
!      end do
   end do Alldim1

   if ( s_dim_1d(i) == 0 ) then
      write(string1, '(A,I3,A)') ' state 1d dimension(',i,') was not assigned and = 0'
      call error_handler(E_ERR, 'trans_coord',trim(string1), source, revision, revdate)
   end if
end do

end subroutine trans_coord


   subroutine read_cam_2Dreal(file_name, cfield)
!======================================================
!  subroutine read_cam_2Dreal(file_name, cfield)

! Subroutine to read in 2D/horizontal CAM fields, such as PHIS.
! Handles both logicall rectangular arrays (FV/Eul) and irregular (SE-CAM/cubed-sphere).
! Can NOT be called with cfield = a record variable  (time,[dim2,]dim1).

!------------------------------------------------------
character (len=*), intent(in)  :: file_name
character (len=8), intent(in)  :: cfield

!------------------------------------------------------
integer :: ncfileid, ncfldid      ! NetCDF variables
integer :: field_dimids(3) = MISSING_I    ! Array of dimension IDs for cfield
                                          ! (2 space (FV) and time dimension (CAM .h0. files).
integer :: i_dim1, i_dim2         ! Variables to reference the dimension(s) of cfield
integer :: num_dim1, num_dim2     ! NetCDF file variable dimension sizes, for comparison to file_name's
integer :: slon_index, slat_index, lon_index, lat_index !indices of [s]lon and [s]lat
                                                        ! within the list of dimensions
integer :: n,m
character (len=NF90_MAX_NAME) :: name_dim1,name_dim2    ! Names of dimensions of cfield
real(r8), allocatable         :: var(:,:)               ! Temp array used by nc_get_var

! If cam_phis doesn't exist, try to find a CAM history (h0) file that
! may have the cfield on it.
! ERROR: "This name (system) does not have a type, and must have an explicit type."
if (file_name == cam_phis .and. .not.file_exist(trim(file_name))) then
!   ! CS; add code to scripts and namelists to locate a h0 file and link to it.
!   !
!   m = system('set h0 = `ls *.h0.*`; if ($status == 0) ln -s $h0[1] '//cam_phis)
   write(string1,'(2A)') trim(file_name),  &
        ' is missing; trying to find a CAM history file (h0) to provide '//cfield
   call error_handler(E_WARN, 'read_cam_2Dreal', trim(string1), source, revision, revdate)
end if

! Open the file and get dimension information.
if (file_exist(trim(file_name))) then
   call nc_check(nf90_open(path = trim(file_name), mode = nf90_nowrite, ncid = ncfileid), &
              'static_init_model:read_cam_2Dreal', 'opening '//trim(file_name))
   if (print_details .and. do_out) write(*, *) 'file_name for ',cfield,' is ', trim(file_name)

   ! get field id
   call nc_check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid), &
              'read_cam_2Dreal', 'inq_varid: '//cfield)

   ! get dimension 'id's
   call nc_check(nf90_inquire_variable(ncfileid, ncfldid, dimids = field_dimids), &
              'read_cam_2Dreal', 'inquire_variable: '//cfield)

   ! get dimension sizes
   ! The first spatial dimension is always present.
   call nc_check(nf90_inquire_dimension(ncfileid, field_dimids(1), name_dim1, num_dim1 ), &
                 'read_cam_2Dreal', 'inquire_dimension: '//name_dim1)
   if (field_dimids(2) == MISSING_I)  then
      num_dim2 = 1
      name_dim2 = 'no2ndDim'
   else
      call nc_check(nf90_inquire_dimension(ncfileid, field_dimids(2), name_dim2, num_dim2 ), &
                    'read_cam_2Dreal', 'inquire_dimension: '//name_dim2)
   end if

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
      end if

      if (field_dimids(2) /= MISSING_I) then
         i_dim2 = dim_sizes(find_name(name_dim2,dim_names))
         if ( num_dim2 /= i_dim2 ) then
            write(string1,'(A,2I8,A)') 'i_dim2, num_dim2, name_dim2 =', &
                                          i_dim2, num_dim2, trim(name_dim2)
            call error_handler(E_MSG, 'read_cam_2Dreal', trim(string1), source, revision, revdate)
            write(string1,'(A,4I12)') 'horizontal dimensions mismatch of initial files and topog ', &
                  i_dim2, num_dim2
            call error_handler(E_ERR, 'read_cam_2Dreal', trim(string1), source, revision, revdate)
         end if
      end if
   end if
else
   write(string1,'(2A)') trim(file_name),  &
        ' is missing; I do not know how to find it.'
   call error_handler(E_ERR, 'read_cam_2Dreal', trim(string1), source, revision, revdate)
end if

! Allocate local arrays, based on size of this variable on the file.
allocate (var(num_dim1, num_dim2))

! Read surface geopotential from cam_phis for use in vertical interpolation in height.
! Coordinate order not affected by CAM model version.
   call nc_check(nf90_get_var(ncfileid, ncfldid, var, start=(/1,1/), &
                 count=(/num_dim1, num_dim2/)), 'read_cam_2Dreal', trim(cfield))

! assign values to phis grids for use by the rest of the module.
if (cfield == 'PHIS    ') then

   if (alloc_phis) allocate (phis(num_dim1, num_dim2))
   ! Don't want to set alloc_phis = false yet; there may be staggered phis to set.
   phis(1:num_dim1,1:num_dim2) = var

   ! If needed, generate phis on the staggered grids.
   slon_index = find_name('slon    ',dim_names)
   slat_index = find_name('slat    ',dim_names)
   lat_index  = find_name('lat     ',dim_names)
   lon_index  = find_name('lon     ',dim_names)

   ! CS these sections will be skipped for grids with no staggered grids, e.g. cubed sphere
   if (slon_index /= 0) then
      if (alloc_phis) allocate (phis_stagr_lon (dim_sizes(slon_index), dim_sizes( lat_index)))
      do n=1,dim_sizes( lat_index)
         ! CS Would a better interpolation (e.g. splines) help anything in a meaningful way?
         ! ? Better interpolation using quad_interp?  But it's 'bilinear', so no differenc?
         phis_stagr_lon(1,n) = .5 * (phis(1,n) + phis(dim_sizes(lon_index),n))
         do m=2,dim_sizes(slon_index)
            phis_stagr_lon(m,n) = .5 * (phis(m-1,n) + phis(m,n))
         end do
      end do
   end if

   if (slat_index /= 0) then
      if (alloc_phis) allocate (phis_stagr_lat (dim_sizes( lon_index), dim_sizes(slat_index)))
      do n=1,dim_sizes(slat_index)
         do m=1,dim_sizes( lon_index)
            phis_stagr_lat(m,n) = .5 * (phis(m,n) + phis(m,n+1))
         end do
      end do
   end if
   alloc_phis = .false.

end if

call nc_check(nf90_close(ncfileid), 'read_cam_2Dreal', 'closing '//trim(file_name))

deallocate (var)

end subroutine read_cam_2Dreal

   subroutine read_cam_2Dint(file_name, cfield, field, num_dim1, num_dim2)
!======================================================
!  subroutine read_cam_2Dint(file_name, cfield, field, num_dim1, num_dim2)
!
! Read 2d integer field from, e.g., HommeMapping.nc
! Called by create_cs_grid_arrays (from static_init_model).

!------------------------------------------------------
character (len=*),    intent(in)  :: file_name
character (len=*),    intent(in)  :: cfield
integer, allocatable, intent(out) :: field(:,:)

!------------------------------------------------------
integer :: ncfileid, ncfldid                            !NetCDF variables
integer :: field_dimids(2) = MISSING_I                  !Array of dimension IDs for cfield
integer                       :: num_dim1, num_dim2     !The dimension(s) of cfield
character (len=NF90_MAX_NAME) :: name_dim1,name_dim2    !Names of dimensions of cfield

if (file_exist(file_name)) then
   call nc_check(nf90_open(path = trim(file_name), mode = nf90_nowrite, ncid = ncfileid), &
              'read_cam_2Dint', 'opening '//trim(file_name))
   if (print_details .and. do_out) write(*, *) 'file_name for ',cfield,' is ', trim(file_name)

   ! get field id
   call nc_check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid), &
              'read_cam_2Dint', 'inq_varid: '//cfield)

   ! get dimension 'id's
   call nc_check(nf90_inquire_variable(ncfileid, ncfldid, dimids = field_dimids), &
              'read_cam_2Dint', 'inquire_variable: '//cfield)

   ! get dimension sizes
   ! The first spatial dimension is always present.
   call nc_check(nf90_inquire_dimension(ncfileid, field_dimids(1), name_dim1, num_dim1 ), &
                 'read_cam_2Dint', 'inquire_dimension: '//name_dim1)
   if (field_dimids(2) /= MISSING_I)  then
      call nc_check(nf90_inquire_dimension(ncfileid, field_dimids(2), name_dim2, num_dim2 ), &
                    'read_cam_2Dint', 'inquire_dimension: '//name_dim2)
   else
      num_dim2 = 1
      name_dim2 = 'no2ndDim'
   end if

   if (print_details .and. do_out) &
      write(*,*) cfield,' dimensions num_dim1, num_dim2 = ',num_dim1, num_dim2
else
   write(string1,'(2A)') trim(file_name),  &
        ' is missing; I do not know how to find it.'
   call error_handler(E_ERR, 'read_cam_2Dint', trim(string1), source, revision, revdate)
end if

! Allocate array, based on size of this variable on the file.
allocate (field(num_dim1,num_dim2))

if (field_dimids(2) /= MISSING_I)  then
   call nc_check(nf90_get_var(ncfileid, ncfldid, field, start=(/1,1/), &
                 count=(/num_dim1, num_dim2/)), 'read_cam_2Dint', trim(cfield))
else
   call nc_check(nf90_get_var(ncfileid, ncfldid, field),  &
                  'read_cam_2Dint', trim(cfield))
   !call nc_check(nf90_get_var(ncfileid, ncfldid, field, start=(/1/), &
   !              count=(/num_dim1/)), 'read_cam_2Dint', trim(cfield))
end if

call nc_check(nf90_close(ncfileid), 'read_cam_2Dint', 'closing '//trim(file_name))

end subroutine read_cam_2Dint


   subroutine create_cs_grid_arrays
!=======================================================================
! subroutine create_cs_grid_arrays
!
! Subroutine to create arrays of relationships between cubed sphere nodes (corners)
! and cell centers, including bearings between nodes.
! These will be used to identify the cell containing an observation.
! The relationships read from HommeMapping.nc will be augmented.
! All will be stored in global storage, and written to a new file for
! subsequent use.

! Local variables
!  ncenters, ncorners
integer, dimension(4) :: sh_corn, n
integer :: col, nbr, c, cent
! integer :: itemp, mask(max_neighbors), num_n, min_ind(1)
! real(r8) :: temp,lon1, lat1, lon2, lat2 , dist, angle
integer :: num_n, min_ind(1)
real(r8) :: dist, angle
!real(r8), dimension(3) :: bearings, x_planar, y_planar
real(r8), dimension(3) :: bearings

! Get array of corner nodes/columns which define the cells (identified by 'center').
if (file_exist(homme_map_file)) then
   call read_cam_2Dint(homme_map_file, 'element_corners', corners,ncenters,ncorners)
   ! done in static_init_mod:  ncol = dim_sizes(find_name('ncol    ',dim_names))
   if (ncenters /= (ncol -2) ) then
      write(string1, *) trim(homme_map_file),' ncenters inconsistent with ncol-2 ', ncenters, ncol
      call error_handler(E_ERR,'create_cs_grid_arrays',string1,source,revision,revdate)
   else
      allocate ( num_nghbrs           (ncol), &
                 centers(max_neighbors,ncol), &
                 a          (3,ncorners,ncenters),   &
                 b          (2,ncorners,ncenters),   &
                 x_ax_bearings(ncorners,ncenters),   &
                 x_planar   (3,ncorners,ncenters),   &
                 y_planar   (3,ncorners,ncenters)    )
      ! Allocated by nc_read_cs_grid           corners  (ncorners,ncenters),   &
                 !diagonals(max_neighbors,ncol),  &
                 !neighbors(max_neighbors,ncol),  &
      ! Initialize the grid variables
      num_nghbrs = 0;
      centers = MISSING_I
      a        = MISSING_R8;        b = MISSING_R8
      x_planar = MISSING_R8; y_planar = MISSING_R8
      x_ax_bearings = MISSING_R8
   end if
else
   write(string1, *) trim(homme_map_file),' can not be found '
   call error_handler(E_ERR,'create_cs_grid_arrays',string1,source,revision,revdate)
end if

! Invert the element_corners array to compile all of the neighbors of each node (corner).
! Loop over HommeMapping cell centers.
Quads: do cent = 1,ncenters
   Corns: do c = 1,4
      ! Get the node numbers that define this cell
      ! and shift them to create a separate mapping for each corner/node.
      ! Shift the section of corners 1 place to the 'left'/lower for the first corner,
      ! 2 for the 2nd, etc.  This will put the node closest to the ob in position 4
      ! (of the shifted corners). Then the (x,y) origin will the the closest node,
      ! and the indexing of the a,b,x_ax_bearing arrays will be easy.
      ! Shifting preserves the order of the corners as we go around the cell (clockwise).
      sh_corn = cshift(corners(cent,:), c)
      if (sh_corn(4) < 10) write(*,'(A,5I7,/31X,4I7)') 'c, 4 corners, shifted = ', &
         c,(corners(cent,nbr),nbr=1,4),(sh_corn(nbr),nbr=1,4)
      ! Increment the number of neighbors (and centers ) this corner(node) has.
      ! sh_corn(4) is used for all cases because the corner we're working on always
      ! ends up in that position, when c is incremented, then the corners are shifted.
      n(c) = num_nghbrs(sh_corn(4)) + 1
      ! Update the number of neighbors of each corner of the cell,
      num_nghbrs(    sh_corn(4)) = n(c)
      ! Store the info that this center is associated with 4 node->neighbor pairs.
      centers(n(c),sh_corn(4)) = cent

      ! Define the planar coordinates for this center/cell and this corner/node.
      ! The 4th corner is the origin, and the cell side from the 4th to the 3rd is
      ! the x-axis of this cell's coordinate system for this corner.
      ! This is established in the definition of bearings().
      ! This choice makes mapping coefficients a(0) and b(0) = 0.
      ! It also helps make the indexing of bearings easy to use and store.

      ! Descend through neighbors so that bearings(3) is defined when needed at loop end.
      if (sh_corn(4) < 10) then
         write(*,'(A,3F10.6)') 'create_cs_grid:    lon1, lat1 = ', &
              lon_rad(sh_corn(4)), lat_rad(sh_corn(4))
      end if
      do nbr = 3,1,-1
         ! Bearings from the current origin node of the cell to the other 3.
         bearings(nbr) = bearing(lon_rad(sh_corn(4)),   lat_rad(sh_corn(4)),  &
                                 lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)) )

         dist = get_dist(cs_locs(sh_corn(4)), cs_locs(sh_corn(nbr)), 0, 0, .true.)

         if (sh_corn(4) < 10) then
            write(*,'(A,3F10.6)') 'create_cs_grid:    lon2, lat2, bearing = ', &
                 lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)), bearings(nbr)
         end if

         ! This order looks wrong, but we need to change the sign of angles from the
         ! clockwise direction used by bearings to the counterclockwise direction used by
         ! trig functions.
         angle = bearings(3) - bearings(nbr)
         ! Normalize to -PI < angle <= PI.
         angle = mod(angle,PI) - PI*int(angle/PI)
         ! Set the planar location of this corner/node.
         x_planar(nbr,c,cent) = dist * cos(angle)
         y_planar(nbr,c,cent) = dist * sin(angle)
      end do

      ! Store the baseline for use when interpolating to an ob location.
      x_ax_bearings(c,cent) = bearings(3)
      ! Define another bearings array to allow coord_ind_cs to find the right
      ! cell around the closest node by using a search through bearings,
      ! rather than a call to unit_square_location.
      !   Propagate sort_bearings to writing and reading of HommeMapping_cs_grid.nc file.
      !   real(r8), allocatable, :: sort_bearings(max_neighbors,ncol)
      !   sort_bearings(n(c),sh_corn(4)) = bearings(3)

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
      a(3,c,cent) = x_planar(3,c,cent)
      a(2,c,cent) = x_planar(1,c,cent)
      a(1,c,cent) = x_planar(2,c,cent) - x_planar(1,c,cent) - x_planar(3,c,cent)
      b(2,c,cent) = y_planar(1,c,cent)
      b(1,c,cent) = y_planar(2,c,cent) - y_planar(1,c,cent)

      if (cent < 10) then
         write(*,'(A,1p4E12.4)') 'create_cs_grid_arrays: a = ',(a(nbr,c,cent),nbr=1,3)
         write(*,'(A,1p4E12.4)') 'create_cs_grid_arrays: b = ',(b(nbr,c,cent),nbr=1,2)
      end if

      if (a(3,c,cent)* a(2,c,cent) *a(1,c,cent) == 0._r8) then
         write(*,'(A,2I8,A,1p3E12.4)') 'a(:,',c,cent,') = ',(a(nbr,c,cent),nbr=1,3)
         write(*,'(A,(6X,2F10.6))') 'create_cs_grid:    lon, lat for nghr=1-3  = ', &
              (lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)), nbr=1,3)
         write(*,'(A,5I7,/31X,4I7)') 'c, 4 corners, shifted = ', &
               c,(corners(cent,nbr),nbr=1,4),(sh_corn(nbr),nbr=1,4)
      end if
      if (b(2,c,cent) *b(1,c,cent) == 0._r8) then
         write(*,'(A,2I8,A,1p3E12.4)') 'b(:,',c,cent,') = ',(b(nbr,c,cent),nbr=1,2)
         write(*,'(A,(6X,2F10.6))') 'create_cs_grid:    lon, lat for nghr=1-3  = ', &
              (lon_rad(sh_corn(nbr)), lat_rad(sh_corn(nbr)), nbr=1,3)
         write(*,'(A,5I7,/31X,4I7)') 'c, 4 corners, shifted = ', &
               c,(corners(cent,nbr),nbr=1,4),(sh_corn(nbr),nbr=1,4)
      end if

      ! Don't bother for now:
      ! I could save parts of the coefficients of the equation for m
      ! which maps the planar coordinates to the unit square.
      ! The other parts of the coefficients depend on the (x,y) of the ob,
      ! and aren't available at this point (see unit_square_location()).  That is:
      ! m^2 *axb(1)  + m*[axb(2) +f(x,y,a,b)]  + [aXb(3) +g(x,y,a,b)] = 0

      !   aXb(1,quad) = a(1)*b(2) - a(2)*b(1)
      !   aXb(2,quad) = a(3)*b(2)
      !   !aXb(3,quad) = 0.  ! Terms eliminated by choice of origin and x-axis.
      !   if (quad < 10) then
      !      write(*,'(A,I8,A,1p4E12.4)') 'quad_interp: aXb(1:3,',quad,') = ',aXb(1:3,quad)
      !   end if

   end do Corns
end do Quads


! Check that all nodes have at least 3 neighbors and no more than 6.
do col = 1,ncol
   if (num_nghbrs(col) < 3 .or. num_nghbrs(col) > max_neighbors) then
      write(*,'(A,I6,A,6I8)') 'ERROR: num_nghbrs(',col,') different than expected = ', &
           num_nghbrs(col)
   end if
end do

! Reorder the neighbors so that they are sequential around each node.
! This will be clockwise, from the neighbor with the least bearing, to the greatest.
! -PI < bearing <= PI.

!do col = 1,ncol
!   num_n = num_nghbrs(col)
!   bear(1:num_n) = sort_bearings(1:num_n,col)
!   if (col < 10) then
!      write(*,'(A,I6,A,6F10.6)') 'bearings(:,',col,') = ',(bearings(nbr,col),nbr=1,num_n)
!   end if
!   ! Create a mask of indices to use for reordering all of the neighbor arrays.
!   mask(1:num_n) = (/ (nbr,nbr=1,num_n) /)
!   ! Successively move the remaining neighbor with smallest bearing into the next spot.
!   ! The num_n-1 iteration puts the last 2 elements in the right positions.
!   do nbr = 1,num_n-1
!      ! Find index of the minimum bearing among the ones that haven't been sorted.
!      ! minloc returns the index within the subarray that's passed to it,
!      ! so the bear indices that have been excluded must be added back on.
!      min_ind = minloc(bear(nbr:num_n)) + nbr - 1
!      ! Move the current first index and value out of the way.
!      itemp = mask(nbr)
!      temp  = bear(nbr)
!      ! Put the smallest remaining value and its index in the current first position.
!      mask(nbr) = mask(min_ind(1))
!      bear(nbr) = bear(min_ind(1))
!      ! Put the out of the way value and its index in the 'vacated position' for further sorting.
!      mask(min_ind(1)) = itemp
!      bear(min_ind(1)) = temp
!      if (col < 10) then
!         write(*,'(A,I6,6F10.6)') '   nbr, bear(:) = ',nbr,bear(:)
!         write(*,'(A,6I4)')       '        mask(:) = ',mask(:)
!      end if
!   end do
!
!   ! Use the mask to re-arrange the neighbors in all the arrays.
!   if (col < 10) then
!      write(*,'(A,6I8)') 'neighbors before mask = ',neighbors(1:num_n,col)
!   end if
!   sort_bearings (1:num_n,col) = bear(1:num_n)
!   centers       (1:num_n,col) = centers  ( (/ (mask(nbr),nbr=1,num_n) /), col)
!   if (col < 10) then
!      write(*,'(A,6I8)') 'neighbors after mask = ',neighbors(1:num_n,col)
!      write(*,*) '----------------- '
!   end if
!end do

return

end subroutine create_cs_grid_arrays


   function nc_write_cs_grid_file( cs_grid_file, homme_map_file ) result (ierr)
!=======================================================================
!  function nc_write_cs_grid_file( cs_grid_file, homme_map_file ) result (ierr)

! Write out number of neighbors, neighbors, corners, centers, and bearings
! to a netCDF file once for this grid at the beginning of the assimilation.
! Called by static_init_model.

character(len=*), intent(in) :: cs_grid_file, homme_map_file
integer :: ierr, ncFileID

integer ::                               &
        ncentersID,       centersVarID,  &
        ncornersID,       cornersVarID,  &
               aID,             aVarID,  &
               bID,             bVarID,  &
   x_ax_bearingsID, x_ax_bearingsVarID,  &
            ncolID,                      &
   max_neighborsID,    num_nghbrsVarID

ierr = 0

! Create the file
call nc_check(nf90_create(path = trim(cs_grid_file), cmode = NF90_SHARE, ncid = ncFileID), &
              'nc_write_cs_grid_file', 'create '//trim(cs_grid_file))

write(string1,*) trim(cs_grid_file),' is ncFileID ',ncFileID
call error_handler(E_MSG,'nc_write_cs_grid_file',string1,source,revision,revdate)

! Define the dimensions
call nc_check(nf90_def_dim(ncid=ncFileID,                                          &
              name="ncenters",      len = ncenters,      dimid = ncentersID), &
              'nc_write_cs_grid_file', 'def_dim ncenters '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=ncFileID,                                          &
              name="ncorners",      len = ncorners,      dimid = ncornersID), &
              'nc_write_cs_grid_file', 'def_dim ncorners '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=ncFileID,                                  &
              name="max_neighbors", len = max_neighbors, dimid = max_neighborsID), &
              'nc_write_cs_grid_file', 'def_dim max_neighbors'//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=ncFileID,                                  &
              name="ncol",          len = ncol,          dimid = ncolID), &
              'nc_write_cs_grid_file', 'def_dim ncol '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=ncFileID,                                  &
              name="ncoef_a",          len = 3,          dimid = aID), &
              'nc_write_cs_grid_file', 'def_dim a '//trim(cs_grid_file))
call nc_check(nf90_def_dim(ncid=ncFileID,                                  &
              name="ncoef_b",          len = 2,          dimid = bID), &
              'nc_write_cs_grid_file', 'def_dim b '//trim(cs_grid_file))

! Write Global Attributes
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "title", trim(cs_grid_file)), &
              'nc_write_cs_grid_file',   'put_att title '//trim(cs_grid_file))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_mod_source", source ), &
              'nc_write_cs_grid_file',   'put_att model_mod_source '//trim(cs_grid_file))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_mod_revision", revision ), &
              'nc_write_cs_grid_file',   'put_att model_mod_revision '//trim(cs_grid_file))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_mod_revdate", revdate ), &
              'nc_write_cs_grid_file',   'put_att model_mod_revdate '//trim(cs_grid_file))

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "elements_per_cube_edge", ne ), &
              'nc_write_cs_grid_file',   'put_att elements_per_cube_edge '//trim(cs_grid_file))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "nodes_per_element_edge", np ), &
              'nc_write_cs_grid_file',   'put_att nodes_per_elements_edge '//trim(cs_grid_file))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "HommeMapping_file", homme_map_file ), &
              'nc_write_cs_grid_file',   'put_att HommeMapping_file '//trim(cs_grid_file))

! Create variables and attributes.
call nc_check(nf90_def_var(ncid=ncFileID, name="num_nghbrs", xtype=nf90_int, &
              dimids=(/ ncolID /), varid=num_nghbrsVarID),  &
              'nc_write_cs_grid_file', 'def_var num_nghbrs')
call nc_check(nf90_put_att(ncFileID, num_nghbrsVarID, "long_name", &
              "number of neighbors of each node/column"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(ncFileID, num_nghbrsVarID, "units",     "none"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(ncFileID, num_nghbrsVarID, "valid_range", &
              (/ 1,max_neighbors /)), 'nc_write_cs_grid_file', 'put_att valid_range')

call nc_check(nf90_def_var(ncid=ncFileID, name="centers", xtype=nf90_int, &
              dimids=(/ max_neighborsID,ncolID /), varid=centersVarID),  &
              'nc_write_cs_grid_file', 'def_var centers')
call nc_check(nf90_put_att(ncFileID, centersVarID, "long_name", &
              "cells which use node/column as a corner"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(ncFileID, centersVarID, "units",     "none"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(ncFileID, centersVarID, "valid_range", &
              (/ 1, ncenters /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(ncFileID, centersVarID, "missing_value", &
              (/ MISSING_I /)), 'nc_write_cs_grid_file', 'put_att missing_value')

call nc_check(nf90_def_var(ncid=ncFileID, name="corners", xtype=nf90_int, &
              dimids=(/ ncentersID,ncornersID /), varid=cornersVarID),  &
              'nc_write_cs_grid_file', 'def_var corners')
call nc_check(nf90_put_att(ncFileID, cornersVarID, "long_name", &
              "corners/nodes of each cell "), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(ncFileID, cornersVarID, "units",     "none"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(ncFileID, cornersVarID, "valid_range", &
              (/ 1, ncol /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(ncFileID, cornersVarID, "missing_value", &
              (/ MISSING_I /)), 'nc_write_cs_grid_file', 'put_att missing_value')

call nc_check(nf90_def_var(ncid=ncFileID, name="a", xtype=nf90_double, &
              dimids=(/ aID,ncornersID,ncentersID /), varid=aVarID),  &
              'nc_write_cs_grid_file', 'def_var a')
call nc_check(nf90_put_att(ncFileID, aVarID, "long_name",  &
              "Coefficients of mapping from planar x coord to unit square"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(ncFileID, aVarID, "units",     "none"), &
              'nc_write_cs_grid_file', 'units')
!call nc_check(nf90_put_att(ncFileID, aVarID, "valid_range", &
!              (/ 1, ncol /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(ncFileID, aVarID, "missing_value", &
              (/ MISSING_R8 /)), 'nc_write_cs_grid_file', 'put_att missing_value')


call nc_check(nf90_def_var(ncid=ncFileID, name="b", xtype=nf90_double, &
              dimids=(/ bID,ncornersID,ncentersID /), varid=bVarID),  &
              'nc_write_cs_grid_file', 'def_var b')
call nc_check(nf90_put_att(ncFileID, bVarID, "long_name", &
              "Coefficients of mapping from planar y coord to unit square"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(ncFileID, bVarID, "units",     "none"), &
              'nc_write_cs_grid_file', 'units')
!call nc_check(nf90_put_att(ncFileID, bVarID, "valid_range", &
!              (/ 1, ncol /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(ncFileID, bVarID, "missing_value", &
              (/ MISSING_R8 /)), 'nc_write_cs_grid_file', 'put_att missing_value')

call nc_check(nf90_def_var(ncid=ncFileID, name="x_ax_bearings", xtype=nf90_double, &
              dimids=(/ ncornersID,ncentersID /), varid=x_ax_bearingsVarID),  &
              'nc_write_cs_grid_file', 'def_var x_ax_bearings')
call nc_check(nf90_put_att(ncFileID, x_ax_bearingsVarID, "long_name", &
              "bearing (clockwise from North) from origin node(corner 4) of each mapping to corner 3"), &
              'nc_write_cs_grid_file', 'long_name')
call nc_check(nf90_put_att(ncFileID, x_ax_bearingsVarID, "units",     "radians"), &
              'nc_write_cs_grid_file', 'units')
call nc_check(nf90_put_att(ncFileID, x_ax_bearingsVarID, "valid_range", &
              (/ -PI, PI /)), 'nc_write_cs_grid_file', 'put_att valid_range')
call nc_check(nf90_put_att(ncFileID, x_ax_bearingsVarID, "missing_value", &
              (/ MISSING_R8 /)), 'nc_write_cs_grid_file', 'put_att missing_value')

! Leave define mode so we can fill
call nc_check(nf90_enddef(ncFileID), 'nc_write_cs_grid_file', 'enddef '//trim(cs_grid_file))

! sync to disk, but leave open
call nc_check(nf90_sync(ncFileID), 'nc_write_cs_grid_file', 'sync '//trim(cs_grid_file))

! Fill the variables
call nc_check(nf90_put_var(ncFileID, num_nghbrsVarID, num_nghbrs),  &
              'nc_write_model_vars ','put_var num_nghbrs ')
call nc_check(nf90_put_var(ncFileID, centersVarID, centers),        &
              'nc_write_model_vars ','put_var centers ')
call nc_check(nf90_put_var(ncFileID, cornersVarID, corners),        &
              'nc_write_model_vars ','put_var centers ')
call nc_check(nf90_put_var(ncFileID, aVarID, a),    &
              'nc_write_model_vars ','put_var a ')
call nc_check(nf90_put_var(ncFileID, bVarID, b),    &
              'nc_write_model_vars ','put_var b ')
call nc_check(nf90_put_var(ncFileID, x_ax_bearingsVarID, x_ax_bearings),      &
              'nc_write_model_vars ','put_var x_ax_bearings ')

end function nc_write_cs_grid_file


   function nc_read_cs_grid_file() result (ierr)
!=======================================================================
!  function nc_read_cs_grid_file() result (ierr)

! Read the number of neighbors, corners, centers, a and b coefficients, and x_ax_bearings
! from a netCDF file once for this grid at the beginning of the assimilation.

integer :: ierr, ncFileID, ncFldID, nc_size, num_dims, max_nghbrs, shp(2)
character(len=NF90_MAX_NAME) :: nc_name

ierr = 0

! Open the cubed sphere grid relationships file
call nc_check(nf90_open(path = trim(cs_grid_file), mode = nf90_nowrite, ncid = ncFileID), &
      'nc_read_cs_grid_file', 'opening '//trim(cs_grid_file))

! learn how many dimensions are defined in this file.
call nc_check(nf90_inquire(ncFileID, num_dims), 'nc_read_cs_grid_file_size', 'inquire num_dims')

! Dimensions written out:
!              name="ncenters",      len = ncenters,      dimid = ncentersID), &
!              name="ncorners",      len = ncorners,      dimid = ncornersID), &
!              name="max_neighbors", len = max_neighbors, dimid = max_neighborsID), &
!              name="ncol",          len = ncol,          dimid = ncolID), &
!              name="ncoef_a",       len = 3,             dimid = aID), &
!              name="ncoef_b",       len = 2,             dimid = bID), &
call nc_check(nf90_inquire_dimension(ncFileID, 1, nc_name, ncenters), &
              'nc_read_cs_grid_file_size', 'inquire for '//trim(nc_name))

call nc_check(nf90_inquire_dimension(ncFileID, 2, nc_name, ncorners), &
              'nc_read_cs_grid_file_size', 'inquire for '//trim(nc_name))

call nc_check(nf90_inquire_dimension(ncFileID, 3, nc_name, max_nghbrs), &
              'nc_read_cs_grid_file_size', 'inquire for '//trim(nc_name))
if (trim(nc_name) /= 'max_neighbors') then
   ierr = 10
   write(string1, *) trim(cs_grid_file),' max_nghbrs does not match ', trim(model_config_file)
   call error_handler(E_ERR,'nc_read_cs_grid_file',string1,source,revision,revdate)
end if
! Check value against the namelist/parameter value.
if (max_nghbrs /= max_neighbors) then
   ierr = 15
   write(string1, *) trim(cs_grid_file),' max_nghbrs does not match max_neighbors', &
         max_nghbrs,max_neighbors
   call error_handler(E_ERR,'nc_read_cs_grid_file',string1,source,revision,revdate)
end if

call nc_check(nf90_inquire_dimension(ncFileID, 4, nc_name, nc_size), &
              'nc_read_cs_grid_file_size', 'inquire for '//trim(nc_name))
if (nc_size == ncol .and. trim(nc_name) == 'ncol') then
   allocate ( corners(ncenters,ncorners),  &
              num_nghbrs           (ncol), &
              centers(max_neighbors,ncol), &
              x_ax_bearings  (ncorners,ncenters), &
              a            (3,ncorners,ncenters), &
              b            (2,ncorners,ncenters), &
              x_planar     (3,ncorners,ncenters), &
              y_planar     (3,ncorners,ncenters)  )
! x_planar only dimensioned here for convenience.
! After debugging, it will only have dimension(3) and be a local variable where it's used.
   ! Initialize the grid variables
   num_nghbrs = MISSING_I; centers = MISSING_I;
   a = MISSING_R8;         b = MISSING_R8
   x_ax_bearings = MISSING_R8
   if (allocated(centers) .and. do_out) then
      shp = shape(centers)
      PRINT*,'Shape of centers = ',shp
   end if
   if (allocated(corners) .and. do_out) then
      shp = shape(corners)
      PRINT*,'Shape of corners = ',shp
   end if
else
   ierr = 20
   write(string1, *) trim(cs_grid_file),' ncol does not match ', trim(model_config_file)
   call error_handler(E_ERR,'nc_read_cs_grid_file',string1,source,revision,revdate)
end if

!if (print_details .and. do_out) write(*,*) 'Dims info = ',i, trim(dim_names(i)), dim_sizes(i)

call nc_check(nf90_inq_varid(ncFileID, 'num_nghbrs', ncFldID), &
                'nc_read_cs_grid_file', 'inq_varid num_nghbrs')
! if (print_details .and. do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
call nc_check(nf90_get_var(ncFileID, ncFldID, num_nghbrs ), &
                'nc_read_cs_grid_file', 'get_var num_nghbrs')

call nc_check(nf90_inq_varid(ncFileID, 'centers', ncFldID), &
                'nc_read_cs_grid_file', 'inq_varid centers')
! if (print_details .and. do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
call nc_check(nf90_get_var(ncFileID, ncFldID, centers ), &
                'nc_read_cs_grid_file', 'get_var centers')

call nc_check(nf90_inq_varid(ncFileID, 'corners', ncFldID), &
                'nc_read_cs_grid_file', 'inq_varid corners')
! if (print_details .and. do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
call nc_check(nf90_get_var(ncFileID, ncFldID, corners, &
                           start=(/1,1/),count=(/ncenters,ncorners/) ), &
                'nc_read_cs_grid_file', 'get_var corners')

call nc_check(nf90_inq_varid(ncFileID, 'a', ncFldID), &
                'nc_read_cs_grid_file', 'inq_varid a')
! if (print_details .and. do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
call nc_check(nf90_get_var(ncFileID, ncFldID, a ), &
                'nc_read_cs_grid_file', 'get_var a')

call nc_check(nf90_inq_varid(ncFileID, 'b', ncFldID), &
                'nc_read_cs_grid_file', 'inq_varid b')
! if (print_details .and. do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
call nc_check(nf90_get_var(ncFileID, ncFldID, b ), &
                'nc_read_cs_grid_file', 'get_var b')

call nc_check(nf90_inq_varid(ncFileID, 'x_ax_bearings', ncFldID), &
                'nc_read_cs_grid_file', 'inq_varid x_ax_bearings')
! if (print_details .and. do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
call nc_check(nf90_get_var(ncFileID, ncFldID, x_ax_bearings ), &
                'nc_read_cs_grid_file', 'get_var x_ax_bearings')

! Close the file
call nc_check(nf90_close(ncFileID), 'nc_read_cs_grid_file', 'closing '//trim(cs_grid_file))

end function nc_read_cs_grid_file


   real function bearing(lon1,lat1,lon2,lat2)
!=======================================================================
!  real function bearing(lon1,lat1,lon2,lat2)
!
! Calculate the direction along the great circle from point 1 on a sphere
! to point 2, relative to north.
! All inputs should have units of radians.
! Output is radians.
! From http://www.movable-type.co.uk/scripts/latlong.html

real(r8), intent(in)    :: lon1,lat1, lon2,lat2

real(r8) :: lon1c,lon2c, cos_lat2, del_lon

! Make sure the poles are handled consistently:
! If the longitude of the pole point is defined as 0.0,
! then the bearing to a nearby point will = the longitude of the point.
! This is consistent/continuous with the bearing from points extremely near
! the pole.
if (half_PI - (abs(lat1)) < epsilon(lat1)) then
   lon1c = 0._r8
else
   lon1c = lon1
end if
if (half_PI - (abs(lat2)) < epsilon(lat2)) then
   lon2c = 0._r8
else
   lon2c = lon2
end if

cos_lat2 = cos(lat2)
del_lon  = lon2c - lon1c
! Normalize del_lon to -pi<angle<pi.
del_lon = mod(del_lon,PI) - PI*int(del_lon/PI)
bearing = atan2(                                cos_lat2*sin(del_lon),  &
                cos(lat1)*sin(lat2) - sin(lat1)*cos_lat2*cos(del_lon) )
return

end function bearing


   subroutine nc_read_model_atts(att, att_vals, nflds, ncFileID)
!=======================================================================
!  subroutine nc_read_model_atts(att, att_vals, nflds, ncFileID)
!
! reads the value of an attribute for each of the fields in cflds.
!
! should be called with att = one of the attributes from the program variable
! input file, which will be written to the Posterior and Prior.nc files

!----------------------------------------------------------------------
! Local workspace
integer :: i, nchars, ierr
integer :: ncFileID, ncFldID, ncAttID, att_type

!----------------------------------------------------------------------
integer,                                intent(in)  :: nflds
character (len=*),                      intent(in)  :: att
character (len=128), dimension(nflds), intent(out)  :: att_vals

! open CAM 'initial' file
! DONE ALREADY in static_init_model
! call check(nf90_open(path = trim(model_config_file), mode = nf90_nowrite, &
!           ncid = ncFileID))

! read CAM 'initial' file attribute desired
if (print_details .and. do_out) PRINT*,'nc_read_model_atts: reading ',trim(att)
do i = 1,nflds
   call nc_check(nf90_inq_varid(ncFileID, trim(cflds(i)), ncFldID), 'nc_read_model_atts', &
                 'inq_varid '//trim(cflds(i)))

! could be inquire_attribute
!
   ierr = nf90_inquire_attribute(ncFileID, ncFldID, trim(att), att_type, nchars, ncAttID)

   if (ierr == nf90_noerr) then
      call nc_check(nf90_get_att(ncFileID, ncFldID, trim(att) ,att_vals(i) ), &
                    'nc_read_model_atts', 'get_att '//trim(att))
      att_vals(i)(nchars+1:128) = ' '
      if (print_details .and. do_out) WRITE(*,'(A,1X,I6,I6,1X,A,1X,A)') &
         att, ncFldID, nchars, cflds(i), trim(att_vals(i))
   else
      WRITE(*,*) ncFldID, cflds(i), 'NOT AVAILABLE'
   end if
end do

end subroutine nc_read_model_atts


   subroutine nc_read_global_att(ncFileID, att, att_val)
!=======================================================================
!  subroutine nc_read_global_att(ncFileID, att, att_val)
!
! Reads the value of a global attribute.

!----------------------------------------------------------------------
integer,           intent(in)  :: ncFileID
character (len=*), intent(in)  :: att
integer,           intent(out) :: att_val

!----------------------------------------------------------------------
! Local workspace
integer :: ierr
integer :: ncAttID, att_type, nchars

! NF90_GLOBAL is the psuedo-variable name used for global attributes.
ierr = nf90_inquire_attribute(ncFileID, NF90_GLOBAL, trim(att), &
       xtype=att_type, len=nchars, attnum=ncAttID)
!if (print_details .and. do_out)  &
!   PRINT*,'nc_read_global_att: passed nf90_inquire_attribute.  type, ID = ',att_type,ncAttID


if (ierr == nf90_noerr) then
   call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, trim(att) ,att_val), &
                 'nc_read_global_att', 'get_att '//trim(att))
   if (print_details .and. do_out) WRITE(*,'(A,I5,2A, I6)') &
      'nc_read_global_att for file ',ncFileID,' attribute and value = ',att, att_val
end if

return

end subroutine nc_read_global_att

   subroutine read_cam_coord(ncfileid, var, cfield)
!=======================================================================
! subroutine read_cam_coord(ncfileid, var , cfield)
!
! read CAM 'initial' file coordinate, i.e. 'lat     ','lon     ','gw      '
!           ,'hyai    ',...

!----------------------------------------------------------------------
integer,            intent(in)  :: ncfileid   ! file and field IDs
character (len=8),  intent(in)  :: cfield
type(grid_1d_type), intent(out) :: var

!----------------------------------------------------------------------
! Local workspace
integer :: i, coord_size             ! grid/array indices
integer :: ncfldid    ! file and field IDs
integer :: fld_exist  ! some grid fields don't exist on some CAM initial files (slat, slon, ...?)
integer :: ncerr      ! other nc errors; don't abort
integer, dimension(nf90_max_var_dims) :: coord_dims
! Some attributes are _Fillvalue (real) which I'll ignore for now.
! The following are used to repack the attributes I want into a compact form
integer :: num_atts, keep_atts, alen
integer                                     :: att_type
character (len=nf90_max_name)               :: att_name
character (len=nf90_max_name), allocatable  :: att_names(:)
character (len=128),           allocatable  :: att_vals(:)
real(r8)                                    :: resol, resol_1, resol_n


! read CAM 'initial' file field desired
fld_exist = nf90_inq_varid(ncfileid, trim(cfield), ncfldid)
if (fld_exist /= nf90_noerr ) then
   var%label = '        '
   return
end if

ncerr = nf90_inquire_variable(ncfileid, ncfldid, dimids = coord_dims, nAtts = num_atts)
if (ncerr /= nf90_noerr ) then
   write(*,*) 'Variable ',cfield,' dimids = ',coord_dims(1)
   write(*,*) nf90_strerror(ncerr)
   var%label = '        '
   var%dim_id = 0
   return
end if
if (coord_dims(1) == 0) then
   coord_size = 1                 ! to handle P0
else
   coord_size = dim_sizes(coord_dims(1))
end if

allocate(att_names(num_atts), att_vals(num_atts))

! get attributes
keep_atts = 0
do i=1,num_atts
   call nc_check(nf90_inq_attname(ncfileid, ncfldid, i, att_name), &
                 'read_cam_coord', 'inq_attname '//trim(att_name))

! CAM FV initial files have coordinates with attributes that are numerical, not character.
! (_FillValue).  These are not used because the coordinates are dimensioned exactly
! the right size.  I'll test for the type of att, and if it's not char, I'll ignore it
! and reduce the num_atts by 1.

! Otherwise I need a var@atts_type and separate var%atts_vals_YYY for each NetCDF
! external type (6 of them) I might run into.

   call nc_check(nf90_inquire_attribute(ncfileid, ncfldid, att_name, xtype=att_type,len=alen ), &
                 'read_cam_coord', 'inquire_attribute '//trim(att_name))

   if (att_type == nf90_char) then
      keep_atts = keep_atts + 1
      att_names(keep_atts) = att_name
      call nc_check(nf90_get_att(ncfileid, ncfldid, att_name, att_vals(keep_atts)), &
                    'read_cam_coord', 'get_att '//trim(att_name) )
      att_vals(keep_atts)(alen+1:128) = ' '

   else
      if (do_out) write(*,*) '                ignoring attribute ',trim(att_name),    &
                    ' because it is not a character type'
   end if
end do

! allocate space for this grid structure, and put vector length and # attributes in structure
call init_grid_1d_instance(var, coord_size, keep_atts)

var%label = cfield
var%dim_id = coord_dims(1)

do i = 1,keep_atts
   var%atts_names(i) = att_names(i)
   var%atts_vals(i)  = att_vals(i)
end do

! call check(nf90_get_var(ncfileid, ncfldid, var%vals(1:coord_size) ,start=(/1/) &
call nc_check(nf90_get_var(ncfileid, ncfldid, var%vals, start=(/1/) &
    ,count=(/coord_size/) ), 'read_cam_coord' ,'get_var '//cfield)

! Determine whether coordinate is regularly spaced,
! If so, store the coordinate resolution in the grid_1d_type.
if (cfield(1:2) == 'hy') then
   var%resolution = MISSING_R8
else
   resol_1 = var%vals(2) - var%vals(1)
   ! CS; resol_1 can be 0.
   if (resol_1 /= 0._r8) then
      resol = 1.0_r8/resol_1
      var%resolution = resol_1

      ! Test all other coordinate spacings.  If any of them differ from the first
      ! by more than epsilon (smallest meaningful number relative to the coordinate spacings)
      ! then spacing is irregular.
      Res: do i = 3,coord_size
         resol_n = var%vals(i) - var%vals(i-1)
         if (((resol_n - resol_1) *resol) > epsilon(resol_n)) then
            var%resolution = -1.0_r8
            exit Res
         end if
      end do Res
   else
      var%resolution = -1.0_r8
   end if
   ! CS end resol_1 change
end if

if (print_details .and. do_out) then
   write(*,'(3A,I6,A,I8,A,1pE12.4)')  'reading ',cfield,' using id ',ncfldid,  &
          ' size ',coord_size,' resolution ', var%resolution
   WRITE(*,*) 'first, last val: ', var%vals(1),var%vals(coord_size)
   ! to get entire array, in case of debugging, use this instead:
   ! WRITE(*,*) (var%vals(i),i=1,coord_size)
end if

deallocate(att_names, att_vals)

end subroutine read_cam_coord



   subroutine init_grid_1d_instance(var, length, num_atts)
!=======================================================================
! subroutine init_grid_1d_instance(var)
!
! Initializes an instance of a cam grid variable

type(grid_1d_type), intent(out) :: var
integer,            intent(in ) :: length, num_atts

! Initialize the storage space and return
allocate( var%vals      (length))
var%vals = 0.0_r8
allocate( var%atts_names(num_atts))
allocate( var%atts_vals (num_atts))

var%length = length
var%num_atts = num_atts

end subroutine init_grid_1d_instance


   subroutine end_grid_1d_instance(var)
!=======================================================================
! subroutine end_grid_1d_instance(var)
!
! Ends an instance of a cam grid_1d variable

type(grid_1d_type), intent(inout) :: var

deallocate(var%vals, var%atts_names, var%atts_vals)

end subroutine end_grid_1d_instance


   subroutine order_state_fields(cflds,nflds)
!=======================================================================
! subroutine order_state_fields(cflds,nflds)
!
! fills cflds with state_names for use in I/O of caminput.nc
! Could eventually tally the number of each kind of field; 2D,3D
! and compare each entry against a master list.
! Sort by class of variable too? So user could provide one unordered list?
! Also assigns TYPE_s for use various routines.

integer :: i, i1, nfld
integer, intent(in) :: nflds
character (len = *), dimension(nflds), intent(out) :: cflds

allocate (TYPE_1D(state_num_1d),TYPE_2D(state_num_2d),TYPE_3D(state_num_3d))
! kdr where should these be deallocated?


nfld = 0

! 0D fields
do i=1,state_num_0d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_0d(i)
end do

! 1D fields
do i=1,state_num_1d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_1d(i)
   TYPE_1D(i) = nfld
!  CS
   if (state_names_1d(i) == 'PS      ') TYPE_PS    = nfld
   if (state_names_1d(i) == 'EFGWORO ') TYPE_EFGWORO = nfld
   if (state_names_1d(i) == 'FRACLDV ') TYPE_FRACLDV = nfld
   if (state_names_1d(i) == 'TBOT    ') TYPE_TBOT  = nfld
   if (state_names_1d(i) == 'TS      ') TYPE_TS    = nfld
   if (state_names_1d(i) == 'TSOCN   ') TYPE_TSOCN = nfld
end do

! 2D fields
do i=1,state_num_2d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_2d(i)
   TYPE_2D(i) = nfld
   if (state_names_2d(i) == 'PS      ') TYPE_PS    = nfld
   if (state_names_2d(i) == 'EFGWORO ') TYPE_EFGWORO = nfld
   if (state_names_2d(i) == 'FRACLDV ') TYPE_FRACLDV = nfld
   if (state_names_2d(i) == 'TBOT    ') TYPE_TBOT  = nfld
   if (state_names_2d(i) == 'TS      ') TYPE_TS    = nfld
   if (state_names_2d(i) == 'TSOCN   ') TYPE_TSOCN = nfld
!  cubed sphere
   if (state_names_2d(i) == 'T       ') TYPE_T = nfld
   if (state_names_2d(i) == 'U       ') TYPE_U = nfld
   if (state_names_2d(i) == 'V       ') TYPE_V = nfld
   ! if (state_names_2d(i) == 'US      ') TYPE_U = nfld
   ! if (state_names_2d(i) == 'VS      ') TYPE_V = nfld
   if (state_names_2d(i) == 'Q       ') TYPE_Q = nfld
   if (state_names_2d(i) == 'CLDICE  ') TYPE_CLDICE = nfld
   if (state_names_2d(i) == 'CLDLIQ  ') TYPE_CLDLIQ = nfld
   if (state_names_2d(i) == 'LCWAT   ') TYPE_LCWAT  = nfld
   if (state_names_2d(i) == 'QCWAT   ') TYPE_QCWAT  = nfld
   ! chem
   if (state_names_2d(i) == 'CO      ') TYPE_CO  = nfld
   if (state_names_2d(i) == 'CO2     ') TYPE_CO2 = nfld
   if (state_names_2d(i) == 'NO      ') TYPE_NO  = nfld
   if (state_names_2d(i) == 'NO2     ') TYPE_NO2 = nfld
   if (state_names_2d(i) == 'CH4     ') TYPE_CH4 = nfld
   if (state_names_2d(i) == 'NH3     ') TYPE_NH3 = nfld
   if (state_names_2d(i) == 'O       ') TYPE_O   = nfld
   if (state_names_2d(i) == 'O3      ') TYPE_O3  = nfld
end do

! 3D fields (including q)
do i=1,state_num_3d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_3d(i)
   TYPE_3D(i) = nfld
   if (state_names_3d(i) == 'T       ') TYPE_T = nfld
   if (state_names_3d(i) == 'U       ') TYPE_U = nfld
   if (state_names_3d(i) == 'V       ') TYPE_V = nfld
   if (state_names_3d(i) == 'US      ') TYPE_U = nfld
   if (state_names_3d(i) == 'VS      ') TYPE_V = nfld
   if (state_names_3d(i) == 'Q       ') TYPE_Q = nfld
   if (state_names_3d(i) == 'CLDICE  ') TYPE_CLDICE = nfld
   if (state_names_3d(i) == 'CLDLIQ  ') TYPE_CLDLIQ = nfld
   if (state_names_3d(i) == 'LCWAT   ') TYPE_LCWAT  = nfld
   if (state_names_3d(i) == 'QCWAT   ') TYPE_QCWAT  = nfld
   ! chem
   if (state_names_3d(i) == 'CO      ') TYPE_CO  = nfld
   if (state_names_3d(i) == 'CO2     ') TYPE_CO2 = nfld
   if (state_names_3d(i) == 'NO      ') TYPE_NO  = nfld
   if (state_names_3d(i) == 'NO2     ') TYPE_NO2 = nfld
   if (state_names_3d(i) == 'CH4     ') TYPE_CH4 = nfld
   if (state_names_3d(i) == 'NH3     ') TYPE_NH3 = nfld
   if (state_names_3d(i) == 'O       ') TYPE_O   = nfld
   if (state_names_3d(i) == 'O3      ') TYPE_O3  = nfld

end do

if (nfld .ne. nflds) then
   write(string1, *) 'nfld = ',nfld,', nflds = ',nflds,' must be equal '
   call error_handler(E_ERR, 'order_state_fields', string1, source, revision, revdate)
else if (do_out) then
   if (print_details) then
      write(logfileunit,'(/A/)') 'State vector is composed of these fields: '
   !   write(logfileunit,'((8(A8,1X)))') (cflds(i),i=1,nflds)
      do i=1,state_num_0d
         write(logfileunit,'(/A,I4)') cflds(i), TYPE_1D(i)
      end do
      i1 = state_num_0d
      do i=1,state_num_1d
         write(logfileunit,'(/A,I4)') cflds(i1+i), TYPE_1D(i)
      end do
      i1 = i1 + state_num_1d
      do i=1,state_num_2d
         write(logfileunit,'(/A,I4)') cflds(i1+i), TYPE_2D(i)
      end do
      i1 = i1 + state_num_2d
      do i=1,state_num_3d
         write(logfileunit,'(/A,I4)') cflds(i1+i), TYPE_3D(i)
      end do
      write(logfileunit,'(/A)')        'TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q, TYPE_CLDLIQ, TYPE_CLDICE = '
      write(logfileunit,'((8(I8,1X)))') TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q, TYPE_CLDLIQ, TYPE_CLDICE
! chem; add species 'CO', 'CO2', 'NO', 'NO2', 'CH4', 'NH3', 'O', 'O3 '
      write(logfileunit,'(/A)')        'TYPE_CO, TYPE_CO2, TYPE_NO, TYPE_NO2, TYPE_CH4, TYPE_NH3, TYPE_O, TYPE_O3 = '
      write(logfileunit,'((8(I8,1X)))') TYPE_CO, TYPE_CO2, TYPE_NO, TYPE_NO2, TYPE_CH4, TYPE_NH3, TYPE_O, TYPE_O3
   else
      call error_handler(E_MSG, 'order_state_fields', 'State vector is composed of these fields: ')
      do i = 1,nflds
         call error_handler(E_MSG, 'order_state_fields', trim(cflds(i)))
      end do
   end if
end if

return

end subroutine order_state_fields


   subroutine map_kinds()
!=======================================================================
! subroutine map_kinds()

! ? Should this be a function instead; removes need to dimension obs_loc_in arbitrarily
!   and wastefully.  But then it's called millions of times, instead of accessing an
!   array that's defined once.

! Makes an array of 'locations within the state vector' of the obs kinds
! that come from obs_kind_mod, which we anticipate CAM's model_mod will need.
! The obs kind that's needed will be the index into this array,
! the corresponding value will be the position of that field (not individual variable)
! within the state vector according to state_name_Xd.
! This subroutine will be called from static_init_model, so it will not have to be
! recomputed for every obs.
! Also maps the local model_mod TYPE_s onto the DART QTY_s by the same mechanism.

! other QTY_ possibilities are listed after the 'use obs_kind_mod' statement

integer :: i

! Physically 2D fields
dart_to_cam_kinds(QTY_SURFACE_PRESSURE) = TYPE_PS
if (TYPE_PS /= MISSING_I) cam_to_dart_kinds(TYPE_PS) = QTY_SURFACE_PRESSURE

dart_to_cam_kinds(QTY_GRAV_WAVE_DRAG_EFFIC) = TYPE_EFGWORO
if (TYPE_EFGWORO /= MISSING_I) &
   cam_to_dart_kinds(TYPE_EFGWORO) = QTY_GRAV_WAVE_DRAG_EFFIC

dart_to_cam_kinds(QTY_GRAV_WAVE_STRESS_FRACTION) = TYPE_FRACLDV
if (TYPE_FRACLDV /= MISSING_I) &
   cam_to_dart_kinds(TYPE_FRACLDV) = QTY_GRAV_WAVE_STRESS_FRACTION

! dart_to_cam_kinds(QTY_SURFACE_TEMPERATURE  ?  ) = TYPE_TS
! dart_to_cam_kinds(QTY_SEA_SURFACE_TEMPERATURE  ?  ) = TYPE_TSOCN

! Physically 3D fields
dart_to_cam_kinds(QTY_TEMPERATURE)        = TYPE_T
dart_to_cam_kinds(QTY_U_WIND_COMPONENT)   = TYPE_U
dart_to_cam_kinds(QTY_V_WIND_COMPONENT)   = TYPE_V
dart_to_cam_kinds(QTY_SPECIFIC_HUMIDITY)  = TYPE_Q
dart_to_cam_kinds(QTY_CLOUD_LIQUID_WATER) = TYPE_CLDLIQ
dart_to_cam_kinds(QTY_CLOUD_ICE)          = TYPE_CLDICE
! dart_to_cam_kinds(QTY_CLOUD_WATER  ?  ) = TYPE_LCWAT
dart_to_cam_kinds(QTY_CO)  = TYPE_CO
dart_to_cam_kinds(QTY_CO2) = TYPE_CO2
dart_to_cam_kinds(QTY_NO)  = TYPE_NO
dart_to_cam_kinds(QTY_NO2) = TYPE_NO2
dart_to_cam_kinds(QTY_CH4) = TYPE_CH4
dart_to_cam_kinds(QTY_NH3) = TYPE_NH3
! dart_to_cam_kinds(QTY_O)   = TYPE_O
dart_to_cam_kinds(QTY_O3)  = TYPE_O3

if (TYPE_T      /= MISSING_I) cam_to_dart_kinds(TYPE_T)      = QTY_TEMPERATURE
if (TYPE_U      /= MISSING_I) cam_to_dart_kinds(TYPE_U)      = QTY_U_WIND_COMPONENT
if (TYPE_V      /= MISSING_I) cam_to_dart_kinds(TYPE_V)      = QTY_V_WIND_COMPONENT
if (TYPE_Q      /= MISSING_I) cam_to_dart_kinds(TYPE_Q)      = QTY_SPECIFIC_HUMIDITY
if (TYPE_CLDLIQ /= MISSING_I) cam_to_dart_kinds(TYPE_CLDLIQ) = QTY_CLOUD_LIQUID_WATER
if (TYPE_CLDICE /= MISSING_I) cam_to_dart_kinds(TYPE_CLDICE) = QTY_CLOUD_ICE
! cam_to_dart_kinds(TYPE_LCWAT) = QTY_CLOUD_WATER  ?
if (TYPE_CO     /= MISSING_I) cam_to_dart_kinds(TYPE_CO)  = QTY_CO
if (TYPE_CO2    /= MISSING_I) cam_to_dart_kinds(TYPE_CO2) = QTY_CO2
if (TYPE_NO     /= MISSING_I) cam_to_dart_kinds(TYPE_NO)  = QTY_NO
if (TYPE_NO2    /= MISSING_I) cam_to_dart_kinds(TYPE_NO2) = QTY_NO2
if (TYPE_CH4    /= MISSING_I) cam_to_dart_kinds(TYPE_CH4) = QTY_CH4
if (TYPE_NH3    /= MISSING_I) cam_to_dart_kinds(TYPE_NH3) = QTY_NH3
! if (TYPE_O  /= MISSING_I) cam_to_dart_kinds(TYPE_O)   = QTY_O
if (TYPE_O3     /= MISSING_I) cam_to_dart_kinds(TYPE_O3)  = QTY_O3


if (print_details .and. do_out) then
   write(*,*) 'OBS_KIND   FIELD_TYPE'
   do i=1,100
      if (dart_to_cam_kinds(i) /= MISSING_I) write(*,'(2I8)') i, dart_to_cam_kinds(i)
   end do
end if

! In the future, if fields are not ordered nicely, or if users are specifying
! correspondence of obs fields with state fields, I may want code like:
! The max size of QTY_ should come from obs_kind_mod
! do i=1,state_num_3d
!    if (state_names_3d(i)(1:1) == 'T' .and. &
!        QTY_TEMPERATURE <= 100) ) dart_to_cam_kinds(QTY_TEMPERATURE) = TYPE_3D(i)
! end do

return

end subroutine map_kinds


   subroutine fill_gc
!=======================================================================
!  subroutine fill_gc
!
! Subroutine to generate location_types of the cubed sphere grid
! and put them into get_close_type cs_gc, with other derived components.

integer :: c

allocate (cs_locs(ncol), cs_kinds(ncol))

! CS inputs in degrees.
do c=1,ncol
   cs_locs(c)  = set_location(lon%vals(c), lat%vals(c), MISSING_R8, VERTISUNDEF)
   cs_kinds(c) = 0
end do
! if (print_details) then
!    write(*,'(A)') 'cs_locs(1:ncol:10000) passed to get_close_obs_init'
!    do c=1,ncol,10000
!       write(*,'(1p3e15.7)') get_location(cs_locs(c))
!    end do
! end if

! Initialize cs_gc%maxdist using the maximum grid spacing.
! cs_gc is only output.  coarse_grid is only input.
! Using the nominal grid resolution is plenty large.
! There will always be at least 2 nodes within 1 coarse_grid in all directions.
call get_close_maxdist_init(cs_gc, coarse_grid)

! Use cs_gc%maxdist and node locations to define the rest of cs_gc.
call get_close_obs_init(cs_gc, ncol, cs_locs)

return

end subroutine fill_gc

! End of static_init_model section
!#######################################################################

! Module I/O to/from DART and files

   subroutine read_cam_init(file_name, var, model_time)
!=======================================================================
! subroutine read_cam_init(file_name, var, model_time)
!

character(len = *),        intent(in)    :: file_name
type(model_type),          intent(inout) :: var
type(time_type), optional, intent(out)   :: model_time

! Local workspace
integer :: i, k, n, m, ifld  ! grid and constituent indices
integer :: ncfileid, ncfldid, dimid, varid, dimlen
real(r8), allocatable :: temp_3d(:,:,:), temp_2d(:,:)

integer :: iyear, imonth, iday, ihour, imin, isec, rem
integer, allocatable, dimension(:) :: datetmp, datesec

!----------------------------------------------------------------------
! read CAM 'initial' file domain info
call nc_check(nf90_open(path = trim(file_name), mode = nf90_nowrite, ncid = ncfileid), &
      'read_cam_init', 'opening '//trim(file_name))

! The temporary arrays into which fields are read are dimensioned by the largest values of
! the sizes of the dimensions listed in f_dim_RANKd
! f_dim_max contents assume that time is always the last dimension on NetCDF files,
! so f_dim_max(4,3) and f_dim_max(3,2) are the non-spatial dimensions to ignore here.
allocate (temp_3d(f_dim_max(1,3),f_dim_max(2,3),f_dim_max(3,3)) &
         ,temp_2d(f_dim_max(1,2),f_dim_max(2,2)) )

! read CAM 'initial' file fields desired

ifld = 0
!0d fields; scalars are recognized and handled differently than vectors
!           by NetCDF
! kdr; nf90_put_var was probably bombing because
!      netcdf recognized that it was writing a scalar, called an f77
!      scalar put_var, and then choked when it saw the count = ...
do i= 1, state_num_0d
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid), &
                'read_cam_init', 'inq_varid '//trim(cflds(ifld)))
   if (print_details .and. do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 1D; TIME(=1)
   call nc_check(nf90_get_var(ncfileid, ncfldid, var%vars_0d(i) ), &
                'read_cam_init', 'get_var '//trim(cflds(ifld)))
end do


!1d fields
do i= 1, state_num_1d
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid), &
                'read_cam_init', 'inq_varid '//trim(cflds(ifld)))
! debug; remove next two lines when I'm sure I don't need them.  Also for 2d and 3d
! done in trans_coord already
!   call check(nf90_inquire_variable(ncfileid, ncfldid, dimids=f_dimid_1d(1:2,i)))
   if (print_details .and. do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid

   ! s_dim_1d should = f_dim_1d
   call nc_check(nf90_get_var(ncfileid, ncfldid, var%vars_1d(1:s_dim_1d(i), i) &
             ,start=(/1,1/) ,count=(/f_dim_1d(1,i), 1/) ), &
             'read_cam_init', 'get_var '//trim(cflds(ifld)))
end do

! CS; simplify by removing old CAM coordinate order?
!2d fields
do i= 1, state_num_2d
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid), &
                'read_cam_init', 'inq_varid '//trim(cflds(ifld)))
! done in trans_coord already
!   call check(nf90_inquire_variable(ncfileid, ncfldid, dimids=f_dimid_2d(1:3,i)))
   if (print_details .and. do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
   ! fields on file are 3D; lon, lat, (usually) and then  TIME(=1)
   ! Need to use temp_Nd here too; I am coding for not knowing what 2 of the 3 dimensions the
   ! 2d fields will have.
   call nc_check(nf90_get_var(ncfileid, ncfldid, temp_2d(1:f_dim_2d(1,i), 1:f_dim_2d(2,i))  &
                              ,start=(/1,1,1/)    ,count=(/f_dim_2d(1,i),   f_dim_2d(2,i), 1/) ), &
                'read_cam_init', 'get_var '//trim(cflds(ifld)))
   if (s_dim_2d(1,i) == f_dim_2d(1,i)) then
      var%vars_2d(1:s_dim_2d(1,i), 1:s_dim_2d(2,i),i) = &
          temp_2d(1:f_dim_2d(1,i), 1:f_dim_2d(2,i)  )
   else if (s_dim_2d(1,i) == f_dim_2d(2,i)) then
      do k=1,s_dim_2d(1,i)
      do m=1,s_dim_2d(2,i)   ! first temp dim is inner loop for faster reads
         var%vars_2d(k,m,i) = temp_2d(m,k)
      end do
      end do
   else
      ! error
   end if
end do

! CS; simplify by removing old CAM coordinate order?
! 3d fields
do i=1, state_num_3d
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid), &
                'read_cam_init', 'inq_varid '//trim(cflds(ifld)))
! done in trans_coord already
!   call nc_check(nf90_inquire_variable(ncfileid, ncfldid, dimids=f_dimid_3d(1:4,i)), &
!                 'read_cam_init')
   if (print_details .and. do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 4D; lon, lev, lat, TIME(=1)
!                     or; lon, lat, lev, TIME
   call nc_check(nf90_get_var(ncfileid, ncfldid                                                 &
             ,temp_3d(1:f_dim_3d(1,i), 1:f_dim_3d(2,i), 1:f_dim_3d(3,i)) ,start=(/1,1,1,1/)  &
             ,count=(/  f_dim_3d(1,i),   f_dim_3d(2,i),   f_dim_3d(3,i),1 /) ), &
             'read_cam_init', 'get_var '//trim(cflds(ifld)))
! read into dummy variable, then repackage depending on coord_order
! put it into lev,lon,lat order.
   if (coord_order == 1) then
!     lon,lev,lat as in original CAM
      do n=1,s_dim_3d(3,i)  ! lats
      do k=1,s_dim_3d(1,i)  ! levs
      do m=1,s_dim_3d(2,i)  ! lons/first file dim   put in inner loop for faster reads
         var%vars_3d(k,m,n, i) = temp_3d(m,k,n)
      end do
      end do
      end do
   else if (coord_order == 2) then
!     lon,lat,lev as in  CAM > 3.0.3
!      do n=1,s_dim_3d(2,i)  ! lats
!      do k=1,s_dim_3d(1,i)  ! levs
!      do m=1,s_dim_3d(3,i)  ! lons/first file dim   put in inner loop for faster reads
      do n=1,s_dim_3d(3,i)  ! lats
      do k=1,s_dim_3d(1,i)  ! levs
      do m=1,s_dim_3d(2,i)  ! lons/first file dim   put in inner loop for faster reads
         var%vars_3d(k,m,n, i) = temp_3d(m,n,k)
      end do
      end do
      end do
   end if
end do

! Read the time of the current state.
! All the caminput.nc files I have seen have two variables of
! length 'time' (the unlimited dimension): date, datesec
! The rest of the routine presumes there is but one time in the file -
! print warning message if this is not the case.

if (present( model_time)) then

   call nc_check(nf90_inq_dimid(ncfileid, 'time', dimid), &
          'read_cam_init', 'inq_dimid time '//trim(file_name))
   call nc_check(nf90_inquire_dimension(ncfileid, dimid, len=dimlen), &
          'read_cam_init', 'inquire_dimension time '//trim(file_name))

   if (dimlen /= 1) then
       write(string1,*)'UNUSUAL - ',trim(file_name),' has',dimlen,'times. Expected 1.'
       call error_handler(E_MSG, 'read_cam_init', string1, source, revision, revdate)
   end if

   allocate(datetmp(dimlen), datesec(dimlen))

   call nc_check(nf90_inq_varid(ncfileid, 'date', varid), &
          'read_cam_init', 'inq_varid date '//trim(file_name))
   call nc_check(nf90_get_var(ncfileid, varid, values=datetmp), &
          'read_cam_init', 'get_var date '//trim(file_name))

   call nc_check(nf90_inq_varid(ncfileid, 'datesec', varid), &
          'read_cam_init', 'inq_varid datesec '//trim(file_name))
   call nc_check(nf90_get_var(ncfileid, varid, values=datesec), &
          'read_cam_init', 'get_var datesec '//trim(file_name))

   ! The 'date' is YYYYMMDD ... datesec is 'current seconds of current day'
   iyear  = datetmp(dimlen) / 10000
   rem    = datetmp(dimlen) - iyear*10000
   imonth = rem / 100
   iday   = rem - imonth*100

   ihour  = datesec(dimlen) / 3600
   rem    = datesec(dimlen) - ihour*3600
   imin   = rem / 60
   isec   = rem - imin*60

   ! some cam files are from before the start of the gregorian calendar.
   ! since these are 'arbitrary' years, just change the offset.

   if (iyear < 1601) then
      write(logfileunit,*)' '
      write(     *     ,*)' '
      write(string1,*)'WARNING - ',trim(file_name),' changing year from ',iyear,'to',iyear+1601
      call error_handler(E_MSG, 'read_cam_init', string1, source, revision, &
                   revdate, text2='to make it a valid Gregorian date.')
      write(logfileunit,*)' '
      write(     *     ,*)' '
      iyear = iyear + 1601
   end if

   model_time = set_date(iyear,imonth,iday,ihour,imin,isec)

   if (do_out) then
      call print_date(model_time,' read_cam_init ... input date')
      call print_time(model_time,' read_cam_init ... input time')
      call print_date(model_time,' read_cam_init ... input date',logfileunit)
      call print_time(model_time,' read_cam_init ... input time',logfileunit)
   end if

   deallocate(datetmp, datesec)

end if

call nc_check(nf90_close(ncfileid), 'read_cam_init', 'closing '//trim(file_name))

deallocate (temp_3d,temp_2d)

end subroutine read_cam_init



   subroutine write_cam_coord_def(ncFileID, c_name, coord, dim_id, c_id)
!=======================================================================================
! subroutine write_cam_coord_def(ncFileID, c_name, coord, dim_id, c_id)

character (len=8),  intent(in)  :: c_name
integer,            intent(in)  :: ncFileID, dim_id
type(grid_1d_type), intent(in)  :: coord
integer,            intent(out) :: c_id

integer  :: i
!integer  :: nch

call nc_check(nf90_def_var(ncFileID, name=c_name, xtype=nf90_double, dimids=dim_id, &
                        varid=c_id), 'write_cam_coord_def', 'def_var '//trim(c_name))
!if (print_details .and. do_out) write(*,'(/A,A)') 'write_cam_coord_def;  ', trim(c_name)

do i=1,coord%num_atts
!   if (print_details .and. do_out) then
!!      nch = len_trim(coord%atts_vals(i))
!!                 i,trim(coord%atts_names(i)),' ', coord%atts_vals(i)(1:nch)
!      write(*,*) '   i, att_name, att_val', &
!                 i,trim(coord%atts_names(i)),' ', trim(coord%atts_vals(i))
!   end if
   call nc_check(nf90_put_att(ncFileID, c_id, coord%atts_names(i), coord%atts_vals(i)), &
                 'write_cam_coord_def', 'put_att '//trim(coord%atts_names(i)))
end do

return

end subroutine write_cam_coord_def

   subroutine write_cam_init(file_name, var, model_time)
!=======================================================================
! subroutine write_cam_init(file_name, var, model_time)

! write CAM 'initial' file fields that have been updated

character (len = *), intent(in)           :: file_name
type(model_type),    intent(inout)        :: var
type(time_type),     intent(in), optional :: model_time

integer               :: i, k, n, m, ifld, ncfileid, ncfldid, f_dim1, f_dim2
integer               :: iyear, imonth, iday, ihour, imin, isec
integer               :: dimid, dimlen, varid
integer, allocatable, dimension(:) :: datetmp, datesec
real(r8), allocatable :: temp_3d(:,:,:), temp_2d(:,:)

! Read CAM 'initial' file domain info
call nc_check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid), &
           'write_cam_init', 'opening '//trim(file_name))

! The temporary arrays into which fields are read are dimensioned by the largest values of
! the sizes of the dimensions listed in coord_RANKd
! debug
! ! ! This may not work for writing fields that have smaller sizes!
!     Unless array section is explicitly given in array argument
allocate (temp_3d(f_dim_max(1,3),f_dim_max(2,3),f_dim_max(3,3)))
allocate (temp_2d(f_dim_max(1,2),f_dim_max(2,2)))


if (print_details .and. do_out) write(*,*) 'write_cam_init; f_dim_max(:2) = ',f_dim_max(1,2),f_dim_max(2,2)

ifld = 0
! 0d fields are first
! kdr; Tim says that nf90_put_var was probably bombing because
!      netcdf recognized that it was writing a scalar, called an f77
!      scalar put_var, and then choked when it saw the count = ...
!      So, this padding is probably not necessary.
do i = 1, state_num_0d
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid), &
                 'write_cam_init', 'put_var '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(ncfileid, ncfldid, var%vars_0d(i) ), &
                 'write_cam_init', 'put_var '//trim(cflds(ifld)))
end do

! 1d fields
do i = 1, state_num_1d
   ! CS added this from 2d loop below.
   ! special code:  check and error out if the PS field has gone negative
   if (state_names_1d(i) == 'PS') then
      if (minval(var%vars_1d(:,i)) < 0._r8) then
         write(string1, *)'PS has negative values; should not happen'
         call error_handler(E_ERR, 'write_cam_init', string1, source, revision, revdate)
      end if
   end if
   ifld = ifld + 1
   call nc_check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid), &
                 'write_cam_init', 'inq_var '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(ncfileid, ncfldid, var%vars_1d(1:s_dim_1d(i), i), &
                                     start=(/1, 1/), count = (/s_dim_1d(i), 1/)),  &
                 'write_cam_init', 'inq_var '//trim(cflds(ifld)))
end do

! CS; simplify by removing old CAM coordinate order?

! 2d fields ; tricky because coordinates may have been rearranged to put them in the order
! of (lev, lon, lat) choosing only 2.  The original coordinate order is in f_dimid_2d.
! debug
! if (print_details .and. do_out) write(*,'(A/A)') 'write_cam_init 2D',' i f_dim1 f_dim2'

do i = 1, state_num_2d
   ! special code:  check and error out if the PS field has gone negative
   if (state_names_2d(i) == 'PS') then
      if (minval(var%vars_2d(:,:,i)) < 0._r8) then
         write(string1, *)'PS has negative values; should not happen'
         call error_handler(E_ERR, 'write_cam_init', string1, source, revision, revdate)
      end if
   ! CS: this section copied from 3d below, and converted to 2d.
   else if (state_names_2d(i) == 'Q') then
      where (var%vars_2d(:,:,i) < 1.e-12_r8) var%vars_2d(:,:,i) = 1.e-12_r8
   else if (state_names_2d(i) == 'CLDLIQ' .or. &
            state_names_2d(i) == 'CLDICE') then
      where (var%vars_2d(:,:,i) < 0._r8)     var%vars_2d(:,:,i) = 0._r8
   else if (state_names_2d(i) == 'T') then
      if (minval(var%vars_2d(:,:,i)) < 0._r8) then
         write(string1, *)'T has negative values; should not happen'
         call error_handler(E_ERR, 'write_cam_init', string1, source, revision, revdate)
      end if
   end if
   if (f_dimid_2d(1,i) == s_dimid_2d(1,i)) then
      ! model_mod and caminput store this variable with the same coordinate order.

      f_dim1 = s_dim_2d(1,i)
      f_dim2 = s_dim_2d(2,i)
      do n=1,f_dim2
      do m=1,f_dim1
         temp_2d(m,n) = var%vars_2d(m,n,i)
      end do
      end do
! debug
!      if (print_details .and. do_out) write(*,'(3I3,2F14.6)') i, f_dim1, f_dim2  , temp_2d(1,1), temp_2d(f_dim1, f_dim2)
   else if (f_dimid_2d(1,i) == s_dimid_2d(2,i)) then
      ! model_mod and caminput store this variable with transposed coordinate order.
      f_dim1 = s_dim_2d(2,i)
      f_dim2 = s_dim_2d(1,i)
      do m=1,f_dim1
      do n=1,f_dim2
         temp_2d(m,n) = var%vars_2d(n,m,i)
      end do
      end do
   else
      ! There aren't any more to try
   end if

   ifld = ifld + 1
   call nc_check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid),           &
                 'write_cam_init','inq_varid '//trim(cflds(ifld)))
!   writing part of temp_2d defined by start and count; don't need indices.
!   call nc_check(nf90_put_var(ncfileid, ncfldid, temp_2d(1:f_dim1, 1:f_dim2),    &
   call nc_check(nf90_put_var(ncfileid, ncfldid, temp_2d(1:f_dim1,1:f_dim2),     &
                           start=(/1, 1, 1/), count = (/f_dim1, f_dim2, 1/)),    &
                 'write_cam_init','put_var '//trim(cflds(ifld)))
!   call nc_check(nf90_put_var(ncfileid, ncfldid, var%vars_2d(1:s_dim_2d(1,i), 1:s_dim_2d(1,i), i),&
!                               start=(/1, 1, 1/), count = (/s_dim_2d(1,i),   s_dim_2d(2,i), 1/)), &
!                 'write_cam_init')
end do

! 3d fields; all 3 coordinates are present, and the order for model_mod fields is always the same.
do i = 1, state_num_3d
   ! special code:  set a minimum threshold for certain variables
   if (state_names_3d(i) == 'Q') then
      where (var%vars_3d(:,:,:,i) < 1.e-12_r8) var%vars_3d(:,:,:,i) = 1.e-12_r8
   else if (state_names_3d(i) == 'CLDLIQ' .or. &
            state_names_3d(i) == 'CLDICE') then
      where (var%vars_3d(:,:,:,i) < 0._r8)     var%vars_3d(:,:,:,i) = 0._r8
   else if (state_names_3d(i) == 'T') then
      if (minval(var%vars_3d(:,:,:,i)) < 0._r8) then
         write(string1, *)'T has negative values; should not happen'
         call error_handler(E_ERR, 'write_cam_init', string1, source, revision, revdate)
      end if
   end if

! CS; simplify by removing old CAM coordinate order?
   !  Repackage depending on coord_order then write the dummy variable.
   if (coord_order == 1) then
      ! lon,lev,lat as in original CAM
      do n=1,s_dim_3d(3,i)      ! lats
      do m=1,s_dim_3d(2,i)      ! lons
      do k=1,s_dim_3d(1,i)      ! levs
         temp_3d(m,k,n) = var%vars_3d(k,m,n,i)
      end do
      end do
      end do
      ! temp_3d         (1:s_dim_3d(1,i), 1:s_dim_3d(2,i), 1:s_dim_3d(3,i))      &
      !   = var%vars_3d(1:s_dim_3d(1,i), 1:s_dim_3d(2,i), 1:s_dim_3d(3,i), i)
   else if (coord_order == 2) then
      ! lon,lat,lev as in new CAM
      do n=1,s_dim_3d(3,i)
      do m=1,s_dim_3d(2,i)
      do k=1,s_dim_3d(1,i)
         temp_3d(m,n,k) = var%vars_3d(k,m,n,i)
         ! Meagher added this for quality control when assim cloud properties
         ! if((cflds(ifld)=='CLOUD   ').and.(temp_3d(m,n,k)<0)) temp_3d(m,n,k)=0
         ! if((cflds(ifld)=='CLOUD   ').and.(temp_3d(m,n,k)>1)) temp_3d(m,n,k)=1
      end do
      end do
      end do
   end if

   ifld = ifld + 1
   call nc_check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid),   &
                 'write_cam_init', 'inq_varid '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(ncfileid, ncfldid                          &
             ,temp_3d(1:f_dim_3d(1,i), 1:f_dim_3d(2,i), 1:f_dim_3d(3,i)) &
             ,start=(/1,1,1,1/) ,count=(/f_dim_3d(1,i), f_dim_3d(2,i), f_dim_3d(3,i), 1/) ), &
                 'write_cam_init', 'put_var '//trim(cflds(ifld)))
end do

! If model_time specified, write the time of the current state.
! All the caminput.nc files I have seen have two variables of
! length 'time' (the unlimited dimension): date, datesec
! The rest of the routine presumes there is but one time in the file -
! print warning message if this is not the case.

if (present( model_time)) then

   call nc_check(nf90_inq_dimid(ncfileid, 'time', dimid), &
          'write_cam_init', 'inq_dimid time '//trim(file_name))
   call nc_check(nf90_inquire_dimension(ncfileid, dimid, len=dimlen), &
          'write_cam_init', 'inquire_dimension time '//trim(file_name))

   if (dimlen /= 1) then
       write(string1,*)'UNUSUAL - ',trim(file_name),' has',dimlen,'times. Expected 1.'
       call error_handler(E_MSG, 'write_cam_init', string1, source, revision, revdate)
   end if

   allocate(datetmp(dimlen), datesec(dimlen))


   call nc_check(nf90_inq_varid(ncfileid, 'date', varid), &
          'write_cam_init', 'inq_varid date '//trim(file_name))
   call nc_check(nf90_get_var(ncfileid, varid, values=datetmp), &
          'write_cam_init', 'get_var date '//trim(file_name))

   call nc_check(nf90_inq_varid(ncfileid, 'datesec', varid), &
          'write_cam_init', 'inq_varid datesec '//trim(file_name))
   call nc_check(nf90_get_var(ncfileid, varid, values=datesec), &
          'write_cam_init', 'get_var datesec '//trim(file_name))

   call get_date(model_time, iyear, imonth, iday, ihour, imin, isec)


   ! The 'date' is YYYYMMDD ... datesec is 'current seconds of current day'
   datetmp(dimlen) = iyear*10000 + imonth*100 + iday
   datesec(dimlen) = ihour*3600 + imin*60 + isec

   call nc_check(nf90_put_var(ncfileid, varid, values=datetmp), &
          'write_cam_init', 'put_var date '//trim(file_name))

   call nc_check(nf90_put_var(ncfileid, varid, values=datesec), &
          'write_cam_init', 'put_var datesec '//trim(file_name))

   deallocate(datetmp, datesec)

end if

call nc_check(nf90_close(ncfileid), 'write_cam_init', 'close cam initial file')

deallocate (temp_3d, temp_2d)

end subroutine write_cam_init


!CS Not needed in CESM+DART framework?
   subroutine write_cam_times(model_time, adv_time)
!=======================================================================
! subroutine write_cam_times(model_time, adv_time)

! writes model time and advance time into a file called 'times'

type(time_type), intent(in) :: model_time, adv_time

integer :: tfile_unit, cam_date, cam_tod, nhtfrq
integer :: year, month, day, hour, minute, second
type(time_type) :: forecast_length


tfile_unit = open_file("times", "formatted", "write")

! end time is first, then beginning time
!  -namelist "&camexp START_YMD=$times[3] START_TOD=$times[4]
!                     STOP_YMD=$times[1] STOP_TOD=$times[2] NHTFRQ=$times[5] "

call get_date(adv_time, year, month, day, hour, minute, second)

cam_date = year*10000 + month*100 + day
cam_tod  = hour*3600  + minute*60 + second

write (tfile_unit,'(I8.8,1X,I8)') cam_date, cam_tod

call get_date(model_time, year, month, day, hour, minute, second)

cam_date = year*10000 + month*100 + day
cam_tod  = hour*3600  + minute*60 + second

write (tfile_unit,'(I8.8,1X,I8)') cam_date, cam_tod

! calculate number of hours in forecast, and pass to history tape
! write frequency

forecast_length = adv_time - model_time

call get_time(forecast_length, second, day)

hour = second/3600
minute = mod(second,3600)
if (minute.ne.0) &
   call error_handler(E_ERR, 'write_cam_times', &
      ' not integer number of hours; nhtfrq error', source, revision, revdate);

! convert to hours, and negative to signal units are hours

! nhtfrq = -1*((((year-1)*365 + (month-1))*30 + (day-1))*24 + hour)
nhtfrq = -1*(day*24 + hour)
write (tfile_unit,'(I8)') nhtfrq

close(tfile_unit)


end subroutine write_cam_times


!> This is just renamed to it with assim_tools_mod because I renamed wrf, get_state_meta_data_distrib
!> It does not do anything different yet.
   subroutine get_state_meta_data_distrib(state_ens_handle, index_in, location, var_kind)
!=======================================================================
! subroutine get_state_meta_data(index_in, location, var_kind, set_loc)
!
! Given an integer index into the state vector structure, returns the
! associated location.
! The location may have components that are MISSING_R8 values, since some fields
! don't have locations in all three dimensions, i.e. PS has no vertical level,
! and other fiendish fields to be devised by parameterization studies may not
! have a longitude, or latitude.  The which_vert should take care of the vertical
! coordinate (it will be ignored), but the others will require more interesting  fixes.
! See order_state_fields for the QTY_s (and corresponding model_mod TYPE_s).
!
! This is not a function because the more general form of the call has a second
! intent(out) optional argument var_kind.  Maybe a functional form should be added?

type(ensemble_type), intent(in)  :: state_ens_handle
integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_kind

integer  :: which_vert
integer  :: i, indx, index_1, index_2, index_3, nfld
integer  :: box, slice
logical  :: lfound

! CS replaced 'goto 10' mechanism with
! if (.not.lfound) then
!    ...
!    lfound = .true.
!    exit
!    ...
! end if
! I could get rid of 'exit' by putting if-test inside loops,
! but then loops would be entered when not necessary.

! Return vertical coordinate as one of the following.
! integer, parameter :: VERTISUNDEF    = -2 ! has no vertical location (undefined)
! integer, parameter :: VERTISSURFACE  = -1 ! surface value
! integer, parameter :: VERTISLEVEL    =  1 ! by level
! integer, parameter :: VERTISPRESSURE =  2 ! by pressure
! integer, parameter :: VERTISHEIGHT   =  3 ! by height


real(r8) :: lon_val, lat_val, lev_val

character(len=8)   :: dim_name

lfound    = .false.

! In order to find what variable this is, and its location, I must subtract the individual
! components of the state vector, since they may have varying sizes.
! Save the original index.
! index_in will be < 0 if it's an identity obs (called from convert_vert)
indx = abs(index_in)
which_vert = MISSING_I
index_1 = 0; index_2 = 0; index_3 = 0; nfld = 0
lon_val = MISSING_R8; lat_val = MISSING_R8; lev_val = MISSING_R8

! Cycle through 0d state variables
do i=1,state_num_0d
   nfld = nfld + 1
   if (indx == 1 ) then
      which_vert = VERTISUNDEF
      !goto 10
      lfound = .true.
      exit
   else
      indx = indx - 1
   end if
end do

! Cycle through 1d state variables
! Note that indices of fields can have varying dimensions.
if (.not.lfound) then
do i=1,state_num_1d
   nfld = nfld + 1
   if (indx > s_dim_1d(i) ) then
      indx = indx - s_dim_1d(i)
   else
      ! We've found the desired field; now find lat, lon or lev of indx
      dim_name = dim_names(s_dimid_1d(i))
      !  Don't return the value of lev(index_1), which = 1000(A+B)
      !  dim_name = dim_names(s_dimid_1d(1,i))
      !  call coord_val(dim_name, index_1, lon_val, lat_val, lev_val)
      if (dim_name == 'lev     ') then
         lev_val = real(indx)
      else
         call coord_val(dim_name, indx, lon_val, lat_val, lev_val)
      end if

      which_vert = which_vert_1d(i)
      !goto 10
      lfound = .true.
      exit
   end if
end do
end if

! Cycle through 2d state variables.
! Note that indices of fields can have varying dimensions.
if (.not.lfound) then
do i=1,state_num_2d
   nfld = nfld + 1
   slice = s_dim_2d(1,i) * s_dim_2d(2,i)
   if (indx > slice ) then
      indx = indx - slice
   else
      ! We've found the desired field.
      ! Now find lat and/or lon and/or lev of indx if called by assim_tools_mod:filter_assim

      ! # second dimension rows to subtract off; temporary value for index_2
      index_2 = (indx -1) / s_dim_2d(1,i)
      index_1 = indx - (index_2 * s_dim_2d(1,i))
      dim_name = dim_names(s_dimid_2d(1,i))
      ! CS Check that this works for dimensions (nlev,ncol).
      if (do_out .and. indx == 1) then
         write(*,'(A,3I7,A)') 'get_state_meta_data: index_in, index_1, index_2, dim_name', &
                                                    index_in, index_1, index_2, dim_name
         write(*,'(A,I7,A,2I7)') '   s_dim_2d(1:2,',i,') = ',s_dim_2d(1,i), s_dim_2d(2,i)
      end if

      ! Find the coordinate value (i.e. 270.5) of the first dimension index (i.e. 54)
      if (dim_name == 'lev     ') then
         lev_val = real(index_1)
      else
         call coord_val(dim_name, index_1, lon_val, lat_val, lev_val)
      end if

      ! index_2 of the variable in question is 1 more than the # subtracted off to get index_1
      index_2 = index_2 + 1
      dim_name = dim_names(s_dimid_2d(2,i))
      !  Don't return the value of lev(index_1), which = 1000(A+B)
      !  dim_name = dim_names(s_dimid_2d(1,i))
      !  call coord_val(dim_name, index_1, lon_val, lat_val, lev_val)
      if (dim_name == 'lev     ') then
         lev_val = real(index_2)
      else
         call coord_val(dim_name, index_2, lon_val, lat_val, lev_val)
      end if

      which_vert = which_vert_2d(i)

      !goto 10
      lfound = .true.
      exit
   end if
end do
end if

! Cycle through 3d state variables
! Note that indices of fields can have varying dimensions.
if (.not.lfound) then
do i=1,state_num_3d
   nfld = nfld + 1
   box = s_dim_3d(1,i) * s_dim_3d(2,i) * s_dim_3d(3,i)
   if (indx > box ) then
      indx = indx - box
   else
      ! We've found the desired field.
      ! Now find lat and/or lon and/or lev of indx if called by assim_tools_mod:filter_assim

         ! # of (first x second dimension) slices to subtract off in order to find the current slice
         ! debug; try printing out all this index info when there's just 1 obs to be processed.
         slice = s_dim_3d(1,i) * s_dim_3d(2,i)
         index_3 = (indx -1) / slice         ! temporary value used to find index_2 and index_1
         index_2 = (indx -1 - (index_3 * slice)) / s_dim_3d(1,i)    ! same for index_2 to find index_1
         index_1 = (indx - (index_3 * slice) - (index_2 * s_dim_3d(1,i)))

   !     Don't return the value of lev(index_1), which = 1000(A+B)
   !      dim_name = dim_names(s_dimid_3d(1,i))
   !      call coord_val(dim_name, index_1, lon_val, lat_val, lev_val)
   !     Return index_1 as the vertical location
         lev_val = real(index_1)

         ! index_2 of the variable in question is one more than the number subtracted off to get index_1
         index_2 = index_2 + 1
         dim_name = dim_names(s_dimid_3d(2,i))
         call coord_val(dim_name, index_2, lon_val, lat_val, lev_val)

         ! index_3 of the variable in question is one more than the number subtracted off to get index_1
         index_3 = index_3 + 1
         dim_name = dim_names(s_dimid_3d(3,i))
         call coord_val(dim_name, index_3, lon_val, lat_val, lev_val)

         which_vert = which_vert_3d(i)

      !goto 10
      lfound = .true.
      exit
   end if
end do
end if

! 10 continue

! This will malfunction for fields that are filled with MISSING_R8 for lat_val or lon_val.
if (lon_val == MISSING_R8 .or. lat_val == MISSING_R8 ) then
   write(string1, *) 'Field ',cflds(nfld),' has no lon or lat dimension.  ', &
         'What should be specified for it in the call to location?'
   call error_handler(E_ERR, 'get_state_meta_data', string1, source, revision, revdate)
else
   if (lat_val <= -90.0_r8) lat_val = -89.9999999_r8
   if (lat_val >=  90.0_r8) lat_val =  89.9999999_r8
   ! CS; Set arguments must be in degrees.
   location = set_location(lon_val, lat_val, lev_val, which_vert)
end if

! If the type is wanted, return it
if (present(var_kind)) then
   if (index_in < 0) then
      ! used by convert_vert which wants the CAM field index, not the DART QTY_
      var_kind = nfld
   else if (index_in > 0) then
      ! used by call from assim_tools_mod:filter_assim, which wants the DART QTY_
      var_kind = cam_to_dart_kinds(nfld)
   end if
end if

end subroutine get_state_meta_data_distrib




   subroutine ens_mean_for_model(filter_ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: filter_ens_mean(:)

call error_handler(E_ERR, 'ens_mean_for_model', 'not allowed in distributed version')

end subroutine ens_mean_for_model


   function get_model_size()
!=======================================================================
! function get_model_size()
!

integer :: get_model_size

get_model_size = model_size

end function get_model_size




   function get_model_time_step()
!=======================================================================
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be replaced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

! Time_step_atmos is global static storage
get_model_time_step =  Time_step_atmos

end function get_model_time_step



   function nc_write_model_atts( ncFileID ) result (ierr)
!=======================================================================
! function nc_write_model_atts( ncFileID ) result (ierr)
!
! Writes the model-specific attributes to a netCDF file.
! TJH Fri Aug 29 MDT 2003
!

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!-----------------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: MemberDimID, StateVarDimID, TimeDimID,ScalarDimID
integer :: xVarID,StateVarID, StateVarVarID
integer :: P_id(num_dims)
integer :: i, ifld, dim_id, g_id
integer :: grid_id(grid_num_1d)
character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

ierr = 0     ! assume normal termination

!-------------------------------------------------------------------------------
! Make sure ncFileID refers to an open netCDF file,
! and then put into define mode.
! nf90_Inquire  returns all but the ncFileID; these were defined in the calling routine.
!    More dimensions, variables and attributes will be added in this routine.
!-------------------------------------------------------------------------------

write(string1,*) 'ncFileID', ncFileID
call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID), &
              'nc_write_model_atts', 'Inquire '//trim(string1))
call nc_check(nf90_Redef(ncFileID), 'nc_write_model_atts', 'Redef '//trim(string1))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID), &
              'nc_write_model_atts', 'inq_dimid copy')
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID), &
              'nc_write_model_atts', 'inq_dimid time')

if ( TimeDimID /= unlimitedDimId ) then
  write(string1,*)'Time dimension ID ',TimeDimID,'must match Unlimited Dimension ID ',unlimitedDimId
  call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
end if

!-------------------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                        len=model_size, dimid = StateVarDimID),  &
              'nc_write_model_atts', 'def_dim StateVariable')

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
! I think this format is screwing up the next statement
! write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
write(str1,'("YYYY MM DD HH MM SS = ",i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1),        &
              'nc_write_model_atts', 'put_att creation_date'//trim(str1))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision),   &
              'nc_write_model_atts', 'put_att model_revision'//trim(revision))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate),     &
              'nc_write_model_atts', 'put_att model_revdate'//trim(revdate))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","CAM"),               &
              'nc_write_model_atts','put_att model CAM')

! how about namelist input? might be nice to save ...

!-------------------------------------------------------------------------------
! Define the new dimensions IDs
!-------------------------------------------------------------------------------

! All the dim_ids on caminput.nc were read in, and will be written out here,
! and some will be used to define the coordinate variables below.
! They have different dimids for this file than they had for caminput.nc
! P_id serves as a map between the 2 sets.
if (print_details .and. do_out) then
   PRINT*,'num_dims = ',num_dims
   write(*,*) ' dimens,       name,  size, cam dim_id, P[oste]rior id'
end if
do i = 1,num_dims
   if (trim(dim_names(i)) /= 'time')  then
      call nc_check(nf90_def_dim (ncid=ncFileID, name=trim(dim_names(i)), len=dim_sizes(i),  &
                    dimid=P_id(i)), 'nc_write_model_atts','def_dim '//trim(dim_names(i)))
   else
     P_id(i) = 0
   end if
   if (print_details .and. do_out) write(*,'(I5,1X,A13,1X,3(I7,2X))') &
      i,trim(dim_names(i)),dim_sizes(i), dim_ids(i), P_id(i)
end do
call nc_check(nf90_def_dim(ncid=ncFileID, name="scalar",   len = 1,   dimid = ScalarDimID) &
             ,'nc_write_model_atts', 'def_dim scalar')

!-------------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and their attributes
!-------------------------------------------------------------------------------

! grid longitudes, latitudes, levels, and other coordinates.
! grid_id() is filled here; it's the dimid of the desired coordinate *on this P_Diag.nc file*.
! It's used to write coordinates.  dim_ids keeps track of the dimensions of coordinates
! and fields on caminput.nc files.  There's some overlap of names, unfortunately.
! The argument after the 'xxx    ' label is a structure with all the relevant info in it.
! The structures are defined in "Grid fields" and filled by calls to init_grid_1d_instance
! in read_cam_coord.
! CS; ncol doesn't belong here because it's just a dimension, not a coordinate variable.

grid_id = MISSING_I

if (lon%label /= '        ')  then
   dim_id = P_id(lon%dim_id)
   g_id   = find_name('lon     ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'lon     ',lon , dim_id, grid_id(g_id))
end if
if (lat%label /= '        ')  then
   dim_id = P_id(lat%dim_id)
   g_id   = find_name('lat     ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'lat     ',lat , dim_id, grid_id(g_id))
end if
if (lev%label /= '        ')  then
   dim_id = P_id(lev%dim_id)
   g_id   = find_name('lev     ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'lev     ',lev , dim_id, grid_id(g_id))
! Gaussian weights -- because they're there.
end if
if (gw%label /= '        ')  then
   dim_id = P_id(gw%dim_id)
   g_id   = find_name('gw      ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'gw      ',gw  , dim_id, grid_id(g_id))
! Hybrid grid level coefficients, parameters
end if
if (hyam%label /= '        ')  then
   dim_id = P_id(hyam%dim_id)
   g_id   = find_name('hyam    ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'hyam    ',hyam, dim_id, grid_id(g_id))
end if
if (hybm%label /= '        ')  then
   dim_id = P_id(hybm%dim_id)
   g_id   = find_name('hybm    ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'hybm    ',hybm, dim_id, grid_id(g_id))
end if
if (hyai%label /= '        ')  then
   dim_id = P_id(hyai%dim_id)
   g_id   = find_name('hyai    ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'hyai    ',hyai, dim_id, grid_id(g_id))
end if
if (hybi%label /= '        ')  then
   dim_id = P_id(hybi%dim_id)
   g_id   = find_name('hybi    ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'hybi    ',hybi, dim_id, grid_id(g_id))
end if
if (slon%label /= '        ')  then
   dim_id = P_id(slon%dim_id)
   g_id   = find_name('slon    ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'slon    ',slon, dim_id, grid_id(g_id))
end if
if (slat%label /= '        ')  then
   dim_id = P_id(slat%dim_id)
   g_id   = find_name('slat    ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'slat    ',slat, dim_id, grid_id(g_id))
end if
if (ilev%label /= '        ')  then
   dim_id = P_id(ilev%dim_id)
   g_id   = find_name('ilev    ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'ilev    ',ilev, dim_id, grid_id(g_id))
end if
if (P0%label /= '        ')  then
   dim_id = P_id(P0%dim_id)
   ! dim_id here is 0; will that work?  It's what's read in from caminput.nc
   ! If not, then I'll need to (re)define grid_0d_kind, etc.
   g_id   = find_name('P0      ',grid_names_1d)
   call write_cam_coord_def(ncFileID,'P0      ',P0  , dim_id, grid_id(g_id))
end if

if (print_details .and. do_out) then
   write(*,*) '1d field#, grid_id, grid_names_1d'
   do i=1,grid_num_1d
      write(*,*) 'grid_ = ', i, grid_id(i), trim(grid_names_1d(i))
   end do
end if

if ( output_state_vector ) then

! CAM3; need to adapt this to state_long_names, state_units, etc?


   !----------------------------------------------------------------------------
   ! Create attributes for the state vector
   !----------------------------------------------------------------------------

   ! Define the state vector coordinate variable
   call nc_check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int,           &
              dimids=StateVarDimID, varid=StateVarVarID),                                   &
                 'nc_write_model_atts','def_var  state vector')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"),   &
                 'nc_write_model_atts','put_att long_name state vector ')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical"),           &
                 'nc_write_model_atts','put_att units state vector ' )
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)), &
                 'nc_write_model_atts','put_att valid range state vector ')
   ! Define the actual state vector
   call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real,                 &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), varid=StateVarID), &
                 'nc_write_model_atts','def_var state vector')
   call nc_check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"),   &
                 'nc_write_model_atts','put_att long_name model state or fcopy ')
   call nc_check(nf90_put_att(ncFileID, StateVarId, "vector_to_prog_var","CAM"),            &
                 'nc_write_model_atts','put_att vector_to_prog_var CAM ')

   ! Leave define mode so we can fill
   call nc_check(nf90_enddef(ncfileID), 'nc_write_model_atts','enddef ')

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ),         &
                 'nc_write_model_atts','put_var StateVar ')

else

   !----------------------------------------------------------------------------
   ! Create the (empty) Variables and the Attributes
   !----------------------------------------------------------------------------

   ! 0-d fields
   ifld = 0
   do i = 1,state_num_0d
      ifld = ifld + 1
      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ MemberDimID, unlimitedDimID /),                             &
                 varid  = xVarID),                                                       &
                 'nc_write_model_atts','def_var 0d '//trim(cflds(ifld)))
      call nc_check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)), &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)),          &
                 'nc_write_model_atts','put_att units ')
!       call nc_check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)),                         &
!                 'nc_write_model_atts','def_var ')
   end do

   ! 1-d fields
   do i = 1,state_num_1d
      ifld = ifld + 1
      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ P_id(s_dimid_1d(1)), MemberDimID, unlimitedDimID /),        &
                 varid  = xVarID),                                                       &
                 'nc_write_model_atts','def_var 1d '//trim(cflds(ifld)))
      call nc_check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)), &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)),          &
                 'nc_write_model_atts','put_att units ')
!       call nc_check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)),                         &
!                 'nc_write_model_atts','def_var ')
   end do

   ! 2-d fields
   do i = 1,state_num_2d
      ifld = ifld + 1
      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ P_id(s_dimid_2d(1,i)), P_id(s_dimid_2d(2,i)),               &
                             MemberDimID, unlimitedDimID /),                             &
                 varid  = xVarID),                                                       &
                 'nc_write_model_atts','def_var 2d '//trim(cflds(ifld)))
      call nc_check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)), &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)),          &
                 'nc_write_model_atts','put_att units ')
!       call nc_check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)),                         &
!                 'nc_write_model_atts','def_var ')
   end do

   ! 3-d fields
   do i = 1,state_num_3d
      ifld = ifld + 1
      call nc_check(nf90_def_var                                                              &
           (ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real,                           &
            dimids = (/ P_id(s_dimid_3d(1,i)), P_id(s_dimid_3d(2,i)), P_id(s_dimid_3d(3,i)),  &
                        MemberDimID, unlimitedDimID /),                                       &
            varid  = xVarID),                                                                 &
                 'nc_write_model_atts','def_var 3d'//trim(cflds(ifld)))
      call nc_check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)),      &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)),               &
                 'nc_write_model_atts','put_att units ')
!       call nc_check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)),                         &
!                 'nc_write_model_atts','def_var ')
   end do

   ! Leave define mode so we can fill variables
   call nc_check(nf90_enddef(ncfileID), &
                 'nc_write_model_atts','enddef ')

!-------------------------------------------------------------------------------
! Fill the coordinate variables
! Each 'vals' vector has been dimensioned to the right size for its coordinate.
! The default values of 'start' and 'count'  write out the whole thing.
! CS; ncol doesn't belong here because it's just a dimension, not a coordinate variable.
!-------------------------------------------------------------------------------

if (lon%label  /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('lon     ',grid_names_1d)),  lon%vals) &
                 ,'nc_write_model_atts', 'put_var lon')
if (lat%label  /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('lat     ',grid_names_1d)),  lat%vals) &
                 ,'nc_write_model_atts', 'put_var lat')
if (lev%label  /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('lev     ',grid_names_1d)),  lev%vals) &
                 ,'nc_write_model_atts', 'put_var lev')
if (gw%label   /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('gw      ',grid_names_1d)),   gw%vals) &
                 ,'nc_write_model_atts', 'put_var gw')
if (hyam%label /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('hyam    ',grid_names_1d)), hyam%vals) &
                 ,'nc_write_model_atts', 'put_var hyam')
if (hybm%label /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('hybm    ',grid_names_1d)), hybm%vals) &
                 ,'nc_write_model_atts', 'put_var hybm')
if (hyai%label /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('hyai    ',grid_names_1d)), hyai%vals) &
                 ,'nc_write_model_atts', 'put_var hyai')
if (hybi%label /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('hybi    ',grid_names_1d)), hybi%vals) &
                 ,'nc_write_model_atts', 'put_var hybi')
if (slon%label /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('slon    ',grid_names_1d)), slon%vals) &
                 ,'nc_write_model_atts', 'put_var slon')
if (slat%label /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('slat    ',grid_names_1d)), slat%vals) &
                 ,'nc_write_model_atts', 'put_var slat')
if (ilev%label /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('ilev    ',grid_names_1d)), ilev%vals) &
                 ,'nc_write_model_atts', 'put_var ilev')
if (P0%label   /= '        ')                                                                     &
    call nc_check(nf90_put_var(ncFileID, grid_id(find_name('P0      ',grid_names_1d)),   P0%vals) &
                 ,'nc_write_model_atts', 'put_var P0')

end if

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID),'nc_write_model_atts', 'sync ')

write (*,*)'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...'

end function nc_write_model_atts



   function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!=======================================================================
! function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!
! Writes the model-specific variables to a netCDF file
! TJH 25 June 2003
!

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!-----------------------------------------------------------------------------------------
type(model_type) :: Var

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, ncfldid
integer :: ifld, i

character (len=8) :: cfield

ierr = 0     ! assume normal termination

!-------------------------------------------------------------------------------
! CAM storage bounds are 'tight' -- no indices needed
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
! then get all the Variable ID's we need.
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID), &
              'nc_write_model_vars','Inquire ')

if ( output_state_vector ) then

   call nc_check(nf90_inq_varid(ncFileID, "state", StateVarID),'nc_write_model_vars ','inq_varid state' )
   call nc_check(nf90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)),'nc_write_model_vars ','put_var state')

else

   !----------------------------------------------------------------------------
   ! Fill the variables
   !----------------------------------------------------------------------------

   call init_model_instance(Var)     ! Explicity released at end of routine.

   call vector_to_prog_var(statevec,  Var)

   ifld = 0
   ZeroDVars : do i = 1, state_num_0d
      ifld = ifld + 1
      cfield = trim(cflds(ifld))
      call nc_check(nf90_inq_varid(ncFileID, cfield, ncfldid),       &
                    'nc_write_model_vars ','inq_varid 0d '//cfield)
      call nc_check(nf90_put_var(ncFileID, ncfldid, Var%vars_0d(i),             &
                                 start=(/ copyindex, timeindex /) ),            &
                    'nc_write_model_vars ','put_var 0d '//cfield)
!, &
!                 count=(/1, 1/) ))
   end do ZeroDVars

   ! 'start' and 'count' are needed here because Var%vars_Nd are dimensioned by the largest
   ! values for the dimensions of the rank N fields, but some of the fields are smaller than that.
   OneDVars : do i = 1, state_num_1d
      ifld = ifld + 1
      cfield = trim(cflds(ifld))
      call nc_check(nf90_inq_varid(ncFileID, cfield, ncfldid),                                    &
                    'nc_write_model_vars ','inq_varid 1d '//cfield)
      call nc_check(nf90_put_var(ncFileID, ncfldid,                                               &
                    Var%vars_1d(1:s_dim_1d(  i), i),                                              &
                    start=   (/ 1              ,copyindex, timeindex /),                          &
                    count=   (/   s_dim_1d(  i),1        , 1/) ),                                 &
                    'nc_write_model_vars ','put_var 1d '//cfield)
   end do OneDVars

   ! Write out 2D variables as 2 of (lev,lon,lat), in that order, regardless of caminput.nc
   ! coordinate order
   ! The sizes can be taken from s_dim_2d, even though the s_dimid_2d don't pertain to this
   ! P_Diag.nc file, because the dimids were mapped correctly in nc_write_model_atts.
   TwoDVars : do i = 1, state_num_2d
      ifld = ifld + 1
      cfield = trim(cflds(ifld))
      call nc_check(nf90_inq_varid(ncFileID, cfield, ncfldid),                                    &
                    'nc_write_model_vars ','inq_varid 2d '//cfield)
      call nc_check(nf90_put_var(ncFileID, ncfldid,                                               &
                    Var%vars_2d(1:s_dim_2d(1,i),1:s_dim_2d(2,i), i),                              &
                    start=   (/ 1              ,1              , copyindex, timeindex /),         &
                    count=   (/   s_dim_2d(1,i),  s_dim_2d(2,i), 1        , 1/) ),                &
                    'nc_write_model_vars ','put_var 2d '//cfield)
   end do TwoDVars

   ! Write out 3D variables as (lev,lon,lat) regardless of caminput.nc coordinate order
   ThreeDVars : do i = 1,state_num_3d
      ifld = ifld + 1
      cfield = trim(cflds(ifld))
      call nc_check(nf90_inq_varid(ncFileID, cfield, ncfldid),                                    &
                    'nc_write_model_vars ','inq_varid 3d '//cfield)
      call nc_check(nf90_put_var(ncFileID, ncfldid,                                               &
                 Var%vars_3d(1:s_dim_3d(1,i),1:s_dim_3d(2,i),1:s_dim_3d(3,i),i)                   &
                 ,start=   (/1              ,1              ,1              ,copyindex,timeindex/)&
                 ,count=   (/  s_dim_3d(1,i),  s_dim_3d(2,i),  s_dim_3d(3,i),1        ,1/) ),     &
                    'nc_write_model_vars ','put_var 3d '//cfield)
   end do ThreeDVars

end if

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

!write (*,*)'Finished filling variables ...'
call nc_check(nf90_sync(ncFileID),'nc_write_model_vars ','sync ')
!write (*,*)'netCDF file is synched ...'

! temporary output
!print*,'num_calced, num_searched = ',num_calced, num_searched

call end_model_instance(Var)   ! should avoid any memory leaking

end function nc_write_model_vars


! End of Module I/O

!#######################################################################

! model_interpolate section

!> Distributed version of model_interpolate
subroutine model_interpolate_distrib(state_ens_handle, location, obs_kind, istatus, interp_val)
!=======================================================================
! CS; this subroutine is now a short routine that calls
!     either a rectangular grid version (the old model_interpolate for eul/FV)
!     or new non-rectangular (cubed-sphere) code.
! C2D: This does get KINDs from filter, not specific obs TYPEs.

type(ensemble_type), intent(in) :: state_ens_handle
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_kind
integer,            intent(out) :: istatus(:)
real(r8),           intent(out) :: interp_val(:)

if (l_rectang) then
   call interp_lonlat_distrib(state_ens_handle, location, obs_kind, interp_val, istatus)
else
   call interp_cubed_sphere_distrib(state_ens_handle, location, obs_kind, interp_val, istatus)
end if

end subroutine model_interpolate_distrib


   subroutine interp_cubed_sphere_distrib(state_ens_handle, obs_loc, obs_kind, interp_val, istatus)
!=======================================================================
!  subroutine interp_cubed_sphere(st_vec, obs_loc, obs_kind, interp_val, istatus)
!
! Find the cell that encloses an ob at 'obs_loc'
! and interpolate the values of obs_kind from the cell's corners to that location.

type(ensemble_type), intent(in) :: state_ens_handle
type(location_type), intent(in) :: obs_loc
integer,             intent(in) :: obs_kind
integer,            intent(out) :: istatus(:)
real(r8),           intent(out) :: interp_val(:)

integer               :: i, closest
integer, allocatable  :: track_vstatus(:), vstatus(:)

integer               :: s_type, s_type_1d, s_type_2d
integer, dimension(4) :: quad_corners   ! node numbers of the corners of the enclosing cell
character (len=8)     :: col_name, lev_name
real(r8)              :: lon_lat_lev(3), l, m
real(r8), allocatable :: vals(:, :)
integer :: ens_size, e

ens_size = state_ens_handle%num_copies -5
allocate(track_vstatus(ens_size), vstatus(ens_size))
allocate(vals(ens_size, 4))

! Start with no errors in
istatus(:) = 0
vstatus(:) = 0
track_vstatus(:) = 0

! Get the observation (horizontal) position, in degrees
! CS; truly 'degrees'?  Yes, in location_mod location%{lon,lat} are stored as radians,
!     but get_location converts them.
lon_lat_lev = get_location(obs_loc)

! check whether model_mod can interpolate the requested variable
! Pressure (3d) can't be specified as a state vector field (so s_type will = MISSING_I),
! but can be calculated for CAM, so obs_kind = QTY_PRESSURE is acceptable.
s_type = dart_to_cam_kinds(obs_kind)
if (s_type == MISSING_I .and. &
   (obs_kind .ne. QTY_PRESSURE) .and.  (obs_kind .ne. QTY_SURFACE_ELEVATION)) then
   istatus = 3
   interp_val = MISSING_R8
   write(*,*) 'Wrong type of obs = ', obs_kind
   return
end if

! Get horizontal grid specs

! Set [ncol,lev] names to defaults, which may be overwritten for variables in the state vector,
! but not for other acceptable variables (3D pressure, surface elevation, ...?)
col_name = 'ncol    '

! ? How to separate the 3D 'other' variables from 2D 'other' variables?
!   Can't do it automatically/generically because they're not part of state vector
!   and that info isn't coming from DART.
if (obs_kind .eq. QTY_SURFACE_ELEVATION) then
   lev_name = 'none    '
else if (obs_kind .eq. QTY_PRESSURE) then
   lev_name = 'lev     '
end if

! See interp_lonlat for original code.
! There can't be any 0d ob fields.
!     CS; what about earth rotation obs?  Are they converted into 2D fields?
! And PS is a 1d field in SE-CAM.
! Positions within the rank 2 and 3 fields.
s_type_1d = s_type - state_num_0d
s_type_2d = s_type_1d - state_num_1d

if (s_type == MISSING_I .and. &
   (obs_kind .eq. QTY_PRESSURE) .or. (obs_kind .eq. QTY_SURFACE_ELEVATION)) then
   ! use defaults col_name set above
else if (s_type <= state_num_0d ) then
   ! error; can't deal with observed variables that are 0D in model_mod.
   istatus(:) = 3
   interp_val(:) = MISSING_R8
   write(*,*) 'Cannot handle 0D state vars, s_type = ', s_type
   return
else if (s_type_1d > 0 .and. s_type_1d <= state_num_1d) then
   col_name = dim_names(s_dimid_1d(s_type_1d))
else if (s_type_2d > 0 .and. s_type_2d <= state_num_2d) then
   col_name = dim_names(s_dimid_2d(1,s_type_2d))
   lev_name = dim_names(s_dimid_2d(2,s_type_2d))
else
   istatus(:) = 3
   interp_val(:) = MISSING_R8
   write(*,*) 'Unexpected state type value, s_type = ', s_type
   return
end if


! Need the node names of the corners of the cell which contains the observation.
! Also return the location (l,m) of the ob in the unit square space,
! since it needs to be calculated to find the cell, and is needed for the interpolation.
! It looks silly to pass cs_locs and cs_kinds (globally available) here,
! but other calls to coord_ind_cs need different arguments there,
! so these arguments need to be in the list.
! convert_vert:
if (print_details .and. do_out) then
   write(*,'(2A,I3)') 'interp_cubed_sphere; col_name, obs_kind, s_type ',col_name, obs_kind, s_type
   write(*,'(A,1p3E20.12)') 'lon,lat,lev of obs = ',get_location(obs_loc)
   write(*,*) 'Calling coord_ind_cs with closest_only = false'
end if

call coord_ind_cs(obs_loc, obs_kind, cs_locs, cs_kinds, quad_corners, l, m, closest, .false.)
! obs_which not needed?                  .false., nint(query_location(obs_loc)))
if (print_details .and. do_out) then
   write(string1,'(A,4I8)') 'quad_corners = ',(quad_corners(i), i=1,4)
   call error_handler(E_MSG, 'interp_cubed_sphere', string1)
   write(*,'(A,I8,A,1p2E20.12) ')'{lon,lat}(',quad_corners(4),' = ',lon%vals(quad_corners(4)), lat%vals(quad_corners(4))
!    write(string1,'(A,4I8)') 'obs_kind  before get_val_pressure = ',&
!         obs_kind
!    call error_handler(E_MSG, 'interp_cubed_sphere', string1)
end if

! The interpolation.
! First interpolate the field in the vertical at each of the 4 corners
! to the height of the ob.
! The subroutines and arrays appear to want indices for the lon and lat dimensions,
! while the cubed sphere has only the ncol horizontal dimension.
! This is handled by passing the ncol index as 'lon_index', and 1 as 'lat_index'.
! Then get_val (the bottom of the calling trees) uses these correctly for the cubed sphere.


if (obs_kind == QTY_SURFACE_ELEVATION) then
   ! Acceptable KIND that's not in the state vector
   ! convert from geopotential height to real height in meters
   call error_handler(E_ERR, 'model_interpolate ', 'QTY_SURFACE_ELEVATION not done')

! Move this to the end of the block?  It's no good here; short circuits GPS
! which asks for pressures on heights
! else if (obs_kind == QTY_PRESSURE) then
!    ! Calculate pressures from surface pressures and A and B coeffs.

else if (vert_is_level(obs_loc)) then
   ! Case 1: model level specified in vertical
   ! Pobs

  call error_handler(E_ERR, 'model_interpolate ', 'vert_is_level')
   ! Pobs end

else if (vert_is_pressure(obs_loc)) then
   ! which_vert is pressure for this obs

   call get_val_pressure_distrib(vals(:, 1), state_ens_handle, quad_corners(1), 1, lon_lat_lev(3), obs_kind, vstatus)
   track_vstatus = vstatus

   call get_val_pressure_distrib(vals(:, 2), state_ens_handle, quad_corners(2), 1, lon_lat_lev(3), obs_kind, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   call get_val_pressure_distrib(vals(:, 3), state_ens_handle, quad_corners(3), 1, lon_lat_lev(3), obs_kind, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   call get_val_pressure_distrib(vals(:, 4), state_ens_handle, quad_corners(4), 1, lon_lat_lev(3), obs_kind, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo
   vstatus = track_vstatus

else if (vert_is_height(obs_loc)) then
   ! which_vert is height for this obs
   ! HK I don't understand why you are passing lon_lat_lev(3) if you are passing obs_loc
   call get_val_height_distrib(vals(:, 1), state_ens_handle, quad_corners(1), 1, lon_lat_lev(3), obs_loc, obs_kind, vstatus)
   track_vstatus = vstatus

   call get_val_height_distrib(vals(:, 2), state_ens_handle, quad_corners(2), 1, lon_lat_lev(3), obs_loc, obs_kind, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   call get_val_height_distrib(vals(:, 3), state_ens_handle, quad_corners(3), 1, lon_lat_lev(3), obs_loc, obs_kind, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   call get_val_height_distrib(vals(:, 4), state_ens_handle, quad_corners(4), 1, lon_lat_lev(3), obs_loc, obs_kind, vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo
   vstatus = track_vstatus

else if (vert_is_surface(obs_loc)) then
   ! location_mod:interactive_location asks for surface obs to have vertical coord = ps(hPa)
   ! The 'lev' argument is set to 1 because there is no level for these types, and 'lev' will be
   ! ignored.
   ! HK get_val_distrib only fails if the index is out of the grid, that is why there is only one return code
   ! not ens_size of them.
   call get_val_distrib(vals(:, 1), state_ens_handle, ens_size, quad_corners(1), 1, 1, obs_kind, vstatus(1))
   if (vstatus(1) /= 1) call get_val_distrib(vals(:, 2), state_ens_handle, ens_size, quad_corners(2), 1, 1, obs_kind, vstatus(1))
   if (vstatus(1) /= 1) call get_val_distrib(vals(:, 3), state_ens_handle, ens_size, quad_corners(3), 1, 1, obs_kind, vstatus(1))
   if (vstatus(1) /= 1) call get_val_distrib(vals(:, 4), state_ens_handle, ens_size, quad_corners(4), 1, 1, obs_kind, vstatus(1))
   vstatus = track_vstatus(1)

! Need option for vert_is_undefined
else
   write(*,*) '   No vert option chosen!'
end if

! HK loop around ensemble
do e = 1, ens_size

   ! lat is already converted to degrees by get_location
   if (abs(lon_lat_lev(2)) > max_obs_lat_degree .and. vstatus(e) /= 1) then
     istatus(e) = 4
   else
      istatus(e) = vstatus(e)
   end if

enddo

! Then interpolate horizontally to the (lon,lat) of the ob.
! The following uses Jeff's recommended 'generalized quadrilateral interpolation', as in
! http://www.particleincell.com/2012/quad-interpolation/.
! Most of the work is done in create_cs_grid_arrays() and coord_ind_cs().

! Interpolate from the cell's corners to the ob location on the unit square.
! This is done by weighting the field at each corner by the rectangular area
! ((l,m) space) diagonally across the ob location from the corner.
! AKA 'linear area weighting'.

! The vals indices are consistent with how mapping of corners was done,
! and how quad_corners was assigned.

! HK loop around ensemble
do e = 1, ens_size

   if (istatus(e) /= 1) then
      interp_val(e) = vals(e, 2) *        l *       m  &
                 + vals(e, 1) * (1._r8-l)*       m  &
                 + vals(e, 4) * (1._r8-l)*(1._r8-m) &
                 + vals(e, 3) *        l *(1._r8-m)
   !Exp8 debug; more output   if (print_details .and. do_out) then
      if (print_details ) then
         ! write(string1,'(A,2F10.6,1pE20.12)') 'quad_interp; l,m, interpolated val = ',l,m,interp_val
         write(string1,'(A,2F10.6,1pE20.12)') ' l,m, interpolated val = ', &
               l,m,interp_val
         call error_handler(E_MSG, 'interp_cubed_sphere', string1)
      end if
   else
      interp_val(e) = MISSING_R8 !HK do you want this in here?
      if (print_details) then
      write(string1,'(A,2F10.6,1pE20.12)') 'istatus = 1, no interpolation'
         call error_handler(E_MSG, 'interp_cubed_sphere', string1)
      end if
   endif

enddo


end subroutine interp_cubed_sphere_distrib


   subroutine unit_square_location(quad, closest, location, lon_o,lat_o, found_quad, l,m, origin, closest_only)
!=======================================================================
!  subroutine unit_square_location(quad, closest, location, lon_o,lat_o , found_quad, l,m, origin)
!
! Subroutine based on http://www.particleincell.com/2012/quad-interpolation/.
! The idea is to derive a mapping from any convex quadrilateral(x,y) onto a unit square (l,m).
! Also map the location of the ob onto that square.
! This is a bilinear interpolation,
! x = a0 + a1*l*m + a2*m + a3*l
! y = b0 + b1*l*m + b2*m + b3*l
! so (?) does not take into account the curvature of the quadrilateral.
! That has been handled by the intermediate mapping from (lon,lat) to (x,y).
!
! A higher order method exists (Nagata 2005: Simple Local Interpolation of Surfaces
! Using Normal Vectors) to map curved quadrilaterals onto the unit square,
! but the inverse map cannot be done analytically(?), so is not developed here.
! Instead another mapping is made before this one: the lon,lat coordinates
! of the corners/nodes are converted to a flat planar coordinate system by using
! the distances and directions from one node to the other three.  See create_cs_grid_arrays.
! The relative error in the distances appears to be <~ 2.E-4, which is avoided
! by defining a planar coordinates system for each corner of each quad.
! Then the ob is never near the 'far edges', where distance distortion could be a problem.

integer,             intent(in)    :: quad, closest
type(location_type), intent(in)    :: location
real(r8),            intent(in)    :: lon_o, lat_o
logical,             intent(in)    :: closest_only           !long debug temp arg
real(r8),            intent(out)   :: l,m
integer,             intent(inout) :: found_quad
integer,             intent(out)   :: origin

! Observation location in the planar space.
real(r8) :: x_o, y_o
! Locations in unit square space, and coefficients of quadratic equation for m.
real(r8) :: aa, bb, cc, det, angle, d, bearing_o
integer  :: oc(1)

l = MISSING_R8;    m = MISSING_R8

! Observation:
! Map the location of the ob into the planar space

! Figure out which corner of cell is the closest to the ob.
! Used to get the correct x_ax_bearing and a and b coeffs.
oc = minloc(corners(quad,:), mask = corners(quad,:) == closest)
origin = oc(1)

! The bearing of the observation.
bearing_o = bearing(lon_rad(closest),lat_rad(closest),lon_o*DEG2RAD,lat_o*DEG2RAD )

! Calculate the difference of the ob bearing from x_axis of this quad.
! The order is opposite of what might be expected because bearings are measured clockwise,
! while angles are measured counterclockwise.
angle = x_ax_bearings(origin,quad) - bearing_o
! Normalize angle to -pi<angle<pi.
angle = mod(angle,PI) - PI*int(angle/PI)
! Calculate the distance from this quad's origin.
d = get_dist(cs_locs(closest), location, no_vert=.true.)
x_o = d * cos(angle)
y_o = d * sin(angle)

! Coefficients of the quadratic equation for m for this quad.
aa = a(1,origin,quad)*b(2,origin,quad) - a(2,origin,quad)*b(1,origin,quad)
bb = a(3,origin,quad)*b(2,origin,quad) - a(1,origin,quad)*y_o             + b(1,origin,quad)*x_o
cc =                                   - a(3,origin,quad)*y_o
if (aa .eq.0._r8 ) then
   write(string1,'(A,I6,1x,1p4E12.4)') 'quad, Cannot divide by a(1)*b(2) - a(2)*b(1) : ',quad,  &
                            a(1,origin,quad),b(2,origin,quad), &
                            a(2,origin,quad),b(1,origin,quad)
   call error_handler(E_ERR, 'unit_square_location', string1,source,revision,revdate)
end if
if ( a(3,origin,quad) .eq. 0._r8 .and. do_out) then
   write(string1,'(A,I6,1x,1p4E12.4)') 'quad, cc=0: a(3)*y_o : ',quad, a(3,origin,quad)
   call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
end if

! Calculate m from the binomial equation, given the quadratic equation coefficients.
det = bb*bb - 4._r8*aa*cc
if (det >= 0._r8 ) then
   m = (-bb + sqrt(det))/(2._r8*aa)
else
   write(string1,'(A,I6,1X,1p4E12.4)') 'quad, angle, d, x_o, y_o',quad, angle, d, x_o, y_o
   call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
   write(string1,'(A,1p4E12.4)') 'Cannot divide by a(1)*b(2) - a(2)*b(1) : ',         &
                            a(1,origin,quad),b(2,origin,quad), &
                            a(2,origin,quad),b(1,origin,quad)
   call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
   write(string1,'(A,1p4E12.4)') 'det<0: a(3)*y_o : ', a(3,origin,quad)
   call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
   write(string1, '(A,1p2E11.3,A,I8)') 'b^2-4ac < 0: can not interpolate ',x_o,y_o,' for cell ',quad
   call error_handler(E_ERR, 'unit_square_location', string1, source, revision, revdate)
end if

if ( m < 0._r8 .or. m > 1._r8) then
   if (found_quad == 0) then
      ! Return with found_quad still = failure (0) to test the next quad.
      return
   else if (found_quad == 1) then
      ! Exit with error if m is outside valid range.
      write(string1, *) 'location of ob in unit square is out of bounds m = [0,1] ',m
      call error_handler(E_ERR, 'unit_square_location', string1, source, revision, revdate)
   else
      write(string1, *) 'found_quad does not have a valid input value of 0 or 1'
      call error_handler(E_ERR, 'unit_square_location', string1, source, revision, revdate)
   end if
end if

! Use m to calculate the 2nd unit square coordinate value.
l = (x_o -a(2,origin,quad)*m) / (a(3,origin,quad) + a(1,origin,quad)*m)
! ? closest_only not needed in this subroutine?
if (l .eq. 0._r8 .and. do_out .and. closest_only) then
   write(string1,'(A,I6,1X,1p4E12.4)') 'quad, x_o - a(2)*m = ',quad, x_o ,a(2,origin,quad),m
   call error_handler(E_MSG, 'unit_square_location', string1,source,revision,revdate)
end if

! Exit with error if l is outside valid range.
if (l < 0._r8 .or. l > 1._r8) then
   if (found_quad == 0) then
      return
   else if (found_quad == 1) then
      ! Exit with error if l is outside valid range.
      write(string1, *) 'location of ob in unit square is out of bounds l = [0,1] ',l
      call error_handler(E_ERR, 'unit_square_location', string1, source, revision, revdate)
   else
      write(string1, *) 'found_quad does not have a valid input value of 0 or 1'
      call error_handler(E_ERR, 'unit_square_location', string1, source, revision, revdate)
   end if
end if

! If we get this far, then this quad contains the ob.
! Return with found_quad = 1; success
found_quad = 1

return

end subroutine unit_square_location


   subroutine interp_lonlat_distrib(state_ens_handle, location, obs_kind, interp_val, istatus)
!=======================================================================
!
! Find the 4 corners of the lon-lat grid cell that encloses an ob at 'location'
! and interpolate the values of obs_kind to that location.

type(ensemble_type), intent(in) :: state_ens_handle
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_kind
integer,            intent(out) :: istatus(:)
real(r8),           intent(out) :: interp_val(:)

integer  :: i, e
real(r8) :: bot_lon, top_lon, delta_lon,                                &
            lon_below, lat_below, lat_above, lev_below,                 &
            lon_fract, lat_fract, temp_lon,            &
            lon_lat_lev(3), level, pressure, height

real(r8), allocatable :: val_11(:), val_12(:), val_21(:), val_22(:), a(:, :)
integer,  allocatable :: vstatus(:), istatus_distrib(:)

integer  :: s_type, s_type_01d,s_type_2d,s_type_3d,   &
            lon_ind_below, lon_ind_above, lat_ind_below, lat_ind_above, &
            num_lons
character (len=8)   :: lon_name, lat_name, lev_name
integer   :: ens_size
integer, allocatable :: track_vstatus(:)

! Would it be better to pass state as prog_var_type (model state type) to here?
! As opposed to the stripped state vector. YES. This would give time interp.
! capability here; do we really want this or should it be pushed up?

! istatus   meaning                  return expected obs?   assimilate?
! 0         obs and model are fine;  yes                    yes
! 1         fatal problem;           no                     no
! 2         exclude valid obs        yes                    no
! 3         unfamiliar obs type      no                     no

! These are fields which were observed, and will have 3d locations, but the
! corresponding state-vector component could, conceivably, be missing one of the dimensions.
! The only use for such fields I have thought of is parameterization tuning.
! Such fields would not have observations associated with them.
! for now I will assume that observed fields are not missing any dimensions.
! PS is missing lev, although we have never assimilated those obs.

! model_interpolate will continue to use state passed to it;
!    recalc p_col and model_h columns as needed.
!    no need to convert to a standard vert coord; no distance calc involved.

ens_size = state_ens_handle%num_copies -5 ! Now calculating mean copy also
allocate(val_11(ens_size),val_12(ens_size), val_21(ens_size), val_22(ens_size))
allocate(a(ens_size, 2))
allocate(vstatus(ens_size), istatus_distrib(ens_size))
allocate(track_vstatus(ens_size))

! Start with no errors in 
istatus(:) = 0
vstatus(:) = 0
istatus_distrib(:) = 0
val_11 = MISSING_R8
val_12 = MISSING_R8
val_21 = MISSING_R8
val_22 = MISSING_R8

! Get the observation (horizontal) position, in degrees
! In location_mod location%{lon,lat} are stored as radians, but get_location converts them.
lon_lat_lev = get_location(location)

! check whether model_mod can interpolate the requested variable
! Pressure (3d) can't be specified as a state vector field (so s_type will = MISSING_I),
! but can be calculated for CAM, so obs_kind = QTY_PRESSURE is acceptable.
s_type = dart_to_cam_kinds(obs_kind)

if (s_type == MISSING_I .and. &
   (obs_kind .ne. QTY_PRESSURE) .and.  (obs_kind .ne. QTY_SURFACE_ELEVATION)) then
   istatus = 3
   interp_val = MISSING_R8
   write(*,*) 'Wrong type of obs = ', obs_kind
   return
end if

! Get lon and lat grid specs

! Set [lon,lat,lev] names to a default, which will be overwritten for variables
! in the state vector, but not for other acceptable variables (3D pressure, surface
! elevation, ...?)
lon_name = 'lon     '
lat_name = 'lat     '
! ? How to separate the 3D from 2D 'other' variables?
!   Can't do it automatically/generically because they're not part of state vector
!   and that info isn't coming from DART.
if (obs_kind .eq. QTY_SURFACE_ELEVATION) then
   lev_name = 'none    '
else if (obs_kind .eq. QTY_PRESSURE) then
   lev_name = 'lev     '
end if

! There can't be any 0d or 1d ob fields, so lump them together for elimination in this search.
! CS; what about earth rotation obs?  Are they converted into 2D fields?
!     And PS is a 1d field in SE-CAM.
s_type_01d = state_num_0d + state_num_1d
! Positions within the rank 2 and 3 fields
s_type_2d = s_type - s_type_01d
s_type_3d = s_type_2d - state_num_2d

if (s_type == MISSING_I .and. &
   (obs_kind .eq. QTY_PRESSURE) .or.  (obs_kind .eq. QTY_SURFACE_ELEVATION)) then
   ! use defaults lon_name and lat_name set above
else if (s_type <= state_num_0d + state_num_1d) then
   ! error; can't deal with observed variables that are 0 or 1D in model_mod.
   istatus = 3
   interp_val = MISSING_R8
   write(*,*) 'Cannot handle 0 or 1d state vars, s_type = ', s_type
   return
else if (s_type_2d > 0 .and. s_type_2d <= state_num_2d) then
   lon_name = dim_names(s_dimid_2d(1,s_type_2d))
   lat_name = dim_names(s_dimid_2d(2,s_type_2d))
   lev_name = 'none    '
else if (s_type_3d > 0 .and. s_type_3d <= state_num_3d) then
   lon_name = dim_names(s_dimid_3d(2,s_type_3d))
   lat_name = dim_names(s_dimid_3d(3,s_type_3d))
   lev_name = dim_names(s_dimid_3d(1,s_type_3d))
else
   istatus = 3
   interp_val = MISSING_R8
   write(*,*) 'Unexpected state type value, s_type = ', s_type
   return
end if

!------------------------------------------------------------------------------
! Gack!
! staggered longitudes; slon (4x5 fv grid) = [-2.5, 2.5,...,352.5]  !
!                        lon ( "         ) = [    0.,  5.,...,  355.]
!! This is a problem for lon = 359, for example.  It's not in the range of slon.
!------------------------------------------------------------------------------

! Compute bracketing lon indices
! Define a local longitude to deal with CAM's weird staggered longitude grid.
temp_lon = lon_lat_lev(1)

if (lon_name == 'lon     ') then
   num_lons  = lon%length
   bot_lon   = lon%vals(1)
   top_lon   = lon%vals(num_lons)
   delta_lon = lon%vals(2) - lon%vals(1)
else if (lon_name == 'slon    ') then
   num_lons  = slon%length
   bot_lon   = slon%vals(1)
   top_lon   = slon%vals(num_lons)
   delta_lon = slon%vals(2) - slon%vals(1)
   ! Make certain longitudes conform to the weird CAM staggered grid.
   if ((lon_lat_lev(1) - top_lon) >= delta_lon) temp_lon = lon_lat_lev(1) - 360._r8
end if

! Remove = (top_lon) to prevent lon_above from being out of bounds
! if (temp_lon >= bot_lon .and. temp_lon <= top_lon) then
if (temp_lon >= bot_lon .and. temp_lon   <  top_lon) then
   ! adding the 1 makes up for subtracting the bot_lon.
   lon_ind_below = int((temp_lon - bot_lon) / delta_lon) + 1
   lon_ind_above = lon_ind_below + 1
   lon_fract = (temp_lon - ((lon_ind_below - 1) * delta_lon + bot_lon)) / delta_lon
else
! At wraparound point
   lon_ind_above = 1
   lon_ind_below = num_lons
! never happens; lon starts with 0, slon starts with -2.5 for 4x5 FV grid
!   if (temp_lon < bot_lon) temp_lon = temp_lon + 360.0_r8
   lon_fract = (temp_lon - top_lon) / delta_lon
end if


! Next, compute neighboring lat rows
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE MAY BE WRONG
! Inefficient search used for latitudes in Gaussian grid. Might want to speed up.
! CAM-FV; lat = -90., ...   ,90.
!        slat =   -88.,...,88.

call coord_index(lat_name, lon_lat_lev(2), lat_ind_above, lat_ind_below)

if (lat_ind_above == lat_ind_below) then
   if (lat_ind_above == 1) then
      lat_fract = 0.0_r8
   else                     !both must be equal to the max (s)lat index
      lat_fract = 1.0_r8
   end if
else
   if (lat_ind_above < lat_ind_below) then
      ! switch order
      i = lat_ind_above
      lat_ind_above = lat_ind_below
      lat_ind_below = i
   end if
   ! only lat_xxx is changed by these calls
   call coord_val(lat_name, lat_ind_below, lon_below, lat_below, lev_below)
   call coord_val(lat_name, lat_ind_above, lon_below, lat_above, lev_below)
   lat_fract = (lon_lat_lev(2) - lat_below) / (lat_above - lat_below)
end if

! Now, need to find the values for the four corners
! determine the vertical coordinate: model level, pressure, or height
! Future?; this assumes that obs with a vertical location have 2 horizontal locations too.
!          The state vector may have fields for which this isn't true, but no obs we've seen
!          so far violate this assumption.  It would have to be a synthetic obs, like some
!          sort of average.
if (obs_kind == QTY_SURFACE_ELEVATION) then
   call error_handler(E_ERR, 'model_interpolate ', 'QTY_SURFACE_ELEVATION not done')

else if (vert_is_level(location)) then
   call error_handler(E_ERR, 'model_interpolate ', 'vert_is_level')

else if (vert_is_pressure(location)) then
   ! which_vert is pressure for this obs
   pressure = lon_lat_lev(3)
   call get_val_pressure_distrib(val_11, state_ens_handle, lon_ind_below, lat_ind_below, pressure, obs_kind,vstatus)
   track_vstatus = vstatus

   call get_val_pressure_distrib(val_12, state_ens_handle, lon_ind_below, lat_ind_above, pressure, obs_kind,vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 )  track_vstatus(e) = vstatus(e)
   enddo

   call get_val_pressure_distrib(val_21, state_ens_handle, lon_ind_above, lat_ind_below, pressure, obs_kind,vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 ) track_vstatus(e) = vstatus(e)
   enddo

   call get_val_pressure_distrib(val_22, state_ens_handle, lon_ind_above, lat_ind_above, pressure, obs_kind,vstatus)
   do e = 1, ens_size
      if (vstatus(e) /= 0 ) track_vstatus(e) = vstatus(e)
   enddo
   vstatus = track_vstatus

else if (vert_is_surface(location)) then
   call error_handler(E_ERR, 'model_interpolate ', 'vert_is_surface')
! Need option for vert_is_undefined
else
   write(*,*) '   No vert option chosen!'

end if

! HK loop around ensembles
do e = 1, ens_size
   ! lat is already converted to degrees by get_location
   if (abs(lon_lat_lev(2)) > max_obs_lat_degree .and. vstatus(e) /= 1) then
      istatus_distrib(e) = 4
   else
      istatus_distrib(e) = vstatus(e)
   end if

   ! indices of vals are (longitude, latitude)
   if (istatus_distrib(e) /= 1) then
      a(e, 1) = lon_fract * val_21(e) + (1.0_r8 - lon_fract) * val_11(e)
      a(e, 2) = lon_fract * val_22(e) + (1.0_r8 - lon_fract) * val_12(e)

      interp_val(e) = lat_fract * a(e, 2) + (1.0_r8 - lat_fract) * a(e, 1)

   else
      interp_val(e) = MISSING_R8
   end if

enddo

! Set the element of ps that's tested elsewhere back to MISSING_R8, to signal
! other routines to calculate the ps arrays for themselves
! Currently (10/26/06) this flag is not used.
! ps(1,1) = MISSING_R8

!> @todo Sort out istatus
istatus = istatus_distrib

end subroutine interp_lonlat_distrib


! Pobs
   subroutine get_val_level(val, st_vec, lon_index, lat_index, level, obs_kind, istatus)
!=======================================================================
!   subroutine get_val_level(val, st_vec, lon_index, lat_index, level, obs_kind, istatus)
!
! Gets the value on level for variable obs_kind
! at lon_index, lat_index horizontal grid point
!
! written by Kevin Raeder, based on code from Hui Liu 4/28/2006 and get_val_pressure
!         from Jeff Anderson
!
! This version excludes observations below lowest level and above
! highest level.

real(r8), intent(out) :: val
real(r8), intent(in)  :: st_vec(:)
integer,  intent(in)  :: lon_index, lat_index, level, obs_kind
integer,  intent(out) :: istatus

integer               :: vstatus, num_levs
real(r8)              :: p_surf, threshold

! No errors to start with
istatus = 0
vstatus = 0

! This assumes that all variables are defined on model levels, not on interface levels.
num_levs = dim_sizes(find_name('lev     ',dim_names))

! Interpolate in vertical to get two bounding levels
if (level > num_levs .or. level < 1) then
   ! Exclude obs below the model's lowest level and above the highest level
   istatus = 1
   val = MISSING_R8
else
   if (highest_obs_level == MISSING_R8) then
      ! To do this completely right; p_surf would depend on whether obs_kind was on a staggered grid
      ! and highest_obs_level would be recalculated for each location passed in.
      ! Since the hybrid coord system is pure pressure at the top levels, I'll ignore these for now.
      ! 3Dp:  No, st_vec should be used, not ens_mean.  And p has not(?) been
      ! allocated or initialized to the st_vec values by the model_interpolate call.
      call get_val(p_surf, st_vec, lon_index, lat_index, -1, QTY_SURFACE_PRESSURE, vstatus)
      call plevs_cam (p_surf, num_levs, p_col)
      ! p_col(1:num_levs) = p(1:num_levs,lon_index,lat_index)

      highest_obs_level = 1.0_r8
      threshold = highest_obs_pressure_Pa
      do while ((p_col(nint(highest_obs_level))) < threshold )
         highest_obs_level = highest_obs_level + 1.0_r8
      end do
   end if

   if (level < nint(highest_obs_level)) then
      ! Exclude from assimilation the obs above a user specified level
      ! but still calculate the expected obs.
      istatus = 2
   else
      istatus = 0
   end if

   if (obs_kind == QTY_PRESSURE) then
      ! Can't get the value from get_val because 3d pressure is not a model variable.
      ! Can calculate it from ps.

      ! ps is on A-grid, so no need to check for staggered grids
      ! 3Dp: See 'No' comment, above.
      call get_val(p_surf, st_vec, lon_index, lat_index, -1, QTY_SURFACE_PRESSURE, vstatus)
      if (vstatus > 0) then
         val = MISSING_R8
         istatus = 1
         return
      end if
      ! Next, get the values on the levels for this ps
      call plevs_cam (p_surf, num_levs, p_col)

      val = p_col(level)
      ! val = p(level,lon_index,lat_index)
      ! 3Dp end
   else
       call get_val(val, st_vec, lon_index, lat_index, level, obs_kind, vstatus)
   end if

   if (vstatus /= 0) then
      istatus = 1
      val = MISSING_R8
   end if
end if

end subroutine get_val_level
! Pobs end


!> Distributed version of get_val_pressure
subroutine get_val_pressure_distrib(val, state_ens_handle, lon_index, lat_index, pressure, obs_kind, istatus)
!=======================================================================
!
! Gets the vertically interpolated value on pressure for variable obs_kind
! at lon_index, lat_index horizontal grid point
!
! This version excludes observations below lowest level pressure and above
! highest level pressure.

real(r8),            intent(out) :: val(:)
type(ensemble_type), intent(in)  :: state_ens_handle
real(r8),            intent(in)  :: pressure
integer,             intent(in)  :: lon_index, lat_index, obs_kind
integer,             intent(out) :: istatus(:)

real(r8), allocatable :: bot_val(:), top_val(:), p_surf(:), frac(:)
real(r8), allocatable :: ps_local(:, :) ! HK What is this for?
integer, allocatable  :: bot_lev(:), top_lev(:)
integer               :: i, vstatus, num_levs, slon_index
integer               :: fld_index !HK What is this for?
integer               :: ens_size, e

ens_size = state_ens_handle%num_copies -5
slon_index = find_name('slon    ',dim_names)

allocate(bot_val(ens_size), top_val(ens_size), p_surf(ens_size), frac(ens_size))
allocate(ps_local(ens_size, 2))
!> @todo HK I don't know why you need two values, one is just + 
allocate(bot_lev(ens_size), top_lev(ens_size))

! No errors to start with
istatus = 0
vstatus = 0

! Need to get the surface pressure at this point.
! Find out whether the observed field is a staggered in CAM.
! Check whether the state vector has wind components on staggered grids, i.e. whether CAM is FV.
! find_name returns 0 if the field name is not found in the cflds list.
! Add more staggered variables later?
! Can I make a more generic test; loop over all QTY_s, then check whether any of the
!    associated dimensions are staggered?   Sounds too expensive to be worth it. . .?
! ? Implement 3Dp here?  or should/can it not use the ens mean PS field?
fld_index   = find_name('PS      ',cflds)
i = index_from_grid(1,lon_index,lat_index,  fld_index)
call get_state(ps_local(:, 1), i, state_ens_handle)


if (obs_kind == QTY_U_WIND_COMPONENT .and. find_name('US      ', cflds) /= 0) then
   ! ps defined on lat grid (-90...90, nlat = nslat + 1),
   !    need it on staggered lat grid, which starts half a grid spacing north
   ! What about poles?
   !    It's OK; ps is defined on the 'lat' grid, which runs [-90...90],
   !    while staggered US is on the 'slat' grid, defined only *inside* this range.
   !p_surf = ps_stagr_lat(lon_index, lat_index)

   i = index_from_grid(1,lon_index,lat_index+1,fld_index)
   call get_state(ps_local(:, 2), i, state_ens_handle)
   p_surf      = (ps_local(:, 1) + ps_local(2, :))* 0.5_r8

else if (obs_kind == QTY_V_WIND_COMPONENT .and. find_name('VS      ', cflds) /= 0) then
   ! lon =     0...     255 (for 5 degree grid)
   !slon = -2.5 ... 252.5
   if (lon_index == slon%length) then
      i = index_from_grid(1,1,          lat_index ,fld_index)
   else
      i = index_from_grid(1,lon_index+1,lat_index ,fld_index)
   end if
   call get_state(ps_local(:, 2), i, state_ens_handle)
   p_surf      = (ps_local(:, 1) + ps_local(2, :))* 0.5_r8
else
   ! A-grid ps can be retrieved from state vector, which was used to define ps on entry to
   ! model_interpolate.
   ! PSx4 Change this to use PS passed in.
   ! p_surf = ps(lon_index, lat_index)
   p_surf     = ps_local(:, 1)
end if
!   if (vstatus > 0) then
!      val = MISSING_R8
!      istatus = 1
!      return
!   end if

! Next, get the pressures on the levels for this ps
! Assuming we'll only need pressures on model mid-point levels, not interface levels.
! This pressure column will be for the correct grid for obs_kind, since p_surf was taken
!     from the grid-correct ps[_stagr] grid

num_levs = dim_sizes(find_name('lev     ',dim_names))
allocate(p_col_distrib(ens_size, num_levs))
call plevs_cam_distrib(p_surf, num_levs, p_col_distrib, ens_size)

do e = 1, ens_size

   if (pressure <= p_col_distrib(e,1) .or. pressure >= p_col_distrib(e,num_levs)) then
      ! Exclude obs below the model's lowest level and above the highest level
      ! We *could* possibly use ps and p(num_levs) to interpolate for points below the lowest level.
      !    if (print_details) PRINT*,'get_val_pressure: pressure is out of range'
      istatus(e) = 1
      val(e) = MISSING_R8

   else

      if (pressure < highest_obs_pressure_Pa) then
         ! Exclude from assimilation the obs above a user specified level
         ! but still calculate the expected obs.
         istatus(e) = 2
      endif

      do i =2, num_levs
         if (pressure < p_col_distrib(e,i)) then
            top_lev(e) = i -1
            bot_lev(e) = i
            frac(e) = (p_col_distrib(e,i) - pressure) / (p_col_distrib(e,i) - p_col_distrib(e, i-1))
            exit
         endif
      enddo
   endif
enddo


! Pobs
if (obs_kind == QTY_PRESSURE) then
   ! can't get pressure on levels from state vector; get it from p_col instead
   ! get_val_pressure is called for 4 different columns, which will have different p_cols
   ! for each ps is on A-grid, so no need to check for staggered grids
   do e = 1, ens_size
      if(istatus(e) == 0 .or. istatus(e) == 2 ) then
         bot_val(e) = p_col_distrib(e, bot_lev(e))
         top_val(e) = p_col_distrib(e, top_lev(e))
         val(e) = (1.0_r8 - frac(e)) * bot_val(e) + frac(e) * top_val(e)
      endif
   enddo
!      if (abs((val - pressure)/val) > 1.0E-12) then
!         ! We're looking for a pressure on a model level, which is exactly what p_col provides,
!! NOT HERE; that happens in get_val_level
!         write(string1, *) 'val /= pressure = ',val,pressure,' when val is a P obs '
!         call error_handler(E_WARN, 'get_val_pressure', string1, source, revision, revdate)
!      end if
else
! Pobs end

  ! need to grab values for each bot_val
  do e = 1, ens_size ! HK you only need to do this for distinct bot_vals
     if(istatus(e) == 0  .or. istatus(e) == 2) then
        call get_val_distrib(bot_val, state_ens_handle, ens_size, lon_index, lat_index, bot_lev(e), obs_kind, vstatus)
        if (vstatus == 0) call get_val_distrib(top_val, state_ens_handle, ens_size, lon_index, lat_index, top_lev(e), obs_kind, vstatus)
        if (vstatus == 0) then
            val(e) = (1.0_r8 - frac(e)) * bot_val(e) + frac(e) * top_val(e)
        else
           istatus(e) = 1
           val(e) = MISSING_R8
        end if

     endif

  enddo

end if


deallocate(bot_val, top_val, p_surf, frac, p_col_distrib)
deallocate(bot_lev, top_lev)

end subroutine get_val_pressure_distrib



!> Distributed version of get height
subroutine get_val_height_distrib(val, state_ens_handle, lon_index, lat_index, height, location, obs_kind, istatus)
!=======================================================================
!
! Gets the vertically interpolated value on height for variable obs_kind
! at lon_index, lat_index horizontal grid point
!
! written by Kevin Raeder, based on code from Hui Liu 4/28/2006 and get_val_pressure
!         from Jeff Anderson
!
! This version excludes observations below lowest level height and above
! highest level height.
! 

real(r8),            intent(out) :: val(:)
type(ensemble_type), intent(in)  :: state_ens_handle
real(r8),            intent(in)  :: height !HK Just one height?
integer,             intent(in)  :: lon_index, lat_index, obs_kind
type(location_type), intent(in)  :: location
integer,             intent(out) :: istatus(:)

real(r8), allocatable :: bot_val(:), top_val(:), p_surf(:), frac(:)
integer,  allocatable :: top_lev(:), bot_lev(:)
real(r8), allocatable :: ps_local(:,:)
integer,  allocatable :: vstatus(:)
integer               :: i, num_levs
integer               :: fld_index, ind
integer               :: ens_size, e

logical  :: stagr_lon, stagr_lat

ens_size = state_ens_handle%num_copies -5

allocate(model_h_distrib(ens_size, num_levs))
allocate(bot_val(ens_size), top_val(ens_size), p_surf(ens_size), frac(ens_size))
allocate(p_col_distrib(ens_size, num_levs))
!> @todo HK I don't know why you need two values, one is just + 1 to the other
allocate(bot_lev(ens_size), top_lev(ens_size))
allocate(ps_local(ens_size, 2))
allocate(vstatus(ens_size))

! No errors to start with
istatus(:) = 0
vstatus(:) = 0
stagr_lon = .false.
stagr_lat = .false.

! Need to get the surface pressure at this point for dcz2. 
! Assuming we'll only need pressures on model mid-point levels, not interface levels.
num_levs = dim_sizes(find_name('lev     ',dim_names))

print*, '**** HEIGHT HEIGHT HIEGHT ****'

! Need to get the surface pressure at this point. 
! Check whether the state vector has wind components on staggered grids, i.e. whether CAM is FV.
! find_name returns 0 is the field name is not found in the cflds list.
! See get_val_press for more documentation.
fld_index   = find_name('PS      ',cflds)
ind         = index_from_grid(1,lon_index,lat_index,  fld_index)
call get_state(ps_local(:, 1), i, state_ens_handle)

if     (obs_kind == QTY_U_WIND_COMPONENT .and. find_name('US      ', cflds) /= 0) then
   !p_surf = ps_stagr_lat(lon_index, lat_index)
   stagr_lat = .true.
   ind         = index_from_grid(1,lon_index,lat_index+1,fld_index)
   call get_state(ps_local(:, 2), i, state_ens_handle)
   p_surf      = (ps_local(:, 1) + ps_local(:, 2))* 0.5_r8
else if (obs_kind == QTY_V_WIND_COMPONENT .and. find_name('VS      ', cflds) /= 0) then
   !p_surf = ps_stagr_lon(lon_index, lat_index)
   stagr_lon = .true.
   if (lon_index == slon%length) then
      ind       = index_from_grid(1,1,          lat_index ,fld_index)
   else
      ind       = index_from_grid(1,lon_index+1,lat_index ,fld_index)
   end if
   call get_state(ps_local(:, 2), i, state_ens_handle)
   p_surf      = (ps_local(:, 1) + ps_local(:, 2))* 0.5_r8
else
   p_surf = ps_local(:, 1)
end if

! Next, get the heights on the levels for this ps

! merge/MPI
! We want to use the new vec for each new ob on height because the state was updated 
! for all previous obs, and we want to use the most up to date state to get the best location.

call model_heights_distrib(state_ens_handle, ens_size, p_surf, num_levs, location, model_h_distrib, vstatus)

!HK I don't think the above vstatus is used.

! debug
! write(logfileunit,'(A,6F7.0,/(10F7.0))') 'heights = ',(model_h(i), i=1,num_levs)

! Interpolate in vertical to get two bounding levels

do e = 1, ens_size

   if (height >= model_h_distrib(e, 1) .or. height <= model_h_distrib(e, num_levs)) then
      ! Exclude obs below the model's lowest level and above the highest level
      istatus(e) = 1
      val(e) = MISSING_R8

   else
      ! This should be redefined every time(?), not just for the first (arbitrary) entry.
      !   if (highest_obs_height_m == MISSING_R8) then
      call plevs_cam_distrib(p_surf, num_levs, p_col_distrib, ens_size) ! does this need to be the distribued version
      do i=1,num_levs
         if (p_col_distrib(e,i) > highest_obs_pressure_Pa) then
            highest_obs_height_m = model_h_distrib(e, i)
            goto 10
         endif
      end do

10    if (height > highest_obs_height_m ) then
         ! Exclude from assimilation the obs above a user specified level
         ! but still calculate the expected obs.
         istatus(e) = 2
      else
         istatus(e) = 0
      endif

      ! Search down through heights
      do i = 2, num_levs
         if (height > model_h_distrib(e, i)) then
            top_lev(e) = i -1
            bot_lev(e) = i
            frac(e) = (model_h_distrib(e, i) - height      ) / &
                   (model_h_distrib(e, i) - model_h_distrib(e, i-1))
            goto 21
         endif
      end do

      21 continue

    endif

enddo

! Pobs HK what is the point of all these Pobs comments?
if (obs_kind == QTY_PRESSURE) then
   ! Observing a pressure on a height surface sounds silly.  But for completeness:
   ! get_val_height is called for 4 different columns, which will have different p_cols for each.
   ! It's also requested by obs_def_gps_mod.

   ! Next, get the values on the levels for this ps
   ! ps is on A-grid, so no need to check for staggered grids
   print*, '****** I do not think you have to call plevs_cam again'
   !call plevs_cam_distrib(p_surf, num_levs, p_col, ens_size) !HK Why are you calling this again?

   do e = 1, ens_size
      bot_val(e) = p_col_distrib(e, bot_lev(e))
      top_val(e) = p_col_distrib(e, top_lev(e))
      val(e) = (1.0_r8 - frac(e)) * bot_val(e) + frac(e) * top_val(e)
   enddo

else

   do e = 1, ens_size
      if(istatus(e) == 0  .or. istatus(e) == 2) then

         call get_val_distrib(bot_val, state_ens_handle, ens_size, lon_index, lat_index, bot_lev(e), obs_kind, vstatus(1))
         if (vstatus(1) == 0) call get_val_distrib(bot_val, state_ens_handle, ens_size, lon_index, lat_index, bot_lev(e), obs_kind, vstatus(1))
         if (vstatus(1) == 0) then
            val(e) = (1.0_r8 - frac(e)) * bot_val(e) + frac(e) * top_val(e)
         else
            istatus(e) = 1
            val(e) = MISSING_R8
         endif

      endif
   enddo

end if

deallocate(bot_val, top_val, p_surf, frac, p_col_distrib)
deallocate(bot_lev, top_lev, vstatus)
deallocate(model_h_distrib) !HK I don't think you want to do this.

end subroutine get_val_height_distrib


   subroutine get_val(val, st_vec, lon_index, lat_index, level, obs_kind, istatus)
!=======================================================================
! subroutine get_val(val, st_vec, lon_index, lat_index, level, obs_kind, istatus)
!

real(r8), intent(out) :: val
real(r8), intent(in)  :: st_vec(:)
integer, intent(in)   :: lon_index, lat_index, level, obs_kind
integer, intent(out)  :: istatus

integer :: indx, field_type

! No errors to start with
istatus = 0

field_type = dart_to_cam_kinds(obs_kind)
if (field_type <= 0 .or. &
    field_type > state_num_0d + state_num_1d + state_num_2d + state_num_3d) then
   istatus = 1
   val = 0.0_r8
   return
end if
indx = index_from_grid(level, lon_index, lat_index, field_type)
!PRINT*,'  get_val; indx = ',indx,' for lev,lon_i,lat_i,field_type ',  &
!       level, lon_index, lat_index, field_type

val = st_vec(indx)

return

end subroutine get_val

!=======================================================================
subroutine get_val_distrib(val, state_ens_handle, ens_size, lon_index, lat_index, level, obs_kind, istatus)

!> how may pieces of state to grab
integer,             intent(in)  :: ens_size
real(r8),            intent(out) :: val(ens_size)
type(ensemble_type), intent(in)  :: state_ens_handle
integer,             intent(in)  :: lon_index, lat_index, level, obs_kind
integer,             intent(out) :: istatus

integer :: indx, field_type

! HK get_val_distrib only fails if the index is out of the grid, that is why there is only one return code
! not ens_size of them.

! No errors to start with
istatus = 0

field_type = dart_to_cam_kinds(obs_kind)
if (field_type <= 0 .or. &
    field_type > state_num_0d + state_num_1d + state_num_2d + state_num_3d) then
   istatus = 1
   val = 0.0_r8
   return
end if

indx = index_from_grid(level, lon_index, lat_index, field_type)
!val = x(indx)
call get_state(val, indx, state_ens_handle)

return

end subroutine get_val_distrib

! End of model_interpolate section

!#######################################################################

! Vector-field translations


   subroutine prog_var_to_vector(var, st_vec)
!=======================================================================
! subroutine prog_var_to_vector(var, st_vec)
!

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

! Start copying fields to straight vector
indx = 0

!  0d variables
do nf = 1, state_num_0d
   indx = indx + 1
   st_vec(indx) = var%vars_0d(nf)
end do

!  1d variables
do nf = 1, state_num_1d
   do i=1,s_dim_1d(nf)
      indx = indx + 1
      st_vec(indx) = var%vars_1d(i, nf)
   end do
end do

!  Surface pressure and other 2d ; load by first dimension;
!      lon and/or lat and/or lev, in that precedence, as loaded in to vars_2d in trans_coord
do nf = 1, state_num_2d
   do j=1,s_dim_2d(2,nf)
   do i=1,s_dim_2d(1,nf)  !levs or lons  both reads and writes will be contiguous in this case
      indx = indx + 1
      st_vec(indx) = var%vars_2d(i, j, nf)
   end do
   end do
end do

!  3D fields, loaded by columns (note the coordinate order).
!  This is also looping over latitude values in the outer loop,
!  so if there is some spatial searching; all the pole points will be closer together
!  in memory than in the previous/Bgrid structure.
do nf= 1, state_num_3d

!   if (print_details .and. do_out) then
!      write(string1, '(A,4I5)') 'fld, nlons, nlats, nlevs ',nf &
!                          ,s_dim_3d(2,nf),s_dim_3d(3,nf),s_dim_3d(1,nf)
!      call error_handler(E_MSG, 'prog_var_to_vector', string1, source, revision, revdate)
!   end if

   do i=1,s_dim_3d(3,nf)   !lats
   do j=1,s_dim_3d(2,nf)   !lons
   do k=1,s_dim_3d(1,nf)   !levs  both reads and writes will be contiguous in this case
      indx = indx + 1
      st_vec(indx) = var%vars_3d(k,j,i, nf)
   end do
   end do
   end do
end do

! Temporary check
if (indx /= model_size) then
   write(string1, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'prog_var_to_vector', string1, source, revision, revdate)
end if

end subroutine prog_var_to_vector




   subroutine vector_to_prog_var(st_vec, var)
!=======================================================================
! subroutine vector_to_prog_var(st_vec, var)
!

real(r8),         intent(in)  :: st_vec(:)
type(model_type), intent(out) :: var

integer :: i, j, k, nf, indx

! Start copying fields from straight vector
indx = 0

! 0d arrays
do nf = 1, state_num_0d
   indx = indx + 1
   var%vars_0d(nf) = st_vec(indx)
end do

!  1d arrays
do nf = 1, state_num_1d
   do i=1,s_dim_1d(nf)
      indx = indx + 1
      var%vars_1d(i, nf) = st_vec(indx)
   end do
end do

!  Surface pressure and other 2d fields
do nf = 1, state_num_2d
   do j = 1, s_dim_2d(2,nf)
   do i = 1, s_dim_2d(1,nf)
      indx = indx + 1
      var%vars_2d(i, j, nf) = st_vec(indx)
   end do
   end do
end do

! 3D fields; see comments in prog_var_to_vect
do nf = 1, state_num_3d
!   if (print_details .and. do_out) then
!      write(string1, '(A,4I5)') 'fld, nlons, nlats, nlevs ',nf &
!                       ,s_dim_3d(2,nf),s_dim_3d(3,nf),s_dim_3d(1,nf)
!      call error_handler(E_MSG, 'vector_to_prog_var', string1, source, revision, revdate)
!   end if
   do i = 1, s_dim_3d(3,nf)
   do j = 1, s_dim_3d(2,nf)
   do k = 1, s_dim_3d(1,nf)
      indx = indx + 1
      var%vars_3d(k,j,i, nf) = st_vec(indx)
   end do
   end do
   end do
end do

! Temporary check
if (indx /= model_size) then
   write(string1, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'vector_to_prog_var', string1, source, revision, revdate)
end if

end subroutine vector_to_prog_var



! End of Vector-field translations

!#######################################################################

subroutine get_close_obs_distrib(filt_gc, base_obs_loc, base_obs_type, locs, kinds, &
                            num_close, close_indices, distances, state_ens_handle)
!----------------------------------------------------------------------------
!
! get_close_obs takes as input an observation location, a DART TYPE (not KIND),
! and a list of all potential locations and KINDS on this task.
!
! get_close_obs 
!    *) converts vertical coordinates as needed,
!    *) calls location_mod/threed_sphere:get_close_obs,
!       to which it sends this (converted) array of locations,
!    *) gets back the distances and indices of those locations that are
!       "close" to the base ob.
!    *) tests for being above the highest_obs_pressure_Pa threshold, 
!       and increases the vertical distance based on height above highest_*.
!
! get_close_obs will use the ensemble average to convert the obs and/or state 
!               vertical location(s) to a standard (pressure) vertical location
!
! Possible Ideas/Notes for Future Development.
!    3D model_h would be useful here; calc once and use over and over.
!       Reinstall height/lon slice code for model_heights to facilitate that.
!          CS; height/lon slice code wouldn't be helpful for cubed sphere grid.
!    3D pressures also useful here;
!       Reinstall height/lon slice code for plevs_cam to facilitate that.
!          CS; height/lon slice code wouldn't be helpful for cubed sphere grid.
!    throw away ens_mean after it's been used (or don't worry about it for now).

type(ensemble_type) :: state_ens_handle

type(get_close_type),         intent(in)    :: filt_gc
type(location_type),          intent(in)    :: base_obs_loc
type(location_type),          intent(inout) :: locs(:)
integer,                      intent(in)    :: base_obs_type, kinds(:)
integer,                      intent(out)   :: num_close, close_indices(:)
real(r8),                     intent(out)   :: distances(:)

! remove some (unused) variables?
integer                :: k, t_ind
integer                :: base_which, local_base_which, obs_which, local_obs_which
real(r8), dimension(3) :: base_array, local_base_array, obs_array, local_obs_array
real(r8)               :: increment, threshold, thresh_wght
type(location_type)    :: local_base_obs_loc, local_locs
integer                :: base_obs_kind

! If base_obs vert type is not pressure; convert it to pressure
base_which    = nint(query_location(base_obs_loc))
base_array    = get_location(base_obs_loc)
base_obs_kind = get_quantity_for_type_of_obs(base_obs_type)

if (base_which == VERTISPRESSURE) then
   if (vert_coord == 'pressure') then
      local_base_obs_loc = base_obs_loc
      local_base_array   = get_location(base_obs_loc)  ! needed in num_close loop
      local_base_which   = base_which
   else if (vert_coord == 'log_invP') then
      ! change so that ln pressure is the vertical coordinate:
      call convert_vert_distrib(state_ens_handle, base_array, base_which, local_base_array, local_base_which, base_obs_loc, base_obs_kind)
      local_base_obs_loc = set_location(base_array(1), base_array(2), local_base_array(3), &
                                        local_base_which)
   end if
else

   !Upgrading convert_vert to use field profiles at the actual ob location is
   !probably not worthwhile: that approx horiz location of the obs is used only to
   !convert its vert coord to pressure (if necessary),
   !which, in turn, is used to modify the distance if the ob or model variable
   !is higher than highest_XXX_mb.  That modification tapers to 0,
   !so any errors introduced by this approx will be continuous and random,
   !introducing no bias.

   call convert_vert_distrib(state_ens_handle, base_array, base_which, local_base_array, local_base_which, base_obs_loc, base_obs_kind)
   ! CS; Set arguments must be in degrees.
   local_base_obs_loc = set_location(base_array(1), base_array(2), local_base_array(3), &
                                     local_base_which)
end if

! Get all the potentially close obs but no distances (optional argument dist(:) is not present)
call loc_get_close_obs(filt_gc, local_base_obs_loc, base_obs_type, locs, kinds, &
                       num_close, close_indices)

if (vert_coord == 'pressure') then
   threshold = highest_state_pressure_Pa
else if (vert_coord == 'log_invP') then
   threshold = log(P0%vals(1)/highest_state_pressure_Pa)
end if

if (threshold > 0.0_r8) thresh_wght = 1._r8/(threshold * threshold)

do k = 1, num_close

   ! CS; not the 2D index into the nodes/corners.
   ! The indices in close_obs refer to the subset of (state) vars or obs ON 1 TASK.
   ! That subset is (re)labeled 1...num_vars_task#, where num_vars_task# ~ state_vec_size/ntasks.
   ! So those indices can't tell me which state vector element I'm dealing with.
   ! I need to use the location of each returned close_indices to learn anything about it.
   ! So how did SE30_GMexHiRes_Exp13 work so well when I assumed that close_indices
   ! has global indices in it.

   t_ind = close_indices(k)
   obs_array = get_location(locs(t_ind))
! query_location returns location%which_vert, if not 'attr' argument is given.
   obs_which = nint(query_location(locs(t_ind)))

   if (obs_which == VERTISPRESSURE ) then
      if (vert_coord == 'pressure') then
         ! put the vertical (pressure) of the state/ob in local storage
         local_obs_array(3) = obs_array(3)
         local_obs_which    = obs_which
      else if (vert_coord == 'log_invP') then

         call convert_vert_distrib(state_ens_handle, obs_array, obs_which, local_obs_array, local_obs_which,  locs(t_ind), kinds(t_ind))

         ! save the converted location back into the original list.
         ! huge improvement in speed since we only do the vertical convert
         ! once per location, instead of num_close * nobs times.
         locs(t_ind) = set_location( local_obs_array(1), local_obs_array(2), &
                                     local_obs_array(3), local_obs_which)

      end if
   else
      ! Convert vertical coordinate of locs to pressure.
      ! If horiz_dist_only is true, the vertical location and which_vert aren't used by get_dist,
      ! but need to be defined for set_loc and are used in the damping section below no matter what.
      ! CS 2nd line of args added for coord_ind_cs.
      ! D2C changed order of args for easier reading.
      call convert_vert_distrib(state_ens_handle, obs_array, obs_which, local_obs_array, local_obs_which, &
                        locs(t_ind), kinds(t_ind))

      ! save the converted location back into the original list.
      ! huge improvement in speed since we only do the vertical convert
      ! once per location, instead of num_close * nobs times.
      locs(t_ind) = set_location( local_obs_array(1), local_obs_array(2), &
                                  local_obs_array(3), local_obs_which)

      ! obs_which = -2 (VERTISUNDEF) means this ob is vertically close to base_obs, no matter what.
      if (local_obs_array(3) == MISSING_R8) then
         local_obs_array(3) = local_base_array(3)
         local_obs_which = local_base_which
      end if
   end if

   ! CS; Set arguments must be in degrees.
   local_locs = set_location(obs_array(1), obs_array(2), local_obs_array(3), &
                                   local_obs_which)

!  nsc fri, 13mar09
!  allow a namelist specified kind string to restrict the impact of those
!  obs kinds to only other obs and state vars of the same kind.
   if ((impact_kind_index >= 0)                .and. &
       (impact_kind_index == base_obs_type)    .and. &
       (impact_kind_index /= kinds(t_ind))) then
      distances(k) = 999999._r8     ! arbitrary very large distance
   else if (local_base_which == VERTISUNDEF) then
      ! The last argument, no_vert = .true., makes get_dist calculate horizontal distance only.
      distances(k) = get_dist(local_base_obs_loc, local_locs, base_obs_type, kinds(t_ind),.true.)
      ! Then no damping can be done since vertical distance is undefined.
      ! ? Is this routine called *both* to get model points close to a real obs,
      !   AND ob close to a model point?  I want damping in the latter case,
      !   even if ob has which_vert = VERTISUNDEF.
      !   I think that testing on local_base_which will do that.
   else
      distances(k) = get_dist(local_base_obs_loc, local_locs, base_obs_type, kinds(t_ind))

      ! Damp the influence of obs (which are below the namelist variable highest_obs_pressure_Pa)
      ! on variables above highest_state_pressure_Pa.
      ! This section could also change the distance based on the QTY_s of the base_obs and obs.

      ! distances = 0 for some for synthetic obs.
      ! Additive increase, based on height above threshold, works better than multiplicative

      ! See model_mod circa 1/1/2007 for other damping algorithms.

      ! for vertical = ln(psurf/p) the obs is ABOVE the cutoff if
      ! local_obs_array(3)-threshold > 0. So reverse the difference.
      if (vert_coord == 'pressure') then
         increment = threshold - local_obs_array(3)
      else if (vert_coord == 'log_invP') then
         increment = local_obs_array(3) - threshold
      end if
      ! This if-test handles the case where no damping is performed, i.e.
      ! highest_state_pressure_Pa = 0 and threshold = 0.
      if (increment > 0) then
         distances(k) = distances(k) + increment * increment * thresh_wght
         ! too sharp      distances(k) = distances(k) + increment / threshold
      end if
   end if

end do

end subroutine get_close_obs_distrib


!=======================================================================
   subroutine convert_vert_distrib(state_ens_handle, old_array, old_which, new_array, new_which, old_loc, old_kind)
!=======================================================================
!
! Uses model information and subroutines to convert the vertical location of an ob
! (prior, model state variable, or actual ob) into the standard vertical coordinate (pressure).
! Called by model_mod:get_close_obs.
! Kevin Raeder 10/26/2006

! CS I need to find the closest horizontal grid point to an ob
! for convert_vert, before the 3D call to get_close_obs.
! This is done for FV/Eul using the monotonicity of the lon and lat arrays,
! which we do not have for CS.
! So I need a call to get_close_obs (do not worry about a small search distance
! around the ob specified) and ask for distances.
! get_close_obs will compare the which_vert of the ob with that of the old_loc
! and do horizontal distance if they are different, 3D distance if they are the same.

type(ensemble_type),    intent(in)    :: state_ens_handle
real(r8), dimension(3), intent(in)    :: old_array
integer,                intent(in)    :: old_which
real(r8), dimension(3), intent(inout) :: new_array
integer,                intent(out)   :: new_which
type(location_type),    intent(in)    :: old_loc
integer,                intent(in)    :: old_kind

! "Pass through" variables needed by coord_ind_cs

integer   :: num_levs, top_lev, bot_lev
integer   :: quad_corners(4), closest
integer   :: lon_ind, lat_ind, kind_indx
real(r8)  :: p_surf, frac, l, m, p_surf_corner(4)

character(len=8) :: cam_varname

! this code does not alter the lat/lon, only the vertical.
! but still return a full location for subsequent use.
new_array(1) = old_array(1)
new_array(2) = old_array(2)

! these should be set by the code below; it's an error if not.
new_which    = MISSING_I
new_array(3) = MISSING_R8
num_levs     = lev%length

if (old_which == VERTISPRESSURE .or. old_which == VERTISHEIGHT  .or. &
    old_which == VERTISLEVEL    .or. old_which == VERTISSURFACE .or. &
    old_which == VERTISUNDEF   ) then
   ! Proceed
else
   ! fatal error - there should be no other options for vert.
   write(string1,'(''obs at '',3(F9.5,1x),I2,'' has bad vertical type'')') &
                   old_array, old_which
   call error_handler(E_ERR, 'convert_vert', string1,source,revision,revdate)
end if

! Need lon and lat indices to select ps for calc of p_col for vertical conversion.
! This is no longer an approximation; P is interpolated to the actual ob location.

if (old_which == VERTISLEVEL ) then
   ! Need the PS at the model grid point.  Get from location.
   kind_indx = dart_to_cam_kinds(old_kind)
   if ( kind_indx < 0 ) then
      write(string1,*)'old_kind  is ',old_kind,' | kind_indx is ',kind_indx
      write(string2,*)'get_name_for_quantity of old_kind ', trim(get_name_for_quantity(old_kind))
      call error_handler(E_ERR,'convert_vert',string1,source,revision,revdate,text2=string2)
   end if

   if (l_rectang) then
      ! This is not as general as in earlier model_mods; assumes 2D obs ! locations are (lon, lat)
      ! (and 3D are (lev,lon,lat) but what else could they be?).

      cam_varname = trim(cflds(kind_indx))
      if (cam_varname == 'US') then
         call coord_index('lon     ', old_array(1), lon_ind)
         call coord_index('slat    ', old_array(2), lat_ind)
         !p_surf = ps_stagr_lat(lon_index, lat_index)
         p_surf = 0.5*(get_surface_pressure(state_ens_handle, lon_ind, lat_ind) + &
                     get_surface_pressure(state_ens_handle, lon_ind, lat_ind +1) )
         ! WHAT ABOUT FIELDS THAT MIGHT COME ON ilevS ?   have lev_which_dimid from above;
         !     test = ilev%dim_id or lev%dim_id
         ! Next get the values on the levels for this ps
      ! unnecessary complication?   num_levs = dim_sizes(find_name('lev     ',dim_names))
         call plevs_cam (p_surf, num_levs, p_col) !HK I think this can just be plevs_cam
      else if (cam_varname == 'VS') then
         call coord_index('slon    ', old_array(1), lon_ind)
         call coord_index('lat     ', old_array(2), lat_ind)
         print*, 'can this work 1'
         p_surf = ps_stagr_lon(lon_ind,lat_ind)
         call plevs_cam(p_surf, num_levs, p_col)
      else
         call coord_index('lon     ', old_array(1), lon_ind)
         call coord_index('lat     ', old_array(2), lat_ind)
         print*, 'can this work 2'
         p_surf = ps(lon_ind,lat_ind)
         ! 3Dp
         print*, 'can this work 3'
         p_col(1:num_levs) = p(1:num_levs,lon_ind,lat_ind)

      end if
   else
      ! Cubed sphere; more complicated search for indices of this location.
      ! 3Dp Is this coord_ind_cs call necessary?
      ! The 3D index into the state vector is NOT known. HK - Yes it is.
      ! We have the 3D index into the subset of state vars ON 1 TASK.
      ! That subset has its own indexing, which is what's available here.
      ! The relationship to global indices is stuck back in filter.
      ! pmo-assim divergence debug
! Reduced output bug?  filter no longer worked after I commented these (and others) out
      ! if (do_out) write(*,'(A)') '       convert_vert; calling coord_ind_cs for VERTISLEVEL'
      call coord_ind_cs(old_loc, old_kind, cs_locs, cs_kinds, quad_corners, l, m, closest, .true.)
! obs_which not needed?
!                         .true., old_which)
      !p_surf = ps(closest,1)
      p_surf = get_surface_pressure(state_ens_handle, closest, 1)
      !HK p_col(1:num_levs) = p(1:num_levs,closest,1)
      !>@todo Is this a 3D pressure array?
      call plevs_cam(p_surf, num_levs, p_col)
   end if
else
   !>@todo What is this code doing? Can you do this without calling interp?
   ! Find ps at the ob point.  Need to interpolate.
   if (l_rectang) then
       call error_handler(E_ERR, 'not done for l_rectang', 'yet. Why do you need to call interploate?')
   !   ! Only interested in P (columns), so don't need to worry about staggered grids here.
   !   call interp_lonlat_distrib(state_ens_handle, old_loc, QTY_SURFACE_PRESSURE, all_psurf, istatus)
   !   p_surf = all_psurf(ens_size) !>@todo Sort this out only need the mean
   else
      !call interp_cubed_sphere_distrib(state_ens_handle, old_loc, QTY_SURFACE_PRESSURE, all_psurf, istatus)
      !p_surf = all_psurf(ens_size) !>@todo Sort this out only need the mean

      !HK idea: get the psurf at each of the indices and interpolate?
      ! Stolen from the interpolate part of model_interpolate
      call coord_ind_cs(old_loc, old_kind, cs_locs, cs_kinds, quad_corners, l, m, closest, .false.)
      p_surf_corner(1) = get_surface_pressure(state_ens_handle, quad_corners(1), 1)
      p_surf_corner(2) = get_surface_pressure(state_ens_handle, quad_corners(2), 1)
      p_surf_corner(3) = get_surface_pressure(state_ens_handle, quad_corners(3), 1)
      p_surf_corner(4) = get_surface_pressure(state_ens_handle, quad_corners(4), 1)
      p_surf  = p_surf_corner(2) *        l *       m  &
              + p_surf_corner(1) * (1._r8-l)*       m  &
              + p_surf_corner(4) * (1._r8-l)*(1._r8-m) &
              + p_surf_corner(3) *        l *(1._r8-m)

   end if
   ! 4 means that ob is beyond lat_max from namelist, so go ahead, ignoring it.
   ! HK convert_vert uses the mean only - this is messy
   !>@todo No status to check now, what should you do?
   !if (istatus(ens_size) /= 0 .and. istatus(ens_size) /= 4) then
   !   write(string1,'(A,I8)') 'inter_X failed for QTY_SURFACE_PRESSURE.  istatus = ',istatus
   !   call error_handler(E_ERR, 'convert_vert', string1,source,revision,revdate)
   !end if


end if

! Need the vertical pressure structure for this column

! Convert vertical coordinate from one of the following to pressure.
! integer, parameter :: VERTISUNDEF    = -2 ! has no vertical location (undefined)
! integer, parameter :: VERTISSURFACE  = -1 ! surface value
! integer, parameter :: VERTISLEVEL    =  1 ! by level
! integer, parameter :: VERTISPRESSURE =  2 ! by pressure
! integer, parameter :: VERTISHEIGHT   =  3 ! by height

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
else if (old_which == VERTISSURFACE ) then
   ! surface field; change which_vert for the distance calculation
   if (vert_coord == 'pressure') then
      new_array(3) =  p_surf
   else if (vert_coord == 'log_invP') then
    ! new_array(3) = log(p_surf/p_surf)
      new_array(3) = 0._r8
   end if
   new_which = VERTISPRESSURE
elseif (old_which == VERTISPRESSURE) then
   if (vert_coord == 'log_invP') then
      new_array(3) = log(p_surf/old_array(3))
      new_which = VERTISPRESSURE
   end if
else if (old_which == VERTISLEVEL ) then
   ! WHAT ABOUT FIELDS THAT MIGHT COME ON ilevS ?   have lev_which_dimid from above;
   !     test = ilev%dim_id or lev%dim_id
   ! OR do this for all columns in static_init_model_dist, which would make PS (and P)
   ! globally available for all regions?
   if (vert_coord == 'pressure') then
      new_array(3) =            p_col(nint(old_array(3)))
   else if (vert_coord == 'log_invP') then
      new_array(3) = log(p_surf/p_col(nint(old_array(3))))
   end if
   new_which = VERTISPRESSURE
else if (old_which == VERTISHEIGHT) then
   call error_handler(E_ERR, 'convert_vert_distrib ', 'height not done')

end if

return

end subroutine convert_vert_distrib

!------------------------------------------------------------
! Subroutines from mpas_atm/model_mod.f90, for using cartesian coordinates to
! find closest node to an ob.

subroutine init_closest_center()

! use nCells(ncol), latCell, lonCell to initialize a GC structure
! to be used later in find_closest_cell_center().

! set up a GC in the locations mod

integer :: i

allocate(cs_locs_xyz(ncol))

! real(r8), allocatable :: lonCell(:) ! cell center longitudes (degrees, original radians in file)
! real(r8), allocatable :: latCell(:) ! cell center latitudes  (degrees, original radians in file)
do i=1, ncol
   cs_locs_xyz(i) = xyz_set_location(lon%vals(i), lat%vals(i), 0.0_r8, radius)
end do

! the width (2nd arg of ...init) really isn't used anymore, but it's part of the
! interface so we have to pass some number in.
call xyz_get_close_maxdist_init(cs_gc_xyz, 1.0_r8)
call xyz_get_close_obs_init    (cs_gc_xyz, ncol, cs_locs_xyz)

end subroutine init_closest_center

!------------------------------------------------------------

function find_closest_cell_center(lat, lon)

! Determine the cell index for the closest center to the given point
! 2D calculation only.

real(r8), intent(in)  :: lat, lon
integer               :: find_closest_cell_center

type(xyz_location_type) :: pointloc
integer ::  closest_cell, rc
! real(r8) :: closest_dist
logical, save :: search_initialized = .false.

! do this exactly once.
if (.not. search_initialized) then
   call init_closest_center()
   search_initialized = .true.
end if

pointloc = xyz_set_location(lon, lat, 0.0_r8, radius)

call xyz_find_nearest(cs_gc_xyz, pointloc, cs_locs_xyz, closest_cell, rc)

! decide what to do if we don't find anything.
if (rc /= 0 .or. closest_cell < 0) then
   if (do_out) print *, 'cannot find nearest cell to lon, lat: ', lon, lat
   find_closest_cell_center = -1
   return
end if

! this is the cell index for the closest center
find_closest_cell_center = closest_cell

end function find_closest_cell_center

!------------------------------------------------------------

subroutine finalize_closest_center()

! get rid of storage associated with GC for cell centers.

call xyz_get_close_obs_destroy(cs_gc_xyz)

end subroutine finalize_closest_center


! End of get_close_obs section

!#######################################################################

! Initial conditions for DART

  subroutine pert_model_state(state, pert_state, interf_provided)
!=======================================================================
! subroutine pert_model_state(state, pert_state, interf_provided)
!
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

type(random_seq_type)   :: random_seq
type(model_type)        :: var_temp
integer                 :: i, j, k, m, pert_fld, mode, field_num
integer                 :: dim1, dim2, dim3, member
real(r8)                :: pert_val
integer, save           :: seed
logical                 :: perturbed

! FIX for 1D 0D  fields?

! perturb model state vector values for the filter initial conditions.
! the input is a single model state vector that has (different) gaussian
! noise added to each member to generate an initial ensemble.

interf_provided = .true.

call init_model_instance(var_temp)
call vector_to_prog_var(state,var_temp)

! If first call, then initialize a seed to use for initializing random sequences.
if (first_pert_call) then
   ! this line generates a unique base number, and subsequent calls add 1
   ! each time (which happens if there are multiple ensemble members/task).
   ! it is assuming there are no more than 1000 ensembles/task, which seems safe
   ! given the current sizes of state vecs and hardware memory.  this will make
   ! the results reproduce for runs with the same number of MPI tasks.  it will
   ! NOT give the same random sequence if you change the task count.
   seed = (my_task_id()+1) * 1000
   first_pert_call = .false.
end if

! After using the seed, increment by one so if this routine is called again
! for a different ensemble member it will generate a different set of nums.
call init_random_seq(random_seq, seed)
seed = seed + 1

pert_fld = 1

Vars2Perturb : do while (pert_names(pert_fld) /= '        ')

   ! Keep track of whether or not this field is matched and was perturbed.
   perturbed = .false.

   ExistingVars : do m=1,nflds

      if (pert_names(pert_fld) /= cflds(m)) cycle ExistingVars

      perturbed = .true.

      call error_handler(E_MSG,'pert_model_state', 'Perturbing '//trim(pert_names(pert_fld)))

      ! FIXME : not robust. ens_member is always 0 in CESM context.
      !  Choose mode of perturbations/resets;
      if (pert_sd(pert_fld) <= 0.0_r8 ) then
         ! Set each ensemble member to it's own constant value,
         ! as found in pert_base_vals(ens_member).
         ! This only works when setting a single field = to a different constant value
         ! for each ensemble member.
         ! Could add more fields by overloading pert_base_vals and
         ! adding code to find those values.
         mode = ens_member
      else
         ! Set each *field* to it's own pert_base_val +/- pert_sd
         mode = pert_fld
      end if

      if (m <= state_num_2d + state_num_1d + state_num_0d) then
         field_num = m - state_num_1d - state_num_0d
         dim1 = dim_sizes(s_dimid_2d(1,field_num))
         dim2 = dim_sizes(s_dimid_2d(2,field_num))

         ! reset base values to value provided in namelist.
         if ( pert_base_vals(mode) /= -888888.0d0 ) then
            if (print_details) then
               WRITE(     *     ,*) '   around new base value ',pert_base_vals(mode)
               WRITE(logfileunit,*) '   around new base value ',pert_base_vals(mode)
            endif
            var_temp%vars_2d(1:dim1,1:dim2,field_num) = pert_base_vals(mode)
         end if

         if (print_details) then
             WRITE(     *     ,'(A,A8,A3,1x,2E24.15)') 'org first and last state for ',&
                    cflds(m),' = ', var_temp%vars_2d(   1,   1,field_num), &
                                    var_temp%vars_2d(dim1,dim2,field_num)
             WRITE(logfileunit,'(A,A8,A3,1x,2E24.15)') 'org first and last state for ',&
                    cflds(m),' = ', var_temp%vars_2d(   1,   1,field_num), &
                                    var_temp%vars_2d(dim1,dim2,field_num)
         endif

         ! randomly perturb each point around its base value.
         if (pert_sd(pert_fld) > 0.0_r8 ) then
            do j = 1, dim2
            do i = 1, dim1
               pert_val = random_gaussian(random_seq, var_temp%vars_2d(i,j,field_num), &
                                          pert_sd(mode))
               var_temp%vars_2d(i,j,field_num) = pert_val
            end do
            end do
         end if

         if (print_details) then
            WRITE(     *     ,'(A,A8,A3,1x,2E24.15)') 'new first and last state for ',&
                    cflds(m),' = ', var_temp%vars_2d(   1,   1,field_num), &
                                    var_temp%vars_2d(dim1,dim2,field_num)
            WRITE(logfileunit,'(A,A8,A3,1x,2E24.15)') 'new first and last state for ',&
                    cflds(m),' = ', var_temp%vars_2d(   1,   1,field_num), &
                                    var_temp%vars_2d(dim1,dim2,field_num)
         endif

      else ! do the 3D fields

         field_num = m - state_num_2d - state_num_1d - state_num_0d
         dim1 = dim_sizes(s_dimid_3d(1,field_num))
         dim2 = dim_sizes(s_dimid_3d(2,field_num))
         dim3 = dim_sizes(s_dimid_3d(3,field_num))

         if (print_details) then
            WRITE(     *     ,'(A,A8,A3,1x,2E24.15)') 'org first and last state for ',&
                    cflds(m), ' = ', var_temp%vars_3d(   1,   1,   1,field_num),  &
                                     var_temp%vars_3d(dim1,dim2,dim3,field_num)
            WRITE(logfileunit,'(A,A8,A3,1x,2E24.15)') 'org first and last state for ',&
                    cflds(m), ' = ', var_temp%vars_3d(   1,   1,   1,field_num),  &
                                     var_temp%vars_3d(dim1,dim2,dim3,field_num)
         endif

         ! reset base values to value provided in namelist.
         if ( pert_base_vals(mode) /= -888888.0d0 ) then
            if (print_details) then
               WRITE(     *     ,*) '  perturbed around new base value ',pert_base_vals(mode)
               WRITE(logfileunit,*) '  perturbed around new base value ',pert_base_vals(mode)
            endif
            var_temp%vars_3d(1:dim1,1:dim2,1:dim3,field_num) = pert_base_vals(mode)
         end if

         ! randomly perturb each point around its base value.
         if (pert_sd(pert_fld) > 0.0_r8 ) then
            if (print_details) then
               WRITE(     *     ,*) 'Perturbing base value of ',cflds(m),' by st dev ',pert_sd(mode)
               WRITE(logfileunit,*) 'Perturbing base value of ',cflds(m),' by st dev ',pert_sd(mode)
            endif
            do j = 1, dim3
            do i = 1, dim2
            do k = 1, dim1
               ! pert_val = rand#(O(0-1)) * standard dev  + mean
               pert_val = random_gaussian(random_seq, var_temp%vars_3d(k,i,j,field_num), &
                                          pert_sd(mode))
               var_temp%vars_3d(k,i,j,field_num) = pert_val
            end do
            end do
            end do
         end if

         if (print_details) then
            WRITE(     *     ,'(A,A8,A3,1x,2E24.15)') 'new first and last state for ',&
                    cflds(m), ' = ', var_temp%vars_3d(   1,   1,   1,field_num),  &
                                     var_temp%vars_3d(dim1,dim2,dim3,field_num)
            WRITE(logfileunit,'(A,A8,A3,1x,2E24.15)') 'new first and last state for ',&
                    cflds(m), ' = ', var_temp%vars_3d(   1,   1,   1,field_num),  &
                                     var_temp%vars_3d(dim1,dim2,dim3,field_num)
         endif

      end if

   end do ExistingVars

   if (.not. perturbed) then
      write(string1,*)trim(pert_names(pert_fld)),' not found in list of state variables.'
      write(string2,*)'but was supposed to be used to perturb.'
      call error_handler(E_ERR,'pert_model_state', string1, source, revision, revdate, text2=string2)
   end if

   pert_fld = pert_fld + 1

end do Vars2Perturb

call prog_var_to_vector(var_temp,pert_state)
call end_model_instance(var_temp)

end subroutine pert_model_state


   subroutine init_conditions(st_vec)
!=======================================================================
! subroutine init_conditions(st_vec)
!
! Reads in restart initial conditions  -- noop for CAM

real(r8), intent(inout) :: st_vec(:)

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

   function index_from_grid(lev_ind, lon_ind, lat_ind, ifld)
!=======================================================================
!
! This could be more efficient by calculating and globally storing the beginning index
! of each ifld and/or the size of each ifld.

integer, intent(in) :: lev_ind, lon_ind, lat_ind, ifld

integer :: index_from_grid
integer :: i, j, fld_ind(2), done

index_from_grid = 0
done = 0

! Cycle through 0d state variables
do i=1,state_num_0d
   index_from_grid = index_from_grid + 1
   if (ifld == i) then
      return
   end if
end do
done = done + state_num_0d

! Cycle through 1d state variables
! Note that indices of fields can have varying dimensions.
do i=1,state_num_1d
   if (ifld - done == i) then
      if (dim_names(s_dimid_1d(i)) == 'lon     ' .or. &
          dim_names(s_dimid_1d(i)) == 'slon    ') index_from_grid = index_from_grid + lon_ind
      if (dim_names(s_dimid_1d(i)) == 'lat     ' .or. &
          dim_names(s_dimid_1d(i)) == 'slat    ') index_from_grid = index_from_grid + lat_ind
      if (dim_names(s_dimid_1d(i)) == 'lev     ' .or. &
          dim_names(s_dimid_1d(i)) == 'ilev    ') index_from_grid = index_from_grid + lev_ind
      ! CS lon_ind has been pirated by ncol.
      if (dim_names(s_dimid_1d(i)) == 'ncol     ') index_from_grid = index_from_grid + lon_ind
      return
   else
      index_from_grid = index_from_grid + s_dim_1d(i)
   end if
end do
done = done + state_num_1d

! Cycle through 2d state variables.
! Note that indices of fields can have varying dimensions.
do i=1,state_num_2d
   if (ifld - done == i) then
      ! We've found the desired field; now find index of lev and/or lon and/or lat
      do j=1,2
         if (dim_names(s_dimid_2d(j,i)) == 'lon     ' .or. &
             dim_names(s_dimid_2d(j,i)) == 'slon    '     ) fld_ind(j) = lon_ind
         if (dim_names(s_dimid_2d(j,i)) == 'lat     ' .or. &
             dim_names(s_dimid_2d(j,i)) == 'slat    '     ) fld_ind(j) = lat_ind
         if (dim_names(s_dimid_2d(j,i)) == 'lev     ' .or. &
             dim_names(s_dimid_2d(j,i)) == 'ilev    '     ) fld_ind(j) = lev_ind
         ! CS lon_ind has been pirated by ncol.
         if (dim_names(s_dimid_2d(j,i)) == 'ncol     ')     fld_ind(j) = lon_ind
      end do

      index_from_grid = index_from_grid + (fld_ind(2)-1) * s_dim_2d(1,i) + fld_ind(1)
      return
   else
      index_from_grid = index_from_grid +  s_dim_2d(2,i) * s_dim_2d(1,i)
   end if
end do
done = done + state_num_2d


! Cycle through 3d state variables
! Note that indices of fields can have varying dimensions.
! CS There won't be any 3d fields if ncol is a dimension.
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
   end if
end do

end function index_from_grid



   function find_name(nam, list)
!=======================================================================
! function find_nam(nam, list)

integer                                    :: find_name
character (len=*),              intent(in) :: nam
character (len=*), dimension(:),intent(in) :: list

integer :: i

find_name = 0
do i = 1,size(list,1)
if (nam == list(i)) then
   find_name = i
   return
end if
end do

end function find_name



   subroutine coord_val(dim_name, indx, lon_val, lat_val, lev_val)
!==========================================================================
! subroutine coord_val(dim_name, indx, lon_val, lat_val, lev_val)

! given the name of the coordinate to be searched and the index into that array,
! returns the coordinate value  in either lon_val, lat_val, or lev_val.
! All 3 _val arguments are present so that this routine can return the value
! in the coordinate that the calling routine wants it to be, and search/placement doesn't
! have to be done there.
! Used by get_state_meta_data and model_interpolate.

character (len=8), intent(in)    :: dim_name
integer,           intent(in)    :: indx
real(r8),          intent(inout) :: lon_val, lat_val, lev_val

! Check for acceptable value of indx?

if (dim_name == 'lon     ') lon_val = lon%vals(indx)
if (dim_name == 'lat     ') lat_val = lat%vals(indx)
if (dim_name == 'lev     ') lev_val = lev%vals(indx)
if (dim_name == 'slon    ') then
   lon_val = slon%vals(indx)
   ! CAM staggered longitude grid -2.5, ..., 352.5 (FV4x5)
   ! but DART wants to see 0.,...,360.  only.
   if (lon_val < 0._r8) lon_val = lon_val + 360._r8
end if
if (dim_name == 'slat    ') lat_val = slat%vals(indx)
if (dim_name == 'ilev    ') lev_val = ilev%vals(indx)
if (dim_name == 'ncol    ') then
   lon_val = lon%vals(indx)
   lat_val = lat%vals(indx)
end if
! Add more for other coords?  hyam...?  Not for now; never referenced indirectly

return

end subroutine coord_val


   subroutine coord_ind_cs(obs_loc, obs_kind, loc_list, kinds_list, &
                           quad_corners, l, m, closest, closest_only)
!==========================================================================
!  subroutine coord_ind_cs(obs_loc, obs_kind, loc_list, kinds_list, &
!                          quad_corners, l, m, closest, closest_only)

! Primary input/output variables
integer,              intent(out) :: quad_corners(4)
real(r8),             intent(out) :: l, m
integer,              intent(out) :: closest

! Variables needed by loc_get_close_obs:
type(location_type),  intent(in)  :: obs_loc,  loc_list(:)
integer,              intent(in)  :: obs_kind, kinds_list(:)
logical,              intent(in)  :: closest_only

! Output from loc_get_close_obs
integer  :: num_close
!
! It would be nice if these could be smaller, but I don't know what number would work.
! It has to be large enough to accommodate all of the grid points that might lie
! within 2xcutoff; resolution and location dependent.
! The size must be specified here; (:) yields an error, and 'allocatable' doesn't help.
integer  :: close_ind(ncol)
real(r8) :: dist(ncol)

! Local Variables
! dist_# in radians (can't be initialized here or it will get the 'save' property,
! and will not be reset during subsequent entries to this subroutine.
real(r8) :: dist_1, dist_2
real(r8) :: lon_lat_lev(3)
integer  :: k, k1, k2, found_quad, closest2, origin

dist_1 = 10._r8;     dist_2 = 10._r8

!? Is this an endless loop?:
!filter calls get_close_obs (the one in model_mod, not location), giving it gc
!model_mod takes the input, does its modifications and calls
!  location_mod:get_close_obs via 'use location: loc_get_close_obs => get_close_obs'
!-> get_close_obs
!-> convert_vert
!-> coord_ind_cs
!-> loc_get_close_obs  passes the cs_gc it got from here to the location_mod version of get_close_obs.
!-> result of location_mod:get_close_obs is used/passed back up,
!location_mod:get_close_obs does NOT use model_mod:get_close_obs.
!so no endless loop.

! If which_vert of obs_loc and loc_list are different,
! then only the horizontal distance will be calculated, which is what we want.
! When called from interp_cubed_sphere loc_list()%which_vert is UNDEFINED,
! while obs_loc will not be UNDEFINED.
! Check the values.
lon_lat_lev = get_location(obs_loc)

k1 = MISSING_I
closest = MISSING_I

! See whether this obs_ is a state variable
! long debug
! write(*,'(2(A,I8))')'num_close for ',obs_kind,' = ',num_close
! This could be done by 2 calls to minloc(dist), with the 2nd call using a mask
! to prevent finding the closest, which was found in the first call.
! But would those 2 intrinsic searches through dist be faster than my 1 explicit search?

if (closest_only) then
   ! Use xyz/cartesian coordinates to quickly find the closest node.
   closest = find_closest_cell_center(lon_lat_lev(2), lon_lat_lev(1))
   ! xyz check
   if (do_out) then
! Reduced output bug?  filter no longer worked after I commented these (and others) out
      ! write(*,'(A,2F15.10)') 'lon,lat of ob   = ',lon_lat_lev(1),   lon_lat_lev(2)
      ! write(*,'(A,I8,2F15.10)') 'closest, lon,lat = ',closest,lon%vals(closest),lat%vals(closest)
   end if
   !
   ! mla = minloc(dist(1:num_close))
   ! closest = close_ind(mla(1))
   ! If convert_vert only needs the closest node, don't find the l,m weights.
   return
else
   ! Look for the 2 closest nodes, using slower way of getting all of the close obs
   ! and searching for the 2 closest.
   ! Could I have an array of max_dist(ncol), which could be used to reset max_dist
   ! in cs_gc based on lon_lat_lev?   NO.
   ! So loc_get_close_obs is going to return lists that are 64x larger in the
   ! refined region than in the coarse region.
   call loc_get_close_obs(cs_gc, obs_loc, obs_kind, loc_list, kinds_list, &
                       num_close, close_ind, dist)

   do k = 1,num_close
      if (dist(k) < dist_2) then
         ! Replace 2nd with new one.
         k2 = k
         closest2 = close_ind(k)
         dist_2 = dist(k)
         if (dist_2 < dist_1) then
            ! Switch 1st and new 2nd.
            k2 = k1
            k1 = k
            dist_2 = dist_1
            dist_1 = dist(k)
            closest2 = closest
            closest  = close_ind(k)
         end if
      end if
   end do
end if

if (closest == MISSING_I) then
!    write(string1,'(A,3F10.2,2I4,1p2E12.4)')                  &
!         'lon_lat_lev, obs_kind, num_close, dist_1, dist_2 = ', &
!          lon_lat_lev, obs_kind, num_close, dist_1, dist_2
!    call error_handler(E_ERR, 'coord_ind_cs', string1,source,revision,revdate)
end if

! Find the quad which contains the ob.
! First search the quads around closest.
! If that fails, search the quads around closest2.
! The search consists of calling quad_interp with the ob location
! and letting it determine whether the ob location maps into the unit square.
! (Re-use num_close; we're done with old value.)
num_close = num_nghbrs(closest)

! long debug
! if (do_out) write(*,'(2(A,I8))')'coord_ind_cs: num_nghbrs(',closest,') = ',num_close
! found_quad signals a failure(0)/success(1) signal.
found_quad = 0
do k=1,num_close
! sort_bearings; use these to figure out which k is the right quad to send to unit_square_location.
! move unit_s out of loop.
! ? closest_only not used this far into the subroutine?
   call unit_square_location(centers(k,closest), closest, obs_loc,          &
                             lon_lat_lev(1),lon_lat_lev(2), found_quad, l,m, origin, closest_only)
   if (found_quad == 1) then
      exit
   end if
end do

! Try the 2nd closest point, if the first failed.
if (found_quad == 0 .and. closest2 /= MISSING_I) then
   num_close = num_nghbrs(closest2)
   ! long debug
!    if (do_out) write(*,'(2(A,I8))')'coord_ind_cs: num_nghbrs2(',closest,') = ',num_close
   do k=1,num_close
! sort_bearings; use these to figure out which k is the right quad to send to unit_square_location.
! move unit_s out of loop.
      call unit_square_location(centers(k,closest2), closest2, obs_loc,         &
                                lon_lat_lev(1),lon_lat_lev(2), found_quad, l,m, origin, closest_only)
      if (found_quad == 1) then
         ! Put '2nd closest' information into 'closest'.
         dist_1 = dist_2
         !dist_1 = temp
         k1 = closest
         closest  = closest2
         !cs_index = close_ind(k)
!          write(string1,'(A,2F10.7,2I8,1p2E12.4)') 'l, m, closest2, origin2 = ', l, m, closest, origin
!          call error_handler(E_MSG, 'coord_ind_cs', string1,source,revision,revdate)
         exit
      end if
   end do
end if

if (found_quad == 1) then
   ! Need to shift corners according to which was chosen as the origin corner
   ! in num_close loop, above.  The weighted interp calculation assumes, as in
   ! the create_cs_grid_arrays mapping scheme, that the origin node is corner 4.
   quad_corners(1:4) = cshift(corners(centers(k,closest),1:4), origin)
else
   ! Both closest nodes failed; abort
   write(string1, '(A,2I8,A,2F10.4)') &
         'Neither of the 2 closest nodes ',  k1,closest2, &
         ' is a corner of the quad containing ob at ', lon_lat_lev(1),lon_lat_lev(2)
   call error_handler(E_ERR, 'coord_ind_cs', string1,source,revision,revdate)
end if

return

end subroutine coord_ind_cs


   subroutine coord_index(dim_name, val, indx, other_indx)
!==========================================================================
! subroutine coord_index(dim_name, indx, val, indx, other_indx)

! Given the name of the (Eulerian or FV) coordinate to be searched and the value,
! Returns the index of the closest coordinate value.
! Optionally returns the next closest index too, which may be < or > the closest.
! Used by get_state_meta_data.
! CS See coord_ind_cs for more expensive search through non-rectangular coordinate.

character (len=8), intent(in)  :: dim_name
real(r8),          intent(in)  :: val
integer,           intent(out) :: indx
integer, optional, intent(out) :: other_indx

real(r8), pointer :: coord(:)
real(r8)          :: diff_upper, diff_lower, val_local, resol
integer           :: coord_len, i

val_local = val

if (dim_name == 'lon     ') then
   coord     => lon%vals
   coord_len =  lon%length
   resol     =  lon%resolution
else if (dim_name == 'lat     ') then
   coord     => lat%vals
   coord_len =  lat%length
   resol     =  lat%resolution
else if (dim_name == 'lev     ') then
   coord     => lev%vals
   coord_len =  lev%length
   resol     =  lev%resolution
else if (dim_name == 'slon    ') then
   coord     => slon%vals
   coord_len =  slon%length
   resol     =  slon%resolution
   ! Make sure longitudes conform to the weird CAM staggered grid.
   if ((val - coord(coord_len)) >= (coord(coord_len)-coord(coord_len-1))) &
      val_local = val_local - 360._r8
else if (dim_name == 'slat    ') then
   coord     => slat%vals
   coord_len =  slat%length
   resol     =  slat%resolution
else if (dim_name == 'ilev    ') then
   coord     => ilev%vals
   coord_len =  ilev%length
   resol     =  ilev%resolution
else
   ! should not happen; fatal error.
   write(string1, *) 'unexpected dim_name, ', trim(dim_name)
   call error_handler(E_ERR, 'coord_index', string1,source,revision,revdate)
end if

! further check?  for blunders check that coord(1) - val is smaller than coord(2) - coord(1), etc.
! Assumes that coordinates are monotonic; not true for hyam, hyai.  But we don't reference them.
! The first 2 if blocks work for latitudes and levels.  Longitudes must be handled in the calling
!    routine.
if (val_local <= coord(1)) then
   indx = 1
   if (present(other_indx)) other_indx = 1
   return
else if (val_local >= coord(coord_len)) then
   indx = coord_len
   if (present(other_indx)) other_indx = coord_len
   return
else
   if (resol > 0.0_r8) then
      ! temp output
      num_calced = num_calced + 1
      ! regularly spaced; calculate the index
      indx = NINT((val_local - coord(1))/resol) + 1

      if (present(other_indx)) then
         if (val_local > coord(indx)) then
            other_indx = indx + 1
         else
            other_indx = indx - 1
         end if
      end if
   else
      ! temp output
      num_searched = num_searched + 1
      ! IRregularly spaced (but monotonically increasing); search for the index
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
            end if
            return
         end if
      end do
   end if
end if

end subroutine coord_index



   subroutine set_ps_arrays (vec)
!=====================================================================

real(r8), intent(in) :: vec(:)

integer :: ind, m,n, slon_index, slat_index, lon_index, lat_index, lev_index, &
           fld_index_2d, fld_index, dim1, dim2, dim_lev, ncol_index

! CS  could use l_rectang instead
ncol_index = find_name('ncol    ',dim_names)

! Assign values to ps grids for use by the rest of the module.
! CS Find PS among state_names_1d.
! Assuming ps is the first 2D field in state_num_2d

if (ncol_index /= 0) then
   ! Non-rectangular grid; 1 horizontal dimension.

   !CS PS is a 1d field (cubed sphere).
   dim2    = 1
   fld_index_2d = find_name('PS      ',state_names_1d)
   dim1    = s_dim_1d(fld_index_2d)
   lev_index  = find_name('lev     ',dim_names)
   dim_lev = dim_sizes(lev_index)

   if (alloc_ps) then
      allocate (ps         (dim1, dim2))
      allocate (p (dim_lev, dim1, dim2))
      alloc_ps = .false.
   end if

   fld_index = find_name('PS      ',cflds)
   ind       = index_from_grid(1,1,1,fld_index) -1
!    write(*,'(A,3I8)') 'Filling p with dimensions ',dim_lev,dim1,dim2
   do m=1,dim1
      ind = ind + 1
      ps(m,1) = vec(ind)
      call plevs_cam(ps(m,1), dim_lev, p(1:dim_lev,m,1))
   end do
   if (do_out) then
      PRINT*,'Finished assignment of PS for ',dim1,' elements'
      do n=1,dim_lev
         write(*,'(3X,1p5E12.4)') (p(n,m,1),m=1,dim1,10000)
      end do
   end if
else
   ! Rectangular grid; 2 horizontal dimensions.
   lon_index  = find_name('lon     ',dim_names)
   lat_index  = find_name('lat     ',dim_names)
   slon_index = find_name('slon    ',dim_names)
   slat_index = find_name('slat    ',dim_names)
   lev_index  = find_name('lev     ',dim_names)

   fld_index_2d = find_name('PS      ',state_names_2d)
   dim2 = s_dim_2d(2,fld_index_2d)
   dim1 = s_dim_2d(1,fld_index_2d)
   dim_lev = dim_sizes(lev_index)

   if (alloc_ps) then
      allocate (ps(         dim_sizes(lon_index), dim_sizes(lat_index)))
      allocate (p (dim_lev, dim_sizes(lon_index), dim_sizes(lat_index)))
      if (slon_index /= 0) &
         allocate( ps_stagr_lon (dim_sizes(slon_index), dim_sizes( lat_index)))
      if (slat_index /= 0) &
         allocate( ps_stagr_lat (dim_sizes( lon_index), dim_sizes(slat_index)))
      alloc_ps = .false.
   end if

   ps = reshape(vec,(/dim1,dim2/))

   !fld_index    = find_name('PS      ',cflds)
   !ind = index_from_grid(1,1,1,fld_index) -1
   ! 3Dp; fill the p(:,:,:) array.
!    write(*,'(A,3I8)') 'Filling p with dimensions ',dim_lev,dim1,dim2
   do n=1,dim2
   do m=1,dim1
      call plevs_cam(ps(dim1,dim2), dim_lev, p(1:dim_lev,m,n))
   !    p(1:dim_lev,m,n) = p_col(1:dim_lev)
   end do
   end do
   ! write(*,'(A,1pe12.4)') 'p(1,      1,   1) = ',p(1,      1,   1)
   ! write(*,'(A,1pe12.4)') 'p(dim_lev,1,   1) = ',p(dim_lev,1,   1)
   ! write(*,'(A,1pe12.4)') 'p(1,      dim1,1) = ',p(1,      dim1,1)
   ! write(*,'(A,1pe12.4)') 'p(dim_lev,dim1,1) = ',p(dim_lev,dim1,1)
   ! write(*,'(A,1pe12.4)') 'p(1,      1,   dim2) = ',p(1,      1,   dim2)
   ! write(*,'(A,1pe12.4)') 'p(dim_lev,1,   dim2) = ',p(dim_lev,1,   dim2)
   ! write(*,'(A,1pe12.4)') 'p(1,      dim1,dim2) = ',p(1,      dim1,dim2)
   ! write(*,'(A,1pe12.4)') 'p(dim_lev,dim1,dim2) = ',p(dim_lev,dim1,dim2)

   ! ps(1:s_dim_2d(1,1), 1:s_dim_2d(2,1)) = var%vars_2d(1:s_dim_2d(1,1),1:s_dim_2d(2,1), 1)

   if (slon_index /= 0) then
      do n=1,dim_sizes( lat_index)
         ! The first element of slon is 1/2 grid obs *west* of lon(1) (which
         ! = 0. longitude)
         ! Make p[hi]s_stagr_lon line up with the slon array.
         ! The index to the second ps can be either slon or lon because num_slon = num_lon
         ps_stagr_lon(1,n) = .5 * (ps(1,n) + ps(dim_sizes(slon_index),n))
         do m=2,dim_sizes(slon_index)
            ps_stagr_lon(m,n) = .5 * (ps(m-1,n) + ps(m,n))
         end do
      end do
   end if

   if (slat_index /= 0) then
      do n=1,dim_sizes(slat_index)
         do m=1,dim_sizes(lon_index)
            ps_stagr_lat(m,n) = .5 * (ps(m,n) + ps(m,n+1))
         end do
      end do
   end if
end if

end subroutine set_ps_arrays


   subroutine plevs_cam (p_surf, num_levs, pmid )
!=======================================================================
!
! Purpose:
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure.
!
! Method:
!
! Author: B. Boville (plevs0),
!         Kevin Raeder modified  8/1/03 to use hy[ab][im] from within module
!         rather than in a common block, for use in DART,
!
!-----------------------------------------------------------------------

! coef; commented these out
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use pmgrid

! coef; commented these out
! #include <comhyb.h>

!-----------------------------------------------------------------------
real(r8), intent(in)  :: p_surf        ! Surface pressure (pascals)
integer,  intent(in)  :: num_levs
real(r8), intent(out) :: pmid(:)   ! Pressure at model levels
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer k             ! Longitude, level indices
!-----------------------------------------------------------------------
!
! Set midpoint pressures and layer thicknesses
!
! coef
do k=1,num_levs
   pmid(k) = hyam%vals(k)*P0%vals(1) + hybm%vals(k)*p_surf
end do

return

end subroutine plevs_cam

!> Distributed version of plevs_cam, nothing special, just works on an array
subroutine plevs_cam_distrib(p_surf, num_levs, pmid, ens_size)
!=======================================================================
!
! Purpose:
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure.
!
! Method:
!
! Author: B. Boville (plevs0), 
!         Kevin Raeder modified  8/1/03 to use hy[ab][im] from within module
!         rather than in a common block, for use in DART,
!
!-----------------------------------------------------------------------

! coef; commented these out
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use pmgrid

! coef; commented these out
! #include <comhyb.h>

!-----------------------------------------------------------------------
integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: p_surf(ens_size)        ! Surface pressure (pascals)
integer,  intent(in)  :: num_levs
real(r8), intent(out) :: pmid(ens_size, num_levs)   ! Pressure at model levels
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer k             ! Longitude, level indices
!-----------------------------------------------------------------------
!
! Set midpoint pressures and layer thicknesses
!
! coef
do k=1,num_levs
   pmid(:,k) = hyam%vals(k)*P0%vals(1) + hybm%vals(k)*p_surf
end do

return

end subroutine plevs_cam_distrib

!===============================================================================
   subroutine model_heights_distrib(state_ens_handle, ens_size, p_surf, num_levs, base_obs_loc, model_h_distrib, istatus)

!===============================================================================
! This routine calculates geometrical height (m) at mid-layers of the CAM model
!
! was Hui's dcz2ccm1
!    has globally defined inputs:
!          hyam(num_levs),hybm(num_levs),hyai(num_levs),hybi(num_levs) -
!          hybrid vertical coefficients, top to bottom. (P = P0*hyam + ps*hybm)
!             P0 - Hybrid base pressure (pascals)
! Kevin Raeder converted to single column version 4/28/2006
!              removed longitude dimension entirely and extra arrays 10/2006
!   5/31/2013; Rewritten to adapt to convert_vert handling obs TYPEs,
!              not obs KINDS, and to handle lonlat and cubed sphere
!              grids/interpolations.
!HK I think you are doing too much work with interp_cubed sphere - isn't it the same horizonal location
! each time?

! type(model_type), intent(in) :: Var
type(ensemble_type), intent(in) :: state_ens_handle
integer,             intent(in) :: ens_size
real(r8),            intent(in) :: p_surf(ens_size)
integer,             intent(in) :: num_levs
type(location_type), intent(in) :: base_obs_loc

! OUTPUT: geometrical height at midlayer (m)  hui liu /03/29/2004 added.
! model_h is already allocated, and in global storage, but must be listed here. !HK Why?
real(r8), dimension(:, :), intent(out) :: model_h_distrib ! HK This is the name of a global variable
integer,                intent(out) :: istatus(:)

! local variables; ps must be dimensioned as an array because dcz2 has it that way
real(r8), dimension(ens_size, num_levs) :: phi, tv, q, t, pterm
real(r8) :: pmln(num_levs+1), hyba(num_levs+1,2), hybb(num_levs+1,2)
real(r8) :: phi_surf(ens_size), ht_tmp, rd, rv, rr_factor
real(r8) :: lon_lat_lev(3)

integer :: k, i, e, vstatus(ens_size), track_vstatus(ens_size)

istatus = 0
vstatus = 0

! CS Should these come from common_mod?
! That might be inconsistent with how levels, etc were defined in CAM originally.
! DART's values are 287.0_r8 and 461.6_r8.
rd = 287.05_r8
rv = 461.51_r8
rr_factor = (rv/rd) - 1.0_r8

DO k = 1,num_levs ! HK why is this a loop?
   model_h_distrib(:, k) = 0.0_r8
   phi(:, k)     = 0.0_r8
   pterm(:, k)   = 0.0_r8
END DO

! copy to temporary arrays

!    All arrays except hyba, hybb are oriented top to bottom
!    modified to be consistent with CAM3.0 where hyai, hyam are top to bottom
!    H Liu, 04/05/2004

! interface=1 extra level at model bottom
k = num_levs +1
hyba(1,1) = hyai%vals(k)
hybb(1,1) = hybi%vals(k)
!   hyam(26) = 0 -> hyba(2,2) = 0, so it would be safe to set hyba(1,2) = 0.
!   This element is referenced below, but not ultimately used.
hyba(1,2) = 0.0_r8
hybb(1,2) = 1.0_r8
! hybX go from bottom to top; b coeffs multiply sigma, and coord is pure sigma
!      at the bottom, so hybb = 1.0 there.

! mid-points=2; note that hyXm(num_levs + 1) is not defined (= MISSING_R8)
do k = 2,num_levs +1
   i = num_levs +2 - k
   hyba(k,1) = hyai%vals(i)
   hybb(k,1) = hybi%vals(i)
   hyba(k,2) = hyam%vals(i)
   hybb(k,2) = hybm%vals(i)
! orig
! do k = 1,num_levs
!   hyba(k+1,2) = hyam(num_levs+1 - k)
!   hybb(k+1,2) = hybm(num_levs+1 - k)
end do

! Calculate phi_surf and tv for this column, for use by dcz2.
if (l_rectang) then
   ! ens_mean can be passed when interpolating phis because QTY_SURFACE_ELEVATION
   ! is handled by explicit references to phis, not ens_mean.
   call interp_lonlat_distrib(state_ens_handle, base_obs_loc, QTY_SURFACE_ELEVATION, phi_surf, vstatus)
!  call interp_lonlat(ens_mean, base_obs_loc, QTY_SURFACE_ELEVATION, phi_surf, vstatus)
   track_vstatus = vstatus

   do k = 1, num_levs
      call interp_lonlat_distrib(state_ens_handle, base_obs_loc, QTY_TEMPERATURE,       t(:, k), vstatus)
!     call interp_lonlat(ens_mean, base_obs_loc, QTY_TEMPERATURE,       t(k), vstatus)
      do e = 1, ens_size
         if (vstatus(e) > 0) then
            track_vstatus(e) = k
         endif
      enddo

      call interp_lonlat_distrib(state_ens_handle, base_obs_loc, QTY_SPECIFIC_HUMIDITY, q(:, k), vstatus)
!     call interp_lonlat(ens_mean, base_obs_loc, QTY_SPECIFIC_HUMIDITY, q(k), vstatus)
      do e = 1, ens_size
         if (vstatus(e) > 0) then
            track_vstatus(e) = num_levs+k  !HK I think you are going to have to deal with this - it could be overwritten
         endif
      enddo
      tv(:, k) = t(:, k)*(1.0_r8 + rr_factor*q(:, k))
   end do
else
   call interp_cubed_sphere_distrib(state_ens_handle, base_obs_loc, QTY_SURFACE_ELEVATION, phi_surf, vstatus)
!  call interp_cubed_sphere(ens_mean, base_obs_loc, QTY_SURFACE_ELEVATION, phi_surf, vstatus)
      track_vstatus = vstatus
   ! ? Use vstatus to tell interp_cubed_sphere to use the same location as last call?
   do k = 1, num_levs
      call interp_cubed_sphere_distrib(state_ens_handle, base_obs_loc, QTY_TEMPERATURE,       t(:, k), vstatus)
!     call interp_cubed_sphere(ens_mean, base_obs_loc, QTY_TEMPERATURE,       t(k), vstatus)
      do e = 1, ens_size
         if (vstatus(e) > 0) then
            track_vstatus(e) = k
         endif
      enddo
      call interp_cubed_sphere_distrib(state_ens_handle, base_obs_loc, QTY_SPECIFIC_HUMIDITY, q(:, k), vstatus)
!     call interp_cubed_sphere(ens_mean, base_obs_loc, QTY_SPECIFIC_HUMIDITY, q(k), vstatus)
      do e = 1, ens_size
         if (vstatus(e) > 0) then
            track_vstatus(e) = num_levs+k
         endif
      enddo
      tv(:, k) = t(:, k)*(1.0_r8 + rr_factor*q(:, k))
   end do
end if

do e = 1, ens_size
   if (track_vstatus(e) == 0 ) then ! only do this calculation for members that did not fail
      call dcz2(p_surf(e), phi_surf(e), tv, P0%vals(1) ,hyba, hybb, num_levs, pmln, pterm, phi)
   endif
enddo

istatus = track_vstatus

! used; hybb, hyba, hprb
! calced in dcz2;  pmln, pterm , zslice

! Need the latitude of the ob location.
lon_lat_lev = get_location(base_obs_loc)
do k = 1,num_levs
   ht_tmp = phi(e, k) * 0.001_r8        ! convert to km for following call only
   model_h_distrib(e, k) = gph2gmh (ht_tmp, lon_lat_lev(2)) * 1000.0_r8           ! convert back to m
end do

end subroutine  model_heights_distrib


   subroutine get_interp_prof (prof, vec, num_levs, lon_index, lat_index, &
                               stagr_lon, stagr_lat, kind_cam, vstatus)
!=====================================================================

real(r8), intent(out) :: prof(num_levs)
integer,  intent(out) :: vstatus
! type(model_type), intent(in) :: state
real(r8), intent(in)  :: vec(:)
integer,  intent(in)  :: kind_cam, num_levs, lon_index, lat_index
logical,  intent(in)  :: stagr_lon, stagr_lat

real(r8)  :: var(num_levs,0:1,0:1), weight
integer   :: k, lon_index_local, vec_ind
!integer :: fld_index
!character :: cfld

! So far only called from model_heights to get T and Q profiles for Tv calculation.
! If the point at which we need T and Q is a staggered US or VS point,
! then T and Q must be interpolated to that point.
! lon_index and lat_index are the indices of the point at which we need the profiles,
!   not necessarily one of the points where the profiles exist.
! stagr_xx tells where to look for profiles and what interpolating to do.
! staggered longitudes have indices the same as the next eastward  (+) unstaggered longitudes
! Staggered latitudes  have indices the same as the next southward (-) unstaggered lats.
! The indices called for here look weird, but think of it this way;
!    a point which is staggered from the usual grid has the same indices as the unstaggered point
!    to it's southeast.  The northwest corner of the box is then one higher latitude,
!    and one less (westward) longitude.
! Pole points should be handled in the calling routine by passing the correct stagr_xx,
! so that this program can count on having values for all the lat and lon indices referenced.

vstatus = 0
var = 0._r8

! Find the profile with the same lat and lon index

!  For now this only handles interpolations of 3D fields.
! fld_index = find_name(cflds(kind_cam), state_names_3d)
! if (fld_index == 0)  fld_index = find_name(cflds(dart_to_cam_kinds(kind_dart)), state_names_2d)
! if (fld_index == 0)  fld_index = find_name(cflds(dart_to_cam_kinds(kind_dart)), state_names_1d)

! Always do the nearest point
!   var(1:num_levs,0,0) = state%vars_3d(1:num_levs ,lon_index,lat_index, fld_index)
! index_from_grid wants kind within whole state vector, not 3d fields
!     vec_ind = index_from_grid(1,lon_index   ,lat_index   , fld_index)
vec_ind = index_from_grid(1,lon_index   ,lat_index   , kind_cam)
var(1:num_levs,0,0) = vec(vec_ind:vec_ind +num_levs -1)
weight = 1.0_r8

if (stagr_lat) then
   ! Find the profile to the north
   ! var(1:num_levs,0,1) = state%vars_3d(1:num_levs ,lon_index  ,lat_index +1, fld_index)
   vec_ind = index_from_grid(1,lon_index   ,lat_index +1, kind_cam)
   var(1:num_levs,0,1) = vec(vec_ind:vec_ind +num_levs -1)
   ! This weighting is not strictly correct for Gaussian latitudes,
   ! but is correct for CAM-FV, which has equally spaced latitudes.
   weight = weight / 2.0_r8
end if

if (stagr_lon) then
   ! Find the profile to the west
   if (lon_index == 1) then
      ! handle the Greenwich meridion points.
      ! Yes, the index should equal the dimension size in this case.
      ! The dimension size is for the unstaggered grid, since that's where T and q are.
      lon_index_local = dim_sizes(find_name('lon     ',dim_names))
      ! var(1:num_levs,1,0) = state%vars_3d(1:num_levs ,lon_index_local,lat_index, fld_index)
   else
      lon_index_local = lon_index -1
      ! var(1:num_levs,1,0) = state%vars_3d(1:num_levs ,lon_index_local,lat_index, fld_index)
   end if
   vec_ind = index_from_grid(1,lon_index_local ,lat_index, kind_cam)
   var(1:num_levs,1,0) = vec(vec_ind:vec_ind +num_levs -1)
   weight = weight / 2.0_r8
end if

if (stagr_lon .and. stagr_lat) then
   ! Find the profile to the northwest
   ! var(1:num_levs,1,1) = state%vars_3d(1:num_levs ,lon_index_local,lat_index +1, fld_index)
   vec_ind = index_from_grid(1,lon_index_local ,lat_index +1, kind_cam)
   var(1:num_levs,1,1) = vec(vec_ind:vec_ind +num_levs -1)
end if

do k = 1,num_levs
   prof(k) = weight * (var(k,0,0) + var(k,1,0) + var(k,0,1) + var(k,1,1))
end do

end subroutine get_interp_prof

!
! height
!=====================================================================
   subroutine dcz2(p_surf,phis0,tv,hprb,hyba,hybb,kmax,pmln, pterm,z2)
!=====================================================================
!       Purpose:
!         To compute geopotential height for a CCM2 hybrid coordinate
!         vertical slice.  Since the vertical integration matrix is a
!         function of latitude and longitude, it is not explicitly
!         computed as for sigma coordinates.  The integration algorithm
!         is derived from Boville's mods in the ibm file hybrid 1mods
!         (6/17/88).  All vertical slice arrays are oriented top to
!         bottom as in CCM2.  This field is on full model levels (aka
!         "midpoints") not half levels.
!
!       Equation references are to "Hybrid Coordinates for CCM1"
!
!----------------------------Code History-------------------------------
!       Original: Jun 25 1991  L Buja
!       Upgrades: Apr 21 1993  L Buja
!                     Truesdale reports difference from CCM2 Z2.
!                     - Make pterm a scratch array (was automatic)
!                     - Make g0=9.80616 (was 9.81).  This appears to
!                        affect only the lowest layer.
!       Change  : Feb    1999: D Shea
!                     NCL changes + error
!                     - The "Invert vertical loop" has to have
!                        the argument checked for >0.0
!-----------------------------------------------------------------------
!-----------------------------Parameters--------------------------------
      real(r8) :: r,g0,rbyg
      parameter (r=287.04_r8,g0=9.80616_r8,rbyg=r/g0)
!-----------------------------Arguments---------------------------------
! Input
!
! Number of vertical levels
      integer kmax

! Surface pressure           (pascals)
      real(r8) :: p_surf

! Surface geopotential
      real(r8) :: phis0

! Virtual temperature, top to bottom
      real(r8) TV(KMAX)

! Hybrid base pressure       (pascals)
      real(r8) :: HPRB

! Hybrid coord coeffs for base pressure
      real(r8) :: HYBA(KMAX+1,2)

!       All arrays except hyba, hybb are oriented top to bottom
!  ground to top, first subscript:

!  = 1 for layer interfaces
!  = 2 for layer midpoints
!  Lowest level is ground for both layer locations

! Hybrid coord coeffs for surf pressure (in same format as hyba)
      real(r8) ::  hybb(kmax+1,2)

! vertical slice scratch space used to
!   hold logs of midpoint pressures
      real(r8) ::  pmln(kmax+1)

! Note: These scratch vertical slices are used to improve computaional efficiency

! Vertical scratch space.
      real(r8)::   pterm(kmax)
! temporary
      real(r8)::   arg
!
! Output ---------------------------------
!
! Geopotential height, top to bottom
      real(r8)::   Z2(KMAX)
!
!--------------------------Local variables------------------------------
! indexes
      integer i,k,l,num
!-----------------------------------------------------------------------
!
      DATA NUM/0/

      NUM = NUM + 1
!

!       Compute intermediate quantities using scratch space

!       Invert vertical loop
!       Compute top only if top interface pressure is nonzero.
!       Implemented by setting loop limit klim
!
!       hyba, hybb are bottom to top, starting at ground.
!       pmln(i,k) is the mid-point pressure of layer k.
!       SHEA MODIFICATION

      DO K = KMAX + 1,1,-1
         i = KMAX-K+2
         ARG = HPRB*HYBA(i,2) + P_surf *HYBB(i,2)
         IF (ARG.GT.0.0_r8) THEN
             PMLN(K) = DLOG(ARG)
         ELSE
             PMLN(K) = 0.0_r8
         END IF
      END DO
!
!       Initialize Z2 to sum of ground height and thickness of
!        top half-layer  (i.e. (phi)sfc in equation 1.14)
!       (Z2(i,1)=top  ->  Z2(i,kmax)=bottom
!       Eq 3.a.109.2  where l=K,k<K  h(k,l) = 1/2 * ln [  p(k+1) / p(k) ]

      DO K = 2,KMAX - 1
         pterm(k) = rbyg*tv(k)*0.5_r8* (pmln(k+1)-pmln(k-1))
      END DO

!
      DO K = 1,KMAX - 1
         z2(k) = phis0/g0 + rbyg*tv(k)*0.5_r8* (PMLN(K+1)-PMLN(K))
      END DO

!       Eq 3.a.109.5  where l=K,k=K  h(k,l) = ln [ pi / (p(k)) ]

      K = KMAX
!
      z2(K) = phis0/g0 + rbyg*tv(k)* (dlog(p_surf*hybb(1,1))-pmln(k))

!       Eq 3.a.109.4  where l=K,k<K  h(k,l) = 1/2*ln[pi*pi/(p(k-1)*p(k))

!
      do k = 1,kmax - 1
          l = kmax
          z2(k) = z2(k) + rbyg*tv(l)* (dlog(p_surf*hybb(1,1))-0.5_r8* &
                                       (pmln(l-1)+pmln(l)))
      end do

!       Add thickness of the remaining full layers
!        (i.e., integrate from ground to highest layer interface)

!       Eqs 1.14 & 3.a.109.3 where l>K, k<K
!                                h(k,l) = 1/2 * ln [ p(l+1)/p(l-1) ]

!
      DO K = 1,KMAX - 2
          DO L = K + 1,KMAX - 1
             Z2(K) = Z2(K) + pterm(L)
          END DO
      END DO

      RETURN
 end subroutine dcz2

! height
!/**----------------------------------------------------------------------
! name       gph2gmh
!  Convert a list of geopotential altitudes to mean sea level altitude.
!
!  input:    h   -- geopotential altitude (in km)
!            lat -- latitude  of profile in degrees.
!  output:   z   -- MSL altitude, in km.
! -----------------------------------------------------------------------*/
   function gph2gmh (h, lat)

      real(r8) ::  h, lat, gph2gmh
      real(r8) ::  be, ae, pi, G, g0, r0, latr

      be = 6356.7516_r8             ! min earth radius, km
      ae = 6378.1363_r8             ! max earth radius, km

      pi = 3.14159265358979_r8
      latr = lat * (pi/180.0_r8)           ! in radians

  ! These are the equations for g0 and r0 in Chris Rocken's paper.
  ! I am not using them because they imply a standard ellipsoid which
  ! is quite different from our standard ellipsoid. D. Hunt 10/28/96
  !G = 0.0098
  !g0 = 0.001 * 9.780356 * (1+0.0052885 * (sin(lat))**2 - 5.9e-6 * (sin(2*lat))**2)
  !r0 = (2*g0)/(3.085462e-6 + 2.27e-9 * cos(2*lat) - 2e-12*cos(4*lat))

      G = 0.00980665_r8          ! WMO reference g value, km/s**2, at 45.542N(S)

      g0 = 0.0_r8
      call gravity (latr, 0.0_r8, g0)
! liu    g0 = g0 * 0.00001_r8             ! convert to km/s**2

! compute local earth's radius using ellipse equation
!
      r0 = dsqrt ( ae**2 * dcos(latr)**2 + be**2 * dsin(latr)**2)

!     if (h.eq.-999.0_r8) then
!        z = -999.0_r8
!     else
! Compute altitude above sea level

         gph2gmh = (r0 * h) / (((g0*r0)/G) - h)
!     end if

end function gph2gmh

!=============================================================
   subroutine gravity(xlat,alt,galt)
!=============================================================
! This subroutine computes the Earth's gravity at any altitude
! and latitude.  The model assumes the Earth is an oblate
! spheriod rotating at a the Earth's spin rate.  The model
! was taken from "Geophysical Geodesy, Kurt Lambeck, 1988".
!
!  input:    xlat, latitude in radians
!            alt,  altitude above the reference ellipsiod, km
!  output:   galt, gravity at the given lat and alt, cm/sec
!
! Compute acceleration due to the Earth's gravity at any latitude/altitude
! author     Bill Schreiner   5/95
! ------------------------------------------------------------------------

  real(r8) :: xmu, ae, f, w, xm, f2, f4, ge, g, galt, xlat,alt
!
      xmu = 398600.4415_r8       ! km^3/s^2
      ae = 6378.1363_r8          ! km
      f = 1.0_r8/298.2564_r8
      w = 7.292115d-05          ! rad/s
      xm = 0.003468_r8           !
!     f2 = -f + 5.0* 0.50*xm - 17.0/14.0*f*xm + 15.0/4.0*xm**2
!     f4 = -f**2* 0.50 + 5.0* 0.50*f*xm
      f2 = 5.3481622134089D-03
      f4 = 2.3448248012911D-05
!
! compute gravity at the equator, km/s2
!
      ge = xmu/ae**2/(1.0_r8 - f + 1.5_r8*xm - 15.0/14.0*xm*f)
!
! compute gravity at any latitude, km/s2
!
      g = ge*(1.0_r8 + f2*(dsin(xlat))**2 - 1.0/4.0*f4*(dsin(2.0*xlat))**2)
!
! compute gravity at any latitude and at any height, km/s2
!
      galt = g - 2.0*ge*alt/ae*(1.0 + f + xm + (-3.0*f + 5.0* 0.50*xm)*  &
                             (dsin(xlat))**2) + 3.0*ge*alt**2/ae**2
!
!liu     galt = galt*1.0d5              ! convert from km/s2 to cm/s2
!
end subroutine gravity



   subroutine init_model_instance(var)
!=======================================================================
! subroutine init_model_instance(var)
!
! Initializes an instance of a cam model state variable

type(model_type), intent(out) :: var

! Initialize the storage space and return

! The temporary arrays into which fields are read are dimensioned by the largest values of
! the sizes of the dimensions listed in s_dim_RANKd

allocate(  var%vars_0d (                                              state_num_0d))
allocate(  var%vars_1d (s_dim_max(1,1),                               state_num_1d))
allocate(  var%vars_2d (s_dim_max(1,2),s_dim_max(2,2),                state_num_2d))
allocate(  var%vars_3d (s_dim_max(1,3),s_dim_max(2,3),s_dim_max(3,3), state_num_3d))

end subroutine init_model_instance



   subroutine end_model_instance(var)
!=======================================================================
! subroutine end_model_instance(var)
!
! Ends an instance of a cam model state variable


type(model_type), intent(inout) :: var

deallocate(var%vars_0d, var%vars_1d, var%vars_2d, var%vars_3d)

end subroutine end_model_instance


! End of utility routines

!#######################################################################

! Stubs not used by cam/model_mod (this is not all of them)

   subroutine adv_1step(st_vec, Time)
!=======================================================================
! subroutine adv_1step(st_vec, Time)
!

real(r8), intent(inout) :: st_vec(:)

! Time is needed for more general models like this; need to add in to
! low-order models
type(time_type), intent(in) :: Time

! This is a no-op for CAM; only asynch integration
! Can be used to test the assim capabilities with a null advance

! make it an error by default; comment these calls out to actually
! test assimilations with null advance.

call error_handler(E_ERR,'adv_1step', &
                  'CAM model cannot be called as a subroutine; async cannot = 0', &
                  source, revision, revdate)

end subroutine adv_1step




   subroutine end_model()
!=======================================================================
! subroutine end_model()
!
! At some point, this stub should coordinate with atmosphere_end but
! that requires an instance variable.

! Deallocate other variables?

end subroutine end_model



   subroutine init_time(time)
!=======================================================================
! subroutine init_time(time)
!
! For now returns value of Time_init which is set in initialization routines.

type(time_type), intent(out) :: time

! Where should initial time come from here?
! WARNING: CURRENTLY SET TO 0
time = set_time(0, 0)

end subroutine init_time


!--------------------------------------------------------------------------
!> This is supposed to replace set_ps_arrays_distrib
function get_surface_pressure_state(state_ens_handle, ens_size, lon_ind, lat_ind)

integer,             intent(in)  :: ens_size
type(ensemble_type), intent(in)  :: state_ens_handle
integer,             intent(in)  :: lon_ind
integer,             intent(in)  :: lat_ind

real(r8) :: get_surface_pressure_state(ens_size)
!> ifld ... pressure field index
integer  :: ifld
!> ind ... index into state vector
integer  :: ind

ifld = find_name('PS      ',cflds)

! find index into state
ind = index_from_grid(1, lon_ind, lat_ind, ifld)

! get correct piece of state
call get_state(get_surface_pressure_state, ind, state_ens_handle)

end function get_surface_pressure_state

!-------------------------------------------------------------------------
!> Non-array version of get_surface_pressure
function get_surface_pressure_mean(state_ens_handle, lon_ind, lat_ind)

type(ensemble_type), intent(in) :: state_ens_handle
integer,             intent(in) :: lon_ind
integer,             intent(in) :: lat_ind

real(r8) :: get_surface_pressure_mean
integer  :: ifld
integer  :: ind

ifld = find_name('PS      ',cflds)

! find index into state
ind = index_from_grid(1, lon_ind, lat_ind, ifld)

! get correct piece of state
call get_state(get_surface_pressure_mean, ind, state_ens_handle)

end function get_surface_pressure_mean


!--------------------------------------------------------------------
!> This returns the vertical coordinate of an observation in the
!> requested vertical localization coordinate. 
!> Aim: to have only the process who owns the observation do this calulation, 
!> rather than all processeses doing the same calculation in get_close_obs_distrib
subroutine convert_base_obs_location(obs_loc, state_ens_handle, vert_coord, istatus)

type(location_type), intent(inout) :: obs_loc
type(ensemble_type),    intent(in) :: state_ens_handle
real(r8),              intent(out) :: vert_coord
integer,               intent(out) :: istatus

real(r8), dimension(3) :: base_array
!> @todo Should check for identity obs
integer                :: base_obs_kind
integer                :: base_which ! vertical coorardiate
integer                :: istatus_v

base_obs_kind = 1 ! dummy for now, should check for identity obs

base_which = nint(query_location(obs_loc))

!if (base_which /= wrf%dom(1)%localization_coord) then *** WRF specific! ***
   !call convert_vert_distrib(state_ens_handle, obs_loc, base_obs_kind, istatus_v)
!endif
if (base_which /= VERTISPRESSURE) then
   call error_handler(E_ERR, 'convert_base_obs_location ', 'broken here fix this')
endif

istatus = istatus_v

base_array = get_location(obs_loc)
vert_coord = base_array(3)

!> @todot set location so you don't redo this calculation in get_close_obs NOPE  they are two different structures
!obs_loc = set_location(base_array(1), base_array(2), base_array(3), wrf%dom(1)%localization_coord )

end subroutine convert_base_obs_location


!#######################################################################
! end of cam model_mod
!#######################################################################

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source: /home/thoar/CVS.REPOS/DART/models/cam/model_mod.f90,v $
! $Revision$
! $Date$
! $Author$
! $Name:  $
 
!--------------------------------------------------------------------
! DISTRIBUTION;
! > instructions for use
! > new filter_ics for the cases I've distributed
! > init2ud script for them to create their own.
!      list of files that script/model_mod needs; 
!          trans_pv_sv_time0  (the executable)
!          input.nml:perfect_model_nml (with the date+time to put in the filter_ics files)
!          topog_file.nc  (for static_init_model to get PHIS)
!          caminput_#.nc  (to convert into filter_ic.####)

!--------------------------------------------------------------------
! DISTANCE CONVERSION;
! 
! model_interpolate will continue to use state passed to it;
!    recalc p_col and model_h columns as needed.
!    no need to convert to a standard vert coord; no distance calc involved.

! model_mod;get_close_obs will be getting an ob, with its location, and its horizontal distances 
!    to an array of other locations (and the locations).
!       These locations were picked out based on the efficient search/box algorithm.
!    It converts vertical coordinates as needed, 
!    It tests for being above the highest_obs_pressure_mb threshold, and increases the
!       vertical distance based on height above highest_.
!    It calls location_mod/threed_sphere:get_close_obs, 
!       to which it sends this (converted) array of locations.
!    It gets back the new total distances and arrays of those locations that are "close"
!       to the base ob.
! get_close_obs will use the ensemble average passed through new interface; ens_mean_for_model
!    Convert the obs and/or state vertical location(s) to a standard (pressure) vertical location
!    3D model_h would be useful here; calc once and use over and over.
!       reinstall height/lon slice code for model_heights to facilitate that
!    3D pressures also useful here; 
!       reinstall height/lon slice code for plevs_cam to facilitate that
!    throw away ens_mean after it's been used (or don't worry about it for now).
! 
!
!--------------------------------------------------------------------
! 
! I assume that staggered grids (US and VS) are staggered only 1 direction (each), 
! so that surface pressure interpolations to get staggered ps use only 2 A-grid ps values.
! The interpolations for columns of heights are more general, but will do a 2 point interp
!     if the staggering is only in one direction.
! 
! Describe: 
!    - The coordinate orders and translations; CAM initial file, model_mod, and DART _Diag.nc. 
!      Motivations
!    - ! I think I need to have 2 sets of arrays for dimensions and dimids; one describing the 
!        caminput file ;
!          (f_...) and one for the state (s_...) (storage in this module).  
!             Call them f_dim_Nd , f_dimid_Nd
!                       s_dim_Nd , s_dimid_Nd
! 



! 'merge' marks changes for merging fv and eul core model_mods,
!         although there are so many that most are not marked.
!    x  > I'll reorganize state vector into columns of single variables.
!    x  > There can be different numbers of columns for different variables.
!    x     (Do I need to keep corresponding columns of different variables close to each other?)
!    x     I'm not, but this would be done in prog_var_to_vect and vect_to_prog_var.
!    x  > init_model_instance allocates space for a CAM format state vector; 
!    x    I would need to use the max size for each dimension among the variables with the same 
!    x    dimension
!    x  > Check every place where init_model_instance is used
!    x  > converted read_cam_coord (check),  
!   0   > Convert read_cam_horiz too?   convert phis to a grid_2d_type?
!    x  > convert prog_xxx_yyy_
!    x  > get_state_meta_data needs rewriting; new state vector org, and coord dimension names
!  ? x  > model_get_close_states; the get_close functionality will need rewriting if it will be
!    x    used under MPI
!    x  > Put in MPI adaptations (no(?) model_get_close_states, but add model_convert_vert
!    x    and model_alter_dist)
!    x  > FV, other grids; must do interpolations of ps to staggered grid in order to do vertical
!    x    interpolations on staggered grids.
!    x  > model_interpolate and routines it calls not done yet (?!)
!    x  > look at all FIX sites for things to change for this merge - pert_model_state
!    x  > cycle through all the grid arrays; find where they're used and change their references
!         to structure components.
!    x  > num_lxxs are used all over the place; search and replace (can't do a universal 
!         substitution because num_lats has several values now, num_lats num_slats.
!    x  > read_cam_scalar no longer needed?
!    x  > look at all slon usage and be sure that I've accounted for wierd CAM staggering
!    x  > Add interpolations of t and q for calculating heights at staggered grid points
!    x  > Nancy or F90 book:
!    x    grid_1d_type ! Should vals be a pointer?  see type model_type
!    x    Should the contents be private?  components are used within this module.
!            does private restrict reference to within this module, or this type?
!    x    coord_index pointers in to coordinate structures
!    x  > ? standard vertical coordinate?  Ask Jeff.
!           For get_close_obs, but not for model_interpolate
!    x  > Derive highest_obs_level and highest_obs_height from highest_obs_pressure_mb
!            in subroutines that need them.
!    x  > Should I pass state vector(which is what comes from filter_assim) or cam_state
!           (easier to use, but slower?) to model_heights? 
!    x  > See 'premature' below for something else to fix
!    x  > clean out old code from alter_dist, after I'm sure that it works
!    x  > clean out maxval section, when I'm convinced it works.
!    x  > compile with options to tell me about unused variables, etc.
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

! EFFICIENCY
!       + index_from_grid (and others?) could be more efficient by calculating and 
!         globally storing the beginning index of each cfld and/or the size of each cfld.
!       + searchs through dimensions are inefficient
!     x + simplify plevs_cam; eliminate a loop and argument(s?)
!       + add dimid to grid_1d_type; to reduce need to call find_names?
!       + pointers may slow down computations by up to 75%.  Dereference them by passing
!         the pointer array to a subroutine, which will then handle the data directly
!         instead of through a pointer.
!       + global storage of height fields?  but need them on staggered grids (only sometimes)

! ISSUE; In P[oste]rior_Diag.nc ensemble members are written out *between* the field mean/spread
!        pair and the inflation mean/sd pair.  Would it make more sense to put members after
!        both pairs?  Easy to do?

! ISSUE?; model_interpolate assumes that obs with a vertical location have 2 horizontal locations 
!          too.  The state vector may have fields for which this isn't true, but no obs we've seen
!          so far violate this assumption.  It would have to be a synthetic/perfect_model obs, 
!          like some sort of average or parameter value.  

! ISSUE; The KIND_ list from obs_def_mod must be updated when new fields are added to state vector.
!        This could be done by the preprocessor when it inserts the code bits corresponding to the
!        lists of observation types, but it currently (10/06) doesn not.  Document accordingly.

! ISSUE; In convert_vert, if a 2D field has dimensions (lev, lat) then how is p_surf defined?
!        I've set it to P0, but is this correct or meaningful?

! ISSUE; is there a machine precision defined somewhere, which I could use for testing differences
!        betwen real numbers?  search for E-12 to see instances.
!   FIX; Fortran intrinsic function; not implemented yet

! ISSUE: In model_get_close_states the P on levels is calculated for each column within
!        radius of the ob location.  So P is recalculated for the same columns many times
!        for the many obs.  Same for obs on heights.
!        Model_get_close_states will disappear with MPI, but the functionality is still needed
!        by new DART routine get_close_obs.  Functionality will be provided by new routine
!        convert_vert.
!   FIX; It's not strictly correct, but the calculations of pfull and heights on model levels
!        could be moved up front (static_init_model?) and done only once for all columns.
!        Jeff says do this for PHIS, but not for pfull; just make ps a globally accessible
!        variable, to remove the ps = 100000. kluge.

! ISSUE: Should get_val_xxxx return MISSING_VALUE instead of 0 for istat = 1 cases?

! ISSUE: The CCM code (and Hui's packaging) for geopotentials and heights  use different
!        values of the physical constants than DARTS.  In one case Shea changed g from
!        9.81 to 9.80616, to get agreement with CCM(?...), so it may be important.
!        Also, matching with Hui's tests may require using his values;  change to DART
!        after verifying?

! "Pobs" marks changes for providing expected obs of P
!        break from past philosophy; P is not a native CAM variable (but is already calced here)

! NOVERT marks modifications for fields with no vertical location,
! i.e. GWD parameters.


! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$
! 
!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!              Translate to/from state_vector and caminput.nc,
!              Initialize model,
!              Write out CAM fields to Prior and Posterior_Diag.nc,
!              Generate expected obs from model state (model_interpolate)
!              Find state variables (or obs) that are close to a given base observation.
!                  (get_close_obs)
!              Provide randomly perturbed fields for initial ensemble in filter.(pert_model_state)

!
! method: trans_pv_sv reads CAM 'initial' file (netCDF format).
!         Reform fields into a state vector.
!         DART assimilates obs, modifying the values of the state vector.
!         This requires get_close_obs to find state variables to be modified by each ob.
!         Reform state vector back into CAM fields.
!         Replace those fields on the CAM initial file with the new values,
!         preserving all other information on the file.
!         Also read hybrid coordinate coefficients from CAM input file (for plevs_cam)
!
! author: Kevin Raeder 2/14/03  and 8/1/03
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
!
! modified: Tim Hoar 02 Sep 03 
!         nc_write_model_atts, nc_write_model_vars now write out "prognostic"
!         files instead of a nondescript state variable vector glom
!
! augmented; Kevin Raeder 7/1/04 for CAM3.0 and to use namelist input for
!         lists of fields to include in state vector.
!         'CAM3' marks changes
!         Later;?
!         Also; read field attributes from netcdf file and write them out 
!         from nc_write_model_atts, instead of hard-coded there.
!
! augmented; Kevin Raeder (code from Hui Liu and CAM) 4/2006 to add vertical interpolation
!         in height.  Also add checks of requested fields for interpolation.
!
! mostly rewritten: Kevin Raeder  9-10/2006 

!     OVERVIEW; 
!     
!     This module interfaces with the MPI version of DART.  
!
!        model_get_close_states has been replaced with new interfaces, including get_close_obs, 
!     which uses location_mod:get_close_obs.  This allows model_mod to use the 
!     "generic efficient search algorithm" in DART to find the obs/state variables that
!     are close to the given base obs, and to modify the distances according to any desired criteria
!     (obs_kinds, height above a threshold, etc).  
!
!     During the assimilation stage, only a piece of the state vector is available to each 
!     process, and each process calls parts of model_mod.  In order to handle the conversion 
!     of vertical coordinates of obs and/or state variables into a consistent coordinate, 
!     an entire state vector is needed, so the ensemble mean is passed to model_mod before 
!     the assimilation starts.  All locations are now converted to a standard coordinate 
!     (pressure), instead of always converting the state vertical location to that of the ob.
!     The highest_obs_level and ..._height_m variables are derived from highest_obs_pressure_mb
!     namelist variable.
!     
!     This module has been rewritten to handle both the eulerian and finite volume core versions 
!     of CAM (they have different grids), and hopefully the semi-lagrangian dynamics core, 
!     and even lay some groundwork for future dynamical cores (such as HOMME) which are column 
!     oriented, with irregular horizontal grids.
!     
!     The coordinate orders of fields stored in various forms have also been simplified.
!     For example; various vintages of CAM 3D fields may be read in with (lon, lat, lev) or 
!     (lon, lev, lat).  These are uniformly converted to (lev, lon, lat) for use in model_mod.
!     This latter form is different than pre MPI model_mods.  Then such fields are stored in
!     the state vector with the same coordinate order.  They are converted back to the modern
!     CAM coordinate order when written to caminput.nc files.
!     
!     The organization of the state vector has been changed from the B-grid style to a field 
!     and column oriented form; see prog_var_to_vector and vector_to_prog_var.  This allows 
!     each field to have different dimensions than all the others.  The fields are still grouped 
!     by the rank of the array used to store them in CAM form; 0d, 1d, 2d, and 3d.  In the 3d case 
!     all the values of a field for one column are stored continguously, to speed up vertical 
!     interpolations, and calculations of heights on model levels.  The dimensions of each field 
!     are stored in globally available arrays and there are functions and subroutines for accessing 
!     the needed information.
!     
!     Observations on pressure levels or heights require various pieces of the state vector.
!     The surface pressure for the former, and that plus the temperature and moisture profiles
!     for the heights.  These may be needed on the regular A-grid (thermodynamic variables) and
!     grids staggered relative to the A-grid.   Currently, PS for the A-grid and the 2 staggered
!     grids is stored for global access and pressures and heights on model leves are (re)calculated
!     as needed.  In the future it may be deemed worthwhile to store the
!     3d pressures and heights on the 3 grids, but for now that seemed like too much memory to 
!     make that worthwhile.
!
!     It also corrects a misuse of TYPE_s from pre MPI versions; those model_mod specific 
!     identifiers are no longer passed back to filter through get_state_meta_data.  Instead, 
!     DART KIND_ identifiers are used.  If a user wants to add more TYPE_s to the state vector, 
!     then more KIND_s will be needed in the obs_kind_mod and the 'use obs_kind_mod' statement.

!     The coordinates of CAM (lats, lons, etc.) and their dimensions  and attributes are 
!     now read into globally accessible data structures (see grid_1d_type).
!     
!     MODULE ORGANIZATION
!    
!         'use' statements
!          Global storage for describing cam model class
!          Namelist variables with default values follow
!          derived parameters
!          static_init_model section
!          module I/O to/from DART and files
!          model_interpolate section
!          vector-field translations
!          get_close_obs section
!          End of get_close_obs section
!          Utility routines; called by several main subroutines
!          Stubs not used by cam/model_mod (this is not all of them)


!==============================================================================================
!  USE statements

use        netcdf
use        typeSizes

use        types_mod, only : r8, MISSING_I, MISSING_R8
!          add after verification against Hui's tests;  gas_constant_v,gas_constant,ps0,PI,DEG2RAD

use time_manager_mod, only : time_type, set_time, print_time, set_calendar_type,                &
                             THIRTY_DAY_MONTHS, JULIAN, GREGORIAN, NOLEAP, NO_CALENDAR
use    utilities_mod, only : open_file, close_file, find_namelist_in_file, check_namelist_read, &
                             register_module, error_handler, file_exist, E_ERR, E_WARN, E_MSG,  &
                             logfileunit, do_output

!-------------------------------------------------------------------------
use     location_mod, only : location_type, get_location, set_location, query_location,            &
                             LocationDims, LocationName, LocationLName,                            &
                             vert_is_level, vert_is_pressure, vert_is_height, vert_is_surface,     &
                             VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, VERTISPRESSURE, VERTISHEIGHT,&
                             get_close_type, get_close_maxdist_init, get_close_obs_init, get_dist, &
                             loc_get_close_obs => get_close_obs

! get_close_maxdist_init, get_close_obs_init, can be modified here (i.e. to add vertical information
! to the initial distance calcs), but will need subroutine pointers like get_close_obs.

!-----------------------------------------------------------------------------
use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT,                    &
                             KIND_SURFACE_PRESSURE, KIND_TEMPERATURE, KIND_SPECIFIC_HUMIDITY, &
                             KIND_PRESSURE, KIND_CLOUD_LIQUID_WATER, KIND_CLOUD_ICE
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

use    random_nr_mod, only : init_ran1
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

! end of use statements
!==============================================================================================
!
! CAM global/module declarations

implicit none
private

! The last three lines are for interfaces to CAM specific programs; transXXX, etc.
public ::                                                            &
   static_init_model, get_model_size, get_model_time_step,           &
   pert_model_state, get_state_meta_data, model_interpolate,         &
   nc_write_model_atts, nc_write_model_vars,                         &
   init_conditions, init_time, adv_1step, end_model,                 &
   model_type, prog_var_to_vector, vector_to_prog_var,               &
   read_cam_init, read_cam_init_size,                                &
   init_model_instance, end_model_instance, write_cam_init,          &
   get_close_maxdist_init, get_close_obs_init, get_close_obs, ens_mean_for_model


!-----------------------------------------------------------------------
! CVS Generated file description for error handling, do not edit
character(len=128) :: version = "$Id$"
character(len=128) :: tag = "$Name:  $"
character(len=128) :: &
source   = "$Source: /home/thoar/CVS.REPOS/DART/models/cam/model_mod.f90,v $", &
revision = "$Revision$", &
revdate  = "$Date$"
!-----------------------------------------------------------------------

! merge/MPI; 
! DART form of ensemble mean, global storage for use by get_close_obs:convert_vert 
! Ensemble mean is used so that the same "state" will be used for the height calculations
! on all processors, for all ensemble members.
! Nancy; need allocatable here if it's an argument passed into a subroutine, with space
! allocated in another module?
real(r8), allocatable :: ens_mean(:)      

!----------------------------------------------------------------------
! Global storage for describing cam model class
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Definition of variable types
! Values will be defined in order_state_fields
! All fields will have entries in the TYPE_xD corresponding to their orders
!   in state_names_Xd.  The explicitly named TYPE_s are for convenience
integer :: TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q, TYPE_CLDICE, TYPE_CLDLIQ,     &
           TYPE_PHIS, TYPE_SGH, TYPE_PBLH, TYPE_TBOT, TYPE_TS, TYPE_TSOCN,        &
           TYPE_LCWAT, TYPE_QCWAT 
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
type grid_1d_type 
   private
   character (len=8)            :: label
   integer                      :: dim_id
   integer                      :: length
   real(r8), pointer            :: vals(:)
   integer                      :: num_atts
   character (len=128), pointer :: atts_names(:)
   character (len=128), pointer :: atts_vals(:)
end type grid_1d_type

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

! File where basic info about model configuration can be found
character(len = 128) :: model_config_file = 'caminput.nc', &
                        model_version     = '3.0',         &
                        topog_file        = 'topog_file.nc'


! Define location restrictions on which observations are assimilated
! (values are calculated anyway, but istatus is set to 2)
real(r8) :: max_obs_lat_degree        = 90.0_r8
real(r8) :: highest_obs_pressure_mb   = 30.0_r8
real(r8) :: highest_state_pressure_mb = 30.0_r8
! These are not namelist variables, but are related, and calculated from highest_obs_pressure_mb
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

! make these allocatable?  Where would I define the default values?
character (len=8),dimension(1000) :: state_names_0d = (/('        ',iii=1,1000)/)
character (len=8),dimension(1000) :: state_names_1d = (/('        ',iii=1,1000)/)
character (len=8),dimension(1000) :: state_names_2d = (/'PS      ',('        ',iii=1,999)/)

character (len=8),dimension(1000) :: state_names_3d =  &
          (/'T       ','U       ','V       ', &
            'Q       ', ('        ',iii=1,996)/)

! NOVERT  need (a) namelist parameter(s) to define which_vert for each 2D (xD?) field
!         There's a danger of having a mismatch with the state_names_Xd; should this
!         definition be part of state_names_Xd, which is parsed into a name and which_vert
!         after being read?  Not for now.

integer , dimension(1000) :: which_vert_1d = (/(-2,iii=1,1000)/)
integer , dimension(1000) :: which_vert_2d = (/(-1,iii=1,1000)/)
integer , dimension(1000) :: which_vert_3d = (/( 1,iii=1,1000)/)


! Is there a way to exclude state_nums from namelist and have those filled in 
! the  subroutine which sorts state_names?
! Yes, use two namelists model_nml_1 and model_nml_2 at future date

! list of fields which trans_pv_sv_pert0 needs to perturb because they're
! constant valued model parameters and show no spread when start_from_restart = .true.
character (len=8),dimension(1000) :: state_names_pert = (/('        ',iii=1,1000)/)
real(r8), dimension(1000) :: state_names_sd = (/(0.,iii=1,1000)/)


! Specify shortest time step that the model will support
! This is limited below by CAMs fixed time step but is also impacted
! by numerical stability concerns for repeated restarting in leapfrog.
integer :: Time_step_seconds = 21600, Time_step_days = 0

namelist /model_nml/ output_state_vector , model_version , model_config_file & 
                     ,state_num_0d   ,state_num_1d   ,state_num_2d   ,state_num_3d &
                     ,state_names_0d ,state_names_1d ,state_names_2d ,state_names_3d &
                     ,                which_vert_1d  ,which_vert_2d  ,which_vert_3d &
                     ,state_names_pert ,state_names_sd &
                     ,highest_obs_pressure_mb , highest_state_pressure_mb  &
                     ,max_obs_lat_degree ,Time_step_seconds ,Time_step_days 
                     

!---- end of namelist (found in file input.nml) ----
!----------------------------------------------------------------------
! derived parameters

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type(time_type) :: Time_step_atmos

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Random sequence and init for pert_model_state
logical                 :: first_pert_call = .true.
type(random_seq_type)   :: random_seq
! Variable for keeping track of which ensemble member is to be perturbed
! by pert_model_state, which is called by filter for each ensemble member
! for a cold start.
integer                 :: ens_member
data ens_member /0/
logical                 :: do_out


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
! Need more arrays if any future fields are doubly staggered.
! ? merge; should this be a grid_2d_type, specified above?
logical               :: alloc_ps=.true.    ! Flag whether to alloc space for ps[_stagr] 
real(r8), allocatable :: ps(:, :)           ! surface pressure used to calc P and height profiles.
real(r8), allocatable :: ps_stagr_lon(:, :) ! ps used to calc P profiles & heights on grid staggered
                                            !    East-West (i.e. for VS) relative to ps
real(r8), allocatable :: ps_stagr_lat(:, :) ! ps used to calc P profiles & heights on grid staggered
                                            !    North-South (i.e. for US) relative to ps
! height
! Surface potential; used for calculation of geometric heights.
logical               :: alloc_phis=.true.    ! Flag whether to allocate space for phis[_stagr] 
real(r8), allocatable :: phis(:, :)           ! surface geopotential
real(r8), allocatable :: phis_stagr_lon(:, :) ! surface geopotential staggered as for ps
real(r8), allocatable :: phis_stagr_lat(:, :) ! surface geopotential staggered as for ps

! I'd need 3 of these; 1 for A-grid and 2 for C-grids, so don't keep model_h laying around.
! real(r8), allocatable :: model_h(:, :, :) ! cartesian heights of model levels

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

! array for the linking of obs_kinds (KIND_) to model field TYPE_s
! It's filled in map_kinds
! The max size of KIND_ should come from obs_kind_mod
! These should be dimensioned the same size as the total of state_names_Nd.
integer, dimension(100) :: dart_to_cam_kinds = (/(MISSING_I,iii=1,100)/)
integer, dimension(100) :: cam_to_dart_kinds = (/(MISSING_I,iii=1,100)/)
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
! Initializes class data for CAM model (all the stuff that needs to
! be done once. For now, does this by reading info from a fixed
! name netcdf file. 

integer            :: iunit, io, topog_lons, topog_lats, i, num_lons, num_lats, ncfileid
! calendar types listed in time_manager_mod.f90
integer            :: calendar_type = GREGORIAN
character(len=129) :: errstring
integer            :: lon_index, slon_index, lat_index, slat_index


! Register the module
call register_module(source, revision, revdate)

! setting calendar type
! this information is NOT passed to CAM; it must be set in the CAM namelist
call set_calendar_type(calendar_type)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! set the printed output logical variable to reduce printed output;
! depends on whether this is being called by trans_... (read ens member # from file 'element' )
! or by filter (multiple processes, printout controlled by do_output())
if (file_exist('element')) then
! debug; fix this ugliness
   open(unit = 99, file='element', form = 'formatted')
   read(99,*) ens_member
   close(99)
   do_out = .false.
   if (ens_member == 1) do_out = .true.
else
   do_out = do_output()
end if

! Record the namelist values 
if (do_out) write(logfileunit, nml=model_nml)

! Set the model minimum time step from the namelist seconds and days input
Time_step_atmos = set_time(Time_step_seconds, Time_step_days)
! kdr debug
if (do_out) call print_time(Time_step_atmos)

! read CAM 'initial' file domain info
call check(nf90_open(path = trim(model_config_file), mode = nf90_write, ncid = ncfileid))

! Get sizes of dimensions/coordinates from netcdf and put in global storage.
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
if (do_out) write(*, *) 'CAM size initialized as ', model_size

! Allocate space for longitude and latitude global arrays
! and Allocate space for hybrid vertical coord coef arrays
! height; phis
! allocate(lons(num_lons), lats(num_lats), gw(num_lats), hyai(num_levs+1), &
!          hybi(num_levs+1), hyam(num_levs), hybm(num_levs), &
!          phis(num_lons, num_lats) )


! There's a query of caminput.nc within read_cam_coord for the existence of the field.
! The second argument is a grid_1d_type structure
call read_cam_coord (ncfileid, lon,   'lon     ')
call read_cam_coord (ncfileid, lat,   'lat     ')
call read_cam_coord (ncfileid, gw,    'gw      ')
call read_cam_coord (ncfileid, slon,  'slon    ')
call read_cam_coord (ncfileid, slat,  'slat    ')
call read_cam_coord (ncfileid, lev,   'lev     ')
call read_cam_coord (ncfileid, ilev,  'ilev    ')

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
if (do_out) write(*, *) '# of fields in state vector =  ', nflds,' = sum of ', &
                 state_num_0d ,state_num_1d ,state_num_2d ,state_num_3d

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

call check(nf90_close(ncfileid))

!------------------------------------------------------------------------
! height
! Get lons and _lats from a new netcdf file and test for consistency.
! This subroutines also opens the file for reading fields.
! read CAM 'initial' file domain info
if (file_exist(topog_file)) then
   call check(nf90_open(path = trim(topog_file), mode = nf90_nowrite, ncid = ncfileid))
   if (do_out) write(*, *) 'file_name for surface geopotential height is ', trim(topog_file)

   call read_topog_size(ncfileid, topog_lons, topog_lats)
   ! debug
   if (do_out) write(*,*) 'topog_lons, _lats = ',topog_lons, topog_lats

   num_lons = dim_sizes(find_name('lon     ',dim_names))
   num_lats = dim_sizes(find_name('lat     ',dim_names))

   if (topog_lons /= num_lons .or. topog_lats /= num_lats) then
      write(errstring,'(A,4I4)') 'horizontal dimensions mismatch of initial files and topog ' &
            ,num_lons, topog_lons, num_lats, topog_lats
      call error_handler(E_ERR, 'static_init_model', trim(errstring), source, revision, revdate)
   end if
else
   write(errstring,'(A)') 'topog_file.nc is missing; check for bnd_topo in namelistin' 
   call error_handler(E_ERR, 'static_init_model', trim(errstring), source, revision, revdate)
end if

! Read surface geopotential from topog_file for use in vertical interpolation in height.
! Coordinate order not affected by CAM model version.

if (alloc_phis) allocate (phis(topog_lons, topog_lats))

! Make local space to hold the means
allocate(ens_mean(model_size))

call read_cam_horiz (ncfileid, phis , topog_lons, topog_lats, 'PHIS    ')

call check(nf90_close(ncfileid))

!------------------------------------------------------------------------
! arrays for the linking of obs_kinds (KIND_) to model field TYPE_s; 
!    dart_to_cam_kinds and cam_to_dart_kinds
call map_kinds()

contains 
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)

   integer, intent ( in) :: istatus

   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'static_init_model', &
          trim(nf90_strerror(istatus)), source, revision, revdate)

   end subroutine check

end subroutine static_init_model


   subroutine read_cam_init_size(ncfileid)
!=======================================================================
! subroutine read_cam_init_size(file_name, num_dims, dim_sizes, dim_names)
!
! Gets the number, names, and sizes of field dimensions from a CAM init netcdf file
! in file_name (regardless of dynamical core).

integer,  intent(in)  :: ncfileid

integer :: i,j

!------------------------------------
! learn how many dimensions are defined in this file.
call check(nf90_inquire(ncfileid, num_dims))

! where to deallocate?
allocate (dim_ids(num_dims), dim_names(num_dims), dim_sizes(num_dims))

! Cycle through dimids until there aren't any more
! Dimension ids are sequential integers on the NetCDF file.
do i = 1,num_dims
   dim_ids(i) = i
   call check(nf90_inquire_dimension(ncfileid, i, dim_names(i), dim_sizes(i)))
   if (do_out) write(*,*) 'Dims info = ',i, trim(dim_names(i)), dim_sizes(i)
end do

! Find and store shapes of all the state vector fields.  Grouped by rank of fields into 
! separate s_dim_RANKd arrays.
! Fields with the same rank can have different shapes and still be handled; efgworo(lat,lon) and
! frac(lat,lev) will both have their shapes stored in X_dim_2d.
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

!debug   where will this be written?
if (do_out) then
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

contains 

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus) 
   integer, intent ( in) :: istatus 
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'read_cam_init_size', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine read_cam_init_size

   subroutine trans_coord(ncfileid)
!=======================================================================
! subroutine trans_coord(ncfileid)
!

! Figure out which coordinates are lon, lat, lev, based on CAM version
! from the namelist, which has form #.#[.#[.#]]
! merge; add state_num_#d , and 
integer,            intent(in) :: ncfileid

! local workspace
character (len=4)              :: form_version = '(I0)'
character (len=4)              :: char_version
character (len=129)            :: err_string
integer                        :: part, nchars, tot_chars, i, j, k, varid, next
integer, dimension(4)          :: int_version = (/(0,i=1,4)/)    


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
if (do_out) WRITE(*,'(A,A10,4(I3,2X))') 'model_version, version(1:4) = ' &
                            ,model_version,(int_version(i),i=1,4)
   
! assume cam3.0.7 format to start
! test on version cam3.0.3
coord_order = 2
if (int_version(1) < 3) then
   coord_order = 1
elseif (int_version(1) == 3 .and. int_version(2) == 0 .and. int_version(3) < 3) then
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
   call check(nf90_inq_varid(ncfileid, state_names_3d(i), varid), state_names_3d(i))
   ! Get dimension ids for the dimensions of the field
   call check(nf90_inquire_variable(ncfileid, varid, dimids=f_dimid_3d(1:4,i)))

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
      elseif (dim_names(k) == 'lon' .or. dim_names(k) == 'slon') then
         s_dim_3d  (2,i) = dim_sizes(k)
         s_dimid_3d(2,i) = k
      elseif (dim_names(k) == 'lat' .or. dim_names(k) == 'slat') then
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
   call check(nf90_inq_varid(ncfileid, state_names_2d(i), varid), state_names_2d(i))
   call check(nf90_inquire_variable(ncfileid, varid, dimids=f_dimid_2d(1:3,i)))

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
      elseif (dim_names(k) == 'lon' .or. dim_names(k) == 'slon') then
         ! longitude always comes first on CAM initial files.
         ! Otherwise, I'll need a test like for levs, but more complicated.
         s_dim_2d  (next,i) = dim_sizes(k)
         s_dimid_2d(next,i) = k
         next = 2
      elseif (dim_names(k) == 'lat' .or. dim_names(k) == 'slat' ) then
         s_dim_2d  (next,i) = dim_sizes(k)
         s_dimid_2d(next,i) = k
      end if
!      cycle Alldim2
!         end if
!      end do
   end do Alldim2
   if (   s_dim_2d(1,i) == 0 .or.  s_dim_2d(2,i) == 0 ) then
      call error_handler(E_ERR, 'trans_coord', &
          'num_[lons,lats,levs] was not assigned and = 0' , source, revision, revdate)
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
   call check(nf90_inq_varid       (ncfileid, state_names_1d(i), varid), state_names_1d(i))
   call check(nf90_inquire_variable(ncfileid, varid, dimids=f_dimid_1d(1:2,i)))

   Alldim1: do j = 1,2       ! time and 1 space
      k = f_dimid_1d(j,i)
      f_dim_1d(j,i) = dim_sizes(k)
!      do k = 1,num_dims
!         if (dim_ids(k) == f_dimid_1d(j,i)) then
      if (dim_names(k) == 'lon' .or. dim_names(k) == 'slon' .or. &
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
      write(err_string, '(A,I3,A)') ' state 1d dimension(',i,') was not assigned and = 0' 
      call error_handler(E_ERR, 'trans_coord',trim(err_string), source, revision, revdate) 
   end if
end do

contains 

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus, string) 
   integer, intent (in) :: istatus 
   character(len=*), optional, intent(in) :: string
   character(len=255) :: errstring

   if (present(string)) then
      write(errstring, *)trim(adjustl(string))//' '//trim(adjustl(nf90_strerror(istatus)))
   else
      write(errstring, *)                            trim(adjustl(nf90_strerror(istatus)))
   end if
      
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'trans_coord', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine trans_coord


   subroutine read_topog_size(ncfileid, num_lons, num_lats)
                                           
!=======================================================================
! subroutine read_topog_size(file_name, num_lons, num_lats)
!
! Gets the number of lons and lats from a CAM surface netcdf file
! in file_name

integer, intent(out)          :: num_lons, num_lats

integer                       :: londimid, levdimid, latdimid, ncfileid, ncfldid
integer                       :: phis_dimids(3)
character (len=NF90_MAX_NAME) :: clon,clat

! get field id
call check(nf90_inq_varid(ncfileid, 'PHIS', ncfldid))
! get dimension 'id's
call check(nf90_inquire_variable(ncfileid, ncfldid, dimids = phis_dimids))
! call check(nf90_inq_dimid(ncfileid, 'lon', londimid))
! call check(nf90_inq_dimid(ncfileid, 'lat', latdimid))

! get dimension sizes
call check(nf90_inquire_dimension(ncfileid, phis_dimids(1), clon , num_lons ))
call check(nf90_inquire_dimension(ncfileid, phis_dimids(2), clat , num_lats )) 

! check for correct order will be done (implicitly) in calling routine.

contains

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
   integer, intent ( in) :: istatus
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'read_topog_size', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine read_topog_size



   subroutine read_cam_horiz(ncfileid, var, dim1, dim2, cfield)
!======================================================
! should be called with cfield = a 2D record variable  (time,lat,lon):

implicit none                                                                                         
!------------------------------------------------------
integer,                         intent(in)  :: ncfileid, dim1, dim2
real(r8), dimension(dim1, dim2), intent(out) :: var
character (len=8),               intent(in)  :: cfield

!------------------------------------------------------
integer :: i, j, ifld             ! grid indices
integer :: ncfldid
integer :: n,m, slon_index, slat_index, lat_index, lon_index

if (do_out) PRINT*,'read_cam_horiz; reading ',cfield
call check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid), cfield)
call check(nf90_get_var(ncfileid, ncfldid, var, start=(/1,1,1/), &
           count=(/dim1, dim2, 1/)), trim(cfield))

if (do_out) PRINT*,'read_cam_horiz; reading ',cfield,' using id ',ncfldid, dim1, dim2

! assign values to phis grids for use by the rest of the module.
if (cfield == 'PHIS    ') then

   slon_index = find_name('slon    ',dim_names)
   slat_index = find_name('slat    ',dim_names)
   lat_index  = find_name('lat     ',dim_names)
   lon_index  = find_name('lon     ',dim_names)

   phis = var

   if (slon_index /= 0) then
      if (alloc_phis) allocate (phis_stagr_lon (dim_sizes(slon_index), dim_sizes( lat_index)))
      do n=1,dim_sizes( lat_index)
         phis_stagr_lon(1,n) = .5 * (phis(1,n) + phis(dim_sizes( lon_index),n))
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

contains

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus,string)
   integer,                    intent(in) :: istatus
   character(len=*), optional, intent(in) :: string
   character(len=255) :: errstring

   if (present(string)) then
      write(errstring, *)trim(adjustl(string))//' '//trim(adjustl(nf90_strerror(istatus)))
   else
      write(errstring, *)                            trim(adjustl(nf90_strerror(istatus)))
   end if
      
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'read_cam_horiz', &
          errstring, source, revision, revdate)
   end subroutine check

end subroutine read_cam_horiz



   subroutine nc_read_model_atts(att, att_vals, nflds, ncfileid)
!=======================================================================
! subroutine nc_read_model_atts(att, att_vals, nflds)
!
! reads the value of an attribute for each of the fields in cflds.
!
! should be called with att = one of the attributes from the program variable
! input file, which will be written to the Posterior and Prior.nc files

!----------------------------------------------------------------------
! Local workspace
integer :: i, nchars, ierr
integer :: ncfileid, ncfldid, ncattid, att_type

!----------------------------------------------------------------------
integer,                                intent(in)  :: nflds
character (len=*),                      intent(in)  :: att 
character (len=128), dimension(nflds), intent(out)  :: att_vals 

! open CAM 'initial' file 
! DONE ALREADY in static_init_model
! call check(nf90_open(path = trim(model_config_file), mode = nf90_nowrite, &
!           ncid = ncfileid))

! read CAM 'initial' file attribute desired
if (do_out) PRINT*,'reading ',att
do i = 1,nflds
   call check(nf90_inq_varid(ncfileid, trim(cflds(i)), ncfldid))
! could be inquire_attribute
! 
   ierr = nf90_inquire_attribute(ncfileid, ncfldid, trim(att), att_type, nchars, ncattid) 

   if (ierr == nf90_noerr) then
      att_vals(i)(1:64) = '                                                                '
      att_vals(i)(65:128) = '                                                                '
      call check(nf90_get_att(ncfileid, ncfldid, trim(att) ,att_vals(i) ))
      if (do_out) WRITE(*,'(I6,2A)') ncfldid, cflds(i), trim(att_vals(i))
   else
      WRITE(*,*) ncfldid, cflds(i), 'NOT AVAILABLE'
   end if
end do

contains 

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus) 
   integer, intent ( in) :: istatus 
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_read_model_atts', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine nc_read_model_atts

   subroutine read_cam_coord(ncfileid, var, cfield)
!=======================================================================
! subroutine read_cam_coord(ncfileid, var , cfield)
!
! read CAM 'initial' file domain info
! should be called with cfield = one of :
!          (/'lat     ','lon     ','gw      '
!           ,'hyai    ','hybi    ','hyam    ','hybm    '/)

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
integer :: num_atts, keep_atts  
integer                                     :: att_type 
character (len=nf90_max_name)               :: att_name
character (len=nf90_max_name), allocatable  :: att_names(:)
character (len=128),           allocatable  :: att_vals(:)


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
   call check(nf90_inq_attname(ncfileid, ncfldid, i, att_name))

! FV initial files have coordinates with attributes that are numerical, not character.
! (_FillValue).  These are not used because the coordinates are dimensioned exactly
! the right size.  I'll test for the type of att, and if it's not char, I'll ignore it
! and reduce the num_atts by 1.

! Otherwise I need a var@atts_type and separate var%atts_vals_YYY for each NetCDF 
! external type (6 of them) I might run into.

   call check(nf90_inquire_attribute(ncfileid, ncfldid, att_name, xtype=att_type) )

   if (att_type == nf90_char) then
      keep_atts = keep_atts + 1
      att_names(keep_atts) = att_name
      call check(nf90_get_att(ncfileid, ncfldid, att_name, att_vals(keep_atts)))

   else
      if (do_out) write(*,*) '                ignoring attribute ',att_name,    &
                    ' because it is not a character type'
   endif
end do

! allocate space for this grid structure, and put vector length and # attributes in structure
call init_grid_1d_instance(var, coord_size, keep_atts)

var%label = cfield
var%dim_id = coord_dims(1)

do i = 1,keep_atts
   var%atts_names(i) = att_names(i)
   var%atts_vals(i)  = att_vals(i)
enddo

! call check(nf90_get_var(ncfileid, ncfldid, var%vals(1:coord_size) ,start=(/1/) &
call check(nf90_get_var(ncfileid, ncfldid, var%vals, start=(/1/) &
    ,count=(/coord_size/) ))

if (do_out) then
   PRINT*,'reading ',cfield,' using id ',ncfldid,' size ',coord_size
   WRITE(*,*) (var%vals(i),i=1,coord_size)
end if

deallocate(att_names, att_vals)

contains 

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus) 
   integer, intent ( in) :: istatus 
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'read_cam_coord', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

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
character (len = 129) :: errstring

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
end do

! 2D fields
do i=1,state_num_2d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_2d(i)
   TYPE_2D(i) = nfld
   if (state_names_2d(i) == 'PS      ') TYPE_PS    = nfld
   if (state_names_2d(i) == 'TBOT    ') TYPE_TBOT  = nfld
   if (state_names_2d(i) == 'TS      ') TYPE_TS    = nfld
   if (state_names_2d(i) == 'TSOCN   ') TYPE_TSOCN = nfld
! These are on topog file, not initial file anymore
!   if (state_names_2d(i) == 'PHIS    ') TYPE_PHIS  = nfld
!   if (state_names_2d(i) == 'SGH     ') TYPE_SGH   = nfld
!   if (state_names_2d(i) == 'PBLH    ') TYPE_PBLH  = nfld
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

end do

if (nfld .ne. nflds) then
   write(errstring, *) 'nfld = ',nfld,', nflds = ',nflds,' must be equal '
   call error_handler(E_ERR, 'order_state_fields', errstring, source, revision, revdate)
elseif (do_out) then
   write(logfileunit,'(/A/)') 'State vector is composed of '
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
   write(logfileunit,'(/A)') 'TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q = ' 
   write(logfileunit,'((8(I8,1X)))') TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q
end if

return

end subroutine order_state_fields


   subroutine map_kinds()
!=======================================================================
! subroutine map_kinds()

! ? Should this be a function instead; removes need to dimension obs_loc_in arbitrarily
!   and wastefully.  But then it's called millions of times, instead of accessing an
!   array that's defined once.

! Makes an array of 'locations within the state vector'
! of  all the available obs kinds that come from obs_kind_mod. 
! The obs kind that's needed will be the index into this array,
! the corresponding value will be the position of that field (not individual variable) 
! within the state vector according to state_name_Xd.  
! This subroutine will be called from static_init_model, so it will not have to be 
! recomputed for every obs.
! Also maps the local model_mod TYPE_s onto the DART KIND_s by the same mechanism.

! other KIND_ possibilities are listed after the 'use obs_kind_mod' statement

integer :: i

! 2D fields
dart_to_cam_kinds(KIND_SURFACE_PRESSURE) = TYPE_PS
cam_to_dart_kinds(TYPE_PS) = KIND_SURFACE_PRESSURE
! dart_to_cam_kinds(KIND_SURFACE_TEMPERATURE  ?  ) = TYPE_TS
! dart_to_cam_kinds(KIND_SEA_SURFACE_TEMPERATURE  ?  ) = TYPE_TSOCN

! 3D fields
dart_to_cam_kinds(KIND_TEMPERATURE)        = TYPE_T
dart_to_cam_kinds(KIND_U_WIND_COMPONENT)   = TYPE_U
dart_to_cam_kinds(KIND_V_WIND_COMPONENT)   = TYPE_V
dart_to_cam_kinds(KIND_SPECIFIC_HUMIDITY)  = TYPE_Q
dart_to_cam_kinds(KIND_CLOUD_LIQUID_WATER) = TYPE_CLDLIQ
dart_to_cam_kinds(KIND_CLOUD_ICE)          = TYPE_CLDICE
! dart_to_cam_kinds(KIND_CLOUD_WATER  ?  ) = TYPE_LCWAT

cam_to_dart_kinds(TYPE_T)      = KIND_TEMPERATURE
cam_to_dart_kinds(TYPE_U)      = KIND_U_WIND_COMPONENT
cam_to_dart_kinds(TYPE_V)      = KIND_V_WIND_COMPONENT
cam_to_dart_kinds(TYPE_Q)      = KIND_SPECIFIC_HUMIDITY
cam_to_dart_kinds(TYPE_CLDLIQ) = KIND_CLOUD_LIQUID_WATER
cam_to_dart_kinds(TYPE_CLDICE) = KIND_CLOUD_ICE
! cam_to_dart_kinds(TYPE_LCWAT) = KIND_CLOUD_WATER  ?  


if (do_out) then
   write(*,*) 'OBS_KIND   FIELD_TYPE'
   do i=1,100
      if (dart_to_cam_kinds(i) /= MISSING_I) write(*,'(2I8)') i, dart_to_cam_kinds(i)
   end do
end if

! In the future, if fields are not ordered nicely, or if users are specifying
! correspondence of obs fields with state fields, I may want code like:
! The max size of KIND_ should come from obs_kind_mod
! do i=1,state_num_3d
!    if (state_names_3d(i)(1:1) == 'T' .and. &
!        KIND_TEMPERATURE <= 100) ) dart_to_cam_kinds(KIND_TEMPERATURE) = TYPE_3D(i)
! end do 

return

end subroutine map_kinds

! End of static_init_model section
!#######################################################################

! module I/O to/from DART and files

   subroutine read_cam_init(file_name, var)
!=======================================================================
! subroutine read_cam_init(file_name, var)
!

character(len = *), intent(in)  :: file_name
type(model_type),   intent(out) :: var

! Local workspace
integer :: i, k, n, m, ifld  ! grid and constituent indices
integer :: ncfileid, ncfldid
real(r8), allocatable :: temp_3d(:,:,:), temp_2d(:,:)

!----------------------------------------------------------------------
! read CAM 'initial' file domain info
call check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid))

! The temporary arrays into which fields are read are dimensioned by the largest values of
! the sizes of the dimensions listed in f_dim_RANKd
! f_dim_max contents assume that time is always the last dimension on NetCDF files,
! so f_dim_max(4,3) and f_dim_max(3,2) are the non-spatial dimensions to ignore here.
allocate (temp_3d(f_dim_max(1,3),f_dim_max(2,3),f_dim_max(3,3)) &
         ,temp_2d(f_dim_max(1,2),f_dim_max(2,2)) )

! Imported from trans_pv_sv.  Nancy; is this a better place for it?
! call init_model_instance(var)

! read CAM 'initial' file fields desired

ifld = 0
!0d fields; scalars are recognized and handled differently than vectors
!           by NetCDF
! kdr; nf90_put_var was probably bombing because
!      netcdf recognized that it was writing a scalar, called an f77
!      scalar put_var, and then choked when it saw the count = ...
do i= 1, state_num_0d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   if (do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 1D; TIME(=1)
   call check(nf90_get_var(ncfileid, ncfldid, var%vars_0d(i) ))
end do
             

!1d fields
do i= 1, state_num_1d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
! debug; remove next two lines when I'm sure I don't need them.  Also for 2d and 3d
! done in trans_coord already
!   call check(nf90_inquire_variable(ncfileid, ncfldid, dimids=f_dimid_1d(1:2,i)))
   if (do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
   
   ! s_dim_1d should = f_dim_1d
   call check(nf90_get_var(ncfileid, ncfldid, var%vars_1d(1:s_dim_1d(i), i) &
             ,start=(/1,1/) ,count=(/f_dim_1d(1,i), 1/) ))
end do

!2d fields
do i= 1, state_num_2d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
! done in trans_coord already
!   call check(nf90_inquire_variable(ncfileid, ncfldid, dimids=f_dimid_2d(1:3,i)))
   if (do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
   ! fields on file are 3D; lon, lat, (ususally) and then  TIME(=1)
   ! Need to use temp_Nd here too; I am coding for not knowing what 2 of the 3 dimensions the
   ! 2d fields will have.
   call check(nf90_get_var(ncfileid, ncfldid, temp_2d(1:f_dim_2d(1,i), 1:f_dim_2d(2,i))  &
                              ,start=(/1,1,1/) ,count=(/f_dim_2d(1,i),   f_dim_2d(2,i), 1/) ))
   if (s_dim_2d(1,i) == f_dim_2d(1,i)) then
      var%vars_2d(1:s_dim_2d(1,i), 1:s_dim_2d(2,i),i) = &
          temp_2d(1:f_dim_2d(1,i), 1:f_dim_2d(2,i)  )
   elseif (s_dim_2d(1,i) == f_dim_2d(2,i)) then
      do k=1,s_dim_2d(1,i)
      do m=1,s_dim_2d(2,i)   ! first temp dim is inner loop for faster reads
         var%vars_2d(k,m,i) = temp_2d(m,k)
      end do
      end do
   else
      ! error
   end if
end do

! 3d fields
do i=1, state_num_3d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
! done in trans_coord already
!   call check(nf90_inquire_variable(ncfileid, ncfldid, dimids=f_dimid_3d(1:4,i)))
   if (do_out) PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 4D; lon, lev, lat, TIME(=1) 
!                     or; lon, lat, lev, TIME
   call check(nf90_get_var(ncfileid, ncfldid                                                 &
             ,temp_3d(1:f_dim_3d(1,i), 1:f_dim_3d(2,i), 1:f_dim_3d(3,i)) ,start=(/1,1,1,1/)  &
             ,count=(/  f_dim_3d(1,i),   f_dim_3d(2,i),   f_dim_3d(3,i),1 /) ))
! read into dummy variable, then repackage depending on coord_order
! put it into lon,lev,lat order for continuity with Tim's/historical.
   if (coord_order == 1) then
!     lon,lev,lat as in original CAM
      do n=1,s_dim_3d(3,i)  ! lats
      do k=1,s_dim_3d(1,i)  ! levs 
      do m=1,s_dim_3d(2,i)  ! lons/first file dim   put in inner loop for faster reads
         var%vars_3d(k,m,n, i) = temp_3d(m,k,n)
      end do
      end do
      end do
   elseif (coord_order == 2) then
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

call check(nf90_close(ncfileid))

deallocate (temp_3d,temp_2d)

contains 

   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus) 
   integer, intent ( in) :: istatus 
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'read_cam_init', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine read_cam_init



   subroutine write_cam_coord_def(ncFileID, c_name, coord, dim_id, c_id)
!=======================================================================================
! subroutine write_cam_coord_def(ncFileID, c_name, coord, dim_id, c_id)

character (len=8),  intent(in)  :: c_name
integer,            intent(in)  :: ncFileID, dim_id
type(grid_1d_type), intent(in)  :: coord
integer,            intent(out) :: c_id

integer  :: i, nch

call check(nf90_def_var(ncFileID, name=c_name, xtype=nf90_double, dimids=dim_id, &
                        varid=c_id) )
if (do_out) write(*,'(/A,A)') 'write_cam_coord_def;  ', trim(c_name)

do i=1,coord%num_atts
   if (do_out) then
!      nch = len_trim(coord%atts_vals(i))
!                 i,trim(coord%atts_names(i)),' ', coord%atts_vals(i)(1:nch)
      write(*,*) '   i, att_name, att_val', &
                 i,trim(coord%atts_names(i)),' ', trim(coord%atts_vals(i))
   endif
   call check(nf90_put_att(ncFileID, c_id, coord%atts_names(i), coord%atts_vals(i)))
end do

return
contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
   integer, intent ( in) :: istatus
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'write_cam_coord_def', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end subroutine write_cam_coord_def

   subroutine write_cam_init(file_name, var)
!=======================================================================
! subroutine write_cam_init(file_name, var)

! write CAM 'initial' file fields that have been updated

character (len = *), intent(in) :: file_name
type(model_type),    intent(in) :: var

integer               :: i, k, n, m, ifld, ncfileid, ncfldid, f_dim1, f_dim2
real(r8), allocatable :: temp_3d(:,:,:), temp_2d(:,:)

! Read CAM 'initial' file domain info
call check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid))
 
! The temporary arrays into which fields are read are dimensioned by the largest values of
! the sizes of the dimensions listed in coord_RANKd
! debug
! ! ! This may not work for writing fields that have smaller sizes!
!     Unless array section is explicitly given in array argument
allocate (temp_3d(f_dim_max(1,3),f_dim_max(2,3),f_dim_max(3,3)))
allocate (temp_2d(f_dim_max(1,2),f_dim_max(2,2)))


if (do_out) write(*,*) 'write_cam_init; f_dim_max(:2) = ',f_dim_max(1,2),f_dim_max(2,2)

ifld = 0
! 0d fields are first
! kdr; Tim says that nf90_put_var was probably bombing because
!      netcdf recognized that it was writing a scalar, called an f77
!      scalar put_var, and then choked when it saw the count = ...
!      So, this padding is probably not necessary.
do i = 1, state_num_0d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid, var%vars_0d(i) ))
end do 

! 1d fields 
do i = 1, state_num_1d
   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid, var%vars_1d(1:s_dim_1d(i), i), &
                                  start=(/1, 1/), count = (/s_dim_1d(i), 1/)))
end do 

! 2d fields ; tricky because coordinates may have been rearranged to put them in the order
! of (lev, lon, lat) choosing only 2.  The original coordinate order is in f_dimid_2d.
! debug 
! if (do_out) write(*,'(A/A)') 'write_cam_init 2D',' i f_dim1 f_dim2'

do i = 1, state_num_2d
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
!      if (do_out) write(*,'(3I3,2F14.6)') i, f_dim1, f_dim2  , temp_2d(1,1), temp_2d(f_dim1, f_dim2)
   elseif (f_dimid_2d(1,i) == s_dimid_2d(2,i)) then
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
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
!   writing part of temp_2d defined by start and count; don't need indices.
!   call check(nf90_put_var(ncfileid, ncfldid, temp_2d(1:f_dim1, 1:f_dim2), &
   call check(nf90_put_var(ncfileid, ncfldid, temp_2d(1:f_dim1,1:f_dim2), &   
                           start=(/1, 1, 1/), count = (/f_dim1, f_dim2, 1/)))
!   call check(nf90_put_var(ncfileid, ncfldid, var%vars_2d(1:s_dim_2d(1,i), 1:s_dim_2d(1,i), i), &
!                               start=(/1, 1, 1/), count = (/s_dim_2d(1,i),   s_dim_2d(2,i), 1/)))
end do 

! 3d fields; all 3 coordinates are present, and the order for model_mod fields is always the same.
do i = 1, state_num_3d
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
   elseif (coord_order == 2) then
      ! lon,lat,lev as in new CAM
      do n=1,s_dim_3d(3,i)
      do m=1,s_dim_3d(2,i)
      do k=1,s_dim_3d(1,i)
         temp_3d(m,n,k) = var%vars_3d(k,m,n,i)
      end do
      end do
      end do
   end if

   ifld = ifld + 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid                                               &
             ,temp_3d(1:f_dim_3d(1,i), 1:f_dim_3d(2,i), 1:f_dim_3d(3,i)) &
             ,start=(/1,1,1,1/) ,count=(/f_dim_3d(1,i), f_dim_3d(2,i), f_dim_3d(3,i), 1/) ))
end do

call check(nf90_close(ncfileid))

deallocate (temp_3d, temp_2d)

contains 
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)

   integer, intent ( in) :: istatus

   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'write_cam_init', &
          trim(nf90_strerror(istatus)), source, revision, revdate)

   end subroutine check

end subroutine write_cam_init


   subroutine get_state_meta_data(index_in, location, var_kind)
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
! See order_state_fields for the KIND_s (and corresponding model_mod TYPE_s).
!
! This is not a function because the more general form of the call has a second 
! intent(out) optional argument var_kind.  Maybe a functional form should be added?

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_kind 

integer  :: which_vert 
integer  :: i, indx, index_1, index_2, index_3, nfld
integer  :: box, slice

real(r8) :: lon_val, lat_val, lev_val

character(len=129) :: errstring
character(len=8)   :: dim_name

! In order to find what variable this is, and its location, I must subtract the individual 
! components of the state vector, since they may have varying sizes.
! Save the original index.   index_in will be < 0 if it's an identity obs (called from convert_vert)
indx = abs(index_in)
which_vert = MISSING_I
index_1 = 0; index_2 = 0; index_3 = 0; nfld = 0
lon_val = MISSING_R8; lat_val = MISSING_R8; lev_val = MISSING_R8

! Cycle through 0d state variables
do i=1,state_num_0d
   nfld = nfld + 1
   if (indx == 1 ) then
      which_vert = -2
      goto 10
   else
      indx = indx - 1
   end if
end do

! Cycle through 1d state variables
! Note that indices of fields can have varying dimensions.
do i=1,state_num_1d
   nfld = nfld + 1
   if (indx > s_dim_1d(i) ) then
      indx = indx - s_dim_1d(i)
   else
      ! We've found the desired field; now find lat, lon or lev of indx
      if (index_in > 0) then
         dim_name = dim_names(s_dimid_1d(i))
         if (dim_name == 'lev     ') then
            lev_val = real(indx)
         else
            call coord_val(dim_name, indx, lon_val, lat_val, lev_val)
         endif
   
         which_vert = which_vert_1d(i)
      end if
      goto 10
   end if
end do

! Cycle through 2d state variables.
! Note that indices of fields can have varying dimensions.
do i=1,state_num_2d
   nfld = nfld + 1
   slice = s_dim_2d(1,i) * s_dim_2d(2,i)
   if (indx > slice ) then
      indx = indx - slice
   else
      ! We've found the desired field. 
      ! Now find lat and/or lon and/or lev of indx if called by assim_tools_mod:filter_assim
      if (index_in > 0) then
         ! # second dimension rows to subtract off; temporary value for index_2
         index_2 = (indx -1) / s_dim_2d(1,i) 
         index_1 = indx - (index_2 * s_dim_2d(1,i))
         dim_name = dim_names(s_dimid_2d(1,i))
         ! Find the coordinate value (i.e. 270.5) of the first dimension index (i.e. 54)
         if (dim_name == 'lev     ') then
            lev_val = real(index_1)
         else
            call coord_val(dim_name, index_1, lon_val, lat_val, lev_val)
         endif

         ! index_2 of the variable in question is 1 more than the # subtracted of to get index_1
         index_2 = index_2 + 1
         dim_name = dim_names(s_dimid_2d(2,i))
         if (dim_name == 'lev     ') then
            lev_val = real(index_2)
         else
            call coord_val(dim_name, index_2, lon_val, lat_val, lev_val)
         endif

         which_vert = which_vert_2d(i)

      end if
      goto 10
   end if
end do

! Cycle through 3d state variables
! Note that indices of fields can have varying dimensions.
do i=1,state_num_3d
   nfld = nfld + 1
   box = s_dim_3d(1,i) * s_dim_3d(2,i) * s_dim_3d(3,i)
   if (indx > box ) then
      indx = indx - box
   else
      ! We've found the desired field. 
      ! Now find lat and/or lon and/or lev of indx if called by assim_tools_mod:filter_assim
      if (index_in > 0) then
         ! # (first x second dimension) slices to subtract off in order to find the current slice
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

         ! index_2 of the variable in question is one more than the number subtracted of to get index_1
         index_2 = index_2 + 1
         dim_name = dim_names(s_dimid_3d(2,i))
         call coord_val(dim_name, index_2, lon_val, lat_val, lev_val)

         ! index_3 of the variable in question is one more than the number subtracted of to get index_1
         index_3 = index_3 + 1
         dim_name = dim_names(s_dimid_3d(3,i))
         call coord_val(dim_name, index_3, lon_val, lat_val, lev_val)
   
         which_vert = which_vert_3d(i)
      end if
      goto 10
   end if
end do


10 continue

! If the type is wanted, return it
if (present(var_kind)) then
   if (index_in < 0) then
      ! used by 
      var_kind = nfld
   else if (index_in > 0) then
      ! used by call from assim_tools_mod:filter_assim
      var_kind = cam_to_dart_kinds(nfld)
      location = set_location(lon_val, lat_val, lev_val, which_vert)  
   end if
end if

end subroutine get_state_meta_data




   subroutine ens_mean_for_model(filter_ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: filter_ens_mean(:)

ens_mean = filter_ens_mean

! Fill ps, ps_stagr_lxx if not filled yet.
! WATCH OUT that it's not still filled with something other than ens_mean
call set_ps_arrays(ens_mean)


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
! Returns the the time step of the model. In the long run should be repalced
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
! merge
integer :: MemberDimID, StateVarDimID, TimeDimID,ScalarDimID
integer :: xVarID,StateVarID, StateVarVarID 
integer :: P_id(num_dims)
integer :: i, ifld, dim_id, g_id
integer :: grid_id(grid_num_1d) 
character(len=129)    :: errstring 
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

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_Redef(ncFileID))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies
!-------------------------------------------------------------------------------

call check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID))
call check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID))

if ( TimeDimID /= unlimitedDimId ) then
  write(errstring,*)'Time dimension ID ',TimeDimID,'must match Unlimited Dimension ID ',unlimitedDimId
  call error_handler(E_ERR,'nc_write_model_atts', errstring, source, revision, revdate)
end if

!-------------------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-------------------------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                        len=model_size, dimid = StateVarDimID))

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","CAM"))

! how about namelist input? might be nice to save ...

!-------------------------------------------------------------------------------
! Define the new dimensions IDs
!-------------------------------------------------------------------------------

! All the dim_ids on caminput.nc were read in, and will be written out here,
! and some will be used to define the coordinate variables below.
! They have different dimids for this file than they had for caminput.nc
! P_id serves as a map between the 2 sets.
if (do_out) write(*,*) 'dimens ,name      ,size, cam dim_id, P[oste]rior id'
do i = 1,num_dims
   if (trim(dim_names(i)) /= 'time')  call check(nf90_def_dim                    &
      (ncid=ncFileID, name=trim(dim_names(i)), len=dim_sizes(i), dimid=P_id(i)))
   if (do_out) write(*,'(A7,A8,3(I4,5X))') i,trim(dim_names(i)),dim_sizes(i), dim_ids(i), P_id(i)
end do
call check(nf90_def_dim(ncid=ncFileID, name="scalar",   len = 1,   dimid = ScalarDimID))

!-------------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and their attributes
!-------------------------------------------------------------------------------

! grid longitudes, latitudes, levels, and other coordinates.
! grid_id() is filled here; it's the dimid of the desired coordinate *on this P_Diag.nc file*.  
! It's used to write coordinates.  dim_ids keeps track of the dimensions of coordinates 
! and fields on caminput.nc files.  There's some overlap of names, unfortunately.
! The argument after the 'xxx    ' label is a structure with all the relevant info in it.

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

if (do_out) then
do i=1,grid_num_1d
   write(*,*) 'grid_ = ', i, grid_id(i), grid_names_1d(i)
end do
end if

if ( output_state_vector ) then

! CAM3; need to adapt this to state_long_names, state_units, etc?


   !----------------------------------------------------------------------------
   ! Create attributes for the state vector
   !----------------------------------------------------------------------------

   ! Define the state vector coordinate variable
   call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID))
   call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"))
   call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical") )
   call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)))

   ! Define the actual state vector
   call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), varid=StateVarID))
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"))
   call check(nf90_put_att(ncFileID, StateVarId, "vector_to_prog_var","CAM"))

   ! Leave define mode so we can fill 
   call check(nf90_enddef(ncfileID))

   ! Fill the state variable coordinate variable
   call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ))

else

   !----------------------------------------------------------------------------
   ! Create the (empty) Variables and the Attributes
   !----------------------------------------------------------------------------

   ! 0-d fields
   ifld = 0
   do i = 1,state_num_0d
      ifld = ifld + 1
      call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ MemberDimID, unlimitedDimID /), &
                 varid  = xVarID))
      call check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)))
      call check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)))
!       call check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)))
   end do

   ! 1-d fields
   do i = 1,state_num_1d
      ifld = ifld + 1
      call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ P_id(s_dimid_1d(1)), MemberDimID, unlimitedDimID /), &
                 varid  = xVarID))
      call check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)))
      call check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)))
!       call check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)))
   end do

   ! 2-d fields
   do i = 1,state_num_2d
      ifld = ifld + 1
      call check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ P_id(s_dimid_2d(1,i)), P_id(s_dimid_2d(2,i)),            &
                             MemberDimID, unlimitedDimID /),                          &
                 varid  = xVarID))
      call check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)))
      call check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)))
!       call check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)))
   end do

   ! 3-d fields
   do i = 1,state_num_3d
      ifld = ifld + 1
      call check(nf90_def_var                                                                 &
           (ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real,                           &
            dimids = (/ P_id(s_dimid_3d(1,i)), P_id(s_dimid_3d(2,i)), P_id(s_dimid_3d(3,i)),  &
                        MemberDimID, unlimitedDimID /),                                       &
            varid  = xVarID))
      call check(nf90_put_att(ncFileID, xVarID, "long_name", state_long_names(ifld)))
      call check(nf90_put_att(ncFileID, xVarID, "units", state_units(ifld)))
!       call check(nf90_put_att(ncFileID, xVarID, "units_long_name", state_units_long_names(ifld)))
   end do

   ! Leave define mode so we can fill variables
   call check(nf90_enddef(ncfileID))

!-------------------------------------------------------------------------------
! Fill the coordinate variables
! Each 'vals' vector has been dimensioned to the right size for its coordinate.  
! The default values of 'start' and 'count'  write out the whole thing.
!-------------------------------------------------------------------------------

if (do_out) write(*,*) 'nc_write_model_atts; filling coords'

if (lon%label  /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('lon     ',grid_names_1d)),  lon%vals ))
if (lat%label  /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('lat     ',grid_names_1d)),  lat%vals ))
if (lev%label  /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('lev     ',grid_names_1d)),  lev%vals ))
if (gw%label   /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('gw      ',grid_names_1d)),   gw%vals ))
if (hyam%label /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('hyam    ',grid_names_1d)), hyam%vals ))
if (hybm%label /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('hybm    ',grid_names_1d)), hybm%vals ))
if (hyai%label /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('hyai    ',grid_names_1d)), hyai%vals ))
if (hybi%label /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('hybi    ',grid_names_1d)), hybi%vals ))
if (slon%label /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('slon    ',grid_names_1d)), slon%vals ))
if (slat%label /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('slat    ',grid_names_1d)), slat%vals ))
if (ilev%label /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('ilev    ',grid_names_1d)), ilev%vals ))
if (P0%label   /= '        ')  &
    call check(nf90_put_var(ncFileID, grid_id(find_name('P0      ',grid_names_1d)),   P0%vals ))

end if

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call check(nf90_sync(ncFileID))

write (*,*)'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...' 

contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
   integer, intent ( in) :: istatus
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_atts', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

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
integer                            :: dim_lev, dim_lat  ! old/new CAM coord order accounting

!-----------------------------------------------------------------------------------------
type(model_type) :: Var

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, ncfldid
integer :: ifld, i

ierr = 0     ! assume normal termination

!-------------------------------------------------------------------------------
! CAM storage bounds are 'tight' -- no indices needed
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! then get all the Variable ID's we need.
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

if ( output_state_vector ) then

   call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)))                               

else

   !----------------------------------------------------------------------------
   ! Fill the variables
   !----------------------------------------------------------------------------

   call init_model_instance(Var)     ! Explicity released at end of routine. 

   call vector_to_prog_var(statevec,  Var)
   
   ifld = 0
   ZeroDVars : do i = 1, state_num_0d    
      ifld = ifld + 1
      call check(NF90_inq_varid(ncFileID, trim(cflds(ifld)), ncfldid))
      call check(nf90_put_var( ncFileID, ncfldid, Var%vars_0d(i), &
                 start=(/ copyindex, timeindex /) ))
!, &
!                 count=(/1, 1/) ))
   end do ZeroDVars

   ! 'start' and 'count' are needed here because Var%vars_Nd are dimensioned by the largest
   ! values for the dimensions of the rank N fields, but some of the fields are smaller than that.
   OneDVars : do i = 1, state_num_1d    
      ifld = ifld + 1
      call check(NF90_inq_varid(ncFileID, trim(cflds(ifld)), ncfldid))
      call check(nf90_put_var( ncFileID, ncfldid,               &
           Var%vars_1d(1:s_dim_1d(  i), i),                      &
              start=(/ 1              ,copyindex, timeindex /), &
              count=(/   s_dim_1d(  i),1        , 1/) ))
   end do OneDVars

   ! Write out 2D variables as 2 of (lev,lon,lat), in that order, regardless of caminput.nc 
   ! coordinate order
   ! The sizes can be taken from s_dim_2d, even though the s_dimid_2d don't pertain to this 
   ! P_Diag.nc file, because the dimids were mapped correctly in nc_write_model_atts.
   TwoDVars : do i = 1, state_num_2d    
      ifld = ifld + 1
      call check(NF90_inq_varid(ncFileID, trim(cflds(ifld)), ncfldid))
      call check(nf90_put_var( ncFileID, ncfldid,                                &
           Var%vars_2d(1:s_dim_2d(1,i),1:s_dim_2d(2,i), i),                       &
              start=(/ 1              ,1              , copyindex, timeindex /), &
              count=(/   s_dim_2d(1,i),  s_dim_2d(2,i), 1        , 1/) ))
   end do TwoDVars

   ! Write out 3D variables as (lev,lon,lat) regardless of caminput.nc coordinate order
   ThreeDVars : do i = 1,state_num_3d
      ifld = ifld + 1
      call check(NF90_inq_varid(ncFileID, trim(cflds(ifld)), ncfldid))
      call check(nf90_put_var( ncFileID, ncfldid,                                            &
           Var%vars_3d(1:s_dim_3d(1,i),1:s_dim_3d(2,i),1:s_dim_3d(3,i),i)                    &
              ,start=(/1              ,1              ,1              ,copyindex,timeindex /)&
              ,count=(/  s_dim_3d(1,i),  s_dim_3d(2,i),  s_dim_3d(3,i),1        ,1/) ))
   end do ThreeDVars

end if

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
write (*,*)'netCDF file is synched ...'

call end_model_instance(Var)   ! should avoid any memory leaking

! Reset highest_obs_level and highest_obs_height_m to MISSING_R8 so that next
! assimilation can use its own values; it might be a new obs_seq file with new
! input.nml and new highest_obs_pressure_mb.
! It's here because this is the last, easily identifiable, thing done in
! model_mod before filter moves on to the next model advance/assimilation.
! Nah, there's no mechanism to change highest_obs_pressure in the midst of a
! obs_seq file run, and any values defined here will be lost when a new obs_seq
! file is processed (meaning a new filter_assim is run).
! highest_obs_level    = MISSING_R8
! highest_obs_height_m = MISSING_R8


contains
   ! Internal subroutine - checks error status after each netcdf, prints 
   !                       text message each time an error code is returned. 
   subroutine check(istatus)
   integer, intent ( in) :: istatus
   if (istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_vars', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end function nc_write_model_vars


! End of module I/O

!#######################################################################

! model_interpolate section


   subroutine model_interpolate(x, location, obs_type, interp_val, istatus)
!=======================================================================
!

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
integer,            intent(out) :: istatus
real(r8),           intent(out) :: interp_val

integer  :: i, vstatus
real(r8) :: bot_lon, top_lon, delta_lon,                                       &
            lon_below, lon_above, lat_below, lat_above, lev_below, lev_above,  &
            lon_fract, lat_fract, vals(2, 2), temp_lon, a(2),                  &
            lon_lat_lev(3), level, pressure, height

integer  :: s_type, s_type_01d,s_type_2d,s_type_3d,   &
            lon_ind_below, lon_ind_above, lat_ind_below, lat_ind_above, &
            num_lons, num_lats
character (len=8) :: lon_name, lat_name, lev_name

! Would it be better to pass state as prog_var_type (model state type) to here?
! As opposed to the stripped state vector. YES. This would give time interp.
! capability here; do we really want this or should it be pushed up?

! istatus   meaning                  return expected obs?   assimilate?
! 0         obs and model are fine;  yes                    yes
! 1         fatal problem;           no                     no
! 2         exclude valid obs        yes                    no
! 3         unfamiliar obs type      no                     no

! These are fields which were observed, and will have 3d locations, but the 
! corresponding state-vector component could, concievably, be missing one of the dimensions.
! The only use for such fields I have thought of is parameterization tuning.
! Such fields would not have observations associated with them.
! for now I will assume that observed fields are not missing any dimensions.
! PS is missing lev, although we have never assimilated those obs.

! Start with no errors in 
istatus = 0
vstatus = 0
vals = 0.0_r8

! Always fill the ps arrays with the state vector here, since most obs and vertical locations
!    will need that info.  "Always" allows ens_mean_for_model to set ps arrays once for all
!    of the get_close_obs calls, without having to remove the ps arrays contents at the end
!    of the get_close_obs calls, which is hard to identify.
call set_ps_arrays(x)

! Get the observation (horizontal) position, in degrees 
lon_lat_lev = get_location(location)

! check whether model_mod can interpolate the requested variable
! Pressure (3d) can't be specified as a state vector field (so s_type will = MISSING_I), 
! but can be calculated for CAM, so obs_type = KIND_PRESSURE is acceptable.
s_type = dart_to_cam_kinds(obs_type) 

if (s_type == MISSING_I .and. obs_type .ne. KIND_PRESSURE) then
   istatus = 3
   interp_val = 0._r8
! check
   write(*,*) 'Wrong type of obs = ', obs_type
   return
end if

! Get lon and lat grid specs
! There can't be any 0d or 1d ob fields, so lump them together for elimination in this search.
s_type_01d = state_num_0d + state_num_1d
! Positions within the rank 2 and 3 fields
s_type_2d = s_type - s_type_01d
s_type_3d = s_type_2d - state_num_2d 
if (s_type <= state_num_0d + state_num_1d) then
   ! error; can't deal with have observed variables that are 0 or 1D in model_mod.
elseif (s_type_2d > 0 .and. s_type_2d <= state_num_2d) then
   lon_name = dim_names(s_dimid_2d(1,s_type_2d))
   lat_name = dim_names(s_dimid_2d(2,s_type_2d))
   lev_name = 'none    '
elseif (s_type_3d > 0 .and. s_type_3d <= state_num_3d) then
   lon_name = dim_names(s_dimid_3d(2,s_type_3d))
   lat_name = dim_names(s_dimid_3d(3,s_type_3d))
   lev_name = dim_names(s_dimid_3d(1,s_type_3d))
end if

!------------------------------------------------------------------------------
! Gack!
! staggered longitudes; slon (4x5 fv grid) = [-2.5, 2.5,...,352.5]  !
!                        lon ( "         ) = [    0.,  5.,...,  355.]
!! This is a problem for lon = 359, for example.  It's not in the range of slon. 
!------------------------------------------------------------------------------

! Compute bracketing lon indices
! Define a local longitude to deal with CAM's wierd staggered longitude grid.
temp_lon = lon_lat_lev(1)

if (lon_name == 'lon     ') then
   num_lons  = lon%length
   bot_lon   = lon%vals(1)
   top_lon   = lon%vals(num_lons)
   delta_lon = lon%vals(2) - lon%vals(1)
elseif (lon_name == 'slon    ') then
   num_lons  = slon%length
   bot_lon   = slon%vals(1)
   top_lon   = slon%vals(num_lons)
   delta_lon = slon%vals(2) - slon%vals(1)
   ! Make certain longitudes conform to the wierd CAM staggered grid.
   if ((lon_lat_lev(1) - top_lon) >= delta_lon) temp_lon = lon_lat_lev(1) - 360._r8
end if

if (temp_lon >= bot_lon .and. temp_lon <= top_lon) then
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
if (vert_is_level(location)) then
   ! Case 1: model level specified in vertical
   ! Pobs
   level = lon_lat_lev(3)
      call get_val_level                  &
      (vals(1, 1), x, lon_ind_below, lat_ind_below, nint(level), obs_type, vstatus)
   if (vstatus /= 1) call get_val_level   &
      (vals(1, 2), x, lon_ind_below, lat_ind_above, nint(level), obs_type, vstatus)
   if (vstatus /= 1) call get_val_level   &
      (vals(2, 1), x, lon_ind_above, lat_ind_below, nint(level), obs_type, vstatus)
   if (vstatus /= 1) call get_val_level   &
      (vals(2, 2), x, lon_ind_above, lat_ind_above, nint(level), obs_type, vstatus)
   ! Pobs end

elseif (vert_is_pressure(location)) then
   ! which_vert is pressure for this obs
   pressure = lon_lat_lev(3)
      call get_val_pressure                  &
      (vals(1,1),x,lon_ind_below,lat_ind_below,pressure,obs_type,vstatus)
   if (vstatus /= 1) call get_val_pressure   &
      (vals(1,2),x,lon_ind_below,lat_ind_above,pressure,obs_type,vstatus)
   if (vstatus /= 1) call get_val_pressure   &
      (vals(2,1),x,lon_ind_above,lat_ind_below,pressure,obs_type,vstatus)
   if (vstatus /= 1) call get_val_pressure   &
      (vals(2,2),x,lon_ind_above,lat_ind_above,pressure,obs_type,vstatus)

elseif (vert_is_height(location)) then
   ! which_vert is height for this obs
   height = lon_lat_lev(3)
      call get_val_height                  &
      (vals(1, 1), x, lon_ind_below, lat_ind_below, height, obs_type, vstatus)
   if (vstatus /= 1) call get_val_height   &
      (vals(1, 2), x, lon_ind_below, lat_ind_above, height, obs_type, vstatus)
   if (vstatus /= 1) call get_val_height   &
      (vals(2, 1), x, lon_ind_above, lat_ind_below, height, obs_type, vstatus)
   if (vstatus /= 1) call get_val_height   &
      (vals(2, 2), x, lon_ind_above, lat_ind_above, height, obs_type, vstatus)

elseif (vert_is_surface(location)) then
   ! location_mod:interactive_location asks for surface obs to have vertical coord = ps(hPa)
   ! The 'lev' argument is set to 1 because there is no level for these types, and 'lev' will be
   ! ignored.
                     call get_val(vals(1,1),x, lon_ind_below, lat_ind_below, 1, obs_type, vstatus)
   if (vstatus /= 1) call get_val(vals(1,2),x, lon_ind_below, lat_ind_above, 1, obs_type, vstatus)
   if (vstatus /= 1) call get_val(vals(2,1),x, lon_ind_above, lat_ind_below, 1, obs_type, vstatus)
   if (vstatus /= 1) call get_val(vals(2,2),x, lon_ind_above, lat_ind_above, 1, obs_type, vstatus)

end if

! lat is already converted to degrees by get_location
if (abs(lon_lat_lev(2)) > max_obs_lat_degree .and. vstatus /= 1) then
   istatus = 2
else
   istatus = vstatus
end if

! indices of vals are (longitude, latitude)
if (istatus /= 1) then
   do i = 1, 2
      a(i) = lon_fract * vals(2, i) + (1.0_r8 - lon_fract) * vals(1, i)
   end do
   interp_val = lat_fract * a(2) + (1.0_r8 - lat_fract) * a(1)
else
   ! ? Should this return MISSING_R8?
   interp_val = 0.0_r8
end if

! Set the element of ps that's tested elsewhere back to MISSING_R8, to signal
! other routines to calculate the ps arrays for themselves
! Currently (10/26/06) this flag is not used.
! ps(1,1) = MISSING_R8

end subroutine model_interpolate


! Pobs
   subroutine get_val_level(val, x, lon_index, lat_index, level, obs_kind, istatus)
!=======================================================================
!   subroutine get_val_level(val, x, lon_index, lat_index, level, obs_kind, istatus)
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
real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: lon_index, lat_index, level, obs_kind
integer,  intent(out) :: istatus

integer               :: top_lev, bot_lev, i, vstatus, num_levs
real(r8)              :: bot_val, top_val, frac, p_surf, threshold
real(r8), allocatable :: pfull(:)

! No errors to start with
istatus = 0
vstatus = 0

! This assumes that all variables are defined on model levels, not on interface levels.
num_levs = dim_sizes(find_name('lev     ',dim_names))

! Interpolate in vertical to get two bounding levels
if (level > num_levs .or. level < 1) then
   ! Exclude obs below the model's lowest level and above the highest level
   istatus = 1
   ! Should this be MISSING_r8?
   val = 0._r8
else 
   if (highest_obs_level == MISSING_R8) then
      allocate(pfull(num_levs))
      ! To do this completely right; p_surf would depend on whether obs_kind was on a staggered grid
      ! and highest_obs_level would be recalculated for each location passed in.
      ! Since the hybrid coord system is pure pressure at the top levels, I'll ignore these for now.
      call get_val(p_surf, x, lon_index, lat_index, -1, KIND_SURFACE_PRESSURE, vstatus)
      call plevs_cam (p_surf, num_levs, pfull)
      highest_obs_level = 1
      threshold = highest_obs_pressure_mb*100.0_r8
      do while ((pfull(nint(highest_obs_level))) < threshold )
         highest_obs_level = highest_obs_level +1
      end do
      deallocate(pfull)
   end if

   if (level < highest_obs_level) then
      ! Exclude from assimilation the obs above a user specified level
      ! but still calculate the expected obs.
      istatus = 2
   else
      istatus = 0
   end if

   if (obs_kind == KIND_PRESSURE) then
      ! Can't get the value from get_val because 3d pressure is not a model variable.
      ! Can calculate it from ps.

      ! ps is on A-grid, so no need to check for staggered grids
      call get_val(p_surf, x, lon_index, lat_index, -1, KIND_SURFACE_PRESSURE, vstatus)
      if (vstatus > 0) then
         val = 0._r8
         istatus = 1
         return
      end if
      ! Next, get the values on the levels for this ps
      allocate (pfull(num_levs))
      call plevs_cam (p_surf, num_levs, pfull)

      val = pfull(level)
      deallocate(pfull)
   else 
       call get_val(val, x, lon_index, lat_index, level, obs_kind, vstatus)
   end if

   if (vstatus /= 0) then
      istatus = 1
      val = 0._r8
   end if
end if

end subroutine get_val_level
! Pobs end


   subroutine get_val_pressure(val, x, lon_index, lat_index, pressure, obs_kind, istatus)
!=======================================================================
!
! Gets the vertically interpolated value on pressure for variable obs_kind
! at lon_index, lat_index horizontal grid point
!
! This version excludes observations below lowest level pressure and above
! highest level pressure.


real(r8), intent(out) :: val
real(r8), intent(in)  :: x(:), pressure
integer,  intent(in)  :: lon_index, lat_index, obs_kind 
integer,  intent(out) :: istatus

real(r8)              :: bot_val, top_val, p_surf, frac
real(r8), allocatable :: pfull(:)
integer               :: top_lev, bot_lev, i, vstatus, num_levs
character(len=129)    :: errstring

! No errors to start with
istatus = 0
vstatus = 0

! Need to get the surface pressure at this point. 
! Find out whether the observed field is a staggered in CAM.
! Check whether the state vector has wind components on staggered grids, i.e. whether CAM is FV.
! find_name returns 0 if the field name is not found in the cflds list.
! Add more staggered variables later?
! Can I make a more generic test; loop over all KIND_s, then check whether any of the 
!    associated dimensions are staggered?   Sounds too expensive to be worth it. . .?
if     (obs_kind == KIND_U_WIND_COMPONENT .and. find_name('US      ', cflds) /= 0) then 
   ! ps defined on lat grid (-90...90, nlat = nslat + 1), 
   !    need it on staggered lat grid, which starts half a grid spacing north
   ! What about poles?  
   !    It's OK; ps is defined on the 'lat' grid, which runs [-90...90], 
   !    while staggered US is on the 'slat' grid, defined only *inside* this range.
   p_surf = ps_stagr_lat(lon_index, lat_index)
elseif (obs_kind == KIND_V_WIND_COMPONENT .and. find_name('VS      ', cflds) /= 0) then
   ! lon =     0...     255 (for 5 degree grid)
   !slon = -2.5 ... 252.5
   p_surf = ps_stagr_lon(lon_index, lat_index)
else
   ! A-grid ps can be retrieved from state vector, which was used to define ps on entry to 
   ! model_interpolate.
   p_surf = ps(lon_index, lat_index)
end if
!   if (vstatus > 0) then
!      val = 0.0_r8
!      istatus = 1
!      return
!   end if

! Next, get the pressures on the levels for this ps
! Assuming we'll only need pressures on model mid-point levels, not interface levels.
! This pressure column will be for the correct grid for obs_kind, since p_surf was taken
!     from the grid-correct ps[_stagr] grid
num_levs = dim_sizes(find_name('lev     ',dim_names))
allocate (pfull(num_levs))
call plevs_cam (p_surf, num_levs, pfull)

if (pressure <= pfull(1) .or. pressure >= pfull(num_levs)) then
   ! Exclude obs below the model's lowest level and above the highest level
   ! We *could* possibly use ps and p(num_levs) to interpolate for points below the lowest level.
   istatus = 1
   val = 0._r8
else 
   ! Interpolate in vertical to get two bounding levels
   if (pressure < highest_obs_pressure_mb * 100.0_r8) then
      ! Exclude from assimilation the obs above a user specified level
      ! but still calculate the expected obs.
      istatus = 2
   else
      istatus = 0
   end if
   ! Search down through pressures
   do i = 2, num_levs 
      if (pressure < pfull(i)) then
         top_lev = i -1
         bot_lev = i
         frac = (pfull(i) - pressure) / &
                (pfull(i) - pfull(i - 1))
         goto 21
      end if
   end do

   21 continue

   ! Pobs
   if (obs_kind == KIND_PRESSURE) then
      ! can't get pressure on levels from state vector; get it from pfull instead
      ! get_val_pressure is called for 4 different columns, which will have different pfulls 
      ! for each ps is on A-grid, so no need to check for staggered grids
      bot_val = pfull(bot_lev)
      top_val = pfull(top_lev)
      val = (1.0_r8 - frac) * bot_val + frac * top_val
!      if (abs((val - pressure)/val) > 1.0E-12) then
!         ! We're looking for a pressure on a model level, which is exactly what pfull provides,
!! NOT HERE; that happens in get_val_level
!         write(errstring, *) 'val /= pressure = ',val,pressure,' when val is a P obs '
!         call error_handler(E_WARN, 'get_val_pressure', errstring, source, revision, revdate)
!      end if
   else 
   ! Pobs end
                        call get_val(bot_val, x, lon_index, lat_index, bot_lev, obs_kind, vstatus)
      if (vstatus == 0) call get_val(top_val, x, lon_index, lat_index, top_lev, obs_kind, vstatus)
      if (vstatus == 0) then
         val = (1.0_r8 - frac) * bot_val + frac * top_val
      else
         istatus = 1
         val = 0._r8
      end if
   end if
   ! Pobs

end if

deallocate (pfull)

end subroutine get_val_pressure



   subroutine get_val_height(val, vec, lon_index, lat_index, height, obs_kind, istatus)
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

real(r8), intent(out) :: val
real(r8), intent(in)  :: vec(:), height
integer,  intent(in)  :: lon_index, lat_index, obs_kind 
integer,  intent(out) :: istatus

real(r8), allocatable :: pfull(:), model_h(:)
integer  :: top_lev, bot_lev, i, vstatus, num_levs
real(r8) :: bot_val, top_val, frac
real(r8) :: p_surf
logical  :: stagr_lon, stagr_lat

! No errors to start with
istatus = 0
vstatus = 0
stagr_lon = .false.
stagr_lat = .false.

! Need to get the surface pressure at this point for dcz2. 
! Assuming we'll only need pressures on model mid-point levels, not interface levels.
num_levs = dim_sizes(find_name('lev     ',dim_names))
allocate (model_h(num_levs))

! Need to get the surface pressure at this point. 
! Check whether the state vector has wind components on staggered grids, i.e. whether CAM is FV.
! find_name returns 0 is the field name is not found in the cflds list.
! See get_val_press for more documentation.
if     (obs_kind == KIND_U_WIND_COMPONENT .and. find_name('US      ', cflds) /= 0) then  
   p_surf = ps_stagr_lat(lon_index, lat_index)
   stagr_lat = .true.
elseif (obs_kind == KIND_V_WIND_COMPONENT .and. find_name('VS      ', cflds) /= 0) then
   p_surf = ps_stagr_lon(lon_index, lat_index)
   stagr_lon = .true.
else
   p_surf = ps(lon_index, lat_index)
end if

! Next, get the heights on the levels for this ps

! merge/MPI
! We want to use the new vec for each new ob on height because the state was updated 
! for all previous obs, and we want to use the most up to date state to get the best location.
call model_heights(vec, p_surf, lon_index, lat_index, num_levs, stagr_lon, stagr_lat, & 
                   model_h, vstatus)

! check

! Interpolate in vertical to get two bounding levels
if (height >= model_h(1) .or. height <= model_h(num_levs)) then
   ! Exclude obs below the model's lowest level and above the highest level
   istatus = 1
   val = 0._r8
! check
else 
   if (highest_obs_height_m == MISSING_R8) then
      allocate(pfull(num_levs))
      call plevs_cam (p_surf, num_levs, pfull)
      do i=1,num_levs
         if (pfull(i) > highest_obs_pressure_mb*100.0_r8) then
            highest_obs_height_m = model_h(i)
            go to 10
         end if
      end do
      deallocate(pfull)
   end if

10 if (height > highest_obs_height_m ) then
      ! Exclude from assimilation the obs above a user specified level
      ! but still calculate the expected obs.
      istatus = 2
! check
      if (do_out) &
      write(*,*) 'get_val_height; height, highest_obs_height_m = ',height, highest_obs_height_m
   else
      istatus = 0
   end if

   ! Search down through heights
   do i = 2, num_levs 
      if (height > model_h(i)) then
         top_lev = i -1
         bot_lev = i
         frac = (model_h(i) - height      ) / &
                (model_h(i) - model_h(i-1))
! check
         goto 21
      end if
   end do

   21 continue


   ! Pobs
   if (obs_kind == KIND_PRESSURE) then
      ! Observing a pressure on a height surface sounds silly.  But for completeness:
      ! get_val_height is called for 4 different columns, which will have different pfulls for each.

      ! Next, get the values on the levels for this ps
      ! ps is on A-grid, so no need to check for staggered grids
      allocate (pfull(num_levs))
      call plevs_cam (p_surf, num_levs, pfull)

      bot_val = pfull(bot_lev)
      top_val = pfull(top_lev)
      deallocate (pfull)
   else 
   ! Pobs end
                        call get_val(bot_val, vec, lon_index, lat_index, bot_lev, obs_kind, vstatus)
      if (vstatus == 0) call get_val(top_val, vec, lon_index, lat_index, top_lev, obs_kind, vstatus)
   ! Pobs
   end if

   if (vstatus == 0) then
      val = (1.0_r8 - frac) * bot_val + frac * top_val
   else
      istatus = 1
      val = 0._r8
   end if
end if

deallocate (model_h)

end subroutine get_val_height



   subroutine get_val(val, x, lon_index, lat_index, level, obs_kind, istatus)
!=======================================================================
! function get_val(x, lon_index, lat_index, level, obs_kind, istatus)
!

real(r8), intent(out) :: val
real(r8), intent(in) :: x(:)
integer, intent(in) :: lon_index, lat_index, level, obs_kind
integer, intent(out) :: istatus

integer :: indx, field_type
character(len=129) :: errstring

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

val = x(indx)

return

end subroutine get_val



! End of model_interpolate section

!#######################################################################

! vector-field translations


   subroutine prog_var_to_vector(var, x)
!=======================================================================
! subroutine prog_var_to_vector(var, x)
!

type(model_type), intent(in)  :: var
real(r8),         intent(out) :: x(:)

integer :: i, j, k, nf, indx
character(len=129) :: errstring

! Load components of state vector, starting with scalars (0D) and finishing with 3D
! A whole field will be loaded (by columns for 3D) before the next field is started.
! This is completely different than the B-grid organization, which loaded all the fields
! at a point before moving on to the next point.  The motivations for this change are:
! 1) This easily allows fields with the same rank, but different sizes to be loaded into
!    the vector (i.e. Ustaggered  and T in the cam-fv)
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
   x(indx) = var%vars_0d(nf)
end do

!  1d variables
do nf = 1, state_num_1d
   do i=1,s_dim_1d(nf)
      indx = indx + 1
      x(indx) = var%vars_1d(i, nf)
   end do
end do

!  Surface pressure and other 2d ; load by first dimension; 
!      lon and/or lat and/or lev, in that precedence, as loaded in to vars_2d in trans_coord
do nf = 1, state_num_2d
   do j=1,s_dim_2d(2,nf)
   do i=1,s_dim_2d(1,nf)  !levs or lons  both reads and writes will be contiguous in this case
      indx = indx + 1
      x(indx) = var%vars_2d(i, j, nf)
   end do
   end do
end do

!  3D fields, loaded by columns (note the coordinate order).
!  This is also looping over latitude values in the outer loop,
!  so if there is some spatial searching; all the pole points will be closer together
!  in memory than in the previous/Bgrid structure.
do nf= 1, state_num_3d

   if (do_out) then
      write(errstring, '(A,4I5)') 'fld, nlons, nlats, nlevs ',nf &
                          ,s_dim_3d(2,nf),s_dim_3d(3,nf),s_dim_3d(1,nf)
      call error_handler(E_MSG, 'prog_var_to_vector', errstring, source, revision, revdate)
   endif

   do i=1,s_dim_3d(3,nf)   !lats
   do j=1,s_dim_3d(2,nf)   !lons
   do k=1,s_dim_3d(1,nf)   !levs  both reads and writes will be contiguous in this case
      indx = indx + 1
      x(indx) = var%vars_3d(k,j,i, nf)
   end do
   end do
   end do
end do

! Temporary check
if (indx /= model_size) then
   write(errstring, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'prog_var_to_vector', errstring, source, revision, revdate)
end if

end subroutine prog_var_to_vector




   subroutine vector_to_prog_var(x, var) 
!=======================================================================
! subroutine vector_to_prog_var(x, var) 
!

real(r8),         intent(in)  :: x(:)
type(model_type), intent(out) :: var

integer :: i, j, k, nf, indx
character(len=129) :: errstring

! Start copying fields from straight vector
indx = 0

! 0d arrays
do nf = 1, state_num_0d
   indx = indx + 1
   var%vars_0d(nf) = x(indx)
end do

!  1d arrays
do nf = 1, state_num_1d
   do i=1,s_dim_1d(nf)
      indx = indx + 1
      var%vars_1d(i, nf) = x(indx)
   end do
end do

!  Surface pressure and other 2d fields 
do nf = 1, state_num_2d
   do j = 1, s_dim_2d(2,nf)
   do i = 1, s_dim_2d(1,nf)
      indx = indx + 1
      var%vars_2d(i, j, nf) = x(indx)
   end do
   end do
end do

! 3D fields; see comments in prog_var_to_vect
do nf = 1, state_num_3d
   if (do_out) then
      write(errstring, '(A,4I5)') 'fld, nlons, nlats, nlevs ',nf &
                       ,s_dim_3d(2,nf),s_dim_3d(3,nf),s_dim_3d(1,nf)
      call error_handler(E_MSG, 'vector_to_prog_var', errstring, source, revision, revdate)
   end if
   do i = 1, s_dim_3d(3,nf)
   do j = 1, s_dim_3d(2,nf)
   do k = 1, s_dim_3d(1,nf)
      indx = indx + 1
      var%vars_3d(k,j,i, nf) = x(indx)
   end do 
   end do
   end do
end do

! Temporary check
if (indx /= model_size) then
   write(errstring, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'vector_to_prog_var', errstring, source, revision, revdate)
end if

end subroutine vector_to_prog_var



! End of vector-field translations

!#######################################################################

! get_close_obs section

   subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                            num_close, close_ind, dist)
!----------------------------------------------------------------------------

!!!ADD IN SOMETHING TO USE EFFICIENTLY IF IT"S AT SAME LOCATION AS PREVIOUS OB!!!

! The kinds are available to do more sophisticated distance computations if needed

implicit none


type(get_close_type), intent(in)  :: gc
type(location_type),  intent(in)  :: base_obs_loc, obs_loc(:)
integer,              intent(in)  :: base_obs_kind, obs_kind(:)
integer,              intent(out) :: num_close, close_ind(:)
real(r8),             intent(out) :: dist(:)

! remove some (unused) variables?
integer                :: i, k, t_ind
integer                :: base_which, local_base_which, obs_which, local_obs_which
real(r8), dimension(3) :: base_array, local_base_array, obs_array, local_obs_array
real(r8)               :: increment, threshold, thresh_wght, signed_magnifier
type(location_type)    :: local_base_obs_loc, local_obs_loc

! If base_obs vert type is not pressure; convert it to pressure
base_which = nint(query_location(base_obs_loc))
if (base_which == VERTISPRESSURE) then
   local_base_obs_loc = base_obs_loc
   local_base_array   = get_location(base_obs_loc)  ! needed in num_close loop
   local_base_which   = base_which
else
   base_array = get_location(base_obs_loc)
   call convert_vert(base_array, base_which, local_base_array, local_base_which, base_obs_kind)
   local_base_obs_loc = set_location(base_array(1), base_array(2), local_base_array(3), &
                                     local_base_which)
end if

! Get all the potentially close obs but no dist (optional argument dist(:) is not present)
call loc_get_close_obs(gc, local_base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                       num_close, close_ind)

threshold = highest_state_pressure_mb *100._r8
thresh_wght = 1._r8/(threshold * threshold)

do k = 1, num_close

   t_ind = close_ind(k)
   obs_array = get_location(obs_loc(t_ind))
   obs_which = nint(query_location(obs_loc(t_ind)))

   if (obs_which == VERTISPRESSURE) then
      ! put the vertical (pressure) of the state/ob in local storage
      local_obs_array(3) = obs_array(3)
      local_obs_which    = obs_which
   else
      ! convert vertical coordinate of obs_loc to pressure.
      call convert_vert(obs_array, obs_which, local_obs_array, local_obs_which, obs_kind(t_ind))

      ! obs_which = -2 (VERTISUNDEF) mean this ob is vertically close to base_obs, no matter what.
      if (local_obs_array(3) == MISSING_R8) then
         local_obs_array(3) = local_base_array(3)
         local_obs_which = local_base_which
      end if

   end if

   local_obs_loc = set_location(obs_array(1), obs_array(2), local_obs_array(3), &
                                   local_obs_which)

   dist(k) = get_dist(local_base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))

   ! Damp the influence of obs below the namelist variable highest_state_pressure_mb on
   ! variables above that level.  
   ! This could also change the distance based on the KIND_s of the base_obs and obs.

   ! dist = 0 for some for synthetic obs.
   ! Multiplicative increase of distance doesn't work for those cases.
   ! I should use additive increase, based on height above threshold.

   ! Push the vertical location of the ob farther from the base_ob.
   ! The 'sign' function gives the magnitude of A times the sign of B.
   ! So it gives the direction, and the magnification of the additive part.
   ! Formula; new_vert_loc = 
   !          old (+/-) (threshold/old) (threshold - old)
   !                     magnifier       something to magnify (in case of co-location)
   ! WATCH OUT: pushing too high will yield press < 0. !
   !    That's OK; get_dist uses abs(loc1%vloc - loc2%vloc), so dist will be > 0.
                                   
! This was much too abrupt; increment is in Pa, but it's being added like radians.
!   increment = threshold - local_obs_array(3)
!   if (increment > 0) then
!      signed_magnifier = sign((threshold/local_obs_array(3)),            &
!                              (local_obs_array(3)-local_base_array(3)))
!      dist(k) = dist(k) + increment * signed_magnifier 
!!      local_obs_array(3) = local_obs_array(3) + increment * signed_magnifier 
!   end if

   increment = threshold - local_obs_array(3)
   if (increment > 0) then
!     This output can be seen in Pre-J/DART/models/cam/work/filter_Tprof_damp_out.out;
!     it all looks fine.
!      write(*,'(A,3(1x,F10.2)))') 'damping state at location ',  &
!           obs_array(1), obs_array(2),local_obs_array(3)
! too sharp      dist(k) = dist(k) + increment / threshold
      dist(k) = dist(k) + increment * increment * thresh_wght
   end if

end do

end subroutine get_close_obs


!=======================================================================
   subroutine convert_vert (old_array, old_which, new_array, new_which, dart_kind)
!=======================================================================
! subroutine convert_vert(old_loc, new_loc, dart_kind)
!
! Uses model information and subroutines to convert the vertical location of an ob 
! (prior, model state variable, or actual ob) into the standard vertical coordinate (pressure).
! Called by model_mod:get_close_obs.
! Kevin Raeder 10/26/2006

integer,                intent(in)    :: dart_kind, old_which
integer,                intent(out)   :: new_which
real(r8), dimension(3), intent(in)    :: old_array
real(r8), dimension(3), intent(inout) :: new_array 

integer   :: i, dim_index, num_levs, top_lev, bot_lev
integer   :: lon_which_dimid=0, lat_which_dimid=0, lon_index=0, lat_index=0
integer   :: rank_kind, cam_kind, istatus
real(r8)  :: p_surf,   phi_surf,  frac
logical   :: stagr_lon = .false., stagr_lat = .false.
real(r8), allocatable :: p_col(:), model_h(:) 
type(location_type)   :: dum_loc

character(len=129)    :: errstring
character(len=8)      :: dim_name

new_array = MISSING_R8

if (old_which == VERTISPRESSURE .or. old_which == VERTISHEIGHT  .or. &
    old_which == VERTISLEVEL    .or. old_which == VERTISSURFACE .or. &
    old_which == VERTISUNDEF   ) then
!  proceed
else
   write(errstring,'(''skipping obs at '',3(F9.5,1x),I2,'' vertical problem'')') &
                   old_array, old_which
   call error_handler(E_MSG, 'convert_vert', errstring,source,revision,revdate)
   return
end if

! Find the nfld of this dart-KIND
if (dart_kind > 0) then
   ! non-identity obs
   cam_kind = dart_to_cam_kinds(dart_kind)
else if (dart_kind < 0) then
   ! identity obs; dart_kind = -1*state_vector_index
   ! Value returned in cam_kind will be the nfld value of this field, not the usual dart_kind.
   call get_state_meta_data(dart_kind, dum_loc, cam_kind)
end if
! Find the index of this kind within its group of same-rank fields
rank_kind = cam_kind

! Figure out what rank CAM field this corresponds to, so that vertical coordinate can be determined
! Also need lon and lat indices to select ps for calc of p_col for vertical conversion.
if (rank_kind <= state_num_0d) then
   call coord_index('lon     ', old_array(1), lon_index)
   call coord_index('lat     ', old_array(2), lat_index)
   go to 10
else
   rank_kind = rank_kind - state_num_0d
end if

if (rank_kind <= state_num_1d) then
   dim_name = dim_names(s_dimid_1d(rank_kind))
   if (dim_name .eq.'lon     ' .or. dim_name .eq.'slon    ' ) then   
! s_dimid_1d holds the single CAM dimension ids of the dimensions of the 1D fields
      lon_which_dimid = 1                             
      call coord_index(dim_name, old_array(1), lon_index)
   elseif (dim_name .eq.'lat     ' .or. dim_name .eq.'slat    ' ) then
      lat_which_dimid = 1
      call coord_index(dim_name, old_array(2), lat_index)

! This may be premature; we haven not converted the 3rd dimension yet!
! This may be spurious; we may not use lev_which_dimid.
! Similarly with other ranks
! Also; coordinate 'lev' is filled with 1000*(A+B), not levels 1,2,...

!   elseif (dim_name .eq.'lev     ' .or. dim_name .eq.'ilev    ' ) then
!      lev_which_dimid = 1
!      call coord_index(dim_name, old_array(3), lev_index)
   end if
   go to 10
else
   rank_kind = rank_kind - state_num_1d
end if
   
if (rank_kind <= state_num_2d) then
   do i=1,2
      dim_name = dim_names(s_dimid_2d(i,rank_kind))
      if (dim_name .eq.'lon     ' .or. dim_name .eq.'slon    ' ) then
         lon_which_dimid = s_dimid_2d(i,rank_kind)   ! assign the CAM longitude dimension, if present
         call coord_index(dim_name, old_array(1), lon_index)
      elseif (dim_name .eq.'lat     ' .or. dim_name .eq.'slat    ' ) then   
         lat_which_dimid = s_dimid_2d(i,rank_kind)
         call coord_index(dim_name, old_array(2), lat_index)
!      elseif (dim_name .eq.'lev     ' .or. dim_name .eq.'ilev    ' ) then
!         lev_which_dimid = s_dimid_2d(i,rank_kind)
!         call coord_index(dim_name, old_array(3), lev_index)
      end if
   end do
   go to 10
else
   rank_kind = rank_kind - state_num_2d
end if

if (rank_kind <= state_num_3d) then
   do i=1,3
      dim_name = dim_names(s_dimid_3d(i,rank_kind))
      if (dim_name .eq.'lon     ' .or. dim_name .eq.'slon    ' ) then
         lon_which_dimid = s_dimid_3d(i,rank_kind)
         call coord_index(dim_name, old_array(1), lon_index)
      elseif (dim_name .eq.'lat     ' .or. dim_name .eq.'slat    ' ) then
         lat_which_dimid = s_dimid_3d(i,rank_kind)
         call coord_index(dim_name, old_array(2), lat_index)
!      elseif (dim_name .eq.'lev     ' .or. dim_name .eq.'ilev    ' ) then
!         lev_which_dimid = s_dimid_3d(i,rank_kind)
!         call coord_index(dim_name, old_array(3), lev_index)
      end if
   end do

   go to 10
else
   ! print error; field not found for dart_kind = ...
end if

10 continue


! doubly staggered not handled here
if     (lat_which_dimid == find_name('slat    ', dim_names)) then 
   stagr_lat = .true.
   p_surf = ps_stagr_lat(lon_index, lat_index)
elseif (lon_which_dimid == find_name('slon    ', dim_names)) then
   stagr_lon = .true.
   p_surf = ps_stagr_lon(lon_index, lat_index)
elseif (lon_which_dimid == 0 .or. lat_which_dimid == 0) then
   ! one of these dimensions is missing from this variable
   p_surf = P0%vals(1)
   ! Or should this be the average of ps over the undefined dimension?
   ! Or is this meaningless?
else
   p_surf = ps(lon_index, lat_index)
end if

! Need the vertical pressure structure for this column
! This routine will be called for :
!   model grid points (from get_close_obs) (just one column of the state vector is the correct one),
!   expected obs (4 times from model_interpolate) (the correct column is an interpolation of the 
!      surrounding 4 columns).  

! Convert vertical coordinate from one of the following to pressure.
! integer, parameter :: VERTISUNDEF    = -2 ! has no vertical location (undefined)
! integer, parameter :: VERTISSURFACE  = -1 ! surface value
! integer, parameter :: VERTISLEVEL    =  1 ! by level
! integer, parameter :: VERTISPRESSURE =  2 ! by pressure
! integer, parameter :: VERTISHEIGHT   =  3 ! by height

if (old_which == VERTISUNDEF) then
   ! Field with no vertical location; get_dist will calculate horiz dist only unless this case
   ! is handled by the calling routine.

   ! If a parameter/state variable is supposed to be close to everything, 
   ! then I would need to have the/an other location to set it to,
   ! Send back new_array empty and test for that in the calling routine, 
   ! where the other location exists.
   ! For model variables user specifies which_vert for each state field, 
   ! so when user specifies undefined, then this should return;
   new_array(3) = MISSING_R8
   new_which    = old_which
elseif (old_which == VERTISSURFACE ) then       
   ! surface field; change which_vert for the distance calculation
   new_array(3) =  p_surf
   new_which    = 2
elseif (old_which == VERTISLEVEL ) then
   ! Next get the values on the levels for this ps
   num_levs = dim_sizes(find_name('lev     ',dim_names))
   allocate (p_col(num_levs))
   call plevs_cam (p_surf, num_levs, p_col)

   ! OR do this for all columns in static_init_model_dist, which would make PS (and P) globally 
   ! available for all regions?
   new_array(3) = p_col(nint(old_array(3)))
   new_which    = 2
   deallocate (p_col)
elseif (old_which == VERTISHEIGHT) then
   num_levs = dim_sizes(find_name('lev     ',dim_names))
   allocate (model_h(num_levs), p_col(num_levs))
   call plevs_cam (p_surf, num_levs, p_col)
   ! ens_mean is global storage that should have been filled 
   ! by a call from filter_assim to ens_mean_for_model.
   call model_heights(ens_mean, p_surf, lon_index, lat_index, num_levs, stagr_lon, stagr_lat,  &
                      model_h, istatus)
   ! Search down through heights
   ! This assumes linear relationship of pressure to height over each model layer, 
   ! when really it's exponential.  How bad is that?
   do i = 2, num_levs 
      if (old_array(3) > model_h(i)) then
         top_lev = i -1
         bot_lev = i
         frac = (old_array(3) - model_h(i)) / (model_h(i-1) - model_h(i))
         goto 21
      end if
   end do

!  write error message if not found within model level heights.

   21 continue

   new_array(3) = (1.0_r8 - frac) * p_col(bot_lev) + frac * p_col(top_lev)
   new_which    = 2

   deallocate (model_h, p_col)
else
   write(errstring, *) 'model which_vert = ',old_which,' not handled in convert_vert '
   call error_handler(E_ERR, 'convert_vert', errstring,source,revision,revdate)
end if

return

end subroutine convert_vert


! End of get_close_obs section

!#######################################################################

! Initial conditions for DART 

   subroutine pert_model_state(state, pert_state, interf_provided)
!       call pert_model_state(ens_mean, temp_ens, interf_provided)
!=======================================================================
! subroutine pert_model_state(state, pert_state, interf_provided)
!
! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

! added to give each ens member a different sequence when perturbing model parameter fields

real(r8), intent(in)    :: state(:)
real(r8), intent(out)   :: pert_state(:)
logical,  intent(out)   :: interf_provided

type(random_seq_type)   :: random_seq
type(model_type)        :: var_temp
integer                 :: i, ipert, j, k, m, field_num, iunit, io
integer                 :: dim1, dim2, dim3
real(r8)                :: pert_val

! FIX for 1D 0D  fields?

! trans_pv_sv_pert0.f90 needs to perturb model parameters for the filter_ics.
! Use the (single) state value as the "ens_mean" here.

interf_provided = .true.

call init_model_instance(var_temp)
call vector_to_prog_var(state,var_temp)

! If first call, then initialize random sequence for perturbations.
! init_random_seq only needed for documentation; initializing the module.
if (first_pert_call) then 
   call init_random_seq(random_seq)
   first_pert_call = .false.
end if

! init_random_seq calls init_ran1, but I need to call init_ran1 with a different seed/temp 
! for each ens_member.  Get a new seed by keeping track of the previous seed.

ens_member = ens_member + 1
call init_ran1(random_seq,-1*ens_member)

ipert = 1
do while (state_names_pert(ipert) /= '        ')
   WRITE(*,*) 'Perturbing ',state_names_pert(ipert)
   do m=1,nflds
      if (state_names_pert(ipert) == cflds(m)) then
         WRITE(*,*) '   Found match  ',cflds(m)
           
         if (m <= state_num_2d + state_num_1d + state_num_0d) then
            field_num = m - state_num_1d - state_num_0d
            dim1 = dim_sizes(s_dimid_2d(1,field_num))
            dim2 = dim_sizes(s_dimid_2d(2,field_num))
            WRITE(*,'(A,1p,2E12.4)') '    first and last state = ', &
                   var_temp%vars_2d(1,1,field_num),var_temp%vars_2d(dim1,dim2,field_num)
            do j = 1, dim2
            do i = 1, dim1
               pert_val = random_gaussian(random_seq, var_temp%vars_2d(i,j,field_num), &
                                          state_names_sd(ipert)) 
               var_temp%vars_2d(i,j,field_num) = pert_val
            end do
            end do
            WRITE(*,'(A,1p,2E12.4)') ' new first and last state = ', &
                   var_temp%vars_2d(1,1,field_num),var_temp%vars_2d(dim1,dim2,field_num)
         else  
            field_num = m - state_num_2d - state_num_1d - state_num_0d
            dim1 = dim_sizes(s_dimid_3d(1,field_num))
            dim2 = dim_sizes(s_dimid_3d(2,field_num))
            dim3 = dim_sizes(s_dimid_3d(3,field_num))

            do j = 1, dim3
            do i = 1, dim2
            do k = 1, dim1
!              pert_val = rand#(O(0-1)) * standard dev  + mean
               pert_val = random_gaussian(random_seq, var_temp%vars_3d(k,i,j,field_num), &
                                          state_names_sd(ipert)) 
               var_temp%vars_3d(k,i,j,field_num) = pert_val
            end do
            end do
            end do
         end if
      end if
   end do
   ipert = ipert + 1
end do
call prog_var_to_vector(var_temp,pert_state)
call end_model_instance(var_temp)

end subroutine pert_model_state


   subroutine init_conditions(x)
!=======================================================================
! subroutine init_conditions(x)
!
! Reads in restart initial conditions  -- noop for CAM

real(r8), intent(inout) :: x(:)

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
integer :: i, j, fld_ind(2), done, indx

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
      ! We've found the desired field; now find indx of lev and/or lon and/or lat 
      do j=1,2
         if (dim_names(s_dimid_2d(j,i)) == 'lon     ' .or. & 
             dim_names(s_dimid_2d(j,i)) == 'slon    '     ) fld_ind(j) = lon_ind
         if (dim_names(s_dimid_2d(j,i)) == 'lat     ' .or. & 
             dim_names(s_dimid_2d(j,i)) == 'slat    '     ) fld_ind(j) = lat_ind
         if (dim_names(s_dimid_2d(j,i)) == 'lev     ' .or. & 
             dim_names(s_dimid_2d(j,i)) == 'ilev    '     ) fld_ind(j) = lev_ind
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
do i=1,state_num_3d
   if (ifld - done == i) then
      ! We've found the desired field; now find indx of lat, lon, lev
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
endif
if (dim_name == 'slat    ') lat_val = slat%vals(indx)
if (dim_name == 'ilev    ') lev_val = ilev%vals(indx)
! Add more for other coords?  hyam...?  Not for now; never referenced indirectly

return

end subroutine coord_val


   subroutine coord_index(dim_name, val, indx, other_indx)
!==========================================================================
! subroutine coord_index(dim_name, indx, val, indx, other_indx)

! Given the name of the coordinate to be searched and the value, 
! Returns the index of the closest coordinate value.  
! Optionally returns the next closest index too, which may be < or > the closest.
! Used by get_state_meta_data.

character (len=8), intent(in)  :: dim_name
real(r8),          intent(in)  :: val
integer,           intent(out) :: indx
integer, optional, intent(out) :: other_indx

real(r8), pointer :: coord(:)
real(r8)          :: diff_upper, diff_lower, val_local 
integer           :: coord_len, i

val_local = val

if (dim_name == 'lon     ') then
   coord     => lon%vals
   coord_len =  lon%length
elseif (dim_name == 'lat     ') then
   coord     => lat%vals
   coord_len =  lat%length
elseif (dim_name == 'lev     ') then
   coord     => lev%vals
   coord_len =  lev%length
elseif (dim_name == 'slon    ') then
   coord     => slon%vals
   coord_len =  slon%length
   ! Make sure longitudes conform to the wierd CAM staggered grid.
   if ((val - coord(coord_len)) >= (coord(coord_len)-coord(coord_len-1))) &
      val_local = val_local - 360._r8
elseif (dim_name == 'slat    ') then
   coord     => slat%vals
   coord_len =  slat%length
elseif (dim_name == 'ilev    ') then
   coord     => ilev%vals
   coord_len =  ilev%length
else
   indx = MISSING_I
!  PRINT error?
end if

! further check?  for blunders check that coord(1) - val is smaller than coord(2) - coord(1), etc.
! Assumes that coordinates are monotonic; not true for hyam, hyai.  But we don't reference them.
! The first 2 if blocks work for latitudes and levels.  Longitudes must be handled in the calling
!    routine.
if (val_local < coord(1)) then
   indx = 1
   if (present(other_indx)) other_indx = 1
   return
elseif (val_local > coord(coord_len)) then
   indx = coord_len
   if (present(other_indx)) other_indx = coord_len
   return
else
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

end subroutine coord_index



   subroutine set_ps_arrays (vec)
!=====================================================================

real(r8), intent(in) :: vec(:)

integer :: ind, m,n, slon_index, slat_index, lon_index, lat_index, &
           fld_index_2d, fld_index

lon_index  = find_name('lon     ',dim_names)
lat_index  = find_name('lat     ',dim_names)
slon_index = find_name('slon    ',dim_names)
slat_index = find_name('slat    ',dim_names)

if (alloc_ps) then
   allocate (ps(dim_sizes(lon_index), dim_sizes(lat_index)))              
   if (slon_index /= 0) &
      allocate( ps_stagr_lon (dim_sizes(slon_index), dim_sizes( lat_index)))
   if (slat_index /= 0) &
      allocate( ps_stagr_lat (dim_sizes( lon_index), dim_sizes(slat_index)))
   alloc_ps = .false.
end if

! assign values to ps grids for use by the rest of the module.
! assuming ps is the first 2D field in state_num_2d
fld_index_2d = find_name('PS      ',state_names_2d)
fld_index    = find_name('PS      ',cflds)
ind = index_from_grid(1,1,1,fld_index) -1
do n=1,s_dim_2d(2,fld_index_2d)
do m=1,s_dim_2d(1,fld_index_2d)
   ind = ind + 1
   ps(m,n) = vec(ind)
end do
end do
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
  integer i,k             ! Longitude, level indices
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


!===============================================================================
   subroutine model_heights(vec, p_surf, lon_index, lat_index, num_levs, &
                            stagr_lon, stagr_lat, model_h, istatus)
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

implicit none

! type(model_type), intent(in) :: Var
real(r8),         intent(in)  :: vec(:)
real(r8),         intent(in)  :: p_surf
integer,          intent(in)  :: lon_index, lat_index, num_levs
logical,          intent(in)  :: stagr_lon, stagr_lat
integer,          intent(out) :: istatus

! OUTPUT: geometrical height at midlayer (m)  hui liu /03/29/2004 added.
real(r8),      intent(out) ::  model_h(num_levs)

! local variables; ps must be dimensioned as an array because dcz2 has it that way
real (r8):: phi(num_levs), tv(num_levs), q(num_levs), t(num_levs)
real (r8):: pmln(num_levs+1), hyba(num_levs+1,2), hybb(num_levs+1,2), pterm(num_levs)
real (r8):: phi_surf, gmh, ht_tmp, rd, rv, rr_factor, local_lat
            
integer :: k, i, vstatus
logical :: stagr_lon_local, stagr_lat_local

istatus = 0
vstatus = 0
! Scratch arrays

rd = 287.05_r8
rv = 461.51_r8
rr_factor = (rv/rd) - 1.0_r8

DO k = 1,num_levs
   model_h(k) = 0.0_r8
   phi(k)     = 0.0_r8
   PTERM(k)   = 0.0_r8
END DO

! copy to temporary arrays

!    All arrays except hyba, hybb are oriented top to bottom
!    modified to be consistent with CAM3.0 where hyai, hyam are top to bottom
!    H Liu, 04/05/2004

! interface=1 extra level at model bottom
k = num_levs +1
hyba(1,1) = hyai%vals(k)
hybb(1,1) = hybi%vals(k)

! mid-points=2
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

stagr_lat_local = stagr_lat
if     (stagr_lat) then
   phi_surf = phis_stagr_lat(lon_index, lat_index)
   local_lat = slat%vals(lat_index)
   ! Don't look for points farther north for interpolation
   if (lat_index == dim_sizes(find_name('slat    ',dim_names))) stagr_lat_local = .false.
elseif (stagr_lon) then
   phi_surf = phis_stagr_lon(lon_index, lat_index)
   local_lat = lat%vals(lat_index)
else
   phi_surf = phis(lon_index, lat_index)
   local_lat = lat%vals(lat_index)
end if

! merge/MPI
call get_interp_prof (q,vec,num_levs, lon_index, lat_index, stagr_lon, stagr_lat_local, &
                      TYPE_Q, vstatus)
call get_interp_prof (t,vec,num_levs, lon_index, lat_index, stagr_lon, stagr_lat_local, &
                      TYPE_T, vstatus)
! if (vstatus == 0) call get_val(q(k), x, lon_index, lat_index, k, KIND_SPECIFIC_HUMIDITY, vstatus)
! if (vstatus == 0) call get_val(t(k), x, lon_index, lat_index, k, KIND_TEMPERATURE      , vstatus)

! Calculate tv for this column, for use by dcz2
if (vstatus == 0) then
   do k = 1, num_levs
      tv(k) = t(k)*(1.0_r8 + rr_factor*q(k))
   end do
elseif (vstatus > 0) then
   istatus = 1
   return
end if

call dcz2(p_surf, phi_surf, tv, P0%vals(1) ,hyba, hybb, num_levs, pmln, pterm, phi)

! used; hybb, hyba, hprb
! calced in dcz2;  pmln, pterm , zslice

do k = 1,num_levs
   ht_tmp = phi(k) * 0.001        ! convert to km for following call only
   model_h(k) = gph2gmh (ht_tmp, local_lat) * 1000.0           ! convert back to m
end do

end subroutine  model_heights


   subroutine get_interp_prof (prof, vec, num_levs, lon_index, lat_index, stagr_lon, stagr_lat, &
                            kind_cam, vstatus)
!=====================================================================

real(r8), intent(out) :: prof(num_levs) 
integer,  intent(out) :: vstatus
! type(model_type), intent(in) :: state
real(r8), intent(in)  :: vec(:)
integer,  intent(in)  :: kind_cam, num_levs, lon_index, lat_index
logical,  intent(in)  :: stagr_lon, stagr_lat

real(r8)  :: var(num_levs,0:1,0:1), weight 
integer   :: k, m, n, fld_index, lon_index_local, vec_ind
character :: cfld

! lon_index and lat_index are the indices of the point at which we need the profiles,
!   not necessarily one of the points where the profiles exist.
! stagr_xx tells where to look for profiles and what interpolating to do.
! staggered longitudes have indices the same as the next eastward  (+) unstaggered longitudes
! Staggered latitudes  have indices the same as the next southward (-) unstaggered lats.
! The indices called for here look wierd, but think of it this way;
!    a point which is staggered from the usual grid has the same indices as the unstaggered point
!    to it's southeast.  The northwest corner is then one higher latitude, and one less (westward)
!    longitude.
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
      implicit none
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

! Surface geoptential
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

! Num of longitude points to compute
      integer imax
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
!       pmln(k+1) is pressure at the midpoint of layer k

          PMLN(1) = DLOG(HPRB*HYBA(KMAX,2)+P_surf *HYBB(KMAX,1))
          PMLN(KMAX+1) = DLOG(HPRB*HYBA(1,2)+P_surf *HYBB(1,1))

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
             Z2(K) = Z2(K) + PTERM(L)
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
      implicit none
      
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

  implicit none
  real (r8) :: xmu, ae, f, w, xm, f2, f4, ge, g, galt, xlat,alt
!
      xmu = 398600.4415_r8       ! km^3/s^2
      ae = 6378.1363_r8          ! km
      f = 1.0_r8/298.2564_r8
      w = 7.292115d-05          ! rad/s
      xm = 0.003468_r8           !
!     f2 = -f + 5.0/2.0*xm - 17.0/14.0*f*xm + 15.0/4.0*xm**2
!     f4 = -f**2/2.0 + 5.0/2.0*f*xm
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
      galt = g - 2.0*ge*alt/ae*(1.0 + f + xm + (-3.0*f + 5.0/2.0*xm)*  &
                             (dsin(xlat))**2) + 3.0*ge*alt**2/ae**2
!
!liu     galt = galt*1.0d5		! convert from km/s2 to cm/s2
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

   subroutine adv_1step(x, Time)
!=======================================================================
! subroutine adv_1step(x, Time)
!
! Does single time-step advance for B-grid model with vector state as
! input and output. This is a modified version of subroutine atmosphere
! in original bgrid_solo_atmosphere driver.

real(r8), intent(inout) :: x(:)

! Time is needed for more general models like this; need to add in to
! low-order models
type(time_type), intent(in) :: Time

! This is a no-op for CAM; only asynch integration
! Can be used to test the assim capabilities with a null advance

end subroutine adv_1step




   subroutine end_model()
!=======================================================================
! subroutine end_model()
!
! At some point, this stub should coordinate with atmosphere_end but
! that requires an instance variable.

! release the local copy of the ensemble means.
deallocate(ens_mean)

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



!#######################################################################
! end of cam model_mod
!#######################################################################

end module model_mod


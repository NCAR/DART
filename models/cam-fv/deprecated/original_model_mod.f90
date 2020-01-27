! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! $Id$

! > > > This version has NOT been updated to describe RMA changes.
!       See RMA-KR for changes to Helen's original RMA version.
!       Some comments in here are meant to connect with comments in 
!       the trunk (non-RMA) version as of 2016-7.  These comments
!       and the sections in the trunk may be helpful in tracing the
!       development of the RMA for FV, and help with the development
!       of the RMA SE version.


!>  This is the interface module between remote memory access capable DART (RMA) 
!>  and the atmospheric components of CESM; CAM, WACCM, CAM-Chem (, ...?).  
!>  It contains the required 16 interface procedures, as specified by DART.  
!>  It also contains several utility routines which help translate between CAM and DART 
!>  formats, and deal with time.
!>  It is used by filter and perfect_model_obs.
!>
!>  This module handles the finite volume dynamical core version of CAM.
!>  A separate model_mod will handle CAM-SE, the spectral element dycore.
!>  CAM-FV uses a logically rectangular grid,
!>  while CAM-SE uses the cubed sphere (non-rectangular) horizontal grid.
!>
!>  There is a perturburbation routine for generating and initial ensemble.
!>  This routine is activated by the filter namelist logical perturb_from_single_instance
!>  and the model_mod namelist variable pert_names.
!>  The module does not provide adv_1step or init_conditions because CAM 
!>  is a separate executable and cannot be called as a subroutine.
!>
!>  This module intercepts the get_close_obs() calls and can alter the distances
!>  for obs near the top of the model to reduce the impact on the state near the
!>  top.
!>
!>  The coordinate orders of fields are preserved from the CAM initial file order.
!>
!>  The RMA model_mod does not refer to TYPE_s, since they were replaced by association
!>  with CAM variables and use of find_name. 
!>  In the future, DART QTYs will be associated with CAM variables by the ${comp}_variables
!>  mechanism as in models/clm.
!>  If a user wants to add new CAM variables to the state vector,
!>  then more QTY_s may be needed in the 'use obs_kind_mod' statement and maybe the obs_kind_mod.
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
!>       get_close section
!>       Utility routines; called by several main subroutines
!>       Stubs not used by cam/model_mod (this is not all of them)
!>
!>  See the subversion code logs for history of this module.
!>  There is an html description of this module in ./model_mod.html.
!>

module model_mod

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! CONTRIBUTORS  (aside from the DART team)

!  Ave Arellano did the first work with CAM-Chem, assimilating MOPPITT CO observations
!  into CAM-Chem using the FV core.  Jerome Barre and Benjamin Gaubert took up the development
!  work from Ave, and prompted several additions to DART, as well as cam/model_mod.

!  Nick Pedatella developed the first vert_coord = 'log_invP' capability
!  to enable assimilation using WACCM and scale height vertical locations.

! NOTES about the module.

!  This module no longer (RMA) keeps a copy of the ensemble mean in module global storage.
!  That was needed for convert_vert to transform the vertical coordinate of something
!  passed from filter into the coordinate used in model_mod.  But now convert_vert is
!  called directly by filter, where the model states and or mean are available, 
!  so ens_mean is not needed.
!  All locations are now converted to a standard coordinate
!  (pressure or log(P0/pressure), aka scale height), instead of always converting the state 
!  vertical location to that of the ob.  The highest_obs_level and ..._height_m variables 
!  are derived from highest_obs_pressure_Pa namelist variable.
!
!  Surface pressure may be needed on the A-grid (thermodynamic variables) and grids staggered
!  relative to the A-grid (winds).   These are retrieved (A-grid) and/or calculated (staggered)
!  as needed from filter, rather than being stored globally in this module.

!  The coordinates of CAM (lats, lons, etc.) and their dimensions  and attributes are
!  read into globally accessible data structures (see grid_1d_type).
!
!     MODULE ORGANIZATION (search for the following strings to find the corresponding section)
!
!          USE statements
!          Global storage for describing cam model class
!          Namelist variables with default values
!          Derived parameters
!          static_init_model section
!          Module I/O to/from DART and files
!          model_interpolate section
!          Vector-field translations
!          get_close section
!          Utility routines; called by several main subroutines
!          Stubs not used by cam/model_mod (this is not all of them)

!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

!  USE statements

use netcdf
use typeSizes

use types_mod,         only : r8, MISSING_I, MISSING_R8, gravity_const => gravity, &
                              PI, DEG2RAD, RAD2DEG, obstypelength, earth_radius, i8
! FIXME; these constants should be consistent with CESM, not necessarily with DART.
!          add after verification against Hui's tests;  gas_constant_v,gas_constant,ps0,PI,DEG2RAD

use time_manager_mod,  only : time_type, set_time, set_date, print_time, print_date,    &
                              set_calendar_type, get_calendar_type, get_time, get_date, &
                              operator(-), operator(==)

use utilities_mod,     only : open_file, close_file, find_namelist_in_file, check_namelist_read, &
                              register_module, error_handler, file_exist, E_ERR, E_WARN, E_MSG,  &
                              logfileunit, nmlfileunit, do_output, get_unit, do_nml_file, &
                              do_nml_term

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_sync, nc_check, &
                                 nc_add_global_creation_time, nc_redef, nc_enddef

use mpi_utilities_mod, only : my_task_id, task_count

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
use location_mod,      only : location_type, get_location, set_location, query_location,         &
                              is_vertical,                                                       &
                              VERTISUNDEF, VERTISSURFACE, VERTISLEVEL,                           &
                              VERTISPRESSURE, VERTISHEIGHT, VERTISSCALEHEIGHT, write_location,   &
                              get_close_type, get_dist, loc_get_close => get_close

! READ THIS SYNTAX as:
!   There's a subroutine in location_mod named 'get_close'.
!   If I want to use that one in this module then refer to it as 'loc_get_close'.
!   If I call 'get_close', then I'll get the one in this module,
!   which does some stuff I need, AND ALSO CALLS 'loc_get_close'

! FIXME
! I've put a copy of solve_quadratic in this model_mod.
! Eventually it should go into a utilities module.
! use utilities_YYY, only : solve_quadratic

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
use     obs_kind_mod, only : QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, QTY_PRESSURE,       &
                             QTY_SURFACE_PRESSURE, QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY,   &
                             QTY_CLOUD_LIQUID_WATER, QTY_CLOUD_ICE, QTY_CLOUD_FRACTION,      &
                             QTY_GRAV_WAVE_DRAG_EFFIC, QTY_GRAV_WAVE_STRESS_FRACTION,        &
                             QTY_SURFACE_ELEVATION,                                          &
                             QTY_CO, QTY_CO2, QTY_NO, QTY_NO2, QTY_CH4, QTY_NH3, QTY_O3,     &
                             QTY_AOD, QTY_CO01, QTY_CO02, QTY_CO03,                          &
                             QTY_SFCO, QTY_SFCO01, QTY_SFCO02, QTY_SFCO03,                   &
                             QTY_CB1, QTY_CB2, QTY_OC1, QTY_OC2,                             &
                             QTY_SFCB1, QTY_SFCB2, QTY_SFOC1, QTY_SFOC2,                     &
                             QTY_CB102, QTY_CB202, QTY_OC102, QTY_OC202,                     &
                             QTY_SFCB102, QTY_SFCB202, QTY_SFOC102, QTY_SFOC202,             &
                             get_index_for_quantity, get_name_for_quantity, get_quantity_for_type_of_obs


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
! and chemical species from WACCM and CAM-Chem.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

use random_seq_mod,        only : random_seq_type, init_random_seq, random_gaussian

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state, get_state_array

use state_structure_mod,   only : add_domain, get_model_variable_indices, get_dim_name, &
                                  get_num_dims, get_domain_size, get_dart_vector_index, &
                                  get_index_start, get_index_end

use default_model_mod,    only : adv_1step, init_time, init_conditions, nc_write_model_vars

! end of use statements
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

! CAM global/module declarations

implicit none
private

! The first block are the 16 required interfaces.  The following block
! are additional useful interfaces that utility programs can call.
public ::                                                      &
   static_init_model, get_model_size,                          &
   shortest_time_between_assimilations,                        &
   pert_model_copies, get_state_meta_data, model_interpolate,  &
   nc_write_model_atts, nc_write_model_vars,                   &
   init_conditions, init_time, adv_1step, end_model,           &
   get_close_obs, get_close_state,                             &
   convert_vertical_obs, convert_vertical_state,               &
   query_vert_localization_coord, read_model_time, write_model_time

public ::                                                   &
   model_type, prog_var_to_vector, vector_to_prog_var,      &
   read_cam_init,                                           &
   init_model_instance, end_model_instance, write_cam_init, &
   write_cam_times

!-----------------------------------------------------------------------
! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
!-----------------------------------------------------------------------

integer :: component_id ! for add_domain.

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Global storage for describing cam model class
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

!>@todo FIXME: this should be an i8 to handle really large
!> state vectors, but that change ripples through several layers
!> of code.  nevertheless it should be done.
integer :: model_size

! This list of dimensions used to define fields will be ordered as they are on the caminput.nc file.
integer                                   :: num_dims
integer,                      allocatable :: dim_sizes(:)
character(len=NF90_MAX_NAME), allocatable :: dim_names(:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Grid fields
! These structures are used by nc_write_model_atts.
! They are dimensioned in create_grid_1d_instance and filled in read_cam_coord.

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

integer, parameter :: no_lev = MISSING_I  ! constant to tell get_val_level there are no levels.
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

! grid_2d_type ?
! integer :: grid_num_2d = 0              ! # of 2d grid fields to read from file
! ? should phis be in grid_names_2d?
! character (len=8),dimension(100) :: grid_names_2d = (/(' ',iii=1,100)/)

! CAM-chem 
! These lists were provided by Jerome Barre' and/or Avelino Arellano.
! They implemented the unit conversion in subroutine read_cam_init in Lanai (and/or earlier DARTs).
! The Manhattan implementation at the end of model_interpolate was by Kevin Raeder.
! FIXME It would be better if the following 2 vectors were read from an external file....
! If meteorological variables (including PRESSURE), or SURFACE_ELEVATION need to have 
! their units converted, their names and conversion factors could be entered in these lists.

integer, parameter :: chemical_list=128
! Names of chemical species. 
character(len=16)  :: solsym(chemical_list) = &
(/'O3              ','O               ','O1D             ','N2O             ','NO              ', &
  'NO2             ','NO3             ','HNO3            ','HO2NO2          ','N2O5            ', &
  'H2              ','OH              ','HO2             ','H2O2            ','CH4             ', &
  'CO              ','CH3O2           ','CH3OOH          ','CH2O            ','CH3OH           ', &
  'C2H5OH          ','C2H4            ','EO              ','EO2             ','CH3COOH         ', &
  'GLYALD          ','C2H6            ','C2H5O2          ','C2H5OOH         ','CH3CHO          ', &
  'CH3CO3          ','CH3COOOH        ','C3H6            ','C3H8            ','C3H7O2          ', &
  'C3H7OOH         ','PO2             ','POOH            ','CH3COCH3        ','RO2             ', &
  'ROOH            ','BIGENE          ','ENEO2           ','MEK             ','MEKO2           ', &
  'MEKOOH          ','BIGALK          ','ALKO2           ','ALKOOH          ','ISOP            ', &
  'ISOPO2          ','ISOPOOH         ','MVK             ','MACR            ','MACRO2          ', &
  'MACROOH         ','MCO3            ','HYDRALD         ','HYAC            ','CH3COCHO        ', &
  'XO2             ','XOOH            ','C10H16          ','TERPO2          ','TERPOOH         ', &
  'TOLUENE         ','CRESOL          ','TOLO2           ','TOLOOH          ','XOH             ', &
  'BIGALD          ','GLYOXAL         ','PAN             ','ONIT            ','MPAN            ', &
  'ISOPNO3         ','ONITR           ','SOA             ','SO2             ','DMS             ', &
  'NH3             ','NH4             ','NH4NO3          ','Rn              ','Pb              ', &
  'HCN             ','CH3CN           ','C2H2            ','HCOOH           ','HOCH2OO         ', &
  'H2SO4           ','SOAG            ','so4_a1          ','pom_a1          ','soa_a1          ', &
  'bc_a1           ','dst_a1          ','ncl_a1          ','num_a1          ','so4_a2          ', &
  'soa_a2          ','ncl_a2          ','num_a2          ','dst_a3          ','ncl_a3          ', &
  'so4_a3          ','num_a3          ','CO01            ','CO02            ','CO03            ', &
  'CO04            ','CO05            ','CO06            ','CO07            ','CO08            ', &
  'CO09            ','CB1             ','CB2             ','OC1             ','OC2             ', &
  'CB101           ','CB201           ','OC101           ','OC201           ', &
  'CB102           ','CB202           ','OC102           ','OC202           '  &
  /)

! The molar mass of each chemical species
real(r8) :: adv_mass(chemical_list) =  &
(/47.9982_r8,     15.9994_r8,     15.9994_r8,     44.01288_r8,  30.00614_r8,    &
  46.00554_r8,    62.00494_r8,    63.01234_r8,    79.01174_r8,  108.01048_r8,   &
  2.0148_r8,      17.0068_r8,     33.0062_r8,     34.0136_r8,   16.0406_r8,     &
  28.0104_r8,     47.032_r8,      48.0394_r8,     30.0252_r8,   32.04_r8,       &
  46.0658_r8,     28.0516_r8,     61.0578_r8,     77.0572_r8,   60.0504_r8,     &
  60.0504_r8,     30.0664_r8,     61.0578_r8,     62.0652_r8,   44.051_r8,      &
  75.0424_r8,     76.0498_r8,     42.0774_r8,     44.0922_r8,   75.0836_r8,     &
  76.091_r8,      91.083_r8,      92.0904_r8,     58.0768_r8,   89.0682_r8,     &
  90.0756_r8,     56.1032_r8,     105.1088_r8,    72.1026_r8,   103.094_r8,     &
  104.1014_r8,    72.1438_r8,     103.1352_r8,    104.1426_r8,  68.1142_r8,     &
  117.1198_r8,    118.1272_r8,    70.0878_r8,     70.0878_r8,   119.0934_r8,    &
  120.1008_r8,    101.0792_r8,    100.113_r8,     74.0762_r8,   72.0614_r8,     &
  149.1186_r8,    150.126_r8,     136.2284_r8,    185.234_r8,   186.2414_r8,    &
  92.1362_r8,     108.1356_r8,    173.1406_r8,    174.148_r8,   190.1474_r8,    &
  98.0982_r8,     58.0356_r8,     121.04794_r8,   119.07434_r8, 147.08474_r8,   &
  162.11794_r8,   147.12594_r8,   144.132_r8,     64.0648_r8,   62.1324_r8,     &
  17.02894_r8,    18.03634_r8,    80.04128_r8,    222.0_r8,     207.2_r8,       &
  27.02514_r8,    41.05094_r8,    26.0368_r8,     46.0246_r8,   63.0314_r8,     &
  98.0784_r8,     12.011_r8,      115.10734_r8,   12.011_r8,    12.011_r8,      &
  12.011_r8,      135.064039_r8,  58.442468_r8,   1.0074_r8,    115.10734_r8,   &
  12.011_r8,      58.442468_r8,   1.0074_r8,      135.064039_r8,58.442468_r8,   &
  115.10734_r8,   1.0074_r8,      28.0104_r8,     28.0104_r8,   28.0104_r8,     &
  28.0104_r8,     28.0104_r8,     28.0104_r8,     28.0104_r8,   28.0104_r8,     &
  28.0104_r8,     12.011_r8,      12.011_r8,      12.011_r8,    12.011_r8,      &
  12.011_r8,      12.011_r8,      12.011_r8,    12.011_r8,      &
  12.011_r8,      12.011_r8,      12.011_r8,    12.011_r8       &
/)

! 2 unit conversion arrays derived from adv_mass will be filled in map_qtys.
real(r8), parameter :: molar_mass_dry_air = 28.9644_r8

! CAM-chem end

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Namelist variables with default values follow

! Files where basic info about model configuration can be found
character(len=128) :: &
   model_config_file = 'caminput.nc',             & ! An example cam initial file.
   cam_phis          = 'cam_phis.nc',             & ! Separate source of PHIS/topography.
   model_version     = '6.0'


! Define location restrictions on which observations are assimilated
! (values are calculated anyway, but istatus is set to 2)
! RMA-KR; This would be a good time to change 'log_invP' to 'scale_ht' or 'scaled_h?
character(len=8) :: vert_coord = 'pressure'            ! or 'log_invP'
real(r8) :: max_obs_lat_degree        = 90.0_r8
real(r8) :: highest_obs_pressure_Pa   = 1000.0_r8
real(r8) :: highest_state_pressure_Pa = 9400.0_r8

! Namelist variables and default values for defining state vector.

integer :: state_num_0d = 0              ! # of scalars fields to read from file
integer :: state_num_1d = 0              ! # of 1d fields to read from file
integer :: state_num_2d = 0              ! # of 2d fields to read from file
integer :: state_num_3d = 0              ! # of 3d fields to read from file

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

! Special for an experiment.  Specify one string kind e.g QTY_CLOUD_LIQUID and
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


namelist /model_nml/ vert_coord, model_version, cam_phis,        &
                       state_num_0d,   state_num_1d,   state_num_2d,   state_num_3d,  &
                     state_names_0d, state_names_1d, state_names_2d, state_names_3d,  &
                                      which_vert_1d,  which_vert_2d,  which_vert_3d,  &
                     pert_names, pert_sd, pert_base_vals,                             &
                     highest_obs_pressure_Pa, highest_state_pressure_Pa,              &
                     max_obs_lat_degree, Time_step_seconds, Time_step_days,           &
                     impact_only_same_kind, print_details,                            &
                     model_config_file

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
! Random sequence and init for pert_model_copies
type(random_seq_type)   :: random_seq
integer                 :: ens_member = 0
logical                 :: output_task0

! common message string used by many subroutines
character(len=512) :: string1, string2, string3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
integer :: nflds         ! # fields to read

! f_dim_#d are the sizes of the coordinates of each variable as found on caminput file.
! RMA-KR
! s_dim_#d and s_dimid_#d are no longer needed, because this model mod is specialized 
!    for a single CAM; dynamical core, coordinate orders in the initial file, etc.
integer, allocatable :: f_dim_3d(:,:), f_dim_2d(:,:), f_dim_1d(:,:),  &
                        f_dimid_3d(:,:), f_dimid_2d(:,:), f_dimid_1d(:,:)

integer :: f_dim_max(4,3)


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Surface pressures, used by vertical interpolation routines.
!
! I assume that staggered grids (US and VS) are staggered only 1 direction (each),
! so that surface pressure interpolations to get staggered ps use only 2 A-grid ps values.
! The interpolations for columns of heights are more general, but will do a 2 point interp
!     if the staggering is only in one direction.

! height
! Surface potential; used for calculation of geometric heights.
logical               :: alloc_phis=.true.    ! Flag whether to allocate space for phis
real(r8), allocatable :: phis(:, :)           ! surface geopotential

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! RMA-KR; cubed sphere (CAM-SE) section removed from here.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Array 'cflds' is filled with simple loops over state_names_xxx.
! I could replace that with code which orders namelist input field names
! into cflds, regardless of their original order, and tallies how many of each.
! Is there a way to exclude state_nums from namelist and have those filled in
! the same subroutine?
! RMA-KR; this may/will be replaced by the ${comp}_variables mechanism.

character(len=8), allocatable :: cflds(:)

! Attribute values for the fields which comprise the state vector.
! These are filled by nc_read_model_atts.
character(len=nf90_max_name), allocatable :: state_long_names(:)
character(len=nf90_max_name), allocatable :: state_units(:)

! Arrays for linking obs_qtys(QTY_) and model variable names are filled in map_qtys.
! The max size of QTY_ should come from obs_kind_mod
! These should be dimensioned the same size as the total of state_names_Nd.
character(len=8) :: dart_to_cam_types(300) = ''
integer          :: cam_to_dart_qtys(300) = MISSING_I
! Strategy; array elements are only changed for conversion factors that are != 1.0.
!           Then convert_mmr2vmr = MISSING_R8 triggers a convert_units of 1.0 in interp_lonlat.
! So far, the conversion from obs units back to state units is no needed.
! If it becomes needed: 
! 1) define array convert_vmr2mmr(MAX_STATE_NAMES) = MISSING_R8 
! 2) Add lines to function map_qtys similar to the convert_mmr2vmr lines:
!       convert_vmr2mmr(i) = 1.0_r8/convert_mmr2vmr(i)
real(r8)         :: convert_mmr2vmr(MAX_STATE_NAMES) = MISSING_R8

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
! No TYPE_s and SE code

subroutine static_init_model()

! Initializes class data for CAM model (all the stuff that needs to be done once).
! For now, does this by reading info from a fixed name netcdf file.

integer  :: iunit, io, i, nc_file_ID
integer  :: max_levs
real(r8), allocatable :: clampfield(:,:)
! RMA-KR; clampfield added to assist restricting the range of some variable values.

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

! Set the printed output logical variable to reduce printed output;
output_task0 = do_output()

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
call read_cam_init_size(nc_file_ID)

! RMA-KR; model size is now calculated in state_structure_mod/get_domain_size

! Allocate space for global coordinate arrays and read them in.
! There's a query of caminput.nc within read_cam_coord for the existence of the field.
! The second argument is a grid_1d_type structure
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
call read_cam_coord(nc_file_ID, 'P0', P0)    ! thats a p-zero

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
   endif
else if (vert_coord == 'log_invP') then
   highest_state_scale_h = scale_height(p_surface=P0%vals(1), p_above=highest_state_pressure_Pa)
   model_top             = scale_height(p_surface=P0%vals(1), p_above=(hyai%vals(1)*P0%vals(1)) )
   highest_state_loc = set_location(1.0_r8,1.0_r8,highest_state_scale_h,VERTISSCALEHEIGHT)
   model_top_loc     = set_location(1.0_r8,1.0_r8,model_top,            VERTISSCALEHEIGHT)
   if (highest_state_scale_h /= model_top) then
      damp_wght = 1.0_r8/get_dist(highest_state_loc,model_top_loc,no_vert=.false.)
   endif
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
! Construct array of variables to be clamped - used in filter netcdf write not in write_cam_init
allocate(clampfield(nflds, 2))
call set_clamp_fields(clampfield)

! Add a component to the state vector
component_id = add_domain('caminput.nc', nflds, cflds, clamp_vals = clampfield)
deallocate(clampfield)

! Compute overall model size and put in global storage
model_size = get_domain_size(component_id)
if (output_task0) then
   write(string1, '(A,I9)') 'CAM state vector size: ', model_size
   call error_handler(E_MSG, 'static_init_model', string1)
endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Get field attributes needed by nc_write_model_atts from caminput.nc.
allocate(state_long_names(nflds), state_units(nflds))    
call nc_read_model_atts(nc_file_ID, 'long_name', state_long_names)
call nc_read_model_atts(nc_file_ID, 'units', state_units)

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

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! RMA-KR; 
!    p_col is now a local variable, allocated when/where it's needed.
! Fills arrays for the linking of obs_qtys (QTY_) to model field names
call map_qtys()

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! If restricting impact of a particular kind to only obs and state vars
! of the same kind, look up and set the kind index.
! RMA-KR This will/may be replaced by Nancy's more general code for restricting
!        the influence of obs on listed variables.
if (len_trim(impact_only_same_kind) > 0) then
   impact_kind_index = get_index_for_quantity(impact_only_same_kind)
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

enddo

! Find and store shapes of all the state vector fields.  Grouped by rank of fields into
! separate f_dim_RANKd arrays.
call read_coord(nc_file_ID)


! The arrays into which CAM fields are put are dimensioned by the largest values of
! the sizes of the dimensions listed in Y_dim_RANKd, Y=[sf], RANK=[1-3] .
! The second dimension denotes the rank of the array for which the first dim
! gives the max size(s).
if (state_num_1d > 0) then
   f_dim_max(1:2, 1) = maxval(f_dim_1d, dim=2)   ! gets the max value of f_dim_1d (1:2, :)
else
   f_dim_max(1:2, 1) = 0
endif

if (state_num_2d > 0) then
   f_dim_max(1:3, 2) = maxval(f_dim_2d, dim=2)   ! gets the max values of f_dim_2d (1:3, :)
else
   f_dim_max(1:3, 2) = 0
endif

if (state_num_3d > 0) then
   f_dim_max(1:4, 3) = maxval(f_dim_3d, dim=2)   ! gets the max values of f_dim_3d (1:4, :)
else
   f_dim_max(1:4, 3) = 0
endif

if (print_details .and. output_task0 ) then
   if (state_num_1d > 0) then
      write(string1,*) 'f_dim_1d = ',f_dim_1d
      write(string2,*) (f_dim_max(i,1),i=1,3)
      call error_handler(E_MSG, 'read_cam_init_size', string1,source,revision,revdate, text2=string2)
   endif

   do i=1,2
      write(string1,*) 'f_dim_2d = ',(f_dim_2d(i,j),j=1,state_num_2d),'f_dim_max = ',f_dim_max(i,2)
      call error_handler(E_MSG, 'read_cam_init_size', string1,source,revision,revdate)
   enddo

   do i=1,3
      write(string1,'(A,(10I4))') 'f_dim_3d = ',(f_dim_3d(i,j),j=1,state_num_3d)
      write(string2,'(A,(10I4))') 'f_dim_max = ',f_dim_max(i,3)
      call error_handler(E_MSG, 'read_cam_init_size', string1,source,revision,revdate, text2=string2)
   enddo
endif

end subroutine read_cam_init_size

!-----------------------------------------------------------------------

subroutine read_coord(nc_file_ID)

! Figure out which coordinates are lon, lat, lev, based on CAM version
! from the namelist, which has form #.#[.#[.#]].

integer,            intent(in) :: nc_file_ID

! local workspace
character(len=4) :: form_version = '(I0)'
character(len=4) :: char_version
integer          :: part, nchars, tot_chars, i, j, k, varid, next
integer          :: int_version(4)

int_version = (/ (0,i=1,4) /)

! Choose order . . .  no longer needed because this model_mod is specialized to 
! CAM-FV in CESM1.x and later.

! Cycle through each field's dimension IDs.
! Pick the dimensions needed out of dim_sizes, using the dimension names in dim_names.
! Fill the state dimids according to the order model_mod wants to see.  (lev, lon, lat).

! 3D
if (state_num_3d > 0) then
   allocate(f_dim_3d(4,state_num_3d), f_dimid_3d(4,state_num_3d))
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
   enddo Alldim3

   if (   f_dim_3d(1,i) == 0 .or.  f_dim_3d(2,i) == 0 .or.  f_dim_3d(3,i) == 0 ) then
      call error_handler(E_ERR, 'trans_coord', &
          'num_[lons,lats,levs] was not assigned and = 0' , source, revision, revdate)
   endif
enddo

! 2D
if (state_num_2d > 0) then
   allocate(f_dim_2d(3,state_num_2d), f_dimid_2d(3,state_num_2d))
   f_dim_2d = 0;  f_dimid_2d  = 0;
endif

do i = 1,state_num_2d
   call nc_check(nf90_inq_varid(nc_file_ID, state_names_2d(i), varid), &
              'trans_coord', 'inq_varid '//trim(state_names_2d(i)))
   call nc_check(nf90_inquire_variable(nc_file_ID, varid, dimids=f_dimid_2d(1:3,i)), &
              'trans_coord', 'inquire_variable '//trim(state_names_2d(i)))

   ! Extract spatial dimids from the fields dimids
   Alldim2: do j = 1,3      ! time and 2 space
      k = f_dimid_2d(j,i)
      f_dim_2d(j,i) = dim_sizes(k)
   enddo Alldim2
   if (   f_dim_2d(1,i) == 0 .or.  f_dim_2d(2,i) == 0 ) then
      call error_handler(E_ERR, 'trans_coord', &
          'num_[lons,lats,levs] was not assigned and = 0' , source, revision, revdate)
   endif
enddo

! 1D
if (state_num_1d > 0) then
   allocate(f_dim_1d(2,state_num_1d), f_dimid_1d(2,state_num_1d))
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
   enddo Alldim1

   if ( f_dim_1d(1, i) == 0 ) then
      write(string1, '(A,I3,A)') ' state 1d dimension(',i,') was not assigned and = 0'
      call error_handler(E_ERR, 'trans_coord',trim(string1), source, revision, revdate)
   endif
enddo

end subroutine read_coord

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
      name_dim2 = 'no2ndDim'
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
         if (num_dim2 /= i_dim2) then
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
      name_dim2 = 'no2ndDim'
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
                 count=(/ num_dim1, num_dim2 /)), 'read_cam_2Dint', trim(cfield))
else
   call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, field),  &
                  'read_cam_2Dint', trim(cfield))
endif

call nc_check(nf90_close(nc_file_ID), 'read_cam_2Dint', 'closing '//trim(file_name))

end subroutine read_cam_2Dint

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

integer,            intent(in)  :: nc_file_ID
character(len=*),   intent(in)  :: cfield
type(grid_1d_type), intent(out) :: var

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

! Moving this from the specification statement to here caused it to
! be initialized every time read_cam_coord is called, and 'broke' it.  
! Previously, P0 may have ended up using the value left over from the last call,
! which was for one of the initial file dimension variables, which was wrong,
! but seems to have worked.  
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

if (print_details .and. output_task0) then
   write(string1,*) 'After inquire_variable for ',cfield,' coord_dimid = ',coord_dimid(1)
   call error_handler(E_MSG, 'read_cam_coord', string1,source,revision,revdate)
endif

if (coord_dimid(1) == MISSING_I) then
   ! to handle P0
   coord_size = 1                 
   coord_dimid(1) = 0      ! This is the dimid for time, which has length 1,
                           ! But time is the record/unlimited dimension, so this may not work.
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
   call error_handler(E_MSG, 'read_cam_coord', string1,source,revision,revdate, text2=string2)
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

integer :: i, i1, nfld

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
enddo

! 2D fields
do i=1,state_num_2d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_2d(i)
enddo

! 3D fields (including q)
do i=1,state_num_3d
   nfld = nfld + 1
   cflds(nfld)(:) = state_names_3d(i)
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
         write(string1,'(A,I4)') cflds(i)
         call error_handler(E_MSG, 'order_state_fields', string1, source, revision, revdate)
      end do
      i1 = state_num_0d
      do i=1,state_num_1d
         write(string1,'(A,I4)') cflds(i1+i)
         call error_handler(E_MSG, 'order_state_fields', string1, source, revision, revdate)
      end do
      i1 = i1 + state_num_1d
      do i=1,state_num_2d
         write(string1,'(A,I4)') cflds(i1+i)
         call error_handler(E_MSG, 'order_state_fields', string1, source, revision, revdate)
      end do
      i1 = i1 + state_num_2d
      do i=1,state_num_3d
         write(string1,'(A,I4)') cflds(i1+i)
         call error_handler(E_MSG, 'order_state_fields', string1, source, revision, revdate)
      end do
   else
      call error_handler(E_MSG, 'order_state_fields', 'State vector is composed of these fields: ')
      do i = 1,nflds
         call error_handler(E_MSG, 'order_state_fields', trim(cflds(i)))
      enddo
   endif
endif

end subroutine order_state_fields

!-----------------------------------------------------------------------

subroutine map_qtys()

! ? Should this be a function instead; removes need to dimension obs_loc_in arbitrarily
!   and wastefully.  But then it's called millions of times, instead of accessing a small
!   array that's defined once.

! Makes an array of 'locations within the state vector' of the obs kinds
! that come from obs_kind_mod, which we anticipate CAM's model_mod will need.
! The obs kind that's needed will be the index into this array,
! the corresponding value will be the name of that field.
! This name will be used with find_name.
! This subroutine will be called from static_init_model, so it will not have to be
! recomputed for every ob.
! Also maps the model variable names onto the DART QTY_s by the same mechanism.

! other QTY_ possibilities are listed after the 'use obs_kind_mod' statement

integer :: i

! Physically 2D fields

i = find_name('PS',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SURFACE_PRESSURE) = 'PS'
   cam_to_dart_qtys(i) = QTY_SURFACE_PRESSURE
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('AEROD_v',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_AOD) = 'AEROD_v'
   cam_to_dart_qtys(i) = QTY_AOD
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFCO',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFCO) = 'SFCO'
   cam_to_dart_qtys(i) = QTY_SFCO
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFCO01',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFCO01) = 'SFCO01'
   cam_to_dart_qtys(i) = QTY_SFCO01
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFCO02',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFCO02) = 'SFCO02'
   cam_to_dart_qtys(i) = QTY_SFCO02
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFCO03',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFCO03) = 'SFCO03'
   cam_to_dart_qtys(i) = QTY_SFCO03
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFOC1',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFOC1) = 'SFOC1'
   cam_to_dart_qtys(i) = QTY_SFOC1
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFOC2',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFOC2) = 'SFOC2'
   cam_to_dart_qtys(i) = QTY_SFOC2
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFCB1',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFCB1) = 'SFCB1'
   cam_to_dart_qtys(i) = QTY_SFCB1
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFCB2',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFCB2) = 'SFCB2'
   cam_to_dart_qtys(i) = QTY_SFCB2
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFOC102',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFOC102) = 'SFOC102'
   cam_to_dart_qtys(i) = QTY_SFOC102
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFOC202',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFOC202) = 'SFOC202'
   cam_to_dart_qtys(i) = QTY_SFOC202
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFCB102',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFCB102) = 'SFCB102'
   cam_to_dart_qtys(i) = QTY_SFCB102
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('SFCB202',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SFCB202) = 'SFCB202'
   cam_to_dart_qtys(i) = QTY_SFCB202
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('EFGWORO',cflds)
if (i/= MISSING_I) then
   dart_to_cam_types(    QTY_GRAV_WAVE_DRAG_EFFIC) = 'EFGWORO'
   cam_to_dart_qtys(i) = QTY_GRAV_WAVE_DRAG_EFFIC
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('FRACLDV',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_GRAV_WAVE_STRESS_FRACTION) = 'FRACLDV'
   cam_to_dart_qtys(i) = QTY_GRAV_WAVE_STRESS_FRACTION
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

! dart_to_cam_types(QTY_SURFACE_TEMPERATURE  ?  ) = TYPE_TS
! dart_to_cam_types(QTY_SEA_SURFACE_TEMPERATURE  ?  ) = TYPE_TSOCN
!    convert_mmr2vmr(i) = mmr2vmr(i)

! Physically 3D fields
i = find_name('T',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_TEMPERATURE)        = 'T'
   cam_to_dart_qtys(i)      = QTY_TEMPERATURE
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('US',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_U_WIND_COMPONENT)   = 'US'
   cam_to_dart_qtys(i)      = QTY_U_WIND_COMPONENT
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('VS',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_V_WIND_COMPONENT)   = 'VS'
   cam_to_dart_qtys(i)      = QTY_V_WIND_COMPONENT
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('Q',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_SPECIFIC_HUMIDITY)  = 'Q'
   cam_to_dart_qtys(i)      = QTY_SPECIFIC_HUMIDITY
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('CLDLIQ',cflds)
if (i /= MISSING_I)  then
   dart_to_cam_types(    QTY_CLOUD_LIQUID_WATER) = 'CLDLIQ'
  cam_to_dart_qtys(i) = QTY_CLOUD_LIQUID_WATER
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('CLDICE',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_CLOUD_ICE)          = 'CLDICE'
   cam_to_dart_qtys(i) = QTY_CLOUD_ICE
   convert_mmr2vmr(i) = mmr2vmr(i)
endif
!    dart_to_cam_types(QTY_CLOUD_WATER  ?  ) = 'LCWAT'
! cam_to_dart_qtys(i) = QTY_CLOUD_WATER  ?
!    convert_mmr2vmr(i) = mmr2vmr(i)

i = find_name('CO',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_CO)  = 'CO'
   cam_to_dart_qtys(i) = QTY_CO
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('CO01',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_CO01)  = 'CO01'
   cam_to_dart_qtys(i) = QTY_CO01
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('CO02',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_CO02)  = 'CO02'
   cam_to_dart_qtys(i) = QTY_CO02
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('CO03',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_CO03)  = 'CO03'
   cam_to_dart_qtys(i) = QTY_CO03
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('OC1',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_OC1)  = 'OC1'
   cam_to_dart_qtys(i) = QTY_OC1
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('OC2',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_OC2)  = 'OC2'
   cam_to_dart_qtys(i) = QTY_OC2
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('CB1',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_CB1)  = 'CB1'
   cam_to_dart_qtys(i) = QTY_CB1
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('CB2',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_CB2)  = 'CB2'
   cam_to_dart_qtys(i) = QTY_CB2
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('OC102',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_OC102)  = 'OC102'
   cam_to_dart_qtys(i) = QTY_OC102
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('OC202',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_OC202)  = 'OC202'
   cam_to_dart_qtys(i) = QTY_OC202
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('CB102',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_CB102)  = 'CB102'
   cam_to_dart_qtys(i) = QTY_CB102
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('CB202',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_CB202)  = 'CB202'
   cam_to_dart_qtys(i) = QTY_CB202
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('CO2',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_CO2) = 'CO2'
   cam_to_dart_qtys(i) = QTY_CO2
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('NO',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_NO)  = 'NO'
   cam_to_dart_qtys(i) = QTY_NO
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('NO2',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_NO2) = 'NO2'
   cam_to_dart_qtys(i) = QTY_NO2
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('CH4',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_CH4) = 'CH4'
   cam_to_dart_qtys(i) = QTY_CH4
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

i = find_name('NH3',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_NH3) = 'NH3'
   cam_to_dart_qtys(i) = QTY_NH3
   convert_mmr2vmr(i) = mmr2vmr(i)
endif

! i = find_name('O',cflds)
! if (i /= MISSING_I) then
!    dart_to_cam_types(   QTY_O)   = 'O'
!    cam_to_dart_qtys(i) = QTY_O
! endif

i = find_name('O3',cflds)
if (i /= MISSING_I) then
   dart_to_cam_types(    QTY_O3)  = 'O3'
   cam_to_dart_qtys(i) = QTY_O3
   convert_mmr2vmr(i) = mmr2vmr(i)
endif


if (print_details .and. output_task0) then
   write(string1,*) 'OBS_QTY   FIELD_TYPE'
   call error_handler(E_MSG, 'map_qtys', string1,source,revision,revdate)
   do i=1,300
      if (dart_to_cam_types(i) /= '') then
         write(string1,'(I8,A)') i, dart_to_cam_types(i)
         call error_handler(E_MSG, 'map_qtys', string1,source,revision,revdate)
      end if
   end do
end if

end subroutine map_qtys

!-----------------------------------------------------------------------
! CAM-chem 3))
! Function to calculate the unit conversion factors, which make 
! estimated obs have units consistent with actual obs in model_interpolate.

function mmr2vmr(var_index)

integer, intent(in) :: var_index

real(r8) :: mmr2vmr
integer  :: chem_index

mmr2vmr = 1.0_r8
do chem_index=1,chemical_list
   if ( cflds(var_index) .eq. solsym(chem_index) ) then
      mmr2vmr = molar_mass_dry_air/adv_mass(chem_index)
      write(string1,'(2A,I4)') 'State field(= chemical name), mmr2vmr = ', &
             solsym(chem_index), mmr2vmr
      call error_handler(E_MSG, 'mmr2vmr', string1,source,revision,revdate)
      exit
   endif
enddo

end function mmr2vmr

!-----------------------------------------------------------------------

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
   call nc_check(nf90_get_var(nc_file_ID, nc_var_ID, var%vars_1d(1:f_dim_1d(1, i), i), &
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

   var%vars_2d(1:f_dim_2d(1,i),1:f_dim_2d(2,i),i) = &
       temp_2d(1:f_dim_2d(1,i),1:f_dim_2d(2,i))

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

   var%vars_3d(1:f_dim_3d(1,i), 1:f_dim_3d(2,i), 1:f_dim_3d(3,i),i) = &
       temp_3d(1:f_dim_3d(1,i), 1:f_dim_3d(2,i), 1:f_dim_3d(3,i))

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
integer               :: nc_file_ID, nc_var_ID
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
   call nc_check(nf90_put_var(nc_file_ID, nc_var_ID, var%vars_1d(1:f_dim_1d(1, i),i),   &
                             start=(/ 1, timeindex /), count = (/ f_dim_1d(1, i), 1 /)), &
                 'write_cam_init', 'put_var '//trim(cflds(ifld)))
enddo

do i = 1, state_num_2d
   ! special code:  check and error out if the PS field has gone negative
   if (state_names_2d(i) == 'PS') then
      if (minval(var%vars_2d(:,:,i)) < 0.0_r8) then
         write(string1, *)'PS has negative values; should not happen'
         call error_handler(E_ERR, 'write_cam_init', string1, source, revision, revdate)
      endif
   endif

   ! 2d fields ; tricky because coordinates may have been rearranged.

       temp_2d(1:f_dim_2d(1, i),1:f_dim_2d(2,i)) = &
   var%vars_2d(1:f_dim_2d(1, i),1:f_dim_2d(2,i), i )

   ifld = ifld + 1
   call nc_check(nf90_inq_varid(nc_file_ID, trim(cflds(ifld)), nc_var_ID),           &
                 'write_cam_init','inq_varid '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(nc_file_ID, nc_var_ID,    temp_2d(1:f_dim_2d(1, i),1:f_dim_2d(2,i)),     &
                      start=(/ 1, 1, timeindex /), count = (/ f_dim_2d(1, i),  f_dim_2d(2,i), 1/)),    &
                 'write_cam_init','put_var '//trim(cflds(ifld)))
enddo

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

       temp_3d(1:f_dim_3d(1,i), 1:f_dim_3d(2,i), 1:f_dim_3d(3,i)) = &
   var%vars_3d(1:f_dim_3d(1,i), 1:f_dim_3d(2,i), 1:f_dim_3d(3,i),i)

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
!> Optionally returns the DART QTY of the variable.
!> 
!> @param[in]    index_in
!> The 'index' of a variable in the state vector, whose physical location 
!> and possibly variable kind are needed,
!>
!> @param[inout] location
!> The DART location_type location of the variable denoted by 'index'
!> 
!> @param[out]   var_kind
!> The optional argument which can return the DART QTY of the variable.


subroutine get_state_meta_data(index_in, location, var_kind)

! Given an integer index into the state vector structure, returns the
! associated location.
! The location may have components that are MISSING_R8 values, since some fields
! don't have locations in all three dimensions, i.e. PS has no vertical level,
! and other fiendish fields to be devised by parameterization studies may not
! have a longitude, or latitude.  The which_vert should take care of the vertical
! coordinate (it will be ignored), but the others will require more interesting  fixes.
! See order_state_fields for the QTY_s (and corresponding model variable names).

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_kind

integer  :: which_vert
integer  :: i, indx, index_1, index_2, index_3, nfld
integer  :: box, slice
logical  :: lfound

real(r8) :: lon_val, lat_val, lev_val
integer  :: ip, jp, kp, dom_id
integer  :: ndims

if (.not. module_initialized) call static_init_model()

lon_val = MISSING_R8
lat_val = MISSING_R8
lev_val = MISSING_R8

! get the state indices from dart index
! RMA-KR; Will this work for cubed sphere or other 'unstructured' grids?
!         I think so; ip, jp, and kp are interpreted according to the dimension
!         name in the coord_val calls, next.
call get_model_variable_indices(index_in, ip ,jp ,kp ,var_id=nfld, dom_id=dom_id)

! convert to lat, lon, lev coordinates
call coord_val(get_dim_name(dom_id,nfld,3), kp, lon_val, lat_val, lev_val)
call coord_val(get_dim_name(dom_id,nfld,2), jp, lon_val, lat_val, lev_val)
call coord_val(get_dim_name(dom_id,nfld,1), ip, lon_val, lat_val, lev_val)

ndims = get_num_dims(dom_id, nfld)

! RMA-KR; This will need to be changed for CAM-SE; 1d and 2d
if( ndims == 2 ) then
   which_vert = which_vert_2d(nfld)
else
   which_vert = which_vert_3d(nfld-state_num_2d)
endif

! This routine should error out for fields that have MISSING_R8 in lat_val or lon_val.
if (lon_val == MISSING_R8 .or. lat_val == MISSING_R8 ) then
   write(string1, *) 'Field ',cflds(nfld),' has no lon or lat dimension.  ', &
         'What should be specified for it in the call to location?'
   call error_handler(E_ERR, 'get_state_meta_data', string1, source, revision, revdate)
else
   location = set_location(lon_val, lat_val, lev_val, which_vert)
endif

! If the type is wanted, return it
if (present(var_kind)) then
   ! used by call from assim_tools_mod:filter_assim, which wants the DART QTY_
   var_kind = cam_to_dart_qtys(nfld)
endif

end subroutine get_state_meta_data

!-----------------------------------------------------------------------
!>
!> Function get_model_size assigns the 'model_size' calculated in static_init_model
!> to the function result 'get_model_size'.

function get_model_size()

integer(i8) :: get_model_size

if (.not. module_initialized) call static_init_model()

get_model_size = model_size

end function get_model_size

!-----------------------------------------------------------------------
!>
!> Function shortest_time_between_assimilations assigns the 'Time_step_atmos' calculated in 
!> static_init_model to the function result 'shortest_time_between_assimilations'.

function shortest_time_between_assimilations()

! Returns the shortest time you want to ask the model to
! advance in a single step

type(time_type) :: shortest_time_between_assimilations

if (.not. module_initialized) call static_init_model()

shortest_time_between_assimilations =  Time_step_atmos

end function shortest_time_between_assimilations

!-----------------------------------------------------------------------
!>
!> nc_write_model_atts
!> writes the model-specific attributes to a netCDF file.
!> 
!> @param[in] ncid      
!>  netCDF file identifier
!>
!> @param[in] domain_id      
!>  domain identifier (CAM has only 1 domain).

subroutine nc_write_model_atts( ncid, domain_id ) 


integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer :: n_dims, n_vars, n_attribs, unlimited_dim_ID
integer :: member_dim_ID, state_var_dim_ID, time_dim_ID,scalar_dim_ID
integer :: x_var_ID,state_var_ID, state_var_var_ID
! Add 1 to num_dims, for P0.
! This hard-wiring should be replaced if more D0 'coordinates' are added.
integer :: P_id(num_dims+1)
integer :: i, ifld, dim_id, g_id
integer :: grid_id(grid_num_1d)

if (.not. module_initialized) call static_init_model()

! Write Global Attributes

call nc_redef(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source)
call nc_add_global_attribute(ncid, "model_revision", revision)
call nc_add_global_attribute(ncid, "model_revdate", revdate)

call nc_add_global_attribute(ncid, "model", "CAM")

! Define the new dimensions IDs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! They have different dimids for this file than they had for caminput.nc
! P_id serves as a map between the 2 sets.
if (print_details .and. output_task0) then
   write(string1,*) 'num_dims = ',num_dims
   write(string2,*) ' dimens,       name,  size, cam dim_id, P[oste]rior id'
   call error_handler(E_MSG, 'nc_write_model_atts', string1,source,revision,revdate, text2=string2)
endif

! P_id debug
! This loops over the number of DIMENSIONS/COORDINATES on the file, not including P0.
! So P_id needs to be defined for P0 after this loop.
do i = 1,num_dims
   if (trim(dim_names(i)) /= 'time')  then
      call nc_check(nf90_def_dim(ncid, name=trim(dim_names(i)), len=dim_sizes(i),  &
                    dimid=P_id(i)), 'nc_write_model_atts','def_dim '//trim(dim_names(i)))
   else
     ! time, not P0
     P_id(i) = 0
   endif
   if (print_details .and. output_task0) then
      write(string1,'(I5,1X,A13,1X,2(I7,2X))') i,trim(dim_names(i)),dim_sizes(i), P_id(num_dims)
      call error_handler(E_MSG, 'nc_write_model_atts', string1,source,revision,revdate)
   endif
enddo

call nc_check(nf90_def_dim(ncid, name="scalar",   len=1, dimid=scalar_dim_ID) &
             ,'nc_write_model_atts', 'def_dim scalar')
call nc_check(nf90_def_dim(ncid, name="P0",   len=1, dimid=P_id(num_dims+1)) &
             ,'nc_write_model_atts', 'def_dim scalar')
if (print_details .and. output_task0) then
   write(string1,'(I5,1X,A13,1X,2(I7,2X))') i,'P0',P0%length, P_id(i)
   call error_handler(E_MSG, 'nc_write_model_atts', string1,source,revision,revdate)
endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Create the (empty) Coordinate Variables and their attributes
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! grid longitudes, latitudes, levels, and other coordinates.
! grid_id() is filled here; it's the dimid of the desired coordinate *on this P_Diag.nc file*.
! It's used to write coordinates.  ! There's some overlap of names, unfortunately.
! The argument after the 'xxx    ' label is a structure with all the relevant info in it.
! The structures are defined in "Grid fields" and filled by calls to create_grid_1d_instance
! in read_cam_coord.

grid_id = MISSING_I

if (lon%label /= ' ')  then
   dim_id = P_id(lon%dim_id)
   g_id   = find_name('lon',grid_names_1d)
   call write_cam_coord_def(ncid,'lon',lon , dim_id, grid_id(g_id))
endif
if (lat%label /= ' ')  then
   dim_id = P_id(lat%dim_id)
   g_id   = find_name('lat',grid_names_1d)
   call write_cam_coord_def(ncid,'lat',lat , dim_id, grid_id(g_id))
endif
if (lev%label /= ' ')  then
   dim_id = P_id(lev%dim_id)
   g_id   = find_name('lev',grid_names_1d)
   call write_cam_coord_def(ncid,'lev',lev , dim_id, grid_id(g_id))
! Gaussian weights -- because they're there.
endif
if (gw%label /= ' ')  then
   dim_id = P_id(gw%dim_id)
   g_id   = find_name('gw',grid_names_1d)
   call write_cam_coord_def(ncid,'gw',gw  , dim_id, grid_id(g_id))
! Hybrid grid level coefficients, parameters
endif
if (hyam%label /= ' ')  then
   dim_id = P_id(hyam%dim_id)
   g_id   = find_name('hyam',grid_names_1d)
   call write_cam_coord_def(ncid,'hyam',hyam, dim_id, grid_id(g_id))
endif
if (hybm%label /= ' ')  then
   dim_id = P_id(hybm%dim_id)
   g_id   = find_name('hybm',grid_names_1d)
   call write_cam_coord_def(ncid,'hybm',hybm, dim_id, grid_id(g_id))
endif
if (hyai%label /= ' ')  then
   dim_id = P_id(hyai%dim_id)
   g_id   = find_name('hyai',grid_names_1d)
   call write_cam_coord_def(ncid,'hyai',hyai, dim_id, grid_id(g_id))
endif
if (hybi%label /= ' ')  then
   dim_id = P_id(hybi%dim_id)
   g_id   = find_name('hybi',grid_names_1d)
   call write_cam_coord_def(ncid,'hybi',hybi, dim_id, grid_id(g_id))
endif
if (slon%label /= ' ')  then
   dim_id = P_id(slon%dim_id)
   g_id   = find_name('slon',grid_names_1d)
   call write_cam_coord_def(ncid,'slon',slon, dim_id, grid_id(g_id))
endif
if (slat%label /= ' ')  then
   dim_id = P_id(slat%dim_id)
   g_id   = find_name('slat',grid_names_1d)
   call write_cam_coord_def(ncid,'slat',slat, dim_id, grid_id(g_id))
endif
if (ilev%label /= ' ')  then
   dim_id = P_id(ilev%dim_id)
   g_id   = find_name('ilev',grid_names_1d)
   call write_cam_coord_def(ncid,'ilev',ilev, dim_id, grid_id(g_id))
endif
if (P0%label /= ' ')  then
   dim_id = P_id(num_dims+1)
   ! At some point, replace the kluge of putting P0 in with 'coordinates' 
   ! by defining grid_0d_kind, etc.
   g_id   = find_name('P0',grid_names_1d)
   call write_cam_coord_def(ncid,'P0',P0  , dim_id, grid_id(g_id))
endif

if (print_details .and. output_task0) then
   write(string1,*) '1d field#, grid_id, grid_names_1d'
   call error_handler(E_MSG, 'nc_write_model_atts', string1,source,revision,revdate)
   do i=1,grid_num_1d
      write(string1,*) 'grid_ = ', i, grid_id(i), trim(grid_names_1d(i))
      call error_handler(E_MSG, 'nc_write_model_atts', string1,source,revision,revdate)
   enddo
endif

! Leave define mode so we can fill variables
call nc_enddef(ncid)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Fill the coordinate variables
! Each 'vals' vector has been dimensioned to the right size for its coordinate.
! The default values of 'start' and 'count'  write out the whole thing.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (lon%label  /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('lon',grid_names_1d)),  lon%vals) &
                 ,'nc_write_model_atts', 'put_var lon')
if (lat%label  /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('lat',grid_names_1d)),  lat%vals) &
                 ,'nc_write_model_atts', 'put_var lat')
if (lev%label  /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('lev',grid_names_1d)),  lev%vals) &
                 ,'nc_write_model_atts', 'put_var lev')
if (gw%label   /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('gw',grid_names_1d)),   gw%vals) &
                 ,'nc_write_model_atts', 'put_var gw')
if (hyam%label /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('hyam',grid_names_1d)), hyam%vals) &
                 ,'nc_write_model_atts', 'put_var hyam')
if (hybm%label /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('hybm',grid_names_1d)), hybm%vals) &
                 ,'nc_write_model_atts', 'put_var hybm')
if (hyai%label /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('hyai',grid_names_1d)), hyai%vals) &
                 ,'nc_write_model_atts', 'put_var hyai')
if (hybi%label /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('hybi',grid_names_1d)), hybi%vals) &
                 ,'nc_write_model_atts', 'put_var hybi')
if (slon%label /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('slon',grid_names_1d)), slon%vals) &
                 ,'nc_write_model_atts', 'put_var slon')
if (slat%label /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('slat',grid_names_1d)), slat%vals) &
                 ,'nc_write_model_atts', 'put_var slat')
if (ilev%label /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('ilev',grid_names_1d)), ilev%vals) &
                 ,'nc_write_model_atts', 'put_var ilev')
if (P0%label   /= ' ') &
    call nc_check(nf90_put_var(ncid, grid_id(find_name('P0',grid_names_1d)),   P0%vals) &
                 ,'nc_write_model_atts', 'put_var P0')

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Flush the buffer and leave netCDF file open
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call nc_sync(ncid)

end subroutine nc_write_model_atts

! End of Module I/O

!#######################################################################

! model_interpolate section

!-----------------------------------------------------------------------
!>
!> Subroutine model_interpolate
!> Interpolates the provided state vector (on model grid points) to an arbitrary
!> location in the atmosphere (e.g. where an observation is).
!> 
!> @param[in]    state_handle
!> The DART ensemble_type structure which gives access to the ensemble of model states.
!>
!> @param[in] :: ens_size
!> The size of the ensemble.
!> 
!> @param[in] :: location
!> The DART location_type 'location' of the desired state estimate.
!> 
!> @param[in] :: obs_kind
!> The DART QTY of the variable being estimated.
!> 
!> @param[out] :: expected_obs
!> The ensemble state estimate of the 'obs_kind' at 'location'.
!> 
!> @param[out] :: istatus
!> A flag to signal the success of the interpolation.


subroutine model_interpolate(state_handle, ens_size, location, obs_kind, expected_obs, istatus)

! This subroutine is now a short routine that calls
! either a rectangular grid version for eul/FV
! or non-rectangular for cubed-sphere code.
! This does get QTYs from filter, not specific obs TYPEs.

! Model_interpolate must return a positive value for istatus for a failure.
! 0 means success, negative values are reserved for DART internal use.

type(ensemble_type),   intent(in)  :: state_handle
integer,               intent(in)  :: ens_size
type(location_type),   intent(in)  :: location ! The DART location_type 'location' of the desired state estimate.
integer,               intent(in)  :: obs_kind ! The DART QTY of the variable being estimated.
real(r8),              intent(out) :: expected_obs(ens_size) ! The state estimate of the 'obs_kind' at 'location'
integer,               intent(out) :: istatus(ens_size) ! A flag to signal the success of the interpolation.

! FIXME; In future DARTs it may be useful to return the DART QTY too.
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

call interp_lonlat(state_handle, ens_size, location, obs_kind, expected_obs, istatus)

end subroutine model_interpolate

!-----------------------------------------------------------------------

recursive subroutine interp_lonlat(state_handle, ens_size, obs_loc, obs_kind, interp_val, istatus)

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

type(ensemble_type),   intent(in)  :: state_handle
integer,               intent(in)  :: ens_size
type(location_type),   intent(in)  :: obs_loc
integer,               intent(in)  :: obs_kind
real(r8),              intent(out) :: interp_val(ens_size)
integer,               intent(out) :: istatus(ens_size)

! FIXME; In future DARTs it may be useful to return the DART QTY too.
!        also convert to a field name (DART subroutine (get_raw_...?)).

integer  :: i, vstatus(ens_size), cur_vstatus(ens_size)
real(r8) :: bot_lon, top_lon, delta_lon,                                &
            lon_below, lat_below, lat_above, lev_below,                 &
            lon_fract, lat_fract, temp_lon, a(ens_size, 2),           &
            lon_lat_lev(3), convert_units
real(r8), dimension(ens_size) :: val_11, val_12, val_21, val_22

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
istatus(:)     = 1
cur_vstatus(:) = 1
vstatus(:)     = 0 ! Assume good so you can keep track of vstatus
val_11(:)      = MISSING_R8
val_12(:)      = MISSING_R8
val_21(:)      = MISSING_R8
val_22(:)      = MISSING_R8
interp_val(:)  = MISSING_R8
! Get the observation (horizontal) position, in degrees
lon_lat_lev = get_location(obs_loc)

! Check whether model_mod can interpolate the requested variable.
! Pressure (3d) can't be specified as a state vector field (so s_type will = MISSING_I),
! but can be calculated for CAM, so obs_kind = QTY_PRESSURE is acceptable.
! obs_kind truly is a DART QTY variable, generally passed from
! obs_def/obs_def_XXX.f90: call interpolate.
! HK I think s_type is the index in cflds
! RMA-KR; use a new mechanism to define s_type (as in 'clm_variables')
!   > > > Just loop through cflds until I find it.
!            Need the state_name of this obs_kind
!         CLM uses obs_kind.  Is there a 1 to 1 match of CAM variables and DART QTYs?
!    > > >It requires hard-wiring all of the potential QTYs in the 'select case (obs_kind)' structure.
!         Could still have dart_to_cam?
!         Does paradigm of separating vars into 0d, 1d, 2d, and 3d make sense?
s_type = find_name(dart_to_cam_types(obs_kind),cflds)

if (s_type == MISSING_I) then
   if (obs_kind /= QTY_PRESSURE .and. obs_kind /= QTY_SURFACE_ELEVATION) then
      write(string1,*) 'Wrong type of obs = ', obs_kind
      call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
      return
   else
      ! CAM-chem 5))
      !          This will be used when interp_val is calculated,
      !          but define it here, as soon as it can be.
      ! Define for the non-chemical, non-state QTYs.
      convert_units = 1.0_r8
   endif
else
   ! CAM-chem Define it here for state variables
   convert_units = convert_mmr2vmr(s_type) 
endif

! Get lon and lat dimension names.

! Set [lon,lat,lev] names to a default, which will be overwritten for variables
! in the state vector, but not for other acceptable variables (3D pressure, surface
! elevation, ...?)
lon_name = 'lon'
lat_name = 'lat'
if (obs_kind == QTY_SURFACE_ELEVATION) then
   lev_name = 'none'
elseif (obs_kind == QTY_PRESSURE) then
   lev_name = 'lev'
! else
!    set below
endif

! Need to get lon, lat, lev dimension names for this field


! DART can't handle any 0d or 1d ob fields, so lump them together for elimination in this search.
s_type_01d = state_num_0d + state_num_1d
s_type_2d = s_type - s_type_01d
s_type_3d = s_type_2d - state_num_2d

! HK This if statement is just finding the rank of the variable (0D, 1D, 2D, 3D).
if (s_type == MISSING_I .and. &
   (obs_kind == QTY_PRESSURE) .or.  (obs_kind == QTY_SURFACE_ELEVATION)) then
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
   lon_name = get_lon_name(s_type)
   lat_name = get_lat_name(s_type)
   lev_name = 'none'
elseif (s_type_3d > 0 .and. s_type_3d <= state_num_3d) then
   lon_name = get_lon_name(s_type)
   lat_name = get_lat_name(s_type)
   lev_name = get_lev_name(s_type)
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
if (obs_kind == QTY_SURFACE_ELEVATION) then
   ! Acceptable field that's not in the state vector, same across the ensemble
   ! convert from geopotential height to real height in meters
   val_11(:) = phis(lon_ind_below, lat_ind_below) / gravity_const
   val_12(:) = phis(lon_ind_below, lat_ind_above) / gravity_const
   val_21(:) = phis(lon_ind_above, lat_ind_below) / gravity_const
   val_22(:) = phis(lon_ind_above, lat_ind_above) / gravity_const
   if (val_11(1) == MISSING_R8 .or. &
       val_12(1) == MISSING_R8 .or. &
       val_21(1) == MISSING_R8 .or. &
       val_22(1) == MISSING_R8 ) then
      vstatus(:) = 1 
      write(string1,*) 'interp_lonlat: val_##(mem1) = MISSING_R* for ',&
                       'lon, lat near ',lon_ind_above, lat_ind_above
      call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
   endif

elseif (is_vertical(obs_loc, "LEVEL")) then
   ! Pobs
   ! FIXME; I may want to change get_val_level to accept REAL level, not INT.
   !        What's the benefit?
   !        But it would be inconsistent with lon_ and lat_ indices,
   !           and I'd have to create an integer level anyway.
   !        May also want to handle staggered vertical grid (ilev).
   call get_val_level(state_handle, ens_size, lon_ind_below, lat_ind_below, nint(lon_lat_lev(3)), obs_kind, val_11, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val_level(state_handle, ens_size, lon_ind_below, lat_ind_above, nint(lon_lat_lev(3)), obs_kind, val_12, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val_level(state_handle, ens_size, lon_ind_above, lat_ind_below, nint(lon_lat_lev(3)), obs_kind, val_21, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val_level(state_handle, ens_size, lon_ind_above, lat_ind_above, nint(lon_lat_lev(3)), obs_kind, val_22, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)


elseif (is_vertical(obs_loc, "PRESSURE")) then
   call get_val_pressure(state_handle, ens_size,lon_ind_below,lat_ind_below,lon_lat_lev(3),obs_kind,val_11,cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val_pressure(state_handle, ens_size,lon_ind_below,lat_ind_above,lon_lat_lev(3),obs_kind,val_12,cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val_pressure(state_handle, ens_size,lon_ind_above,lat_ind_below,lon_lat_lev(3),obs_kind,val_21,cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val_pressure(state_handle, ens_size,lon_ind_above,lat_ind_above,lon_lat_lev(3),obs_kind,val_22,cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)

elseif (is_vertical(obs_loc, "HEIGHT")) then
   call get_val_height(state_handle, ens_size, lon_ind_below, lat_ind_below, lon_lat_lev(3), obs_loc, obs_kind, val_11, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val_height(state_handle, ens_size, lon_ind_below, lat_ind_above, lon_lat_lev(3), obs_loc, obs_kind, val_12, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val_height(state_handle,  ens_size,lon_ind_above, lat_ind_below, lon_lat_lev(3), obs_loc, obs_kind, val_21, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val_height(state_handle,  ens_size,lon_ind_above, lat_ind_above, lon_lat_lev(3), obs_loc, obs_kind, val_22, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)

elseif (is_vertical(obs_loc, "SURFACE")) then
   ! The 'lev' argument is set to 1 because there is no level for these types, and 'lev' will be
   ! ignored.
   call get_val(state_handle, ens_size, lon_ind_below, lat_ind_below, 1, obs_kind, val_11, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val(state_handle, ens_size, lon_ind_below, lat_ind_above, 1, obs_kind, val_12, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val(state_handle, ens_size, lon_ind_above, lat_ind_below, 1, obs_kind, val_21, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val(state_handle, ens_size, lon_ind_above, lat_ind_above, 1, obs_kind, val_22, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)

! This needs to be at the end of the block.  Otherwise, it short circuits GPS
! which asks for pressures on heights.
! elseif (obs_kind == QTY_PRESSURE) then
!    ! Calculate pressures from surface pressures and A and B coeffs.
!     write(string1,'(A)') 'No code available yet for obs_kind = QTY_PRESSURE '
!     call error_handler(E_ERR, 'interp_lon_lat', string1)

elseif (is_vertical(obs_loc, "SCALE_HEIGHT")) then
   ! Need option for this case
   write(string1,*)'Scale height is not an acceptable vert coord yet.  Skipping observation'
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
   return

! Need option for is_vertical("UNDEFINED")
else
   write(string1,*) '   No vert option chosen!'
   call error_handler(E_WARN, 'interp_lonlat', string1,source,revision,revdate)
   return

endif

! Conundrum (but unimportant for now): an ob could be excluded for > 1 reason.
! E.g. it's too far north and it's above the highest_obs_pressure_Pa.
! What istatus to return? a 2 (or more) digit number?  Like vstatus*10 + 4?
! RMA-KR; Note that there's no early return based on an interpolation failure.
!         The interpolation is done for those members for whom it's possible
!         and the others get 'failed' istatus, which is returned to the calling routine.

if (abs(lon_lat_lev(2)) > max_obs_lat_degree) then
   ! Define istatus to be a combination of vstatus (= 0 or 2 (for higher than highest_obs...))
   ! and whether the ob is poleward of the limits set in the namelist (+ 4).
   ! Too confusing for now; 
   ! istatus(:) = 10*vstatus + 4
   istatus(:) = 2
else
   istatus(:) = vstatus(:)
endif

where (istatus == 0 .or. istatus == 2) ! These are success codes
   ! indices of vals are (longitude, latitude)
   a(:, 1) = lon_fract * val_21 + (1.0_r8 - lon_fract) * val_11
   a(:, 2) = lon_fract * val_22 + (1.0_r8 - lon_fract) * val_12

   ! CAM-chem 6)); multiply the result by the unit conversion factor
   interp_val(:) = (lat_fract * a(:, 2) + (1.0_r8 - lat_fract) * a(:, 1)) * convert_units
endwhere

end subroutine interp_lonlat

!-----------------------------------------------------------------------

! Pobs
subroutine get_val_level(state_handle, ens_size, lon_index, lat_index, level, obs_kind, val, istatus)

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

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
integer,             intent(in)  :: level
integer,             intent(in)  :: obs_kind
real(r8),            intent(out) :: val(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer  :: vstatus(ens_size), i, indx
real(r8) :: p_surf(ens_size), threshold
integer  :: imem
real(r8), allocatable :: p_col(:)

! Start with failure condition
istatus(:) = 1
vstatus(:) = 1
val(:) = MISSING_R8

! This assumes that all variables are defined on model levels, not on interface levels.
! Exclude obs below the model's lowest level and above the highest level,
! but go ahead with surface fields (level = no_lev).
if (level /= no_lev .and. (level > lev%length .or. level < 1)) return
allocate(p_col(lev%length))

! Interpolate in vertical to get two bounding levels, but treat pressure
! specially since it has to be computed from PS instead of interpolated.

if (obs_kind == QTY_PRESSURE) then

   ! p_surf is returned in pascals, which is the right units for plevs_cam() below.
   ! RMA-KR; level is irrelevant for PS, and should not cause a failure even now that 
   !         io/state_structure_mod.f90:get_dart_vector_index is the eventual recipient of that index.
   !         Only lon and lat dimensions will be used to find the index into the state vector;
   !         'level' will not be used.  Same for the pre-RMA trunk version.
   call get_val(state_handle, ens_size, lon_index, lat_index, no_lev, QTY_SURFACE_PRESSURE, p_surf, vstatus)
   if (all(vstatus /= 0)) then
      deallocate(p_col)
      return
   endif
   ! Next, get the values on the levels for this PS.
   do imem = 1, ens_size
      if (vstatus(imem) == 0) then
         call plevs_cam (p_surf(imem), lev%length, p_col)
         val(imem) = p_col(level)
      endif
   enddo

else

   call get_val(state_handle, ens_size, lon_index, lat_index, level, obs_kind, val, vstatus)

endif

! if this routine is called with a location that has a vertical level above
! the pressure cutoff, go ahead and compute the value but return an istatus=2
! (unless some other error occurs later in this routine).  note that smaller
! level numbers are higher up in the atmosphere; level 1 is at the top.

if (level < highest_obs_level) then
   istatus(:) = 2
else
   istatus(:) = vstatus
endif

deallocate(p_col)

end subroutine get_val_level

!-----------------------------------------------------------------------

subroutine get_val_pressure(state_handle, ens_size, lon_index, lat_index, pressure, obs_qty, val, istatus)

! Gets the vertically interpolated value on pressure for variable obs_qty
! at lon_index, lat_index horizontal grid point
!
! This routine indicates things with the return code:
!   istatus 0 - success
!   istatus 1 - failure (e.g. above or below highest/lowest level, or can't
!                          interpolate the value)
!   istatus 2 - val is set successfully, but vert is above highest_obs_pressure
!
! Excludes observations below lowest level pressure and above highest level pressure.

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
real(r8),            intent(in)  :: pressure
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
integer,             intent(in)  :: obs_qty
real(r8),            intent(out) :: val(ens_size)
integer,             intent(out) :: istatus(ens_size)

real(r8), dimension(ens_size) :: bot_val, top_val, p_surf, frac
real(r8), dimension(ens_size) :: ps_local(ens_size, 2)
integer,  dimension(ens_size) :: top_lev, bot_lev, vstatus, cur_vstatus
! RMA-KR; cur_vstatus was explicitly dimensioned (ens_size), which was redundant.
integer               :: fld_index
integer(i8)           :: i, imem
real(r8), allocatable :: p_col(:,:)

! Start with error condition.
istatus(:) = 1
cur_vstatus(:) = 1
vstatus(:) = 0 ! so you can track statuses
val(:)     = MISSING_R8
p_surf(:) = MISSING_R8

! Need to get the surface pressure at this point.
! Find out whether the observed field is a staggered field in CAM.
! Check whether the state vector has wind components on staggered grids, i.e. whether CAM is FV.
! find_name returns 0 if the field name is not found in the cflds list.

fld_index   = find_name('PS',cflds)
i = index_from_grid(1,lon_index,lat_index,  fld_index)
ps_local(:, 1) = get_state(i, state_handle)

if (obs_qty == QTY_U_WIND_COMPONENT .and. find_name('US', cflds) /= 0) then
   ! ps defined on lat grid (-90...90, nlat = nslat + 1),
   !    need it on staggered lat grid, which starts half a grid spacing north.

   i = index_from_grid(1,lon_index,lat_index+1,fld_index)
   ps_local(:, 2) = get_state(i, state_handle)
   p_surf(:)      = (ps_local(:, 1) + ps_local(:, 2))* 0.5_r8
elseif (obs_qty == QTY_V_WIND_COMPONENT .and. find_name('VS', cflds) /= 0) then
   ! lon =     0...     255 (for 5 degree grid)
   !slon = -2.5 ... 252.5
   if (lon_index == slon%length) then
      i = index_from_grid(1,1,          lat_index ,fld_index)
   else
      i = index_from_grid(1,lon_index+1,lat_index ,fld_index)
   endif
   ps_local(:, 2) = get_state(i, state_handle)
   p_surf(:)     = (ps_local(:, 1) + ps_local(:, 2))* 0.5_r8
else
   ! A-grid ps can be retrieved from state vector, which was used to define ps on entry to
   ! model_interpolate.
   p_surf(:)     = ps_local(:, 1)
endif

! Next, get the pressures on the levels for this ps
! Assuming we'll only need pressures on model mid-point levels, not interface levels.
! This pressure column will be for the correct grid for obs_qty, since p_surf was taken
!     from the grid-correct ps[_stagr] grid
allocate(p_col(lev%length, ens_size))
p_col(:,:) = MISSING_R8
do imem = 1, ens_size
   call plevs_cam(p_surf(imem), lev%length, p_col(:, imem))
enddo

do imem = 1, ens_size
   if (pressure <= p_col(1, imem) .or. pressure >= p_col(lev%length, imem)) then
      vstatus(imem) = 1
      ! Exclude obs below the model's lowest level and above the highest level
      ! We *could* possibly use ps and p(lev%length) to interpolate for points below the lowest level.
      !return
   endif
enddo

! Interpolate in vertical to get two bounding levels

! Search down through pressures for each ensemble member
do imem = 1, ens_size
   if (vstatus(imem) == 0) then
      levloop: do i = 2, lev%length
         if (pressure < p_col(i, imem)) then
            top_lev(imem) = i -1
            bot_lev(imem) = i
            frac(imem) = (p_col(i, imem) - pressure) / &
                  (p_col(i, imem) - p_col(i - 1, imem))
            exit levloop
         endif
      enddo levloop
   else
      ! This is to avoid top_lev and bot_lev getting nonsense values
      top_lev(imem) = 1
      bot_lev(imem) = 2
   endif
enddo

if (obs_qty == QTY_PRESSURE) then
   ! can't get pressure on levels from state vector; get it from p_col instead
   do imem = 1, ens_size
      bot_val(imem) = p_col(bot_lev(imem), imem)
      top_val(imem) = p_col(top_lev(imem), imem)
   enddo
else
   call get_val_array_of_levels(state_handle, ens_size, lon_index, lat_index, bot_lev, obs_qty, bot_val, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val_array_of_levels(state_handle, ens_size, lon_index, lat_index, top_lev, obs_qty, top_val, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
endif


! Failed to get value for interpolation; return istatus = 1
where (vstatus == 0)
   istatus = 0
   val = (1.0_r8 - frac) * bot_val + frac * top_val
elsewhere
   istatus = 1
   val = MISSING_R8
endwhere

! if this routine is called with a location that has a vertical pressure above
! the pressure cutoff, go ahead and compute the value but return an istatus=2
! (unless some other error occurs later in this routine).
if (pressure < highest_obs_pressure_Pa) then
   where (istatus == 0) istatus = 2
endif

deallocate(p_col)

end subroutine get_val_pressure

!-----------------------------------------------------------------------

subroutine get_val_height(state_handle, ens_size, lon_index, lat_index, height, location, obs_kind, val, istatus)

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

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
real(r8),            intent(in)  :: height
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
real(r8),            intent(out) :: val(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer     :: i, fld_index
integer,  dimension(ens_size) :: top_lev, bot_lev, vstatus, cur_vstatus
real(r8), dimension(ens_size) :: bot_val, top_val, frac
integer(i8) :: ind
real(r8)    :: p_surf(ens_size), ps_local(ens_size, 2)
logical     :: stagr_lon, stagr_lat
real(r8), allocatable :: p_col(:, :), model_h(:, :) !(lev%length, ens_size)
integer :: imem

! Start with error condition.
! RMA-KR; should vstatus start with 1?  Then change comment to 'start with error condition'.
!         vstatus is first passed to model_heights, which sets it to 1, so this is irrelevant.
istatus(:)   = 1
vstatus(:)   = 1
cur_vstatus(:) = 1
val(:)       = MISSING_R8
stagr_lon = .false.
stagr_lat = .false.

! Assuming we'll only need pressures on model mid-point levels, not interface levels.
allocate(p_col(lev%length, ens_size))
allocate(model_h(lev%length, ens_size))
! Need to get the surface pressure at this point.
! Check whether the state vector has wind components on staggered grids, i.e. whether CAM is FV.
! See get_val_pressure for more documentation.
fld_index   = find_name('PS',cflds)
ind         = index_from_grid(1,lon_index,lat_index,  fld_index)
ps_local(:, 1) = get_state(ind, state_handle)

! find_name returns 0 if the field name is not found in the cflds list.
if (obs_kind == QTY_U_WIND_COMPONENT .and. find_name('US', cflds) /= 0) then
   stagr_lat = .true.
   ind = index_from_grid(1,lon_index,lat_index+1,fld_index)
   ps_local(2, :) = get_state(ind, state_handle)
   p_surf(:) = (ps_local(1, :) + ps_local(2, :))* 0.5_r8
elseif (obs_kind == QTY_V_WIND_COMPONENT .and. find_name('VS', cflds) /= 0) then
   stagr_lon = .true.
   if (lon_index == slon%length) then
      ind = index_from_grid(1,1,          lat_index ,fld_index)
   else
      ind = index_from_grid(1,lon_index+1,lat_index ,fld_index)
   endif
   ps_local(:, 2) = get_state(ind, state_handle)
   p_surf(:) = (ps_local(:, 1) + ps_local(:, 2))* 0.5_r8
else
   p_surf(:) = ps_local(:, 1)
endif

! Next, get the heights on the levels for this ps

! We want to use the new vec for each new ob on height because the state was updated
! for all previous obs, and we want to use the most up to date state to get the best location.
! The above comment is untrue - the state is not updated, either it is the forward operator
! before assimilation, or it is the mean (not updated during assimilation)
call model_heights(state_handle, ens_size, lev%length, p_surf, location, model_h, vstatus)
if (all(vstatus == 1)) return    ! Failed to get model heights; return istatus = 1

! Exclude obs below the model's lowest level and above the highest level
do imem = 1, ens_size
  if (height >= model_h(1, imem) .or. height <= model_h(lev%length, imem)) vstatus(imem) = 1 ! Fail 
enddo

! ? Implement 3Dp here?  or should/can it not use the ens mean PS field?
do imem = 1, ens_size
   call plevs_cam(p_surf(imem), lev%length, p_col(:, imem))
enddo

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
! HK You have a highest_obs_height_m for each ensemble member. Is this what you want?
! HK The trunk will ens up with highest_obs_height_m equal to its first ensemble
! member at the first obseravation location.
!> @todo The location used in the distributed forward operator will be different 
!> on each task for the highest_obs_height_calculation
if (highest_obs_height_m == MISSING_R8) then
   ! Search until we find a good member
   memloop: do imem = 1, ens_size 
      if (vstatus(imem) == 0) then
         levloop: do i=2,lev%length
            if (p_col(i, imem) > highest_obs_pressure_Pa) then
               ! highest_obs_height_m = model_h(i)
               highest_obs_height_m = model_h(i, imem) + (model_h(i-1, imem)-model_h(i, imem))*  &
                                             ((p_col(i, imem)-highest_obs_pressure_Pa) / &
                                              (p_col(i, imem)-p_col(i-1, imem)))
               write(string1, *) 'highest_obs_height_m = ',highest_obs_height_m
               call error_handler(E_MSG,'get_val_height', string1, &
                                  source, revision, revdate)
               exit memloop
            endif
         enddo levloop
      endif
   enddo memloop
endif


! Interpolate in vertical to get two bounding levels.
! Search down through heights and set the enclosing level numbers
! and the fraction between them.  There has already been a test to
! ensure the height is between the levels (and has discarded values
! exactly equal to the limits), so this will always succeed.
do imem = 1, ens_size
   if (vstatus(imem) == 0) then
      lev2loop: do i = 2, lev%length
         if (height > model_h(i, imem)) then
            top_lev(imem) = i -1
            bot_lev(imem) = i
            frac(imem) = (model_h(i, imem) - height      ) / &
                  (model_h(i, imem) - model_h(i-1, imem))
            exit lev2loop
         endif
      enddo lev2loop
   else ! This is so you can make a call to get_val_array_of_levels without
        ! looking at the vstatus of each ensemble member.
      bot_lev(imem) = 2
      top_lev(imem) = 1
   endif
enddo

if (obs_kind == QTY_PRESSURE) then
   do imem = 1, ens_size
      bot_val(imem) = p_col(bot_lev(imem), imem)
      top_val(imem) = p_col(top_lev(imem), imem)
   enddo
else
   call get_val_array_of_levels(state_handle, ens_size, lon_index, lat_index, bot_lev, obs_kind, bot_val, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   call get_val_array_of_levels(state_handle, ens_size, lon_index, lat_index, top_lev, obs_kind, top_val, cur_vstatus)
   call update_vstatus(ens_size, cur_vstatus, vstatus)
   ! Failed to get a value to use in interpolation
   !if (vstatus == 1) return
endif

istatus(:) = vstatus(:)

where (istatus == 0) 
   val = (1.0_r8 - frac) * bot_val + frac * top_val
endwhere

if (height > highest_obs_height_m ) then
   ! if this routine is called with a location that has a vertical height above
   ! the pressure cutoff, pass back the value but return an istatus=2
   ! (Only for successful forward operators)
   where(istatus == 0) istatus = 2
endif


deallocate(p_col, model_h)

end subroutine get_val_height

!-----------------------------------------------------------------------

subroutine get_val(state_handle, ens_size, lon_index, lat_index, level, obs_kind, val, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
integer,             intent(in)  :: level
integer,             intent(in)  :: obs_kind
real(r8),            intent(out) :: val(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer(i8) :: indx
integer     :: field_type

! Start with error condition.
istatus(:) = 1
val(:) = MISSING_R8

field_type = find_name(dart_to_cam_types(obs_kind),cflds)
if (field_type <= 0 .or. field_type > nflds) return

indx = index_from_grid(level, lon_index, lat_index, field_type)
!> @todo pull this check out or error
! HK: This check is not correct for XCESM 
! RMA-KR; is this check related to synthetic obs, which have state indices < 0?
!if (indx > 0 .and. indx <= model_size) then
   istatus(:) = 0
   val = get_state(indx, state_handle)
!endif

end subroutine get_val


!-----------------------------------------------------------------------
!> Same as get val but accepts an array of levels
subroutine get_val_array_of_levels(state_handle, ens_size, lon_index, lat_index, levels, obs_kind, val, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
integer,             intent(in)  :: levels(ens_size)
integer,             intent(in)  :: obs_kind
real(r8),            intent(out) :: val(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer(i8) :: indx(ens_size)
integer     :: field_type
integer     :: imem

! Start with error condition.
istatus(:) = 1
val(:) = MISSING_R8

field_type = find_name(dart_to_cam_types(obs_kind),cflds)
if (field_type <= 0 .or. field_type > nflds) return

do imem = 1, ens_size
   indx(imem) = index_from_grid(levels(imem), lon_index, lat_index, field_type)
enddo
! HK: This check is not correct for XCESM.
! RMA-KR; is this check related to synthetic obs, which have state indices < 0?
!if (indx > 0 .and. indx <= model_size) then
   istatus(:) = 0
   call get_state_array(val, indx, state_handle)
!endif

end subroutine get_val_array_of_levels


!-----------------------------------------------------------------------

subroutine set_highest_obs_limit()

! Verify that the value for highest_obs_pressure_Pa in the namelist is ok.
!
! If this routine detects an error it calls the error handler with a
! fatal error.  If it returns, the namelist value is ok.
!
! Sets the module global variable 'highest_obs_level', and references
! the hybm array.

integer  :: i, lowest_ok
real(r8) :: p_surf, top
real(r8), allocatable :: p_col(:)
! This assumes that all variables are defined on model levels, not on interface levels.
allocate(p_col(lev%length))

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
call plevs_cam(P0%vals(1), lev%length, p_col)

! Loop downwards through pressure values (1 = model top, lev%length = bottom).
! The level is set to the highest level which is below the given threshold.


! RMA-KR; this lev%length condition was added to ensure that highest_obs_level
!         doesn't end up with value lev%length+1 due to the loop running all the way through.
High: do highest_obs_level=1,lev%length
   if (p_col(highest_obs_level) > highest_obs_pressure_Pa .or. &
       highest_obs_level == lev%length) exit High
enddo High

! Test whether user has set highest_obs_pressure_Pa to be uselessly small (high),
! which causes problems for setting highest_obs_height_m.
! highest model level (mid-layer) pressure:
top = hyam%vals(1)*P0%vals(1)
if (highest_obs_pressure_Pa < top) then
   write(string1, '(2A)') 'Namelist variable "highest_obs_pressure_Pa" is too small', &
                          ' (located above the model atmosphere)'
   write(string2, '(A,1pe15.5)') '   Reset to at least ',top
   call error_handler(E_ERR, 'set_highest_obs_limit', string1, source, revision, revdate, text2=string2)
endif

! Test to be sure user hasn't set level so low that contributions from
! terrain are going to be an issue.  If they do, tell them in the error
! message what the valid limit is.
if (hybm%vals(highest_obs_level) > 0.0_r8) then
   lowest_ok = 1
   findzero: do i=2, lev%length
      if (hybm%vals(i) > 0.0_r8) then
         lowest_ok = i-1
         exit findzero
      endif
   enddo findzero
   write(string1, '(A)') 'invalid value for namelist "highest_obs_pressure_Pa"'
   write(string2, '(A)') 'value is too large (located out of the pure pressure levels of the atmosphere)'
   write(string3, '(A,F9.3,A)') 'must specify a value located above pressure ', p_col(lowest_ok), ' Pascals'
   call error_handler(E_ERR, 'set_highest_obs_limit', string1, source, revision, revdate, &
                      text2=string2, text3=string3)
endif

deallocate(p_col)

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
   do i=1,f_dim_1d(1,nf)
      indx = indx + 1
      st_vec(indx) = var%vars_1d(i, nf)
   enddo
enddo

!  2d variables
do nf = 1, state_num_2d
   do j=1,f_dim_2d(2,nf)
   do i=1,f_dim_2d(1,nf)
      indx = indx + 1
      st_vec(indx) = var%vars_2d(i, j, nf)
   enddo
   enddo
enddo

!  3D fields
!  This section is only entered for models with logically rectangular grids,
!  which will have dimensions level, longitude, and latitude.
! RMA-KR; the indices in vars_3d are reversed compared to the non-RMA trunk,
!         due to no longer re-ordering dimensions to a standard order.
do nf= 1, state_num_3d
   do k=1,f_dim_3d(3,nf)
   do j=1,f_dim_3d(2,nf)
   do i=1,f_dim_3d(1,nf)
      indx = indx + 1
      st_vec(indx) = var%vars_3d(i, j, k, nf)
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
   do i=1,f_dim_1d(1, nf)
      indx = indx + 1
      var%vars_1d(i, nf) = st_vec(indx)
   enddo
enddo

!  2d fields
do nf = 1, state_num_2d
   do j = 1, f_dim_2d(2,nf)
   do i = 1, f_dim_2d(1,nf)
      indx = indx + 1
      var%vars_2d(i, j, nf) = st_vec(indx)
   enddo
   enddo
enddo

! 3D fields; see comments in prog_var_to_vect
! RMA-KR; the indices in vars_3d are reversed compared to the non-RMA trunk,
!         due to no longer re-ordering dimensions to a standard order.
do nf = 1, state_num_3d
   do k = 1, f_dim_3d(3,nf)
   do j = 1, f_dim_3d(2,nf)
   do i = 1, f_dim_3d(1,nf)
      indx = indx + 1
      var%vars_3d(i, j, k, nf) = st_vec(indx)
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
! get_close section

!-----------------------------------------------------------------------
!>
!> Subroutine get_close_obs
!> 
!> get_close_obs takes as input an "observation" location, a DART TYPE (not QTY),
!> and a list of all potentially close locations and QTYs on this task.
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
!> @param[in]    base_loc
!> The DART location_type location of the observation, which is the target of *get_close_obs*
!> 
!> @param[in]    base_type 
!> The DART TYPE (not QTY) of the observation
!> 
!> @param[inout] locs(:)
!> The DART location_type locations of the potentially close state variables
!> 
!> @param[in]    kinds(:)
!> The DART QTYs of the potentially close state variables
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
!>
!> @param[in]    state_handle
!> The DART ensemble_type structure which gives access to the ensemble of model states.

subroutine get_close_obs(filt_gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_indices, distances, state_handle)

type(get_close_type), intent(in)    :: filt_gc
type(location_type),  intent(in)    :: base_loc
integer,              intent(in)    :: base_type
type(location_type),  intent(inout) :: locs(:)
integer,              intent(in)    :: loc_qtys(:)
integer,              intent(in)    :: loc_types(:)
integer,              intent(out)   :: num_close
integer,              intent(out)   :: close_indices(:)
real(r8),             intent(out), optional :: distances(:)
type(ensemble_type),  intent(in),  optional :: state_handle

call get_close(filt_gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_indices, distances, state_handle)

end subroutine get_close_obs

!-----------------------------------------------------------------------

subroutine get_close_state(filt_gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                         num_close, close_indices, distances, state_handle)

type(get_close_type), intent(in)    :: filt_gc
type(location_type),  intent(in)    :: base_loc
integer,              intent(in)    :: base_type
type(location_type),  intent(inout) :: locs(:)
integer,              intent(in)    :: loc_qtys(:)
integer(i8),          intent(in)    :: loc_indx(:)
integer,              intent(out)   :: num_close
integer,              intent(out)   :: close_indices(:)
real(r8),             intent(out), optional :: distances(:)
type(ensemble_type),  intent(in),  optional :: state_handle

call get_close(filt_gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_indices, distances, state_handle)

end subroutine get_close_state

!-----------------------------------------------------------------------

subroutine get_close(filt_gc, base_loc, base_type, locs, loc_qtys, &
                     num_close, close_indices, distances, state_handle)

! get_close_obs takes as input an "observation" location, a DART TYPE (not QTY),
! and a list of all potentially close locations and QTYs on this task.
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
type(location_type),  intent(in)    :: base_loc
integer,              intent(in)    :: base_type
type(location_type),  intent(inout) :: locs(:)
integer,              intent(in)    :: loc_qtys(:)
integer,              intent(out)   :: num_close
integer,              intent(out)   :: close_indices(:)
real(r8),             intent(out), optional :: distances(:)
type(ensemble_type),  intent(in),  optional :: state_handle

! FIXME remove some (unused) variables?
integer  :: k, t_ind
integer  :: base_which, local_base_which, obs_which, local_obs_which
integer  :: base_obs_kind
real(r8) :: base_array(3), local_base_array(3), obs_array(3), local_obs_array(3)
real(r8) :: damping_dist, threshold, thresh_wght
type(location_type) :: local_base_loc, local_loc, vert_only_loc

if (.not. module_initialized) call static_init_model()

! If base_obs vert type is not pressure; convert it to pressure
base_which    = nint(query_location(base_loc))
base_array    = get_location(base_loc)
base_obs_kind = get_quantity_for_type_of_obs(base_type)

! Upgrading convert_vert to use field profiles at the actual ob location is
! probably not worthwhile: that approx horiz location of the obs is used only to
! convert its vert coord to pressure (if necessary),
! which, in turn, is used to modify the distance if the ob or model variable
! is higher than highest_XXX_Pa.  That modification tapers to 0,
! so any errors introduced by this approx will be continuous and random,
! introducing no bias.
if (base_which == VERTISPRESSURE .and. vert_coord == 'pressure') then
   local_base_loc = base_loc
   local_base_array   = get_location(base_loc)  ! needed in num_close loop
   local_base_which   = base_which
else
   call convert_vert(state_handle, base_array, base_which, base_loc, base_obs_kind, &
                     local_base_array, local_base_which)
   local_base_loc = set_location(base_array(1), base_array(2), local_base_array(3), &
                                     local_base_which)
endif

! Get all the potentially close obs but no distances. 
call loc_get_close(filt_gc, local_base_loc, base_type, locs, loc_qtys, &
                   num_close, close_indices)

do k = 1, num_close

   ! The indices in close_obs refer to the subset of (state) vars or obs ON 1 TASK.
   ! That subset is (re)labeled 1...num_vars_task#, where num_vars_task# ~ state_vec_size/ntasks.
   ! So those indices can't tell me which state vector element I'm dealing with.
   ! I need to use the location of each returned close_indices to learn anything about it.

   t_ind = close_indices(k)
   obs_array = get_location(locs(t_ind))
   ! query_location returns location%which_vert, if no 'attr' argument is given.
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
      call convert_vert(state_handle, obs_array, obs_which, locs(t_ind), loc_qtys(t_ind), &
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

!>@todo FIXME this should be removed and replaced by calls to obs_impact
!> in the assim_tools module.
!  allow a namelist specified kind string to restrict the impact of those
!  obs kinds to only other obs and state vars of the same kind.
   if ((impact_kind_index >= 0)                .and. &
       (impact_kind_index == base_obs_kind)    .and. &
       (impact_kind_index /= loc_qtys(t_ind))) then
      if(present(distances)) distances(k) = 999999.0_r8     ! arbitrary very large distance

   else
      ! Need to damp the influence of all obs (VERTISUNDEF, VERTISSURFACE too) on model state vars
      ! above highest_state_pressure_Pa.

      ! The which vert of local_base_loc determines how vertical distance to local_loc is calculated.
      ! It can be VERTISSCALEHEIGHT.
      if(present(distances)) distances(k) = get_dist(local_base_loc, local_loc, base_type, loc_qtys(t_ind))

      ! Damp the influence of obs, which are below the namelist variable highest_OBS_pressure_Pa,
      ! on variables above highest_STATE_pressure_Pa.
      ! This section could also change the distance based on the QTY_s of the base_obs and obs.

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
   
         if(present(distances)) distances(k) = distances(k) + damping_dist * damping_dist * damp_wght

      endif
   endif

enddo

end subroutine get_close

!-----------------------------------------------------------------------
!> wrapper for convert_vert so it can be called from assim_tools
!>
!> @param[in]    state_handle
!> The DART ensemble_type structure which gives access to the ensemble of model states.
!>
!> @param[inout] obs_loc
!> The DART location_type location of the observation.
!> 
!> @param[in]    obs_kind 
!> The DART QTY of the observation.
!>
!> @param[out]   vstatus 
!> The status of the conversion from one vertical location to another.
!>
!--------------------------------------------------------------------

subroutine convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, status)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:), loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: status(:)

real(r8) :: old_array(3)
integer  :: old_which, wanted_vert
type(location_type) :: old_loc

real(r8) :: new_array(3)
integer  :: new_which, i


status(:) = 0 ! I don't think cam has a return status for vertical conversion
wanted_vert = query_vert_localization_coord() 

do i=1, num
   old_which = query_location(locs(i), 'which_vert')
   if (old_which == wanted_vert) cycle

   old_loc = locs(i)
   old_array = get_location(locs(i))

   call convert_vert(state_handle, old_array, old_which, old_loc, loc_qtys(i), new_array, new_which)

   if(new_which == MISSING_I) then
      status(i) = 1
   else
      locs(i) = set_location(new_array(1), new_array(2), new_array(3), new_which)
   endif
enddo


end subroutine convert_vertical_obs

!--------------------------------------------------------------------
!>@todo FIXME there should be a more efficient way to convert
!>state locations - no interpolation in the horizontal is needed.

subroutine convert_vertical_state(state_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

real(r8) :: old_array(3)
integer  :: old_which, wanted_vert
type(location_type) :: old_loc

real(r8) :: new_array(3)
integer  :: new_which, i

wanted_vert = query_vert_localization_coord() 

do i=1, num
   old_which = query_location(locs(i), 'which_vert')
   if (old_which == wanted_vert) cycle

   old_loc = locs(i)
   old_array = get_location(locs(i))

   old_which = query_location(locs(i), 'which_vert')
   if (old_which == wanted_vert) cycle

   call convert_vert(state_handle, old_array, old_which, old_loc, loc_qtys(i), new_array, new_which)

   ! this is converting state locations.  it shouldn't fail.
   if(new_which == MISSING_I) then
      istatus = 1
      return
   else
      locs(i) = set_location(new_array(1), new_array(2), new_array(3), new_which)
   endif
enddo

istatus = 0

end subroutine convert_vertical_state


!-----------------------------------------------------------------------

subroutine convert_vert(state_handle, old_array, old_which, old_loc, old_kind, new_array, new_which)

! Uses model information and subroutines to convert the vertical location of an ob
! (prior, model state variable, or actual ob) into the standard vertical coordinate
! (pressure or log_invP = log(P0/ps)).
! Kevin Raeder 10/26/2006
! updated 2014 for WACCM use; log_invP vertical coordinate.

type(ensemble_type), intent(in)    :: state_handle
real(r8),            intent(in)    :: old_array(3)
integer,             intent(in)    :: old_which
type(location_type), intent(in)    :: old_loc
integer,             intent(in)    :: old_kind
real(r8),            intent(inout) :: new_array(3)
integer,             intent(out)   :: new_which

integer  :: top_lev, bot_lev
integer  :: istatus(1), closest
integer  :: lon_ind, lat_ind, cam_type
! p_surf dimensioned (1) because it's input to interp_lonlat, 
! which needs it to be an array because of RMA.
real(r8)              :: p_surf(1), frac, l, m, new_pressure
type(location_type)   :: temp_loc
integer               :: slon_index

character(len=8) :: cam_varname
integer :: ens_size ! To call interp_lonlat with ens_size of 1
real(r8), allocatable :: p_col(:)
real(r8), allocatable :: model_h(:)

ens_size = 1

!HK not building ps arrays.
slon_index = find_name('slon',dim_names)

! this code does not alter the lat/lon, only the vertical.
! but still return a full location for subsequent use.
new_array(1) = old_array(1)
new_array(2) = old_array(2)

! these should be set by the code below; it's an error if not.
new_which    = MISSING_I
new_array(3) = MISSING_R8
allocate(p_col(lev%length))

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
   cam_type = find_name(dart_to_cam_types(old_kind),cflds)
   if (cam_type < 0) then
      write(string1,*)'old_kind  is ',old_kind,' | cam_type is ',cam_type
      write(string2,*)'get_name_for_quantity of old_kind ', trim(get_name_for_quantity(old_kind))
      call error_handler(E_ERR,'convert_vert',string1,source,revision,revdate,text2=string2)
   endif

   ! Assumes 2D obs locations are (lon, lat) and 3D are (lev,lon,lat).

   ! Get the column of pressures at this location, from the ensemble mean.

   cam_varname = trim(cflds(cam_type))
   if (cam_varname == 'US') then
      call coord_index('lon', old_array(1), lon_ind)
      call coord_index('slat', old_array(2), lat_ind)
      !p_surf = ps_stagr_lat(lon_ind,lat_ind)
      p_surf = 0.5*(get_surface_pressure(state_handle, ens_size, lon_ind, lat_ind) + &
                    get_surface_pressure(state_handle, ens_size, lon_ind, lat_ind +1) )
      ! WHAT ABOUT FIELDS THAT MIGHT COME ON ilevS ?   have lev_which_dimid from above;
      !     test = ilev%dim_id or lev%dim_id
      call plevs_cam(p_surf(1), lev%length, p_col)
   elseif (cam_varname == 'VS') then
      call coord_index('slon', old_array(1), lon_ind)
      call coord_index('lat', old_array(2), lat_ind)
      !p_surf = ps_stagr_lon(lon_ind,lat_ind)
      if ( lon_ind == 1 ) then
         p_surf = 0.5*(get_surface_pressure(state_handle, ens_size, lon_ind,               lat_ind) + &
                       get_surface_pressure(state_handle, ens_size, dim_sizes(slon_index), lat_ind) )
      else
         p_surf = 0.5*(get_surface_pressure(state_handle, ens_size, lon_ind -1, lat_ind) + &
                       get_surface_pressure(state_handle, ens_size, lon_ind,    lat_ind) )
      endif
      call plevs_cam(p_surf(1), lev%length, p_col)
   else
      call coord_index('lon', old_array(1), lon_ind)
      call coord_index('lat', old_array(2), lat_ind)
      !p_surf = ps(lon_ind,lat_ind)
      p_surf = get_surface_pressure(state_handle, ens_size, lon_ind, lat_ind)
      call plevs_cam(p_surf(1), lev%length, p_col)
   endif
else
   ! Make a vertical location that has a vert type of surface.
   ! Don't need lon_lat_vert array because old_array is passed in,
   ! which is get_location(old_loc)
   temp_loc = set_location(old_array(1), old_array(2), 0.0_r8, VERTISSURFACE)
   ! Find ps at the ob point.  Need to interpolate.
   ! Only interested in P (columns), so don't need to worry about staggered grids here.
   call interp_lonlat(state_handle, ens_size, temp_loc, QTY_SURFACE_PRESSURE, p_surf, istatus)
   if (istatus(1) == 1) then
      write(string1,'(A,I8)') 'interp_X failed for QTY_SURFACE_PRESSURE.'
      call write_location(0, old_loc, charstring=string2)
      call error_handler(E_ERR, 'convert_vert', string1,source,revision,revdate, text2=string2)
   endif

   call plevs_cam(p_surf(1), lev%length, p_col)

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
      new_array(3) =  p_surf(1)
      new_which = VERTISPRESSURE
   elseif (vert_coord == 'log_invP') then
    ! Scale height at the surface is 0.0_r8 by definition [log(p_surf/p_surf)]
      new_array(3) = 0.0_r8
      new_which = VERTISSCALEHEIGHT
   endif

elseif (old_which == VERTISPRESSURE) then
   if (vert_coord == 'pressure') then
      new_array(3) =  old_array(3)
      new_which = VERTISPRESSURE
   elseif (vert_coord == 'log_invP') then
      new_array(3) = scale_height(p_surface=p_surf(1), p_above=old_array(3))
      new_which = VERTISSCALEHEIGHT
   endif

elseif (old_which == VERTISSCALEHEIGHT) then
   if (vert_coord == 'pressure') then
      new_array(3) = p_surf(1) / exp(old_array(3))
      new_which = VERTISPRESSURE
   elseif (vert_coord == 'log_invP') then
      new_array(3) = old_array(3)
      new_which = old_which
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
      new_array(3) = scale_height(p_surface=p_surf(1), p_above=p_col(nint(old_array(3))))
      new_which = VERTISSCALEHEIGHT
   endif

elseif (old_which == VERTISHEIGHT) then

   allocate(model_h(lev%length))
   call model_heights(state_handle, ens_size, lev%length, p_surf(1), old_loc,  model_h, istatus(1))
   if (istatus(1) == 1) then
      write(string1, *) 'model_heights failed'
      call error_handler(E_ERR, 'convert_vert', string1)
      ! return
   endif

   ! Search down through heights
   ! This assumes linear relationship of pressure to height over each model layer,
   ! when really it's exponential.  How bad is that?
!   bot_lev = 2
!   do while (old_array(3) <= model_h(bot_lev) .and. bot_lev <= lev%length)
!      bot_lev = bot_lev + 1
!   end do
   Bottom: do bot_lev = 2,lev%length
      if (old_array(3) > model_h(bot_lev) .or. bot_lev == lev%length) exit Bottom
   end do Bottom
   if (bot_lev > lev%length) bot_lev = lev%length
   top_lev = bot_lev - 1

   ! Write warning message if not found within model level heights.
   ! Maybe this should return failure somehow?
   if (top_lev == 1 .and. old_array(3) > model_h(1)) then
      ! above top of model
      frac = 1.0_r8
      write(string1, *) 'ob height ',old_array(3),' above CAM levels at ' &
                          ,old_array(1) ,old_array(2) ,' for kind',old_kind
      call error_handler(E_MSG, 'convert_vert', string1,source,revision,revdate)
   elseif (bot_lev <= lev%length) then
      ! within model levels
      frac = (old_array(3) - model_h(bot_lev)) / (model_h(top_lev) - model_h(bot_lev))
   else
      ! below bottom of model
      frac = 0.0_r8
      write(string1, *) 'ob height ',old_array(3),' below CAM levels at ' &
                          ,old_array(1) ,old_array(2) ,' for kind',old_kind
      call error_handler(E_MSG, 'convert_vert', string1,source,revision,revdate)
   endif

   new_pressure = (1.0_r8 - frac) * p_col(bot_lev) + frac * p_col(top_lev)

   if (vert_coord == 'pressure') then
      new_array(3) = new_pressure
      new_which = VERTISPRESSURE
   else if (vert_coord == 'log_invP') then
      new_array(3) = scale_height(p_surface=p_surf(1), p_above=new_pressure)
      new_which = VERTISSCALEHEIGHT
   endif

   deallocate(model_h)

else
   write(string1, *) 'model which_vert = ',old_which,' not handled in convert_vert '
   call error_handler(E_ERR, 'convert_vert', string1,source,revision,revdate)
endif

deallocate(p_col)

return

end subroutine convert_vert

! End of get_close section

!#######################################################################

! Initial conditions for DART

!------------------------------------------------------------------
!> Perturbs a model state copy for generating initial ensembles.
!> Routine which could provide a custom perturbation routine to
!> generate initial ensembles.  The default (if interface is not
!> provided) is to add gaussian noise to each item in the state vector.
!>
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

subroutine pert_model_copies(state_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_handle
integer,   intent(in) :: ens_size
real(r8),  intent(in) :: pert_amp
logical,  intent(out) :: interf_provided

type(model_type)        :: var_temp
integer                 :: j, k, m, pert_fld, mode, field_num
integer                 :: dim1, dim2, dim3, member
integer, save           :: seed
logical                 :: perturbed
integer(i8)             :: start_index, end_index, i ! for a variable
integer                 :: copy
real(r8)                :: random_number

! The input is a single model state vector that has (different) gaussian
! noise added to each member to generate an initial ensemble.

if (.not. module_initialized) call static_init_model()

interf_provided = .true.

! This will make the results reproduce for runs with the same number of MPI tasks.
! It will NOT give the same random sequence if you change the task count.
k = (my_task_id()+1) * 1000
call init_random_seq(random_seq, k)

pert_fld = 1

Vars2Perturb : do pert_fld=1,100
   if (pert_names(pert_fld) == ' ') exit Vars2Perturb

   ! Keep track of whether or not this field is matched and was perturbed.
   perturbed = .false.

   ExistingVars : do m=1,nflds

      if (pert_names(pert_fld) /= cflds(m)) cycle ExistingVars

      perturbed = .true.

      start_index = get_index_start(component_id, m)
      end_index = get_index_end(component_id, m)
      if (output_task0) then
         write(string1,'(3A,2I8,A,I8)') 'Perturbing ',trim(pert_names(pert_fld)), &
               ' start,stop = ',start_index,end_index,' seed=', k
         call error_handler(E_MSG,'pert_model_copies', string1)
      endif

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

      if (print_details .and. output_task0) then
         write(string1,'(2A,I8,A,1pE12.3)') &
              '   WARNING: filter_nml:perturbation_amplitude is not being used. ', &
              '   INSTEAD: model_nml:pert_sd(',mode,') = ',pert_sd(mode) 
         call error_handler(E_WARN,'pert_model_copies', string1)
      endif

      ! Handle the fields

      ! reset base values to value provided in namelist.
      if (pert_base_vals(mode) /= MISSING_R8) then
         if (print_details) then
            write(string1,*) 'Using a new base value ',pert_base_vals(mode), 'for ',cflds(m)
            call error_handler(E_MSG, 'pert_model_copies', string1, source, revision, revdate)
         endif
         where (state_handle%my_vars > start_index .and. state_handle%my_vars < end_index)
            state_handle%copies(copy, :) = pert_base_vals(mode)
         endwhere
      endif

      ! randomly perturb each point around its base value.
      if (pert_sd(pert_fld) > 0.0_r8 ) then
         do i = 1, state_handle%my_num_vars
            if (state_handle%my_vars(i) >= start_index .and. state_handle%my_vars(i) <= end_index) then
               do copy = 1, ens_size
! RMA-KR This looks like it gives the same seed to each member.
!        But this is how it was done in the trunk, which worked.
                  state_handle%copies(copy, i) = random_gaussian(random_seq, state_handle%copies(copy, i), pert_sd(mode))
               enddo
            endif
         enddo
      endif

   enddo ExistingVars

   if (.not. perturbed) then
      write(string1,*)trim(pert_names(pert_fld)),' not found in list of state variables.'
      write(string2,*)'but was supposed to be used to perturb.'
      call error_handler(E_ERR,'pert_model_copies', string1, source, revision, revdate, text2=string2)
   endif

enddo Vars2Perturb

end subroutine pert_model_copies

! End of initial model state section

!#######################################################################

! Utility routines; called by several main subroutines

!-----------------------------------------------------------------------

function index_from_grid(lev_ind, lon_ind, lat_ind, ifld)

! Calculate the index into the state vector, given the coordinate indices
! and the field number (out of nflds state vector fields).

integer, intent(in) :: lev_ind
integer, intent(in) :: lon_ind
integer, intent(in) :: lat_ind
integer, intent(in) :: ifld
integer(i8)         :: index_from_grid

integer :: i, j, k

i = -1
j = -1
k = -1

! Need to convert from lev_ind, lon_ind, lat_ind to i, j, k
!>  @todo Should just store staggared info for each ifld in static_init_model_mod
!RMA-KR; these sections could be condensed into 1, inside a loop over the 3 dimensions,
!        by defining ijk(3)
if (get_dim_name(component_id, ifld, 1) == 'lev')        i = lev_ind
if (get_dim_name(component_id, ifld, 1) == 'lon' .or. &
    get_dim_name(component_id, ifld, 1) == 'slon')       i = lon_ind
if (get_dim_name(component_id, ifld, 1) == 'lat' .or. &
    get_dim_name(component_id, ifld, 1) == 'slat')       i = lat_ind

if (get_dim_name(component_id, ifld, 2) == 'lev')        j = lev_ind
if (get_dim_name(component_id, ifld, 2) == 'lon' .or. &
    get_dim_name(component_id, ifld, 2) == 'slon')       j = lon_ind
if (get_dim_name(component_id, ifld, 2) == 'lat' .or. &
    get_dim_name(component_id, ifld, 2) == 'slat')       j = lat_ind


if (get_dim_name(component_id, ifld, 3) == 'lev')        k = lev_ind
if (get_dim_name(component_id, ifld, 3) == 'lon' .or. &
    get_dim_name(component_id, ifld, 3) == 'slon')       k = lon_ind
if (get_dim_name(component_id, ifld, 3) == 'lat' .or. &
    get_dim_name(component_id, ifld, 3) == 'slat')       k = lat_ind

index_from_grid = get_dart_vector_index(i, j, k, component_id, ifld)


end function index_from_grid

!-----------------------------------------------------------------------

function find_name(nam, list)

character(len=*), intent(in) :: nam
character(len=*), intent(in) :: list(:)
integer                      :: find_name

integer :: i

! find_name = 0
find_name = MISSING_I
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
! RMA-KR; the ncol section was removed for this CAM-FV model_mod.
!         Will be needed for the CAM-SE version.

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

subroutine plevs_cam (p_surf, n_levels, pmid )

! Define the pressures of the layer midpoints from the
! coordinate definitions and the surface pressure.

real(r8), intent(in)  :: p_surf    ! Surface pressure (pascals)
integer,  intent(in)  :: n_levels
real(r8), intent(out) :: pmid(lev%length)   ! Pressure at model levels

integer :: k

! Set midpoint pressures and layer thicknesses

do k=1,n_levels
   pmid(k) = hyam%vals(k)*P0%vals(1) + hybm%vals(k)*p_surf
enddo

end subroutine plevs_cam

!-----------------------------------------------------------------------

subroutine model_heights(state_handle, ens_size, n_levels, p_surf, base_obs_loc, model_h, istatus)

! This routine calculates geometrical height (m) at mid-layers of the CAM model
!
! was Hui's dcz2ccm1
!    has globally defined inputs:
!          hyam(lev%length),hybm(lev%length),hyai(lev%length),hybi(lev%length) =
!          hybrid vertical coefficients, top to bottom.
!          (P = P0*hyam + ps*hybm)
!          P0 - Hybrid base pressure (pascals)
!
! Kevin Raeder converted to single column version 4/28/2006
!              removed longitude dimension entirely and extra arrays 10/2006
!   5/31/2013; Rewritten to adapt to convert_vert handling obs TYPEs,
!              not obs QTYs, and to handle lonlat and cubed sphere
!              grids/interpolations.

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: n_levels
real(r8),            intent(in) :: p_surf(ens_size)
type(location_type), intent(in) :: base_obs_loc

real(r8), intent(out) :: model_h(n_levels, ens_size)
integer,  intent(out) :: istatus(ens_size)

! local variables; 
real(r8), dimension(ens_size, n_levels) :: phi, tv, q, t, mmr_o1, mmr_o2, mmr_h1, mmr_n2
real(r8), dimension(ens_size)           :: h_surf, ht_tmp
real(r8), dimension(n_levels+1,2)       :: hybrid_As, hybrid_Bs

! CS Should these come from common_mod?
! That might be inconsistent with how levels, etc were defined in CAM originally.
! DART's values are 287.0_r8 and 461.6_r8.
real(r8), parameter :: rd = 287.05_r8
real(r8), parameter :: rv = 461.51_r8
real(r8), parameter :: rr_factor = (rv/rd) - 1.0_r8

real(r8) :: lon_lat_lev(3)
type(location_type) :: temp_obs_loc

integer :: k, i, imem
integer :: vstatus(ens_size)
istatus(:) = 1
model_h(:,:) = MISSING_R8

! RMA-KR; CAM-SE section was removed from here.

! lat, lon and vertical in height
lon_lat_lev = get_location(base_obs_loc)

!> @todo I don't think hybrid_As and hybrid_Bs change thoughout a run of filter
! RMA-KR; That's true.  They could be put in global storage and initialized in static_init_mod
!         after hy[ab][im] have been read in.
! copy hybrid_As, hybrid_Bs to temporary arrays to pass to dcz2
! All arrays except hybrid_As, hybrid_Bs are oriented top to bottom.

! The 'interface' levels have an 'extra' level at model bottom, compared to the midpoint levels.
! Initialize this extra level, before filling the rest in a loop.
k = n_levels +1
hybrid_As(1,1) = hyai%vals(k)
hybrid_Bs(1,1) = hybi%vals(k)

!   hyam(n_levels) = 0 -> hybrid_As(2,2) = 0, so it
!   would be safe to set  hybrid_As(1,2) = 0.
!   It's safe because this element is used to set pmln in dcz2, but that element of pmln is never used.
hybrid_As(1,2) = 0.0_r8

! hyb[im]  have non-0 values at the bottom, 0s at the top;
!      hyb[im] coeffs multiply sigma in the calculation of pressure on levels,
!      and CAM's vertical coord is pure sigma at the bottom, so hybrid_Bs = 1.0 there.
hybrid_Bs(1,2) = 1.0_r8

! mid-points: 2nd dimension of hybrid_[AB]s = 2
! note that hyXm(n_levels + 1) is not defined (= MISSING_R8)
do k = 2,n_levels +1
   i = n_levels +2 - k
   hybrid_As(k,1) = hyai%vals(i)
   hybrid_Bs(k,1) = hybi%vals(i)
   hybrid_As(k,2) = hyam%vals(i)
   hybrid_Bs(k,2) = hybm%vals(i)
enddo

! Calculate h_surf and tv for this column, for use by dcz2.
call interp_lonlat(state_handle, ens_size, base_obs_loc, QTY_SURFACE_ELEVATION, h_surf, vstatus)
if (any(vstatus == 1)) then
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
do k = 1, n_levels
   ! construct a location with the same lat/lon but cycle though the model levels
   temp_obs_loc = set_location(lon_lat_lev(1), lon_lat_lev(2), real(k,r8), VERTISLEVEL)

   call interp_lonlat(state_handle, ens_size, temp_obs_loc, QTY_TEMPERATURE, t(:, k), vstatus)
   if (any(vstatus == 1)) then
      write(string1,'(A,I2,A)') 'Temperature level ',k, &
           ' could not be interpolated in interp_lonlat'
      call error_handler(E_WARN, 'model_heights', string1)
      return
   endif
   call interp_lonlat(state_handle, ens_size, temp_obs_loc, QTY_SPECIFIC_HUMIDITY, q(:, k), vstatus)
   if (any(vstatus == 1)) then
      write(string1,'(A,I2,A)') 'specific humidity level ',k, &
           ' could not be interpolated in interp_lonlat'
      call error_handler(E_WARN, 'model_heights', string1)
      return
   endif

   tv(:, k) = t(:, k)*(1.0_r8 + rr_factor*q(:, k))
enddo

do imem = 1, ens_size
   call dcz2(n_levels, p_surf(imem), h_surf(imem), tv(imem,:), P0%vals(1) , &
             hybrid_As, hybrid_Bs, phi(imem,:))
enddo

! used; hybrid_Bs, hybrid_As, hprb
! output from dcz2;  phi

! Conversion from geopotential height to geometric height depends on latitude
! Convert to kilometers for gph2gmh call, then back to meters for return value.
do k = 1,n_levels
   ht_tmp(:) = phi(:, k) * 0.001_r8        ! convert to km for following call only
   do imem = 1, ens_size
      model_h(k, imem) = gph2gmh(ht_tmp(imem), lon_lat_lev(2)) * 1000.0_r8
   enddo
enddo

! model_heights returns only istatus 0 or 1
! RMA-KR; model_heights uses a somewhat different status strategy than get_val_...
!         It uses 'return's (istatus(:)(all) set to 1) if it fails along the way,
!         rather than continuing on with calculations for those members that don't fail.
!         So if we arrived here, all is well, and return 'success' in all values of istatus.
!         
istatus = 0

end subroutine  model_heights

!-----------------------------------------------------------------------

subroutine dcz2(kmax,p_surf,h_surf,tv,hprb,hybrid_As,hybrid_Bs,z2)

! Compute geopotential height for a CESM hybrid coordinate column.
! All arrays except hybrid_As, hybrid_Bs are oriented top to bottom.
! hybrid_[AB]s first subscript:
!  = 1 for layer interfaces
!  = 2 for layer midpoints
! hybrid_As coord coeffs for P0 reference pressure term in plevs_cam
! hybrid_Bs coord coeffs for surf pressure term in plevs_cam (in same format as hybrid_As)

integer,  intent(in)  :: kmax                ! Number of vertical levels
real(r8), intent(in)  :: p_surf              ! Surface pressure           (pascals)
real(r8), intent(in)  :: h_surf               ! Surface height (m)
real(r8), intent(in)  :: tv(kmax)            ! Virtual temperature, top to bottom
real(r8), intent(in)  :: hprb                ! Hybrid base pressure       (pascals)
real(r8), intent(in)  :: hybrid_As(kmax+1,2)
real(r8), intent(in)  :: hybrid_Bs(kmax+1,2)
real(r8), intent(out) :: z2(kmax)            ! Geopotential height, top to bottom

! Local variables
real(r8), parameter :: r = 287.04_r8    ! Different than model_heights !
real(r8), parameter :: g0 = 9.80616_r8  ! Different than model_heights:gph2gmh:G !
real(r8), parameter :: rbyg=r/g0
real(r8) :: pterm(kmax)         ! pressure profile
real(r8) :: pmln(kmax+1)        ! logs of midpoint pressures

integer  :: i,k,l
real(r8) :: ARG

! Compute intermediate quantities using scratch space

! DEBUG: z2 was unassigned in previous code.
z2(:) = MISSING_R8

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

! Initialize z2 to sum of ground height and thickness of top half layer
! DEBUG; this is NOT adding the thickness of the 'top' half layer.
!        it's adding the thickness of the half layer at level K,
do K = 1,kmax - 1
   z2(k) = h_surf + rbyg*tv(k)*0.5_r8* (pmln(K+1)-pmln(K))
enddo
z2(kmax) = h_surf + rbyg*tv(kmax)* (log(p_surf*hybrid_Bs(1,1))-pmln(kmax))

! DEBUG; THIS is adding the half layer at the BOTTOM.
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
! the sizes of the dimensions listed in f_dim_RANKd.  Those are stored in f_dim_max.

allocate(var%vars_0d(                                              state_num_0d))
allocate(var%vars_1d(f_dim_max(1,1),                               state_num_1d))
allocate(var%vars_2d(f_dim_max(1,2),f_dim_max(2,2),                state_num_2d))
allocate(var%vars_3d(f_dim_max(1,3),f_dim_max(2,3),f_dim_max(3,3), state_num_3d))

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
!> Subroutine end_model
!> deallocates arrays that are in module global storage.

subroutine end_model()

deallocate(dim_names, dim_sizes, phis)
deallocate(state_long_names, state_units)
deallocate(cflds)

if (allocated(f_dim_3d)) then
   deallocate(f_dim_3d, f_dimid_3d)
endif
if (allocated(f_dim_2d)) then
   deallocate(f_dim_2d, f_dimid_2d)
endif
if (allocated(f_dim_1d)) then
   deallocate(f_dim_1d, f_dimid_1d)
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

end subroutine end_model

!-------------------------------------------------------------------------
!> This replaces set_ps_arrays.  It handles the whole ensemble,
!> when needed, as required by RMA.
function get_surface_pressure(state_handle, ens_size, lon_ind, lat_ind)

integer,             intent(in)  :: ens_size
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: lon_ind
integer,             intent(in)  :: lat_ind

real(r8)    :: get_surface_pressure(ens_size)
integer     :: ifld ! pressure field index
integer(i8) :: ind  ! index into state vector

ifld = find_name('PS      ',cflds)

! find index into state
ind = index_from_grid(1, lon_ind, lat_ind, ifld)

! get correct piece of state
get_surface_pressure = get_state(ind, state_handle)

end function get_surface_pressure

!-----------------------------------------------------------------------
subroutine update_vstatus(ens_size, current_vstatus, vstatus)

integer, intent(in)  :: ens_size
integer, intent(in)  :: current_vstatus(ens_size)
integer, intent(out) :: vstatus(ens_size)
logical :: bail_out ! quit because all the ensemble members have failed

! RMA-KR; Is this bail_out code commented out because it's handled in the calling routines?
!bail_out = .false.
! only update if there are no previous failures
where(vstatus == 0) vstatus = current_vstatus
!if(all(vstatus /= 0)) bail_out = .true.  ! Every ensemble member has reached failure

end subroutine update_vstatus
!-----------------------------------------------------------------------

! RMA-KR; set_print_details is not used in this module.
subroutine set_print_details(how)

! reset the print_details module global variable to control
! how much output there is

logical, intent(in) :: how

print_details = how

end subroutine set_print_details

!--------------------------------------------------------------------
!> construct restart file name for reading
!> model time for CESM format?
function construct_file_name_in(stub, domain, copy)

character(len=512), intent(in) :: stub
integer,            intent(in) :: domain
integer,            intent(in) :: copy
character(len=256) :: construct_file_name_in

! fv_testcase.cam_0003.i.2004-01-15-00000.nc
! RMA-KR; Why is the file type (i) and date hard-wired?
!         Where is this used?
!            io/io_filenames_mod.f90;  when restart name can't be read from rpointer, build a name.
write(construct_file_name_in, '(A, i4.4, A)') TRIM(stub), copy, '.i.2004-01-15-00000.nc'

end function construct_file_name_in

!--------------------------------------------------------------------
!> pass the vertical localization coordinate to assim_tools_mod
function query_vert_localization_coord()

integer :: query_vert_localization_coord

query_vert_localization_coord = VERTISUNDEF

if (vert_coord == 'pressure') query_vert_localization_coord = VERTISPRESSURE
if (vert_coord == 'log_invP') query_vert_localization_coord = VERTISSCALEHEIGHT

end function query_vert_localization_coord

!--------------------------------------------------------------------
!> read the time from the input file
function read_model_time(file_name)

character(len=256), intent(in) :: file_name

type(time_type) :: read_model_time

integer :: i, k, n, m, ifld  
integer :: nc_file_ID, nc_var_ID, dimid, varid, dimlen
integer :: iyear, imonth, iday, ihour, imin, isec, rem
integer :: timestep
integer,  allocatable :: datetmp(:), datesec(:)

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

read_model_time = set_date(iyear,imonth,iday,ihour,imin,isec)

end function read_model_time

!-----------------------------------------------------------------------
!>@todo this routine should write the model time when 
!>      creating files from scratch
subroutine write_model_time(ncid, dart_time)

integer,             intent(in) :: ncid !< netcdf file handle
type(time_type),     intent(in) :: dart_time

call error_handler(E_MSG, 'write_model_time', 'no routine for cam-fv write model time')

end subroutine write_model_time


!--------------------------------------------------------------------
!> Construct an arry to pass to add_domain that contains the clamping info
!> for each variable. Note for non-netcdf read this is done in write_cam_init
subroutine set_clamp_fields(clampfield)

real(r8), intent(out) :: clampfield(nflds, 2) ! min, max for each field

integer :: i

clampfield(:, :) = MISSING_R8 ! initalize to no clamping

do i = 1, nflds
   if(cflds(i) == 'Q')      clampfield(i, 1) = 1.0e-12_r8
   if(cflds(i) == 'CLDLIQ') clampfield(i, 1) =  0.0_r8
   if(cflds(i) == 'CLDICE') clampfield(i, 1) =  0.0_r8
enddo

end subroutine

!--------------------------------------------------------------------
function get_lon_name(var)

integer, intent(in) :: var ! s_type - order in state vectors
character(len=8) :: get_lon_name

integer :: i

get_lon_name = 'lon'  ! default to not staggered

do i = 1, get_num_dims(component_id, var)
   if (get_dim_name(component_id, var, i)=='slon') then
      get_lon_name = 'slon'
      exit
   endif
enddo

end function get_lon_name

!--------------------------------------------------------------------
function get_lat_name(var)

integer, intent(in) :: var ! s_type - order in state vectors
character(len=8) :: get_lat_name

integer :: i

get_lat_name = 'lat'  ! default to not staggered

do i = 1, get_num_dims(component_id, var)
   if (get_dim_name(component_id, var, i)=='slat') then
      get_lat_name = 'slat'
      exit
   endif
enddo

end function get_lat_name

!--------------------------------------------------------------------
function get_lev_name(var)

integer, intent(in) :: var ! s_type - order in state vectors
character(len=8) :: get_lev_name

integer :: i

get_lev_name = 'lev'  ! default to not staggered

do i = 1, get_num_dims(component_id, var)
   if (get_dim_name(component_id, var, i)=='ilev') then
      get_lev_name = 'ilev'
      exit
   endif
enddo

end function get_lev_name

!#######################################################################
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> @mainpage Remote Memory Access version of DART code.
!>
!> Forward operator, vertical conversion
!>
!> WRF
!>
!> For distributed phb array see distributed_phb_model_mod.f90
!> \todo To do list
!> @author dart@ucar.edu
!>
!> \subpage test
!> 
!> \page test

!> WRF model mod

module model_mod

! Assimilation interface for WRF model
!> \defgroup wrf model_mod
!> Model mod
!> 
!> Model mod for WRF
!> @{
!-----------------------------------------------------------------------
!
!     interface for WRF
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------
!-----------------------------------------------------------------------

use         types_mod,   only : r8, i8, deg2rad, missing_r8, ps0, earth_radius, &
                                gas_constant, gas_constant_v, gravity, pi,      &
                                digits12

use  time_manager_mod,   only : time_type, set_time, set_calendar_type, GREGORIAN,&
                                set_date, get_date

use      location_mod,   only : location_type, get_location, set_location, &
                                query_location, VERTISUNDEF, VERTISSURFACE, &
                                VERTISLEVEL, VERTISPRESSURE, VERTISHEIGHT, &
                                VERTISSCALEHEIGHT, vertical_localization_on, &
                                set_vertical_localization_coord, &
                                get_close_type, get_dist, is_vertical, &
                                loc_get_close => get_close_obs

use     utilities_mod,  only  : file_exist, open_file, close_file, &
                                register_module, error_handler, E_ERR, E_WARN, &
                                E_MSG, nmlfileunit, logfileunit, do_output, &
                                find_namelist_in_file, check_namelist_read, &
                                find_textfile_dims, file_to_text, &
                                do_nml_file, do_nml_term, scalar

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, nc_check, &
                                 nc_begin_define_mode, nc_end_define_mode

use  mpi_utilities_mod,  only : my_task_id, task_count, all_reduce_min_max

use     random_seq_mod,  only : random_seq_type, init_random_seq, random_gaussian

use      obs_kind_mod,   only : QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, &
                                QTY_SURFACE_PRESSURE, QTY_TEMPERATURE, &
                                QTY_SPECIFIC_HUMIDITY, QTY_SURFACE_ELEVATION, &
                                QTY_PRESSURE, QTY_VERTICAL_VELOCITY, &
                                QTY_DENSITY, QTY_FLASH_RATE_2D, &
                                QTY_RAINWATER_MIXING_RATIO, QTY_HAIL_MIXING_RATIO, &
                                QTY_GRAUPEL_MIXING_RATIO, QTY_SNOW_MIXING_RATIO, &
                                QTY_CLOUDWATER_MIXING_RATIO, QTY_ICE_MIXING_RATIO, &
                                QTY_CLOUD_FRACTION, QTY_CONDENSATIONAL_HEATING,  &
                                QTY_VAPOR_MIXING_RATIO, QTY_2M_TEMPERATURE, &
                                QTY_2M_SPECIFIC_HUMIDITY, QTY_10M_U_WIND_COMPONENT, &
                                QTY_10M_V_WIND_COMPONENT, &
                                QTY_ICE_NUMBER_CONCENTRATION, QTY_GEOPOTENTIAL_HEIGHT, &
                                QTY_POTENTIAL_TEMPERATURE, QTY_SOIL_MOISTURE, &
                                QTY_DROPLET_NUMBER_CONCENTR, QTY_SNOW_NUMBER_CONCENTR, &
                                QTY_RAIN_NUMBER_CONCENTR, QTY_GRAUPEL_NUMBER_CONCENTR, &
                                QTY_HAIL_NUMBER_CONCENTR, QTY_HAIL_VOLUME, &
                                QTY_GRAUPEL_VOLUME, QTY_DIFFERENTIAL_REFLECTIVITY, &
                                QTY_RADAR_REFLECTIVITY, QTY_POWER_WEIGHTED_FALL_SPEED, &
                                QTY_SPECIFIC_DIFFERENTIAL_PHASE, &
                                QTY_VORTEX_LAT, QTY_VORTEX_LON, &
                                QTY_VORTEX_PMIN, QTY_VORTEX_WMAX, &
                                QTY_SKIN_TEMPERATURE, QTY_LANDMASK, &
                                QTY_SURFACE_TYPE, &
                                get_index_for_quantity, get_num_quantities, &
                                get_name_for_quantity

use ensemble_manager_mod,  only : ensemble_type, get_my_num_vars, get_my_vars

use sort_mod,              only : sort

use distributed_state_mod, only : get_state

use default_model_mod,   only : adv_1step, nc_write_model_vars, &
                                init_conditions => fail_init_conditions, &
                                init_time => fail_init_time

use state_structure_mod, only : add_domain, get_model_variable_indices, &
                                state_structure_info, &
                                get_index_start, get_index_end, &
                                get_dart_vector_index

! FIXME:
! the kinds 
! QTY_ICE_NUMBER_CONCENTRATION should be QTY_ICE_NUMBER_CONCENTR
! to be consistent with the other concentration names.

!nc -- module_map_utils split the declarations of PROJ_* into a separate module called
!nc --   misc_definitions_module 
use         map_utils, only : proj_info, map_init, map_set, latlon_to_ij, &
                              ij_to_latlon, gridwind_to_truewind

use misc_definitions_module, only : PROJ_LATLON, PROJ_MERC, PROJ_LC, PROJ_PS, PROJ_CASSINI, &
                                    PROJ_CYL

use netcdf
use typesizes

implicit none
private


!-----
! DART requires 18 specific public interfaces from model_mod.f90 

! routines with code in this module
public ::  get_model_size,                &
           get_state_meta_data,           &
           shortest_time_between_assimilations, &
           static_init_model,             &
           model_interpolate,             &
           nc_write_model_atts,           &
           get_close_obs,                 &
           get_close_state,               &
           convert_vertical_obs,          &
           convert_vertical_state,        &
           pert_model_copies,             &
           read_model_time,               &
           write_model_time,              &
           end_model


! required routines where the code is in another module
public ::  adv_1step,          &
           init_time,          &
           init_conditions,    &
           nc_write_model_vars

!-----
! Here is the appropriate place for other users to make additional routines
!   contained within model_mod available for public use:
public ::  get_number_domains,          &
           get_wrf_static_data,         &
           model_pressure_distrib,      &
           model_height_distrib,        &
           pres_to_zk,                  &
           height_to_zk,                &
           get_domain_info,             &
           get_wrf_state_variables,     &
           fill_default_state_table,    &
           read_wrf_dimensions,         &
           get_number_of_wrf_variables, &
           get_variable_bounds,         &
           set_variable_bound_defaults, &
           get_variable_size_from_file, &
           get_wrf_date, set_wrf_date,  &
           vert_convert,                &
           height_diff_check

! public parameters
public :: max_state_variables,     &
          num_state_table_columns, &
          num_bounds_table_columns

! types
public :: wrf_dom, wrf_static_data_for_dart

! Interfaces for array and single value versions of subroutines/functions
! This is because the forward operator works on the whole ensemble, and the
! vertical conversion only uses the mean copy.


! HK ? interp_4pressure interface needed?

!-----------------------------------------------------------------------
! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'wrf/model_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

! miscellaneous
integer, parameter :: max_state_variables = 100
integer, parameter :: num_state_table_columns = 5
integer, parameter :: num_bounds_table_columns = 4

!-----------------------------------------------------------------------
! Model namelist parameters with default values.
!
! center_search_half_length:  half length (in meter) of the searching box to locate 
!                             minimum pressure at a grid point
! center_spline_scale: coarse grid to spline interp. fine grid ratio
!-----------------------------------------------------------------------

logical :: output_state_vector     = .false.  ! output prognostic variables
logical :: default_state_variables = .true.   ! use default state list?
character(len=129) :: wrf_state_variables(num_state_table_columns,max_state_variables) = 'NULL'
character(len=129) :: wrf_state_bounds(num_bounds_table_columns,max_state_variables) = 'NULL'
integer :: num_domains          = 1
integer :: calendar_type        = GREGORIAN
integer :: assimilation_period_seconds = 21600
! Max height a surface obs can be away from the actual model surface
! and still be accepted (in meters)
real (kind=r8) :: sfc_elev_max_diff  = -1.0_r8   ! could be something like 200.0_r8
real (kind=r8) :: center_search_half_length = 500000.0_r8
real(r8) :: circulation_pres_level = 80000.0_r8
real(r8) :: circulation_radius     = 108000.0_r8
integer :: center_spline_grid_scale = 10
integer :: vert_localization_coord = VERTISHEIGHT
! Allow observations above the surface but below the lowest sigma level.
logical :: allow_obs_below_vol = .false.
! Do the interpolation of pressure values only after taking the log (.true.)
! vs doing a linear interpolation directly in pressure units (.false.)
logical :: log_vert_interp  = .true.
logical :: log_horz_interpM = .false.
logical :: log_horz_interpQ = .false.
!nc -- we are adding these to the model.nml until they appear in the NetCDF files
logical :: polar = .false.         ! wrap over the poles
logical :: periodic_x = .false.    ! wrap in longitude or x
logical :: periodic_y = .false.    ! used for single column model, wrap in y
!JPH -- single column model flag 
logical :: scm        = .false.    ! using the single column model
logical :: allow_perturbed_ics = .false.  ! should spin the model up for a while after

! obsolete items; ignored by this code. 
! non-backwards-compatible change. should be removed, 
! but see note below about namelist.
integer :: num_moist_vars
logical :: surf_obs, soil_data, h_diab

! adv_mod_command moved to dart_to_wrf namelist; ignored here.
character(len = 72) :: adv_mod_command = ''

! num_moist_vars, surf_obs, soil_data, h_diab, and adv_mod_command
! are IGNORED no matter what their settings in the namelist are.
! they are obsolete, but removing them here will cause a fatal error
! until users remove them from their input.nml files as well.
namelist /model_nml/ num_moist_vars, &
                     num_domains, calendar_type, surf_obs, soil_data, h_diab, &
                     default_state_variables, wrf_state_variables, &
                     wrf_state_bounds, sfc_elev_max_diff, &
                     adv_mod_command, assimilation_period_seconds, &
                     allow_obs_below_vol, vert_localization_coord, &
                     center_search_half_length, center_spline_grid_scale, &
                     circulation_pres_level, circulation_radius, polar, &
                     periodic_x, periodic_y, scm, allow_perturbed_ics

! if you need to check backwards compatibility, set this to .true.
! otherwise, leave it as false to use the more correct geometric height
logical :: use_geopotential_height = .false.

character(len = 20) :: wrf_nml_file = 'namelist.input'
logical :: have_wrf_nml_file = .false.
integer :: num_obs_kinds = 0
logical, allocatable :: in_state_vector(:)
integer, allocatable  :: domain_id(:) ! Global storage to interface with state_structure_mod.

!-----------------------------------------------------------------------

! Private definition of domain map projection use by WRF

!nc -- added in CASSINI and CYL according to module_map_utils convention
!JPH -- change variable name from map_sphere to more specific map_latlon
integer, parameter :: map_latlon = 0, map_lambert = 1, map_polar_stereo = 2, map_mercator = 3
integer, parameter :: map_cassini = 6, map_cyl = 5

! Private definition of model variable types

real (kind=r8), PARAMETER    :: kappa = 2.0_r8/7.0_r8 ! gas_constant / cp
real (kind=r8), PARAMETER    :: ts0 = 300.0_r8        ! Base potential temperature for all levels.

! hybrid_opt 2 is a hybrid coordinate - terrain following at the
! surface, and straight pressure levels at the top.
integer, parameter :: VERT_HYBRID = 2

!---- private data ----

! Got rid of surf_var as global private variable for model_mod and just defined it locally
!   within model_interpolate

TYPE wrf_static_data_for_dart

   integer  :: bt, bts, sn, sns, we, wes, sls
   real(r8) :: dx, dy, dt, p_top
   integer  :: map_proj
   real(r8) :: cen_lat,cen_lon
   type(proj_info) :: proj

   ! Boundary conditions -- hopefully one day these will be in the global attributes of the
   !   input NetCDF file ("periodic_x" and "polar" are namelist items in the &bdy_control
   !   section of a standard WRF "namelist.input" file), but for now we have included them
   !   in the "model_nml" group of DART's own "input.nml".  Above, their default values are
   !   both set to .false. 
   logical  :: periodic_x
   logical  :: periodic_y
   logical  :: polar
   logical  :: scm

   integer  :: domain_size
   integer  :: localization_coord
   integer  :: hybrid_opt
   real(r8), dimension(:),     pointer :: znu, dn, dnw, zs, znw, c1h, c2h
   real(r8), dimension(:,:),   pointer :: mub, hgt
   real(r8), dimension(:,:),   pointer :: latitude, latitude_u, latitude_v
   real(r8), dimension(:,:),   pointer :: longitude, longitude_u, longitude_v
   real(r8), dimension(:,:,:), pointer :: phb

   ! NEWVAR:  Currently you have to add a new type here if you want to use
   ! NEWVAR:  a WRF variable which is not one of these types.  This will go
   ! NEWVAR:  away eventually, we hope.  Search for NEWVAR for other places
   ! NEWVAR:  the code has to change.

   ! JPH local variables to hold type indices
   integer :: type_u, type_v, type_w, type_t, type_qv, type_qr, type_hdiab, &
              type_qndrp, type_qnsnow, type_qnrain, type_qngraupel, type_qnice, &
              type_qc, type_qg, type_qi, type_qs, type_gz, type_refl, type_fall_spd, &
              type_dref, type_spdp, type_qh, type_qnhail, type_qhvol, type_qgvol, &
              type_cldfra

   integer :: type_u10, type_v10, type_t2, type_th2, type_q2, &
              type_ps, type_mu, type_tsk, type_tslb, type_sh2o, &
              type_smois, type_2dflash

   integer :: number_of_wrf_variables
   integer(i8), dimension(:,:), pointer :: var_index
   integer,     dimension(:,:), pointer :: var_size
   integer,     dimension(:),   pointer :: var_type
   integer,     dimension(:),   pointer :: var_index_list
   logical,     dimension(:),   pointer :: var_update_list
   integer,     dimension(:),   pointer :: dart_kind
   integer,     dimension(:,:), pointer :: land
   real(r8),    dimension(:),   pointer :: lower_bound,upper_bound
   character(len=10), dimension(:),pointer :: clamp_or_fail
   character(len=129),dimension(:),pointer :: description, units, stagger, coordinates


end type wrf_static_data_for_dart

type wrf_dom
   type(wrf_static_data_for_dart), pointer :: dom(:)
   integer(i8) :: model_size
end type wrf_dom

type(wrf_dom) :: wrf

! JPH move map stuff into common (can move back into S/R later?)
real(r8) :: stdlon,truelat1,truelat2 !,latinc,loninc

! have a single, module global error string (rather than 
! replicate it in each subroutine and use up more stack space)
character(len=129) :: errstring, msgstring2, msgstring3

contains

!#######################################################################

subroutine static_init_model()

! Initializes class data for WRF

integer :: ncid
integer :: io, iunit

character (len=1)     :: idom
logical, parameter    :: debug = .false.
integer               :: ind, i, j, k, id
integer(i8)           :: dart_index
integer               :: my_index
integer               :: var_element_list(max_state_variables)
logical               :: var_update_list(max_state_variables)
real(r8)              :: var_bounds_table(max_state_variables,2)
! holds the variable names for a domain when calling add_domain
character(len=129)    :: netcdf_variable_names(max_state_variables)

!----------------------------------------------------------------------

! Register the module
call register_module(source, revision, revdate)

! Begin by reading the namelist input
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Temporary warning until this namelist item is removed.
if (adv_mod_command /= '') then
   msgstring2 = "Set the model advance command in the &dart_to_wrf_nml namelist"
   call error_handler(E_MSG, 'static_init_model:', &
         "WARNING: adv_mod_command ignored in &model_mod namelist", &
          text2=msgstring2)
endif

allocate(wrf%dom(num_domains))
allocate(domain_id(num_domains))

! get default state variable table if asked
if ( default_state_variables ) then
  wrf_state_variables = 'NULL'
  call fill_default_state_table(wrf_state_variables)
  msgstring2 = 'Set "default_state_variables" to .false. in the namelist'
  msgstring3 = 'to use the "wrf_state_variables" list instead.'
  call error_handler(E_MSG, 'static_init_model:', &
                  'Using predefined wrf variable list for dart state vector.', &
                   text2=msgstring2, text3=msgstring3)
endif

if ( debug ) then
  print*,'WRF state vector table'
  print*,'default_state_variables = ',default_state_variables
  print*,wrf_state_variables
endif

!---------------------------
! set this array so we know exactly which obs kinds are
! allowed to be interpolated (and can give a reasonably
! helpful error message if not).
!---------------------------

num_obs_kinds = get_num_quantities()
allocate(in_state_vector(num_obs_kinds))
call fill_dart_kinds_table(wrf_state_variables, in_state_vector)


! set calendar type
call set_calendar_type(calendar_type)

! Store vertical localization coordinate
! Only 4 are allowed: level(1), pressure(2), height(3), or scale height(4)
! Everything else is assumed height
if (vert_localization_coord == VERTISLEVEL) then
   wrf%dom(:)%localization_coord = VERTISLEVEL
elseif (vert_localization_coord == VERTISPRESSURE) then
   wrf%dom(:)%localization_coord = VERTISPRESSURE
elseif (vert_localization_coord == VERTISHEIGHT) then
   wrf%dom(:)%localization_coord = VERTISHEIGHT
elseif (vert_localization_coord == VERTISSCALEHEIGHT) then
   wrf%dom(:)%localization_coord = VERTISSCALEHEIGHT
else
   write(msgstring2,*)'vert_localization_coord must be one of ', &
                     VERTISLEVEL, VERTISPRESSURE, VERTISHEIGHT, VERTISSCALEHEIGHT
   write(errstring,*)'vert_localization_coord is ', vert_localization_coord
   call error_handler(E_ERR,'static_init_model', errstring, source, revision,revdate, &
                      text2=msgstring2)
endif

! the agreement amongst the dart/wrf users was that there was no need to
! read the wrf namelist, since the only thing it was extracting was the
! timestep, which is part of the wrf input netcdf file.
! call read_dt_from_wrf_nml()

dart_index = 1

WRFDomains : do id=1,num_domains

   write( idom , '(I1)') id

   ! only print this once, no matter how many parallel tasks are running
   if (do_output()) then
      write(     *     ,*) '******************'
      write(     *     ,*) '**  DOMAIN # ',idom,'  **'
      write(     *     ,*) '******************'
      write(logfileunit,*) '******************'
      write(logfileunit,*) '**  DOMAIN # ',idom,'  **'
      write(logfileunit,*) '******************'
   endif

   if(file_exist('wrfinput_d0'//idom)) then

      call nc_check( nf90_open('wrfinput_d0'//idom, NF90_NOWRITE, ncid), &
                     'static_init_model','open wrfinput_d0'//idom )

   else

      call error_handler(E_ERR,'static_init_model', &
           'Please put wrfinput_d0'//idom//' in the work directory.', source, revision,revdate)

   endif

   if(debug) write(*,*) ' ncid is ',ncid

!-------------------------------------------------------
! read WRF dimensions
!-------------------------------------------------------
   call read_wrf_dimensions(ncid,wrf%dom(id)%bt, wrf%dom(id)%bts, &
                                 wrf%dom(id)%sn, wrf%dom(id)%sns, &
                                 wrf%dom(id)%we, wrf%dom(id)%wes, &
                                 wrf%dom(id)%sls)

!-------------------------------------------------------
! read WRF file attributes
!-------------------------------------------------------
   call read_wrf_file_attributes(ncid,id)

!-------------------------------------------------------
! assign boundary condition flags
!-------------------------------------------------------

   call assign_boundary_conditions(id)

!-------------------------------------------------------
! read static data
!-------------------------------------------------------

   call read_wrf_static_data(ncid,id)


!-------------------------------------------------------
! next block set up map
!-------------------------------------------------------

   call setup_map_projection(id)

!-------------------------------------------------------
! end block set up map
!-------------------------------------------------------

! get the number of wrf variables wanted in this domain's state
   wrf%dom(id)%number_of_wrf_variables = get_number_of_wrf_variables(id,wrf_state_variables,var_element_list, var_update_list)

! allocate and store the table locations of the variables valid on this domain
   allocate(wrf%dom(id)%var_index_list(wrf%dom(id)%number_of_wrf_variables))
   wrf%dom(id)%var_index_list = var_element_list(1:wrf%dom(id)%number_of_wrf_variables)

! allocation for wrf variable types 
   allocate(wrf%dom(id)%var_type(wrf%dom(id)%number_of_wrf_variables))

! allocation for update/nocopyback/noupdate
   allocate(wrf%dom(id)%var_update_list(wrf%dom(id)%number_of_wrf_variables))
   wrf%dom(id)%var_update_list = var_update_list(1:wrf%dom(id)%number_of_wrf_variables)

! allocation for dart kinds
   allocate(wrf%dom(id)%dart_kind(wrf%dom(id)%number_of_wrf_variables))

! allocation of var size 
   allocate(wrf%dom(id)%var_size(3,wrf%dom(id)%number_of_wrf_variables))

! allocation for wrf variable metadata
   allocate(wrf%dom(id)%stagger(wrf%dom(id)%number_of_wrf_variables))
   allocate(wrf%dom(id)%description(wrf%dom(id)%number_of_wrf_variables))
   allocate(wrf%dom(id)%units(wrf%dom(id)%number_of_wrf_variables))
   allocate(wrf%dom(id)%coordinates(wrf%dom(id)%number_of_wrf_variables))

! set default bounds checking
   allocate(wrf%dom(id)%lower_bound(wrf%dom(id)%number_of_wrf_variables))
   allocate(wrf%dom(id)%upper_bound(wrf%dom(id)%number_of_wrf_variables))
   allocate(wrf%dom(id)%clamp_or_fail(wrf%dom(id)%number_of_wrf_variables))
   call set_variable_bound_defaults(wrf%dom(id)%number_of_wrf_variables, &
                                    wrf%dom(id)%lower_bound, &
                                    wrf%dom(id)%upper_bound, &
                                    wrf%dom(id)%clamp_or_fail)

!  build the variable indices
!  this accounts for the fact that some variables might not be on all domains

   do ind = 1,wrf%dom(id)%number_of_wrf_variables

      ! actual location in state variable table
      my_index =  wrf%dom(id)%var_index_list(ind)

      wrf%dom(id)%var_type(ind) = ind ! types are just the order for this domain
      wrf%dom(id)%dart_kind(ind) = get_index_for_quantity(trim(wrf_state_variables(2,my_index)))

      if ( debug ) then
         print*,'dart kind identified: ',trim(wrf_state_variables(2,my_index)),' ',wrf%dom(id)%dart_kind(ind)
      endif

      ! get stagger and variable size
      call get_variable_size_from_file(ncid,id,  &
                                       wrf_state_variables(1,my_index), &
                                       wrf%dom(id)%bt, wrf%dom(id)%bts, &
                                       wrf%dom(id)%sn, wrf%dom(id)%sns, &
                                       wrf%dom(id)%we, wrf%dom(id)%wes, & 
                                       wrf%dom(id)%stagger(ind),        &
                                       wrf%dom(id)%var_size(:,ind))

      ! get other variable metadata; units, coordinates and description
      call get_variable_metadata_from_file(ncid,id,  &
                                       wrf_state_variables(1,my_index), &
                                       wrf%dom(id)%description(ind),         &
                                       wrf%dom(id)%coordinates(ind),         &
                                       wrf%dom(id)%units(ind) )

      if ( debug ) then
         print*,'variable size ',trim(wrf_state_variables(1,my_index)),' ',wrf%dom(id)%var_size(:,ind)
      endif

      !  add bounds checking information
      call get_variable_bounds(wrf_state_bounds, wrf_state_variables(1,my_index), &
                               wrf%dom(id)%lower_bound(ind), wrf%dom(id)%upper_bound(ind), &
                               wrf%dom(id)%clamp_or_fail(ind))

      if ( debug ) then
         write(*,*) 'Bounds for variable ',  &
         trim(wrf_state_variables(1,my_index)), &
         ' are ',wrf%dom(id)%lower_bound(ind), &
         wrf%dom(id)%upper_bound(ind), &
         wrf%dom(id)%clamp_or_fail(ind)
      endif

      write(errstring, '(A,I4,2A)') 'state vector array ', ind, ' is ', trim(wrf_state_variables(1,my_index))
      call error_handler(E_MSG, 'static_init_model: ', errstring)
   enddo

   if (do_output()) then
      write(     *     ,*)
      write(logfileunit,*)
   endif

! close data file, we have all we need

   call nc_check(nf90_close(ncid),'static_init_model','close wrfinput_d0'//idom)

! indices into 1D array - hopefully this becomes obsolete
! JPH changed last dimension here from num_model_var_types
   !HK allocate(wrf%dom(id)%dart_ind(wrf%dom(id)%wes,wrf%dom(id)%sns,wrf%dom(id)%bts,wrf%dom(id)%number_of_wrf_variables))

! start and stop of each variable in vector
   allocate(wrf%dom(id)%var_index(2,wrf%dom(id)%number_of_wrf_variables))

!---------------------------
! end block to be obsolete
!---------------------------

!---------------------------
! at this point we need all information assigned to each variable
! then just loop thru the table
!---------------------------

   !HK wrf%dom(id)%dart_ind = 0

! NOTE: this could be put into the loop above if wrf%dom(id)%dart_ind is
! eventually not needed
! Here we use ind instead of type as the 4th dimension.  In is linked to the
! specific type via wrf%dom(id)%var_index_list(ind).  This saves some
! space from the previous implementation but I am not yet sure of other
! problems that it might cause.

   do ind = 1,wrf%dom(id)%number_of_wrf_variables

      my_index =  wrf%dom(id)%var_index_list(ind)

      if ( debug ) then
         write(*,*) 'Assigning dart vector indices for var_type ',wrf%dom(id)%var_type(ind)
         write(*,*) 'affiliated with WRF variable ',trim(wrf_state_variables(1,my_index)),' of size ',wrf%dom(id)%var_size(:,ind)
      endif

      wrf%dom(id)%var_index(1,ind) = dart_index
      do k=1,wrf%dom(id)%var_size(3,ind)
         do j=1,wrf%dom(id)%var_size(2,ind)
            do i=1,wrf%dom(id)%var_size(1,ind)
               ! HK wrf%dom(id)%dart_ind(i,j,k,ind) &
               !                               = dart_index
               dart_index = dart_index + 1
            enddo
         enddo
      enddo
      wrf%dom(id)%var_index(2,ind) = dart_index - 1

      if ( debug ) write(*,*) 'assigned start, stop ',wrf%dom(id)%var_index(:,ind)

   enddo ! loop through all viable state variables on this domain
   if ( id == 1 ) then
     wrf%dom(id)%domain_size = dart_index - 1
   else
     wrf%dom(id)%domain_size = dart_index - 1
     do ind = id-1, 1, -1
       wrf%dom(id)%domain_size = wrf%dom(id)%domain_size - wrf%dom(ind)%domain_size
     enddo
   endif


   ! NEWVAR: If you add a new wrf array type which is not yet in this list, currently
   ! NEWVAR: you will have to add it here, and add a type_xx for it, and also add
   ! NEWVAR: a in_state_vector case in the select statement.  search for NEWVAR.

   ! JPH now that we have the domain ID just go ahead and get type indices once
   ! NOTE: this is not strictly necessary - can use only stagger info in the future (???)
   wrf%dom(id)%type_u      = get_type_ind_from_type_string(id,'U')
   wrf%dom(id)%type_v      = get_type_ind_from_type_string(id,'V')
   wrf%dom(id)%type_w      = get_type_ind_from_type_string(id,'W')
   wrf%dom(id)%type_t      = get_type_ind_from_type_string(id,'THM')
   wrf%dom(id)%type_gz     = get_type_ind_from_type_string(id,'PH')
   wrf%dom(id)%type_qv     = get_type_ind_from_type_string(id,'QVAPOR')
   wrf%dom(id)%type_qr     = get_type_ind_from_type_string(id,'QRAIN')
   wrf%dom(id)%type_qc     = get_type_ind_from_type_string(id,'QCLOUD')
   wrf%dom(id)%type_qg     = get_type_ind_from_type_string(id,'QGRAUP')
   wrf%dom(id)%type_qh     = get_type_ind_from_type_string(id,'QHAIL')
   wrf%dom(id)%type_qi     = get_type_ind_from_type_string(id,'QICE')
   wrf%dom(id)%type_qs     = get_type_ind_from_type_string(id,'QSNOW')
   wrf%dom(id)%type_cldfra = get_type_ind_from_type_string(id,'CLDFRA')
   wrf%dom(id)%type_qgvol  = get_type_ind_from_type_string(id,'QVGRAUPEL')
   wrf%dom(id)%type_qhvol  = get_type_ind_from_type_string(id,'QVHAIL')
   wrf%dom(id)%type_qnice  = get_type_ind_from_type_string(id,'QNICE')
   wrf%dom(id)%type_qndrp  = get_type_ind_from_type_string(id,'QNDRP')
   wrf%dom(id)%type_qnsnow = get_type_ind_from_type_string(id,'QNSNOW')
   wrf%dom(id)%type_qnrain = get_type_ind_from_type_string(id,'QNRAIN')
   wrf%dom(id)%type_qngraupel = get_type_ind_from_type_string(id,'QNGRAUPEL')
   wrf%dom(id)%type_qnhail = get_type_ind_from_type_string(id,'QNHAIL')
   wrf%dom(id)%type_u10    = get_type_ind_from_type_string(id,'U10')
   wrf%dom(id)%type_v10    = get_type_ind_from_type_string(id,'V10')
   wrf%dom(id)%type_t2     = get_type_ind_from_type_string(id,'T2')
   wrf%dom(id)%type_th2    = get_type_ind_from_type_string(id,'TH2')
   wrf%dom(id)%type_q2     = get_type_ind_from_type_string(id,'Q2')
   wrf%dom(id)%type_ps     = get_type_ind_from_type_string(id,'PSFC')
   wrf%dom(id)%type_mu     = get_type_ind_from_type_string(id,'MU')
   wrf%dom(id)%type_tsk    = get_type_ind_from_type_string(id,'TSK')
   wrf%dom(id)%type_2dflash = get_type_ind_from_type_string(id,'FLASH_RATE_2D')
   wrf%dom(id)%type_tslb   = get_type_ind_from_type_string(id,'TSLB')
   wrf%dom(id)%type_smois  = get_type_ind_from_type_string(id,'SMOIS')
   wrf%dom(id)%type_sh2o   = get_type_ind_from_type_string(id,'SH2O')
   wrf%dom(id)%type_refl   = get_type_ind_from_type_string(id,'REFL_10CM')
   wrf%dom(id)%type_dref   = get_type_ind_from_type_string(id,'DIFF_REFL_10CM')
   wrf%dom(id)%type_spdp   = get_type_ind_from_type_string(id,'SPEC_DIFF_10CM')
   wrf%dom(id)%type_fall_spd = get_type_ind_from_type_string(id,'FALL_SPD_Z_WEIGHTED')
   !wrf%dom(id)%type_fall_spd = get_type_ind_from_type_string(id,'VT_DBZ_WT')
   wrf%dom(id)%type_hdiab  = get_type_ind_from_type_string(id,'H_DIABATIC')

   ! variable bound table for setting upper and lower bounds of variables 
   var_bounds_table(1:wrf%dom(id)%number_of_wrf_variables,1) = wrf%dom(id)%lower_bound
   var_bounds_table(1:wrf%dom(id)%number_of_wrf_variables,2) = wrf%dom(id)%upper_bound

   ! List of netcdf variable names in the domain
   do i = 1, wrf%dom(id)%number_of_wrf_variables
      my_index =  wrf%dom(id)%var_index_list(i) ! index in wrf_state_variables
      netcdf_variable_names(i) = wrf_state_variables(1, my_index)
   enddo


   ! add domain - not doing anything with domain_id yet so just overwriting it
   domain_id(id) = add_domain( 'wrfinput_d0'//idom, &
                           wrf%dom(id)%number_of_wrf_variables, &
                           var_names   = netcdf_variable_names(1:wrf%dom(id)%number_of_wrf_variables), &
                           kind_list   = wrf%dom(id)%dart_kind, &
                           clamp_vals  = var_bounds_table(1:wrf%dom(id)%number_of_wrf_variables,:), &
                           update_list = var_update_list(1:wrf%dom(id)%number_of_wrf_variables) )

   if (debug) call state_structure_info(domain_id(id))

enddo WRFDomains 

wrf%model_size = dart_index - 1
write(errstring,*) ' wrf model size is ',wrf%model_size
call error_handler(E_MSG, 'static_init_model', errstring)

! tell the location module how we want to localize in the vertical
call set_vertical_localization_coord(vert_localization_coord)

end subroutine static_init_model

!#######################################################################
!> Convert from the state structure id to the wrf domain number.
!> These are the same if there is only WRF involved in the assimilation
!> The state structure id may be different if WRF is coupled with another
!> model.
function get_wrf_domain(state_id)

integer, intent(in) :: state_id
integer :: get_wrf_domain

integer :: i

do i = 1, num_domains
   if (domain_id(i) == state_id) then
      get_wrf_domain = i
      return
   endif
enddo

call error_handler(E_ERR, 'get_wrf_domain', 'not a valid domain')

end function get_wrf_domain

!#######################################################################

function get_model_size()

integer(i8) :: get_model_size

get_model_size = wrf%model_size

end function get_model_size

!#######################################################################

function get_number_domains()

integer :: get_number_domains

get_number_domains = num_domains

end function get_number_domains

!#######################################################################

function get_wrf_static_data(dom_num)

integer, intent(in)  :: dom_num

type(wrf_static_data_for_dart) :: get_wrf_static_data

get_wrf_static_data = wrf%dom(dom_num)

return
end function get_wrf_static_data

!#######################################################################

subroutine get_wrf_state_variables(state_var)

character(len=129), intent(out) :: state_var(num_state_table_columns,max_state_variables)

state_var = wrf_state_variables

end subroutine get_wrf_state_variables

!#######################################################################

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations
integer :: model_dt, assim_dt

! We need to coordinate the desired assimilation window to be a 
! multiple of the model time step (which has no precision past integer seconds).

model_dt = nint(wrf%dom(1)%dt)

! The integer arithmetic does its magic.
assim_dt = (assimilation_period_seconds / model_dt) * model_dt

shortest_time_between_assimilations = set_time(assim_dt)

end function shortest_time_between_assimilations


!#######################################################################


subroutine get_state_meta_data(index_in, location, var_type_out, id_out)

! Given an integer index into the DART state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?

! this version has an optional last argument that is never called by
! any of the dart code, which can return the wrf domain number.

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type_out, id_out

integer     :: var_type, dart_type
integer(i8) :: index
integer     :: ip, jp, kp
integer     :: nz, ny, nx
logical     :: var_found
real(r8)    :: lon, lat, lev
character(len=129) :: string1

integer :: i, id, var_id, state_id
logical, parameter :: debug = .false.

! from the dart index get the local variables indices
call get_model_variable_indices(index_in, ip, jp, kp, var_id=var_id, dom_id=state_id)

! convert from state_structure domain number to wrf.
id = get_wrf_domain(state_id)

! at this point, (ip,jp,kp) refer to indices in the variable's own grid

if(debug) write(*,*) ' ip, jp, kp for index ',ip,jp,kp,index
if(debug) write(*,*) ' Var type: ',var_type

var_type  = wrf%dom(id)%var_type(var_id)
dart_type = wrf%dom(id)%dart_kind(var_id)

! first obtain lat/lon from (ip,jp)
call get_wrf_horizontal_location( ip, jp, var_type, id, lon, lat )

! return staggered levels for vertical types?  check this later.
if( (var_type == wrf%dom(id)%type_w ) .or. (var_type == wrf%dom(id)%type_gz) ) then
  lev = real(kp) - 0.5_r8
else
  lev = real(kp)
endif

if(debug) write(*,*) 'lon, lat, lev: ',lon, lat, lev

location = set_location(lon, lat, lev, VERTISLEVEL)

! return DART variable kind if requested
if(present(var_type_out)) var_type_out = dart_type

! return domain id if requested
if(present(id_out)) id_out = id

end subroutine get_state_meta_data

!--------------------------------------------------------------------
!> Distributed version of model interpolate
!> obs_kind is called as location type in assim model model
!> Storing the mean copy level location as the observation location
!> to save recomputation of v_p in vert_convert in get_close_obs
!
! Should this code be simplified so there is not so much repetition?
! This is the main forward operator subroutine for WRF.
! Given an ob (its DART location and kind), the corresponding model
! value is computed at nearest i,j,k. Thus, first i,j,k is obtained
! from ob lon,lat,z and then the state value that corresponds to
! the ob kind is interpolated.
!
! No location conversions are carried out in this subroutine. See
! get_close_obs, where ob vertical location information is converted
! to the requested vertical coordinate type.
subroutine model_interpolate(state_handle, ens_size, location, obs_kind, expected_obs, istatus)

! x:       Full DART state vector relevant to what's being updated
!          in the filter (mean or individual members).
! istatus: Returned 0 if everything is OK, >0 if error occured.
!                   1 = missing domain id
!                   2 = bad vertical level
!                   3 = unsupported obs kind
!                  10 = polar observation while restrict_polar namelist true
!                  99 = unknown reason (reached end of routine with missing_r8
!                       as obs_val)

! modified 26 June 2006 to accomodate vortex attributes
! modified 13 December 2006 to accomodate changes for the mpi version
! modified 22 October 2007 to accomodate global WRF (3.0)
! modified 18 November 2008 to allow unknown kinds to return without stopping
! modified  5 February 2010 to add circulation calc for vortex obs

! Helen Kershaw - Aim: to not require the whole state vector

! arguments
type(ensemble_type),   intent(in)  :: state_handle
integer,               intent(in)  :: ens_size
type(location_type),   intent(in)  :: location 
integer,               intent(in)  :: obs_kind
real(r8),              intent(out) :: expected_obs(ens_size)
integer,               intent(out) :: istatus(ens_size)

! local
logical, parameter  :: debug = .false.
logical, parameter  :: restrict_polar = .false.
logical, parameter  :: use_old_vortex = .true.   ! set to .false. to use circ code
real(r8), parameter :: drad = pi / 18.0_r8
real(r8)            :: xloc, yloc, xloc_u, yloc_v, xyz_loc(3)
integer             :: i, i_u, j, j_v, k2
real(r8)            :: dx,dy,dxm,dym,dx_u,dxm_u,dy_v,dym_v
integer             :: id
logical             :: surf_var
real(r8) :: a1(ens_size) 
real(r8) :: zloc(ens_size)
integer  :: k(ens_size)
real(r8) :: dz(ens_size), dzm(ens_size)
real(r8) :: utrue(ens_size), vtrue(ens_size)

! from getCorners
integer, dimension(2) :: ll, lr, ul, ur, ll_v, lr_v, ul_v, ur_v
integer            :: rc, i1, i2
integer(i8)        :: ill, ilr, iul, iur

real(r8) :: fld(2,ens_size)
real(r8), allocatable, dimension(:,:) :: v_h, v_p

! local vars, used in finding sea-level pressure and vortex center
real(r8), allocatable, dimension(:, :)    :: t1d, p1d, qv1d, z1d
real(r8), allocatable, dimension(:, :)    :: z1d_1, z1d_2
real(r8), allocatable, dimension(:,:)  :: pp, pd
real(r8), allocatable, dimension(:,:,:)  :: vfld
real(r8), allocatable, dimension(:,:,:)  :: uwnd, vwnd, vort
real(r8), allocatable, dimension(:)    :: x1d, y1d, xx1d, yy1d
integer  :: center_search_half_size, center_track_xmin, center_track_ymin, &
            center_track_xmax, center_track_ymax, circ_half_size, &
            circ_xmin, circ_xmax, circ_ymin, circ_ymax, xlen, ylen, &
            xxlen, yylen, ii1, ii2, cxlen, cylen, imax, jmax
real(r8) :: clat, clon, cxloc, cyloc, vcrit, &
            circ_half_length, asum, distgrid, dgi1, dgi2

real(r8) :: magwnd(ens_size), maxwspd(ens_size), circ(ens_size)

! local vars, used in calculating density, pressure, height
real(r8) :: rho1(ens_size) , rho2(ens_size) , rho3(ens_size), rho4(ens_size)
real(r8) :: pres1(ens_size), pres2(ens_size), pres3(ens_size), pres4(ens_size), pres(ens_size)
logical  :: is_lev0(ens_size) 

! local var for terrain elevation check for surface stations 
real(r8)            :: mod_sfc_elevation


real(r8) :: x_ill(ens_size), x_iul(ens_size), x_ilr(ens_size), x_iur(ens_size), ugrid(ens_size), vgrid(ens_size)
real(r8) :: x_ugrid_1(ens_size), x_ugrid_2(ens_size), x_ugrid_3(ens_size), x_ugrid_4(ens_size)
real(r8) :: x_vgrid_1(ens_size), x_vgrid_2(ens_size), x_vgrid_3(ens_size), x_vgrid_4(ens_size)
integer  :: e, count, uk !< index varibles for loop
real(r8) :: failedcopies(ens_size)
integer, allocatable  :: uniquek(:)
integer  :: ksort(ens_size)

integer(i8) :: ugrid_1, ugrid_2, ugrid_3, ugrid_4, vgrid_1, vgrid_2, vgrid_3, vgrid_4
integer(i8) :: z1d_ind1, z1d_ind2, t1d_ind, qv1d_ind

id = 1

! call static_data_sizes(domain=id)

! Initialize stuff
istatus(:) = 0
fld(:,:) = missing_r8
expected_obs(:) = missing_r8  !< array of obs_vals
failedcopies(:) = 1

! If identity observation (obs_kind < 0), then no need to interpolate
if ( obs_kind < 0 ) then

   ! identity observation -> -(obs_kind)=DART state vector index
   ! obtain state value directly from index
   expected_obs(:) = get_state(int(-1*obs_kind, i8), state_handle)

! Otherwise, we need to do interpolation
else

   ! Is this a valid kind to interpolate?  Set up in the static_init_model code,
   ! based on entries in wrf_state_vector namelist item.
   if (.not. in_state_vector(obs_kind)) then
      write(errstring, *) 'cannot interpolate ' // trim(get_name_for_quantity(obs_kind)) &
                           // ' with the current WRF arrays in state vector'
      call error_handler(E_ERR, 'model_interpolate', errstring, &
                                 source, revision, revdate)
   endif

   ! Unravel location_type information
   xyz_loc = get_location(location)

   !----------------------------------
   ! 0. Prelude to Interpolation
   !----------------------------------
   
   ! 0.a Horizontal stuff

   ! first obtain domain id, and mass points (i,j)
! JPH --- scm is only defined for d1
   if ( .not. scm ) then
      call get_domain_info(xyz_loc(1),xyz_loc(2),id,xloc,yloc)
   else
      id = 1
      xloc = 1.0_r8
      yloc = 1.0_r8
   endif

   ! check that we obtained a valid domain id number
   if (id==0) then
      istatus = 1
      return
   endif

   if ( debug ) then
     write(*,*) 'retreiving obs kind ',obs_kind,' on domain ',id
   endif

   !*****************************************************************************
   ! Check polar-b.c. constraints -- if restrict_polar = .true., then we are not 
   !   processing observations poleward of the 1st or last mass grid points.
   ! If we have tried to pass a polar observation, then exit with istatus = 10
   if ( wrf%dom(id)%polar .and. restrict_polar ) then
      if ( yloc < 1.0_r8 .or. yloc >= real(wrf%dom(id)%sn,r8) ) then

         ! Perhaps write to dart_log.out?
         write(errstring,*)'Obs cannot be polar with restrict_polar on: yloc = ',yloc
         call error_handler(E_WARN,'model_interpolate', errstring, &
              source, revision, revdate)

         istatus = 10  ! istatus 10, if it's not used, will mean the observation is too polar
         print*, 'model_mod.f90 :: model_interpolate :: No polar observations!  istatus = ', istatus
         return
      endif
   endif
   !*****************************************************************************
   
   ! print info if debugging
   if(debug) then
      i = xloc
      j = yloc
      print*,xyz_loc(2), xyz_loc(1), xloc,yloc
      write(*,*) ' corners of lat '
      write(*,*) wrf%dom(id)%latitude(i,j),wrf%dom(id)%latitude(i+1,j),  &
           wrf%dom(id)%latitude(i,j+1), &
           wrf%dom(id)%latitude(i+1,j+1)
      write(*,*) ' corners of long '
      write(*,*) wrf%dom(id)%longitude(i,j),wrf%dom(id)%longitude(i+1,j),  &
           wrf%dom(id)%longitude(i,j+1), &
           wrf%dom(id)%longitude(i+1,j+1)
   endif
   
   ! get integer (west/south) grid point and distances to neighboring grid points
   ! distances are used as weights to carry out horizontal interpolations
   call toGrid(xloc,i,dx,dxm)
   call toGrid(yloc,j,dy,dym)
   
   ! 0.b Vertical stuff

   if ( debug ) then
      write(*,*) 'is_vertical(PRESSURE) ',is_vertical(location,"PRESSURE")
      write(*,*) 'is_vertical(HEIGHT) ',is_vertical(location,"HEIGHT")
      write(*,*) 'is_vertical(LEVEL) ',is_vertical(location,"LEVEL")
      write(*,*) 'is_vertical(SURFACE) ',is_vertical(location,"SURFACE")
      write(*,*) 'is_vertical(UNDEFINED) ',is_vertical(location,"UNDEFINED")
   endif

   ! HK
   ! Allocate both a vertical height and vertical pressure coordinate -- 0:bt x ens_size
   allocate(v_h(0:wrf%dom(id)%bt, ens_size), v_p(0:wrf%dom(id)%bt, ens_size))

   ! Set surf_var to .false. and then change in vert_is_surface section if necessary
   surf_var = .false.

   ! Determine corresponding model level for obs location
   ! This depends on the obs vertical coordinate
   !   From this we get a meaningful z-direction real-valued index number
   if(is_vertical(location,"LEVEL")) then
      ! Ob is by model level
      zloc = xyz_loc(3)

   elseif(is_vertical(location,"PRESSURE")) then
      ! Ob is by pressure: get corresponding mass level zloc from
      ! computed column pressure profile
      call get_model_pressure_profile_distrib(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,id,v_p,state_handle, ens_size)

      !print*, 'v_p distrib ', v_p

! testing get_model_pressure_profile_distrib
!      istatus = 0
!      expected_obs = -9
!      return

      ! get pressure vertical co-ordinate
      call pres_to_zk_distrib(xyz_loc(3), v_p, wrf%dom(id)%bt, ens_size, zloc,is_lev0)
      if(debug .and. obs_kind /= QTY_SURFACE_PRESSURE) &
                print*,' obs is by pressure and zloc,lev0 =',zloc, is_lev0
      if(debug) print*,'model pressure profile'
      if(debug) print*,v_p
      
      ! If location is above model surface but below the lowest sigma level,
      ! the default is to reject it.  But if the namelist value is true, then
      ! accept the observation and later on extrapolate the values from levels
      ! 1 and 2 downward.

      !HK ensemble loop, however do you reject the obs for all ensembles later anyway?

      do e = 1, ens_size

         if (is_lev0(e)) then
            ! the pres_to_zk() routine has returned a valid zloc in case we
            ! want to use it.  the default is to reject the observation and so
            ! we overwrite it with missing -- but, if the namelist value is set
            ! to true, leave zloc alone.
            if (.not. allow_obs_below_vol) zloc(e) = missing_r8
               if (debug .and. .not. allow_obs_below_vol) print*, 'setting zloc missing'

            ! else need to set a qc here?
         endif

      enddo

   elseif(is_vertical(location,"HEIGHT")) then

      ! Ob is by height: get corresponding mass level zloc from
      ! computed column height profile
      call get_model_height_profile_distrib(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,id,v_h, state_handle, ens_size)
      ! get height vertical co-ordinate
      do e = 1, ens_size ! HK should there be a height_to_zk_distrib?
         call height_to_zk(xyz_loc(3), v_h(:, e), wrf%dom(id)%bt,zloc(e),is_lev0(e))
         if(debug) print*,' obs is by height and zloc,lev0 =',zloc(e), is_lev0(e)
         if(debug) print*,'model height profile'
         if(debug) print*,v_h

         ! If location is above model surface but below the lowest sigma level,
         ! the default is to reject it.  But if the namelist value is true, then
         ! accept the observation and later on extrapolate the values from levels
         ! 1 and 2 downward.
         if (is_lev0(e)) then
            ! the height_to_zk() routine has returned a valid zloc in case we
            ! want to use it.  the default is to reject the observation and so
            ! we overwrite it with missing.  but if the namelist value is set
            ! to true, leave zloc alone.
            if (.not. allow_obs_below_vol) zloc(e) = missing_r8
            if (debug .and. .not. allow_obs_below_vol) print*, 'setting zloc missing member ', e
            ! else need to set a qc here?
         endif

      enddo

   elseif(is_vertical(location,"SURFACE") .or. obs_kind == QTY_SURFACE_ELEVATION) then
      zloc = 1.0_r8
      surf_var = .true.
      if(debug) print*,' obs is at the surface = ', xyz_loc(3)

      ! if you want to have a distance check to see if the station height
      ! is too far away from the model surface height, here is the place to
      ! reject the observation.

      ! Elevation check function drawn from Ryan's wrf preprocessing code. 
      ! The elevation is now passed in instead of calling model_interpolate:
      ! sfc_elev_max_diff  - if < 0 routine is skipped.
      ! mod_sfc_elevation  - interpolated model surface height at the lowest
      !                    - model layer.
      ! z_loc - the third array element of xyz_loc
      ! the station elevation against the estimated model surface height at the
      ! station location, the maximum difference in elevation allowed (m), and the
      ! observation location. There is no check for whether the third element for
      ! the xyz_loc array is anything other than station height in meters. The
      ! function returns a logical where .true. means the obs station elevation
      ! 'passed' the height check. Here, if a height check fails, set an istatus
      ! of '1' and bail out of this routine. 
      if ( sfc_elev_max_diff >= 0 ) then
      ! Check to make sure retrieved integer gridpoints are in valid range
         if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
              boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) ) then

            call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
            if ( rc .ne. 0 ) &
                print*, 'model_mod.f90 :: model_interpolate :: getCorners HGT for sfc rc = ', rc

            ! Interpolation for the HGT field -- HGT is NOT part of state vector x, but rather
            !   in the associated domain meta data
            mod_sfc_elevation = dym*( dxm*wrf%dom(id)%hgt(ll(1), ll(2)) + &
                                       dx*wrf%dom(id)%hgt(lr(1), lr(2)) ) + &
                                 dy*( dxm*wrf%dom(id)%hgt(ul(1), ul(2)) + &
                                       dx*wrf%dom(id)%hgt(ur(1), ur(2)) )
         endif
         if ( .not. height_diff_check(sfc_elev_max_diff,xyz_loc(3),mod_sfc_elevation) ) zloc = missing_r8
      endif

   elseif(is_vertical(location,"UNDEFINED")) then
      ! the zloc value should not be used since there is no actual vertical
      ! location for this observation, but give zloc a valid value to avoid
      ! the error checks below for missing_r8
      zloc  = 0.0_r8
      if(debug) print*,' obs height is intentionally undefined'

   else
      write(errstring,*) 'wrong option for which_vert ', &
                         nint(query_location(location,'which_vert'))
      call error_handler(E_ERR,'model_interpolate', errstring, &
           source, revision, revdate)

   endif


   ! Deal with missing vertical coordinates -- return with istatus .ne. 0
   ! HK This is annoying.  Back to earlier question of QC if one ensemble fails do we
   ! discard all the obervations?
   ! For now, to replicate the original code behaviour keep track of the failed copies

   ! ********* bail out early code **********
   !if(any(zloc == missing_r8)) then
   !   print*, 'One of the ensembles has failed the vertical co ordinate rank', my_task_id()
   !   expected_obs = missing_r8 !HK this is redundant? expected_obs is initialized to missing_r8
   !   istatus = 2
   !   deallocate(v_h, v_p)
   !   return
   !endif
   ! ********* endof bail out early code ********

   do e = 1, ens_size
      if(zloc(e) == missing_r8) then
          failedcopies(e) = missing_r8
          istatus(e) = 2
      endif
   enddo

   ! Set a working integer k value -- if (int(zloc) < 1), then k = 1
   k = max(1,int(zloc)) !HK k is now ensemble size


   ! Find the unique k values
   ksort = sort(k)

   count = 1
   do e = 2, ens_size
       if ( ksort(e) /= ksort(e-1) ) count = count + 1
   enddo

   allocate(uniquek(count))
   uniquek(:) = -1 
   uk = 1
   do e = 1, ens_size
      if ( all(uniquek /= k(e)) ) then
         uniquek(uk) = k(e)
         uk = uk + 1
      endif
   enddo

   ! The big horizontal interp loop below computes the data values in the level
   ! below and above the actual location, and then does a separate vertical
   ! interpolation (if the obs is not a 2d surface obs).  The two values are
   ! stored in fld(1:2).  Set them to missing here, and if the code below cannot
   ! compute a value, it can just drop out and not have to explicitly set it to
   ! missing anymore.
   fld(:,:) = missing_r8 !HK

   !----------------------------------
   ! 1. Horizontal Interpolation 
   !----------------------------------

   ! This part is the forward operator -- compute desired model state value for given point.

   ! Strategy is to do the horizontal interpolation on two different levels in the
   !   vertical, and then to do the vertical interpolation afterwards, since it depends on
   !   what the vertical coordinate is

   ! Large if-structure to select on obs_kind of desired field....
   ! Table of Contents:
   ! 1.a Horizontal Winds (U, V, U10, V10)
   ! 1.b Sensible Temperature (T, T2)
   ! 1.c Potential Temperature (Theta, TH2)
   ! 1.d Density (Rho)
   ! 1.e Vertical Wind (W)
   ! 1.f Specific Humidity (SH, SH2)
   ! 1.g Vapor Mixing Ratio (QV, Q2)
   ! 1.h Rainwater Mixing Ratio (QR)
   ! 1.i.1 Graupel Mixing Ratio (QG)
   ! 1.i.2 Hail Mixing Ratio (QH)
   ! 1.j Snow Mixing Ratio (QS)
   ! 1.k Ice Mixing Ratio (QI)
   ! 1.l Cloud Mixing Ratio (QC)
   ! 1.m Droplet Number Concentration (QNDRP)
   ! 1.n Ice Number Concentration (QNICE)
   ! 1.o Snow Number Concentration (QNSNOW)
   ! 1.p Rain Number Concentration (QNRAIN)
   ! 1.q.1 Graupel Number Concentration (QNGRAUPEL) 
   ! 1.q.2 Hail Number Concentration (QNHAIL)
   ! 1.r Previous time step condensational heating (H_DIABATIC)
   ! 1.s Reflectivity weighted precip fall speed (FALL_SPD_Z_WEIGHTED)
   ! 1.t Pressure (P)
   ! 1.u Vortex Center Stuff from Yongsheng
   ! 1.v.1 Radar Reflectivity (REFL_10CM)
   ! 1.v.2 Differential Reflectivity (DIFF_REFL_10CM)
   ! 1.v.3 Specific Differential Phase (SPEC_DIFF_10CM)
   ! 1.w Geopotential Height (GZ)
   ! 1.x Surface Elevation (HGT)
   ! 1.y Surface Skin Temperature (TSK)
   ! 1.z Land Mask (XLAND)

   ! NEWVAR:  Currently you have to add a new case here to tell the code what
   !   field inside the state vector you will be interpolating in.  the eventual
   !   plan is for there to be a default case which all simple interps fall into,
   !   but for now we still have to add code.

   ! NOTE: the previous version of this code checked for surface observations with the syntax:
   !          "if(.not. vert_is_surface(location) .or. .not. surf_var) then"
   !   We identified this as redundant because surf_var is changed from .false. only by
   !     the above code (section 0.b), which must be traced through before one can arrive
   !     at the following forward operator code.  Hence, we can remove the call to 
   !     vert_is_surface.

   if (obs_kind == QTY_RAINWATER_MIXING_RATIO .or. & 
       obs_kind == QTY_GRAUPEL_MIXING_RATIO .or. &
       obs_kind == QTY_HAIL_MIXING_RATIO .or. &
       obs_kind == QTY_SNOW_MIXING_RATIO .or. &
       obs_kind == QTY_ICE_MIXING_RATIO .or. &
       obs_kind == QTY_CLOUDWATER_MIXING_RATIO .or. &
       obs_kind == QTY_DROPLET_NUMBER_CONCENTR .or. &
       obs_kind == QTY_ICE_NUMBER_CONCENTRATION .or. &
       obs_kind == QTY_SNOW_NUMBER_CONCENTR .or. &
       obs_kind == QTY_RAIN_NUMBER_CONCENTR .or. &
       obs_kind == QTY_GRAUPEL_NUMBER_CONCENTR .or. &
       obs_kind == QTY_HAIL_NUMBER_CONCENTR .or. &
       obs_kind == QTY_CONDENSATIONAL_HEATING .or. &
       obs_kind == QTY_POWER_WEIGHTED_FALL_SPEED .or. &
       obs_kind == QTY_RADAR_REFLECTIVITY .or. &
       obs_kind == QTY_DIFFERENTIAL_REFLECTIVITY .or. &
       obs_kind == QTY_SPECIFIC_DIFFERENTIAL_PHASE .or. &
       obs_kind == QTY_CLOUD_FRACTION) then

       call simple_interp_distrib(fld, wrf, id, i, j, k, obs_kind, dxm, dx, dy, dym, uniquek, ens_size, state_handle)
       if (all(fld == missing_r8)) goto 200
      

      ! don't accept negative fld
      if (obs_kind == QTY_RAINWATER_MIXING_RATIO .or. &
          obs_kind == QTY_GRAUPEL_MIXING_RATIO .or. &
          obs_kind == QTY_HAIL_MIXING_RATIO .or. &
          obs_kind == QTY_SNOW_MIXING_RATIO .or. &
          obs_kind == QTY_ICE_MIXING_RATIO .or. &
          obs_kind == QTY_CLOUDWATER_MIXING_RATIO .or. &
          obs_kind == QTY_DROPLET_NUMBER_CONCENTR .or. &
          obs_kind == QTY_ICE_NUMBER_CONCENTRATION .or. &
          obs_kind == QTY_SNOW_NUMBER_CONCENTR .or. &
          obs_kind == QTY_RAIN_NUMBER_CONCENTR .or. &
          obs_kind == QTY_GRAUPEL_NUMBER_CONCENTR .or. &
          obs_kind == QTY_HAIL_NUMBER_CONCENTR ) then

          fld = max(0.0_r8, fld) ! Don't accept negative

      endif

   !-----------------------------------------------------
   ! 1.a Horizontal Winds (U, V, U10, V10)

   ! We need one case structure for both U & V because they comprise a vector which could need
   !   transformation depending on the map projection (hence, the call to gridwind_to_truewind)
   elseif( obs_kind == QTY_U_WIND_COMPONENT .or. obs_kind == QTY_V_WIND_COMPONENT) then   ! U, V

     ! This is for 3D wind fields -- surface winds later
      if(.not. surf_var) then

         if ( ( wrf%dom(id)%type_u >= 0 ) .and. ( wrf%dom(id)%type_v >= 0 ) ) then

            ! xloc and yloc are indices on mass-grid.  If we are on a periodic longitude domain,
            !   then xloc can range from [1 wes).  This means that simply adding 0.5 to xloc has
            !   the potential to render xloc_u out of the valid mass-grid index bounds (>wes).
            !   To remedy this, we can either do periodicity check on xloc_u value, or we can
            !   leave it to a subroutine or function to alter xloc to xloc_u if the observation
            !   type requires it.
            xloc_u = xloc + 0.5_r8
            yloc_v = yloc + 0.5_r8
   
            ! Check periodicity if necessary -- but only subtract 'we' because the U-grid
            !   cannot have an index < 1 (i.e., U(wes) = U(1) ).
            if ( wrf%dom(id)%periodic_x .and. xloc_u > real(wrf%dom(id)%wes,r8) ) &
                 xloc_u = xloc_u - real(wrf%dom(id)%we,r8)
   
            ! Get South West gridpoint indices for xloc_u and yloc_v
            call toGrid(xloc_u,i_u,dx_u,dxm_u)
            call toGrid(yloc_v,j_v,dy_v,dym_v)

            do uk = 1, count ! for the different ks

               ! Check to make sure retrieved integer gridpoints are in valid range
               if ( boundsCheck( i_u, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_u) .and. &
                    boundsCheck( i,   wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_v) .and. &
                    boundsCheck( j,   wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_u) .and. &
                    boundsCheck( j_v, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_v) .and. &
                    boundsCheck( uniquek(uk),   .false.,                id, dim=3, type=wrf%dom(id)%type_u) ) then
   
                  ! Need to get grid cell corners surrounding observation location -- with
                  !   periodicity, this could be non-consecutive (i.e., NOT necessarily i and i+1);
                  !   Furthermore, it could be different for the U-grid and V-grid.  Remember, for
                  !   now, we are disallowing observations to be located poleward of the 1st and
                  !   last mass points.
               
                  call getCorners(i_u, j, id, wrf%dom(id)%type_u, ll, ul, lr, ur, rc )
                  if ( rc .ne. 0 ) &
                       print*, 'model_mod.f90 :: model_interpolate :: getCorners U rc = ', rc
               
                  call getCorners(i, j_v, id, wrf%dom(id)%type_v, ll_v, ul_v, lr_v, ur_v, rc )
                  if ( rc .ne. 0 ) &
                       print*, 'model_mod.f90 :: model_interpolate :: getCorners V rc = ', rc

                  ! Now we want to get the corresponding DART state vector indices, and then
                  !   interpolate horizontally on TWO different vertical levels (so that we can
                  !   do the vertical interpolation properly later)
                  do k2 = 1, 2

                     ! Interpolation for the U field
                     ill = get_dart_vector_index(ll(1), ll(2), uniquek(uk)+k2-1, domain_id(id), wrf%dom(id)%type_u)
                     iul = get_dart_vector_index(ul(1), ul(2), uniquek(uk)+k2-1, domain_id(id), wrf%dom(id)%type_u)
                     ilr = get_dart_vector_index(lr(1), lr(2), uniquek(uk)+k2-1, domain_id(id), wrf%dom(id)%type_u)
                     iur = get_dart_vector_index(ur(1), ur(2), uniquek(uk)+k2-1, domain_id(id), wrf%dom(id)%type_u)

                     x_ill = get_state(ill, state_handle)
                     x_iul = get_state(iul, state_handle)
                     x_ilr = get_state(ilr, state_handle)
                     x_iur = get_state(iur, state_handle)

                     ugrid = dym*( dxm_u*x_ill + dx_u*x_ilr ) + dy*( dxm_u*x_iul + dx_u*x_iur )

                     ! Interpolation for the V field
                     ill = get_dart_vector_index(ll_v(1), ll_v(2), uniquek(uk)+k2-1, domain_id(id), wrf%dom(id)%type_v)
                     iul = get_dart_vector_index(ul_v(1), ul_v(2), uniquek(uk)+k2-1, domain_id(id), wrf%dom(id)%type_v)
                     ilr = get_dart_vector_index(lr_v(1), lr_v(2), uniquek(uk)+k2-1, domain_id(id), wrf%dom(id)%type_v)
                     iur = get_dart_vector_index(ur_v(1), ur_v(2), uniquek(uk)+k2-1, domain_id(id), wrf%dom(id)%type_v)

                     x_ill = get_state(ill, state_handle)
                     x_iul = get_state(iul, state_handle)
                     x_ilr = get_state(ilr, state_handle)
                     x_iur = get_state(iur, state_handle)

                     vgrid = dym_v*( dxm*x_ill + dx*x_ilr ) + dy_v*( dxm*x_iul + dx*x_iur )

                        do e = 1, ens_size

                           if ( k(e) == uniquek(uk) ) then ! interpolate only if is the correct k
                              ! Certain map projections have wind on grid different than true wind (on map)
                              !   subroutine gridwind_to_truewind is in module_map_utils.f90
                              call gridwind_to_truewind(xyz_loc(1), wrf%dom(id)%proj, ugrid(e), vgrid(e), &
                                utrue(e), vtrue(e))
                  
                              ! Figure out which field was the actual desired observation and store that
                              !   field as one of the two elements of "fld" (the other element is the other
                              !   k-level)
                              if( obs_kind == QTY_U_WIND_COMPONENT) then
                                 fld(k2, e) = utrue(e)
                              else   ! must want v
                                 fld(k2, e) = vtrue(e)
                              endif
                           endif
                        enddo
                    enddo
                   endif

               enddo

          endif


      ! This is for surface wind fields -- NOTE: surface winds are on Mass grid 
      ! (therefore, TYPE_T), not U-grid & V-grid.  
      ! Also, because surface winds are at a given single vertical level, 
      !  only fld(1) will be filled.
      else

         if ( ( wrf%dom(id)%type_u10 >= 0 ) .and. ( wrf%dom(id)%type_v10 >= 0 ) ) then

   ! JPH -- should test this for doubly periodic
   ! JPH -- does not pass for SCM config, so just do it below
            ! Check to make sure retrieved integer gridpoints are in valid range
            if ( ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
                   boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) )     &
                   .or. wrf%dom(id)%scm ) then
   
               call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
               if ( rc .ne. 0 ) &
                    print*, 'model_mod.f90 :: model_interpolate :: getCorners U10, V10 rc = ', rc
   
               ! Interpolation for the U10 field
               ill = get_dart_vector_index(ll(1), ll(2), 1, domain_id(id), wrf%dom(id)%type_u10)
               iul = get_dart_vector_index(ul(1), ul(2), 1, domain_id(id), wrf%dom(id)%type_u10)
               ilr = get_dart_vector_index(lr(1), lr(2), 1, domain_id(id), wrf%dom(id)%type_u10)
               iur = get_dart_vector_index(ur(1), ur(2), 1, domain_id(id), wrf%dom(id)%type_u10)

               x_ill = get_state(ill, state_handle)
               x_iul = get_state(iul, state_handle)
               x_ilr = get_state(ilr, state_handle)
               x_iur = get_state(iur, state_handle)

               ugrid = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )
   
               ! Interpolation for the V10 field
               ill = get_dart_vector_index(ll(1), ll(2), 1, domain_id(id), wrf%dom(id)%type_v10)
               iul = get_dart_vector_index(ul(1), ul(2), 1, domain_id(id), wrf%dom(id)%type_v10)
               ilr = get_dart_vector_index(lr(1), lr(2), 1, domain_id(id), wrf%dom(id)%type_v10)
               iur = get_dart_vector_index(ur(1), ur(2), 1, domain_id(id), wrf%dom(id)%type_v10)

               x_ill = get_state(ill, state_handle)
               x_iul = get_state(iul, state_handle)
               x_iur = get_state(iur, state_handle)
               x_ilr = get_state(ilr, state_handle)

               vgrid = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )

               do e = 1, ens_size
                  call gridwind_to_truewind(xyz_loc(1), wrf%dom(id)%proj, ugrid(e), vgrid(e), &
                       utrue(e), vtrue(e))

                  ! U10 (U at 10 meters)
                  if( obs_kind == QTY_U_WIND_COMPONENT) then
                     fld(1, e) = utrue(e)
                  ! V10 (V at 10 meters)
                  else
                     fld(1, e) = vtrue(e)
                  endif

               enddo


            endif
         endif
      endif

   elseif (obs_kind == QTY_10M_U_WIND_COMPONENT .or. obs_kind == QTY_10M_V_WIND_COMPONENT) then
      if ( ( wrf%dom(id)%type_u10 >= 0 ) .and. ( wrf%dom(id)%type_v10 >= 0 ) ) then

   ! JPH -- should test this for doubly periodic
   ! JPH -- does not pass for SCM config, so just do it below
         ! Check to make sure retrieved integer gridpoints are in valid range
         if ( ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
                boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) )     &
                .or. wrf%dom(id)%scm ) then

            call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
            if ( rc .ne. 0 ) &
                 print*, 'model_mod.f90 :: model_interpolate :: getCorners U10, V10 rc = ', rc

            ! Interpolation for the U10 field
            ill = get_dart_vector_index(ll(1), ll(2), 1, domain_id(id), wrf%dom(id)%type_u10)
            iul = get_dart_vector_index(ul(1), ul(2), 1, domain_id(id), wrf%dom(id)%type_u10)
            ilr = get_dart_vector_index(lr(1), lr(2), 1, domain_id(id), wrf%dom(id)%type_u10)
            iur = get_dart_vector_index(ur(1), ur(2), 1, domain_id(id), wrf%dom(id)%type_u10)

            x_ill = get_state(ill, state_handle)
            x_iul = get_state(iul, state_handle)
            x_ilr = get_state(ilr, state_handle)
            x_iur = get_state(iur, state_handle)

            ugrid = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )

            ! Interpolation for the V10 field
            ill = get_dart_vector_index(ll(1), ll(2), 1, domain_id(id), wrf%dom(id)%type_v10)
            iul = get_dart_vector_index(ul(1), ul(2), 1, domain_id(id), wrf%dom(id)%type_v10)
            ilr = get_dart_vector_index(lr(1), lr(2), 1, domain_id(id), wrf%dom(id)%type_v10)
            iur = get_dart_vector_index(ur(1), ur(2), 1, domain_id(id), wrf%dom(id)%type_v10)

            x_ill = get_state(ill, state_handle)
            x_iul = get_state(iul, state_handle)
            x_iur = get_state(iur, state_handle)
            x_ilr = get_state(ilr, state_handle)

            vgrid = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )

            do e = 1, ens_size
               call gridwind_to_truewind(xyz_loc(1), wrf%dom(id)%proj, ugrid(e), vgrid(e), &
                    utrue(e), vtrue(e))

               ! U10 (U at 10 meters)
               if( obs_kind == QTY_10M_U_WIND_COMPONENT) then
                  fld(1, e) = utrue(e)
               ! V10 (V at 10 meters)
               else
                  fld(1, e) = vtrue(e)
               endif
            enddo
         endif
      endif
   !-----------------------------------------------------
   ! 1.b Sensible Temperature (T, T2)

   elseif ( obs_kind == QTY_TEMPERATURE ) then
      ! This is for 3D temperature field -- surface temps later
      !print*, 'k ', k

      if ( wrf%dom(id)%type_t >= 0 ) then

         do uk = 1, count ! for the different ks

            ! Check to make sure retrieved integer gridpoints are in valid range
            if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
                 boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) .and. &
                 boundsCheck( uniquek(uk), .false.,                id, dim=3, type=wrf%dom(id)%type_t ) ) then

               call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
               if ( rc .ne. 0 ) &
                    print*, 'model_mod.f90 :: model_interpolate :: getCorners T rc = ', rc
            
               ! Interpolation for T field at level k
               ill = get_dart_vector_index(ll(1), ll(2), uniquek(uk), domain_id(id), wrf%dom(id)%type_t)
               iul = get_dart_vector_index(ul(1), ul(2), uniquek(uk), domain_id(id), wrf%dom(id)%type_t)
               ilr = get_dart_vector_index(lr(1), lr(2), uniquek(uk), domain_id(id), wrf%dom(id)%type_t)
               iur = get_dart_vector_index(ur(1), ur(2), uniquek(uk), domain_id(id), wrf%dom(id)%type_t)

               x_iul = get_state(iul, state_handle)
               x_ill = get_state(ill, state_handle)
               x_ilr = get_state(ilr, state_handle)
               x_iur = get_state(iur, state_handle)

               ! In terms of perturbation potential temperature
               a1 = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )

               pres1 = model_pressure_t_distrib(ll(1), ll(2), uniquek(uk), id, state_handle, ens_size)
               pres2 = model_pressure_t_distrib(lr(1), lr(2), uniquek(uk), id, state_handle, ens_size)
               pres3 = model_pressure_t_distrib(ul(1), ul(2), uniquek(uk), id, state_handle, ens_size)
               pres4 = model_pressure_t_distrib(ur(1), ur(2), uniquek(uk), id, state_handle, ens_size)

               ! Pressure at location
               pres = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )

               do e = 1, ens_size
                  if ( k(e) == uniquek(uk) ) then
                     ! Full sensible temperature field
                     fld(1, e) = (ts0 + a1(e))*(pres(e)/ps0)**kappa
                  endif
               enddo

               ! Interpolation for T field at level k+1
               ill = get_dart_vector_index(ll(1), ll(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_t)
               iul = get_dart_vector_index(ul(1), ul(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_t)
               ilr = get_dart_vector_index(lr(1), lr(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_t)
               iur = get_dart_vector_index(ur(1), ur(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_t)

               x_ill = get_state(ill, state_handle)
               x_iul = get_state(iul, state_handle)
               x_iur = get_state(iur, state_handle)
               x_ilr = get_state(ilr, state_handle)

               ! In terms of perturbation potential temperature
               a1 = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )

               pres1 = model_pressure_t_distrib(ll(1), ll(2), uniquek(uk)+1, id, state_handle, ens_size)
               pres2 = model_pressure_t_distrib(lr(1), lr(2), uniquek(uk)+1, id, state_handle, ens_size)
               pres3 = model_pressure_t_distrib(ul(1), ul(2), uniquek(uk)+1, id, state_handle, ens_size)
               pres4 = model_pressure_t_distrib(ur(1), ur(2), uniquek(uk)+1, id, state_handle, ens_size)

               ! Pressure at location
               pres = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )

               do e = 1, ens_size
                  if ( k(e) == uniquek(uk) ) then
                  ! Full sensible temperature field
                  fld(2, e) = (ts0 + a1(e))*(pres(e)/ps0)**kappa
                  endif
               enddo
            endif
         enddo
      else
         fld = missing_r8
      end if
   elseif (obs_kind == QTY_2M_TEMPERATURE) then
      ! This is for 2-meter temperature
      if ( wrf%dom(id)%type_t2 >= 0 ) then ! HK is there a better way to do this?
         call surface_interp_distrib(fld, wrf, id, i, j, obs_kind, wrf%dom(id)%type_t2, dxm, dx, dy, dym, ens_size, state_handle)
         if (all(fld == missing_r8)) goto 200
      else
         fld = missing_r8
      endif

   !-----------------------------------------------------
   ! 1.c Potential Temperature (Theta, TH2)

   ! Note:  T is perturbation potential temperature (potential temperature - ts0)
   !   TH2 is potential temperature at 2 m
   elseif ( obs_kind == QTY_POTENTIAL_TEMPERATURE ) then
      ! This is for 3D potential temperature field -- surface pot temps later
      if(.not. surf_var) then

         if ( wrf%dom(id)%type_t >= 0 ) then

            do uk = 1, count

            ! Check to make sure retrieved integer gridpoints are in valid range
            if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
                 boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) .and. &
                 boundsCheck( uniquek(uk), .false.,                id, dim=3, type=wrf%dom(id)%type_t ) ) then
         
               call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
               if ( rc .ne. 0 ) &
                    print*, 'model_mod.f90 :: model_interpolate :: getCorners Theta rc = ', rc
               
               ! Interpolation for Theta field at level k
               ill = get_dart_vector_index(ll(1), ll(2), uniquek(uk),  domain_id(id),wrf%dom(id)%type_t)
               iul = get_dart_vector_index(ul(1), ul(2), uniquek(uk),  domain_id(id),wrf%dom(id)%type_t)
               ilr = get_dart_vector_index(lr(1), lr(2), uniquek(uk),  domain_id(id),wrf%dom(id)%type_t)
               iur = get_dart_vector_index(ur(1), ur(2), uniquek(uk),  domain_id(id),wrf%dom(id)%type_t)

               x_ill = get_state(ill, state_handle)
               x_iul = get_state(iul, state_handle)
               x_ilr = get_state(ilr, state_handle)
               x_iur = get_state(iur, state_handle)

               do e = 1, ens_size
                  if ( k(e) == uniquek(uk) ) then
                     fld(1, e) = ts0 + dym*( dxm*x_ill(e) + dx*x_ilr(e)) + dy*( dxm*x_iul(e) + dx*x_iur(e) )
                  endif
               enddo
   
               ! Interpolation for Theta field at level k+1
               ill = get_dart_vector_index(ll(1), ll(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_t)
               iul = get_dart_vector_index(ul(1), ul(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_t)
               ilr = get_dart_vector_index(lr(1), lr(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_t)
               iur = get_dart_vector_index(ur(1), ur(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_t)

               x_ill = get_state(ill, state_handle)
               x_ill = get_state(ill, state_handle)
               x_ilr = get_state(ilr, state_handle)
               x_iur = get_state(iur, state_handle)

               do e = 1, ens_size
                  if ( k(e) == uniquek(uk) ) then
                    fld(2, e) = ts0 + dym*( dxm*x_ill(e) + dx*x_ilr(e) ) + dy*( dxm*x_iul(e) + dx*x_iur(e) )
                 endif
               enddo
             endif
           enddo
         endif

      ! This is for surface potential temperature (TH2)
      else
         
         if ( wrf%dom(id)%type_th2 >= 0 ) then

            call surface_interp_distrib(fld, wrf, id, i, j, obs_kind, wrf%dom(id)%type_th2, dxm, dx, dy, dym, ens_size, state_handle)
            if (all(fld == missing_r8)) goto 200
   
            endif
      endif

   !-----------------------------------------------------
   ! 1.d Density (Rho)
   elseif (obs_kind == QTY_DENSITY) then

      do uk = 1, count ! for the different ks

      ! Check to make sure retrieved integer gridpoints are in valid range
      if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
           boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) .and. &
           boundsCheck( uniquek(uk), .false.,                id, dim=3, type=wrf%dom(id)%type_t ) ) then
         
         call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
         if ( rc .ne. 0 ) &
              print*, 'model_mod.f90 :: model_interpolate :: getCorners Rho rc = ', rc
      
         ! calculate full rho at corners of interp box
         ! and interpolate to desired horizontal location

         ! Hmmm, it does not appear that Rho is part of the DART state vector, so there
         !   is not a reference to wrf%dom(id)%dart_ind -- we'll have to go right from
         !   the corner indices

         ! Interpolation for the Rho field at level k
         rho1 = model_rho_t_distrib(ll(1), ll(2), uniquek(uk), id, state_handle, ens_size)
         rho2 = model_rho_t_distrib(lr(1), lr(2), uniquek(uk), id, state_handle, ens_size)
         rho3 = model_rho_t_distrib(ul(1), ul(2), uniquek(uk), id, state_handle, ens_size)
         rho4 = model_rho_t_distrib(ur(1), ur(2), uniquek(uk), id, state_handle, ens_size)

         do e = 1, ens_size
            if (k(e) == uniquek(uk) ) then
               fld(1, e) = dym*( dxm*rho1(e) + dx*rho2(e) ) + dy*( dxm*rho3(e) + dx*rho4(e) )
            endif
         enddo

         ! Interpolation for the Rho field at level k+1
         rho1 = model_rho_t_distrib(ll(1), ll(2), uniquek(uk)+1, id, state_handle, ens_size)
         rho2 = model_rho_t_distrib(lr(1), lr(2), uniquek(uk)+1, id, state_handle, ens_size)
         rho3 = model_rho_t_distrib(ul(1), ul(2), uniquek(uk)+1, id, state_handle, ens_size)
         rho4 = model_rho_t_distrib(ur(1), ur(2), uniquek(uk)+1, id, state_handle, ens_size)

         do e = 1, ens_size
            if (k(e) == uniquek(uk) ) then
               fld(2, e) = dym*( dxm*rho1(e) + dx*rho2(e) ) + dy*( dxm*rho3(e) + dx*rho4(e) )
            endif
         enddo

      endif

      enddo

   !-----------------------------------------------------
   ! 1.e Vertical Wind (W)

   elseif ( obs_kind == QTY_VERTICAL_VELOCITY ) then

      ! Adjust zloc for staggered ZNW grid (or W-grid, as compared to ZNU or M-grid)
      zloc = zloc + 0.5_r8

      !>@todo what should you do with this?
      k = max(1,int(zloc))

      deallocate(uniquek)

      ! Re-find the unique k values
      ksort = sort(k)
   
      count = 1
      do e = 2, ens_size
          if ( ksort(e) /= ksort(e-1) ) count = count + 1
      enddo
   
      allocate(uniquek(count))
      uniquek(:) = -1 
    
      uk = 1
      do e = 1, ens_size
         if ( all(uniquek /= k(e)) ) then
            uniquek(uk) = k(e)
            uk = uk + 1
         endif
      enddo

      call simple_interp_distrib(fld, wrf, id, i, j, k, obs_kind, dxm, dx, dy, dym, uniquek, ens_size, state_handle )
      if (all(fld == missing_r8)) goto 200

    !-----------------------------------------------------
   ! 1.f Specific Humidity (SH, SH2)
   ! Look at me
   ! Convert water vapor mixing ratio to specific humidity:
   else if( obs_kind == QTY_SPECIFIC_HUMIDITY ) then

      ! This is for 3D specific humidity -- surface spec humidity later
      if(.not. surf_var) then

         ! First confirm that vapor mixing ratio is in the DART state vector
         if ( wrf%dom(id)%type_qv >= 0 ) then

            UNIQUEK_LOOP: do uk = 1, count ! for the different ks

               ! Check to make sure retrieved integer gridpoints are in valid range
               if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
                    boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) .and. &
                    boundsCheck( uniquek(uk), .false.,                id, dim=3, type=wrf%dom(id)%type_t ) ) then

                  call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc ) ! HK why is this type_t
                  if ( rc .ne. 0 ) &
                       print*, 'model_mod.f90 :: model_interpolate :: getCorners SH rc = ', rc

                  ! Interpolation for SH field at level k
                  ill = get_dart_vector_index(ll(1), ll(2), uniquek(uk), domain_id(id), wrf%dom(id)%type_qv)
                  iul = get_dart_vector_index(ul(1), ul(2), uniquek(uk), domain_id(id), wrf%dom(id)%type_qv)
                  ilr = get_dart_vector_index(lr(1), lr(2), uniquek(uk), domain_id(id), wrf%dom(id)%type_qv)
                  iur = get_dart_vector_index(ur(1), ur(2), uniquek(uk), domain_id(id), wrf%dom(id)%type_qv)

                  x_ill = get_state(ill, state_handle)
                  x_iul = get_state(iul, state_handle)
                  x_ilr = get_state(ilr, state_handle)
                  x_iur = get_state(iur, state_handle)

                  do e = 1, ens_size
                     if ( k(e) == uniquek(uk) ) then ! interpolate only if it is the correct k
                        a1 = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur ) ! think this can go outside the k loop
                        fld(1,e) = a1(e) /(1.0_r8 + a1(e))
                     endif
                  enddo

                  ! Interpolation for SH field at level k+1
                  ill = get_dart_vector_index(ll(1), ll(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_qv)
                  iul = get_dart_vector_index(ul(1), ul(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_qv)
                  ilr = get_dart_vector_index(lr(1), lr(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_qv)
                  iur = get_dart_vector_index(ur(1), ur(2), uniquek(uk)+1, domain_id(id), wrf%dom(id)%type_qv)

                  x_ill = get_state(ill, state_handle)
                  x_ilr = get_state(ilr, state_handle)
                  x_iul = get_state(iul, state_handle)
                  x_iur = get_state(iur, state_handle)

                  do e = 1, ens_size
                     if ( k(e) == uniquek(uk) ) then ! interpolate only if it is the correct
                        a1 = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )
                        fld(2,e) = a1(e) /(1.0_r8 + a1(e))
                     endif
                  enddo
              endif
            enddo UNIQUEK_LOOP

         endif

      ! This is for surface specific humidity (calculated from Q2)
      else
         
         ! confirm that field is in the DART state vector
         if ( wrf%dom(id)%type_q2 >= 0 ) then
            ! Check to make sure retrieved integer gridpoints are in valid range
            if ( ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
                   boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) ) &
                   .or. wrf%dom(id)%scm ) then
               
               call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
               if ( rc .ne. 0 ) &
                    print*, 'model_mod.f90 :: model_interpolate :: getCorners Q2 rc = ', rc

               ! Interpolation for the SH2 field
               ill = get_dart_vector_index(ll(1), ll(2), 1, domain_id(id), wrf%dom(id)%type_q2)
               iul = get_dart_vector_index(ul(1), ul(2), 1, domain_id(id), wrf%dom(id)%type_q2)
               ilr = get_dart_vector_index(lr(1), lr(2), 1, domain_id(id), wrf%dom(id)%type_q2)
               iur = get_dart_vector_index(ur(1), ur(2), 1, domain_id(id), wrf%dom(id)%type_q2)

               x_ill = get_state(ill, state_handle)
               x_iul = get_state(iul, state_handle)
               x_iur = get_state(iur, state_handle)
               x_ilr = get_state(ilr, state_handle)

               a1 = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )
               fld(1,:) = a1 /(1.0_r8 + a1)

            endif
         endif
      endif


   !-----------------------------------------------------
   ! 1.g Vapor Mixing Ratio (QV, Q2)  
   else if( obs_kind == QTY_VAPOR_MIXING_RATIO ) then

      ! This is for 3D vapor mixing ratio -- surface QV later
      if(.not. surf_var) then
         call simple_interp_distrib(fld, wrf, id, i, j, k, obs_kind, dxm, dx, dy, dym, uniquek, ens_size, state_handle )
         if (all(fld == missing_r8)) goto 200
      else ! This is for surface QV (Q2)
         ! Confirm that right field is in the DART state vector
         if ( wrf%dom(id)%type_q2 >= 0 ) then
            !HK I am not sure what the type should be
            call surface_interp_distrib(fld, wrf, id, i, j, obs_kind, wrf%dom(id)%type_q2, dxm, dx, dy, dym, ens_size, state_handle)
            if (all(fld == missing_r8)) goto 200
         endif
      endif

      ! Don't accept negative water vapor amounts (?)
     fld = max(0.0_r8, fld)

   else if (obs_kind == QTY_2M_SPECIFIC_HUMIDITY ) then
      ! FIXME: Q2 is actually a mixing ratio, not a specific humidity
      if ( wrf%dom(id)%type_q2 >= 0 ) then
         ! Check to make sure retrieved integer gridpoints are in valid range
         if ( ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
                boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) ) &
                .or. wrf%dom(id)%scm ) then
            
            call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
            if ( rc .ne. 0 ) &
                 print*, 'model_mod.f90 :: model_interpolate :: getCorners Q2 rc = ', rc

            ! Interpolation for the SH2 field
            ill = get_dart_vector_index(ll(1), ll(2), 1, domain_id(id), wrf%dom(id)%type_q2)
            iul = get_dart_vector_index(ul(1), ul(2), 1, domain_id(id), wrf%dom(id)%type_q2)
            ilr = get_dart_vector_index(lr(1), lr(2), 1, domain_id(id), wrf%dom(id)%type_q2)
            iur = get_dart_vector_index(ur(1), ur(2), 1, domain_id(id), wrf%dom(id)%type_q2)

            x_ill = get_state(ill, state_handle)
            x_iul = get_state(iul, state_handle)
            x_iur = get_state(iur, state_handle)
            x_ilr = get_state(ilr, state_handle)

            a1 = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )
            fld(1,:) = a1
         endif
      endif

   !-----------------------------------------------------
   ! 1.t Pressure (P)
   else if( obs_kind == QTY_PRESSURE .or. obs_kind == QTY_SURFACE_PRESSURE ) then
            ! This is for the 3D pressure field -- surface pressure later
      if(.not. surf_var) then

         do uk = 1, count

         ! Check to make sure retrieved integer gridpoints are in valid range
         if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
              boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) .and. &
              boundsCheck( uniquek(uk), .false.,                id, dim=3, type=wrf%dom(id)%type_t ) ) then
   
            call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
            if ( rc .ne. 0 ) &
                 print*, 'model_mod.f90 :: model_interpolate :: getCorners P rc = ', rc
         
            ! Hmmm, it does not appear that P is part of the DART state vector, so there
            !   is not a reference to wrf%dom(id)%dart_ind -- we'll have to go right from
            !   the corner indices

            ! Interpolation for the P field at level k
            pres1 = model_pressure_t_distrib(ll(1), ll(2), uniquek(uk), id, state_handle, ens_size)
            pres2 = model_pressure_t_distrib(lr(1), lr(2), uniquek(uk), id, state_handle, ens_size)
            pres3 = model_pressure_t_distrib(ul(1), ul(2), uniquek(uk), id, state_handle, ens_size)
            pres4 = model_pressure_t_distrib(ur(1), ur(2), uniquek(uk), id, state_handle, ens_size)

            do e = 1, ens_size
               if ( k(e) == uniquek(uk) ) then ! interpolate only if it is the correct k
                  fld(1, e) = dym*( dxm*pres1(e) + dx*pres2(e) ) + dy*( dxm*pres3(e) + dx*pres4(e) )
               endif
            enddo

   
            ! Interpolation for the P field at level k+1
            pres1 = model_pressure_t_distrib(ll(1), ll(2), uniquek(uk)+1, id, state_handle, ens_size)
            pres2 = model_pressure_t_distrib(lr(1), lr(2), uniquek(uk)+1, id, state_handle, ens_size)
            pres3 = model_pressure_t_distrib(ul(1), ul(2), uniquek(uk)+1, id, state_handle, ens_size)
            pres4 = model_pressure_t_distrib(ur(1), ur(2), uniquek(uk)+1, id, state_handle, ens_size)

            do e = 1, ens_size
               if ( k(e) == uniquek(uk) ) then ! interpolate only if it is the correct k
                  fld(2, e) = dym*( dxm*pres1(e) + dx*pres2(e) ) + dy*( dxm*pres3(e) + dx*pres4(e) )
               endif
            enddo
   
         endif
         enddo

      !  This is for surface pressure (PSFC)
      else

         if ( wrf%dom(id)%type_ps >= 0 ) then

            ! Check to make sure retrieved integer gridpoints are in valid range
            if ( ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
                   boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) ) &
                   .or. wrf%dom(id)%scm ) then
      
               call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
               if ( rc .ne. 0 ) &
                    print*, 'model_mod.f90 :: model_interpolate :: getCorners PS rc = ', rc
      
               ! Interpolation for the PS field
               ill = get_dart_vector_index(ll(1), ll(2), 1, domain_id(id), wrf%dom(id)%type_ps)
               iul = get_dart_vector_index(ul(1), ul(2), 1, domain_id(id), wrf%dom(id)%type_ps)
               ilr = get_dart_vector_index(lr(1), lr(2), 1, domain_id(id), wrf%dom(id)%type_ps)
               iur = get_dart_vector_index(ur(1), ur(2), 1, domain_id(id), wrf%dom(id)%type_ps)

               x_ill = get_state(ill, state_handle)
               x_iul = get_state(iul, state_handle)
               x_ilr = get_state(ilr, state_handle)
               x_iur = get_state(iur, state_handle)

               do e = 1, ens_size
                  ! I'm not quite sure where this comes from, but I will trust them on it....
                  if ( x_ill(e) /= 0.0_r8 .and. x_ilr(e) /= 0.0_r8 .and. x_iul(e) /= 0.0_r8 .and. &
                       x_iur(e) /= 0.0_r8 ) then
      
                  fld(1, e) = dym*( dxm*x_ill(e) + dx*x_ilr(e) ) + dy*( dxm*x_iul(e) + dx*x_iur(e) )
      
                  ! JPH special treatment for scm configuration, where PS is not defined
                  ! on the boundaries and the weights are already 1,0
                  elseif ( wrf%dom(id)%scm ) then !HK is this redunant? What happens if you fail both conditions.
                     fld(1, e) = x_ill(e)
                  endif
               enddo
            endif
         endif
      endif


   !-----------------------------------------------------
   ! 1.u Vortex Center Stuff from Yongsheng
   else if ( obs_kind == QTY_VORTEX_LAT  .or. obs_kind == QTY_VORTEX_LON .or. &
             obs_kind == QTY_VORTEX_PMIN .or. obs_kind == QTY_VORTEX_WMAX ) then

      do uk = 1, count ! for the different ks

      ! Check to make sure retrieved integer gridpoints are in valid range
      if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
           boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) .and. &
           boundsCheck( uniquek(uk), .false.,                id, dim=3, type=wrf%dom(id)%type_t ) ) then

         !!   define a search box bounded by center_track_***
         center_search_half_size = nint(center_search_half_length/wrf%dom(id)%dx)
         if ( wrf%dom(id)%periodic_x ) then
            center_track_xmin       = i-center_search_half_size
            center_track_xmax       = i+center_search_half_size
         else
            center_track_xmin = max(1,i-center_search_half_size)
            center_track_xmax = min(wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu),i+center_search_half_size)
         endif
         if ( wrf%dom(id)%periodic_y ) then
            center_track_ymin       = j-center_search_half_size
            center_track_ymax       = j+center_search_half_size
         else
            center_track_ymin = max(1,j-center_search_half_size)
            center_track_ymax = min(wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu),j+center_search_half_size)
         endif
         if ( center_track_xmin<1 .or. center_track_xmax>wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu) .or. &
              center_track_ymin<1 .or. center_track_ymax>wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu) .or. &
              center_track_xmin >= center_track_xmax .or. center_track_ymin >= center_track_ymax) then

            print*,'i,j,center_search_half_length,center_track_xmin(max),center_track_ymin(max)'
            print*,i,j,center_search_half_length,center_track_xmin,center_track_xmax,center_track_ymin,center_track_ymax
            write(errstring,*)'Wrong setup in center_track_nml'
            call error_handler(E_ERR,'model_interpolate', errstring, source, revision, revdate)

         endif 

         if ( obs_kind == QTY_VORTEX_LAT .or. obs_kind == QTY_VORTEX_LON .or. &
              obs_kind == QTY_VORTEX_PMIN ) then

            !!   define spline interpolation box dimensions
            xlen = center_track_xmax - center_track_xmin + 1
            ylen = center_track_ymax - center_track_ymin + 1
            xxlen = (center_track_xmax - center_track_xmin)*center_spline_grid_scale + 1
            yylen = (center_track_ymax - center_track_ymin)*center_spline_grid_scale + 1

            allocate(x1d(xlen), y1d(ylen))  ;  allocate(xx1d(xxlen), yy1d(yylen))
            allocate(pd(xlen,ylen))         ;  allocate(pp(xxlen,yylen))
            allocate(vfld(xlen,ylen, ens_size))

            do i1 = 1,xlen
               x1d(i1) = (i1-1)+center_track_xmin
            enddo
            do i2 = 1,ylen
               y1d(i2) = (i2-1)+center_track_ymin
            enddo
            do ii1 = 1,xxlen
               xx1d(ii1) = center_track_xmin+real(ii1-1,r8)*1.0_r8/real(center_spline_grid_scale,r8)
            enddo
            do ii2 = 1,yylen
               yy1d(ii2) = center_track_ymin+real(ii2-1,r8)*1.0_r8/real(center_spline_grid_scale,r8)
            enddo

         endif

         if ( (obs_kind == QTY_VORTEX_LAT .or. obs_kind == QTY_VORTEX_LON) .and. (.not. use_old_vortex) ) then

            !  determine window that one would need wind components, thus circulation
            circ_half_size   = nint(circulation_radius/wrf%dom(id)%dx)
            circ_half_length = circulation_radius / wrf%dom(id)%dx

            if ( wrf%dom(id)%periodic_x ) then
               circ_xmin = i-center_search_half_size-circ_half_size
               circ_xmax = i+center_search_half_size+circ_half_size
            else
               circ_xmin = max(1,i-center_search_half_size-circ_half_size)
               circ_xmax = min(wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu),i+center_search_half_size+circ_half_size)
            endif

            if ( wrf%dom(id)%periodic_y ) then
               circ_ymin = j-center_search_half_size-circ_half_size
               circ_ymax = j+center_search_half_size+circ_half_size
            else
               circ_ymin = max(1,j-center_search_half_size-circ_half_size)
               circ_ymax = min(wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu),j+center_search_half_size+circ_half_size)
            endif

            cxlen = circ_xmax-circ_xmin+1
            cylen = circ_ymax-circ_ymin+1
            allocate(uwnd(cxlen+2,cylen+2, ens_size))
            allocate(vwnd(cxlen+2,cylen+2, ens_size))
            allocate(z1d(0:wrf%dom(id)%bt, ens_size))

            do i1 = circ_xmin-1, circ_xmax+1

               ii1 = i1
               if ( wrf%dom(id)%periodic_x ) then
                  if ( i1 > wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu) ) then
                     ii1 = i1 - wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu)
                  elseif ( i1 < 1 ) then
                     ii1 = i1 + wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu)
                  endif
               else
                  if ( i1 > wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu) ) then
                     ii1 = wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu)
                  elseif ( i1 < 1 ) then
                     ii1 = 1
                  endif
               endif

               do i2 = circ_ymin-1, circ_ymax+1

                  ii2 = i2
                  if ( wrf%dom(id)%periodic_y ) then
                     if ( i2 > wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu) ) then
                        ii2 = i2 - wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu)
                     elseif ( i2 < 1 ) then
                        ii2 = i2 + wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu)
                     endif
                  else
                     if ( i2 > wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu) ) then
                        ii2 = i2
                     elseif ( i2 < 1 ) then
                        ii2 = 1
                     endif
                  endif

                  !  calculate the wind components at the desired pressure level
                  do k2 = 1, wrf%dom(id)%var_size(3,wrf%dom(id)%type_t)
                     !z1d(k2) = model_pressure_t(i1,i2,k2,id,x)
                     z1d(k2, :) = model_pressure_t_distrib(i1,i2,k2,id, state_handle, ens_size)
                  enddo
                  z1d(0, :) = z1d(1, :)
                  !call pres_to_zk(circulation_pres_level, z1d, wrf%dom(id)%bt, zloc, is_lev0)
                  call pres_to_zk_distrib(circulation_pres_level, z1d, wrf%dom(id)%bt, ens_size, zloc, is_lev0)

                  do e = 1, ens_size

                     k2 = floor(zloc(e))  ;  dxm = mod(zloc(e),1.0_r8)  ;  dx = 1.0_r8 - dxm

                     if ( zloc(e) >= 1.0_r8 ) then  !  vertically interpolate

!                     ugrid = (dx  * (x(wrf%dom(id)%dart_ind(ii1,  ii2,  k2,  wrf%dom(id)%type_u))  + &
!                                     x(wrf%dom(id)%dart_ind(ii1+1,ii2,  k2,  wrf%dom(id)%type_u))) + &
!                              dxm * (x(wrf%dom(id)%dart_ind(ii1,  ii2,  k2+1,wrf%dom(id)%type_u))  + &
!                                     x(wrf%dom(id)%dart_ind(ii1+1,ii2,  k2+1,wrf%dom(id)%type_u)))) * 0.5_r8

                        ugrid_1 = get_dart_vector_index(ii1,  ii2,  k2,   domain_id(id), wrf%dom(id)%type_u)
                        ugrid_2 = get_dart_vector_index(ii1+1,ii2,  k2,   domain_id(id), wrf%dom(id)%type_u)
                        ugrid_3 = get_dart_vector_index(ii1,  ii2,  k2,   domain_id(id), wrf%dom(id)%type_u)
                        ugrid_4 = get_dart_vector_index(ii1,  ii2,  k2+1, domain_id(id), wrf%dom(id)%type_u )

                        x_ugrid_1 = get_state(ugrid_1, state_handle)
                        x_ugrid_2 = get_state(ugrid_2, state_handle)
                        x_ugrid_3 = get_state(ugrid_3, state_handle)
                        x_ugrid_4 = get_state(ugrid_4, state_handle)

                        ugrid = (dx  * (x_ugrid_1  + x_ugrid_2) + dxm * (x_ugrid_3 + x_ugrid_4)) * 0.5_r8

                        vgrid_1 = get_dart_vector_index(ii1,  ii2,  k2,   domain_id(id),wrf%dom(id)%type_v)
                        vgrid_2 = get_dart_vector_index(ii1,  ii2+1,k2,   domain_id(id),wrf%dom(id)%type_v)
                        vgrid_3 = get_dart_vector_index(ii1,  ii2,  k2+1, domain_id(id),wrf%dom(id)%type_v)
                        vgrid_4 = get_dart_vector_index(ii1,  ii2+1,k2+1, domain_id(id),wrf%dom(id)%type_v)

                        x_vgrid_1 = get_state(vgrid_1, state_handle)
                        x_vgrid_2 = get_state(vgrid_2, state_handle)
                        x_vgrid_3 = get_state(vgrid_3, state_handle)
                        x_vgrid_4 = get_state(vgrid_4, state_handle)

!                     vgrid = (dx  * (x(wrf%dom(id)%dart_ind(ii1,  ii2,  k2,  wrf%dom(id)%type_v))  + &
!                                     x(wrf%dom(id)%dart_ind(ii1,  ii2+1,k2,  wrf%dom(id)%type_v))) + &
!                              dxm * (x(wrf%dom(id)%dart_ind(ii1,  ii2,  k2+1,wrf%dom(id)%type_v))  + &
!                                     x(wrf%dom(id)%dart_ind(ii1,  ii2+1,k2+1,wrf%dom(id)%type_v)))) * 0.5_r8

                        vgrid = (dx  * (vgrid_1 + vgrid_2) + dxm * (vgrid_3 + vgrid_4)) * 0.5_r8

                     else  !  pressure level below ground.  Take model level 1 winds

!                     ugrid = (x(wrf%dom(id)%dart_ind(ii1,  ii2,  1,wrf%dom(id)%type_u)) + &
!                              x(wrf%dom(id)%dart_ind(ii1+1,ii2,  1,wrf%dom(id)%type_u))) * 0.5_r8

                        ugrid_1 = get_dart_vector_index(ii1,  ii2,  1, domain_id(id),wrf%dom(id)%type_u)
                        ugrid_2 = get_dart_vector_index(ii1+1,ii2,  1, domain_id(id),wrf%dom(id)%type_u)

                        x_ugrid_1 = get_state(ugrid_1, state_handle)
                        x_ugrid_2 = get_state(ugrid_2, state_handle)
 
                        ugrid = (x_ugrid_1 + x_ugrid_2) * 0.5_r8

!                     vgrid = (x(wrf%dom(id)%dart_ind(ii1,  ii2,  1,wrf%dom(id)%type_v)) + &
!                              x(wrf%dom(id)%dart_ind(ii1,  ii2+1,1,wrf%dom(id)%type_v))) * 0.5_r8

                        vgrid_1 = get_dart_vector_index(ii1,  ii2,  1, domain_id(id),wrf%dom(id)%type_v)
                        vgrid_2 = get_dart_vector_index(ii1,  ii2+1,1, domain_id(id),wrf%dom(id)%type_v)

                        x_vgrid_1 = get_state(vgrid_1, state_handle)
                        x_vgrid_2 = get_state(vgrid_2, state_handle)

                        vgrid = (x_vgrid_1 + x_vgrid_2) * 0.5_r8

                     endif
                     uwnd(i1-circ_xmin+2,i2-circ_ymin+2, :) = ugrid
                     vwnd(i1-circ_xmin+2,i2-circ_ymin+2, :) = vgrid

                  enddo

               enddo
            enddo

            allocate(vort(cxlen,cylen, ens_size))
            do i1 = circ_xmin, circ_xmax

               dgi1 = 2.0_r8
               if ( .not. wrf%dom(id)%periodic_x ) then
                  if     ( i1 == 1 ) then
                    dgi1 = 1.0_r8
                  elseif ( i1 == wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu) ) then
                    dgi1 = 1.0_r8
                  endif
               endif

               do i2 = circ_ymin, circ_ymax

                  dgi2 = 2.0_r8
                  if ( .not. wrf%dom(id)%periodic_x ) then
                     if     ( i2 == 1 ) then
                       dgi2 = 1.0_r8
                     elseif ( i2 == wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu) ) then
                       dgi2 = 1.0_r8
                     endif
                  endif

                  ii1  = i1-circ_xmin+1
                  ii2  = i2-circ_ymin+1

                  !  compute the vorticity for each point needed to compute circulation
                  vort(ii1,ii2, :) = (vwnd(ii1+2,ii2+1, :) - vwnd(ii1+1,ii2+1, :)) / (wrf%dom(id)%dx * dgi2) + &
                                  (vwnd(ii1+1,ii2+1, :) - vwnd(ii1  ,ii2+1, :)) / (wrf%dom(id)%dx * dgi2) - &
                                  (uwnd(ii1+1,ii2+2, :) - uwnd(ii1+1,ii2+1, :)) / (wrf%dom(id)%dx * dgi1) - &
                                  (uwnd(ii1+1,ii2+1, :) - uwnd(ii1+1,ii2, :  )) / (wrf%dom(id)%dx * dgi1)

               enddo
            enddo

            !  loop over all grid points in search area, compute average vorticity
            !  (circulation) within a certain distance of that grid point
            do i1 = center_track_xmin, center_track_xmax
            do i2 = center_track_ymin, center_track_ymax

               asum = 0.0_r8  ;  circ = 0.0_r8
               do ii1 = max(circ_xmin,i1-circ_half_size), min(circ_xmax,i1+circ_half_size)
               do ii2 = max(circ_ymin,i2-circ_half_size), min(circ_ymax,i2+circ_half_size)

                  distgrid = sqrt(real(ii1-i1,r8) ** 2 + real(ii2-i2,r8) ** 2) * wrf%dom(id)%dx
                  if ( distgrid <= circulation_radius ) then

                     asum = asum + 1.0_r8
                     circ = circ + vort(ii1-circ_xmin+1,ii2-circ_ymin+1, :)

                  endif

               enddo
               enddo

               vfld(i1-center_track_xmin+1,i2-center_track_ymin+1, :) = circ / asum

            enddo
            enddo

            !  find maximum in circulation through spline interpolation
            do e = 1, ens_size
               call splie2(x1d,y1d,vfld(:,:,e),xlen,ylen,pd)
            enddo

            vcrit = -1.0e20_r8
            cxloc = -1
            cyloc = -1
            do ii1 = 1, xxlen
            do ii2 = 1, yylen
               do e = 1, ens_size
                  call splin2(x1d,y1d,vfld(:,:,e),pd,xlen,ylen,xx1d(ii1),yy1d(ii2),pp(ii1,ii2))
               enddo
               if (vcrit < pp(ii1,ii2)) then
                  vcrit = pp(ii1,ii2)
                  cxloc = xx1d(ii1)
                  cyloc = yy1d(ii2)
               endif
            enddo
            enddo

            !  forward operator fails if maximum is at edge of search area
            if ( cxloc-xx1d(1) < 1.0_r8 .or. xx1d(xxlen)-cxloc < 1.0_r8 .or. &
                 cyloc-yy1d(1) < 1.0_r8 .or. yy1d(yylen)-cyloc < 1.0_r8 ) then

               do e = 1, ens_size
                  if ( k(e) == uniquek(uk) ) then ! interpolate only if is the correct k
                     fld(:, e) = missing_r8
                  endif
               enddo

            else

            call ij_to_latlon(wrf%dom(id)%proj, cxloc, cyloc, clat, clon)

            if ( obs_kind == QTY_VORTEX_LAT ) then
               do e = 1, ens_size
                  if ( k(e) == uniquek(uk) ) then ! interpolate only if is the correct k
                     fld(1, e) = clat
                  endif
               enddo

            else
               do e = 1, ens_size
                  if ( k(e) == uniquek(uk) ) then ! interpolate only if is the correct k
                     fld(1, e) = clon
                  endif
               enddo

            endif

            endif

            deallocate(uwnd, vwnd, vort, z1d)
            deallocate(vfld, pd, pp, x1d, y1d, xx1d, yy1d)

         else if ( obs_kind == QTY_VORTEX_PMIN .or. (use_old_vortex .and. & 
                  (obs_kind == QTY_VORTEX_LAT .or. obs_kind == QTY_VORTEX_LON)) ) then

            allocate(p1d(wrf%dom(id)%bt, ens_size),  t1d(wrf%dom(id)%bt, ens_size))
            allocate(qv1d(wrf%dom(id)%bt, ens_size), z1d(wrf%dom(id)%bt, ens_size))
            allocate(z1d_1(wrf%dom(id)%bt, ens_size), z1d_2(wrf%dom(id)%bt, ens_size))

            !  compute SLP for each grid point within the search area
            print*, center_track_xmin, center_track_xmax
            print*, center_track_ymin, center_track_ymax
            do i1 = center_track_xmin, center_track_xmax

               ii1 = i1
               if ( wrf%dom(id)%periodic_x ) then
                  if ( i1 > wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu) ) then
                     ii1 = i1 - wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu)
                  elseif ( i1 < 1 ) then
                     ii1 = i1 + wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu)
                  endif
               endif

               do i2 = center_track_ymin, center_track_ymax

                  ii2 = i2
                  if ( wrf%dom(id)%periodic_y ) then
                     if ( i2 > wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu) ) then
                        ii2 = i2 - wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu)
                     elseif ( i2 < 1 ) then
                        ii2 = i2 + wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu)
                     endif
                  endif

                  do k2 = 1,wrf%dom(id)%var_size(3,wrf%dom(id)%type_t)

!                     p1d(k2) = model_pressure_t(ii1,ii2,k2,id,x)
                      p1d(k2, :) = model_pressure_t_distrib(ii1,ii2,k2, id, state_handle, ens_size)
                      !print*, 'p1d(k2, 1)', p1d(k2, 1)

!                     t1d(k2) = x(wrf%dom(id)%dart_ind(ii1,ii2,k2,wrf%dom(id)%type_t)) + ts0
                     t1d_ind = get_dart_vector_index(ii1,ii2,k2, domain_id(id),wrf%dom(id)%type_t)
                     t1d(k2, :) = get_state( t1d_ind, state_handle)
                     t1d(k2, :) = t1d(k2, :) + ts0
                     !print*, 't1d(k2, 1)', t1d(k2, 1)

!                     qv1d(k2)= x(wrf%dom(id)%dart_ind(ii1,ii2,k2,wrf%dom(id)%type_qv))
                     qv1d_ind = get_dart_vector_index(ii1,ii2,k2, domain_id(id),wrf%dom(id)%type_qv)
                     qv1d(k2, :) = get_state(qv1d_ind, state_handle)
                     !print*, 'qv1d(k2, 1)', qv1d(k2, 1)

!                     z1d(k2) = (x(wrf%dom(id)%dart_ind(ii1,ii2,k2,  wrf%dom(id)%type_gz))+ &
!                                x(wrf%dom(id)%dart_ind(ii1,ii2,k2+1,wrf%dom(id)%type_gz))+ &
!                                wrf%dom(id)%phb(ii1,ii2,k2)+wrf%dom(id)%phb(ii1,ii2,k2+1))*0.5_r8/gravity
                     z1d_ind1 = get_dart_vector_index(ii1,ii2,k2,   domain_id(id),wrf%dom(id)%type_gz)
                     z1d_ind2 = get_dart_vector_index(ii1,ii2,k2+1, domain_id(id),wrf%dom(id)%type_gz)

                     z1d_1(k2, :) = get_state(z1d_ind1, state_handle)
                     z1d_2(k2, :) = get_state(z1d_ind2, state_handle)

                     z1d(k2, :) = ( z1d_1(k2, :)+ z1d_2(k2, :) + &
                                    wrf%dom(id)%phb(ii1,ii2,k2)+wrf%dom(id)%phb(ii1,ii2,k2+1))*0.5_r8/gravity

                     !print*, 'z1d(k2, 1)', z1d(k2, 1)


                  enddo
                  do e = 1, ens_size
                     call compute_seaprs(wrf%dom(id)%bt, z1d(:,e), t1d(:,e), p1d(:,e), qv1d(:,e), &
                                      vfld(i1-center_track_xmin+1,i2-center_track_ymin+1, e), debug)
                  enddo
               enddo

            enddo

            !print*, 'vfld(1, 1)', vfld(1, 1, 1)
            do e = 1, ens_size

               if ( k(e) == uniquek(uk) ) then

                  !  find minimum in MSLP through spline interpolation
                  call splie2(x1d,y1d,vfld(:,:,e),xlen,ylen,pd)

                     vcrit = 1.0e20_r8
                     do ii1=1,xxlen
                     do ii2=1,yylen
                        call splin2(x1d,y1d,vfld(:,:,e),pd,xlen,ylen,xx1d(ii1),yy1d(ii2),pp(ii1,ii2))
                        if ( vcrit > pp(ii1,ii2)) then
                           vcrit = pp(ii1,ii2)
                           cxloc = xx1d(ii1)
                           cyloc = yy1d(ii2)
                        endif
                     enddo
                     enddo

                     !  forward operator fails if maximum is at edge of search area
                     if ( cxloc-xx1d(1) < 1.0_r8 .or. xx1d(xxlen)-cxloc < 1.0_r8 .or. &
                        cyloc-yy1d(1) < 1.0_r8 .or. yy1d(yylen)-cyloc < 1.0_r8 ) then

                        fld(:, e) = missing_r8

                     else

                        call ij_to_latlon(wrf%dom(id)%proj, cxloc, cyloc, clat, clon)

                        if ( obs_kind == QTY_VORTEX_PMIN ) then
                           fld(1, e) = vcrit
                        else if ( obs_kind == QTY_VORTEX_LAT ) then
                           fld(1, e) = clat
                        else
                           fld(1, e) = clon
                        endif

                     endif

               endif

            enddo


            print*, 'fld', fld
            deallocate(p1d, t1d, qv1d, z1d)
            deallocate(vfld, pd, pp, x1d, y1d, xx1d, yy1d)
            if (all(fld == missing_r8)) goto 200

         else if ( obs_kind == QTY_VORTEX_WMAX ) then   !  Maximum wind speed

            maxwspd = 0.0_r8
            do i1 = center_track_xmin, center_track_xmax

               ii1 = i1
               if ( wrf%dom(id)%periodic_x ) then
                  if ( i1 > wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu) ) then
                     ii1 = i1 - wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu)
                  elseif ( i1 < 1 ) then
                     ii1 = i1 + wrf%dom(id)%var_size(1,wrf%dom(id)%type_mu)
                  endif
               endif

               do i2 = center_track_ymin, center_track_ymax

                  ii2 = i2
                  if ( wrf%dom(id)%periodic_y ) then
                     if ( i2 > wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu) ) then
                        ii2 = i2 - wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu)
                     elseif ( i2 < 1 ) then
                        ii2 = i2 + wrf%dom(id)%var_size(2,wrf%dom(id)%type_mu)
                     endif
                  endif

                  if ( ( wrf%dom(id)%type_u10 >= 0 ) .and. ( wrf%dom(id)%type_v10 >= 0 ) ) then
                     !ugrid = x(wrf%dom(id)%dart_ind(ii1,ii2,1,wrf%dom(id)%type_u10))
                     !vgrid = x(wrf%dom(id)%dart_ind(ii1,ii2,1,wrf%dom(id)%type_v10))

                     ugrid_1 = get_dart_vector_index(ii1,ii2,1, domain_id(id),wrf%dom(id)%type_u10)
                     ugrid = get_state(ugrid_1, state_handle)
 
                     vgrid_1 = get_dart_vector_index(ii1,ii2,1, domain_id(id),wrf%dom(id)%type_v10)
                     vgrid = get_state(vgrid_1, state_handle)

                  else

! Same code as above
!                     ugrid = 0.5_r8*(x(wrf%dom(id)%dart_ind(ii1,  ii2,1,wrf%dom(id)%type_u)) + &
!                                     x(wrf%dom(id)%dart_ind(ii1+1,ii2,1,wrf%dom(id)%type_u)))
!                     vgrid = 0.5_r8*(x(wrf%dom(id)%dart_ind(ii1,ii2,  1,wrf%dom(id)%type_v)) + &
!                                     x(wrf%dom(id)%dart_ind(ii1,ii2+1,1,wrf%dom(id)%type_v)))
                     ugrid_1 = get_dart_vector_index(ii1,  ii2,  1, domain_id(id),wrf%dom(id)%type_u)
                     ugrid_2 = get_dart_vector_index(ii1+1,ii2,  1, domain_id(id),wrf%dom(id)%type_u)

                     x_ugrid_1 = get_state(ugrid_1, state_handle)
                     x_ugrid_2 = get_state(ugrid_2, state_handle)
 
                     ugrid = (x_ugrid_1 + x_ugrid_2) * 0.5_r8

                     vgrid_1 = get_dart_vector_index(ii1,  ii2,  1, domain_id(id),wrf%dom(id)%type_v)
                     vgrid_2 = get_dart_vector_index(ii1,  ii2+1,1, domain_id(id),wrf%dom(id)%type_v)

                     x_vgrid_1 = get_state(vgrid_1, state_handle)
                     x_vgrid_2 = get_state(vgrid_2, state_handle)

                     vgrid = (x_vgrid_1 + x_vgrid_2) * 0.5_r8

                  endif

                  magwnd  = sqrt(ugrid * ugrid + vgrid * vgrid)

                  do e = 1, ens_size
                     if ( k(e) == uniquek(uk) ) then ! interpolate only if is the correct k
                        if ( magwnd(e) > maxwspd(e) ) then
                           imax    = i1
                           jmax    = i2
                           maxwspd(e) = magwnd(e)
                        endif
                     endif
                  enddo



               enddo
            enddo

            !  forward operator fails if maximum is at edge of search area
            if ( imax == center_track_xmin .or. jmax == center_track_ymin .or. &
                 imax == center_track_xmax .or. jmax == center_track_ymax ) then
               do e = 1, ens_size
                  if ( k(e) == uniquek(uk) ) then ! interpolate only if is the correct k
                     fld(:, e) = missing_r8
                  endif
               enddo
            else
               do e = 1, ens_size
                  if ( k(e) == uniquek(uk) ) then ! interpolate only if is the correct k
                     fld(1, e) = maxwspd(e)
                  endif
               enddo
            endif

         endif  !  if test on obs_kind

      endif   ! bounds check failed.

      enddo 

!*****************************************************************************
! END OF VERBATIM BIT - what does this mean?
!*****************************************************************************


   !-----------------------------------------------------
   ! 1.w Geopotential Height (GZ)

   !   GZ is on the ZNW grid (bottom_top_stagger), so its bottom-most level is defined to
   !   be at eta = 1 (the surface).  Thus, we have a 3D variable that contains a surface
   !   variable; the same is true for W as well.  If one wants to observe the surface value
   !   of either of these variables, then one can simply operate on the full 3D field 
   !   (toGrid below should return dz ~ 0 and dzm ~ 1) 
   else if( obs_kind == QTY_GEOPOTENTIAL_HEIGHT ) then
      !if( my_task_id() == 0 ) print*, '*** geopotential height forward operator not tested'

      ! make sure vector includes the needed field
      if ( wrf%dom(id)%type_gz >= 0 ) then

         ! Adjust zloc for staggered ZNW grid (or W-grid, as compared to ZNU or M-grid)
         zloc = zloc + 0.5_r8
         k = max(1,int(zloc))  ! Only 1 value of k across the ensemble?

         deallocate(uniquek)
         ! Re-find the unique k values
         ksort = sort(k)
      
         count = 1
         do e = 2, ens_size
             if ( ksort(e) /= ksort(e-1) ) count = count + 1
         enddo
      
         allocate(uniquek(count))
         uniquek(:) = -1 
       
         uk = 1
         do e = 1, ens_size
            if ( all(uniquek /= k(e)) ) then
               uniquek(uk) = k(e)
               uk = uk + 1
            endif
         enddo

         ! Check to make sure retrieved integer gridpoints are in valid range
         if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_gz ) .and. &
              boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_gz ) .and. &
              boundsCheck( k(1), .false.,                id, dim=3, type=wrf%dom(id)%type_gz ) ) then
            
            call getCorners(i, j, id, wrf%dom(id)%type_gz, ll, ul, lr, ur, rc )
            if ( rc .ne. 0 ) &
                 print*, 'model_mod.f90 :: model_interpolate :: getCorners GZ rc = ', rc
            
            ! Interpolation for GZ field at level k
            ill = get_dart_vector_index(ll(1), ll(2), k(1), domain_id(id), wrf%dom(id)%type_gz)
            iul = get_dart_vector_index(ul(1), ul(2), k(1), domain_id(id), wrf%dom(id)%type_gz)
            ilr = get_dart_vector_index(lr(1), lr(2), k(1), domain_id(id), wrf%dom(id)%type_gz)
            iur = get_dart_vector_index(ur(1), ur(2), k(1), domain_id(id), wrf%dom(id)%type_gz)

            x_ill = get_state(ill, state_handle)
            x_iul = get_state(iul, state_handle)
            x_iur = get_state(iur, state_handle)
            x_ilr = get_state(ilr, state_handle)

            fld(1,:) = ( dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur ) + &
                       dym*( dxm*wrf%dom(id)%phb(ll(1), ll(2), k)   + &
                             dx *wrf%dom(id)%phb(lr(1), lr(2), k) ) + &
                       dy *( dxm*wrf%dom(id)%phb(ul(1), ul(2), k)   + &
                             dx *wrf%dom(id)%phb(ur(1), ur(2), k) ) )  / gravity
            
            ! Interpolation for GZ field at level k+1
            ill = get_dart_vector_index(ll(1), ll(2), k(1)+1, domain_id(id), wrf%dom(id)%type_gz)
            iul = get_dart_vector_index(ul(1), ul(2), k(1)+1, domain_id(id), wrf%dom(id)%type_gz)
            ilr = get_dart_vector_index(lr(1), lr(2), k(1)+1, domain_id(id), wrf%dom(id)%type_gz)
            iur = get_dart_vector_index(ur(1), ur(2), k(1)+1, domain_id(id), wrf%dom(id)%type_gz)

            x_ill = get_state(ill, state_handle)
            x_iul = get_state(iul, state_handle)
            x_iur = get_state(iur, state_handle)
            x_ilr = get_state(ilr, state_handle)

            fld(2, :) = ( dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur ) + &
                       dym*( dxm*wrf%dom(id)%phb(ll(1), ll(2), k(1)+1)   + &
                             dx *wrf%dom(id)%phb(lr(1), lr(2), k(1)+1) ) + &
                       dy *( dxm*wrf%dom(id)%phb(ul(1), ul(2), k(1)+1)   + &
                             dx *wrf%dom(id)%phb(ur(1), ur(2), k(1)+1) ) )  / gravity
   
         endif
      endif

     !-----------------------------------------------------
   ! 1.x Surface Elevation (HGT)

   ! Surface Elevation has been added by Ryan Torn to accommodate altimeter observations.
   !   HGT is not in the dart_ind vector, so get it from wrf%dom(id)%hgt.
   else if( obs_kind == QTY_SURFACE_ELEVATION ) then

      if ( debug ) print*,'Getting surface elevation'

      ! Check to make sure retrieved integer gridpoints are in valid range
      if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
           boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) ) then

         call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
         if ( rc .ne. 0 ) &
              print*, 'model_mod.f90 :: model_interpolate :: getCorners HGT rc = ', rc
        
         ! Interpolation for the HGT field -- HGT is NOT part of state vector x, but rather
         !   in the associated domain meta data
         fld(1, :) = dym*( dxm*wrf%dom(id)%hgt(ll(1), ll(2)) + &
                         dx*wrf%dom(id)%hgt(lr(1), lr(2)) ) + &
                   dy*( dxm*wrf%dom(id)%hgt(ul(1), ul(2)) + &
                         dx*wrf%dom(id)%hgt(ur(1), ur(2)) )

      endif


   !-----------------------------------------------------
   ! 1.y Surface Skin Temperature (TSK)

   else if( obs_kind == QTY_SKIN_TEMPERATURE ) then
     ! make sure vector includes the needed field
     if ( wrf%dom(id)%type_tsk >= 0 ) then
        call surface_interp_distrib(fld, wrf, id, i, j, obs_kind, wrf%dom(id)%type_tsk, dxm, dx, dy, dym, ens_size, state_handle)
        if (all(fld == missing_r8)) goto 200
     endif

   !-----------------------------------------------------
   ! 1.z Land Mask (XLAND)

   ! Land Mask has been added to accommodate satellite observations.
   !   XLAND is not in the dart_ind vector, so get it from wrf%dom(id)%land
   else if( obs_kind == QTY_LANDMASK ) then
      if( my_task_id() == 0 ) print*, '*** Land mask forward operator not tested'

      if ( debug ) print*,'Getting land mask'

      ! Check to make sure retrieved integer gridpoints are in valid range
      if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
           boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) ) then
      
         call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
         if ( rc .ne. 0 ) &
              print*, 'model_mod.f90 :: model_interpolate :: getCorners XLAND rc = ', rc
         
         ! Interpolation for the XLAND field -- XLAND is NOT part of state vector x, but rather
         !   in the associated domain meta data
         fld(1, :) = dym*( dxm*real(wrf%dom(id)%land(ll(1), ll(2))) + &
                         dx*real(wrf%dom(id)%land(lr(1), lr(2))) ) + &
                   dy*( dxm*real(wrf%dom(id)%land(ul(1), ul(2))) + &
                         dx*real(wrf%dom(id)%land(ur(1), ur(2))) )

      endif

   else if( obs_kind == QTY_SURFACE_TYPE ) then ! 0 = land, 1 = water, 2 = sea ice (not yet supported)
      if ( debug ) print*,'Getting land mask'

      ! Check to make sure retrieved integer gridpoints are in valid range
      if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
           boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) ) then
      
         call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
         if ( rc .ne. 0 ) &
              print*, 'model_mod.f90 :: model_interpolate :: getCorners XLAND rc = ', rc
         
         ! Interpolation for the XLAND field -- XLAND is NOT part of state vector x, but rather
         !   in the associated domain meta data
         fld(1, :) = dym*( dxm*real(wrf%dom(id)%land(ll(1), ll(2))) + &
                         dx*real(wrf%dom(id)%land(lr(1), lr(2))) ) + &
                   dy*( dxm*real(wrf%dom(id)%land(ul(1), ul(2))) + &
                         dx*real(wrf%dom(id)%land(ur(1), ur(2))) ) - 1

      endif

   !-----------------------------------------------------
   ! If obs_kind is not negative (for identity obs), or if it is not one of the above 15
   !   explicitly checked-for kinds, then set error istatus and missing_r8.
   else

      expected_obs = missing_r8
      istatus = 3
      if (debug) print*, 'unrecognized obs KIND, value = ', obs_kind
      deallocate(v_h, v_p)
      if (allocated(uniquek)) deallocate(uniquek)
      return

   endif


   !----------------------------------
   ! 2. Vertical Interpolation 
   !----------------------------------

   ! Do vertical interpolation -- only for non-surface, non-identity obs.  

   ! The previous section (1. Horizontal Interpolation) has produced a variable called
   !   "fld", which nominally has two entries in it.  3D fields have hopefully produced
   !   2 non-zero entries, whereas surface fields only have filled the first entry.
   ! If a full 3D field, then do vertical interpolation between sandwiched model levels
   !   (k and k+1).

   ! Check to make sure that we did something sensible in the Horizontal Interpolation 
   !   section above.  All valid obs_kinds will have changed fld(1,e).

   !HK I am unsure as to whether this should be done on the array expected_obs or one ensemble
   ! member (e) at a time.
   do e = 1, ens_size

      if ( fld(1,e) == missing_r8 ) then

         expected_obs(e) = missing_r8
   
      else ! We purposefully changed fld(1,e), so continue onward

         ! If a surface variable, or a variable with no particular vertical location
         ! (basically the entire column) then no need to do any vertical interpolation
         if ( surf_var .or. is_vertical(location,"UNDEFINED") ) then

            !obs_val = fld(1)
             expected_obs(e) = fld(1,e) !HK

         ! If an interior variable, then we DO need to do vertical interpolation
         else

            ! First make sure fld(2,:) is no longer a missing value
            if ( fld(2,e) == missing_r8 ) then !HK should be any?
               print *,'fld 2 is missing',surf_var

               expected_obs(e) = missing_r8

            ! Do vertical interpolation -- at this point zloc is >= 1 unless
            ! the namelist value allow_obs_below_vol is true, in which case
            ! it is >= 0, and < 1 is a request to extrapolate.
            else

               ! Get fractional distances between grid points
               call toGrid(zloc(e), k(e), dz(e), dzm(e))
               if (debug) print*, 'zloc(e), k(e), dz(e), dzm(e) = ', zloc(e), k(e), dz(e), dzm(e)
               if (debug) print*, 'fld(1,e), fld(2,e) = ', fld(1,e), fld(2,e)

               ! If you get here and zloc < 1.0, then k will be 0, and
               ! we should extrapolate.  fld(1,:) and fld(2,:) where computed
               ! at levels 1 and 2.

               if (k(e) >= 1) then
                  ! Linearly interpolate between grid points
                  expected_obs(e) = dzm(e)*fld(1,e) + dz(e)*fld(2,e)
                  if (debug) print*, 'interpolated obs_val = ', expected_obs(e)
               else
                  ! Extrapolate below first level.
                  expected_obs(e) = fld(1,e) - (fld(2,e)-fld(1,e))*dzm(e)
                  if (debug) print*, 'extrapolated obs_val = ', expected_obs(e)
               endif

            endif
         endif
      endif

      ! Fill in failed copies
      if ( failedcopies(e) == missing_r8 ) expected_obs(e) = missing_r8

   enddo

endif  ! end of "if ( obs_kind < 0 )"

200  continue

if (all(fld == missing_r8)) then
   expected_obs(:) = missing_r8
   istatus(:) = 99
endif 

! Now that we are done, check to see if a missing value somehow 
! made it through without already having set an error return code.
do e = 1, ens_size
   if ( expected_obs(e) == missing_r8 .and. istatus(e) == 0 ) then
      istatus(e) = 99
   endif
enddo

! Pring the observed value if in debug mode
if(debug) then
  print*,'model_interpolate() return value for obs_kind ',obs_kind, ' = ',expected_obs
endif

! Deallocate variables before exiting
deallocate(v_h, v_p)
deallocate(uniquek)

end subroutine model_interpolate

!#######################################################################
subroutine convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, status)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer,             intent(in)    :: loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: status(:)

integer :: i

do i=1, num
   call vert_convert(state_handle, locs(i), loc_qtys(i), status(i))
enddo

end subroutine convert_vertical_obs


!#######################################################################
subroutine convert_vertical_state(state_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

integer :: i, istat
integer     :: var_type, dart_type
integer(i8) :: index
integer     :: ip, jp, kp
integer     :: nz, ny, nx
logical     :: var_found
real(r8)    :: lon, lat, lev
character(len=129) :: string1

integer :: id, var_id, state_id
logical, parameter :: debug = .false.


istatus = 0

do i=1, num

   ! from the dart index get the local variables indices
   call get_model_variable_indices(loc_indx(i), ip, jp, kp, var_id=var_id, dom_id=state_id)
   
   ! convert from state_structure domain number to wrf.
   id = get_wrf_domain(state_id)
   
   ! at this point, (ip,jp,kp) refer to indices in the variable's own grid
   
   if(debug) write(*,*) ' ip, jp, kp for index ',ip,jp,kp,index
   if(debug) write(*,*) ' Var type: ',var_type
   
   var_type  = wrf%dom(id)%var_type(var_id)
   dart_type = wrf%dom(id)%dart_kind(var_id)
   
   ! first obtain lat/lon from (ip,jp)
   call get_wrf_horizontal_location( ip, jp, var_type, id, lon, lat )
   
   ! now convert to desired vertical coordinate (defined in the namelist)
   if (wrf%dom(id)%localization_coord == VERTISLEVEL) then
      ! here we need level index of mass grid
      if( (var_type == wrf%dom(id)%type_w ) .or. (var_type == wrf%dom(id)%type_gz) ) then
         lev = real(kp) - 0.5_r8
      else
         lev = real(kp)
      endif
   elseif (wrf%dom(id)%localization_coord == VERTISPRESSURE) then
      ! directly convert to pressure
      lev = model_pressure_distrib(ip, jp, kp, id, var_type, state_handle)
   elseif (wrf%dom(id)%localization_coord == VERTISHEIGHT) then
      lev = model_height_distrib(ip, jp, kp, id, var_type, state_handle)
   elseif (wrf%dom(id)%localization_coord == VERTISSCALEHEIGHT) then
      lev = -log(model_pressure_distrib(ip, jp, kp, id, var_type, state_handle) / &
                 model_surface_pressure_distrib(ip, jp, id, var_type, state_handle))
   endif
   
   locs(i) = set_location(lon, lat, lev, wrf%dom(id)%localization_coord)
   
enddo

end subroutine convert_vertical_state


!#######################################################################
!> This is used in the filter_assim. The vertical conversion is done using the 
!> mean state.
!> I think at the moment you are over communicating
subroutine vert_convert(state_handle, location, obs_kind, istatus)

! This subroutine converts a given ob/state vertical coordinate to
! the vertical localization coordinate type requested through the 
! model_mod namelist.
!
! Notes: (1) obs_kind is only necessary to check whether the ob
!            is an identity ob.
!        (2) This subroutine can convert both obs' and state points'
!            vertical coordinates. Remember that state points get
!            their DART location information from get_state_meta_data
!            which is called by filter_assim during the assimilation
!            process.
!        (3) x is the relevant DART state vector for carrying out
!            interpolations necessary for the vertical coordinate
!            transformations. As the vertical coordinate is only used
!            in distance computations, this is actually the "expected"
!            vertical coordinate, so that computed distance is the
!            "expected" distance. Thus, under normal circumstances,
!            x that is supplied to vert_convert should be the
!            ensemble mean. Nevertheless, the subroutine has the
!            functionality to operate on any DART state vector that
!            is supplied to it.

type(ensemble_type),    intent(in)    :: state_handle
type(location_type),    intent(inout) :: location
integer,                intent(in)    :: obs_kind
integer,                intent(out)   :: istatus

! changed zloc to zin and zout, since the point of this routine
! is to convert zloc from one value to another.  ztype{in,out}
! are the vertical types as defined by the 3d sphere locations mod.
real(r8)            :: xloc, yloc, zin, xyz_loc(3), zout
integer             :: id, i, j, rc, ztypein, ztypeout
real(r8)            :: dx,dy,dxm,dym, dz, dzm
integer, dimension(2) :: ll, lr, ul, ur

real(r8), allocatable, dimension(:) :: v_h
real(r8), allocatable, dimension(:) :: v_p !< only need the mean value

! local vars, used in calculating pressure and height
real(r8)            :: pres1, pres2, pres3, pres4
real(r8)            :: presa, presb, psurf
real(r8)            :: hgt1, hgt2, hgt3, hgt4, hgta, hgtb

!HK
real(r8) :: zk
integer  :: k
logical  :: lev0

! assume failure.
istatus = 1

! first off, check if ob is identity ob.  if so get_state_meta_data() will 
! return location information already in the requested vertical type.
!>@todo This in not true anymore if you don't convert all the state variables 
! to the localization coordinate in get_state_meta_data
if (obs_kind < 0) then
   call get_state_meta_data(int(obs_kind,i8),location)
   istatus = 0
   return
endif

! if the existing coord is already in the requested vertical units
! or if the vert is 'undef' which means no specifically defined
! vertical coordinate, return now. 
ztypein  = nint(query_location(location, 'which_vert'))
ztypeout = vert_localization_coord
if ((ztypein == ztypeout) .or. (ztypein == VERTISUNDEF)) then
   istatus = 0
   return
endif

! we do need to convert the vertical.  start by
! extracting the location lat/lon/vert values.
xyz_loc = get_location(location)

!if(my_task_id() == 0) write(10, *) xyz_loc

! the routines below will use zin as the incoming vertical value
! and zout as the new outgoing one.  start out assuming failure
! (zout = missing) and wait to be pleasantly surprised when it works.
zin  = xyz_loc(3)
zout = missing_r8

! if the vertical is missing to start with, return it the same way
! with the requested type as out.
if (zin == missing_r8) then
   location = set_location(xyz_loc(1),xyz_loc(2),missing_r8,ztypeout)
   !print*, 'missing rank', my_task_id()
   return
endif

! first obtain domain id, and where we are in the grid (xloc,yloc)
if ( .not. scm ) then
   call get_domain_info(xyz_loc(1),xyz_loc(2),id,xloc,yloc)
else
   id = 1
   xloc = 1.0_r8
   yloc = 1.0_r8
endif
 
! cannot find domain info, return error.  set location to missing value
! but using requested vertical coord.  istatus already set above.
if (id==0) then
   location = set_location(xyz_loc(1),xyz_loc(2),missing_r8,ztypeout)
   !print*, 'can not find domain rank', my_task_id()

   return
endif

! get integer (west/south) grid point and distances to neighboring grid points
! distances are used as weights to carry out horizontal interpolations
call toGrid(xloc,i,dx,dxm)
call toGrid(yloc,j,dy,dym)

! Check that integer indices of Mass grid are in valid ranges for the given
!   boundary conditions (i.e., periodicity).  if not, bail here.
if ( .not. boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .or. &
     .not. boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) ) then
   location = set_location(xyz_loc(1),xyz_loc(2),missing_r8,ztypeout)
   !print*, 'out of bounds rank', my_task_id()

   return
endif

! Get indices of corners (i,i+1,j,j+1), which depend on periodicities.
! since the boundsCheck routines succeeded, this call should never fail
! so make it a fatal error if it does return an error code.
call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
if ( rc /= 0 ) then
   write(errstring,*) 'for i,j: ', i, j, ' getCorners rc = ', rc
   call error_handler(E_ERR,'model_mod.f90::vert_convert', errstring, &
                      source, revision, revdate)
endif

! at this point we have set:  xloc, yloc, i, j, ll, ul, lr, ur, zin, id,
! dx, dxm, dy, dym  already, and i, j have been checked to verify they 
! are valid values for this grid.  if you need k, dz, dzm below you still 
! need to compute and validate them first.

! Convert the incoming vertical type (ztypein) into the vertical 
! localization coordinate given in the namelist (ztypeout).

! convert from:
select case (ztypein)

! -------------------------------------------------------
! incoming vertical coordinate is 'model level number'
! ---------------------------------------------------- ---
case (VERTISLEVEL)

   ! convert into:
   select case (ztypeout)

   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'model level number'
   ! outgoing vertical coordinate should be 'pressure' in Pa
   ! -------------------------------------------------------
   case (VERTISPRESSURE)

      ! get neighboring mass level indices & compute weights to zin
      call toGrid(zin,k,dz,dzm)

      ! Check that integer height index is in valid range.  if not, bail to end
      if(.not. boundsCheck(k, .false., id, dim=3, type=wrf%dom(id)%type_t)) goto 100

      ! compute pressure at all neighboring mass points and interpolate
      presa = scalar(model_pressure_t_distrib(ll(1), ll(2), k  ,id, state_handle,1))
      presb = scalar(model_pressure_t_distrib(ll(1), ll(2), k+1,id, state_handle,1))
      pres1 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(lr(1), lr(2), k  ,id, state_handle,1))
      presb = scalar(model_pressure_t_distrib(lr(1), lr(2), k+1,id, state_handle,1))
      pres2 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(ul(1), ul(2), k  ,id, state_handle,1))
      presb = scalar(model_pressure_t_distrib(ul(1), ul(2), k+1,id, state_handle,1))
      pres3 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(ur(1), ur(2), k  ,id, state_handle,1))
      presb = scalar(model_pressure_t_distrib(ur(1), ur(2), k+1,id, state_handle,1))
      pres4 = dzm*presa + dz*presb
      zout = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'model level number'
   ! outgoing vertical coordinate should be 'height' in meters
   ! -------------------------------------------------------
   case (VERTISHEIGHT) 

      ! get neighboring mass level indices & compute weights to zin
      ! need to add half a grid to get to staggered vertical coordinate
      call toGrid(zin+0.5_r8,k,dz,dzm)
      !print*, 'k', k, 'zin ', zin

      ! Check that integer height index is in valid range.  if not, bail to end
      if(.not. boundsCheck(k, .false., id, dim=3, type=wrf%dom(id)%type_gz)) then
         ! print*, 'bounds check fail'
          goto 100
      endif

      ! HK if this is a piece of state, I believe you don't need to the four corners,
      ! the location is the lower left corner.
      ! compute height at all neighboring mass points and interpolate
      ! You have already converted the state in get_state_meta_data
      hgta = model_height_w_distrib(ll(1), ll(2), k  ,id,state_handle)
      hgtb = model_height_w_distrib(ll(1), ll(2), k+1,id,state_handle)
      hgt1 = dzm*hgta + dz*hgtb
      hgta = model_height_w_distrib(lr(1), lr(2), k  ,id,state_handle)
      hgtb = model_height_w_distrib(lr(1), lr(2), k+1,id,state_handle)
      hgt2 = dzm*hgta + dz*hgtb
      hgta = model_height_w_distrib(ul(1), ul(2), k  ,id,state_handle)
      hgtb = model_height_w_distrib(ul(1), ul(2), k+1,id,state_handle)
      hgt3 = dzm*hgta + dz*hgtb
      hgta = model_height_w_distrib(ur(1), ur(2), k  ,id,state_handle)
      hgtb = model_height_w_distrib(ur(1), ur(2), k+1,id,state_handle)
      hgt4 = dzm*hgta + dz*hgtb
      zout = dym*( dxm*hgt1 + dx*hgt2 ) + dy*( dxm*hgt3 + dx*hgt4 )

   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'model level number'
   ! outgoing vertical coordinate should be 'scale height' 
   ! -------------------------------------------------------
   case (VERTISSCALEHEIGHT)

      ! get neighboring mass level indices & compute weights to zin
      call toGrid(zin,k,dz,dzm)

      ! Check that integer height index is in valid range.  if not, bail to end
      if(.not. boundsCheck(k, .false., id, dim=3, type=wrf%dom(id)%type_t)) goto 100

      ! pressure at height
      presa = scalar(model_pressure_t_distrib(ll(1), ll(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(ll(1), ll(2), k+1,id,state_handle,1))
      pres1 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(lr(1), lr(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(lr(1), lr(2), k+1,id,state_handle,1))
      pres2 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(ul(1), ul(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(ul(1), ul(2), k+1,id,state_handle,1))
      pres3 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(ur(1), ur(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(ur(1), ur(2), k+1,id,state_handle,1))
      pres4 = dzm*presa + dz*presb
      zout = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )

      ! surface pressure
      pres1 = model_pressure_s_distrib(ll(1), ll(2), id, state_handle)
      pres2 = model_pressure_s_distrib(lr(1), lr(2), id, state_handle)
      pres3 = model_pressure_s_distrib(ul(1), ul(2), id, state_handle) 
      pres4 = model_pressure_s_distrib(ur(1), ur(2), id, state_handle)
      zout = -log(zout / (dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )))


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'model level number'
   ! outgoing vertical coordinate is unrecognized
   ! -------------------------------------------------------
   case default 
      write(errstring,*) 'Requested vertical coordinate not recognized: ', ztypeout
      call error_handler(E_ERR,'vert_convert', errstring, &
                         source, revision, revdate,  &
                         text2='Incoming vertical coordinate was model level.')


   end select   ! incoming vert was model level

! -------------------------------------------------------
! incoming vertical coordinate is 'pressure' in Pa
! -------------------------------------------------------
case (VERTISPRESSURE)

   ! If obs is by pressure: get corresponding mass level zk,
   ! then get neighboring mass level indices and compute weights 

   ! get model pressure profile and
   ! get pressure vertical co-ordinate in model level number
   allocate(v_p(0:wrf%dom(id)%bt)) 
   !HK This has already been called in model interpolate
   ! - not for observations that were not in the assimilate catagory
   call get_model_pressure_profile_distrib(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,id,v_p, state_handle,1)

     !if (my_task_id() == 0) then
     !    write(10, *) v_p
     !endif

   call pres_to_zk(zin, v_p(:), wrf%dom(id)%bt, zk, lev0)
   deallocate(v_p)

     !if (my_task_id() == 0) then
     !    write(10, *) zk
     !endif

   ! if you cannot get a model level out of the pressure profile, bail to end
   if ( zk == missing_r8 ) goto 100 !HK I don't think zk = missing_r8 even if there is an error

   ! convert into:
   select case (ztypeout)

   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'pressure' in Pa
   ! outgoing vertical coordinate should be 'model level'
   ! -------------------------------------------------------
   case (VERTISLEVEL)
      ! pres_to_zk() above converted pressure into a real number
      ! of vertical model levels, including the fraction.
      zout = zk


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'pressure' in Pa
   ! outgoing vertical coordinate should be 'height' in meters
   ! -------------------------------------------------------
   case (VERTISHEIGHT)
      ! adding 0.5 to get to the staggered vertical grid
      ! because height is on staggered vertical grid
      call toGrid(zk+0.5, k, dz, dzm)

      ! Check that integer height index is in valid range.  if not, bail to end
      if(.not. boundsCheck(k, .false., id, dim=3, type=wrf%dom(id)%type_gz)) goto 100

      ! HK should check whether this is state or not. If it is state, there is no need to 
      ! do all four corners, it is the lower left corner.

      ! compute height at all neighboring mass points and interpolate
      hgta = model_height_w_distrib(ll(1), ll(2), k  ,id, state_handle)
      hgtb = model_height_w_distrib(ll(1), ll(2), k+1,id, state_handle)
      hgt1 = dzm*hgta + dz*hgtb
      hgta = model_height_w_distrib(lr(1), lr(2), k  ,id, state_handle)
      hgtb = model_height_w_distrib(lr(1), lr(2), k+1,id, state_handle)
      hgt2 = dzm*hgta + dz*hgtb
      hgta = model_height_w_distrib(ul(1), ul(2), k  ,id, state_handle)
      hgtb = model_height_w_distrib(ul(1), ul(2), k+1,id, state_handle)
      hgt3 = dzm*hgta + dz*hgtb
      hgta = model_height_w_distrib(ur(1), ur(2), k  ,id, state_handle)
      hgtb = model_height_w_distrib(ur(1), ur(2), k+1,id, state_handle)
      hgt4 = dzm*hgta + dz*hgtb
      zout = dym*( dxm*hgt1 + dx*hgt2 ) + dy*( dxm*hgt3 + dx*hgt4 )
      
     !if (my_task_id() == 0) then
     !    write(10, *) zout
     !endif

   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'pressure' in Pa
   ! outgoing vertical coordinate should be 'scale height' 
   ! -------------------------------------------------------
   case (VERTISSCALEHEIGHT)
      call toGrid(zk,k,dz,dzm)

      ! Check that integer height index is in valid range.  if not, bail to end
      if(.not. boundsCheck(k, .false., id, dim=3, type=wrf%dom(id)%type_t)) goto 100

      ! compute surface pressure at all neighboring mass points and interpolate
      pres1 = model_pressure_s_distrib(ll(1), ll(2), id, state_handle)
      pres2 = model_pressure_s_distrib(lr(1), lr(2), id, state_handle)
      pres3 = model_pressure_s_distrib(ul(1), ul(2), id, state_handle)
      pres4 = model_pressure_s_distrib(ur(1), ur(2), id, state_handle)
      zout = -log(zin / (dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )))

   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'pressure'
   ! outgoing vertical coordinate is unrecognized
   ! -------------------------------------------------------
   case default 
      write(errstring,*) 'Requested vertical coordinate not recognized: ', ztypeout
      call error_handler(E_ERR,'vert_convert', errstring, &
                         source, revision, revdate,  &
                         text2='Incoming vertical coordinate was pressure.')


   end select   ! incoming vert was pressure


! -------------------------------------------------------
! incoming vertical coordinate is 'height' in meters
! -------------------------------------------------------
case (VERTISHEIGHT)

   ! If obs is by height: get corresponding mass level zk,
   ! then get neighboring mass level indices and compute weights

   ! get model height profile and
   ! get height vertical co-ordinate in model level number 
   allocate(v_h(0:wrf%dom(id)%bt))
   call get_model_height_profile_distrib(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,id,v_h, state_handle,1)
   call height_to_zk(zin, v_h, wrf%dom(id)%bt,zk,lev0)
   deallocate(v_h)

   ! convert into:
   select case (ztypeout)

   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'height' in meters
   ! outgoing vertical coordinate should be 'model level'
   ! -------------------------------------------------------
   case (VERTISLEVEL)
      ! height_to_zk() above converted pressure into a real number
      ! of vertical model levels, including the fraction.
      zout = zk


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'height' in meters
   ! outgoing vertical coordinate should be 'pressure' in Pa
   ! -------------------------------------------------------
   case (VERTISPRESSURE)
      call toGrid(zk,k,dz,dzm)

      ! Check that integer height index is in valid range.  if not, bail to end
      if(.not. boundsCheck(k, .false., id, dim=3, type=wrf%dom(id)%type_t)) goto 100

      ! compute pressure at all neighboring mass points and interpolate
      presa = scalar(model_pressure_t_distrib(ll(1), ll(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(ll(1), ll(2), k+1,id,state_handle,1))
      pres1 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(lr(1), lr(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(lr(1), lr(2), k+1,id,state_handle,1))
      pres2 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(ul(1), ul(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(ul(1), ul(2), k+1,id,state_handle,1))
      pres3 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(ur(1), ur(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(ur(1), ur(2), k+1,id,state_handle,1))
      pres4 = dzm*presa + dz*presb
      zout = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'height' in meters
   ! outgoing vertical coordinate should be 'scale height'
   ! -------------------------------------------------------
   case (VERTISSCALEHEIGHT)
      call toGrid(zk,k,dz,dzm)

      ! Check that integer height index is in valid range.  if not, bail to end
      if(.not. boundsCheck(k, .false., id, dim=3, type=wrf%dom(id)%type_t)) goto 100

      ! pressure at height
      presa = scalar(model_pressure_t_distrib(ll(1), ll(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(ll(1), ll(2), k+1,id,state_handle,1))
      pres1 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(lr(1), lr(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(lr(1), lr(2), k+1,id,state_handle,1))
      pres2 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(ul(1), ul(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(ul(1), ul(2), k+1,id,state_handle,1))
      pres3 = dzm*presa + dz*presb
      presa = scalar(model_pressure_t_distrib(ur(1), ur(2), k  ,id,state_handle,1))
      presb = scalar(model_pressure_t_distrib(ur(1), ur(2), k+1,id,state_handle,1))
      pres4 = dzm*presa + dz*presb
      zout = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )
 
      ! surface pressure
      pres1 = model_pressure_s_distrib(ll(1), ll(2), id, state_handle)
      pres2 = model_pressure_s_distrib(lr(1), lr(2), id, state_handle)
      pres3 = model_pressure_s_distrib(ul(1), ul(2), id, state_handle) 
      pres4 = model_pressure_s_distrib(ur(1), ur(2), id, state_handle)
      zout = -log(zout / (dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )))


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'height' in meters
   ! outgoing vertical coordinate is unrecognized
   ! -------------------------------------------------------
   case default 
      write(errstring,*) 'Requested vertical coordinate not recognized: ', ztypeout
      call error_handler(E_ERR,'vert_convert', errstring, &
                         source, revision, revdate,  &
                         text2='Incoming vertical coordinate was height.')


   end select   ! incoming vert was height


! -------------------------------------------------------
! incoming vertical coordinate is 'scale height' 
! -------------------------------------------------------
case (VERTISSCALEHEIGHT)

   ! If obs is by scale height: compute the surface pressure, 
   ! get corresponding mass level zk, then get neighboring mass 
   ! level indices and compute weights

   pres1 = model_pressure_s_distrib(ll(1), ll(2), id,state_handle) 
   pres2 = model_pressure_s_distrib(lr(1), lr(2), id,state_handle)
   pres3 = model_pressure_s_distrib(ul(1), ul(2), id,state_handle) 
   pres4 = model_pressure_s_distrib(ur(1), ur(2), id,state_handle) 
   psurf = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )

   ! get model pressure profile and
   ! get pressure vertical co-ordinate in model level number
   allocate(v_p(0:wrf%dom(id)%bt))
   call get_model_pressure_profile_distrib(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,id,v_p, state_handle, 1)
   call pres_to_zk(exp(-zin)*psurf, v_p, wrf%dom(id)%bt,zk,lev0)
   deallocate(v_p)

   ! if you cannot get a model level out of the pressure profile, bail to end
   if (zk == missing_r8) goto 100

   ! convert into:
   select case (ztypeout)

   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'scale height'
   ! outgoing vertical coordinate should be 'model level'
   ! -------------------------------------------------------
   case (VERTISLEVEL)
      ! pres_to_zk() above converted scale height/pressure into 
      ! a real number of vertical model levels, including the fraction.
      zout = zk


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'scale height'
   ! outgoing vertical coordinate should be 'height' in meters
   ! -------------------------------------------------------
   case (VERTISHEIGHT)
      ! adding 0.5 to get to the staggered vertical grid
      ! because height is on staggered vertical grid
      call toGrid(zk+0.5,k,dz,dzm)

      ! Check that integer height index is in valid range.  if not, bail to end
      if(.not. boundsCheck(k, .false., id, dim=3, type=wrf%dom(id)%type_gz)) goto 100

      ! compute height at all neighboring mass points and interpolate
      hgta = model_height_w_distrib(ll(1), ll(2), k  ,id,state_handle)
      hgtb = model_height_w_distrib(ll(1), ll(2), k+1,id,state_handle)
      hgt1 = dzm*hgta + dz*hgtb
      hgta = model_height_w_distrib(lr(1), lr(2), k  ,id,state_handle)
      hgtb = model_height_w_distrib(lr(1), lr(2), k+1,id,state_handle)
      hgt2 = dzm*hgta + dz*hgtb
      hgta = model_height_w_distrib(ul(1), ul(2), k  ,id,state_handle)
      hgtb = model_height_w_distrib(ul(1), ul(2), k+1,id,state_handle)
      hgt3 = dzm*hgta + dz*hgtb
      hgta = model_height_w_distrib(ur(1), ur(2), k  ,id,state_handle)
      hgtb = model_height_w_distrib(ur(1), ur(2), k+1,id,state_handle)
      hgt4 = dzm*hgta + dz*hgtb
      zout = dym*( dxm*hgt1 + dx*hgt2 ) + dy*( dxm*hgt3 + dx*hgt4 )


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'scale height'
   ! outgoing vertical coordinate should be 'pressure' in Pa
   ! -------------------------------------------------------
   case (VERTISPRESSURE)
      zout = exp(-zin)*psurf


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'scale height'
   ! outgoing vertical coordinate is unrecognized
   ! -------------------------------------------------------
   case default 
      write(errstring,*) 'Requested vertical coordinate not recognized: ', ztypeout
      call error_handler(E_ERR,'vert_convert', errstring, &
                         source, revision, revdate,  &
                         text2='Incoming vertical coordinate was scale height.')


   end select   ! incoming vert was scale height

! -------------------------------------------------------
! incoming vertical coordinate is 'surface' (assumes zin is height in meters)
! -------------------------------------------------------
case(VERTISSURFACE)

   ! convert into:
   select case (ztypeout)

   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'surface'
   ! outgoing vertical coordinate should be 'model level'
   ! -------------------------------------------------------
   case (VERTISLEVEL)
      zout = 1.0_r8


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'surface'
   ! outgoing vertical coordinate should be 'pressure' in Pa
   ! -------------------------------------------------------
   case (VERTISPRESSURE)

      ! compute surface pressure at all neighboring mass points
      pres1 = model_pressure_s_distrib(ll(1), ll(2), id, state_handle)
      pres2 = model_pressure_s_distrib(lr(1), lr(2), id, state_handle)
      pres3 = model_pressure_s_distrib(ul(1), ul(2), id, state_handle)
      pres4 = model_pressure_s_distrib(ur(1), ur(2), id, state_handle)
      zout = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'surface'
   ! outgoing vertical coordinate should be 'scale height' 
   ! -------------------------------------------------------
   case (VERTISSCALEHEIGHT)
      zout = -log(1.0_r8)


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'surface'
   ! outgoing vertical coordinate should be 'height' in meters
   ! -------------------------------------------------------
   case (VERTISHEIGHT)
      ! a surface ob is assumed to have height as vertical coordinate.
      ! this code needs to be revised if this is not true 
      ! (in that case uncomment lines below to get terrain height
      ! from model)
      zout = zin
      !! or: directly interpolate terrain height at neighboring mass points
      !zout = dym*( dxm*wrf%dom(id)%hgt(i,  j) + &
      !              dx*wrf%dom(id)%hgt(i+1,j) ) + &
      !        dy*( dxm*wrf%dom(id)%hgt(i,  j+1) + &
      !              dx*wrf%dom(id)%hgt(i+1,j+1) )


   ! -------------------------------------------------------
   ! incoming vertical coordinate is 'surface'
   ! outgoing vertical coordinate is unrecognized
   ! -------------------------------------------------------
   case default 
      write(errstring,*) 'Requested vertical coordinate not recognized: ', ztypeout
      call error_handler(E_ERR,'vert_convert', errstring, &
                         source, revision, revdate,  &
                         text2='Incoming vertical coordinate was surface.')


   end select   ! incoming vert was surface

! -------------------------------------------------------
! incoming vertical coordinate has no case section
! -------------------------------------------------------
case default
   write(errstring,*) 'Incoming vertical coordinate type not recognized: ',ztypein
   call error_handler(E_ERR,'vert_convert', errstring, &
        source, revision, revdate)

end select   ! incoming z vertical type


! on error, come here.  istatus was set to 1 and zout to missing_r8
! so unless they have been reset to good values, things did not work.
100 continue 

! Returned location 
location = set_location(xyz_loc(1),xyz_loc(2),zout,ztypeout)

!if(my_task_id() == 0 ) print*, xyz_loc(1), xyz_loc(2), zout, ztypeout

! Set successful return code only if zout has good value
if(zout /= missing_r8) istatus = 0

end subroutine vert_convert

!#######################################################################

subroutine get_wrf_horizontal_location( i, j, var_type, id, long, lat )

integer,  intent(in)  :: i,j,var_type, id
real(r8), intent(out) :: long, lat

! given i, j indices into the horizontal grid return the lat/long.
! if u or v staggering use the staggered grids, otherwise use the mass 
! grid.  this code has changed -- earlier versions only had the mass 
! grid available and used it to compute cell midpoints and called them 
! the staggered points.  now that all three grids are being read, look 
! up the point locations directly from the appropriate array.

if (var_type == wrf%dom(id)%type_u) then
   long = wrf%dom(id)%longitude_u(i,j)
   lat  = wrf%dom(id)%latitude_u(i,j)
elseif (var_type == wrf%dom(id)%type_v) then
   long = wrf%dom(id)%longitude_v(i,j)
   lat  = wrf%dom(id)%latitude_v(i,j)
else
   long = wrf%dom(id)%longitude(i,j)
   lat  = wrf%dom(id)%latitude(i,j)
endif

do while (long <   0.0_r8)
   long = long + 360.0_r8
end do
do while (long > 360.0_r8)
   long = long - 360.0_r8
end do

end subroutine get_wrf_horizontal_location



!***********************************************************************


subroutine nc_write_model_atts( ncid, id ) 
!-----------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: id

logical, parameter :: write_precip = .false.

integer :: weDimID, weStagDimID, snDimID, snStagDimID, &
     btDimID, btStagDimID, slSDimID

integer :: DXVarID, DYVarID, TRUELAT1VarID, TRUELAT2VarID, STAND_LONVarID
integer :: CEN_LATVarID, CEN_LONVarID, MAP_PROJVarID
integer :: PERIODIC_XVarID, POLARVarID

integer :: DNVarID, ZNUVarID, DNWVarID, phbVarID, &
     MubVarID, LonVarID, LatVarID, ilevVarID, XlandVarID, hgtVarID , LatuVarID, &
     LatvVarID, LonuVarID, LonvVarID, ZNWVarID, C1HVarID, C2HVarID

integer :: TimeDimID

! currently unused, but if needed could be added back in.  these fields
! only appear to be supported in certain projections, so the code should
! test to be sure they exist before trying to read them from the netcdf file.
!integer, dimension(num_domains) :: MapFacMVarID, MapFacUVarID, MapFacVVarID

integer :: var_id
integer :: i, ret, tmp

character(len=129) :: title
character(len=32) :: context = 'nc_write_model_atts'

character (len=1)     :: idom

character(len=129), allocatable, dimension(:) :: textblock
integer :: ind, my_index
character(len=NF90_MAX_NAME) :: attname, varname
character(len=129) :: unitsval, descriptionval, coordinatesval, long_nameval, coordinate_char
logical               :: debug = .false.
character(len=256) :: filename



! use netcdf file id for identification
write(filename,*) 'ncid', ncid

!-------------------------------------------------------------------------------
! Put file into define mode and
! Write Global Attributes 
!-------------------------------------------------------------------------------
call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate", revdate )

call nc_add_global_attribute(ncid, "model", "wrf")

!-----------------------------------------------------------------
! Define the dimensions IDs
!-----------------------------------------------------------------

!>@todo all the wrf files use Time as the first dimension
!> but if we are creating a file from scratch, we need this
!> to be made
ret = nf90_inq_dimid(ncid, "Time", TimeDimID)
if (ret /= NF90_NOERR) then
   call nc_check(nf90_def_dim(ncid=ncid, name="Time", &
                 len = 1,  dimid = TimeDimID), &
                 'nc_write_model_atts','def_dim domain')
endif


!>@todo FIXME all the variables below should have Time as
!> the first dimension.

call nc_check(nf90_def_dim(ncid=ncid, name='west_east', &
              len = wrf%dom(id)%we,  dimid = weDimID), &
              'nc_write_model_atts','def_dim west_east')
call nc_check(nf90_def_dim(ncid=ncid, name='west_east_stag',   &
              len = wrf%dom(id)%wes, dimid = weStagDimID), &
              'nc_write_model_atts','def_dim west_east_stag')
call nc_check(nf90_def_dim(ncid=ncid, name='south_north',      &
              len = wrf%dom(id)%sn,  dimid = snDimID), &
              'nc_write_model_atts','def_dim south_north')
call nc_check(nf90_def_dim(ncid=ncid, name='south_north_stag', &
              len = wrf%dom(id)%sns, dimid = snStagDimID), &
              'nc_write_model_atts','def_dim south_north_stag')
call nc_check(nf90_def_dim(ncid=ncid, name='bottom_top',       &
              len = wrf%dom(id)%bt,  dimid = btDimID), &
              'nc_write_model_atts','def_dim bottom_top')
call nc_check(nf90_def_dim(ncid=ncid, name='bottom_top_stag',  &
              len = wrf%dom(id)%bts, dimid = btStagDimID), &
              'nc_write_model_atts','def_dim bottom_top_stag')
call nc_check(nf90_def_dim(ncid=ncid, name='soil_layers_stag',  &
              len = wrf%dom(id)%sls, dimid = slSDimID), &
              'nc_write_model_atts','def_dim soil_layers_stag')

!-----------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Create the (empty) static variables and their attributes
! Commented block is from wrfinput
!-----------------------------------------------------------------

call nc_check(nf90_def_var(ncid, name='DX', xtype=nf90_real, &
               varid=DXVarID), &
              'nc_write_model_atts','def_var DX')
call nc_check(nf90_put_att(ncid, DXVarID, 'long_name', 'X HORIZONTAL RESOLUTION'), &
              'nc_write_model_atts','put_att DX long_name')
call nc_check(nf90_put_att(ncid, DXVarID, 'description', 'X HORIZONTAL RESOLUTION'), &
              'nc_write_model_atts','put_att DX description')
call nc_check(nf90_put_att(ncid, DXVarID, 'units', 'm'), &
              'nc_write_model_atts','put_att DX units')

call nc_check(nf90_def_var(ncid, name='DY', xtype=nf90_real, &
               varid=DYVarID), &
              'nc_write_model_atts','def_var DY')
call nc_check(nf90_put_att(ncid, DYVarID, 'long_name', 'Y HORIZONTAL RESOLUTION'), &
              'nc_write_model_atts','put_att DY long_name')
call nc_check(nf90_put_att(ncid, DYVarID, 'description', 'Y HORIZONTAL RESOLUTION'), &
              'nc_write_model_atts','put_att DY description')
call nc_check(nf90_put_att(ncid, DYVarID, 'units', 'm'), &
              'nc_write_model_atts','put_att DY units')

call nc_check(nf90_def_var(ncid, name='TRUELAT1', xtype=nf90_real, &
               varid=TRUELAT1VarID), &
              'nc_write_model_atts','def_var TRUELAT1')
call nc_check(nf90_put_att(ncid, TRUELAT1VarID, 'long_name', &
              'first standard parallel'), &
              'nc_write_model_atts','put_att TRUELAT1 long_name')
call nc_check(nf90_put_att(ncid, TRUELAT1VarID, 'description', &
              'first standard parallel'), &
              'nc_write_model_atts','put_att TRUELAT1 description')
call nc_check(nf90_put_att(ncid, TRUELAT1VarID, 'units', &
              'degrees, negative is south'), &
              'nc_write_model_atts','put_att TRUELAT1 units')

call nc_check(nf90_def_var(ncid, name='TRUELAT2', xtype=nf90_real, &
               varid=TRUELAT2VarID), &
              'nc_write_model_atts','def_var TRUELAT2')
call nc_check(nf90_put_att(ncid, TRUELAT2VarID, 'long_name', &
              'second standard parallel'), &
              'nc_write_model_atts','put_att TRUELAT2 long_name')
call nc_check(nf90_put_att(ncid, TRUELAT2VarID, 'description', &
              'second standard parallel'), &
              'nc_write_model_atts','put_att TRUELAT2 description')
call nc_check(nf90_put_att(ncid, TRUELAT2VarID, 'units', &
              'degrees, negative is south'), &
              'nc_write_model_atts','put_att TRUELAT2 units')

call nc_check(nf90_def_var(ncid, name='STAND_LON', xtype=nf90_real, &
               varid=STAND_LONVarID), &
              'nc_write_model_atts','def_var STAND_LON')
call nc_check(nf90_put_att(ncid, STAND_LONVarID, 'long_name', &
              'standard longitude'), &
              'nc_write_model_atts','put_att STAND_LON long_name')
call nc_check(nf90_put_att(ncid, STAND_LONVarID, 'description', &
              'standard longitude'), &
              'nc_write_model_atts','put_att STAND_LON description')
call nc_check(nf90_put_att(ncid, STAND_LONVarID, 'units', &
              'degrees, negative is west'), &
              'nc_write_model_atts','put_att STAND_LON units')

call nc_check(nf90_def_var(ncid, name='CEN_LAT', xtype=nf90_real, &
               varid=CEN_LATVarID), &
              'nc_write_model_atts','def_var CEN_LAT')
call nc_check(nf90_put_att(ncid, CEN_LATVarID, 'long_name', &
              'center latitude'), &
              'nc_write_model_atts','put_att CEN_LAT long_name')
call nc_check(nf90_put_att(ncid, CEN_LATVarID, 'description', &
              'center latitude'), &
              'nc_write_model_atts','put_att CEN_LAT description')
call nc_check(nf90_put_att(ncid, CEN_LATVarID, 'units', &
              'degrees, negative is south'), &
              'nc_write_model_atts','put_att CEN_LAT units')

call nc_check(nf90_def_var(ncid, name='CEN_LON', xtype=nf90_real, &
               varid=CEN_LONVarID), &
              'nc_write_model_atts','def_var CEN_LON')
call nc_check(nf90_put_att(ncid, CEN_LONVarID, 'long_name', &
              'central longitude'), &
              'nc_write_model_atts','put_att CEN_LON long_name')
call nc_check(nf90_put_att(ncid, CEN_LONVarID, 'description', &
              'central longitude'), &
              'nc_write_model_atts','put_att CEN_LON description')
call nc_check(nf90_put_att(ncid, CEN_LONVarID, 'units', &
              'degrees, negative is west'), &
              'nc_write_model_atts','put_att CEN_LON units')

call nc_check(nf90_def_var(ncid, name='MAP_PROJ', xtype=nf90_real, &
               varid=MAP_PROJVarID), &
              'nc_write_model_atts','def_var MAP_PROJ')
call nc_check(nf90_put_att(ncid, MAP_PROJVarID, 'long_name', &
              'domain map projection'), &
              'nc_write_model_atts','put_att MAP_PROJ long_name')
call nc_check(nf90_put_att(ncid, MAP_PROJVarID, 'description', &
              'domain map projection'), &
              'nc_write_model_atts','put_att MAP_PROJ description')
call nc_check(nf90_put_att(ncid, MAP_PROJVarID, 'units', &
              '0=none, 1=Lambert, 2=polar, 3=Mercator, 5=Cylindrical, 6=Cassini'), &
              'nc_write_model_atts','put_att MAP_PROJ units')

!nc -- we need to add in code here to report the domain values for the 
!        boundary condition flags periodic_x and polar.  Since these are
!        carried internally as logicals, they will first need to be 
!        converted back to integers.
call nc_check(nf90_def_var(ncid, name='PERIODIC_X', xtype=nf90_int, &
               varid=PERIODIC_XVarID), &
              'nc_write_model_atts','def_var PERIODIC_X')
call nc_check(nf90_put_att(ncid, PERIODIC_XVarID, 'long_name', &
              'Longitudinal periodic b.c. flag'), &
              'nc_write_model_atts','put_att PERIODIC_X long_name')
call nc_check(nf90_put_att(ncid, PERIODIC_XVarID, 'description', &
              'Longitudinal periodic b.c. flag'), &
              'nc_write_model_atts','put_att PERIODIC_X desciption')
call nc_check(nf90_put_att(ncid, PERIODIC_XVarID, 'units', &
              'logical: 1 = .true., 0 = .false.'), &
              'nc_write_model_atts','put_att PERIODIC_X units')

call nc_check(nf90_def_var(ncid, name='POLAR', xtype=nf90_int, &
               varid=POLARVarID), &
              'nc_write_model_atts','def_var POLAR')
call nc_check(nf90_put_att(ncid, POLARVarID, 'long_name', &
              'Polar periodic b.c. flag'), &
              'nc_write_model_atts','put_att POLAR long_name')
call nc_check(nf90_put_att(ncid, POLARVarID, 'description', &
              'Polar periodic b.c. flag'), &
              'nc_write_model_atts','put_att POLAR description')
call nc_check(nf90_put_att(ncid, POLARVarID, 'units', &
              'logical: 1 = .true., 0 = .false.'), &
              'nc_write_model_atts','put_att POLAR units')



write( idom , '(I1)') id

call nc_check(nf90_def_var(ncid, name='DN', xtype=nf90_real, &
                 dimids= btDimID, varid=DNVarID), &
                 'nc_write_model_atts','def_var DN_do'//idom)
call nc_check(nf90_put_att(ncid, DNVarID, 'long_name', &
                 'dn values on half (mass) levels'), &
                 'nc_write_model_atts','put_att DN_do'//idom//' long_name')
call nc_check(nf90_put_att(ncid, DNVarID, 'description', &
                 'dn values on half (mass) levels'), &
                 'nc_write_model_atts','put_att DN_do'//idom//' description')
call nc_check(nf90_put_att(ncid, DNVarID, 'units', &
                 ''), &
                 'nc_write_model_atts','put_att DN_do'//idom//' units')

call nc_check(nf90_def_var(ncid, name='ZNU', xtype=nf90_real, &
                 dimids= btDimID, varid=ZNUVarID), &
                 'nc_write_model_atts','def_var ZNU')
call nc_check(nf90_put_att(ncid, ZNUVarID, 'long_name', &
                 'eta values on half (mass) levels'), &
                 'nc_write_model_atts','put_att ZNU'//' long_name')
call nc_check(nf90_put_att(ncid, ZNUVarID, 'description', &
                 'eta values on half (mass) levels'), &
                 'nc_write_model_atts','put_att ZNU'//' description')
call nc_check(nf90_put_att(ncid, ZNUVarID, 'units', &
                 ''), &
                 'nc_write_model_atts','put_att ZNU'//' units')

call nc_check(nf90_def_var(ncid, name='ZNW', xtype=nf90_real, &
                 dimids= btStagDimID, varid=ZNWVarID), &
                 'nc_write_model_atts','def_var ZNW')
call nc_check(nf90_put_att(ncid, ZNWVarID, 'long_name', &
                 'eta values on full (mass) levels'), &
                 'nc_write_model_atts','put_att ZNW'//' long_name')
call nc_check(nf90_put_att(ncid, ZNWVarID, 'description', &
                 'eta values on full (mass) levels'), &
                 'nc_write_model_atts','put_att ZNW'//' description')
call nc_check(nf90_put_att(ncid, ZNWVarID, 'units', &
                 ''), &
                 'nc_write_model_atts','put_att ZNW'//' units')

call nc_check(nf90_def_var(ncid, name='DNW', xtype=nf90_real, &
                 dimids= btDimID, varid=DNWVarID), &
                 'nc_write_model_atts','def_var DNW')
call nc_check(nf90_put_att(ncid, DNWVarID, 'long_name', &
                 'dn values on full (w) levels'), &
                 'nc_write_model_atts','def_var DNW'//' long_name')
call nc_check(nf90_put_att(ncid, DNWVarID, 'description', &
                 'dn values on full (w) levels'), &
                 'nc_write_model_atts','def_var DNW'//' description')
call nc_check(nf90_put_att(ncid, DNWVarID, 'units', &
                 ''), &
                 'nc_write_model_atts','def_var DNW'//' units')

if (wrf%dom(id)%hybrid_opt == VERT_HYBRID) then
   call nc_check(nf90_def_var(ncid, name='C1H', xtype=nf90_real, &
                    dimids= btDimID, varid=C1HVarID), &
                    'nc_write_model_atts','def_var C1H')
   call nc_check(nf90_put_att(ncid, C1HVarID, 'long_name', &
                    'dn values on full (w) levels'), &
                    'nc_write_model_atts','def_var C1H'//' long_name')
   call nc_check(nf90_put_att(ncid, C1HVarID, 'description', &
                    'dn values on full (w) levels'), &
                    'nc_write_model_atts','def_var C1H'//' description')
   call nc_check(nf90_put_att(ncid, C1HVarID, 'units', &
                    ''), &
                    'nc_write_model_atts','def_var C1H'//' units')
   
   call nc_check(nf90_def_var(ncid, name='C2H', xtype=nf90_real, &
                    dimids= btDimID, varid=C2HVarID), &
                    'nc_write_model_atts','def_var C2H')
   call nc_check(nf90_put_att(ncid, C2HVarID, 'long_name', &
                    'dn values on full (w) levels'), &
                    'nc_write_model_atts','def_var C2H'//' long_name')
   call nc_check(nf90_put_att(ncid, C2HVarID, 'description', &
                    'dn values on full (w) levels'), &
                    'nc_write_model_atts','def_var C2H'//' description')
   call nc_check(nf90_put_att(ncid, C2HVarID, 'units', &
                    ''), &
                    'nc_write_model_atts','def_var C2H'//' units')
endif

!
!    float MUB(Time, south_north, west_east) ;
!            MUB:FieldType = 104 ;
!            MUB:MemoryOrder = 'XY ' ;
!            MUB:stagger = '' ;
             coordinate_char = 'XLONG'//' XLAT'
call nc_check(nf90_def_var(ncid, name='MUB', xtype=nf90_real, &
                 dimids= (/ weDimID, snDimID /), varid=MubVarID), &
                 'nc_write_model_atts','def_var MUB')
call nc_check(nf90_put_att(ncid, MubVarID, 'long_name', &
                 'base state dry air mass in column'), &
                 'nc_write_model_atts','put_att MUB'//' long_name')
call nc_check(nf90_put_att(ncid, MubVarID, 'description', &
                 'base state dry air mass in column'), &
                 'nc_write_model_atts','put_att MUB'//' description')
call nc_check(nf90_put_att(ncid, MubVarID, 'units', 'Pa'), &
                 'nc_write_model_atts','put_att MUB'//' units')
call nc_check(nf90_put_att(ncid, MubVarID, 'coordinates', &
                 trim(coordinate_char)), &
                 'nc_write_model_atts','put_att MUB'//' coordinates')

! Longitudes
!      float XLONG(Time, south_north, west_east) ;
!         XLONG:FieldType = 104 ;
!         XLONG:MemoryOrder = 'XY ' ;
!         XLONG:stagger = '' ;
   call nc_check(nf90_def_var(ncid, name='XLONG', xtype=nf90_real, &
                 dimids= (/ weDimID, snDimID /), varid=LonVarID),  &
                 'nc_write_model_atts','def_var XLONG')
   call nc_check(nf90_put_att(ncid, LonVarID, 'long_name', &
                 'LONGITUDE, WEST IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLONG'//' long_name')
   call nc_check(nf90_put_att(ncid, LonVarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts','put_att XLONG'//' units')
   call nc_check(nf90_put_att(ncid, LonVarID, 'valid_range', &
                 (/ -180.0_r8, 180.0_r8 /)), &
                 'nc_write_model_atts','put_att XLONG'//' valid_range')
   call nc_check(nf90_put_att(ncid, LonVarID, 'description', &
                 'LONGITUDE, WEST IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLONG'//' description')

!      float XLONG_U(Time, south_north, west_east_stag) ;
!         XLONG:FieldType = 104 ;
!         XLONG:MemoryOrder = 'XY ' ;
!         XLONG:stagger = '' ;
   call nc_check(nf90_def_var(ncid, name='XLONG_U', xtype=nf90_real, &
                 dimids= (/ weStagDimID, snDimID /), varid=LonuVarID),  &
                 'nc_write_model_atts','def_var XLONG_U')
   call nc_check(nf90_put_att(ncid, LonuVarID, 'long_name', &
                 'LONGITUDE, WEST IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLONG_U'//' long_name')
   call nc_check(nf90_put_att(ncid, LonuVarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts','put_att XLONG_U'//' units')
   call nc_check(nf90_put_att(ncid, LonuVarID, 'valid_range', &
                 (/ -180.0_r8, 180.0_r8 /)), &
                 'nc_write_model_atts','put_att XLONG_U'//' valid_range')
   call nc_check(nf90_put_att(ncid, LonuVarID, 'description', &
                 'LONGITUDE, WEST IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLONG_U'//' description')

!      float XLONG_V(Time, south_north_stag, west_east) ;
!         XLONG:FieldType = 104 ;
!         XLONG:MemoryOrder = 'XY ' ;
!         XLONG:stagger = '' ;
call nc_check(nf90_def_var(ncid, name='XLONG_V', xtype=nf90_real, &
                 dimids= (/ weDimID, snStagDimID /), varid=LonvVarID),  &
                 'nc_write_model_atts','def_var XLONG_V')
call nc_check(nf90_put_att(ncid, LonvVarID, 'long_name', &
                 'LONGITUDE, WEST IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLONG_V'//' long_name')
call nc_check(nf90_put_att(ncid, LonvVarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts','put_att XLONG_V'//' units')
call nc_check(nf90_put_att(ncid, LonvVarID, 'valid_range', &
                 (/ -180.0_r8, 180.0_r8 /)), &
                 'nc_write_model_atts','put_att XLONG_V'//' valid_range')
call nc_check(nf90_put_att(ncid, LonvVarID, 'description', &
                 'LONGITUDE, WEST IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLONG_V'//' description')

! Latitudes
!      float XLAT(Time, south_north, west_east) ;
!         XLAT:FieldType = 104 ;
!         XLAT:MemoryOrder = 'XY ' ;
!         XLAT:stagger = '' ;
call nc_check(nf90_def_var(ncid, name='XLAT', xtype=nf90_real, &
                 dimids=(/ weDimID, snDimID /), varid=LatVarID), &
                 'nc_write_model_atts','def_var XLAT') 
call nc_check(nf90_put_att(ncid, LatVarID, 'long_name', &
                 'LATITUDE, SOUTH IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLAT'//' long_name')
call nc_check(nf90_put_att(ncid, LatVarID, 'units', 'degrees_north'), &
                 'nc_write_model_atts','put_att XLAT'//' units')
call nc_check(nf90_put_att(ncid, LatVarID, 'valid_range', &
                 (/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts','put_att XLAT'//' valid_range')
call nc_check(nf90_put_att(ncid, LatVarID, 'description', &
                 'LATITUDE, SOUTH IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLAT'//' description')

!      float XLAT_U(Time, south_north, west_east_stag) ;
!         XLAT_U:FieldType = 104 ;
!         XLAT_U:MemoryOrder = 'XY ' ;
!         XLAT_U:stagger = '' ;
call nc_check(nf90_def_var(ncid, name='XLAT_U', xtype=nf90_real, &
                 dimids=(/ weStagDimID, snDimID /), varid=LatuVarID), &
                 'nc_write_model_atts','def_var XLAT_U') 
call nc_check(nf90_put_att(ncid, LatuVarID, 'long_name', &
                 'LATITUDE, SOUTH IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLAT_U'//' long_name')
call nc_check(nf90_put_att(ncid, LatuVarID, 'units', 'degrees_north'), &
                 'nc_write_model_atts','put_att XLAT_U'//' units')
call nc_check(nf90_put_att(ncid, LatuVarID, 'valid_range', &
                 (/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts','put_att XLAT_U'//' valid_range')
call nc_check(nf90_put_att(ncid, LatuVarID, 'description', &
                 'LATITUDE, SOUTH IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLAT_U'//' description')

!      float XLAT_V(Time, south_north_stag, west_east) ;
!         XLAT_V:FieldType = 104 ;
!         XLAT_V:MemoryOrder = 'XY ' ;
!         XLAT_V:stagger = '' ;
call nc_check(nf90_def_var(ncid, name='XLAT_V', xtype=nf90_real, &
                 dimids=(/ weDimID, snStagDimID /), varid=LatvVarID), &
                 'nc_write_model_atts','def_var XLAT_V') 
call nc_check(nf90_put_att(ncid, LatvVarID, 'long_name', &
                 'LATITUDE, SOUTH IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLAT_V'//' long_name')
call nc_check(nf90_put_att(ncid, LatvVarID, 'units', 'degrees_north'), &
                 'nc_write_model_atts','put_att XLAT_V'//' units')
call nc_check(nf90_put_att(ncid, LatvVarID, 'valid_range', &
                 (/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts','put_att XLAT_V'//' valid_range')
call nc_check(nf90_put_att(ncid, LatvVarID, 'description', &
                 'LATITUDE, SOUTH IS NEGATIVE'), &
                 'nc_write_model_atts','put_att XLAT_V'//' description')

! grid levels
call nc_check(nf90_def_var(ncid, name='level', xtype=nf90_short, &
                 dimids=btDimID, varid=ilevVarID), &
                 'nc_write_model_atts','def_var level')
call nc_check(nf90_put_att(ncid, ilevVarID, 'long_name', &
                 'level index'), &
                 'nc_write_model_atts','put_att level'//' long_name')
call nc_check(nf90_put_att(ncid, ilevVarID, 'description', &
                 'level index'), &
                 'nc_write_model_atts','put_att level'//' description')
call nc_check(nf90_put_att(ncid, ilevVarID, 'units', &
                 ''), &
                 'nc_write_model_atts','put_att level'//' units')

! Land Mask
!    float XLAND(Time, south_north, west_east) ;
!            XLAND:FieldType = 104 ;
!            XLAND:MemoryOrder = 'XY ' ;
!            XLAND:units = 'NA' ;
!            XLAND:stagger = '' ;
             coordinate_char = 'XLONG'//' XLAT'
call nc_check(nf90_def_var(ncid, name='XLAND', xtype=nf90_short, &
                 dimids= (/ weDimID, snDimID /), varid=XlandVarID),  &
                 'nc_write_model_atts','def_var XLAND')
call nc_check(nf90_put_att(ncid, XlandVarID, 'long_name', &
                 'LAND MASK (1 FOR LAND, 2 FOR WATER)'), &
                 'nc_write_model_atts','put_att XLAND'//' long_name')
call nc_check(nf90_put_att(ncid, XlandVarID, 'units', ' '), &
                 'nc_write_model_atts','put_att XLAND'//' units')
call nc_check(nf90_put_att(ncid, XlandVarID, 'coordinates', &
                 trim(coordinate_char)), &
                 'nc_write_model_atts','put_att XLAND'//' coordinates')
call nc_check(nf90_put_att(ncid, XlandVarID, 'valid_range', (/ 1, 2 /)), &
                 'nc_write_model_atts','put_att XLAND'//' valid_range')
call nc_check(nf90_put_att(ncid, XlandVarID, 'description', &
                 'LAND MASK (1 FOR LAND, 2 FOR WATER)'), &
                 'nc_write_model_atts','put_att XLAND'//' description')

! PHB
!    float PHB(Time, bottom_top_stag, south_north, west_east) ;
!            PHB:FieldType = 104 ;
!            PHB:MemoryOrder = 'XYZ' ;
!            PHB:stagger = 'Z' ;
             coordinate_char = 'XLONG'//' XLAT'
call nc_check(nf90_def_var(ncid, name='PHB', xtype=nf90_real, &
        dimids= (/ weDimID, snDimID, btStagDimID /), varid=phbVarId), &
        'nc_write_model_atts','def_var PHB')
call nc_check(nf90_put_att(ncid, phbVarId, 'long_name', &
                 'base-state geopotential'), &
                 'nc_write_model_atts','put_att PHB'//' long_name')
call nc_check(nf90_put_att(ncid, phbVarId, 'description', &
                 'base-state geopotential'), &
                 'nc_write_model_atts','put_att PHB'//' description')
call nc_check(nf90_put_att(ncid, phbVarId, 'units', 'm2/s2'), &
                 'nc_write_model_atts','put_att PHB'//' units')
call nc_check(nf90_put_att(ncid, phbVarId, 'coordinates', &
                 trim(coordinate_char)), &
                 'nc_write_model_atts','put_att PHB'//' coordinates')
call nc_check(nf90_put_att(ncid, phbVarId, 'units_long_name', 'm{2} s{-2}'), &
                 'nc_write_model_atts','put_att PHB'//' units_long_name')

             coordinate_char = 'XLONG'//' XLAT'
call nc_check(nf90_def_var(ncid, name='HGT', xtype=nf90_real, &
                 dimids= (/ weDimID, snDimID /), varid=hgtVarId), &
                 'nc_write_model_atts','def_var HGT')
call nc_check(nf90_put_att(ncid, hgtVarId, 'long_name', 'Terrain Height'), &
                 'nc_write_model_atts','put_att HGT'//' long_name')
call nc_check(nf90_put_att(ncid, hgtVarId, 'description', 'Terrain Height'), &
                 'nc_write_model_atts','put_att HGT'//' description')
call nc_check(nf90_put_att(ncid, hgtVarId, 'units', 'm'), &
                 'nc_write_model_atts','put_att HGT'//' units')
call nc_check(nf90_put_att(ncid, hgtVarId, 'coordinates', &
                 trim(coordinate_char)), &
                 'nc_write_model_atts','put_att HGT'//' coordinates')
call nc_check(nf90_put_att(ncid, hgtVarId, 'units_long_name', 'meters'), &
                 'nc_write_model_atts','put_att HGT'//' units_long_name')

! Leave define mode so we can actually fill the variables.
call nc_end_define_mode(ncid)

!-----------------------------------------------------------------
! Fill the variables we can
!-----------------------------------------------------------------

call nc_check(nf90_put_var(ncid,        DXVarID, wrf%dom(id)%dx), &
              'nc_write_model_atts','put_var dx')
call nc_check(nf90_put_var(ncid,        DYVarID, wrf%dom(id)%dy), &
              'nc_write_model_atts','put_var dy')
call nc_check(nf90_put_var(ncid,  TRUELAT1VarID, wrf%dom(id)%proj%truelat1), &
              'nc_write_model_atts','put_var truelat1')
call nc_check(nf90_put_var(ncid,  TRUELAT2VarID, wrf%dom(id)%proj%truelat2), &
              'nc_write_model_atts','put_var truelat2')
call nc_check(nf90_put_var(ncid, STAND_LONVarID, wrf%dom(id)%proj%stdlon), &
              'nc_write_model_atts','put_var stdlon')
call nc_check(nf90_put_var(ncid,   CEN_LATVarID, wrf%dom(id)%cen_lat), &
              'nc_write_model_atts','put_var cen_lat')
call nc_check(nf90_put_var(ncid,   CEN_LONVarID, wrf%dom(id)%cen_lon), &
              'nc_write_model_atts','put_var cen_lon')
call nc_check(nf90_put_var(ncid,  MAP_PROJVarID, wrf%dom(id)%map_proj), &
              'nc_write_model_atts','put_var map_proj')

!nc -- convert internally logical boundary condition variables into integers before filling
if ( wrf%dom(id)%periodic_x ) then
   tmp = 1
else
   tmp = 0
endif
call nc_check(nf90_put_var(ncid, PERIODIC_XVarID, tmp ), &
              'nc_write_model_atts','put_var PERIODIC_XVarID')

if ( wrf%dom(id)%polar ) then
   tmp = 1
else
   tmp = 0
endif
call nc_check(nf90_put_var(ncid, POLARVarID, tmp ), &
              'nc_write_model_atts','put var POLARVarID')


! defining grid levels
call nc_check(nf90_put_var(ncid,       DNVarID, wrf%dom(id)%dn), &
              'nc_write_model_atts','put_var dn')
call nc_check(nf90_put_var(ncid,      ZNUVarID, wrf%dom(id)%znu), &
              'nc_write_model_atts','put_var znu')
call nc_check(nf90_put_var(ncid,      ZNWVarID, wrf%dom(id)%znw), &
              'nc_write_model_atts','put_var znw')
call nc_check(nf90_put_var(ncid,      DNWVarID, wrf%dom(id)%dnw), &
              'nc_write_model_atts','put_var dnw')
if (wrf%dom(id)%hybrid_opt == VERT_HYBRID) then
   call nc_check(nf90_put_var(ncid,      C1HVarID, wrf%dom(id)%c1h), &
                 'nc_write_model_atts','put_var c1h')
   call nc_check(nf90_put_var(ncid,      C2HVarID, wrf%dom(id)%c2h), &
                 'nc_write_model_atts','put_var c2h')
endif

! defining horizontal
call nc_check(nf90_put_var(ncid,      mubVarID, wrf%dom(id)%mub), &
              'nc_write_model_atts','put_var mub')
call nc_check(nf90_put_var(ncid,      LonVarID, wrf%dom(id)%longitude), &
              'nc_write_model_atts','put_var longitude')
call nc_check(nf90_put_var(ncid,      LonuVarID, wrf%dom(id)%longitude_u), &
              'nc_write_model_atts','put_var longitude_u')
call nc_check(nf90_put_var(ncid,      LonvVarID, wrf%dom(id)%longitude_v), &
              'nc_write_model_atts','put_var longitude_v')
call nc_check(nf90_put_var(ncid,      LatVarID, wrf%dom(id)%latitude), &
              'nc_write_model_atts','put_var latitude')
call nc_check(nf90_put_var(ncid,      LatuVarID, wrf%dom(id)%latitude_u), &
              'nc_write_model_atts','put_var latitude_u')
call nc_check(nf90_put_var(ncid,      LatvVarID, wrf%dom(id)%latitude_v), &
              'nc_write_model_atts','put_var latitude_v')
call nc_check(nf90_put_var(ncid,     ilevVarID, (/ (i,i=1,wrf%dom(id)%bt) /)), &
              'nc_write_model_atts','put_var bt')
call nc_check(nf90_put_var(ncid,    XlandVarID, wrf%dom(id)%land), &
              'nc_write_model_atts','put_var land')
!   call nc_check(nf90_put_var(ncid,  MapFacMVarID(id), wrf%dom(id)%mapfac_m), &
!             'nc_write_model_atts','put_var mapfac_m')
!   call nc_check(nf90_put_var(ncid,  MapFacUVarID(id), wrf%dom(id)%mapfac_u), &
!             'nc_write_model_atts','put_var mapfac_u')
!   call nc_check(nf90_put_var(ncid,  MapFacVVarID(id), wrf%dom(id)%mapfac_v), &
!             'nc_write_model_atts','put_var mapfac_v')
call nc_check(nf90_put_var(ncid,      phbVarID, wrf%dom(id)%phb), &
              'nc_write_model_atts','put_var phb')
call nc_check(nf90_put_var(ncid,      hgtVarID, wrf%dom(id)%hgt), &
              'nc_write_model_atts','put_var hgt')


!-----------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-----------------------------------------------------------------

call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts


!#######################################################################

subroutine end_model()

deallocate(domain_id)

end subroutine end_model

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

!#######################################################################

subroutine toGrid_distrib(x, ens_size, j, dx, dxm)
!HK I don't think this is needed.
!  Transfer obs. x to grid j and calculate its
!  distance to grid j and j+1

  integer,  intent(in)  :: ens_size
  real(r8), intent(in)  :: x(ens_size)
  real(r8), intent(out) :: dx(ens_size), dxm(ens_size)
  integer,  intent(out) :: j(ens_size)

  j = int (x)

  dx = x - real (j)

  dxm= 1.0_r8 - dx

end subroutine toGrid_distrib

!#######################################################################

subroutine pres_to_zk(pres, mdl_v, n3, zk, lev0)

! Calculate the model level "zk" on half (mass) levels,
! corresponding to pressure "pres".

  integer,  intent(in)  :: n3
  real(r8), intent(in)  :: pres
  real(r8), intent(in)  :: mdl_v(0:n3)
  real(r8), intent(out) :: zk
  logical,  intent(out) :: lev0

  integer  :: k

  zk = missing_r8
  lev0 = .false.

  ! if out of range completely, return missing_r8 and lev0 false
  if (pres > mdl_v(0) .or. pres < mdl_v(n3)) return

  ! if above surface but below lowest sigma level, return the
  ! sigma value but set lev0 true
  if(pres <= mdl_v(0) .and. pres > mdl_v(1)) then
    lev0 = .true.
    if (log_vert_interp) then
       zk = (log(mdl_v(0)) - log(pres))/(log(mdl_v(0)) - log(mdl_v(1)))
    else
    zk = (mdl_v(0) - pres)/(mdl_v(0) - mdl_v(1))
    endif
    return
   endif

  ! find the 2 sigma levels the value is between and return that
  ! as a real number, including the fraction between the levels.
  do k = 1,n3-1
     if(pres <= mdl_v(k) .and. pres >= mdl_v(k+1)) then
        if (log_vert_interp) then
           zk = real(k) + (log(mdl_v(k)) - log(pres))/(log(mdl_v(k)) - log(mdl_v(k+1)))
        else
        zk = real(k) + (mdl_v(k) - pres)/(mdl_v(k) - mdl_v(k+1))
        endif
        exit
     endif
  enddo

end subroutine pres_to_zk

!#######################################################################

subroutine pres_to_zk_distrib(pres, mdl_v, n3, ens_size, zk, lev0)
! Calculate the model level "zk" on half (mass) levels,
! corresponding to pressure "pres".

  integer,  intent(in)  :: n3
  real(r8), intent(in)  :: pres
  integer,  intent(in)  :: ens_size
  real(r8), intent(in)  :: mdl_v(0:n3, ens_size)
  real(r8), intent(out) :: zk(ens_size)
  logical,  intent(out) :: lev0(ens_size)

  integer  :: k, e

  zk = missing_r8
  lev0 = .false.

  do e = 1, ens_size

     ! if out of range completely, return missing_r8 and lev0 false
     if (pres > mdl_v(0,e) .or. pres < mdl_v(n3,e)) then
        ! do nothing, keep zk(e) = missing_r8 and lev0(e) = false

     ! if above surface but below lowest sigma level, return the
     ! sigma value but set lev0 true
     elseif(pres <= mdl_v(0,e) .and. pres > mdl_v(1,e)) then
        lev0(e) = .true.
        if (log_vert_interp) then
           zk(e) = (log(mdl_v(0,e)) - log(pres))/(log(mdl_v(0,e)) - log(mdl_v(1,e)))
        else
           zk(e) = (mdl_v(0,e) - pres)/(mdl_v(0,e) - mdl_v(1,e))
        endif

     else

        ! find the 2 sigma levels the value is between and return that
        ! as a real number, including the fraction between the levels.
        do k = 1, n3-1
           if(pres <= mdl_v(k,e) .and. pres >= mdl_v(k+1,e)) then
              if (log_vert_interp) then
                 zk(e) = real(k) + (log(mdl_v(k,e)) - log(pres))/(log(mdl_v(k,e)) - log(mdl_v(k+1,e)))
              else
                 zk(e) = real(k) + (mdl_v(k,e) - pres)/(mdl_v(k,e) - mdl_v(k+1,e))
              endif
              exit
           endif
        enddo

     endif

  enddo

end subroutine pres_to_zk_distrib

!#######################################################################

subroutine height_to_zk(obs_v, mdl_v, n3, zk, lev0)

! Calculate the model level "zk" on half (mass) levels,
! corresponding to height "obs_v".

  real(r8), intent(in)  :: obs_v
  integer,  intent(in)  :: n3
  real(r8), intent(in)  :: mdl_v(0:n3)
  real(r8), intent(out) :: zk
  logical,  intent(out) :: lev0

  integer   :: k

  zk = missing_r8
  lev0 = .false.

  ! if out of range completely, return missing_r8 and lev0 false
  if (obs_v < mdl_v(0) .or. obs_v > mdl_v(n3)) return

  ! if above surface but below lowest 3-d height level, return the
  ! height value but set lev0 true
  if(obs_v >= mdl_v(0) .and. obs_v < mdl_v(1)) then
    lev0 = .true.
    zk = (mdl_v(0) - obs_v)/(mdl_v(0) - mdl_v(1))
    return
  endif

  ! find the 2 height levels the value is between and return that
  ! as a real number, including the fraction between the levels.
  do k = 1,n3-1
     if(obs_v >= mdl_v(k) .and. obs_v <= mdl_v(k+1)) then
        zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
        exit
     endif
  enddo

end subroutine height_to_zk

!#######################################################

subroutine get_model_pressure_profile_distrib(i,j,dx,dy,dxm,dym,n,id,v_p, state_handle, ens_size)

! Calculate the full model pressure profile on half (mass) levels,
! horizontally interpolated at the observation location.

integer,  intent(in)  :: i,j,n,id
real(r8), intent(in)  :: dx,dy,dxm,dym
integer, intent(in)   :: ens_size
real(r8), intent(out) :: v_p(0:n, ens_size)
type(ensemble_type), intent(in)  :: state_handle
integer e !< for ensemble loop

integer, dimension(2) :: ll, lr, ul, ur
integer(i8)           :: ill, ilr, iul, iur
integer               :: k, rc
real(r8), allocatable :: pres1(:), pres2(:), pres3(:), pres4(:)
logical  :: debug = .false.

!HK 
real(r8), allocatable :: x_ill(:), x_ilr(:), x_iul(:), x_iur(:)

allocate(pres1(ens_size), pres2(ens_size), pres3(ens_size), pres4(ens_size))
allocate(x_ill(ens_size), x_ilr(ens_size), x_iul(ens_size), x_iur(ens_size))

if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_t ) .and. &
     boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_t ) ) then

   call getCorners(i, j, id, wrf%dom(id)%type_t, ll, ul, lr, ur, rc )
   if ( rc .ne. 0 ) &
        print*, 'model_mod.f90 :: get_model_pressure_profile :: getCorners rc = ', rc


   do k=1,n
      pres1 = model_pressure_t_distrib(ll(1), ll(2), k,id,state_handle, ens_size)
      pres2 = model_pressure_t_distrib(lr(1), lr(2), k,id,state_handle, ens_size)
      pres3 = model_pressure_t_distrib(ul(1), ul(2), k,id,state_handle, ens_size)
      pres4 = model_pressure_t_distrib(ur(1), ur(2), k,id,state_handle, ens_size)

      v_p(k, :) = interp_4pressure_distrib(pres1, pres2, pres3, pres4, dx, dxm, dy, dym, ens_size)
   enddo


   if (debug) &
        print*, 'model_mod.f90 :: get_model_pressure_profile :: n, v_p() ', n, v_p(1:n, :)

   if ( wrf%dom(id)%type_ps >= 0 ) then

      ill = get_dart_vector_index(ll(1), ll(2), 1, domain_id(id), wrf%dom(id)%type_ps)
      ilr = get_dart_vector_index(lr(1), lr(2), 1, domain_id(id), wrf%dom(id)%type_ps)
      iul = get_dart_vector_index(ul(1), ul(2), 1, domain_id(id), wrf%dom(id)%type_ps)
      iur = get_dart_vector_index(ur(1), ur(2), 1, domain_id(id), wrf%dom(id)%type_ps)

      x_ill = get_state(ill, state_handle)
      x_ilr = get_state(ilr, state_handle)
      x_iul = get_state(iul, state_handle)
      x_iur = get_state(iur, state_handle)

      ! I'm not quite sure where this comes from, but I will trust them on it....
      ! Do you have to do this per ensemble?
      !>@todo This is messy
      do e = 1,ens_size

         if ( x_ill(e) /= 0.0_r8 .and. x_ilr(e) /= 0.0_r8 .and. x_iul(e) /= 0.0_r8 .and. &
              x_iur(e) /= 0.0_r8 ) then

            v_p(0,e:e) = interp_4pressure_distrib(x_ill(e:e), x_ilr(e:e), x_iul(e:e), x_iur(e:e), dx, dxm, dy, dym, 1)

         else

            ! HK I think this is a bug, you are just  going to grab the first copy - is this fixed?
            ! in each iteration of the loop
            call error_handler(E_WARN, 'model_mod.f90 check for correctness', 'Helen')
            pres1 = model_pressure_t_distrib(ll(1), ll(2), 2,id,state_handle, ens_size)
            pres2 = model_pressure_t_distrib(lr(1), lr(2), 2,id,state_handle, ens_size)
            pres3 = model_pressure_t_distrib(ul(1), ul(2), 2,id,state_handle, ens_size)
            pres4 = model_pressure_t_distrib(ur(1), ur(2), 2,id,state_handle, ens_size)

            v_p(0,e:e) = interp_4pressure_distrib(pres1(e:e), pres2(e:e), pres3(e:e), pres4(e:e), dx, dxm, dy, dym, 1, &
                  extrapolate=.true., edgep=v_p(1,e))

         endif

      enddo

   else

      pres1 = model_pressure_t_distrib(ll(1), ll(2), 2,id,state_handle, ens_size)
      pres2 = model_pressure_t_distrib(lr(1), lr(2), 2,id,state_handle, ens_size)
      pres3 = model_pressure_t_distrib(ul(1), ul(2), 2,id,state_handle, ens_size)
      pres4 = model_pressure_t_distrib(ur(1), ur(2), 2,id,state_handle, ens_size)

      v_p(0,:) = interp_4pressure_distrib(pres1, pres2, pres3, pres4, dx, dxm, dy, dym, ens_size, &
              extrapolate=.true., edgep=v_p(1,:))

   endif

   if (debug) &
        print*, 'model_mod.f90 :: get_model_pressure_profile :: v_p(0) ', v_p(0, :)
else
   v_p(:,:) = missing_r8

endif

deallocate(pres1, pres2, pres3, pres4, x_ill, x_ilr, x_iul, x_iur)

end subroutine get_model_pressure_profile_distrib

!#######################################################
!> Only for the mean value.
!> Used in get_state_meta_data for the vertical conversion
function model_pressure_distrib(i, j, k, id, var_type, state_handle)

! Calculate the pressure at grid point (i,j,k), domain id.
! The grid is defined according to var_type.

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: i,j,k,id,var_type
real(r8)                        :: model_pressure_distrib

integer  :: off, type_x
real(r8) :: pres1, pres2
integer  :: ens_size

model_pressure_distrib = missing_r8

! since we are only working with the mean copy
ens_size = 1

! If W-grid (on ZNW levels), then we need to average in vertical, unless
!   we are at the upper or lower boundary in which case we will extrapolate.
if( (var_type == wrf%dom(id)%type_w) .or. (var_type == wrf%dom(id)%type_gz) ) then

   if (var_type == wrf%dom(id)%type_w) then
      type_x = wrf%dom(id)%type_w
   else
      type_x = wrf%dom(id)%type_gz
   endif

   if( k == 1 ) then

      pres1 = scalar(model_pressure_t_distrib(i, j, k,  id, state_handle, ens_size))
      pres2 = scalar(model_pressure_t_distrib(i, j, k+1,id, state_handle, ens_size))
      model_pressure_distrib = interp_pressure(pres1, pres2, extrapolate=.true.)

   elseif( k == wrf%dom(id)%var_size(3,type_x) ) then

      pres1 = scalar(model_pressure_t_distrib(i,j,k-1,id, state_handle, ens_size))
      pres2 = scalar(model_pressure_t_distrib(i,j,k-2,id, state_handle, ens_size))
      model_pressure_distrib = interp_pressure(pres1, pres2, extrapolate=.true.)

   else

      pres1 = scalar(model_pressure_t_distrib(i, j, k,  id, state_handle, ens_size))
      pres2 = scalar(model_pressure_t_distrib(i, j, k-1,id, state_handle, ens_size))
      model_pressure_distrib = interp_pressure(pres1, pres2)

   endif

! If U-grid, then pressure is defined between U points, so average --
!   averaging depends on longitude periodicity
elseif( var_type == wrf%dom(id)%type_u ) then

   if( i == wrf%dom(id)%var_size(1,wrf%dom(id)%type_u) ) then

      ! Check to see if periodic in longitude
      if ( wrf%dom(id)%periodic_x ) then

         ! We are at seam in longitude, take first and last M-grid points
         pres1 = scalar(model_pressure_t_distrib(i-1,j,k,id, state_handle, ens_size))
         pres2 = scalar(model_pressure_t_distrib(1,  j,k,id, state_handle, ens_size))
         model_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)
         
      else

         ! If not periodic, then try extrapolating
         pres1 = scalar(model_pressure_t_distrib(i-1,j,k,id, state_handle, ens_size))
         pres2 = scalar(model_pressure_t_distrib(i-2,j,k,id, state_handle, ens_size))
         model_pressure_distrib = interp_pressure(pres1, pres2, extrapolate=.true., vertical=.false.)

      endif

   elseif( i == 1 ) then

      ! Check to see if periodic in longitude
      if ( wrf%dom(id)%periodic_x ) then

         ! We are at seam in longitude, take first and last M-grid points
         pres1 = scalar(model_pressure_t_distrib(i,             j,k,id, state_handle, ens_size))
         pres2 = scalar(model_pressure_t_distrib(wrf%dom(id)%we,j,k,id, state_handle, ens_size))
         model_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)
         
      else

         ! If not periodic, then try extrapolating
         pres1 = scalar(model_pressure_t_distrib(i,  j,k,id, state_handle, ens_size))
         pres2 = scalar(model_pressure_t_distrib(i+1,j,k,id, state_handle, ens_size))
         model_pressure_distrib = interp_pressure(pres1, pres2, extrapolate=.true., vertical=.false.)

      endif

   else

      pres1 = scalar(model_pressure_t_distrib(i,  j,k,id, state_handle, ens_size))
      pres2 = scalar(model_pressure_t_distrib(i-1,j,k,id, state_handle, ens_size))
      model_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)

   endif

! If V-grid, then pressure is defined between V points, so average --
!   averaging depends on polar periodicity
elseif( var_type == wrf%dom(id)%type_v ) then

   if( j == wrf%dom(id)%var_size(2,wrf%dom(id)%type_v) ) then

      ! Check to see if periodic in latitude (polar)
      if ( wrf%dom(id)%polar ) then

         ! The upper corner is 180 degrees of longitude away
         off = i + wrf%dom(id)%we/2
         if ( off > wrf%dom(id)%we ) off = off - wrf%dom(id)%we

         pres1 = scalar(model_pressure_t_distrib(off,j-1,k,id, state_handle, ens_size))
         pres2 = scalar(model_pressure_t_distrib(i  ,j-1,k,id, state_handle, ens_size))
         model_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)

      ! If not periodic, then try extrapolating
      else

         pres1 = scalar(model_pressure_t_distrib(i,j-1,k,id, state_handle, ens_size))
         pres2 = scalar(model_pressure_t_distrib(i,j-2,k,id, state_handle, ens_size))
         model_pressure_distrib = interp_pressure(pres1, pres2, extrapolate=.true., vertical=.false.)

      endif

   elseif( j == 1 ) then

      ! Check to see if periodic in latitude (polar)
      if ( wrf%dom(id)%polar ) then

         ! The lower corner is 180 degrees of longitude away
         off = i + wrf%dom(id)%we/2
         if ( off > wrf%dom(id)%we ) off = off - wrf%dom(id)%we

         pres1 = scalar(model_pressure_t_distrib(off,j,k,id, state_handle, ens_size))
         pres2 = scalar(model_pressure_t_distrib(i,  j,k,id, state_handle, ens_size))
         model_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)

      ! If not periodic, then try extrapolating
      else

         pres1 = scalar(model_pressure_t_distrib(i,j,  k,id, state_handle, ens_size))
         pres2 = scalar(model_pressure_t_distrib(i,j+1,k,id, state_handle, ens_size))
         model_pressure_distrib = interp_pressure(pres1, pres2, extrapolate=.true., vertical=.false.)

      endif

   else

      pres1 = scalar(model_pressure_t_distrib(i,j,  k,id, state_handle, ens_size))
      pres2 = scalar(model_pressure_t_distrib(i,j-1,k,id, state_handle, ens_size))
      model_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)

   endif

elseif( var_type == wrf%dom(id)%type_mu    .or. var_type == wrf%dom(id)%type_tslb .or. &
        var_type == wrf%dom(id)%type_ps    .or. var_type == wrf%dom(id)%type_u10  .or. &
        var_type == wrf%dom(id)%type_v10   .or. var_type == wrf%dom(id)%type_t2   .or. &
        var_type == wrf%dom(id)%type_th2   .or.                                        &
        var_type == wrf%dom(id)%type_q2    .or. var_type == wrf%dom(id)%type_tsk  .or. &
        var_type == wrf%dom(id)%type_smois .or. var_type == wrf%dom(id)%type_sh2o) then

   model_pressure_distrib = model_pressure_s_distrib(i,j,id, state_handle)
    
else

   pres1 = scalar(model_pressure_t_distrib(i,j,k,id, state_handle, ens_size))
   model_pressure_distrib = pres1

endif

end function model_pressure_distrib

!#######################################################

function model_surface_pressure_distrib(i, j, id, var_type, state_handle)

! Calculate the surface pressure at grid point (i,j), domain id.
! The grid is defined according to var_type.

type(ensemble_type), intent(in) :: state_handle
integer,            intent(in)  :: i,j,id,var_type
real(r8)              :: model_surface_pressure_distrib

integer  :: off
real(r8) :: pres1, pres2

model_surface_pressure_distrib = missing_r8


! If U-grid, then pressure is defined between U points, so average --
!   averaging depends on longitude periodicity
if( var_type == wrf%dom(id)%type_u ) then

   if( i == wrf%dom(id)%var_size(1,wrf%dom(id)%type_u) ) then

      ! Check to see if periodic in longitude
      if ( wrf%dom(id)%periodic_x ) then

         ! We are at seam in longitude, take first and last M-grid points
         pres1 = model_pressure_s_distrib(i-1,j,id, state_handle)
         pres2 = model_pressure_s_distrib(1,  j,id, state_handle)
         model_surface_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)
         
      else

         ! If not periodic, then try extrapolating
         pres1 = model_pressure_s_distrib(i-1,j,id, state_handle)
         pres2 = model_pressure_s_distrib(i-2,j,id, state_handle)
         model_surface_pressure_distrib = interp_pressure(pres1, pres2, extrapolate=.true., vertical=.false.)

      endif

   elseif( i == 1 ) then

      ! Check to see if periodic in longitude
      if ( wrf%dom(id)%periodic_x ) then

         ! We are at seam in longitude, take first and last M-grid points
         pres1 = model_pressure_s_distrib(i,             j,id, state_handle)
         pres2 = model_pressure_s_distrib(wrf%dom(id)%we,j,id, state_handle)
         model_surface_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)
         
      else

         ! If not periodic, then try extrapolating
         pres1 = model_pressure_s_distrib(i,  j,id, state_handle)
         pres2 = model_pressure_s_distrib(i+1,j,id, state_handle)
         model_surface_pressure_distrib = interp_pressure(pres1, pres2, extrapolate=.true., vertical=.false.)

      endif

   else

      pres1 = model_pressure_s_distrib(i,  j,id, state_handle)
      pres2 = model_pressure_s_distrib(i-1,j,id, state_handle)
      model_surface_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)

   endif

! If V-grid, then pressure is defined between V points, so average --
!   averaging depends on polar periodicity
elseif( var_type == wrf%dom(id)%type_v ) then

   if( j == wrf%dom(id)%var_size(2,wrf%dom(id)%type_v) ) then

      ! Check to see if periodic in latitude (polar)
      if ( wrf%dom(id)%polar ) then

         ! The upper corner is 180 degrees of longitude away
         off = i + wrf%dom(id)%we/2
         if ( off > wrf%dom(id)%we ) off = off - wrf%dom(id)%we

         pres1 = model_pressure_s_distrib(off,j-1,id, state_handle)
         pres2 = model_pressure_s_distrib(i  ,j-1,id, state_handle)
         model_surface_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)

      ! If not periodic, then try extrapolating
      else

         pres1 = model_pressure_s_distrib(i,j-1,id, state_handle)
         pres2 = model_pressure_s_distrib(i,j-2,id, state_handle)
         model_surface_pressure_distrib = interp_pressure(pres1, pres2, extrapolate=.true., vertical=.false.)

      endif

   elseif( j == 1 ) then

      ! Check to see if periodic in latitude (polar)
      if ( wrf%dom(id)%polar ) then

         ! The lower corner is 180 degrees of longitude away
         off = i + wrf%dom(id)%we/2
         if ( off > wrf%dom(id)%we ) off = off - wrf%dom(id)%we

         pres1 = model_pressure_s_distrib(off,j,id, state_handle)
         pres2 = model_pressure_s_distrib(i,  j,id, state_handle)
         model_surface_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)

      ! If not periodic, then try extrapolating
      else

         pres1 = model_pressure_s_distrib(i,j,  id, state_handle)
         pres2 = model_pressure_s_distrib(i,j+1,id, state_handle)
         model_surface_pressure_distrib = interp_pressure(pres1, pres2, extrapolate=.true., vertical=.false.)

      endif

   else

      pres1 = model_pressure_s_distrib(i,j,  id, state_handle)
      pres2 = model_pressure_s_distrib(i,j-1,id, state_handle)
      model_surface_pressure_distrib = interp_pressure(pres1, pres2, vertical=.false.)

   endif

else

   model_surface_pressure_distrib = model_pressure_s_distrib(i,j,id, state_handle)

endif

end function model_surface_pressure_distrib

!#######################################################

function model_pressure_t_distrib(i,j,k,id,state_handle, ens_size)

! Calculate total pressure on mass point (half (mass) levels, T-point).

integer,             intent(in) :: ens_size
integer,             intent(in) :: i,j,k,id
type(ensemble_type), intent(in) :: state_handle
real(r8) :: model_pressure_t_distrib(ens_size)

real (kind=r8), PARAMETER    :: rd_over_rv = gas_constant / gas_constant_v
real (kind=r8), PARAMETER    :: cpovcv = 1.4_r8        ! cp / (cp - gas_constant)

integer(i8) :: iqv !< I think this is i for index
integer(i8) :: it !< change to array
real(r8) :: qvf1(ens_size),rho(ens_size), x_iqv(ens_size), x_it(ens_size)

model_pressure_t_distrib = missing_r8

! Adapted the code from WRF module_big_step_utilities_em.F ----
!         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)

! Simplification: alb*mub = (phb(i,j,k+1) - phb(i,j,k))/dnw(k)

if (wrf%dom(id)%type_qv < 0 .or. wrf%dom(id)%type_t < 0) then
  call error_handler(E_ERR, 'model_pressure_t:', &
      'BOTH QVAPOR and T must be in state vector to compute total pressure', &
       source, revision, revdate)
endif

iqv = get_dart_vector_index(i,j,k, domain_id(id), wrf%dom(id)%type_qv)
it  = get_dart_vector_index(i,j,k, domain_id(id), wrf%dom(id)%type_t)

x_iqv = get_state(iqv, state_handle)
x_it  = get_state(it, state_handle)

qvf1(:) = 1.0_r8 + x_iqv(:) / rd_over_rv
!print*, 'qvf1 ', qvf1

rho(:) = model_rho_t_distrib(i,j,k,id,state_handle, ens_size)
!print*, 'rho ', rho

! .. total pressure:
model_pressure_t_distrib = ps0 * ( (gas_constant*(ts0+x_it)*qvf1) / &
     (ps0/rho) )**cpovcv

end function model_pressure_t_distrib

!#######################################################

function model_pressure_s_distrib(i, j, id, state_handle)

! compute pressure at surface at mass point

integer,             intent(in) :: i,j,id
type(ensemble_type), intent(in) :: state_handle

real(r8)              :: model_pressure_s_distrib

integer(i8) :: ips, imu
real(r8)    :: x_imu, x_ips


! make sure one of these is good.
if ( wrf%dom(id)%type_mu < 0 .and. wrf%dom(id)%type_ps < 0 ) then
  call error_handler(E_ERR, 'model_pressure_s:', &
      'One of MU or PSFC must be in state vector to compute surface pressure', &
       source, revision, revdate)
endif

if ( wrf%dom(id)%type_ps >= 0 ) then
   ips = get_dart_vector_index(i,j,1, domain_id(id), wrf%dom(id)%type_ps)
   x_ips = scalar(get_state(ips, state_handle))
   model_pressure_s_distrib = x_ips

else
   imu = get_dart_vector_index(i,j,1, domain_id(id), wrf%dom(id)%type_mu)
   x_imu = minval(get_state(imu, state_handle))
   model_pressure_s_distrib = wrf%dom(id)%p_top + wrf%dom(id)%mub(i,j) + x_imu

endif

end function model_pressure_s_distrib

!#######################################################

function interp_pressure(p1, p2, extrapolate, vertical)
 
! interpolate halfway between pressures 1 and 2 in log units.
! if extrapolate is true, extrapolate where 1 is the edge and
! 2 is the inner value, going 1/2 grid cell out.

real(r8), intent(in)           :: p1, p2
logical,  intent(in), optional :: extrapolate
logical,  intent(in), optional :: vertical
real(r8)                       :: interp_pressure

logical  :: do_interp
logical  :: is_vert
real(r8) :: intermediate

! default is to do interpolation; only extrapolate if the optional
! arg is specified and if it is true.
do_interp = .true.
if (present(extrapolate)) then
   if (extrapolate) do_interp = .false.
endif

! if vert is specified and is false, check log_horz_interpM instead
! of log_vert_interp to decide log vs linear interpolation for the
! Midpoint value.  default is to do vertical interpolation.
is_vert = .true.
if (present(vertical)) then
   is_vert = vertical
endif

! once we like the results, remove the log_vert_interp test.
if (do_interp) then
   if ((      is_vert .and. log_vert_interp )  .or. &
       (.not. is_vert .and. log_horz_interpM)) then
      interp_pressure = exp((log(p1) + log(p2))/2.0_r8)
   else
      interp_pressure = (p1 + p2)/2.0_r8
   endif
else
   if ((      is_vert .and. log_vert_interp )  .or. &
       (.not. is_vert .and. log_horz_interpM)) then
      intermediate = (3.0_r8*log(p1) - log(p2))/2.0_r8
      if (intermediate <= 0.0_r8) then
         interp_pressure = p1
      else
         interp_pressure = exp(intermediate)
      endif
   else
      interp_pressure = (3.0_r8*p1 - p2)/2.0_r8
   endif
endif

end function interp_pressure

!#######################################################

function interp_4pressure_distrib(p1, p2, p3, p4, dx, dxm, dy, dym, ens_size, extrapolate, edgep)
 
! given 4 corners of a quad, where the p1, p2, p3 and p4 points are
! respectively:  lower left, lower right, upper left, upper right
! and dx is the distance in x, dxm is 1.0-dx, dy is distance in y
! and dym is 1.0-dy, interpolate the pressure while converted to log.
! if extrapolate is true, extrapolate where edgep is the edge pressure
! and the 4 points and dx/dy give the location of the inner point.

integer, intent(in)                :: ens_size
real(r8), intent(in)               :: p1(ens_size), p2(ens_size), p3(ens_size), p4(ens_size)
real(r8), intent(in)               :: dx, dxm, dy, dym
logical,  intent(in), optional     :: extrapolate
real(r8), intent(in), optional     :: edgep(ens_size)
real(r8)                           :: interp_4pressure_distrib(ens_size)

logical  :: do_interp
real(r8) :: intermediate(ens_size)
real(r8) :: l1(ens_size), l2(ens_size), l3(ens_size), l4(ens_size)

!HK
integer :: i

! default is to do interpolation; only extrapolate if the optional
! arg is specified and if it is true.  for extrapolation 'edgep' is
! required; it is unused for interpolation.
do_interp = .true.
if (present(extrapolate)) then
   if (extrapolate) do_interp = .false.
endif

if (.not. do_interp .and. .not. present(edgep)) then
  call error_handler(E_ERR, 'interp_4pressure:', &
      'edgep must be specified for extrapolation.  internal error.', &
       source, revision, revdate)
endif

if (log_horz_interpQ) then
   l1 = log(p1)
   l2 = log(p2)
   l3 = log(p3)
   l4 = log(p4)
endif


! once we like the results, remove the log_horz_interpQ test.
if (do_interp) then
   if (log_horz_interpQ) then
      interp_4pressure_distrib = exp(dym*( dxm*l1 + dx*l2 ) + dy*( dxm*l3 + dx*l4 ))
   else
      interp_4pressure_distrib = dym*( dxm*p1 + dx*p2 ) + dy*( dxm*p3 + dx*p4 )
   endif
else
   if (log_horz_interpQ) then
      intermediate = (3.0_r8*log(edgep) - &
                 dym*( dxm*l1 + dx*l2 ) - dy*( dxm*l3 + dx*l4 ))/2.0_r8

      do i = 1, size(intermediate)
         if (intermediate(i) <= 0.0_r8) then
            !>@todo is this correct?
            interp_4pressure_distrib(i) = edgep(i)
         else
            interp_4pressure_distrib(i) = exp(intermediate(i))
         endif
      enddo
   else
      interp_4pressure_distrib = (3.0_r8*edgep - &
                 dym*( dxm*p1 + dx*p2 ) - dy*( dxm*p3 + dx*p4 ))/2.0_r8
   endif
endif

end function interp_4pressure_distrib

!#######################################################
!> Interpolate function
function interp_4pressure(p1, p2, p3, p4, dx, dxm, dy, dym, extrapolate, edgep)
 
! given 4 corners of a quad, where the p1, p2, p3 and p4 points are
! respectively:  lower left, lower right, upper left, upper right
! and dx is the distance in x, dxm is 1.0-dx, dy is distance in y
! and dym is 1.0-dy, interpolate the pressure while converted to log.
! if extrapolate is true, extrapolate where edgep is the edge pressure
! and the 4 points and dx/dy give the location of the inner point.

real(r8), intent(in)               :: p1, p2, p3, p4
real(r8), intent(in)               :: dx, dxm, dy, dym
logical,  intent(in), optional     :: extrapolate
real(r8), intent(in), optional     :: edgep
real(r8)                           :: interp_4pressure

logical  :: do_interp
real(r8) :: intermediate
real(r8) :: l1, l2, l3, l4

! default is to do interpolation; only extrapolate if the optional
! arg is specified and if it is true.  for extrapolation 'edgep' is
! required; it is unused for interpolation.
do_interp = .true.
if (present(extrapolate)) then
   if (extrapolate) do_interp = .false.
endif

if (.not. do_interp .and. .not. present(edgep)) then
  call error_handler(E_ERR, 'interp_4pressure:', &
      'edgep must be specified for extrapolation.  internal error.', &
       source, revision, revdate)
endif

if (log_horz_interpQ) then
   l1 = log(p1)
   l2 = log(p2)
   l3 = log(p3)
   l4 = log(p4)
endif


! once we like the results, remove the log_horz_interpQ test.
if (do_interp) then
   if (log_horz_interpQ) then
      interp_4pressure = exp(dym*( dxm*l1 + dx*l2 ) + dy*( dxm*l3 + dx*l4 ))
   else
      interp_4pressure = dym*( dxm*p1 + dx*p2 ) + dy*( dxm*p3 + dx*p4 )
   endif
else
   if (log_horz_interpQ) then
      intermediate = (3.0_r8*log(edgep) - &
                 dym*( dxm*l1 + dx*l2 ) - dy*( dxm*l3 + dx*l4 ))/2.0_r8
      if (intermediate <= 0.0_r8) then
         interp_4pressure = edgep
      else
         interp_4pressure = exp(intermediate)
      endif
   else
      interp_4pressure = (3.0_r8*edgep - &
                 dym*( dxm*p1 + dx*p2 ) - dy*( dxm*p3 + dx*p4 ))/2.0_r8
   endif
endif

end function interp_4pressure

!#######################################################

function model_rho_t_distrib(i,j,k,id,state_handle, ens_size)

! Calculate the total density on mass point (half (mass) levels, T-point).

integer,             intent(in)  :: ens_size
integer,             intent(in)  :: i,j,k,id
type(ensemble_type), intent(in)  :: state_handle
real(r8) :: model_rho_t_distrib(ens_size)

integer(i8) :: imu,iph,iphp1
real(r8)    :: ph_e(ens_size), x_imu(ens_size), x_iph(ens_size), x_iphp1(ens_size)

model_rho_t_distrib(:) = missing_r8

! Adapted the code from WRF module_big_step_utilities_em.F ----
!         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)

! Simplification: alb*mub = (phb(i,j,k+1) - phb(i,j,k))/dnw(k)

if (wrf%dom(id)%type_mu < 0 .or. wrf%dom(id)%type_gz < 0) then
  call error_handler(E_ERR, 'model_rho_t:', &
      'BOTH MU and PH must be in state vector to compute total density', &
       source, revision, revdate)
endif

imu   = get_dart_vector_index(i,j,1,   domain_id(id), wrf%dom(id)%type_mu)
iph   = get_dart_vector_index(i,j,k,   domain_id(id), wrf%dom(id)%type_gz)
iphp1 = get_dart_vector_index(i,j,k+1, domain_id(id), wrf%dom(id)%type_gz)

x_imu = get_state(imu, state_handle)
x_iph = get_state(iph, state_handle)
x_iphp1 = get_state(iphp1, state_handle)

ph_e = ( (x_iphp1 + wrf%dom(id)%phb(i,j,k+1)) &
       - (x_iph   + wrf%dom(id)%phb(i,j,k  )) ) / wrf%dom(id)%dnw(k)

! now calculate rho = - mu / dphi/deta

if (wrf%dom(id)%hybrid_opt == VERT_HYBRID) then
   model_rho_t_distrib(:) = - (wrf%dom(id)%c1h(k)*(wrf%dom(id)%mub(i,j)+x_imu) + wrf%dom(id)%c2h(k)) / ph_e
else
   model_rho_t_distrib(:) = - (wrf%dom(id)%mub(i,j)+x_imu) / ph_e
endif

end function model_rho_t_distrib

!#######################################################

subroutine get_model_height_profile_distrib(i,j,dx,dy,dxm,dym,n,id,v_h, state_handle, ens_size)

! Calculate the model height profile on half (mass) levels,
! horizontally interpolated at the observation location.
! This routine used to compute geopotential heights; it now
! computes geometric heights.

integer,  intent(in)  :: i,j,n,id
real(r8), intent(in)  :: dx,dy,dxm,dym
integer,  intent(in)  :: ens_size
real(r8), intent(out) :: v_h(0:n, ens_size)
type(ensemble_type), intent(in)  :: state_handle
integer e !< for ensemble loop

real(r8)              :: fll(n+1, ens_size), geop(ens_size), lat(ens_size)
integer(i8)           :: ill, iul, ilr, iur
integer               :: k, rc
integer, dimension(2) :: ll, lr, ul, ur
logical   :: debug = .false.

real(r8), allocatable :: x_ill(:), x_ilr(:), x_iul(:), x_iur(:)

allocate(x_ill(ens_size), x_ilr(ens_size), x_iul(ens_size), x_iur(ens_size))

if (wrf%dom(id)%type_gz < 0) then
  call error_handler(E_ERR, 'get_model_height_profile:', &
      'PH must be in state vector to compute height profile', &
       source, revision, revdate)
endif

if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf%dom(id)%type_gz ) .and. &
     boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf%dom(id)%type_gz ) ) then

   call getCorners(i, j, id, wrf%dom(id)%type_gz, ll, ul, lr, ur, rc )
   if ( rc .ne. 0 ) &
        print*, 'model_mod.f90 :: get_model_height_profile :: getCorners rc = ', rc

   do k = 1, wrf%dom(id)%var_size(3,wrf%dom(id)%type_gz)

      ill = get_dart_vector_index(ll(1), ll(2), k, domain_id(id), wrf%dom(id)%type_gz)
      iul = get_dart_vector_index(ul(1), ul(2), k, domain_id(id), wrf%dom(id)%type_gz)
      ilr = get_dart_vector_index(lr(1), lr(2), k, domain_id(id), wrf%dom(id)%type_gz)
      iur = get_dart_vector_index(ur(1), ur(2), k, domain_id(id), wrf%dom(id)%type_gz)

      x_ill = get_state(ill, state_handle)
      x_ilr = get_state(ilr, state_handle)
      x_iul = get_state(iul, state_handle)
      x_iur = get_state(iur, state_handle)

      geop(:) = ( dym*( dxm*( wrf%dom(id)%phb(ll(1),ll(2),k) + x_ill ) + &
                      dx*( wrf%dom(id)%phb(lr(1),lr(2),k) + x_ilr ) ) + &
                dy*( dxm*( wrf%dom(id)%phb(ul(1),ul(2),k) + x_iul ) + &
                      dx*( wrf%dom(id)%phb(ur(1),ur(2),k) + x_iur ) ) )/gravity

      lat(:) = ( wrf%dom(id)%latitude(ll(1),ll(2)) + &
              wrf%dom(id)%latitude(lr(1),lr(2)) + &
              wrf%dom(id)%latitude(ul(1),ul(2)) + &
              wrf%dom(id)%latitude(ur(1),ur(2)) ) / 4.0_r8

      do e = 1, ens_size
         fll(k, e) = compute_geometric_height(geop(e), lat(e))
      enddo

   end do

   do k=1,n
      v_h(k, :) = 0.5_r8*(fll(k, :) + fll(k+1, :) )
   end do

   v_h(0, :) = dym*( dxm*wrf%dom(id)%hgt(ll(1), ll(2)) + &
                   dx*wrf%dom(id)%hgt(lr(1), lr(2)) ) + &
             dy*( dxm*wrf%dom(id)%hgt(ul(1), ul(2)) + &
                   dx*wrf%dom(id)%hgt(ur(1), ur(2)) )

   if (debug) &
        print*, 'model_mod.f90 :: get_model_height_profile :: n, v_h() ', n, v_h(1:n, :)

   if (debug) &
        print*, 'model_mod.f90 :: get_model_height_profile :: v_h(0) ', v_h(0, :)
 
! If the boundsCheck functions return an unsatisfactory integer index, then set
!   fld as missing data
else

   print*,'Not able the get height_profile'
   print*,i,j,dx,dy,dxm,dym,n,id,wrf%dom(id)%var_size(1,wrf%dom(id)%type_gz), &
        wrf%dom(id)%var_size(2,wrf%dom(id)%type_gz)

   v_h(:, :) =  missing_r8

endif

end subroutine get_model_height_profile_distrib

!#######################################################

function model_height_distrib(i,j,k,id,var_type, state_handle)

! This routine used to compute geopotential heights; it now
! computes geometric heights.

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: i,j,k,id,var_type
real(r8)                        :: model_height_distrib

integer(i8)  :: i1, i2, i3, i4
integer      :: off
real(r8)     :: x_i1(1), x_i2(1), x_i3(1), x_i4(1)
real(r8)     :: geop, lat

model_height_distrib = missing_r8

if (wrf%dom(id)%type_gz < 0) then
  call error_handler(E_ERR, 'model_height:', &
      'PH must be in state vector to compute total pressure', &
       source, revision, revdate)
endif

! If W-grid (on ZNW levels), then we are fine because it is native to GZ
if( (var_type == wrf%dom(id)%type_w) .or. (var_type == wrf%dom(id)%type_gz) ) then

   i1 = get_dart_vector_index(i,j,k, domain_id(id),wrf%dom(id)%type_gz)
   x_i1 = get_state(i1, state_handle)

   geop = minval((wrf%dom(id)%phb(i,j,k)+x_i1)/gravity)
   model_height_distrib = compute_geometric_height(geop, wrf%dom(id)%latitude(i, j))

! If U-grid, then height is defined between U points, both in horizontal 
!   and in vertical, so average -- averaging depends on longitude periodicity
elseif( var_type == wrf%dom(id)%type_u ) then

   if( i == wrf%dom(id)%var_size(1,wrf%dom(id)%type_u) ) then

      ! Check to see if periodic in longitude
      if ( wrf%dom(id)%periodic_x ) then

         ! We are at the seam in longitude, so take first and last mass points
         i1 = get_dart_vector_index(i-1,j,k  , domain_id(id),wrf%dom(id)%type_gz)
         i2 = get_dart_vector_index(i-1,j,k+1, domain_id(id),wrf%dom(id)%type_gz)
         i3 = get_dart_vector_index(1,  j,k  , domain_id(id),wrf%dom(id)%type_gz)
         i4 = get_dart_vector_index(1,  j,k+1, domain_id(id),wrf%dom(id)%type_gz)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)


         geop = minval(( (wrf%dom(id)%phb(i-1,j,k  ) + x_i1) &
                 +(wrf%dom(id)%phb(i-1,j,k+1) + x_i2) &
                 +(wrf%dom(id)%phb(1  ,j,k  ) + x_i3) &
                 +(wrf%dom(id)%phb(1  ,j,k+1) + x_i4) )/(4.0_r8*gravity))
         
         lat = ( wrf%dom(id)%latitude(i-1,j)  &
                +wrf%dom(id)%latitude(i-1,j)  &
                +wrf%dom(id)%latitude(1  ,j)  &
                +wrf%dom(id)%latitude(1  ,j) ) / 4.0_r8

         model_height_distrib = compute_geometric_height(geop, lat)
         
      else

         ! If not periodic, then try extrapolating
         i1 = get_dart_vector_index(i-1,j,k  , domain_id(id),wrf%dom(id)%type_gz)
         i2 = get_dart_vector_index(i-1,j,k+1, domain_id(id),wrf%dom(id)%type_gz)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i1 -1, state_handle)
         x_i4 = get_state(i2 -1, state_handle)


         geop = minval(( 3.0_r8*(wrf%dom(id)%phb(i-1,j,k  )+x_i1) &
                 +3.0_r8*(wrf%dom(id)%phb(i-1,j,k+1)+x_i2) &
                        -(wrf%dom(id)%phb(i-2,j,k  )+x_i3) &
                        -(wrf%dom(id)%phb(i-2,j,k+1)+x_i4) )/(4.0_r8*gravity))

         lat = ( 3.0_r8*wrf%dom(id)%latitude(i-1,j)  &
                +3.0_r8*wrf%dom(id)%latitude(i-1,j)  &
                       -wrf%dom(id)%latitude(i-2,j)  &
                       -wrf%dom(id)%latitude(i-2,j)) / 4.0_r8

         model_height_distrib = compute_geometric_height(geop, lat)

      endif

   elseif( i == 1 ) then

      ! Check to see if periodic in longitude
      if ( wrf%dom(id)%periodic_x ) then

         ! We are at the seam in longitude, so take first and last mass points
         off = wrf%dom(id)%we
         i1 = get_dart_vector_index(i  ,j,k  ,domain_id(id),wrf%dom(id)%type_gz)
         i2 = get_dart_vector_index(i  ,j,k+1,domain_id(id),wrf%dom(id)%type_gz)
         i3 = get_dart_vector_index(off,j,k  ,domain_id(id),wrf%dom(id)%type_gz)
         i4 = get_dart_vector_index(off,j,k+1,domain_id(id),wrf%dom(id)%type_gz)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)

         geop = minval(( (wrf%dom(id)%phb(i  ,j,k  ) + x_i1) &
                 +(wrf%dom(id)%phb(i  ,j,k+1) + x_i2) &
                 +(wrf%dom(id)%phb(off,j,k  ) + x_i3) &
                 +(wrf%dom(id)%phb(off,j,k+1) + x_i4) )/(4.0_r8*gravity))
         
         lat = ( wrf%dom(id)%latitude(i  ,j)  &
                +wrf%dom(id)%latitude(i  ,j)  &
                +wrf%dom(id)%latitude(off,j)  &
                +wrf%dom(id)%latitude(off,j)) / 4.0_r8

         model_height_distrib = compute_geometric_height(geop, lat)

      else

         ! If not periodic, then try extrapolating
         i1 = get_dart_vector_index(i,j,k  ,domain_id(id),wrf%dom(id)%type_gz)
         i2 = get_dart_vector_index(i,j,k+1,domain_id(id),wrf%dom(id)%type_gz)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i1 +1, state_handle)
         x_i4 = get_state(i2 +1, state_handle)


         geop = minval(( 3.0_r8*(wrf%dom(id)%phb(i  ,j,k  )+x_i1)   &
                 +3.0_r8*(wrf%dom(id)%phb(i  ,j,k+1)+x_i2)   &
                        -(wrf%dom(id)%phb(i+1,j,k  )+x_i3) &
                        -(wrf%dom(id)%phb(i+1,j,k+1)+x_i4) )/(4.0_r8*gravity))

         lat = ( 3.0_r8*wrf%dom(id)%latitude(i  ,j)  &
                +3.0_r8*wrf%dom(id)%latitude(i  ,j)  &
                       -wrf%dom(id)%latitude(i+1,j)  &
                       -wrf%dom(id)%latitude(i+1,j)) / 4.0_r8

         model_height_distrib = compute_geometric_height(geop, lat)

      endif

   else

      i1 = get_dart_vector_index(i,j,k  ,domain_id(id),wrf%dom(id)%type_gz)
      i2 = get_dart_vector_index(i,j,k+1,domain_id(id),wrf%dom(id)%type_gz)

      x_i1 = get_state(i1, state_handle)
      x_i2 = get_state(i2, state_handle)
      x_i3 = get_state(i1 -1, state_handle)
      x_i4 = get_state(i2 -1, state_handle)


      geop = minval(( (wrf%dom(id)%phb(i  ,j,k  )+x_i1)   &
              +(wrf%dom(id)%phb(i  ,j,k+1)+x_i2)   &
              +(wrf%dom(id)%phb(i-1,j,k  )+x_i3) &
              +(wrf%dom(id)%phb(i-1,j,k+1)+x_i4) )/(4.0_r8*gravity))

      lat = (  wrf%dom(id)%latitude(i  ,j)  &
              +wrf%dom(id)%latitude(i  ,j)  &
              +wrf%dom(id)%latitude(i-1,j)  &
              +wrf%dom(id)%latitude(i-1,j)) / 4.0_r8

      model_height_distrib = compute_geometric_height(geop, lat)

   endif

! If V-grid, then pressure is defined between V points, both in horizontal 
!   and in vertical, so average -- averaging depends on polar periodicity
elseif( var_type == wrf%dom(id)%type_v ) then

   if( j == wrf%dom(id)%var_size(2,wrf%dom(id)%type_v) ) then

      ! Check to see if periodic in latitude (polar)
      if ( wrf%dom(id)%polar ) then

         ! The upper corner is 180 degrees of longitude away
         off = i + wrf%dom(id)%we/2
         if ( off > wrf%dom(id)%we ) off = off - wrf%dom(id)%we

         i1 = get_dart_vector_index(off,j-1,k  , domain_id(id),wrf%dom(id)%type_gz)
         i2 = get_dart_vector_index(off,j-1,k+1, domain_id(id),wrf%dom(id)%type_gz)
         i3 = get_dart_vector_index(i  ,j-1,k  , domain_id(id),wrf%dom(id)%type_gz)
         i4 = get_dart_vector_index(i  ,j-1,k+1, domain_id(id),wrf%dom(id)%type_gz)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)

         geop = minval(( (wrf%dom(id)%phb(off,j-1,k  )+x_i1) &
                 +(wrf%dom(id)%phb(off,j-1,k+1)+x_i2) &
                 +(wrf%dom(id)%phb(i  ,j-1,k  )+x_i3) &
                 +(wrf%dom(id)%phb(i  ,j-1,k+1)+x_i4) )/(4.0_r8*gravity))
         
         lat = ( wrf%dom(id)%latitude(off,j-1)  &
                +wrf%dom(id)%latitude(off,j-1)  &
                +wrf%dom(id)%latitude(i  ,j-1)  &
                +wrf%dom(id)%latitude(i  ,j-1)) / 4.0_r8

         model_height_distrib = compute_geometric_height(geop, lat)

      else

         ! If not periodic, then try extrapolating
         i1 = get_dart_vector_index(i,j-1,k ,  domain_id(id),wrf%dom(id)%type_gz)
         i2 = get_dart_vector_index(i,j-1,k+1, domain_id(id),wrf%dom(id)%type_gz)
         i3 = get_dart_vector_index(i,j-2,k  , domain_id(id),wrf%dom(id)%type_gz)
         i4 = get_dart_vector_index(i,j-2,k+1, domain_id(id),wrf%dom(id)%type_gz)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)

         geop = minval(( 3.0_r8*(wrf%dom(id)%phb(i,j-1,k  )+x_i1) &
                 +3.0_r8*(wrf%dom(id)%phb(i,j-1,k+1)+x_i2) &
                        -(wrf%dom(id)%phb(i,j-2,k  )+x_i3) &
                        -(wrf%dom(id)%phb(i,j-2,k+1)+x_i4) )/(4.0_r8*gravity))

         lat = ( 3.0_r8*wrf%dom(id)%latitude(i,j-1)  &
                +3.0_r8*wrf%dom(id)%latitude(i,j-1)  &
                       -wrf%dom(id)%latitude(i,j-2)  &
                       -wrf%dom(id)%latitude(i,j-2)) / 4.0_r8

         model_height_distrib = compute_geometric_height(geop, lat)

      endif

   elseif( j == 1 ) then

      ! Check to see if periodic in latitude (polar)
      if ( wrf%dom(id)%polar ) then

         ! The lower corner is 180 degrees of longitude away
         off = i + wrf%dom(id)%we/2
         if ( off > wrf%dom(id)%we ) off = off - wrf%dom(id)%we

         i1 = get_dart_vector_index(off,j,k  , domain_id(id),wrf%dom(id)%type_gz)
         i2 = get_dart_vector_index(off,j,k+1, domain_id(id),wrf%dom(id)%type_gz)
         i3 = get_dart_vector_index(i  ,j,k  , domain_id(id),wrf%dom(id)%type_gz)
         i4 = get_dart_vector_index(i  ,j,k+1, domain_id(id),wrf%dom(id)%type_gz)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)

         geop = minval(( (wrf%dom(id)%phb(off,j,k  )+x_i1) &
                 +(wrf%dom(id)%phb(off,j,k+1)+x_i2) &
                 +(wrf%dom(id)%phb(i  ,j,k  )+x_i3) &
                 +(wrf%dom(id)%phb(i  ,j,k+1)+x_i4) )/(4.0_r8*gravity))
         
         lat = ( wrf%dom(id)%latitude(off,j)  &
                +wrf%dom(id)%latitude(off,j)  &
                +wrf%dom(id)%latitude(i  ,j)  &
                +wrf%dom(id)%latitude(i  ,j)) / 4.0_r8

         model_height_distrib = compute_geometric_height(geop, lat)

      else

         ! If not periodic, then try extrapolating
         i1 = get_dart_vector_index(i,j  ,k  , domain_id(id),wrf%dom(id)%type_gz)
         i2 = get_dart_vector_index(i,j  ,k+1, domain_id(id),wrf%dom(id)%type_gz)
         i3 = get_dart_vector_index(i,j+1,k  , domain_id(id),wrf%dom(id)%type_gz)
         i4 = get_dart_vector_index(i,j+1,k+1, domain_id(id),wrf%dom(id)%type_gz)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)

         geop = minval(( 3.0_r8*(wrf%dom(id)%phb(i,j  ,k  )+x_i1) &
                 +3.0_r8*(wrf%dom(id)%phb(i,j  ,k+1)+x_i2) &
                        -(wrf%dom(id)%phb(i,j+1,k  )+x_i3) &
                        -(wrf%dom(id)%phb(i,j+1,k+1)+x_i4) )/(4.0_r8*gravity))

         lat = ( 3.0_r8*wrf%dom(id)%latitude(i,j  )  &
                +3.0_r8*wrf%dom(id)%latitude(i,j  )  &
                       -wrf%dom(id)%latitude(i,j+1)  &
                       -wrf%dom(id)%latitude(i,j+1)) / 4.0_r8

         model_height_distrib = compute_geometric_height(geop, lat)

      endif

   else

      i1 = get_dart_vector_index(i,j  ,k  , domain_id(id),wrf%dom(id)%type_gz)
      i2 = get_dart_vector_index(i,j  ,k+1, domain_id(id),wrf%dom(id)%type_gz)
      i3 = get_dart_vector_index(i,j-1,k  , domain_id(id),wrf%dom(id)%type_gz)
      i4 = get_dart_vector_index(i,j-1,k+1, domain_id(id),wrf%dom(id)%type_gz)

      x_i1 = get_state(i1, state_handle)
      x_i2 = get_state(i2, state_handle)
      x_i3 = get_state(i3, state_handle)
      x_i4 = get_state(i4, state_handle)

      geop = minval(( (wrf%dom(id)%phb(i,j  ,k  )+x_i1) &
              +(wrf%dom(id)%phb(i,j  ,k+1)+x_i2) &
              +(wrf%dom(id)%phb(i,j-1,k  )+x_i3) &
              +(wrf%dom(id)%phb(i,j-1,k+1)+x_i4) )/(4.0_r8*gravity))

      lat = ( wrf%dom(id)%latitude(i,j  )  &
             +wrf%dom(id)%latitude(i,j  )  &
             +wrf%dom(id)%latitude(i,j-1)  &
             +wrf%dom(id)%latitude(i,j-1)) / 4.0_r8

      model_height_distrib = compute_geometric_height(geop, lat)

   endif

elseif( var_type == wrf%dom(id)%type_mu .or. &
        var_type == wrf%dom(id)%type_ps .or. &
        var_type == wrf%dom(id)%type_tsk) then

   model_height_distrib = wrf%dom(id)%hgt(i,j)

elseif( var_type == wrf%dom(id)%type_tslb  .or. &
        var_type == wrf%dom(id)%type_smois .or. &
        var_type == wrf%dom(id)%type_sh2o ) then

   model_height_distrib = wrf%dom(id)%hgt(i,j) - wrf%dom(id)%zs(k)

elseif( var_type == wrf%dom(id)%type_u10 .or. &
        var_type == wrf%dom(id)%type_v10 ) then

   model_height_distrib = wrf%dom(id)%hgt(i,j) + 10.0_r8

elseif( var_type == wrf%dom(id)%type_t2  .or. &
        var_type == wrf%dom(id)%type_th2 .or. &
        var_type == wrf%dom(id)%type_q2 ) then

   model_height_distrib = wrf%dom(id)%hgt(i,j) + 2.0_r8

else

   i1 = get_dart_vector_index(i,j,k  , domain_id(id),wrf%dom(id)%type_gz)
   i2 = get_dart_vector_index(i,j,k+1, domain_id(id),wrf%dom(id)%type_gz)

   x_i1 = get_state(i1, state_handle)
   x_i2 = get_state(i2, state_handle)

   geop = minval(( (wrf%dom(id)%phb(i,j,k  )+x_i1) &
           +(wrf%dom(id)%phb(i,j,k+1)+x_i2) )/(2.0_r8*gravity))

   lat = wrf%dom(id)%latitude(i,j)

   model_height_distrib = compute_geometric_height(geop, lat)

endif

end function model_height_distrib

!#######################################################

!> Distributed version of model_height_w
!> Only one value, the mean, is used because model_height_w_distrib
!> is only used in the vertical conversion
function model_height_w_distrib(i, j, k, id, state_handle)

! return total height at staggered vertical coordinate
! and horizontal mass coordinates

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: i,j,k,id
real(r8)                        :: x_i1
real(r8)                        :: model_height_w_distrib

integer(i8) :: i1
real(r8)    :: geop

if (wrf%dom(id)%type_gz < 0) then
  call error_handler(E_ERR, 'model_height_w:', &
      'PH must be in state vector to compute staggered model height', &
       source, revision, revdate)
endif

i1 = get_dart_vector_index(i,j,k, domain_id(id),wrf%dom(id)%type_gz)

x_i1 = minval(get_state(i1, state_handle))

geop = (wrf%dom(id)%phb(i,j,k) + x_i1)/gravity
model_height_w_distrib = compute_geometric_height(geop, wrf%dom(id)%latitude(i, j))

end function model_height_w_distrib


!#######################################################

subroutine pert_model_copies(ens_handle, ens_size,  dummy_pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: dummy_pert_amp ! not used
logical,             intent(out)   :: interf_provided

! Perturbs model states for generating initial ensembles.
! Because this requires some care when using - see the comments in the
! code below - you must set a namelist variable to enable this functionality.

! Using this routine is not a substitute for a good set of real initial 
! condition files.  Intended as a last resort, this routine should be used 
! to start a period of free-running to allow the model to spin-up a set of 
! internally consistent states with their own structure before assimilating 
! a set of observations.  A good ensemble of boundary conditions is important
! to evolve the ensemble members differently, which is the goal.
!

real(r8)              :: pert_amount = 0.005   ! 0.5%

real(r8)              :: pert_ampl, range
real(r8)              :: minv, maxv, temp
type(random_seq_type) :: random_seq
integer               :: id, i, j, s, e
logical, allocatable  :: within_range(:)
integer               :: num_variables
real(r8), allocatable :: min_var(:), max_var(:)
integer(i8)           :: start_ind, end_ind
integer(i8), allocatable :: var_list(:)
integer               :: count, copy

! generally you do not want to just perturb a single state to begin an
! experiment, especially for a regional weather model, because the 
! resulting fields will have spread but they won't have organized features.
! we have had good luck with some global atmosphere models where there is
! a lot of model divergence; after a few days of running they evolve into
! plausible conditions that allow assimilation of real obs.
!
! if you really need to start with a single state and proceed, the suggestion
! is to start with small magnitude perturbations and then get a good ensemble
! of boundary conditions and run the model for a while (many days) to let it
! evolve into plausible weather patterns.  then start assimilating real obs.
!
! using this routine requires you to set the new namelist item 
! 'allow_perturbed_ics' to true so you have to read the warnings here or
! in the html docs.
!
! this code will add random noise field by field (T, U, V, etc), and new values
! will not exceed the original max or min values for each field.  this means
! it will not generate illegal values (e.g. negatives for percentages or
! number concentrations) but it also means that if all values in a field are
! identical (e.g. all 0.0) this routine will not change those values.  the code
! can easily be modified to set allowed min and max values here instead of
! using the incoming field min and max values; but you will have to modify
! the code below to enable that functionality.

if (.not. allow_perturbed_ics) then
call error_handler(E_ERR,'pert_model_copies', &
                     'starting WRF model from a single vector requires additional steps', &
                  source, revision, revdate, &
                  text2='see comments in wrf/model_mod.f90::pert_model_copies()')
endif

! NOT REACHED unless allow_perturbed_ics is true in the namelist

! start of pert code
interf_provided = .true.

! Make space for the state vector index numbers that are
! physically located on my task and get the global numbers.

allocate(var_list(get_my_num_vars(ens_handle)))
call get_my_vars(ens_handle, var_list)

! count up the total number of variables across all domains.
num_variables = 0
do id = 1, num_domains
  num_variables = num_variables + wrf%dom(id)%number_of_wrf_variables
enddo

! get the global min/max on a variable by variable basis
allocate(min_var(num_variables), max_var(num_variables))
allocate(within_range(ens_handle%my_num_vars))

count = 1
do id = 1, num_domains
   do i = 1, wrf%dom(id)%number_of_wrf_variables
  
      start_ind = get_index_start(domain_id(id), i)
      end_ind   = get_index_end(domain_id(id), i)

      ! at this point we only have 1 ensemble
      within_range = (var_list >= start_ind .and. var_list <= end_ind)
      min_var(count) = minval(ens_handle%copies(1,:), MASK=within_range)
      max_var(count) = maxval(ens_handle%copies(1,:), MASK=within_range)

      count = count + 1
   enddo
enddo

! find the global min/max values across all tasks.
call all_reduce_min_max(min_var, max_var, num_variables)

deallocate(within_range)

! Now do the perturbing

! using task id as the seed for the random number generator is ok
! because pert_model_copies() is only called once on any single task.  
! it perturbs all ensemble members for the items in the state vector 
! that it owns. because the decomposition will be different with a
! different task count, you will NOT get the same result if you change
! the number of tasks.

call init_random_seq(random_seq, my_task_id()+1)

count = 1 ! min and max are numbered 1 to n, where n is the total number of variables (all domains)
do id = 1, num_domains
   do i = 1, wrf%dom(id)%number_of_wrf_variables

      start_ind = get_index_start(domain_id(id), i)
      end_ind = get_index_end(domain_id(id), i)

      !! Option 1:
      !! make the perturbation amplitude N% of the total
      !! range of this variable.  values could vary a lot
      !! over some of the types, like pressure
      range = max_var(count) - min_var(count)
      pert_ampl = pert_amount * range

      do j=1, ens_handle%my_num_vars
         ! is this state variable index the current variable type we're perturbing?
         if (var_list(j) >= start_ind .and. var_list(j) <= end_ind) then
            do copy = 1, ens_size
               !! Option 2: perturb each value individually
               !! make the perturbation amplitude N% of this value
               !pert_ampl = pert_amount * ens_handle%copies(copy, j)
               ens_handle%copies(copy, j) = random_gaussian(random_seq, ens_handle%copies(copy, j), pert_ampl)
            enddo

            ! keep variable from exceeding the original range
            ens_handle%copies(1:ens_size,j) = max(min_var(count), ens_handle%copies(1:ens_size,j))
            ens_handle%copies(1:ens_size,j) = min(max_var(count), ens_handle%copies(1:ens_size,j))

         endif
      enddo

      count = count + 1

   enddo
enddo

deallocate(var_list, min_var, max_var)

end subroutine pert_model_copies

!-------------------------------------------------------------------------
!>@todo FIXME:
! these routines are only called from the vortex interpolation code, 
! which should be moved to a separate forward operator module.  the forward 
! operator code should define a grid, call model_interpolate() on each
! point of that grid, and then do the same computation it is doing now.


subroutine compute_seaprs ( nz, z, t, p , q ,          &
                            sea_level_pressure, debug)
!-------------------------------------------------------------------------
! compute_seaprs    Estimate sea level pressure.
!
! This routines has been taken "as is" from wrf_user_fortran_util_0.f
!
! This routine assumes
!    index order is (i,j,k)
!    wrf staggering
!    units: pressure (Pa), temperature(K), height (m), mixing ratio (kg kg{-1})
!    availability of 3d p, t, and qv; 2d terrain; 1d half-level zeta string
!    output units of SLP are Pa, but you should divide that by 100 for the
!          weather weenies.
!    virtual effects are included
!
! Dave
!
! cys: change to 1d
! TJH: verified intent() qualifiers, declaration syntax, uses error_handler

      IMPLICIT NONE
      INTEGER,  intent(in)    :: nz
      REAL(r8), intent(in)    :: z(nz), p(nz), q(nz)
      REAL(r8), intent(inout) :: t(nz)
      REAL(r8), intent(out)   :: sea_level_pressure
      LOGICAL,  intent(in)    :: debug

      INTEGER  :: level
      REAL(r8) :: t_surf, t_sea_level

!     Some required physical constants:

      REAL(r8) :: R, G, GAMMA
      PARAMETER (R=287.04_r8, G=9.81_r8, GAMMA=0.0065_r8)

!     Specific constants for assumptions made in this routine:

      REAL(r8) :: TC, PCONST
      PARAMETER (TC=273.16_r8 + 17.5_r8, PCONST = 10000.0_r8)

      LOGICAL  :: ridiculous_mm5_test
      PARAMETER  (ridiculous_mm5_test = .TRUE.)
!     PARAMETER  (ridiculous_mm5_test = .false.)

!     Local variables:

      INTEGER :: k
      INTEGER :: klo, khi

      REAL(r8) :: plo, phi, tlo, thi, zlo, zhi
      REAL(r8) :: p_at_pconst, t_at_pconst, z_at_pconst
      REAL(r8) :: z_half_lowest

      REAL(r8), PARAMETER :: cp           = 7.0_r8*R/2.0_r8
      REAL(r8), PARAMETER :: rcp          = R/cp
      REAL(r8), PARAMETER :: p1000mb      = 100000.0_r8

      LOGICAL ::  l1 , l2 , l3, found

!     Find least zeta level that is PCONST Pa above the surface.  We later use this
!     level to extrapolate a surface pressure and temperature, which is supposed
!     to reduce the effect of the diurnal heating cycle in the pressure field.

      t = t*(p/p1000mb)**rcp

      level = -1

      k = 1
      found = .false.
      do while( (.not. found) .and. (k.le.nz))
         IF ( p(k) .LT. p(1)-PCONST ) THEN
            level = k
            found = .true.
         END IF
         k = k+1
      END DO

      IF ( level .EQ. -1 ) THEN
         PRINT '(A,I4,A)','Troubles finding level ',   &
               NINT(PCONST)/100,' above ground.'
         print*, 'p=',p
         print*, 't=',t
         print*, 'z=',z
         print*, 'q=',q
         write(errstring,*)'Error_in_finding_100_hPa_up'
         call error_handler(E_ERR,'compute_seaprs',errstring,' ',' ',' ')
      END IF


!     Get temperature PCONST Pa above surface.  Use this to extrapolate
!     the temperature at the surface and down to sea level.

      klo = MAX ( level - 1 , 1      )
      khi = MIN ( klo + 1        , nz - 1 )

      IF ( klo .EQ. khi ) THEN
         PRINT '(A)','Trapping levels are weird.'
         PRINT '(A,I3,A,I3,A)','klo = ',klo,', khi = ',khi, &
               ': and they should not be equal.'
         write(errstring,*)'Error_trapping_levels'
         call error_handler(E_ERR,'compute_seaprs',errstring,' ',' ',' ')
      END IF

      plo = p(klo)
      phi = p(khi)
      tlo = t(klo)*(1. + 0.608 * q(klo) )
      thi = t(khi)*(1. + 0.608 * q(khi) )
!     zlo = zetahalf(klo)/ztop*(ztop-terrain(i,j))+terrain(i,j)
!     zhi = zetahalf(khi)/ztop*(ztop-terrain(i,j))+terrain(i,j)
      zlo = z(klo)
      zhi = z(khi)

      p_at_pconst = p(1) - pconst
      t_at_pconst = thi-(thi-tlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)
      z_at_pconst = zhi-(zhi-zlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)

      t_surf = t_at_pconst*(p(1)/p_at_pconst)**(gamma*R/g)
      t_sea_level = t_at_pconst+gamma*z_at_pconst


!     If we follow a traditional computation, there is a correction to the sea level
!     temperature if both the surface and sea level temnperatures are *too* hot.

      IF ( ridiculous_mm5_test ) THEN
         l1 = t_sea_level .LT. TC
         l2 = t_surf      .LE. TC
         l3 = .NOT. l1
         IF ( l2 .AND. l3 ) THEN
            t_sea_level = TC
         ELSE
            t_sea_level = TC - 0.005*(t_surf-TC)**2
         END IF
      END IF

!     The grand finale: ta da!

!     z_half_lowest=zetahalf(1)/ztop*(ztop-terrain(i,j))+terrain(i,j)
      z_half_lowest=z(1)
      sea_level_pressure = p(1) *              &
                           EXP((2.*g*z_half_lowest)/   &
                           (R*(t_sea_level+t_surf)))

!        sea_level_pressure(i,j) = sea_level_pressure(i,j)*0.01

    if (debug) then
      print *,'slp=',sea_level_pressure
    endif
!      print *,'t=',t(10:15,10:15,1),t(10:15,2,1),t(10:15,3,1)
!      print *,'z=',z(10:15,1,1),z(10:15,2,1),z(10:15,3,1)
!      print *,'p=',p(10:15,1,1),p(10:15,2,1),p(10:15,3,1)
!      print *,'slp=',sea_level_pressure(10:15,10:15),     &
!         sea_level_pressure(10:15,10:15),sea_level_pressure(20,10:15)

end subroutine compute_seaprs


      
SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      INTEGER m,n,NN
      REAL(r8), intent(in) :: x1,x2,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      real(r8), intent(out) :: y
      PARAMETER (NN=200)
!     USES spline,splint
      INTEGER j,k
      REAL(r8) y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12    continue
      call spline(x1a,yytmp,m,1.e30_r8,1.e30_r8,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
END subroutine splin2

SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      INTEGER m,n,NN
      REAL(r8), intent(in) :: x1a(m),x2a(n),ya(m,n)
      REAL(r8), intent(out) :: y2a(m,n)
      PARAMETER (NN=200)
!     USES spline
      INTEGER j,k
      REAL(r8) y2tmp(NN),ytmp(NN)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline(x2a,ytmp,n,1.e30_r8,1.e30_r8,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      return
END subroutine splie2

SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL(r8), intent(in) :: yp1,ypn,x(n),y(n)
      REAL(r8), intent(out) :: y2(n)
      PARAMETER (NMAX=800)
      INTEGER i,k
      REAL(r8) p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+ &
      1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
      u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
END subroutine spline


SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL(r8),intent(in) :: x,xa(n),y2a(n),ya(n)
      REAL(r8),intent(out) :: y
      INTEGER k,khi,klo
      REAL(r8) a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
         write(*,*) 'bad xa input in splint'
         stop
      endif 
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
END subroutine splint


!#######################################################################

subroutine get_domain_info(obslon,obslat,id,iloc,jloc,domain_id_start)

real(r8), intent(in)           :: obslon, obslat
integer,  intent(out)          :: id
real(r8), intent(out)          :: iloc, jloc
integer,  intent(in), optional :: domain_id_start

logical               :: dom_found

! given arbitrary lat and lon values, returns closest domain id and
! horizontal mass point grid points (xloc,yloc)

dom_found = .false.

! the default is to start at the innermost domain and stop when
! the location is found.  however if you want to start at a particular
! domain id number, pass it in as the last optional arg.
id = num_domains
if (present(domain_id_start)) then
   if (domain_id_start < 1 .or. domain_id_start > num_domains) then
      write(errstring,  '(A,I1)') 'bad domain_id_start: ', domain_id_start
      write(msgstring2, '(A,I1)') 'must be between 1 and ', num_domains
      call error_handler(E_ERR, 'model_mod', errstring, &
                         source, revision, revdate, text2=msgstring2)
   endif
   id = domain_id_start
endif

do while (.not. dom_found)

   ! Checking for exact equality on real variable types is generally a bad idea.

   if( (wrf%dom(id)%proj%hemi ==  1.0_r8 .and. obslat < -90.0_r8) .or. &
       (wrf%dom(id)%proj%hemi == -1.0_r8 .and. obslat >  90.0_r8) .or. &
       (wrf%dom(id)%proj%code == PROJ_MERC .and. abs(obslat) > 90.0_r8) ) then

      ! catch latitudes that are out of range - ignore them but print out a warning.
      write(errstring, *) 'obs with latitude out of range: ', obslat
      call error_handler(E_MSG, 'model_mod', errstring)

   else
      call latlon_to_ij(wrf%dom(id)%proj,min(max(obslat,-89.9999999_r8),89.9999999_r8),obslon,iloc,jloc)

      ! Array bound checking depends on whether periodic or not -- these are
      !   real-valued indices here, so we cannot use boundsCheck  :( 

      if ( wrf%dom(id)%periodic_x .and. .not. wrf%dom(id)%periodic_y  ) then
         if ( wrf%dom(id)%polar ) then        
            !   Periodic     X & M_grid ==> [1 we+1)    
            !   Periodic     Y & M_grid ==> [0.5 sn+0.5]
            if ( iloc >= 1.0_r8 .and. iloc <  real(wrf%dom(id)%we,r8)+1.0_r8 .and. &
                 jloc >= 0.5_r8 .and. jloc <= real(wrf%dom(id)%sn,r8)+0.5_r8 ) &
                 dom_found = .true.     
         else
            !   Periodic     X & M_grid ==> [1 we+1)    
            !   NOT Periodic Y & M_grid ==> [1 sn]
            if ( iloc >= 1.0_r8 .and. iloc <  real(wrf%dom(id)%we,r8)+1.0_r8 .and. &
                 jloc >= 1.0_r8 .and. jloc <= real(wrf%dom(id)%sn,r8) ) &
                 dom_found = .true.
         endif

      elseif ( wrf%dom(id)%periodic_x .and. wrf%dom(id)%periodic_y ) then
            !   Periodic     X & M_grid ==> [1 we+1)    
            !   Periodic     Y & M_grid ==> [1 sn+1]
            if ( iloc >= 1.0_r8 .and. iloc <  real(wrf%dom(id)%we,r8)+1.0_r8 .and. &
                 jloc >= 1.0_r8 .and. jloc <= real(wrf%dom(id)%sn,r8)+1.0_r8 ) &
                 dom_found = .true.

      else
         if ( wrf%dom(id)%polar ) then        
            !   NOT Periodic X & M_grid ==> [1 we]    
            !   Periodic     Y & M_grid ==> [0.5 sn+0.5]
            if ( iloc >= 1.0_r8 .and. iloc <= real(wrf%dom(id)%we,r8) .and. &
                 jloc >= 0.5_r8 .and. jloc <= real(wrf%dom(id)%sn,r8)+0.5_r8 ) &
                 dom_found = .true.     
         else
            !   NOT Periodic X & M_grid ==> [1 we]    
            !   NOT Periodic Y & M_grid ==> [1 sn]
            if ( iloc >= 1.0_r8 .and. iloc <= real(wrf%dom(id)%we,r8) .and. &
                 jloc >= 1.0_r8 .and. jloc <= real(wrf%dom(id)%sn,r8) ) &
                 dom_found = .true.
         endif 
      endif

   endif

   if (.not. dom_found) then
      id = id - 1
      if (id == 0) return
   endif

end do

end subroutine get_domain_info

!#######################################################################
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, state_handle)

type(get_close_type),          intent(in)     :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)     :: base_type, loc_qtys(:)
integer(i8),                   intent(in)     :: loc_indx(:)
integer,                       intent(out)    :: num_close, close_ind(:)
real(r8),            optional, intent(out)    :: dist(:)
type(ensemble_type), optional, intent(in)     :: state_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist, state_handle)

end subroutine get_close_state

!#######################################################################
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, state_handle)

type(get_close_type),          intent(in)     :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)     :: base_type, loc_qtys(:), loc_types(:)
integer,                       intent(out)    :: num_close, close_ind(:)
real(r8),            optional, intent(out)    :: dist(:)
type(ensemble_type), optional, intent(in)     :: state_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist, state_handle)

end subroutine get_close_obs

!#######################################################################
subroutine get_close(gc, base_loc, base_type, locs, loc_qtys, &
                     num_close, close_ind, dist, state_handle)

! Given a DART ob (referred to as "base") and a set of obs priors or state variables
! (obs_loc, obs_kind), returns the subset of close ones to the "base" ob, their
! indices, and their distances to the "base" ob...

! For vertical distance computations, general philosophy is to convert all vertical
! coordinates to a common coordinate. This coordinate type is defined in the namelist
! with the variable "vert_localization_coord".

! Vertical conversion is carried out by the subroutine vert_convert.

! Note that both base_loc and locs are intent(inout), meaning that these
! locations are possibly modified here and returned as such to the calling routine.
! The calling routine is always filter_assim and these arrays are local arrays
! within filter_assim. In other words, these modifications will only matter within
! filter_assim, but will not propagate backwards to filter.

type(get_close_type),          intent(in)     :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)     :: base_type, loc_qtys(:)
integer,                       intent(out)    :: num_close, close_ind(:)
real(r8),            optional, intent(out)    :: dist(:)
type(ensemble_type), optional, intent(in)     :: state_handle

integer                :: t_ind, istatus1, istatus2, k
integer                :: base_which, local_which
real(r8), dimension(3) :: base_array, local_array
type(location_type)    :: local_loc


! Initialize variables to missing status
num_close = 0
close_ind = -99
if (present(dist)) dist      = 1.0e9

istatus1 = 0
istatus2 = 0

! Convert base vertical coordinate to requested vertical coordinate if necessary

base_array = get_location(base_loc) 
base_which = nint(query_location(base_loc))

if (vertical_localization_on()) then
   if (base_which /= wrf%dom(1)%localization_coord) then
      !print*, 'base_which ', base_which, 'loc coord ', wrf%dom(1)%localization_coord
      call vert_convert(state_handle, base_loc, base_type, istatus1)
      !call error_handler(E_ERR, 'you should not call this ', 'get_close_obs')
   elseif (base_array(3) == missing_r8) then
      istatus1 = 1
   endif
endif

if (istatus1 == 0) then

   ! Get all the potentially close obs but no dist (optional argument dist(:) is not present)
   ! This way, we are decreasing the number of distance computations that will follow.
   ! This is a horizontal-distance operation and we don't need to have the relevant vertical
   ! coordinate information yet (for obs_loc).
   call loc_get_close(gc, base_loc, base_type, locs, loc_qtys, loc_qtys, &
                          num_close, close_ind)

   ! Loop over potentially close subset of obs priors or state variables
   do k = 1, num_close

      t_ind = close_ind(k)
      local_loc   = locs(t_ind)
      local_which = nint(query_location(local_loc))

      if (present(dist)) then
         if (local_which == VERTISUNDEF) then
            dist(k) = get_dist(base_loc, local_loc, base_type, loc_qtys(t_ind))
            cycle
         endif
      endif

      ! Convert local vertical coordinate to requested vertical coordinate if necessary.
      ! This should only be necessary for obs priors, as state location information already
      ! contains the correct vertical coordinate (filter_assim's call to get_state_meta_data).
      if (vertical_localization_on()) then
         if (local_which /= wrf%dom(1)%localization_coord) then
            call vert_convert(state_handle, local_loc, loc_qtys(t_ind), istatus2)
            ! Store the "new" location into the original full local array
            locs(t_ind) = local_loc !HK Overwritting the location
         else
            istatus2 = 0
         endif
      endif

      if (present(dist)) then
         ! Compute distance - set distance to a very large value if vert coordinate is missing
         ! or vert_convert returned error (istatus2=1)
         local_array = get_location(local_loc)
         if (((vertical_localization_on()).and.(local_array(3) == missing_r8)).or.(istatus2 == 1)) then
            dist(k) = 1.0e9
         else
            dist(k) = get_dist(base_loc, local_loc, base_type, loc_qtys(t_ind))
         endif
      endif

      !print*, 'k ', k, 'rank ', my_task_id()

   end do
endif

end subroutine get_close

!#######################################################################
!nc -- additional function from Greg Lawson & Nancy Collins
!
!  logical function boundsCheck determines whether real-valued location indices are
!    within a sensible range based on the assumed (un)staggered grid and based on 
!    whether the domain is assumed to be periodic in a given direction.

function boundsCheck ( ind, periodic, id, dim, type )

  integer,  intent(in)  :: ind, id, dim, type
  logical,  intent(in)  :: periodic

  logical :: boundsCheck  
!  logical, parameter :: restrict_polar = .true.
  logical, parameter :: restrict_polar = .false.

  ! Consider cases in REAL-VALUED indexing:
  !
  ! I. Longitude -- x-direction
  !    A. PERIODIC (period_x = .true.)
  !
  !       Consider Mass-grid (& V-grid) longitude grid with 4 west-east gridpoints
  !         Values  ::  [ -135 -45  45 135 ] .. {225}
  !         Indices ::  [   1   2   3   4  ] .. {1,5}
  !       Complementary U-grid
  !         Values  ::  [ -180 -90  0  90  180 ]
  !         Indices ::  [   1   2   3   4   5  ]
  !
  !       What are the allowable values for a real-valued index on each of these grids?
  !       1. M-grid  --->  [1 5)       ---> [1 we+1)
  !                  --->  [-135 225)  
  !       2. U-grid  --->  [1 5)       ---> [1 wes)
  !                  --->  [-180 180)
  !       [Note that above "allowable values" reflect that one should be able to have
  !        an observation anywhere on a given longitude circle -- the information 
  !        exists in order to successfully interpolate to anywhere over [0 360).]
  !
  !       It is up to the routine calling "boundsCheck" to have handled the 0.5 offset
  !         in indices between the M-grid & U-grid.  Hence, two examples: 
  !          a. If there is an observation location at -165 longitude, then:
  !             * An observation of TYPE_T (on the M-grid) would have ind = 4.667
  !             * An observation of TYPE_U (on the U-grid) would have ind = 1.167
  !          b. If there is an observation location at 0 longitude, then:
  !             * An observation of TYPE_T (on the M-grid) would have ind = 2.5
  !             * An observation of TYPE_U (on the U-grid) would have ind = 3.0
  !
  !    B. NOT periodic (period_x = .false.)
  !
  !       Consider Mass-grid (& V-grid) longitude grid with 4 west-east gridpoints
  !         Values  ::  [  95  105 115 125 ] 
  !         Indices ::  [   1   2   3   4  ] 
  !       Complementary U-grid
  !         Values  ::  [  90  100 110 120 130 ]
  !         Indices ::  [   1   2   3   4   5  ]
  !
  !       What are the allowable values for a real-valued index on each of these grids?
  !       1. M-grid  --->  [1 4]       ---> [1 we]
  !                  --->  [95 125]  
  !       2. U-grid  --->  [1.5 4.5]       ---> [1.5 we+0.5]
  !                  --->  [95 125]
  !       [Note that above "allowable values" reflect that one should only be able to
  !        have an observation within the M-grid, since that is the only way to  
  !        guarantee that the necessary information exists in order to successfully 
  !        interpolate to a specified location.]
  !
  !       It is up to the routine calling "boundsCheck" to have handled the 0.5 offset
  !         in indices between the M-grid & U-grid.  Hence, two examples: 
  !          a. If there is an observation location at 96 longitude, then:
  !             * An observation of TYPE_T (on the M-grid) would have ind = 1.1
  !             * An observation of TYPE_U (on the U-grid) would have ind = 1.6
  !          b. If there is an observation location at 124 longitude, then:
  !             * An observation of TYPE_T (on the M-grid) would have ind = 3.9
  !             * An observation of TYPE_U (on the U-grid) would have ind = 4.4
  !
  ! II. Latitude -- y-direction
  !    A. PERIODIC (polar = .true.)
  !
  !       Consider Mass-grid (& U-Grid) latitude grid with 4 south-north gridpoints
  !         Values  :: [ -67.5 -22.5  22.5  67.5 ] 
  !         Indices :: [   1     2     3     4   ] 
  !       Complementary V-grid 
  !         Values  :: [ -90   -45     0    45    90 ] 
  !         Indices :: [   1     2     3     4     5 ] 
  !
  !       What are the allowable values for a real-valued index on each of these grids?
  !       1. M-grid  --->  [0.5 4.5]   ---> [0.5 sn+0.5]
  !                  --->  [-90 90]  
  !       2. U-grid  --->  [1 5]       ---> [1 sns]
  !                  --->  [-90 90]
  !       [Note that above "allowable values" reflect that one should be able to have
  !        an observation anywhere along a give latitude circle -- the information 
  !        exists in order to successfully interpolate to anywhere over [-90 90]; 
  !        however, in latitude this poses a special challenge since the seams join
  !        two separate columns of data over the pole, as opposed to in longitude
  !        where the seam wraps back on a single row of data.]  
  !
  !       It is up to the routine calling "boundsCheck" to have handled the 0.5 offset
  !         in indices between the M-grid & V-grid.  Hence, two examples: 
  !          a. If there is an observation location at -75 latitude, then:
  !             * An observation of TYPE_T (on the M-grid) would have ind = 0.833
  !             * An observation of TYPE_V (on the V-grid) would have ind = 1.333
  !          b. If there is an observation location at 0 latitude, then:
  !             * An observation of TYPE_T (on the M-grid) would have ind = 2.5
  !             * An observation of TYPE_V (on the V-grid) would have ind = 3.0
  !
  !    B. NOT periodic (polar = .false.)
  !
  !       Consider Mass-grid (& U-Grid) latitude grid with 4 south-north gridpoints
  !         Values  :: [ 10  20  30  40 ] 
  !         Indices :: [  1   2   3   4 ] 
  !       Complementary V-grid 
  !         Values  :: [  5  15  25  35  45 ] 
  !         Indices :: [  1   2   3   4   5 ] 
  !
  !       What are the allowable values for a real-valued index on each of these grids?
  !       1. M-grid  --->  [1 4]   ---> [1 sn]
  !                  --->  [10 40]  
  !       2. U-grid  --->  [1.5 4.5]       ---> [1.5 sn+0.5]
  !                  --->  [10 40]
  !       [Note that above "allowable values" reflect that one should only be able to
  !        have an observation within the M-grid, since that is the only way to  
  !        guarantee that the necessary information exists in order to successfully 
  !        interpolate to a specified location.]
  !
  !       It is up to the routine calling "boundsCheck" to have handled the 0.5 offset
  !         in indices between the M-grid & V-grid.  Hence, two examples: 
  !          a. If there is an observation location at 11 latitude, then:
  !             * An observation of TYPE_T (on the M-grid) would have ind = 1.1
  !             * An observation of TYPE_V (on the V-grid) would have ind = 1.6
  !          b. If there is an observation location at 25 latitude, then:
  !             * An observation of TYPE_T (on the M-grid) would have ind = 2.5
  !             * An observation of TYPE_V (on the V-grid) would have ind = 3.0
  ! 
  ! III. Vertical -- z-direction (periodicity not an issue)
  !    
  !    Consider Mass vertical grid with 4 bottom-top gridpoints
  !      Values  :: [ 0.875 0.625 0.375 0.125 ]
  !      Indices :: [   1     2     3     4   ]
  !    Complementary W-grid
  !      Values  :: [   1   0.75  0.50  0.25    0   ]
  !      Indices :: [   1     2     3     4     5   ]
  !
  !    What are the allowable values for a real-valued index on each of these grids?
  !    1. M-grid  --->  [1 4]           ---> [1 bt]
  !               --->  [0.875 0.125]  
  !    2. W-grid  --->  [1.5 4.5]       ---> [1.5 bt+0.5]
  !               --->  [0.875 0.125]
  !
  !    [Note that above "allowable values" reflect that one should only be able to
  !     have an observation within the M-grid, since that is the only way to  
  !     guarantee that the necessary information exists in order to successfully 
  !     interpolate to a specified location.]
  !

  ! Summary of Allowable REAL-VALUED Index Values ==> INTEGER Index Values 
  !
  ! In longitude (x) direction
  !   Periodic     & M_grid ==> [1 we+1)       ==> [1 wes)
  !   Periodic     & U_grid ==> [1 wes)        ==> [1 wes)
  !   NOT Periodic & M_grid ==> [1 we]         ==> [1 we)
  !   NOT Periodic & U_grid ==> [1.5 we+0.5]   ==> [1 wes)
  ! In latitude (y) direction
  !   Periodic     & M_grid ==> [0.5 sn+0.5]   ==> [0 sns) *though in practice, [1 sn)*
  !   Periodic     & V_grid ==> [1 sns]        ==> [1 sns) *though allowable range, [1.5 sn+.5]*
  !   NOT Periodic & M_grid ==> [1 sn]         ==> [1 sn)
  !   NOT Periodic & V_grid ==> [1.5 sn+0.5]   ==> [1 sns)
  ! In vertical (z) direction
  !                  M_grid ==> [1 bt]         ==> [1 bt)
  !                  W_grid ==> [1.5 bt+0.5]   ==> [1 bts)
  

  ! Assume boundsCheck is false unless we can prove otherwise
  boundsCheck = .false.

  ! First check direction (dimension)
  !   Longitude (x-direction) has dim == 1
  if ( dim == 1 ) then

     ! Next check periodicity
     if ( periodic ) then
        
        ! If periodic in longitude, then no need to check staggering because both
        !   M and U grids allow integer indices from [1 wes)
        if ( ind >= 1 .and. ind < wrf%dom(id)%wes ) boundsCheck = .true.

     else

        ! If NOT periodic in longitude, then we need to check staggering because
        !   M and U grids allow different index ranges

        ! Check staggering by comparing var_size(dim,type) to the staggered dimension 
        !   for dim == 1 stored in wrf%dom(id)
        if ( wrf%dom(id)%var_size(dim,type) == wrf%dom(id)%wes ) then
           ! U-grid allows integer range of [1 wes)
           if ( ind >= 1 .and. ind < wrf%dom(id)%wes ) boundsCheck = .true.
        else  
           ! M & V-grid allow [1 we)
           if ( ind >= 1 .and. ind < wrf%dom(id)%we ) boundsCheck = .true.
        endif

     endif

   !   Latitude (y-direction) has dim == 2
   elseif ( dim == 2 ) then

     ! Next check periodicity
     if ( periodic ) then
        
        ! We need to check staggering because M and V grids allow different indices

!*** NOTE: For now are disallowing observation locations that occur poleward of the 
!            first and last M-grid gridpoints.  This means that this function will 
!            return false for polar observations.  This need not be the case because
!            the information should be available for proper interpolation across the
!            poles, but it will require more clever thinking.  Hopefully this can 
!            be added in later.  

        ! Check staggering by comparing var_size(dim,type) to the staggered dimension 
        !   for dim == 2 stored in wrf%dom(id)
        if ( wrf%dom(id)%var_size(dim,type) == wrf%dom(id)%sns ) then
           ! V-grid allows integer range [1 sns)
           if ( ind >= 1 .and. ind < wrf%dom(id)%sns ) boundsCheck = .true.
        else  
           ! For now we will set a logical flag to more restrictively check the array
           !   bounds under our no-polar-obs assumptions
           if ( restrict_polar ) then
              ! M & U-grid allow integer range [1 sn) in practice (though properly, [0 sns) )
              if ( ind >= 1 .and. ind < wrf%dom(id)%sn ) boundsCheck = .true.
           else
              ! M & U-grid allow integer range [0 sns) in unrestricted circumstances
              if ( ind >= 0 .and. ind < wrf%dom(id)%sns ) boundsCheck = .true.
           endif
        endif
        
     else

        ! We need to check staggering because M and V grids allow different indices
        if ( wrf%dom(id)%var_size(dim,type) == wrf%dom(id)%sns ) then
           ! V-grid allows [1 sns)
           if ( ind >= 1 .and. ind < wrf%dom(id)%sns ) boundsCheck = .true.
        else 
           ! M & U-grid allow [1 sn)
           if ( ind >= 1 .and. ind < wrf%dom(id)%sn ) boundsCheck = .true.
        endif

     endif

  elseif ( dim == 3 ) then

     ! No periodicity to worry about in the vertical!  However, we still need to check
     !   staggering because the ZNU and ZNW grids allow different index ranges
     if ( wrf%dom(id)%var_size(dim,type) == wrf%dom(id)%bts ) then
        ! W vertical grid allows [1 bts)
        if ( ind >= 1 .and. ind < wrf%dom(id)%bts ) boundsCheck = .true.
     else
        ! M vertical grid allows [1 bt)
        if ( ind >= 1 .and. ind < wrf%dom(id)%bt ) boundsCheck = .true.
     endif
  
  else

     print*, 'model_mod.f90 :: function boundsCheck :: dim must equal 1, 2, or 3!'

  endif


end function boundsCheck

!#######################################################################
! get_orners takes in an i and j index, information about domain and grid staggering,
!   and then returns the four cornering gridpoints' 2-element integer indices. 
subroutine getCorners(i, j, id, type, ll, ul, lr, ur, rc)

  integer, intent(in)  :: i, j, id, type
  integer, dimension(2), intent(out) :: ll, ul, lr, ur
  integer, intent(out) :: rc

!  logical, parameter :: restrict_polar = .true.
  logical, parameter :: restrict_polar = .false.

  ! set return code to 0, and change this if necessary
  rc = 0

  !----------------
  ! LOWER LEFT
  !----------------

  ! i and j are the lower left (ll) corner already
  !
  ! NOTE :: once we allow for polar periodicity, the incoming j index could actually 
  !           be 0, which would imply a ll(2) value of 1, with a ll(1) value 180 degrees
  !           of longitude away from the incoming i index!  But we have not included 
  !           this possibility yet.  

  ! As of 22 Oct 2007, this option is not allowed!
  !   Note that j = 0 can only happen if we are on the M (or U) wrt to latitude
  if ( wrf%dom(id)%polar .and. j == 0 .and. .not. restrict_polar ) then

     ! j = 0 should be mapped to j = 1 (ll is on other side of globe)
     ll(2) = 1
     
     ! Need to map i index 180 degrees away
     ll(1) = i + wrf%dom(id)%we/2
     
     ! Check validity of bounds & adjust by periodicity if necessary
     if ( ll(1) > wrf%dom(id)%we ) ll(1) = ll(1) - wrf%dom(id)%we

     ! We shouldn't be able to get this return code if restrict_polar = .true.
!     rc = 1
!     print*, 'model_mod.f90 :: getCorners :: Tried to do polar bc -- rc = ', rc

  else
     
     ll(1) = i
     ll(2) = j

  endif


  !----------------
  ! LOWER RIGHT
  !----------------

  ! Most of the time, the lower right (lr) corner will simply be (i+1,j), but we need to check
  ! Summary of x-direction corners:
  !   Periodic     & M_grid has ind = [1 wes)
  !     ind = [1 we)    ==> ind_p_1 = ind + 1
  !     ind = [we wes)  ==> ind_p_1 = 1
  !   Periodic     & U_grid has ind = [1 wes)
  !     ind = [1 we)    ==> ind_p_1 = ind + 1
  !     ind = [we wes)  ==> ind_p_1 = wes       ( keep in mind that U(1) = U(wes) if periodic )
  !   NOT Periodic & M_grid has ind = [1 we)
  !     ind = [1 we-1)  ==> ind_p_1 = ind + 1
  !     ind = [we-1 we) ==> ind_p_1 = we
  !   NOT Periodic & U_grid has ind = [1 wes)
  !     ind = [1 we)    ==> ind_p_1 = ind + 1
  !     ind = [we wes)  ==> ind_p_1 = wes 

  if ( wrf%dom(id)%periodic_x ) then
    
     ! Check to see what grid we have, M vs. U
     if ( wrf%dom(id)%var_size(1,type) == wrf%dom(id)%wes ) then
        ! U-grid is always i+1 -- do this in reference to already adjusted ll points
        lr(1) = ll(1) + 1
        lr(2) = ll(2)
     else
        ! M-grid is i+1 except if we <= ind < wes, in which case it's 1
        if ( i < wrf%dom(id)%we ) then
           lr(1) = ll(1) + 1
        else
           lr(1) = 1
        endif
        lr(2) = ll(2)
     endif

  else

     ! Regardless of grid, NOT Periodic always has i+1
     lr(1) = ll(1) + 1
     lr(2) = ll(2)

  endif
        

  !----------------
  ! UPPER LEFT
  !----------------

!*** NOTE: For now are disallowing observation locations that occur poleward of the 
!            first and last M-grid gridpoints.  This need not be the case because
!            the information should be available for proper interpolation across the
!            poles, but it will require more clever thinking.  Hopefully this can 
!            be added in later.

  ! Most of the time, the upper left (ul) corner will simply be (i,j+1), but we need to check
  ! Summary of y-direction corners:
  !   Periodic     & M_grid has ind = [0 sns)  *though in practice, [1 sn)*
  !     ind = [1 sn-1)  ==> ind_p_1 = ind + 1
  !     ind = [sn-1 sn) ==> ind_p_1 = sn
  !   Periodic     & V_grid has ind = [1 sns) 
  !     ind = [1 sn)    ==> ind_p_1 = ind + 1
  !     ind = [sn sns)  ==> ind_p_1 = sns  
  !   NOT Periodic & M_grid has ind = [1 sn)
  !     ind = [1 sn-1)  ==> ind_p_1 = ind + 1
  !     ind = [sn-1 sn) ==> ind_p_1 = sn
  !   NOT Periodic & V_grid has ind = [1 sns)
  !     ind = [1 sn)    ==> ind_p_1 = ind + 1
  !     ind = [sn sns)  ==> ind_p_1 = sns 
  !
  ! Hence, with our current polar obs restrictions, all four possible cases DO map into
  !   ul = (i,j+1).  But this will not always be the case.
  
  if ( wrf%dom(id)%polar ) then

     ! Check to see what grid we have, M vs. V
     if ( wrf%dom(id)%var_size(2,type) == wrf%dom(id)%sns ) then
        ! V-grid is always j+1, even if we allow for full [1 sns) range
        ul(1) = ll(1)
        ul(2) = ll(2) + 1
     else
        ! M-grid changes depending on polar restriction
        if ( restrict_polar ) then 
           ! If restricted, then we can simply add 1
           ul(1) = ll(1)
           ul(2) = ll(2) + 1
        else
           ! If not restricted, then we can potentially wrap over the north pole, which 
           !   means that ul(2) is set to sn and ul(1) is shifted by 180 deg.

           if ( j == wrf%dom(id)%sn ) then
              ! j > sn should be mapped to j = sn (ul is on other side of globe)
              ul(2) = wrf%dom(id)%sn
     
              ! Need to map i index 180 degrees away
              ul(1) = ll(1) + wrf%dom(id)%we/2
     
              ! Check validity of bounds & adjust by periodicity if necessary
              if ( ul(1) > wrf%dom(id)%we ) ul(1) = ul(1) - wrf%dom(id)%we

              ! We shouldn't be able to get this return code if restrict_polar = .true.
!              rc = 1
!              print*, 'model_mod.f90 :: getCorners :: Tried to do polar bc -- rc = ', rc

           elseif ( j == 0 ) then
              ! In this case, we have place ll on the other side of the globe, so we 
              !   cannot reference ul to ll
              ul(1) = i
              ul(2) = 1

           else
              ! We can confidently set to j+1
              ul(1) = ll(1)
              ul(2) = ll(2) + 1
           endif

        endif
     endif

  elseif ( wrf%dom(id)%periodic_y ) then

     ! Check to see what grid we have, M vs. U
     if ( wrf%dom(id)%var_size(2,type) == wrf%dom(id)%sns ) then
        ! V-grid is always j+1 -- do this in reference to already adjusted ll points
        ul(1) = ll(1)
        ul(2) = ll(2)+1
     else
        ! M-grid is j+1 except if we <= ind < wes, in which case it's 1
        if ( j < wrf%dom(id)%sn ) then
           ul(2) = ll(2) + 1
        else
           ul(2) = 1
        endif
        ul(1) = ll(1)
     endif

  else

     ! Regardless of grid, NOT Periodic always has j+1
     ul(1) = ll(1) 
     ul(2) = ll(2) + 1

  endif
     

  !----------------
  ! UPPER RIGHT
  !----------------

!*** NOTE: For now are disallowing observation locations that occur poleward of the 
!            first and last M-grid gridpoints.  This need not be the case because
!            the information should be available for proper interpolation across the
!            poles, but it will require more clever thinking.  Hopefully this can 
!            be added in later.

  ! Most of the time, the upper right (ur) corner will simply be (i+1,j+1), but we need to check
  !   In fact, we can largely get away with ur = (lr(1),ul(2)).  Where this will NOT work is
  !   where we have had to re-map the i index to the other side of the globe (180 deg) due to 
  !   the polar boundary condition.  There are no situations where ur(2) will not be equal to
  !   ul(2).  

  ur(2) = ul(2)

  ! Need to check if ur(1) .ne. lr(1)
  if ( wrf%dom(id)%polar .and. .not. restrict_polar ) then

     ! Only if j == 0 or j == sn
     if ( j == 0 .or. j ==  wrf%dom(id)%sn) then
        ! j == 0 means that ll(1) = i + 180 deg, so we cannot use lr(1) -- hence, we will
        !   add 1 to ul(1), unless doing so spans the longitude seam point.
        ! j == sn means that ul(1) = i + 180 deg.  Here we cannot use lr(1) either because
        !   it will be half a domain away from ul(1)+1.  Be careful of longitude seam point.

        !   Here we need to check longitude periodicity and the type of grid
        if ( wrf%dom(id)%periodic_x ) then
    
           ! Check to see what grid we have, M vs. U
           if ( wrf%dom(id)%var_size(1,type) == wrf%dom(id)%wes ) then
              ! U-grid is always i+1 -- do this in reference to already adjusted ll points
              ur(1) = ul(1) + 1
           else
              ! M-grid is i+1 except if we <= ind < wes, in which case it's 1
              if ( ul(1) < wrf%dom(id)%we ) then
                 ur(1) = ul(1) + 1
              else
                 ur(1) = 1
              endif
           endif

        else

           ! Regardless of grid, NOT Periodic always has i+1
           ur(1) = ul(1) + 1

        endif

     ! If not a special j value, then we are set for the ur(1) = lr(1)
     else

        ur(1) = lr(1)

     endif

  ! If not an unrestricted polar periodic domain, then we have nothing to worry about
  else

     ur(1) = lr(1)

  endif

end subroutine getCorners

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine read_wrf_dimensions(ncid,bt,bts,sn,sns,we,wes,sls)

! ncid: input, file handl
! id:   input, domain id

integer, intent(in)            :: ncid
integer, intent(out)           :: bt,bts,sn,sns,we,wes,sls
logical, parameter             :: debug = .false.
integer                        :: var_id 
character (len=NF90_MAX_NAME)  :: name

! get wrf grid dimensions

   call nc_check( nf90_inq_dimid(ncid, "bottom_top", var_id), &
                     'static_init_model','inq_dimid bottom_top')
   call nc_check( nf90_inquire_dimension(ncid, var_id, name, bt), &
                     'static_init_model','inquire_dimension '//trim(name))

   call nc_check( nf90_inq_dimid(ncid, "bottom_top_stag", var_id), &
                     'static_init_model','inq_dimid bottom_top_stag') ! reuse var_id, no harm
   call nc_check( nf90_inquire_dimension(ncid, var_id, name, bts), &
                     'static_init_model','inquire_dimension '//trim(name))

   call nc_check( nf90_inq_dimid(ncid, "south_north", var_id), &
                     'static_init_model','inq_dimid south_north')
   call nc_check( nf90_inquire_dimension(ncid, var_id, name, sn), &
                     'static_init_model','inquire_dimension '//trim(name))

   call nc_check( nf90_inq_dimid(ncid, "south_north_stag", var_id), &
                     'static_init_model','inq_dimid south_north_stag') ! reuse var_id, no harm
   call nc_check( nf90_inquire_dimension(ncid, var_id, name, sns), &
                     'static_init_model','inquire_dimension '//trim(name))

   call nc_check( nf90_inq_dimid(ncid, "west_east", var_id), &
                     'static_init_model','inq_dimid west_east')
   call nc_check( nf90_inquire_dimension(ncid, var_id, name, we), &
                     'static_init_model','inquire_dimension '//trim(name))

   call nc_check( nf90_inq_dimid(ncid, "west_east_stag", var_id), &
                     'static_init_model','inq_dimid west_east_stag')  ! reuse var_id, no harm
   call nc_check( nf90_inquire_dimension(ncid, var_id, name, wes), &
                     'static_init_model','inquire_dimension '//trim(name))

   call nc_check( nf90_inq_dimid(ncid, "soil_layers_stag", var_id), &
                     'static_init_model','inq_dimid soil_layers_stag')  ! reuse var_id, no harm
   call nc_check( nf90_inquire_dimension(ncid, var_id, name, sls), &
                     'static_init_model','inquire_dimension '//trim(name))

   if(debug) then
      write(*,*) ' dimensions bt, sn, we are ',bt, &
           sn, we
      write(*,*) ' staggered  bt, sn, we are ',bts, &
           sns,wes
   endif

   RETURN

end subroutine read_wrf_dimensions

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine read_wrf_file_attributes(ncid,id)

! ncid: input, file handl
! id:   input, domain id

integer, intent(in)   :: ncid, id
logical, parameter    :: debug = .false.
integer :: ret

! get meta data and static data we need

   call nc_check( nf90_get_att(ncid, nf90_global, 'DX', wrf%dom(id)%dx), &
                     'static_init_model', 'get_att DX')
   call nc_check( nf90_get_att(ncid, nf90_global, 'DY', wrf%dom(id)%dy), &
                     'static_init_model', 'get_att DY')
   call nc_check( nf90_get_att(ncid, nf90_global, 'DT', wrf%dom(id)%dt), &
                     'static_init_model', 'get_att DT')
   write(errstring, *) 'dt from wrfinput_d0X file is: ', wrf%dom(id)%dt
   call error_handler(E_MSG, ' ', errstring)
   if(debug) write(*,*) ' dx, dy, dt are ',wrf%dom(id)%dx, &
        wrf%dom(id)%dy, wrf%dom(id)%dt

   call nc_check( nf90_get_att(ncid, nf90_global, 'MAP_PROJ', wrf%dom(id)%map_proj), &
                     'static_init_model', 'get_att MAP_PROJ')
   if(debug) write(*,*) ' map_proj is ',wrf%dom(id)%map_proj

   call nc_check( nf90_get_att(ncid, nf90_global, 'CEN_LAT', wrf%dom(id)%cen_lat), &
                     'static_init_model', 'get_att CEN_LAT')
   if(debug) write(*,*) ' cen_lat is ',wrf%dom(id)%cen_lat

   call nc_check( nf90_get_att(ncid, nf90_global, 'CEN_LON', wrf%dom(id)%cen_lon), &
                     'static_init_model', 'get_att CEN_LON')
   if(debug) write(*,*) ' cen_lon is ',wrf%dom(id)%cen_lon

   call nc_check( nf90_get_att(ncid, nf90_global, 'TRUELAT1', truelat1), &
                     'static_init_model', 'get_att TRUELAT1')
   if(debug) write(*,*) ' truelat1 is ',truelat1

   call nc_check( nf90_get_att(ncid, nf90_global, 'TRUELAT2', truelat2), &
                     'static_init_model', 'get_att TRUELAT2')
   if(debug) write(*,*) ' truelat2 is ',truelat2

   call nc_check( nf90_get_att(ncid, nf90_global, 'STAND_LON', stdlon), &
                     'static_init_model', 'get_att STAND_LON')
   if(debug) write(*,*) ' stdlon is ',stdlon

   ! this attribute is not present in older wrf files, which means it is not
   ! using a hybrid vertical coordinate.  not having this attribute should
   ! not cause an error.
   ret =  nf90_get_att(ncid, nf90_global, 'HYBRID_OPT', wrf%dom(id)%hybrid_opt)
   if (ret /= NF90_NOERR) wrf%dom(id)%hybrid_opt = 0
   if(debug) write(*,*) 'hybrid_opt is ', wrf%dom(id)%hybrid_opt

   RETURN

end subroutine read_wrf_file_attributes

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine assign_boundary_conditions(id)

! id:   input, domain id

integer, intent(in)   :: id
logical, parameter    :: debug = .false.

!nc -- fill in the boundary conditions (periodic_x and polar) here.  This will
!        need to be changed once these are taken from the NetCDF input instead
!        of the model namelist
!      NOTE :: because NetCDF cannot handle logicals, these boundary conditions
!        are likely to be read in as integers.  The agreed upon strategy is to 
!        test whether the integers are equal to 0 (for .false.) or 1 (for .true.)
!        and set the wrf%dom(id)% values to logicals to be used internally within
!        model_mod.f90.
!
!      Jeff Anderson points out that not everyone will convert to wrf3.0 and this
!        global attribute convention may not be backward-compatible.  So we should
!        test for existence of attributes and have defaults (from model_mod 
!        namelist) ready if they do not exist.  Note that defaults are currently 
!        true (as of 24 Oct 2007), but once the attributes arrive, the defaults
!        should be false.
   if ( id == 1 ) then
      wrf%dom(id)%periodic_x = periodic_x
      wrf%dom(id)%periodic_y = periodic_y
      wrf%dom(id)%polar = polar
      wrf%dom(id)%scm = scm
   else
      wrf%dom(id)%periodic_x = .false.
      wrf%dom(id)%periodic_y = .false.
      wrf%dom(id)%polar = .false.      
      wrf%dom(id)%scm = .false.      
   endif
   if(debug) write(*,*) ' periodic_x ',wrf%dom(id)%periodic_x
   if(debug) write(*,*) ' periodic_y ',wrf%dom(id)%periodic_y
   if(debug) write(*,*) ' polar ',wrf%dom(id)%polar
   if(debug) write(*,*) ' scm   ',wrf%dom(id)%scm   

end subroutine assign_boundary_conditions

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine read_wrf_static_data(ncid,id)

! ncid: input, file handle
! id:   input, domain id

integer, intent(in)   :: ncid, id
logical, parameter    :: debug = .false.
integer               :: var_id 

   call nc_check( nf90_inq_varid(ncid, "P_TOP", var_id), &
                     'read_wrf_static_data','inq_varid P_TOP')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%p_top), &
                     'read_wrf_static_data','get_var P_TOP')

!  get 1D (z) static data defining grid levels

   allocate(wrf%dom(id)%dn(1:wrf%dom(id)%bt))
   call nc_check( nf90_inq_varid(ncid, "DN", var_id), &
                     'read_wrf_static_data','inq_varid DN')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%dn), &
                     'read_wrf_static_data','get_var DN')
   if(debug) write(*,*) ' dn ',wrf%dom(id)%dn

   allocate(wrf%dom(id)%znu(1:wrf%dom(id)%bt))
   call nc_check( nf90_inq_varid(ncid, "ZNU", var_id), &
                     'read_wrf_static_data','inq_varid ZNU')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%znu), &
                     'read_wrf_static_data','get_var ZNU')
   if(debug) write(*,*) ' znu is ',wrf%dom(id)%znu

   allocate(wrf%dom(id)%znw(1:wrf%dom(id)%bts))
   call nc_check( nf90_inq_varid(ncid, "ZNW", var_id), &
                     'read_wrf_static_data','inq_varid ZNW')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%znw), &
                     'read_wrf_static_data','get_var ZNW')
   if(debug) write(*,*) ' znw is ',wrf%dom(id)%znw

   allocate(wrf%dom(id)%dnw(1:wrf%dom(id)%bt))
   call nc_check( nf90_inq_varid(ncid, "DNW", var_id), &
                     'read_wrf_static_data','inq_varid DNW')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%dnw), &
                     'read_wrf_static_data','get_var DNW')
   if(debug) write(*,*) ' dnw is ',wrf%dom(id)%dnw

   if (wrf%dom(id)%hybrid_opt == VERT_HYBRID) then
      allocate(wrf%dom(id)%c1h(1:wrf%dom(id)%bt))
      call nc_check( nf90_inq_varid(ncid, "C1H", var_id), &
                        'read_wrf_static_data','inq_varid C1H')
      call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%c1h), &
                        'read_wrf_static_data','get_var C1H')
      if(debug) write(*,*) ' c1h is ',wrf%dom(id)%c1h
   
      allocate(wrf%dom(id)%c2h(1:wrf%dom(id)%bt))
      call nc_check( nf90_inq_varid(ncid, "C2H", var_id), &
                        'read_wrf_static_data','inq_varid C2H')
      call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%c2h), &
                        'read_wrf_static_data','get_var C2H')
      if(debug) write(*,*) ' c2h is ',wrf%dom(id)%c2h
   endif

   allocate(wrf%dom(id)%zs(1:wrf%dom(id)%sls))
   call nc_check( nf90_inq_varid(ncid, "ZS", var_id), &
                     'read_wrf_static_data','inq_varid ZS')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%zs), &
                     'read_wrf_static_data','get_var ZS')

!  get 2D (x,y) base state for mu, latitude, longitude

   allocate(wrf%dom(id)%mub(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncid, "MUB", var_id), &
                     'read_wrf_static_data','inq_varid MUB')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%mub), &
                     'read_wrf_static_data','get_var MUB')
   if(debug) then
      write(*,*) ' corners of mub '
      write(*,*) wrf%dom(id)%mub(1,1),wrf%dom(id)%mub(wrf%dom(id)%we,1),  &
           wrf%dom(id)%mub(1,wrf%dom(id)%sn),wrf%dom(id)%mub(wrf%dom(id)%we, &
           wrf%dom(id)%sn)
   endif

   allocate(wrf%dom(id)%longitude(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncid, "XLONG", var_id), &
                     'read_wrf_static_data','inq_varid XLONG')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%longitude), &
                     'read_wrf_static_data','get_var XLONG')

   allocate(wrf%dom(id)%longitude_u(1:wrf%dom(id)%wes,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncid, "XLONG_U", var_id), &
                     'read_wrf_static_data','inq_varid XLONG_U')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%longitude_u), &
                     'read_wrf_static_data','get_var XLONG_U')

   allocate(wrf%dom(id)%longitude_v(1:wrf%dom(id)%we,1:wrf%dom(id)%sns))
   call nc_check( nf90_inq_varid(ncid, "XLONG_V", var_id), &
                     'read_wrf_static_data','inq_varid XLONG_V')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%longitude_v), &
                     'read_wrf_static_data','get_var XLONG_V')

   allocate(wrf%dom(id)%latitude(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncid, "XLAT", var_id), &
                     'read_wrf_static_data','inq_varid XLAT')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%latitude), &
                     'read_wrf_static_data','get_var XLAT')

   allocate(wrf%dom(id)%latitude_u(1:wrf%dom(id)%wes,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncid, "XLAT_U", var_id), &
                     'read_wrf_static_data','inq_varid XLAT_U')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%latitude_u), &
                     'read_wrf_static_data','get_var XLAT_U')

   allocate(wrf%dom(id)%latitude_v(1:wrf%dom(id)%we,1:wrf%dom(id)%sns))
   call nc_check( nf90_inq_varid(ncid, "XLAT_V", var_id), &
                     'read_wrf_static_data','inq_varid XLAT_V')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%latitude_v), &
                     'read_wrf_static_data','get_var XLAT_V')

   allocate(wrf%dom(id)%land(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncid, "XLAND", var_id), &
                     'read_wrf_static_data','inq_varid XLAND')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%land), &
                     'read_wrf_static_data','get_var XLAND')
   if(debug) then
      write(*,*) ' corners of land '
      write(*,*) wrf%dom(id)%land(1,1),wrf%dom(id)%land(wrf%dom(id)%we,1),  &
           wrf%dom(id)%land(1,wrf%dom(id)%sn),wrf%dom(id)%land(wrf%dom(id)%we, &
           wrf%dom(id)%sn)
   endif

   if(debug) then
      write(*,*) ' corners of lat '
      write(*,*) wrf%dom(id)%latitude(1,1),wrf%dom(id)%latitude(wrf%dom(id)%we,1),  &
           wrf%dom(id)%latitude(1,wrf%dom(id)%sn), &
           wrf%dom(id)%latitude(wrf%dom(id)%we,wrf%dom(id)%sn)
      write(*,*) ' corners of long '
      write(*,*) wrf%dom(id)%longitude(1,1),wrf%dom(id)%longitude(wrf%dom(id)%we,1),  &
           wrf%dom(id)%longitude(1,wrf%dom(id)%sn), &
           wrf%dom(id)%longitude(wrf%dom(id)%we,wrf%dom(id)%sn)
   endif

   allocate(wrf%dom(id)%hgt(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncid, "HGT", var_id), &
                     'read_wrf_static_data','inq_varid HGT')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%hgt), &
                     'read_wrf_static_data','get_var HGT')

! get 3D base state geopotential

   allocate(wrf%dom(id)%phb(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%bts))
   call nc_check( nf90_inq_varid(ncid, "PHB", var_id), &
                     'read_wrf_static_data','inq_varid PHB')
   call nc_check( nf90_get_var(ncid, var_id, wrf%dom(id)%phb), &
                     'read_wrf_static_data','get_var PHB')
   if(debug) then
      write(*,*) ' corners of phb '
      write(*,*) wrf%dom(id)%phb(1,1,1),wrf%dom(id)%phb(wrf%dom(id)%we,1,1),  &
           wrf%dom(id)%phb(1,wrf%dom(id)%sn,1),wrf%dom(id)%phb(wrf%dom(id)%we, &
           wrf%dom(id)%sn,1)
      write(*,*) wrf%dom(id)%phb(1,1,wrf%dom(id)%bts), &
           wrf%dom(id)%phb(wrf%dom(id)%we,1,wrf%dom(id)%bts),  &
           wrf%dom(id)%phb(1,wrf%dom(id)%sn,wrf%dom(id)%bts), &
           wrf%dom(id)%phb(wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bts)
   endif

end subroutine read_wrf_static_data

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine setup_map_projection(id)

! id:   input, domain id

integer, intent(in)   :: id
logical, parameter    :: debug = .false.

integer  :: proj_code
real(r8) :: latinc,loninc

! Initializes the map projection structure to missing values

   call map_init(wrf%dom(id)%proj)

! Populate the map projection structure

!nc -- added in case structures for CASSINI and CYL
!nc -- global wrfinput_d01 has truelat1 = 1.e20, so we need to change it where needed
!nc -- for PROJ_LATLON stdlon and truelat1 have different meanings -- 
!nc --   stdlon --> loninc  and  truelat1 --> latinc
!JPH --- this latinc/loninc calculations are only valid for global domains

   if ( wrf%dom(id)%scm ) then
! JPH -- set to zero which should cause the map utils to return NaNs if called
      latinc = 0.0_r8 
      loninc = 0.0_r8 
   else
      latinc = 180.0_r8/wrf%dom(id)%sn
      loninc = 360.0_r8/wrf%dom(id)%we
   endif

   if(wrf%dom(id)%map_proj == map_latlon) then
      truelat1 = latinc
      stdlon = loninc
      proj_code = PROJ_LATLON
   elseif(wrf%dom(id)%map_proj == map_lambert) then
      proj_code = PROJ_LC
   elseif(wrf%dom(id)%map_proj == map_polar_stereo) then
      proj_code = PROJ_PS
   elseif(wrf%dom(id)%map_proj == map_mercator) then
      proj_code = PROJ_MERC
   elseif(wrf%dom(id)%map_proj == map_cyl) then
      proj_code = PROJ_CYL
   elseif(wrf%dom(id)%map_proj == map_cassini) then
      proj_code = PROJ_CASSINI
   else
      call error_handler(E_ERR,'static_init_model', &
        'Map projection no supported.', source, revision, revdate)
   endif

!nc -- specified inputs to hopefully handle ALL map projections -- hopefully map_set will
!        just ignore the inputs it doesn't need for its map projection of interest (?)
!     
!   NOTE:: We are NOT yet supporting the Gaussian grid or the Rotated Lat/Lon, so we
!            are going to skip the entries:  nlon, nlat, ixdim, jydim, stagger, phi, lambda
!
!      + Gaussian grid uses nlat & nlon
!      + Rotated Lat/Lon uses ixdim, jydim, stagger, phi, & lambda
!
   call map_set( proj_code=proj_code, &
                 proj=wrf%dom(id)%proj, &
                 lat1=wrf%dom(id)%latitude(1,1), &
                 lon1=wrf%dom(id)%longitude(1,1), &
                 lat0=90.0_r8, &
                 lon0=0.0_r8, &
                 knowni=1.0_r8, &
                 knownj=1.0_r8, &
                 dx=wrf%dom(id)%dx, &
                 latinc=latinc, &
                 loninc=loninc, &
                 stdlon=stdlon, &
                 truelat1=truelat1, &
                 truelat2=truelat2  )

end subroutine setup_map_projection

!--------------------------------------------
!--------------------------------------------

subroutine fill_default_state_table(default_table)

character(len=129), intent(out) :: default_table(num_state_table_columns,max_state_variables) 

integer :: row

default_table = 'NULL'
row = 0

! fill default state variable table here.
row = row+1
default_table(:,row) = (/ 'U                         ', &
                          'QTY_U_WIND_COMPONENT      ', &
                          'TYPE_U                    ', &
                          'UPDATE                    ', &
                          '999                       '  /)
row = row+1
default_table(:,row) = (/ 'V                         ', &
                          'QTY_V_WIND_COMPONENT      ', &
                          'TYPE_V                    ', &
                          'UPDATE                    ', &
                          '999                       '  /)
row = row+1
default_table(:,row) = (/ 'W                         ', &
                          'QTY_VERTICAL_VELOCITY     ', &
                          'TYPE_W                    ', &
                          'UPDATE                    ', &
                          '999                       '  /)
row = row+1
default_table(:,row) = (/ 'PH                        ', &
                          'QTY_GEOPOTENTIAL_HEIGHT   ', &
                          'TYPE_GZ                   ', &
                          'UPDATE                    ', &
                          '999                       '  /)
row = row+1
default_table(:,row) = (/ 'THM                       ', &
                          'QTY_POTENTIAL_TEMPERATURE ', &
                          'TYPE_T                    ', &
                          'UPDATE                    ', &
                          '999                       '  /)
row = row+1
default_table(:,row) = (/ 'MU                        ', &
                          'QTY_PRESSURE              ', &
                          'TYPE_MU                   ', &
                          'UPDATE                    ', &
                          '999                       '  /)

row = row+1
default_table(:,row) = (/ 'QVAPOR                    ', &
                          'QTY_VAPOR_MIXING_RATIO    ', &
                          'TYPE_QV                   ', &
                          'UPDATE                    ', &
                          '999                       '  /)

return

end subroutine fill_default_state_table

!--------------------------------------------
!--------------------------------------------

subroutine fill_dart_kinds_table(wrf_state_variables, in_state_vector)

! for each row in the kinds table, tick them off in an array
! of all possible kinds.  then do some simple error checks for
! kinds we know have interactions -- like both wind vectors are
! required to convert directions from projection to lat/lon, etc

character(len=*), intent(in) :: wrf_state_variables(:, :)
logical, intent(inout)       :: in_state_vector(:)

integer :: row, i, nextkind

in_state_vector = .false.

row = size(wrf_state_variables, 2)

! NEWVAR: for a simple new variable you do not need to do anything here.
! NEWVAR: but if a new kind has interactions - like if you have an obs
! NEWVAR: of this kind, you actually interpolate in multiple fields in
! NEWVAR: the state vector (e.g. wind needs both U and V), then you
! NEWVAR: might need to add some code here.  
! NEWVAR: see each of part1, part 2, and part 3 below.

! part 1: mark off all the kinds that the user specifies, plus the
! kinds that are related and can be interpolated from the given kind.

do i = 1, row

   ! end of the list?
   if (wrf_state_variables(2, i) == 'NULL') exit

   nextkind = get_index_for_quantity(trim(wrf_state_variables(2, i)))
   select case(nextkind)

   ! wrf stores potential temperature (temperature perturbations around a 
   ! threshold) but we can interpolate sensible temperature from it
   case (QTY_POTENTIAL_TEMPERATURE)
      in_state_vector(QTY_TEMPERATURE) = .true.
      in_state_vector(QTY_POTENTIAL_TEMPERATURE) = .true.

   ! we use vapor mixing ratio to compute specific humidity
   case (QTY_VAPOR_MIXING_RATIO)
      in_state_vector(QTY_VAPOR_MIXING_RATIO) = .true.
      in_state_vector(QTY_SPECIFIC_HUMIDITY) = .true.

   ! unrecognized kind string in namelist.
   ! 0 is actually QTY_STATE_VARIABLE and not supported here.
   case (-1, 0) 
      write(errstring, *) 'unrecognized KIND string: ' // trim(wrf_state_variables(2, i))
      call error_handler(E_ERR, 'fill_dart_kinds_table', errstring, &
                         source, revision, revdate)
     
   ! everything else, set it to be supported
   case default
      in_state_vector(nextkind) = .true.

   end select

   ! NOTE: PSFC can be labeled either QTY_PRESSURE or QTY_SURFACE_PRESSURE
   ! in the namelist, but make sure however it is labeled that for now we 
   ! allow surface pressure interpolation.  this may go away once we work out
   ! QTY_FOO vs QTY_SURFACE_FOO - are they fundamentally different things
   ! or should the decision be made based on a QTY_FOO and the vertical
   ! location type -- if it is VERTISSURFACE, then you do the 2d calc in the
   ! surface field, otherwise you do the full-up 3d interpolation.
   if ( wrf_state_variables(1, i) == 'PSFC' ) then
      in_state_vector(QTY_SURFACE_PRESSURE) = .true.
   endif

enddo


! part 2: if you specified one kind but the interpolation is going to
! require another, make sure the combinations make sense.  i.e. if you
! have U wind, you also have to have V wind, etc.

do i = 1, size(in_state_vector)

   select case(i)
  
   ! the vortex center computations require wind speeds and phb?
   case (QTY_VORTEX_LAT, QTY_VORTEX_LON, QTY_VORTEX_PMIN, QTY_VORTEX_WMAX)
      if ((.not. in_state_vector(QTY_U_WIND_COMPONENT))   .or. &
          (.not. in_state_vector(QTY_V_WIND_COMPONENT))   .or. &
          (.not. in_state_vector(QTY_TEMPERATURE))        .or. &
          (.not. in_state_vector(QTY_VAPOR_MIXING_RATIO)) .or. &
          (.not. in_state_vector(QTY_PRESSURE))) then
         write(errstring, *) 'VORTEX kinds will require U,V,T,QVAPOR,MU in state vector'
         ! FIXME: not fatal error, just informative at this point.
         call error_handler(E_MSG, 'fill_dart_kinds_table', errstring, &
                            source, revision, revdate)
      endif
 
   ! if you have one wind component you have to have both
   case (QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT)
      if (in_state_vector(QTY_U_WIND_COMPONENT) .neqv. in_state_vector(QTY_V_WIND_COMPONENT)) then
         write(errstring, *) 'WIND kinds will require both U,V in state vector'
         ! FIXME: not fatal error, just informative at this point.
         call error_handler(E_MSG, 'fill_dart_kinds_table', errstring, &
                            source, revision, revdate)
      endif
 
   ! by default anything else is fine

   end select
enddo


! part 3: fields you just have to have, always, and other exceptions
! and things that break the rules.
if (.not. in_state_vector(QTY_GEOPOTENTIAL_HEIGHT)) then
   write(errstring, *) 'PH is always a required field'
   call error_handler(E_ERR, 'fill_dart_kinds_table', errstring, &
                      source, revision, revdate)
endif

! FIXME: is this true?  or is pressure always required, and surface
! pressure required only if you have any of the surface obs?
if ((.not. in_state_vector(QTY_PRESSURE)) .and. &
    (.not. in_state_vector(QTY_SURFACE_PRESSURE))) then
    write(errstring, *) 'One of MU or PSFC is a required field'
    call error_handler(E_ERR, 'fill_dart_kinds_table', errstring, &
                       source, revision, revdate)
endif

! surface elevation and land mask are read outside the state vector mechanism,
! directly from the wrfinput template file, and does not vary from
! one ensemble member to another.
in_state_vector(QTY_SURFACE_ELEVATION) = .true.
in_state_vector(QTY_LANDMASK) = .true.
in_state_vector(QTY_SURFACE_TYPE) = .true.

! there is no field that directly maps to the vortex measurements.
! if you have all the fields it needs, allow them.
if (in_state_vector(QTY_U_WIND_COMPONENT)    .and. &
    in_state_vector(QTY_V_WIND_COMPONENT)    .and. &
    in_state_vector(QTY_TEMPERATURE)         .and. &
    in_state_vector(QTY_VAPOR_MIXING_RATIO)  .and. &
    in_state_vector(QTY_PRESSURE)) then        ! ok to add vortex types
   in_state_vector(QTY_VORTEX_LAT)  = .true.
   in_state_vector(QTY_VORTEX_LON)  = .true.
   in_state_vector(QTY_VORTEX_PMIN) = .true.
   in_state_vector(QTY_VORTEX_WMAX) = .true.
endif
 
! if you have geopotential height and pressure, you can compute
! a density value.
if (in_state_vector(QTY_GEOPOTENTIAL_HEIGHT) .and. &
    in_state_vector(QTY_PRESSURE) ) in_state_vector(QTY_DENSITY) = .true.


! allow reflectivity to be asked for, because the obs_def has an alternative
! way to interpolate it even if it is not in the state vector. if this
! is not allowed it will error exit instead of returning with an invalid istatus
! to indicate to the caller that the interpolation failed.
! ditto for power weighted fall speed.
in_state_vector(QTY_RADAR_REFLECTIVITY) = .true.
in_state_vector(QTY_POWER_WEIGHTED_FALL_SPEED) = .true.



! FIXME:  i was going to suggest nuking this routine all together because it makes
! the default behavior be to exit with an error when requesting to interpolate an
! unknown kind, instead of simply returning with a failed interpolation and letting
! the calling code handle it.  but it does turn out to be amazingly useful when
! trying to interpolate something you do expect to be there, so now i'm unsure what
! to suggest here.    nsc  30jul09

end subroutine fill_dart_kinds_table

!--------------------------------------------
!--------------------------------------------

integer function get_number_of_wrf_variables(id, state_table, var_element_list, var_update_list)

integer, intent(in) :: id
character(len=*), intent(in) :: state_table(num_state_table_columns,max_state_variables) 
integer, intent(out), optional :: var_element_list(max_state_variables)
logical, intent(out), optional :: var_update_list(max_state_variables)

integer :: ivar, num_vars
logical :: debug = .false.

if ( present(var_element_list) ) var_element_list = -1

ivar = 1
num_vars = 0
do while ( trim(state_table(5,ivar)) /= 'NULL' ) 

   if ( variable_is_on_domain(state_table(5,ivar),id) ) then
      num_vars = num_vars + 1
      if ( present(var_element_list) ) var_element_list(num_vars) = ivar

      if (present(var_update_list)) then
         if (state_table(4,ivar) == 'NO_COPY_BACK') then
            var_update_list(num_vars) = .false.
         else
            var_update_list(num_vars) = .true.
         endif
      endif

   endif

   ivar = ivar + 1

enddo ! ivar

if ( debug ) then
  print*,'In function get_number_of_wrf_variables'
  print*,'Found ',num_vars,' state variables for domain id ',id
endif

get_number_of_wrf_variables = num_vars 

return

end function get_number_of_wrf_variables 

!--------------------------------------------
!--------------------------------------------

subroutine set_variable_bound_defaults(nbounds,lb,ub,instructions)

   integer, intent(in)                                :: nbounds
   real(r8), dimension(nbounds), intent(out)          :: lb, ub
   character(len=10), dimension(nbounds), intent(out) :: instructions

   lb(:) = missing_r8
   ub(:) = missing_r8
   instructions(:) = 'NULL'

   return

end subroutine set_variable_bound_defaults

!--------------------------------------------
!--------------------------------------------

subroutine get_variable_bounds(bounds_table,wrf_var_name,lb,ub,instructions)

! matches WRF variable name in bounds table to input name, and assigns
! the bounds and instructions if they exist

   character(len=*), intent(in)    :: bounds_table(num_bounds_table_columns,max_state_variables) 
   character(len=*), intent(in)    :: wrf_var_name
   real(r8),         intent(out)   :: lb,ub
   character(len=10),intent(out)   :: instructions

   character(len=50)               :: wrf_varname_trim, bounds_varname_trim
   character(len=50)               :: bound_trim
   integer :: ivar
   logical :: debug = .false.

   wrf_varname_trim = ''
   wrf_varname_trim = trim(wrf_var_name)

   ivar = 1
   do while ( trim(bounds_table(1,ivar)) /= 'NULL' ) 

     bounds_varname_trim = trim(bounds_table(1,ivar))

     if ( bounds_varname_trim == wrf_varname_trim ) then
 
        bound_trim = trim(bounds_table(2,ivar))
        if ( bound_trim  /= 'NULL' ) then
           read(bound_trim,'(d16.8)') lb 
        else
           lb = missing_r8
        endif

        bound_trim = trim(bounds_table(3,ivar))
        if ( bound_trim  /= 'NULL' ) then
           read(bound_trim,'(d16.8)') ub 
        else
           ub = missing_r8
        endif

        ! instructions are required
        instructions = trim(bounds_table(4,ivar))

        if ( instructions == 'NULL' ) then
          call error_handler(E_ERR,'get_variable_bounds','instructions for bounds on '//wrf_varname_trim//' are required',&
                                    source,revision,revdate)
        endif
 
     endif

     ivar = ivar + 1

   enddo !ivar

   if ( debug ) then
      write(*,*) 'In get_variable_bounds assigned ',wrf_varname_trim
      write(*,*) ' bounds ',lb,ub,instructions
   endif

   return

end subroutine get_variable_bounds

!--------------------------------------------
!--------------------------------------------

logical function variable_is_on_domain(domain_id_string, id)

integer,           intent(in) :: id
character(len=*),  intent(in) :: domain_id_string

integer                       :: domain_int, i

variable_is_on_domain = .false.

   ! if '999' then counts all domains
   if ( trim(domain_id_string) == '999' ) then
      variable_is_on_domain = .true.
   else
      i = 1
      do while ( domain_id_string(i:i) /= ' ' )
         read(domain_id_string(i:i),'(i1)') domain_int
         if ( domain_int == id ) variable_is_on_domain = .true.
         i = i+1
      enddo
   endif

return

end function variable_is_on_domain

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine get_variable_size_from_file(ncid,id,wrf_var_name,bt,bts,sn,sns, &
                                       we,wes,stagger,var_size)

!NOTE: only supports 2D and 3D variables (ignoring time dimension)

! ncid: input, file handle
! id:   input, domain index

integer,           intent(in)   :: ncid, id
integer,           intent(in)   :: bt, bts, sn, sns, we, wes
character(len=*),  intent(in)   :: wrf_var_name
integer,           intent(out)  :: var_size(3)
character(len=129),intent(out)  :: stagger

logical, parameter    :: debug = .false.
integer               :: var_id, ndims, dimids(10) 
integer               :: idim

   stagger = ''

! get variable ID
   call nc_check( nf90_inq_varid(ncid, trim(wrf_var_name), var_id), &
                     'get_variable_size_from_file',                 &
                     'inq_varid '//wrf_var_name)

! get number of dimensions and dimension IDs
   call nc_check( nf90_inquire_variable(ncid, var_id,ndims=ndims,  &
                     dimids=dimids),                               &
                     'get_variable_size_from_file',                &
                     'inquire_variable '//wrf_var_name)

! get dimension length, ignoring first dimension (time)
   do idim = 1,ndims-1
      call nc_check( nf90_inquire_dimension(ncid, dimids(idim),  &
                     len=var_size(idim)),                               &
                     'get_variable_size_from_file',                &
                     'inquire_dimension '//wrf_var_name)
   enddo

! if a 2D variable fill the vertical dimension with 1
   if ( ndims < 4 ) var_size(ndims) = 1

   if ( debug ) then
      print*,'In get_variable_size_from_file got variable size ',var_size
   endif


! get variable attribute stagger
   call nc_check( nf90_get_att(ncid, var_id, 'stagger', stagger), &
                     'get_variable_size_from_file', &
                     'get_att '//wrf_var_name//' '//stagger)

   if ( debug ) then
      print*,'In get_variable_size_from_file got stagger ',trim(stagger),' for variable ',trim(wrf_var_name)
   endif

return

end subroutine get_variable_size_from_file

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine get_variable_metadata_from_file(ncid,id,wrf_var_name,description, &
                                       coordinates,units)

! ncid: input, file handle
! id:   input, domain index

integer, intent(in)               :: ncid, id
character(len=*),   intent(in)    :: wrf_var_name
character(len=129), intent(out)   :: description, coordinates, units

logical, parameter    :: debug = .false.
integer               :: var_id

   call nc_check( nf90_inq_varid(ncid, trim(wrf_var_name), var_id), &
                     'get_variable_metadata_from_file', &
                     'inq_varid '//wrf_var_name)

   description = ''
   call nc_check( nf90_get_att(ncid, var_id, 'description', description), &
                     'get_variable_metadata_from_file', &
                     'get_att '//wrf_var_name//' '//description)

   coordinates = ''
   call nc_check( nf90_get_att(ncid, var_id, 'coordinates', coordinates), &
                     'get_variable_metadata_from_file', &
                     'get_att '//wrf_var_name//' '//coordinates)

   units = ''
   call nc_check( nf90_get_att(ncid, var_id, 'units', units), &
                     'get_variable_metadata_from_file', &
                     'get_att '//wrf_var_name//' '//units)

return

end subroutine get_variable_metadata_from_file

!--------------------------------------------
!--------------------------------------------
! Note get_dart_vector_index depends on this function
integer function get_type_ind_from_type_string(id, wrf_varname)

! simply loop through the state variable table to get the index of the
! type for this domain.  returns -1 if not found

   integer,          intent(in) :: id
   character(len=*), intent(in) :: wrf_varname

   integer                      :: ivar, my_index
   logical                      :: debug = .false.
   character(len=50)            :: wrf_varname_trim, wrf_state_var_trim

   get_type_ind_from_type_string = -1

   do ivar = 1,wrf%dom(id)%number_of_wrf_variables
      
      my_index =  wrf%dom(id)%var_index_list(ivar)

      wrf_state_var_trim = trim(wrf_state_variables(1,my_index))
      wrf_varname_trim   = trim(wrf_varname)
      
      if ( wrf_state_var_trim == wrf_varname_trim ) then

         get_type_ind_from_type_string = ivar

      endif

   enddo ! ivar

   if ( debug ) write(*,*) 'get_type_from_ind ',trim(wrf_varname),' ',get_type_ind_from_type_string

   return

end function get_type_ind_from_type_string

!----------------------------------------------------------------------
! Returns integers taken from tstring
! It is assumed that the tstring char array is as YYYY-MM-DD_hh:mm:ss

subroutine get_wrf_date (tstring, year, month, day, hour, minute, second)

integer,           intent(out) :: year, month, day, hour, minute, second
character(len=19), intent(in)  :: tstring

read(tstring( 1: 4),'(i4)') year
read(tstring( 6: 7),'(i2)') month
read(tstring( 9:10),'(i2)') day
read(tstring(12:13),'(i2)') hour
read(tstring(15:16),'(i2)') minute
read(tstring(18:19),'(i2)') second

return

end subroutine get_wrf_date

!----------------------------------------------------------------------
! Returns integers taken from tstring
! It is assumed that the tstring char array is as YYYY-MM-DD_hh:mm:ss

subroutine set_wrf_date (tstring, year, month, day, hour, minute, second)

integer,           intent(in) :: year, month, day, hour, minute, second
character(len=19), intent(out)  :: tstring

character(len=4)  :: ch_year
character(len=2)  :: ch_month, ch_day, ch_hour, ch_minute, ch_second

write(ch_year,'(i4)') year
write(ch_month,'(i2)') month
if (ch_month(1:1) == " ") ch_month(1:1) = "0"
write(ch_day,'(i2)') day
if (ch_day(1:1) == " ") ch_day(1:1) = "0"
write(ch_hour,'(i2)') hour
if (ch_hour(1:1) == " ") ch_hour(1:1) = "0"
write(ch_minute,'(i2)') minute
if (ch_minute(1:1) == " ") ch_minute(1:1) = "0"
write(ch_second,'(i2)') second
if (ch_second(1:1) == " ") ch_second(1:1) = "0"
tstring(1:4)   = ch_year
tstring(5:5)   = "-"
tstring(6:7)   = ch_month
tstring(8:8)   = "-"
tstring(9:10)  = ch_day
tstring(11:11) = "_"
tstring(12:13) = ch_hour
tstring(14:14) = ":"
tstring(15:16) = ch_minute
tstring(17:17) = ":"
tstring(18:19) = ch_second

end subroutine set_wrf_date


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    height_diff_check - function that determines whether a pair of heights
!                        in meters are closer than the given tolerance.
!                        returns .TRUE. if closer or equal to limit
!
!    max_diff_meters   - maximum difference between 2 elevations (m)
!    height1           - first height (m)
!    height2           - second height (m)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function height_diff_check(max_diff_meters, height1, height2)
 real(r8), intent(in) :: max_diff_meters, height1, height2
 logical              :: height_diff_check

height_diff_check = .true.

if ( abs(height1 - height2) > max_diff_meters ) height_diff_check = .false.

end function height_diff_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    compute_geometric_ht - subroutine the converts geopotential height
!                        into geometric height.
!
!     geopot -- input real value geopotential height
!     lat    -- latitude of input (longitude not needed)
!
!     output is real value geometric height
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function compute_geometric_height(geopot, lat)
 real(r8), intent(in)  :: geopot
 real(r8), intent(in)  :: lat
 real(r8)              :: compute_geometric_height


! -----------------------------------------------------------------------*/
   real(digits12) :: pi2, latr
   real(digits12) :: semi_major_axis, semi_minor_axis, grav_polar, grav_equator
   real(digits12) :: earth_omega, grav_constant, flattening, somigliana
   real(digits12) :: grav_ratio, sin2, termg, termr, grav, eccentricity

!  Parameters below from WGS-84 model software inside GPS receivers.
   parameter(semi_major_axis = 6378.1370d3)    ! (m)
   parameter(semi_minor_axis = 6356.7523142d3) ! (m)
   parameter(grav_polar = 9.8321849378)        ! (m/s2)
   parameter(grav_equator = 9.7803253359)      ! (m/s2)
   parameter(earth_omega = 7.292115d-5)        ! (rad/s)
   parameter(grav_constant = 3.986004418d14)   ! (m3/s2)
   parameter(grav = 9.80665d0)                 ! (m/s2) WMO std g at 45 deg lat
   parameter(eccentricity = 0.081819d0)        ! unitless
   parameter(pi2 = 3.14159265358979d0/180.d0)

!  Derived geophysical constants
   parameter(flattening = (semi_major_axis-semi_minor_axis) / semi_major_axis)

   parameter(somigliana = (semi_minor_axis/semi_major_axis)*(grav_polar/grav_equator)-1.d0)

   parameter(grav_ratio = (earth_omega*earth_omega * &
                semi_major_axis*semi_major_axis * semi_minor_axis)/grav_constant)


   latr = lat * (pi2)        ! in radians
   sin2  = sin(latr) * sin(latr)
   termg = grav_equator * ( (1.d0+somigliana*sin2) / &
           sqrt(1.d0-eccentricity*eccentricity*sin2) )
   termr = semi_major_axis / (1.d0 + flattening + grav_ratio - 2.d0*flattening*sin2)

   ! geometric height conversion.  in the special case that you're trying to
   ! reproduce an old result and want to use the original geopotential height
   ! instead of converting it, set 'use_geopotential_height' to .true. at the
   ! top of this file.  Otherwise the normal conversion to geometric height is
   ! what you want here; we agreed that height was going to be in terms of
   ! geometric height when we did the conversions.

   if (use_geopotential_height) then
      compute_geometric_height = geopot
   else
      compute_geometric_height = (termr*geopot) / ( (termg/grav) * termr - geopot )
   endif

end function compute_geometric_height

!--------------------------------------------------------------------------
!> Perform interpolation across the ensemble
subroutine simple_interp_distrib(fld, wrf, id, i, j, k, obs_kind, dxm, dx, dy, dym, uniquek, ens_size, state_handle)

integer,             intent(in) :: ens_size
type(ensemble_type), intent(in) :: state_handle
type(wrf_dom),       intent(in) :: wrf
integer,             intent(in) :: id
integer,             intent(in) :: i,j
integer,             intent(in) :: uniquek(:)
integer,             intent(in) :: k(ens_size)
integer,             intent(in) :: obs_kind
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8),           intent(out) :: fld(2, ens_size)

integer(i8)           :: ill, iul, ilr, iur
integer               :: rc
integer, dimension(2) :: ll, ul, lr, ur
integer               :: uk, e
logical               :: in_state
integer               :: wrf_type

real(r8), dimension(ens_size) ::x_ill, x_iul, x_ilr, x_iur

! Confirm that the obs kind is in the DART state vector and return the wrf_type
!>@todo should boundsCheck always be temperatue type? This is what it is in the original code
call obs_kind_in_state_vector(in_state, wrf_type, obs_kind, id)

if ( in_state ) then

   ! An observation could be on different levels for each ensemble member. 
   ! But, you don't want to do ens_size*communication, so just do it for the levels you need.
   UNIQUEK_LOOP: do uk = 1, size(uniquek)

   ! Check to make sure retrieved integer gridpoints are in valid range
   if ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf_type ) .and.&
      boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf_type ) .and. &
      boundsCheck( uniquek(uk), .false.,                id, dim=3, type=wrf_type ) ) then
         
         call getCorners(i, j, id, wrf_type, ll, ul, lr, ur, rc )
         if ( rc .ne. 0 ) &
         print*, 'model_mod.f90 :: model_interpolate :: getCorners QNSNOW rc = ', rc
               
         ! Interpolation for QNSNOW field at level k
         ill = get_dart_vector_index(ll(1), ll(2), uniquek(uk), domain_id(id), wrf_type)
         iul = get_dart_vector_index(ul(1), ul(2), uniquek(uk), domain_id(id), wrf_type)
         ilr = get_dart_vector_index(lr(1), lr(2), uniquek(uk), domain_id(id), wrf_type)
         iur = get_dart_vector_index(ur(1), ur(2), uniquek(uk), domain_id(id), wrf_type)

         x_ill = get_state(ill, state_handle)
         x_iul = get_state(iul, state_handle)
         x_ilr = get_state(ilr, state_handle)
         x_iur = get_state(iur, state_handle)

         do e = 1, ens_size
            if ( k(e) == uniquek(uk) ) then
               fld(1,e) = dym*( dxm*x_ill(e) + dx*x_ilr(e) ) + dy*( dxm*x_iul(e) + dx*x_iur(e) )
            endif
         enddo

         ! Interpolation for QNSNOW field at level k+1
         ill = get_dart_vector_index(ll(1), ll(2), uniquek(uk)+1, domain_id(id), wrf_type)
         iul = get_dart_vector_index(ul(1), ul(2), uniquek(uk)+1, domain_id(id), wrf_type)
         ilr = get_dart_vector_index(lr(1), lr(2), uniquek(uk)+1, domain_id(id), wrf_type)
         iur = get_dart_vector_index(ur(1), ur(2), uniquek(uk)+1, domain_id(id), wrf_type)

         x_ill = get_state(ill, state_handle)
         x_iul = get_state(iul, state_handle)
         x_ilr = get_state(ilr, state_handle)
         x_iur = get_state(iur, state_handle)

         do e = 1, ens_size
            if ( k(e) == uniquek(uk) ) then
               fld(2, e) = dym*( dxm*x_ill(e) + dx*x_ilr(e) ) + dy*( dxm*x_iul(e) + dx*x_iur(e) )
            endif
         enddo

      endif

   enddo UNIQUEK_LOOP

else ! not in state

   call error_handler(E_MSG, 'simple_interp_distrib', &
      'obs_kind "'//trim(get_name_for_quantity(obs_kind))//'" is not in state vector 1', &
      source, revision, revdate)
   fld(2, ens_size) = missing_r8
   return

endif

end subroutine simple_interp_distrib

!------------------------------------------------------------------------
!> interpolation for surface fields
subroutine surface_interp_distrib(fld, wrf, id, i, j, obs_kind, wrf_surf_type, dxm, dx, dy, dym, ens_size, state_handle)

integer,             intent(in) :: ens_size
type(ensemble_type), intent(in) :: state_handle
type(wrf_dom),       intent(in) :: wrf
integer,             intent(in) :: id
integer,             intent(in) :: obs_kind
integer,             intent(in) :: wrf_surf_type !< this has to be passed in is beccause obs_kind does not give a unique wrf_type ( could check obs_kind then if surf var?)
integer,             intent(in) :: i,j
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8),           intent(out) :: fld(2, ens_size)

integer(i8)           :: ill, iul, ilr, iur
integer               :: rc
integer, dimension(2) :: ll, ul, lr, ur

logical               :: in_state
integer               :: wrf_type

real(r8), dimension(ens_size) ::x_ill, x_iul, x_ilr, x_iur

fld(:,:) = missing_r8

! Find the wrf_type from the obs kind
! check for in state is performed before surface_interp_distrib is called
!>@todo should boundsCheck always be temperatue type? This is what it is in the original code
call obs_kind_in_state_vector(in_state, wrf_type, obs_kind, id)

! Check to make sure retrieved integer gridpoints are in valid range
if ( ( boundsCheck( i, wrf%dom(id)%periodic_x, id, dim=1, type=wrf_type ) .and. &
       boundsCheck( j, wrf%dom(id)%polar,      id, dim=2, type=wrf_type ) ) &
       .or. wrf%dom(id)%scm ) then

  call getCorners(i, j, id, wrf_type, ll, ul, lr, ur, rc )
  if ( rc .ne. 0 ) &
     print*, 'model_mod.f90 :: model_interpolate :: getCorners T2 rc = ', rc
   
     ! Interpolation for the T2 field
     ill = get_dart_vector_index(ll(1), ll(2), 1, domain_id(id), wrf_surf_type)
     iul = get_dart_vector_index(ul(1), ul(2), 1, domain_id(id), wrf_surf_type)
     ilr = get_dart_vector_index(lr(1), lr(2), 1, domain_id(id), wrf_surf_type)
     iur = get_dart_vector_index(ur(1), ur(2), 1, domain_id(id), wrf_surf_type)

     x_ill = get_state(ill, state_handle)
     x_iul = get_state(iul, state_handle)
     x_iur = get_state(iur, state_handle)
     x_ilr = get_state(ilr, state_handle)

     fld(1, :) = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )

endif

end subroutine surface_interp_distrib

!--------------------------------------------------------------------------
!> test if an obs kind is in the state vector and set wrf_type
subroutine obs_kind_in_state_vector(part_of_state_vector, wrf_type, obs_kind, id)

logical, intent(out) :: part_of_state_vector
integer, intent(out) :: wrf_type !< WRF
integer, intent(in)  :: obs_kind !< DART
integer, intent(in)  :: id


part_of_state_vector = .false. ! assume not in state vector


if    ( ( obs_kind == QTY_VERTICAL_VELOCITY)             .and. ( wrf%dom(id)%type_w >= 0 ) )  then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_w
else if( ( obs_kind == QTY_RAINWATER_MIXING_RATIO )      .and. ( wrf%dom(id)%type_qr >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_qr
else if( ( obs_kind == QTY_GRAUPEL_MIXING_RATIO )        .and. ( wrf%dom(id)%type_qg >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_qg
else if( ( obs_kind == QTY_HAIL_MIXING_RATIO )           .and. ( wrf%dom(id)%type_qh >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_qh
else if( ( obs_kind == QTY_SNOW_MIXING_RATIO )           .and. ( wrf%dom(id)%type_qs >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_qs
else if( ( obs_kind == QTY_ICE_MIXING_RATIO )                   .and. ( wrf%dom(id)%type_qi >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_qi
else if( ( obs_kind == QTY_CLOUDWATER_MIXING_RATIO )          .and. ( wrf%dom(id)%type_qc >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_qc
else if ( ( obs_kind == QTY_CLOUD_FRACTION )             .and. ( wrf%dom(id)%type_cldfra >= 0 ) )then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_cldfra
else if( ( obs_kind == QTY_DROPLET_NUMBER_CONCENTR )     .and. ( wrf%dom(id)%type_qndrp >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_qndrp
else if( ( obs_kind == QTY_ICE_NUMBER_CONCENTRATION )    .and. ( wrf%dom(id)%type_qnice >= 0 ) )then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_qnice
else if( ( obs_kind == QTY_SNOW_NUMBER_CONCENTR )        .and. ( wrf%dom(id)%type_qnsnow >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_qnsnow
else if( ( obs_kind == QTY_RAIN_NUMBER_CONCENTR )        .and. ( wrf%dom(id)%type_qnrain >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_qnrain
else if( ( obs_kind == QTY_GRAUPEL_NUMBER_CONCENTR )     .and. ( wrf%dom(id)%type_qngraupel >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_qngraupel
else if( ( obs_kind == QTY_HAIL_NUMBER_CONCENTR )        .and. ( wrf%dom(id)%type_qnhail >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_qnhail
else if( ( obs_kind == QTY_CONDENSATIONAL_HEATING )      .and. ( wrf%dom(id)%type_hdiab >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_hdiab
else if( ( obs_kind == QTY_POWER_WEIGHTED_FALL_SPEED )   .and. ( wrf%dom(id)%type_fall_spd >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_fall_spd
else if( ( obs_kind == QTY_RADAR_REFLECTIVITY )          .and. ( wrf%dom(id)%type_refl >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_refl
else if( ( obs_kind == QTY_DIFFERENTIAL_REFLECTIVITY )   .and. ( wrf%dom(id)%type_dref >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type =  wrf%dom(id)%type_dref
else if( ( obs_kind == QTY_SPECIFIC_DIFFERENTIAL_PHASE ) .and. ( wrf%dom(id)%type_spdp >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_spdp
else if ( ( obs_kind == QTY_VAPOR_MIXING_RATIO )         .and. ( wrf%dom(id)%type_qv >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_qv
else if ( ( obs_kind == QTY_TEMPERATURE )                  .and. ( wrf%dom(id)%type_t >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_t
else if ( ( obs_kind == QTY_POTENTIAL_TEMPERATURE )        .and. ( wrf%dom(id)%type_t >= 0 ) ) then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_t
else if ( ( obs_kind == QTY_SKIN_TEMPERATURE )              .and. ( wrf%dom(id)%type_tsk >= 0 ) )then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_tsk
else if ( ( obs_kind == QTY_GEOPOTENTIAL_HEIGHT )        .and. ( wrf%dom(id)%type_gz >= 0 ) )then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_gz
else if ( ( obs_kind == QTY_2M_TEMPERATURE )        .and. ( wrf%dom(id)%type_t2 >= 0 ) )then
   part_of_state_vector = .true.
   wrf_type = wrf%dom(id)%type_t2
else
   call error_handler(E_MSG, 'obs_kind_in_state_vector', &
      'obs_kind "'//trim(get_name_for_quantity(obs_kind))//'" is not in state vector 2', &
       source, revision, revdate)
   part_of_state_vector = .false.
   wrf_type = -1
endif

end subroutine obs_kind_in_state_vector

!--------------------------------------------------------------------
!> pass the vertical localization coordinate to assim_tools_mod
function query_vert_localization_coord()

integer :: query_vert_localization_coord

query_vert_localization_coord = vert_localization_coord

end function query_vert_localization_coord

!--------------------------------------------------------------------
!> read the time from the input file
function read_model_time(filename)

character(len=*),  intent(in) :: filename
type(time_type)               :: read_model_time

integer           :: year, month, day, hour, minute, second
integer           :: ret ! netcdf return code
integer           :: ndims, dimids(2), ivtype, ncid, var_id
character(len=80) :: varname
character(len=19) :: timestring
integer           :: i,  idims(2)

call nc_check( nf90_open(filename, NF90_NOWRITE, ncid), &
                  'opening', filename )

call nc_check( nf90_inq_varid(ncid, "Times", var_id), 'read_model_time', &
               'inq_varid Times' )
call nc_check( nf90_inquire_variable(ncid, var_id, varname, xtype=ivtype, &
               ndims=ndims, dimids=dimids), 'read_model_time', &
               'inquire_variable Times' )

do i=1,ndims ! isnt this just 1?
   call nc_check( nf90_inquire_dimension(ncid, dimids(i), &
                   len=idims(i)),'read_model_time','inquire_dimensions Times' )
enddo

call nc_check( nf90_get_var(ncid, var_id, timestring, &
               start = (/ 1, idims(2) /)), 'read_model_time','get_var Times' )

call get_wrf_date(timestring, year, month, day, hour, minute, second)
read_model_time = set_date(year, month, day, hour, minute, second)


call nc_check( nf90_close(ncid) , 'closing', filename)

end function read_model_time

!--------------------------------------------------------------------
!> write the time from the input file
subroutine write_model_time(ncid, dart_time)

use typeSizes
use netcdf

integer,         intent(in) :: ncid
type(time_type), intent(in) :: dart_time

integer :: dim_ids(2), var_id, ret
integer :: year, month, day, hour, minute, second
character(len=19) :: timestring

call get_date(dart_time, year, month, day, hour, minute, second)
call set_wrf_date(timestring, year, month, day, hour, minute, second)

call nc_begin_define_mode(ncid)

! Define Times variable if it does not exist
ret = nf90_inq_varid(ncid, "Times", var_id)
if (ret /= NF90_NOERR) then

   ! check to see if there is a time and date_str_length
   ret = nf90_inq_dimid(ncid, "Time", dim_ids(2))
   ! if Time dimension does not exist create it
   if (ret /= NF90_NOERR) then
      call nc_check(nf90_def_dim(ncid, "Time", nf90_unlimited, dim_ids(2)), &
        "write_model_time def_var dimension Time")
   endif

   ret = nf90_inq_dimid(ncid, "DateStrLen", dim_ids(1))
   if (ret /= NF90_NOERR) then
      ! if DateStrLen dimension does not exist create it.
      call nc_check(nf90_def_dim(ncid, "DateStrLen", len(timestring), dim_ids(1)), &
        "write_model_time def_var dimension dateStrLength")
   endif

   ! use id's to set Times(Time, DateStrLen)
   call nc_check(nf90_def_var(ncid, name="Times", xtype=nf90_char, &
      dimids=dim_ids, varid=var_id), "write_model_time def_var Times")
endif

call nc_end_define_mode(ncid)

call nc_check( nf90_put_var(ncid, var_id, timestring), &
               'write_model_time', 'put_var Times' )

end subroutine write_model_time

!--------------------------------------------------------------------

subroutine static_data_sizes(domain)
integer, intent(in) :: domain

print*
print*, '******** wrf_static_data_for_dart domain ',domain
print*, 'znu, dn, dnw, zs, znw ', size(wrf%dom(domain)%znu), &
                                  size(wrf%dom(domain)%dn ), &
                                  size(wrf%dom(domain)%dnw), &
                                  size(wrf%dom(domain)%zs ), &
                                  size(wrf%dom(domain)%znw)

print*, 'mub, hgt ', size(wrf%dom(domain)%mub), size(wrf%dom(domain)%hgt)

print*, 'latitude, latitude_u, latitude_v ', size(wrf%dom(domain)%latitude), &
                                             size(wrf%dom(domain)%latitude_u), &
                                             size(wrf%dom(domain)%latitude_v)

print*, 'longitude, longitude_u, longitude_v ', size(wrf%dom(domain)%longitude), &
                                                size(wrf%dom(domain)%longitude_u), &
                                                size(wrf%dom(domain)%longitude_v)

print*, 'phb             ', size(wrf%dom(domain)%phb)
print*, 'var_index       ', size(wrf%dom(domain)%var_index)
print*, 'var_size        ', size(wrf%dom(domain)%var_size)
print*, 'var_type        ', size(wrf%dom(domain)%var_type)
print*, 'var_index_list  ', size(wrf%dom(domain)%var_index_list)
print*, 'var_update_list ', size(wrf%dom(domain)%var_update_list)
print*, 'dart_kind       ', size(wrf%dom(domain)%dart_kind)
print*, 'land            ', size(wrf%dom(domain)%land)

print*, 'lower_bound,upper_bound ', size(wrf%dom(domain)%lower_bound), &
                                    size(wrf%dom(domain)%upper_bound)

print*, 'clamp_or_fail   ', size(wrf%dom(domain)%clamp_or_fail)

print*, 'description, units, stagger, coordinates ', size(wrf%dom(domain)%description), &
                                                     size(wrf%dom(domain)%units), &
                                                     size(wrf%dom(domain)%stagger), &
                                                     size(wrf%dom(domain)%coordinates)

print*

end subroutine static_data_sizes

!--------------------------------------------------------------------

end module model_mod

!> @}


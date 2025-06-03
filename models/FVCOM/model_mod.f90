! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! FVCOM interface to DART data assimilation system.
! This was built based on the existing FESOM and ROMS interfaces as templates
! and developed by adopting similar routines

! Units of the variables in state vector:
!
! u   (velocity):  meter / second
! h   (depth)   :  meter
! temperature   :  degrees C
! salinity      :  PSU          
!

! Routines in other modules that are used here.
use types_mod, only : r4, r8, i8, digits12, SECPERDAY, MISSING_R8,       &
                                 rad2deg, deg2rad, PI, MISSING_I, obstypelength

use time_manager_mod, only : time_type, set_time, set_date, get_date, &
                                 get_time, print_time, print_date,        &
                                 set_calendar_type, increment_time,       &
                                 operator(*),  operator(+), operator(-),  &
                                 operator(>),  operator(<), operator(/),  &
                                 operator(/=), operator(<=)

use     location_mod, only : location_type, get_close_type, &
                             loc_get_close_obs => get_close_obs,VERTISHEIGHT, &
                             loc_get_close_state => get_close_state, &
                             set_location, set_location_missing,get_location, &
                             write_location,get_dist


use    utilities_mod, only : error_handler, to_upper, &
                             E_ERR, E_MSG, register_module, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read,LOGFILEUNIT

use    obs_kind_mod, only : get_index_for_quantity

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, nc_check,        &
                                 nc_begin_define_mode, nc_end_define_mode,     &
                                 nc_open_file_readonly, nc_close_file,         &
                                 nc_define_dimension, nc_define_unlimited_dimension

use state_structure_mod, only : add_domain, get_domain_size,get_variable_name, &
                                get_varid_from_kind, state_structure_info, get_kind_index, &
                                get_model_variable_indices,get_dart_vector_index, &
                                get_index_start,get_variable_size
                                

use ensemble_manager_mod, only : ensemble_type

use distributed_state_mod, only : get_state

! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
use default_model_mod, only : pert_model_copies, read_model_time, write_model_time, &
                              init_time => fail_init_time, &
                              init_conditions => fail_init_conditions, &
                              convert_vertical_obs, convert_vertical_state, adv_1step

use netcdf
implicit none
private

! routines required by DART code - will be called from filter and other
! DART executables. 
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          end_model,              &
          static_init_model,      &
          nc_write_model_atts,    &
          get_close_obs,          &
          get_close_state,        &
          pert_model_copies,      &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          adv_1step,              &
          init_time,              &
          init_conditions,        &
          shortest_time_between_assimilations, &
          write_model_time


! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.
integer :: domid ! For state_structure_mod access

! module global storage; maintains values between calls, accessible by
! any subroutine
character(len=512) :: string1, string2, string3
integer :: dom_id ! used to access the state structure
integer :: Ncells,Nnodes,Nlays,Nmaxnodes

real(r8), allocatable, dimension(:) :: latN,lonN,latE,lonE,BDepth
real(r8), allocatable, dimension(:,:) :: sglE,sglN,sfht,WtDepths
integer, allocatable, dimension(:,:) :: NumNdNNs
integer, allocatable, dimension(:) :: NumNNs

type(location_type), allocatable :: cell_locations(:)
integer,             allocatable :: cell_kinds(:)
integer,             allocatable :: close_cell_inds(:)
type(location_type), allocatable :: node_locations(:)
integer,             allocatable :: node_kinds(:)
integer,             allocatable :: close_node_inds(:)
real(r8),            allocatable :: depths(:)

integer, parameter :: MAX_STATE_VARIABLES = 3
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 5
integer, parameter :: VARNAME_INDEX = 1
integer, parameter ::    KIND_INDEX = 2
integer, parameter ::  MINVAL_INDEX = 3
integer, parameter ::  MAXVAL_INDEX = 4
integer, parameter :: REPLACE_INDEX = 5

character(len=obstypelength) :: var_names(max_state_variables)
real(r8) :: var_ranges(max_state_variables,2)
logical  :: var_update(max_state_variables)
integer  :: var_qtys(  max_state_variables)

! variables which are in the module namelist
integer            :: vert_localization_coord = VERTISHEIGHT
integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 21600     ! 86400 | 43200 | 21600
real(r8)           :: model_perturbation_amplitude = 0.0001   ! tiny amounts
logical            :: diagnostic_metadata = .false.
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: model_analysis_filename = 'expno.year.oce.nc'
character(len=256) :: model_grid_filename = 'expno.nc'
character(len=NF90_MAX_NAME) :: variables(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = ' '

namelist /model_nml/             &
   model_analysis_filename,      &
   model_grid_filename,         &
   vert_localization_coord,      &
   diagnostic_metadata,          &
   assimilation_period_days,     &
   assimilation_period_seconds,  &
   model_perturbation_amplitude, &
   calendar,                     &
   variables,                    &
   debug

! Everything needed to describe a variable
integer :: nfields

integer         :: model_size          ! the state vector length
type(time_type) :: model_timestep      ! smallest time to adv model


contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.

subroutine static_init_model()

! Local variables - all the important ones have module scope

      integer :: ncid
      integer :: iunit, io
      integer :: ss, dd

      if ( module_initialized ) return ! only need to do this once.

      ! Print module information to log file and stdout.
      call register_module(source, revision, revdate)

      ! Early initialization before calling other routines
      module_initialized = .true.

      call find_namelist_in_file("input.nml", "model_nml", iunit)
      read(iunit, nml = model_nml, iostat = io)
      call check_namelist_read(iunit, io, "model_nml")

      ! Record the namelist values used for the run 
      if (do_nml_file()) write(nmlfileunit, nml=model_nml)
      if (do_nml_term()) write(     *     , nml=model_nml)

      ! This time is both the minimum time you can ask the model to advance
      ! (for models that can be advanced by filter) and it sets the assimilation
      ! window.  All observations within +/- 1/2 this interval from the current
      ! model time will be assimilated. If this is not settable at runtime 
      ! feel free to hardcode it and remove from the namelist.
      call set_calendar_type( calendar )  ! We can use Julian here

      model_timestep = set_model_time_step()

      call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

      write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
      call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)
      !-------------------------------------------------------------- 
      ! Read the model size, grid and variable dimensions
      ncid = nc_open_file_readonly(model_analysis_filename,'static_init_model')
      ! Get the FVCOM grid -- sizes and variables.
      call parse_variable_input(variables, ncid, model_analysis_filename, nfields)
      call nc_close_file(ncid, 'static_init_model', model_analysis_filename)
      call read_grid()
      ! Define which variables are in the model state
      domid = add_domain(model_analysis_filename, nfields, &
                  var_names  = var_names(1:nfields),       &
                  kind_list  = var_qtys(1:nfields),       &
                  clamp_vals = var_ranges(1:nfields,:),     &
                  update_list= var_update(1:nfields)        )
      
      model_size = get_domain_size(domid)

end subroutine static_init_model
!---------------------------------------------------------------------


!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 
function get_model_size()

      integer(i8) :: get_model_size

      if ( .not. module_initialized ) call static_init_model
      get_model_size = model_size

end function get_model_size
!-------------------------------------------------------------------



!------------------------------------------------------------------
! Given a state handle, a location, and a state quantity,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case a positive istatus should be returned.
!
! For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

subroutine model_interpolate(state_handle, ens_size, location, obs_type, expected_obs, istatus)

!       ERROR codes:
!
!       ISTATUS = 99:  general error in case something terrible goes wrong...
!       ISTATUS = 88:  this kind is not in the state vector
!       ISTATUS = 11:  Could not find a triangle that contains this lat/lon
!       ISTATUS = 12:  Height vertical coordinate out of model range.
!       ISTATUS = 13:  Missing value in interpolation.
!       ISTATUS = 16:  Don't know how to do vertical velocity for now
!       ISTATUS = 17:  Unable to compute pressure values
!       ISTATUS = 18:  altitude illegal
!       ISTATUS = 19:  could not compute u using RBF code
!       ISTATUS = 101: Internal error; reached end of subroutine without
!                      finding an applicable case.

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

! local storage

integer         :: ivar, obs_kind
real(r8)        :: llv(3)    ! lon/lat/vert
integer(i8)     :: index1  ! the DART index of the start of the variable of interest

! local storage

integer(i8)      :: surface_index
integer(i8)      :: ilayer, layer_below, layer_above
real(r8)         :: lon, lat, vert
real(r8)         :: depth_below, depth_above, layer_thick
integer(i8)      :: closest_index_above, closest_index_below
integer(i8)      :: index_above, index_below

type(location_type) :: location_above, location_below

if (.not. module_initialized ) call static_init_model

! This should be the result of the interpolation of a
! given kind (itype) of variable at the given location.
expected_obs(:) = MISSING_R8

! istatus for successful return should be 0. 
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
istatus(:) = 0

! rename for sanity - we can't change the argument names
! to this subroutine, but this really is a kind.
obs_kind = obs_type

! Make sure the DART state has the type (T,S,U,etc.) that we are asking for.
! If we cannot, simply return and 'fail' with an 88

ivar = get_varid_from_kind(domid, obs_kind)
if (ivar < 1) then
   istatus = 88
   return
endif

! Determine the offset into the DART state vector for this variable
index1 = get_index_start(domid,ivar)

! Decode the location into bits for error messages ...
llv  = get_location(location)
lon  = llv(1)    ! degrees East [0,360)
lat  = llv(2)    ! degrees North [-90,90]
vert = llv(3)    ! depth in meters ... even 2D fields have a value of 0.0

surface_index = find_closest_surface_location(location, obs_kind)

if (surface_index < 1) then ! nothing close
   istatus = 11
   return
endif

! If it is a surface variable, we're done.

if (get_variable_size(domid,ivar) == Nnodes) then
   index_above   = index1 + surface_index - 1
   expected_obs  = get_state(index_above, state_handle)
   istatus     = 0
   return
endif
!------------------------------------------------------------------------------ 
! For vertical direction pick the closest layer
layer_below = 0
LAYER: do ilayer = 1,Nlays
   if (WtDepths(ilayer,surface_index) > vert ) then
        layer_below = ilayer
        exit LAYER
   endif
enddo LAYER

if(layer_below == 0) then ! below the deepest level, RETURN
   istatus = 19
   return
elseif (layer_below == 1) then ! too shallow, RETURN
   istatus = 18
   return
endif

layer_above = layer_below - 1

depth_below = WtDepths(layer_below,surface_index) - vert
depth_above = vert - WtDepths(layer_above,surface_index)
layer_thick = WtDepths(layer_below,surface_index) - WtDepths(layer_above,surface_index)

closest_index_above = get_dart_vector_index(surface_index,layer_above,1,domid,ivar) !getting dart index
closest_index_below = get_dart_vector_index(surface_index,layer_below,1,domid,ivar)

!index_above = index1 + closest_index_above - 1  ! May not need this as we get the DART index directly 
!index_below = index1 + closest_index_below - 1

index_above = closest_index_above
index_below = closest_index_below

expected_obs  = ( &
                depth_below*get_state(index_above,state_handle) &
              + depth_above*get_state(index_below,state_handle)) &
              / layer_thick
istatus     = 0

if (do_output() .and. debug > 2) then
   call get_state_meta_data(index_above, location_above)
   call get_state_meta_data(index_below, location_below)

   call write_location(0,location_above,charstring=string1)
   call write_location(0,location      ,charstring=string2)
   call write_location(0,location_below,charstring=string3)

   write(logfileunit,*)
   write(     *     ,*)
   call error_handler(E_MSG,'model_interpolate', '... '//string1, &
              text2=string2, text3=string3)
endif

end subroutine model_interpolate
!----------------------------------------------------------------------------------------



!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = model_timestep

end function shortest_time_between_assimilations

!------------------------------------------------------------------
!>

function set_model_time_step()

! the static_init_model ensures that the model namelists are read.

type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call static_init_model

! these are from the namelist
!FIXME: sanity check these for valid ranges?
set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!-----------------------------------------------------------------------------------
! gets the length of a netCDF dimension given the dimension name.
! This bundles the nf90_inq_dimid and nf90_inquire_dimension routines
! into a slightly easier-to-use function.
function get_dimension_length(ncid, dimension_name, filename) result(dimlen)

integer                      :: dimlen
integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimension_name
character(len=*), intent(in) :: filename
integer :: DimID

      write(string1,*)'inq_dimid '//trim(dimension_name)//' '//trim(filename)
      write(string2,*)'inquire_dimension '//trim(dimension_name)//' '//trim(filename)

      call nc_check(nf90_inq_dimid(ncid, trim(dimension_name), DimID), &
                  'get_dimension_length',string1)
      call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
                  'get_dimension_length', string2)

end function get_dimension_length
!----------------------------------------------------------------------------------


!--------------------------------------------------------------------------------
! Reads the fvcom grid (nodes, elements, and dimensions)
subroutine read_grid()

integer :: ncid,VarID
logical, save :: grid_read = .false.

      if (grid_read) return ! grid reading needs to be done only once
      ! Read the static grid dimensions from the FVCOM grid file.
      call nc_check(nf90_open(trim(model_analysis_filename), nf90_nowrite, ncid), &
                  'get_grid_dimensions', 'open '//trim(model_analysis_filename))
      Ncells = get_dimension_length(ncid, 'nele',model_analysis_filename)
      Nnodes = get_dimension_length(ncid, 'node',model_analysis_filename)
      Nlays = get_dimension_length(ncid, 'siglay',model_analysis_filename)

      if (.not. allocated(latN)) allocate(latN(Nnodes))
      if (.not. allocated(lonN)) allocate(lonN(Nnodes))
      if (.not. allocated(latE)) allocate(latE(Ncells))
      if (.not. allocated(lonE)) allocate(lonE(Ncells))
      if (.not. allocated(sglN)) allocate(sglN(Nlays,Nnodes))
      if (.not. allocated(sglE)) allocate(sglE(Nlays,Ncells))
      if (.not. allocated(BDepth)) allocate(BDepth(Nnodes))
      if (.not. allocated(sfht)) allocate(sfht(1,Nnodes))
      if (.not. allocated(WtDepths)) allocate(WtDepths(Nlays,Nnodes))
      
      call nc_check(nf90_inq_varid(ncid, 'lon', VarID), &
      'get_grid', 'inq_varid lon_node '//trim(model_analysis_filename))
      call nc_check(nf90_get_var( ncid, VarID, lonN), &
            'get_grid', 'get_var lon_node '//trim(model_analysis_filename))

      where (lonN < 0.0_r8) lonN = lonN + 360.0_r8

      call nc_check(nf90_inq_varid(ncid, 'lat', VarID), &
            'get_grid', 'inq_varid lat_node '//trim(model_analysis_filename))
      call nc_check(nf90_get_var( ncid, VarID, latN), &
            'get_grid', 'get_var lat_node '//trim(model_analysis_filename))

      call nc_check(nf90_inq_varid(ncid, 'lonc', VarID), &
      'get_grid', 'inq_varid lon_cell '//trim(model_analysis_filename))
      call nc_check(nf90_get_var( ncid, VarID, lonE), &
            'get_grid', 'get_var lon_cell '//trim(model_analysis_filename))

      where (lonE < 0.0_r8) lonE = lonE + 360.0_r8

      call nc_check(nf90_inq_varid(ncid, 'latc', VarID), &
            'get_grid', 'inq_varid lat_cell '//trim(model_analysis_filename))
      call nc_check(nf90_get_var( ncid, VarID, latE), &
            'get_grid', 'get_var lat_cell '//trim(model_analysis_filename))

      call nc_check(nf90_inq_varid(ncid, 'siglay', VarID), &
            'get_grid', 'inq_varid siglay '//trim(model_analysis_filename))
      call nc_check(nf90_get_var( ncid, VarID, sglN), &
            'get_grid', 'get_var siglay '//trim(model_analysis_filename))

      call nc_check(nf90_inq_varid(ncid, 'siglay_center', VarID), &
            'get_grid', 'inq_varid siglay_center '//trim(model_analysis_filename))
      call nc_check(nf90_get_var( ncid, VarID, sglE), &
            'get_grid', 'get_var siglay_center '//trim(model_analysis_filename))

      call nc_check(nf90_inq_varid(ncid, 'h', VarID), &
            'get_grid', 'inq_varid bathymetry '//trim(model_analysis_filename))
      call nc_check(nf90_get_var( ncid, VarID, BDepth), &
            'get_grid', 'get_var bathymetry '//trim(model_analysis_filename))
      
      call nc_check(nf90_close(ncid), &
                  'get_var','close '//trim(model_analysis_filename))

      call calc_water_depth()

      !----------------------------- other grid properties are read from a solution file -----------------
      call nc_check(nf90_open(trim(model_grid_filename), nf90_nowrite, ncid), &
                  'get_grid_dimensions', 'open '//trim(model_grid_filename))
      Nmaxnodes = get_dimension_length(ncid, 'maxnode',model_grid_filename)
      if (.not. allocated(NumNdNNs)) allocate(NumNdNNs(Nmaxnodes,Nnodes))
      if (.not. allocated(NumNNs)) allocate(NumNNs(Nnodes))

       call nc_check(nf90_inq_varid(ncid, 'ntsn', VarID), &
            'get_grid', 'inq_varid surround_node number '//trim(model_grid_filename))
       call nc_check(nf90_get_var( ncid, VarID, NumNNs), &
            'get_grid', 'get_var surround_node number '//trim(model_grid_filename))

       call nc_check(nf90_inq_varid(ncid, 'nbsn', VarID), &
            'get_grid', 'inq_varid surround_node number '//trim(model_grid_filename))
       call nc_check(nf90_get_var( ncid, VarID, NumNdNNs), &
            'get_grid', 'get_var surround_node number '//trim(model_grid_filename))

      call nc_check(nf90_close(ncid), &
                  'get_var','close '//trim(model_grid_filename))

      if (debug > 5) &
      call error_handler(E_MSG,'read_grid','get close lookup table initialized')

      if (debug > 1) then
            write(string1,*)'nCells      is ', Ncells
            write(string2,*)'nNodes      is ', Nnodes
            write(string3,*)'nVertLevels is ', Nlays
            call error_handler(E_MSG,'read_grid','... '//string1,text2=string2,text3=string3)

            write(string1,*)'latitude  range ',minval(latN(:)),maxval(latN(:))
            write(string2,*)'longitude range ',minval(lonN(:)),maxval(lonN(:))
            write(string3,*)'depths    range ',minval(WtDepths),maxval(WtDepths)
            call error_handler(E_MSG,'read_grid','... '//string1,text2=string2,text3=string3)
      endif

      grid_read = .true.

end subroutine
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

integer  :: nf, myindx
real(r8) :: lon, lat, depth
integer  :: var_id, dom_id, variable_type
integer  ::  iloc, jloc,kloc

if (.not. module_initialized) call static_init_model
iloc = -1

! from the dart index get the local variables indices - can be nodes or cells
! currently DA is done only for T and S which are based on nodes (iloc indices)

!call get_model_variable_indices(index_in, iloc, jloc, kloc, &
!            var_id=var_id, dom_id=dom_id, kind_index=variable_type)

call get_model_variable_indices(index_in, iloc, jloc, kloc, &
            var_id, dom_id, variable_type)


if(iloc == -1) then
     write(string1,*) 'Problem, cannot find base_offst, index_in is: ', index_in
     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
endif

lon   = lonN(iloc)
lat   = latN(iloc)
depth = WtDepths(jloc,iloc)
location = set_location(lon,lat,depth,VERTISHEIGHT)


! should be set to the actual location using set_location()
location = set_location_missing()

! should be set to the physical quantity, e.g. QTY_TEMPERATURE
if (present(var_type)) then
   var_type = variable_type

   if(var_type == MISSING_I) then
      write(string1,*) 'Cannot find DART QTY for indx ', index_in
      write(string2,*) 'variable "'//trim(get_variable_name(dom_id, var_id))//'"'
      call error_handler(E_ERR, 'get_state_meta_data', string1, &
                 source, revision, revdate, text2=string2)
   endif
endif

end subroutine get_state_meta_data
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine calc_water_depth()
integer   :: i,j,ncid,VarID

      call nc_check(nf90_open(trim(model_analysis_filename), nf90_nowrite, ncid), &
                  'get_grid', 'open '//trim(model_analysis_filename))
      call nc_check(nf90_inq_varid(ncid, 'zeta', VarID), &
                  'get_grid', 'inq_varid zeta '//trim(model_analysis_filename))
      call nc_check(nf90_get_var(ncid, VarID, sfht), &
                  'get_grid', 'get_var zeta '//trim(model_analysis_filename))
      call nc_check(nf90_close(ncid), &
                  'get_var','close '//trim(model_analysis_filename))
      
      do i=1,Nnodes
            WtDepths(:,i) = BDepth(i)*sglN(:,i)
      end do

end subroutine calc_water_depth
!-----------------------------------------------------------------------


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc            ! handle to a get_close structure
integer,                       intent(in)    :: base_type     ! observation TYPE
type(location_type),           intent(inout) :: base_loc      ! location of interest
type(location_type),           intent(inout) :: locs(:)       ! obs locations
integer,                       intent(in)    :: loc_qtys(:)   ! QTYS for obs
integer,                       intent(in)    :: loc_types(:)  ! TYPES for obs
integer,                       intent(out)   :: num_close     ! how many are close
integer,                       intent(out)   :: close_ind(:)  ! incidies into the locs array
real(r8),            optional, intent(out)   :: dist(:)       ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_obs'

call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                          num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs
!------------------------------------------------------------------


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc           ! handle to a get_close structure
type(location_type),           intent(inout) :: base_loc     ! location of interest
integer,                       intent(in)    :: base_type    ! observation TYPE
type(location_type),           intent(inout) :: locs(:)      ! state locations
integer,                       intent(in)    :: loc_qtys(:)  ! QTYs for state
integer(i8),                   intent(in)    :: loc_indx(:)  ! indices into DART state vector
integer,                       intent(out)   :: num_close    ! how many are close
integer,                       intent(out)   :: close_ind(:) ! indices into the locs array
real(r8),            optional, intent(out)   :: dist(:)      ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_state'


call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)


end subroutine get_close_state
!-----------------------------------------------------------------


!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()


end subroutine end_model


!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

if ( .not. module_initialized ) call static_init_model

! put file into define mode.

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model", "template")

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts


!----------------------------------------------------------------------
! Read model variables and attributes from model analysis file
subroutine parse_variable_input( state_variables, ncid, filename, ngood )

character(len=*), intent(in)  :: state_variables(:,:)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: filename
integer,          intent(out) :: ngood

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: dimname
integer :: i, j, VarID, dimlen, numdims
logical :: failure

character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 5

real(r8) :: minvalue, maxvalue
integer  :: ios

if ( .not. module_initialized ) call static_init_model

var_names  = 'no_variable_specified'
var_qtys   = MISSING_I
var_ranges = MISSING_R8
var_update = .false.

failure = .FALSE. ! perhaps all with go well

ngood = 0
MyLoop : do i = 1, MAX_STATE_VARIABLES

   if ( variables(1,i) == ' ' .and. variables(2,i) == ' ' ) exit MyLoop ! Found end of list.

   if ( any(state_variables(:,i) == ' ') ) then
      string1 = '...  model_nml:"variables" not fully specified'
      write(string2,*)'failing on line ',i
      call error_handler(E_ERR, 'parse_variable_input', string1, &
                 source, revision, revdate, text2=string2)
   endif

   varname      = trim(state_variables(VARNAME_INDEX,i))
   dartstr      = trim(state_variables(   KIND_INDEX,i))
   minvalstring = trim(state_variables( MINVAL_INDEX,i))
   maxvalstring = trim(state_variables( MAXVAL_INDEX,i))
   state_or_aux = trim(state_variables(REPLACE_INDEX,i))
   call to_upper(state_or_aux)

   var_names(i) = trim(varname)
   var_qtys(i)  = get_index_for_quantity(dartstr)

   read(minvalstring,*,iostat=ios) minvalue
   if (ios == 0) var_ranges(i,1) = minvalue

   read(maxvalstring,*,iostat=ios) maxvalue
   if (ios == 0) var_ranges(i,1) = maxvalue

   if (state_or_aux == 'UPDATE' ) var_update(i) = .true.

   ! Make sure DART kind is valid

   if( var_qtys(i) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') &
            trim(dartstr)
      call error_handler(E_ERR,'parse_variable_input',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in model analysis variable list

   write(string1,'(''variable '',a,'' in '',a)') trim(varname), trim(filename)
   write(string2,'(''there is no '',a)') trim(string1)
   call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
                 'parse_variable_input', trim(string2))

   ! Make sure variable is defined by (Time,nCells) or (Time,nCells,vertical)
   ! unable to support Edges or Vertices at this time.

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
                 'parse_variable_input', 'inquire '//trim(string1))

   ! DimensionLoop : do j = 1,numdims

   !    write(string2,'(''inquire dimension'',i2,'' of '',a)') j,trim(string1)
   !    call nc_check(nf90_inquire_dimension(ncid, dimIDs(j), len=dimlen, name=dimname), &
   !                                        'parse_variable_input', trim(string2))
   !    select case ( trim(dimname) )
   !       case ('T')
   !          ! supported - do nothing
   !       case ('nodes_2d')
   !          ! supported - do nothing
   !       case ('nodes_3d')
   !          ! supported - do nothing
   !       case default
   !          write(string2,'(''unsupported dimension '',a,'' in '',a)') trim(dimname),trim(string1)
   !          call error_handler(E_MSG,'parse_variable_input',string2,source,revision,revdate)
   !          failure = .TRUE.
   !    end select

   ! enddo DimensionLoop

   ! if (failure) then
   !     string2 = 'unsupported dimension(s) are fatal'
   !     call error_handler(E_ERR,'parse_variable_input',string2,source,revision,revdate)
   ! endif

   ! Record the contents of the DART state vector

   if (debug > 0) then
      write(string1,*)'variable ',i,' is ',trim(varname), ' ', trim(dartstr)
      call error_handler(E_MSG,'parse_variable_input',string1)
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == MAX_STATE_VARIABLES) then
   string1 = 'WARNING: There is a possibility you need to increase ''MAX_STATE_VARIABLES'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'parse_variable_input',string1,text2=string2)
endif

end subroutine parse_variable_input
!----------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------
function find_closest_surface_location(location, obs_kind)
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_kind
integer                         :: find_closest_surface_location

type(location_type) :: LocTempNode,LocStNode
integer :: startNd, TotNNs, inn, NxtStrNd, cnt
real(r8) :: lon,lat,vert,DisStoLoc,DisNNtoLoc,MinDis,startNdFL 

find_closest_surface_location = -1

startNdFL = Nnodes/2.0
startNd = floor(startNdFL)
vert = 2.0     ! dummy value to create a location type variable
cnt = 0
outer: do while(cnt < Nnodes)
      lon = lonN(startNd)
      lat = latN(startNd)
      LocStNode = set_location(lon,lat,vert,VERTISHEIGHT)
      DisStoLoc = get_dist(location, LocStNode, no_vert=.true.)
      TotNNs = NumNNs(startNd)
      NxtStrNd = -1
      MinDis = DisStoLoc
      do inn=1,TotNNs      ! No. of NNs
            lon = lonN(NumNdNNs(inn,startNd))  ! lon of NN
            lat = latN(NumNdNNs(inn,startNd))  ! lat of NN
            LocTempNode = set_location(lon,lat,vert,VERTISHEIGHT) ! location of NN
            DisNNtoLoc = get_dist(location, LocTempNode, no_vert=.true.)   ! distance to NN
            if(DisNNtoLoc >= DisStoLoc) then
                  CYCLE
            elseif((DisNNtoLoc < DisStoLoc) .AND. (DisNNtoLoc < MinDis)) then
                  NxtStrNd = NumNdNNs(inn,startNd)
                  MinDis = DisNNtoLoc
            end if
      end do
      cnt = cnt + 1
      if(NxtStrNd == -1) then
            EXIT outer
      else
            startNd = NxtStrNd
      end if
end do outer

find_closest_surface_location = startNd

end function find_closest_surface_location
!----------------------------------------------------------------------------------------

end module model_mod
!========================================================================================

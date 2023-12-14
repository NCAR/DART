! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. Do not change the arguments
! for the public routines.

use        types_mod, only : r8, i8, MISSING_R8, vtablenamelength

use time_manager_mod, only : time_type, set_time

use     location_mod, only : location_type, get_close_type, &
                             loc_get_close_obs => get_close_obs, &
                             loc_get_close_state => get_close_state, &
                             set_location, set_location_missing, &
                             get_location, VERTISLEVEL

use    utilities_mod, only : register_module, error_handler, &
                             E_ERR, E_MSG, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read, &
                             to_upper

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode, &
                                 NF90_MAX_NAME, nc_open_file_readonly, &
                                 nc_get_variable, nc_get_variable_size, nc_close_file

use state_structure_mod, only : add_domain, get_domain_size, get_model_variable_indices, &
                                get_varid_from_kind, get_dart_vector_index

use obs_kind_mod, only : get_index_for_quantity, QTY_U_CURRENT_COMPONENT, &
                         QTY_V_CURRENT_COMPONENT, QTY_DRY_LAND, &
                         QTY_LAYER_THICKNESS

use distributed_state_mod, only: get_state

use ensemble_manager_mod, only : ensemble_type

! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
use default_model_mod, only : pert_model_copies, write_model_time, &
                              init_time => fail_init_time, &
                              init_conditions => fail_init_conditions, &
                              convert_vertical_obs, convert_vertical_state, adv_1step

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


character(len=256), parameter :: source = "model_mod.f90"
character(len=256), parameter :: ocean_geometry = "/glade/work/rarmstrong/cesm/cesm2_3_alpha12b+mom6_marbl/components/mom/standalone/examples/single_column_MARBL/BATS_state_estimation/ensemble/baseline/ocean_geometry.nc"
logical :: module_initialized = .false.
integer :: dom_id  ! used to access the state structure
integer :: nfields ! number of fields in the state vector
integer :: nz      ! the number of vertical layers
integer :: model_size
type(time_type) :: assimilation_time_step
real(r8), parameter :: geolon = 360 - 64.0
real(r8), parameter :: geolat = 31.0
real(r8) :: basin_depth(2,2)

logical, parameter :: debug_interpolation = .false.

! parameters to be used in specifying the DART internal state
integer, parameter :: modelvar_table_height = 13
integer, parameter :: modelvar_table_width = 5

! defining the variables that will be read from the namelist
character(len=256) :: template_file
integer            :: time_step_days
integer            :: time_step_seconds
character(len=vtablenamelength) &
   :: model_state_variables(modelvar_table_height * modelvar_table_width)

namelist /model_nml/ template_file, &
                     time_step_days, &
                     time_step_seconds, &
                     model_state_variables

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.

subroutine static_init_model()

integer :: iunit, io
character(len=vtablenamelength) :: variable_table(modelvar_table_height, modelvar_table_width) 
integer :: state_qty_list(modelvar_table_height)
real(r8):: state_clamp_vals(modelvar_table_height, 2)
logical :: update_var_list(modelvar_table_height)

module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source)

! Read values from the namelist

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
assimilation_time_step = set_time(time_step_seconds, &
                                  time_step_days)

call verify_state_variables(model_state_variables, nfields, variable_table, &
                            state_qty_list, state_clamp_vals, update_var_list)

! print *, "variable_table   = ",variable_table
! print *, "state_qty_list   = ",state_qty_list
! print *, "state_clamp_vals = ",state_clamp_vals
! print *, "update_var_list  = ",update_var_list

dom_id = add_domain(template_file, nfields, &
                    var_names = variable_table(1:nfields, 1), &
                    kind_list = state_qty_list(1:nfields), &
                    clamp_vals = state_clamp_vals, &
                    update_list = update_var_list(1:nfields))

model_size = get_domain_size(dom_id)

call read_num_layers     ! setting the value of nz
call read_ocean_geometry ! determining the basin depth

end subroutine static_init_model

!------------------------------------------------------------------
! Reads the simulation length from a netCDF file.

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type) :: read_model_time

integer :: ncid
character(len=*), parameter :: routine = 'read_model_time'
real(r8) :: days
type(time_type) :: mom6_time
integer :: mom_base_date_in_days, mom_days

mom_base_date_in_days = 139157 ! 1982 1 1 0 0
ncid = nc_open_file_readonly(filename, routine)

call nc_get_variable(ncid, 'Time', days, routine)

call nc_close_file(ncid, routine)
    
mom_days = int(days) + mom_base_date_in_days

read_model_time = set_time(0,mom_days)

end function read_model_time

!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = get_domain_size(dom_id)

end function get_model_size

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

subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: qty
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

integer     :: qty_id, thickness_id, layer_index
integer(i8) :: qty_index, thickness_index, ens_index
real(8)     :: requested_depth, layerdepth_bottom, layerdepth_center, &
               layerdepth_top, depth_below, depth_above, val_below, val_above
real(8)     :: layer_thicknesses(nz, ens_size)
real(8)     :: state_slice(ens_size)
real(8)     :: loc_temp(3)

if ( .not. module_initialized ) call static_init_model

qty_id = get_varid_from_kind(dom_id, qty)

! extracting the requested depth value from `location`.
loc_temp = get_location(location)
requested_depth = loc_temp(3)

if(debug_interpolation) then; print *, "REQUESTED DEPTH = ",requested_depth; end if

! extracting the current layer thicknesses for each ensemble member.
thickness_id = get_varid_from_kind(dom_id, QTY_LAYER_THICKNESS)

do layer_index = 1, nz
    ! layer thicknesses are always extracted from the grid cell at index (1, 1), since the four grid
    ! cells are supposed to represent identical columns.
    thickness_index = get_dart_vector_index(1, 1, layer_index, dom_id, thickness_id)
    layer_thicknesses(layer_index, :) = get_state(thickness_index, state_handle)
end do

! performing the interpolation for each ensemble member individually.
do ens_index = 1, ens_size
    if(debug_interpolation) then
        print *, "BASIN DEPTH     = ",-basin_depth(1,1)
        print *, "----------------------------------------------------"
        print *, "computing centers of layers"
        print *, "----------------------------------------------------"
    end if

    ! locating the layer index to be used as the upper interpolation point for this ensemble member.
    layer_index = nz
    layerdepth_bottom = -basin_depth(1,1)       ! depth at bottom of the layer given by layer_index.
    layerdepth_center = layerdepth_bottom &     ! depth at center of the layer given by layer_index.
                          + .5*layer_thicknesses(layer_index, ens_index)
    layerdepth_top = layerdepth_bottom &        ! depth at top of the layer given by layer_index.
                          + layer_thicknesses(layer_index, ens_index)

    if(debug_interpolation) then
        print *, "layer = ",layer_index,", center = ",layerdepth_center,", bottom = ",layerdepth_bottom,", top = ",layerdepth_top
    end if
    
    do while((layerdepth_center < requested_depth) .and. (layer_index > 1))
        layer_index = layer_index - 1
        layerdepth_bottom = layerdepth_bottom + layer_thicknesses(layer_index, ens_index)
        layerdepth_center = layerdepth_bottom + 0.5 * layer_thicknesses(layer_index, ens_index)
        layerdepth_top = layerdepth_bottom + layer_thicknesses(layer_index, ens_index)

        if(debug_interpolation) then
            print *, "layer = ",layer_index,", bottom = ",layerdepth_bottom,", center = ",layerdepth_center,", top = ",layerdepth_top
        end if

    end do

    ! having located the index of the upper interpolation layer, we now calculate the interpolation.
    if((requested_depth < -basin_depth(1,1)) .or. (layerdepth_top < requested_depth)) then
        ! case where the requested depth is below the ocean floor, or above the ocean surface

        if(debug_interpolation) then
            print *, "----------------------------------------------------"
            print *, "depth is below ocean floor or above ocean surface"
            print *, "----------------------------------------------------"
        end if

        istatus(ens_index) = 1
        expected_obs(ens_index) = MISSING_R8

    else if((layer_index == nz) .or. (layerdepth_center < requested_depth)) then
        ! case where the requested depth is either in the bottom half of the deepest layer,
        ! or the top half of the shallowest layer. In both cases, the "interpolated" value is
        ! simply the current value of that layer in MOM6.

        if(debug_interpolation) then
            print *, "----------------------------------------------------"
            print *, "only using value from layer ",layer_index
            print *, "----------------------------------------------------"
        end if

        istatus(ens_index) = 0
        qty_index = get_dart_vector_index(1, 1, layer_index, dom_id, qty_id)
        state_slice = get_state(qty_index, state_handle)
        expected_obs(ens_index) = state_slice(ens_index)

        if(debug_interpolation) then; print *, "final value = ",expected_obs(ens_index); end if

    else
        ! case where the requested depth is above the center of some layer, and below
        ! the center of another. We interpolate linearly between the nearest layers above
        ! and below.

        if(debug_interpolation) then
            print *, "----------------------------------------------------"
            print *, "interpolating between layers ",layer_index," and ",(layer_index + 1)
        end if

        istatus(ens_index) = 0

        ! computing the depths at the centers of the nearest layers above and below
        depth_above = layerdepth_center
        depth_below = layerdepth_top - layer_thicknesses(layer_index, ens_index) &
                                        - .5*layer_thicknesses(layer_index + 1, ens_index)
        
        ! extracting the quantity values at the layers above and below
        qty_index = get_dart_vector_index(1, 1, layer_index, dom_id, qty_id)
        state_slice = get_state(qty_index, state_handle)
        val_above = state_slice(ens_index)

        qty_index = get_dart_vector_index(1, 1, layer_index + 1, dom_id, qty_id)
        state_slice = get_state(qty_index, state_handle)
        val_below = state_slice(ens_index)

        if(debug_interpolation) then
            print *, "corresponding to the values: ",val_above," and ",val_below
            print *, "----------------------------------------------------"
        end if

        ! linear interpolation
        expected_obs(ens_index) = val_above + (requested_depth - depth_above) * (val_below - val_above) / (depth_below - depth_above)

        if(debug_interpolation) then; print *, "final value = ",expected_obs(ens_index); end if
    end if
end do

end subroutine model_interpolate



!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = assimilation_time_step

end function shortest_time_between_assimilations



!------------------------------------------------------------------
! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.

subroutine get_state_meta_data(index_in, location, qty)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: qty

real(r8) :: lat, lon
integer :: lon_index, lat_index, level, local_qty

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, lon_index, lat_index, level, kind_index=local_qty)

location = set_location(geolon, geolat, real(level,r8), VERTISLEVEL)

if (present(qty)) then
    qty = local_qty
endif


end subroutine get_state_meta_data


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

!------------------------------------------------------------------
! Verify that the namelist was filled in correctly, and check
! that there are valid entries for the dart_kind.
! Returns a table with columns:
!
! netcdf_variable_name ; dart_qty_string ; update_string

subroutine verify_state_variables(state_variables, ngood, table, qty_list, clamp_vals, update_var)

character(len=*),  intent(inout) :: state_variables(:)
integer,           intent(out) :: ngood
character(len=*),  intent(out) :: table(:,:)
integer,           intent(out) :: qty_list(:)   ! kind number
real(r8),          intent(out) :: clamp_vals(:,:)
logical,           intent(out) :: update_var(:) ! logical update

integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname, dartstr, lowerbound, upperbound, update
character(len=256) :: string1, string2

if ( .not. module_initialized ) call static_init_model

nrows = size(table,1)

ngood = 0

MyLoop : do i = 1, nrows

    varname =    trim(state_variables(5*i -4))
    dartstr =    trim(state_variables(5*i -3))
    lowerbound = trim(state_variables(5*i -2))
    upperbound = trim(state_variables(5*i -1))
    update  =    trim(state_variables(5*i   ))

    call to_upper(update)

    table(i,1) = trim(varname)
    table(i,2) = trim(dartstr)
    table(i,3) = trim(lowerbound)
    table(i,4) = trim(upperbound)
    table(i,5) = trim(update)

    if ( table(i,1) == ' ' .and. table(i,2) == ' ' .and. table(i,3) == ' ' .and. table(i,4) == ' ' .and. table(i,5) == ' ') exit MyLoop ! Found end of list.

    if ( table(i,1) == ' ' .or. table(i,2) == ' ' .or. table(i,3) == ' ' .or. table(i,4) == ' ' .or. table(i,5) == ' ') then
        string1 = 'model_nml:model_state_variables not fully specified'
        call error_handler(E_ERR,'verify_state_variables',string1)
    endif

    ! Make sure DART qty is valid

    qty_list(i) = get_index_for_quantity(dartstr)
    if( qty_list(i)  < 0 ) then
        write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
        call error_handler(E_ERR,'verify_state_variables',string1)
    endif
    
    ! Make sure the update variable has a valid name

    select case (update)
        case ('UPDATE')
            update_var(i) = .true.
        case ('NO_COPY_BACK')
            update_var(i) = .false.
        case default
            write(string1,'(A)')  'only UPDATE or NO_COPY_BACK supported in model_state_variable namelist'
            write(string2,'(6A)') 'you provided : ', trim(varname), ', ', trim(dartstr), ', ', trim(update)
            call error_handler(E_ERR,'verify_state_variables',string1, text2=string2)
    end select

    ! reading the clamp values

    if (table(i, 3) /= 'NA') then
        read(table(i,3), '(d16.8)') clamp_vals(i,1)
    else
        clamp_vals(i,1) = MISSING_R8
    endif
    
    if (table(i,4) /= 'NA') then
        read(table(i,4), '(d16.8)') clamp_vals(i,2)
    else
        clamp_vals(i,2) = MISSING_R8
    endif

    ngood = ngood + 1
enddo MyLoop


end subroutine verify_state_variables

!------------------------------------------------------------
function on_v_grid(qty)

integer, intent(in)  :: qty
logical :: on_v_grid

if (qty == QTY_V_CURRENT_COMPONENT) then
    on_v_grid = .true.
else
    on_v_grid = .false.
endif

end function on_v_grid

!----------------------------------------------------------
function on_u_grid(qty)

integer, intent(in)  :: qty
logical :: on_u_grid

if (qty == QTY_U_CURRENT_COMPONENT) then
    on_u_grid = .true.
else
    on_u_grid = .false.
endif

end function on_u_grid

!------------------------------------------------------------
! Read number of vertical layers from mom6 template file
subroutine read_num_layers()

integer :: ncid

character(len=*), parameter :: routine = 'read_num_layers'

ncid = nc_open_file_readonly(template_file)

call nc_get_variable_size(ncid, 'Layer', nz)

call nc_close_file(ncid)

end subroutine read_num_layers

!------------------------------------------------------------
! ocean_geom are 2D state sized static data
! HK Do these arrays become too big in high res cases?
subroutine read_ocean_geometry()

integer :: ncid

character(len=*), parameter :: routine = 'read_ocean_geometry'

! Need nx, ny
if ( .not. module_initialized ) call static_init_model

ncid = nc_open_file_readonly(ocean_geometry)

call nc_get_variable(ncid, 'D', basin_depth, routine)

call nc_close_file(ncid)

end subroutine read_ocean_geometry

!===================================================================
! End of model_mod
!===================================================================
end module model_mod


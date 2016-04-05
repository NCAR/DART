! DART software - Copyright 2004 - 2015 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> This is the interface between the GCOM ocean model and DART.

module model_mod

! Modules that are absolutely required for use are listed

use        types_mod, only : r4, r8, SECPERDAY, MISSING_R4, MISSING_R8, MISSING_I, &
                             rad2deg, PI, obstypelength, metadatalength, earth_radius

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+),  operator(-),          &
                             operator(>),  operator(<),  operator(/),          &
                             operator(/=), operator(<=), operator(==)

use     location_mod, only : location_type, get_dist, get_close_maxdist_init,  &
                             get_close_obs_init, set_location,                 &
                             VERTISHEIGHT, get_location, vert_is_height,       &
                             vert_is_level, vert_is_surface, get_close_type,   &
                             loc_get_close_obs => get_close_obs,               &
                             get_close_obs_destroy

use    utilities_mod, only : register_module, logfileunit, get_unit,     &
                             error_handler, E_ERR, E_WARN, E_MSG,        &
                             nc_check, do_output, to_upper,              &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, file_exist, find_textfile_dims,  &
                             file_to_text, do_output, close_file

use     obs_kind_mod, only : KIND_TEMPERATURE,           &
                             KIND_SALINITY,              &
                             KIND_DRY_LAND,              &
                             KIND_U_CURRENT_COMPONENT,   &
                             KIND_V_CURRENT_COMPONENT,   &
                             KIND_SEA_SURFACE_HEIGHT,    &
                             KIND_SEA_SURFACE_PRESSURE,  &
                             KIND_POTENTIAL_TEMPERATURE, &
                             get_raw_obs_kind_index,     &
                             get_raw_obs_kind_name,      &
                             get_obs_kind_var_type

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use         sort_mod, only: index_sort

use typesizes
use netcdf

implicit none
private

! These routines are required by DART.
! Interfaces to these routines are fixed and cannot be changed in any way.

public :: static_init_model,      &
          get_model_size,         &
          get_model_time_step,    &
          get_state_meta_data,    &
          model_interpolate,      &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          init_conditions,        &
          init_time,              &
          pert_model_state,       &
          adv_1step,              &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model,     &
          end_model

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.
public :: get_gcom_restart_filename, &
          gcom_file_to_dart_vector,  &
          dart_vector_to_gcom_file,  &
          DART_get_var,              &
          test_interpolation,        &
          write_gcom_timeinfo


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2, string3

logical, save :: module_initialized = .false.
logical, save :: interpolation_initialized = .false.

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

! Codes for restricting the range of a variable
integer, parameter :: BOUNDED_NONE  = 0 ! ... unlimited range
integer, parameter :: BOUNDED_BELOW = 1 ! ... minimum, but no maximum
integer, parameter :: BOUNDED_ABOVE = 2 ! ... maximum, but no minimum
integer, parameter :: BOUNDED_BOTH  = 3 ! ... minimum and maximum

! Variables pertaining to the makeup of the DART vector
integer :: nfields
integer, parameter :: max_state_variables = 40
integer, parameter :: num_state_table_columns = 5
character(len=obstypelength) :: variable_table(max_state_variables, num_state_table_columns)

! Codes for interpreting the columns of the variable_table
integer, parameter :: VT_VARNAMEINDX = 1 ! ... variable name
integer, parameter :: VT_KINDINDX    = 2 ! ... DART kind
integer, parameter :: VT_MINVALINDX  = 3 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX  = 4 ! ... maximum value if any
integer, parameter :: VT_STATEINDX   = 5 ! ... update (state) or not

! things which can/should be in the model_nml
logical  :: output_state_vector = .true.
integer  :: assimilation_period_days = 1
integer  :: assimilation_period_seconds = 0
real(r8) :: model_perturbation_amplitude = 0.2
logical  :: update_dry_cell_walls = .false.
integer  :: debug = 0   ! turn up for more and more debug messages
character(len=256) :: gcom_restart_file  = 'gcom_restart.nc'
character(len=256) :: gcom_geometry_file = 'gcom_geometry.nc'
character(len=obstypelength) :: gcom_variables(max_state_variables*num_state_table_columns) = ' '

! FIXME: currently the update_dry_cell_walls namelist value DOES
! NOTHING.  it needs additional code to detect the cells which are
! wet, but within 1 cell of the bottom/sides/etc.

namelist /model_nml/             &
   gcom_restart_file,            &
   gcom_geometry_file,           &
   output_state_vector,          &
   assimilation_period_days,     &
   assimilation_period_seconds,  &
   model_perturbation_amplitude, &
   update_dry_cell_walls,        &
   debug,                        &
   gcom_variables

integer, parameter :: S_index     = 1  ! TJH FIXME remove at some point
integer, parameter :: T_index     = 2  ! TJH FIXME remove at some point
integer, parameter :: U_index     = 3  ! TJH FIXME remove at some point
integer, parameter :: V_index     = 4  ! TJH FIXME remove at some point
integer, parameter :: PSURF_index = 5  ! TJH FIXME remove at some point

integer :: start_index(5) ! TJH FIXME this gets replaced

! Grid parameters - the values will be read from a
! standard GCOM namelist and filled in here.

! nxp1, nyp1 and nzp1 are the maximum size of the grid
integer :: nxp1=-1, nyp1=-1, nzp1=-1  ! number of grid edges
integer :: nx=-1, ny=-1, nz=-1  ! number of grid centers

! locations of cell centers (C) and edges (G) for each axis.
real(r8), allocatable :: ZC(:), ZG(:)

! These arrays store the longitude and latitude of the lower left corner of
! each of the dipole u quadrilaterals and t quadrilaterals.
real(r8), allocatable :: ULAT(:,:,:), ULON(:,:,:), ULEV(:,:,:)
real(r8), allocatable :: VLAT(:,:,:), VLON(:,:,:), VLEV(:,:,:)
real(r8), allocatable :: WLAT(:,:,:), WLON(:,:,:), WLEV(:,:,:)
real(r8), allocatable ::  LAT(:,:,:),  LON(:,:,:),  LEV(:,:,:)

type(location_type), allocatable :: state_locations(:)
integer,             allocatable :: state_kinds(:)
type(get_close_type) :: gc_state

! integer, lowest valid cell number in the vertical
integer, allocatable  :: KMT(:, :), KMU(:, :)

! compute pressure based on depth - can do once upfront.
real(r8), allocatable :: pressure(:)

type(time_type) :: model_time, model_timestep

integer :: model_size    ! the state vector length

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME) :: coordinates
   character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimnames
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer  :: numdims
   integer  :: rank        ! excludes the (singleton) time dimension if present
   integer  :: maxlevels
   integer  :: xtype
   integer  :: varsize     ! prod(dimlens(1:numdims))
   integer  :: index1      ! location in dart state vector of first occurrence
   integer  :: indexN      ! location in dart state vector of last  occurrence
   integer  :: dart_kind
   integer  :: rangeRestricted
   real(r8) :: minvalue
   real(r8) :: maxvalue
   integer  :: spvalINT, missingINT
   real(r4) :: spvalR4, missingR4
   real(r8) :: spvalR8, missingR8
   logical  :: has_fill_value      ! intended for future use
   logical  :: has_missing_value   ! intended for future use
   character(len=obstypelength) :: kind_string
   character(len=512) :: origin    ! the file it came from
   logical  :: update
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

INTERFACE vector_to_variable
      MODULE PROCEDURE vector_to_1d_variable
      MODULE PROCEDURE vector_to_2d_variable
      MODULE PROCEDURE vector_to_3d_variable
END INTERFACE

INTERFACE DART_get_var
      MODULE PROCEDURE get_var_1d
      MODULE PROCEDURE get_var_2d
      MODULE PROCEDURE get_var_3d
END INTERFACE

!------------------------------------------------

! The regular grid used for dipole interpolation divides the sphere into
! a set of regularly spaced lon-lat boxes. The number of boxes in
! longitude and latitude are set by num_reg_x and num_reg_y. Making the
! number of regular boxes smaller decreases the computation required for
! doing each interpolation but increases the static storage requirements
! and the initialization computation (which seems to be pretty small).
integer, parameter :: num_reg_x = 90, num_reg_y = 90

! The max_reg_list_num controls the size of temporary storage used for
! initializing the regular grid. Four arrays
! of size num_reg_x*num_reg_y*max_reg_list_num are needed. The initialization
! fails and returns an error if max_reg_list_num is too small. A value of
! 30 is sufficient for the x3 GCOM grid with 180 regular lon and lat boxes
! and a value of 80 is sufficient for for the x1 grid.
integer, parameter :: max_reg_list_num = 80

! The dipole interpolation keeps a list of how many and which dipole quads
! overlap each regular lon-lat box. The number for the u and t grids are
! stored in u_dipole_num and t_dipole_num. The allocatable arrays
! u_dipole_lon(lat)_list and t_dipole_lon(lat)_list list the longitude
! and latitude indices for the overlapping dipole quads. The entry in
! u_dipole_start and t_dipole_start for a given regular lon-lat box indicates
! where the list of dipole quads begins in the u_dipole_lon(lat)_list and
! t_dipole_lon(lat)_list arrays.

integer :: u_dipole_start(num_reg_x, num_reg_y)
integer :: u_dipole_num  (num_reg_x, num_reg_y) = 0
integer :: t_dipole_start(num_reg_x, num_reg_y)
integer :: t_dipole_num  (num_reg_x, num_reg_y) = 0
integer, allocatable :: u_dipole_lon_list(:), t_dipole_lon_list(:)
integer, allocatable :: u_dipole_lat_list(:), t_dipole_lat_list(:)

! Need to check for pole quads: for now we are not interpolating in them
integer :: pole_x, t_pole_y, u_pole_y

! Have a global variable saying whether this is dipole or regular lon-lat grid
! This should be initialized static_init_model. Code to do this is below.
logical :: dipole_grid

contains

!=======================================================================
! All of the REQUIRED public interfaces come first - by convention.
!=======================================================================


!-----------------------------------------------------------------------
!>
!> Called to do one time initialization of the model.
!> In this case, it reads in the grid information.

subroutine static_init_model()

integer :: iunit, io
integer :: ss, dd

if ( module_initialized ) return ! only need to do this once.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat=io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
call error_handler(E_MSG,'static_init_model','model_nml values are')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output()) write(     *     , nml=model_nml)

call set_calendar_type('Gregorian')

!-----------------------------------------------------------------------
! Set the time step ...
! FIXME ensure model_timestep is multiple of 'param.nml:Writeout_freq'

model_timestep = set_time(assimilation_period_seconds,assimilation_period_days)
call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)
write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1)

!-----------------------------------------------------------------------
! get data dimensions, then allocate space, then open the files
! and actually fill in the arrays.

call get_grid_dims()
call read_grid()

if (do_output() .and. (debug > 20)) call write_grid_interptest()

!---------------------------------------------------------------
! Compile the list of variables to use in the creation
! of the DART state vector.

call parse_variable_table( gcom_variables, nfields, variable_table )

call fill_progvar()

call set_state_locations_kinds()

!-----------------------------------------------------------------------
! compute the offsets into the state vector for the start of each
! different variable type.

if (do_output() .and. (debug > 0)) then
   write(string1,*) '    grid center shape : nx, ny, nz = ', nx, ny, nz
   write(string2,*) 'grid edges  shape : nxp1, nyp1, nzp1 = ', nxp1, nyp1, nzp1
   write(string3,*) 'number of variables = ', nfields, ' model_size = ', model_size
   call error_handler(E_MSG,'static_init_model',string1,text2=string2,text3=string3)
endif

! initialize the pressure array - pressure in bars
! allocate(pressure(nzp1))
! call depth2pressure(nzp1, ZC, pressure)

! Initialize the interpolation routines
call init_interp()

end subroutine static_init_model


!-----------------------------------------------------------------------
!>
!> Returns the size of the model as an integer, i.e.
!> the length of everything you need to pack into a DART vector:
!> the 'state', any parameters you are estimating, etc.
!>
!> Required for all applications.
!>

function get_model_size()
integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!-----------------------------------------------------------------------
!>
!> Returns the minimum amount of time that the model should be
!> advanced in a given implementation.
!>
!> This interface is required for all applications.

function get_model_time_step()
type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

call error_handler(E_MSG,'get_model_time_step','FIXME TJH UNTESTED', &
    text2='should be a multiple of param.dat:&param_nml:Writeout_freq')

end function get_model_time_step


!-----------------------------------------------------------------------
!>
!> Given an integer index into the state vector structure, returns the
!> associated location. A second intent(out) optional argument kind
!> can be returned if the model has more than one type of field (for
!> instance temperature and zonal wind component). This interface is
!> required for all filter applications as it is required for computing
!> the distance between observations and state variables.
!>
!> This interface is required for all applications.

subroutine get_state_meta_data(index_in, location, var_type)

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,  OPTIONAL,  intent(out) :: var_type

integer  :: varindex
integer  :: lon_index, lat_index, lev_index
real(r8) :: mylon, mylat, mylev

if ( .not. module_initialized ) call static_init_model

varindex  = -1 ! if varindex is negative, get_state_indices will calculate it
call get_state_indices(index_in, lon_index, lat_index, lev_index, varindex)

call get_state_lonlatlev(varindex, lon_index, lat_index, lev_index, &
                               mylon, mylat, mylev)

location = set_location(mylon, mylat, mylev, VERTISHEIGHT)

if (present(var_type)) then
   var_type = progvar(varindex)%dart_kind
endif

if (do_output() .and. (debug > 5)) then
   if     ( progvar(varindex)%rank == 1) then
      write(*,100) trim(progvar(varindex)%varname), index_in, lon_index 
   elseif ( progvar(varindex)%rank == 2) then
      write(*,200) trim(progvar(varindex)%varname), index_in, lon_index, lat_index 
   elseif ( progvar(varindex)%rank == 3) then
      if (mylon > 180.0_r8) mylon =  mylon - 360.0_r8
      write(*,300) trim(progvar(varindex)%varname), index_in, lon_index, lat_index , lev_index, mylon, mylat, mylev
   endif
endif

 100 format(A,2(1x,i10))              ! varname, DART index, i
 200 format(A,3(1x,i10))              ! varname, DART index, i, j
 300 format(A,4(1x,i10),3(1x,f18.10)) ! varname, DART index, i, j, k, and the 3 locations

end subroutine get_state_meta_data


!-----------------------------------------------------------------------
!>
!> model_interpolate is the basis for almost all forward observation operators.
!> the given location given a state vector. The type of the variable being
!> interpolated is obs_type since normally this is used to find the expected
!> value of an observation at some location. The interpolated value is
!> returned in interp_val and istatus is 0 for success.
!>
!> This interface is required for all applications.

subroutine model_interpolate(x, location, itype, obs_val, istatus)

real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: itype  ! really a KIND
real(r8),            intent(out) :: obs_val
integer,             intent(out) :: istatus

integer  :: close_ind(model_size) ! undesirable to have something this big ... but ...

! Local storage
real(r8) :: loc_array(3), longitude, latitude, level
integer  :: base_kind
integer  :: i, ivar, varindex, iclose, num_close, indx, num_wanted
integer  :: closest_index, lon_index, lat_index, lev_index
real(r8) :: closest
! real(r8) :: above, below, dist_above, dist_below
! integer  :: istatus1, istatus2

real(r8), allocatable :: distances(:), sorted_distances(:), data_values(:)
integer,  allocatable :: indices(:), sorted_indices(:), close_indices(:)

integer, parameter :: num_neighbors = 12

character(len=obstypelength) :: kind_name

real(r8), parameter :: RAD2M = 40000000.0_r8/(2.0_r8*PI)

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the obs_val will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

obs_val   = MISSING_R8   ! the DART bad value flag
istatus   = 99           ! unknown error
base_kind = itype        ! it really is a DART KIND, so call it a kind
varindex  = -1

! Get the individual location values
loc_array = get_location(location)
longitude = loc_array(1)
latitude  = loc_array(2)
level     = loc_array(3)

if (do_output() .and. (debug > 4))  &
   print *,'want interpolation of DART KIND',base_kind,'at',longitude,latitude,level

kind_name = get_raw_obs_kind_name(base_kind)

! If the requested kind does not exist in our state - just return as a failure
VARLOOP : do ivar = 1,nfields
   if (progvar(ivar)%dart_kind == base_kind) then
      varindex = ivar
      exit VARLOOP
   endif
enddo VARLOOP
if (varindex < 0) return

! Generate the list of indices into the DART vector of the close candidates.

call loc_get_close_obs(gc_state, location, base_kind, state_locations, state_kinds, &
                       num_close, close_ind)

! Sometimes the location is outside the model domain.
! In this case, we cannot interpolate.

if (num_close == 0) then
   istatus = 23
   return
endif

! Loop over close candidates. They come in without regard to what DART KIND they are,
! nor are they sorted by distance. We are only interested in the close locations
! of the right DART KIND.

closest = 1000000.0_r8 ! not very close
closest_index = 0

allocate(distances(num_close),indices(num_close))
num_wanted = 0
CLOSE : do iclose = 1, num_close

   indx = close_ind(iclose)

   if (state_kinds(indx) /= base_kind) cycle CLOSE

   num_wanted = num_wanted + 1
   kind_name  = get_raw_obs_kind_name(state_kinds(indx))
   distances(num_wanted) = get_dist(location, state_locations(indx))
   indices(num_wanted)   = indx

   if (distances(num_wanted) < closest) then
      closest = distances(num_wanted)
      closest_index = indx
   endif

enddo CLOSE

!FIX ! FIXME ... can use the simple inverse_distance below or do this. no need to do both
!FIX !
!FIX ! Given the index of the closest, we can calculate the indices of the neighbors.
!FIX ! Might be faster than sorting several hundred items ...
!FIX ! Might be able to just search the N nearest neighbors.
!FIX
!FIX call get_state_indices(closest_index, lon_index, lat_index, lev_index, varindex)
!FIX
!FIX if (do_output() .and. debug > 0) then
!FIX    write(*,*)'model_interpolate: closest_index, distance',closest_index,closest
!FIX    write(*,*)'model_interpolate: [i,j,k] of closest is',lon_index,lat_index,lev_index
!FIX endif
!FIX
!FIX There are remnants of the POP interpolation that may be useful here.
!FIX lon_lat_interpolate()
!FIX do_interp()
!FIX
!FIX horizontal_interpolation() is simply a stub ...
!FIX call horizontal_interpolation(location, closest_index, num_wanted, &
!FIX        indices(1:num_wanted), distances(1:num_wanted), varindex, x, 'above', &
!FIX        above, dist_above, istatus1)
!FIX call horizontal_interpolation(location, closest_index, num_wanted, &
!FIX        indices(1:num_wanted), distances(1:num_wanted), varindex, x, 'below', &
!FIX        below, dist_below, istatus2)

! The inverse distance weighted scheme should only use a small number of neighbors.
! Now we can sort based on distance.
! Must also keep track of the location into the DART state vector so we can use
! them to get the state values at those locations (that are close enough).

! FIXME ... what is happening is that the vertical levels are closely spaced relative
! to the horizontal separation. Since we have found the closest level, we
! should horizontally interpolate on the levels above and below and then do
! a vertical interpolation.

! TJH I played around with the order of magnitude of vert_normalization_height
! TJH from 6370000 to 6370 and looked at the effect on which were the closest gridcells.
! TJH with 6370000 ... the closest were all the same horizontal i,j and only k changed
! TJH with 637000  ... ditto
! TJH with 63700   ... ditto
! TJH with 6370    ... FINALLY the vertical stayed the same and the nearest neighbors
! TJH were the horizontal neighbors.

allocate(sorted_distances(num_wanted), &
           sorted_indices(num_wanted), &
            close_indices(num_wanted), &
              data_values(num_wanted))

sorted_indices = (/ (i,i=1,num_wanted) /)

call index_sort(distances(1:num_wanted), sorted_indices, num_wanted)

do i=1,min(num_wanted,num_neighbors)
   sorted_distances(i) = distances(sorted_indices(i))
      close_indices(i) =   indices(sorted_indices(i))
        data_values(i) = x(indices(sorted_indices(i)))
enddo

if (do_output() .and. (debug > 4)) then
   write(*,*)'       rank   dartindex   distance(meters)           state' &
             //'_value                  i           j           k'
   do i=1,min(num_wanted,num_neighbors)
      call get_state_indices(close_indices(i), lon_index, lat_index, lev_index, varindex)
      write(*,*)i,close_indices(i),sorted_distances(i)*RAD2M,data_values(i),&
                lon_index, lat_index, lev_index, &
                ULON(lon_index, lat_index, lev_index), &
                ULAT(lon_index, lat_index, lev_index), &
                ULEV(lon_index, lat_index, lev_index)
   enddo
endif

call inverse_distance_interpolation(data_values, sorted_distances, num_neighbors, &
         obs_val, istatus)

if (do_output() .and. (debug > 4)) then
    write(*,*)
    print *, 'model_interpolate: summary'
    print *, 'itype     ', itype, ' is ',trim(kind_name)
    print *, 'num_close ', num_close
    print *, 'closest   ', close_ind(1)
    print *, 'obs_val   ', obs_val
    print *, 'istatus   ', istatus
endif

deallocate(distances, sorted_distances, data_values)
deallocate(indices, sorted_indices, close_indices)

end subroutine model_interpolate


!-----------------------------------------------------------------------
!>
!> Writes the model-specific (static) attributes to a netCDF file.
!> This includes coordinate variables and some metadata, but NOT
!> the model state.
!>
!> assim_model_mod:init_diag_output uses information from the location_mod
!>     to define the location dimension and variable ID. All we need to do
!>     is query, verify, and fill ...
!>
!> This interface is required for all applications. no need to do both

function nc_write_model_atts( ncFileID ) result (ierr)

integer, intent(in) :: ncFileID      ! netCDF file identifier
integer             :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!----------------------------------------------------------------------
! variables if we just blast out one long state vector
!----------------------------------------------------------------------

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

!----------------------------------------------------------------------
! variables if we parse the state vector into prognostic variables.
!----------------------------------------------------------------------

! for the dimensions and coordinate variables
integer :: nxp1_dimid, nyp1_dimid, nzp1_dimid
integer ::   nx_dimid,   ny_dimid,   nz_dimid
integer :: ulonVarID, ulatVarID, ulevVarID
integer :: vlonVarID, vlatVarID, vlevVarID
integer :: wlonVarID, wlatVarID, wlevVarID
integer ::  lonVarID,  latVarID,  levVarID ! [T,S,P]

! for the prognostic variables
integer :: ivar, VarID

!----------------------------------------------------------------------
! variables for the namelist output
!----------------------------------------------------------------------

character(len=129), allocatable :: textblock(:)
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_gcom_namelist

!----------------------------------------------------------------------
! local variables
!----------------------------------------------------------------------

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)  :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10) :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)  :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer           :: values(8)   ! needed by F90 DATE_AND_TIME intrinsic

character(len=NF90_MAX_NAME) :: str1
character(len=NF90_MAX_NAME) :: varname
integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids

integer :: i, myndims

character(len=256) :: filename

! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
! and then put into define mode.
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                   'nc_write_model_atts', 'inquire '//trim(filename))
call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts',   'redef '//trim(filename))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension.
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='NMLlinelen', dimid=LineLenDimID), &
                           'nc_write_model_atts','inq_dimid NMLlinelen')

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
                           'nc_write_model_atts', 'copy dimid '//trim(filename))

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
                           'nc_write_model_atts', 'time dimid '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
             ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
        dimid = StateVarDimID),'nc_write_model_atts', 'state def_dim '//trim(filename))

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,str1    ), &
           'nc_write_model_atts', 'creation put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source  ), &
           'nc_write_model_atts', 'source put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
           'nc_write_model_atts', 'revision put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate ), &
           'nc_write_model_atts', 'revdate put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'GCOM' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

!-------------------------------------------------------------------------------
! Determine shape of most important GCOM control file
!-------------------------------------------------------------------------------

call find_textfile_dims('param.dat', nlines, linelen)
if (nlines > 0) then
  has_gcom_namelist = .true.
else
  has_gcom_namelist = .false.
endif

if (has_gcom_namelist) then
   allocate(textblock(nlines))
   textblock = ''

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nlines', &
                 len = nlines, dimid = nlinesDimID), &
                 'nc_write_model_atts', 'def_dim nlines ')

   call nc_check(nf90_def_var(ncFileID,name='param', xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'nc_write_model_atts', 'def_var param')
   call nc_check(nf90_put_att(ncFileID, nmlVarID, 'long_name',       &
                 'contents of param.dat control file'), 'nc_write_model_atts', 'put_att param')

endif

!-------------------------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.
!-------------------------------------------------------------------------------

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create a variable for the state vector
   !----------------------------------------------------------------------------

  ! Define the state vector coordinate variable and some attributes.
   call nc_check(nf90_def_var(ncid=ncFileID,name='StateVariable', xtype=nf90_int, &
                 dimids=StateVarDimID, varid=StateVarVarID), 'nc_write_model_atts', &
                 'statevariable def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,'long_name','State Variable ID'),&
                 'nc_write_model_atts','statevariable long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'units','indexical'), &
                 'nc_write_model_atts', 'statevariable units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,'valid_range',(/ 1,model_size /)),&
                 'nc_write_model_atts', 'statevariable valid_range '//trim(filename))

   ! Define the actual (3D) state vector, which gets filled as time goes on ...
   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                 dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                 'nc_write_model_atts','state def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'long_name','model state or fcopy'),&
                 'nc_write_model_atts', 'state long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'_FillValue',MISSING_R8),&
                 'nc_write_model_atts', 'state FillValue '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'missing_value',MISSING_R8),&
                 'nc_write_model_atts', 'state missing_value '//trim(filename))

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts','state enddef '//trim(filename))

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                 'nc_write_model_atts', 'state put_var '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to output the prognostic variables.
   !----------------------------------------------------------------------------
   ! Define the new dimensions IDs
   ! FIXME ... do we grab dimension names from all possible vars
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nxp1', &
          len=nxp1,dimid=nxp1_dimid),'nc_write_model_atts','nxp1 def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='nyp1', &
          len=nyp1,dimid=nyp1_dimid),'nc_write_model_atts','nyp1 def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='nzp1', &
          len=nzp1,dimid=nzp1_dimid),'nc_write_model_atts','nzp1 def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nx', &
          len=nx,dimid=nx_dimid),'nc_write_model_atts','nx def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='ny', &
          len=ny,dimid=ny_dimid),'nc_write_model_atts','ny def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='nz', &
          len=nz,dimid=nz_dimid),'nc_write_model_atts','nz def_dim '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! U Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='ulon', xtype=nf90_real, &
                 dimids=(/ nxp1_dimid, ny_dimid, nz_dimid /), varid=ulonVarID),&
                 'nc_write_model_atts', 'ulon def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'long_name', 'longitude'), &
                 'nc_write_model_atts', 'ulon long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'comment', 'longitudes of U grid'), &
                 'nc_write_model_atts', 'ulon comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'axis', 'X'),  &
                 'nc_write_model_atts', 'ulon axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'ulon units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'ulon valid_range '//trim(filename))

   ! U Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='ulat', xtype=nf90_real, &
                 dimids=(/ nxp1_dimid, ny_dimid, nz_dimid /), varid=ulatVarID),&
                 'nc_write_model_atts', 'ulat def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'long_name', 'latitude'), &
                 'nc_write_model_atts', 'ulat long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'comment', 'latitudes of U grid'), &
                 'nc_write_model_atts', 'ulat comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'axis', 'Y'),   &
                 'nc_write_model_atts', 'ulat axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'ulat units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'ulat valid_range '//trim(filename))

   ! U Grid Levels
   call nc_check(nf90_def_var(ncFileID,name='ulev', xtype=nf90_real, &
                 dimids=(/ nxp1_dimid, ny_dimid, nz_dimid /), varid=ulevVarID),&
                 'nc_write_model_atts', 'ulev def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulevVarID, 'long_name', 'depth'), &
                 'nc_write_model_atts', 'ulev long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulevVarID, 'comment', 'depth of U grid'), &
                 'nc_write_model_atts', 'ulev comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulevVarID, 'axis', 'Z'),   &
                 'nc_write_model_atts', 'ulev axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulevVarID, 'units', 'meters'),  &
                 'nc_write_model_atts', 'ulev units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulevVarID, 'positive', 'up'),  &
                 'nc_write_model_atts', 'ulev positive '//trim(filename))

   !----------------------------------------------------------------------------
   ! V Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='vlon', xtype=nf90_real, &
                 dimids=(/ nx_dimid, nyp1_dimid, nz_dimid /), varid=vlonVarID),&
                 'nc_write_model_atts', 'vlon def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'long_name', 'longitude'), &
                 'nc_write_model_atts', 'vlon long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'comment', 'longitudes of V grid'), &
                 'nc_write_model_atts', 'vlon comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'axis', 'X'),  &
                 'nc_write_model_atts', 'vlon axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'vlon units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'vlon valid_range '//trim(filename))

   ! V Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='vlat', xtype=nf90_real, &
                 dimids=(/ nx_dimid, nyp1_dimid, nz_dimid /), varid=vlatVarID),&
                 'nc_write_model_atts', 'vlat def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'long_name', 'latitude'), &
                 'nc_write_model_atts', 'vlat long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'comment', 'latitudes of V grid'), &
                 'nc_write_model_atts', 'vlat comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'axis', 'Y'),   &
                 'nc_write_model_atts', 'vlat axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'vlat units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlatVarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'vlat valid_range '//trim(filename))

   ! V Grid Levels
   call nc_check(nf90_def_var(ncFileID,name='vlev', xtype=nf90_real, &
                 dimids=(/ nx_dimid, nyp1_dimid, nz_dimid /), varid=vlevVarID),&
                 'nc_write_model_atts', 'vlev def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlevVarID, 'long_name', 'depth'), &
                 'nc_write_model_atts', 'vlev long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlevVarID, 'comment', 'depth of V grid'), &
                 'nc_write_model_atts', 'vlev comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlevVarID, 'axis', 'Z'),   &
                 'nc_write_model_atts', 'vlev axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlevVarID, 'units', 'meters'),  &
                 'nc_write_model_atts', 'vlev units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  vlevVarID, 'positive', 'up'),  &
                 'nc_write_model_atts', 'vlev positive '//trim(filename))

   !----------------------------------------------------------------------------
   ! W Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='wlon', xtype=nf90_real, &
                 dimids=(/ nx_dimid, ny_dimid, nzp1_dimid /), varid=wlonVarID),&
                 'nc_write_model_atts', 'wlon def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlonVarID, 'long_name', 'longitude'), &
                 'nc_write_model_atts', 'wlon long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlonVarID, 'comment', 'longitudes of W grid'), &
                 'nc_write_model_atts', 'wlon comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlonVarID, 'axis', 'X'),  &
                 'nc_write_model_atts', 'wlon axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlonVarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'wlon units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'wlon valid_range '//trim(filename))

   ! W Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='wlat', xtype=nf90_real, &
                 dimids=(/ nx_dimid, ny_dimid, nzp1_dimid /), varid=wlatVarID),&
                 'nc_write_model_atts', 'wlat def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlatVarID, 'long_name', 'latitude'), &
                 'nc_write_model_atts', 'wlat long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlatVarID, 'comment', 'latitudes of W grid'), &
                 'nc_write_model_atts', 'wlat comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlatVarID, 'axis', 'Y'),   &
                 'nc_write_model_atts', 'wlat axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlatVarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'wlat units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlatVarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'wlat valid_range '//trim(filename))

   ! W Grid Levels
   call nc_check(nf90_def_var(ncFileID,name='wlev', xtype=nf90_real, &
                 dimids=(/ nx_dimid, ny_dimid, nzp1_dimid /), varid=wlevVarID),&
                 'nc_write_model_atts', 'wlev def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'long_name', 'depth'), &
                 'nc_write_model_atts', 'wlev long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'comment', 'levels of W grid'), &
                 'nc_write_model_atts', 'wlev comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'axis', 'Z'),   &
                 'nc_write_model_atts', 'wlev axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'units', 'meters'),  &
                 'nc_write_model_atts', 'wlev units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'positive', 'up'),  &
                 'nc_write_model_atts', 'wlev positive '//trim(filename))

   !----------------------------------------------------------------------------
   ! [T,S,p] Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lon', xtype=nf90_real, &
                 dimids=(/ nx_dimid, ny_dimid, nz_dimid /), varid=lonVarID),&
                 'nc_write_model_atts', 'lon def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  lonVarID, 'long_name', 'longitude'), &
                 'nc_write_model_atts', 'lon long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  lonVarID, 'comment', 'longitudes of [T,S,p] grid'), &
                 'nc_write_model_atts', 'lon comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  lonVarID, 'axis', 'X'),  &
                 'nc_write_model_atts', 'lon axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  lonVarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lon units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  lonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'lon valid_range '//trim(filename))

   ! [T,S,p] Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='lat', xtype=nf90_real, &
                 dimids=(/ nx_dimid, ny_dimid, nz_dimid /), varid=latVarID),&
                 'nc_write_model_atts', 'lat def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  latVarID, 'long_name', 'latitude'), &
                 'nc_write_model_atts', 'lat long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  latVarID, 'comment', 'latitudes of [T,S,p] grid'), &
                 'nc_write_model_atts', 'lat comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  latVarID, 'axis', 'Y'),   &
                 'nc_write_model_atts', 'lat axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  latVarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'lat units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  latVarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'lat valid_range '//trim(filename))

   ! [T,S,p] Grid Levels
   call nc_check(nf90_def_var(ncFileID,name='lev', xtype=nf90_real, &
                 dimids=(/ nx_dimid, ny_dimid, nz_dimid /), varid=levVarID),&
                 'nc_write_model_atts', 'lev def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  levVarID, 'long_name', 'depth'), &
                 'nc_write_model_atts', 'lev long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  levVarID, 'comment', 'levels of [T,S,p] grid'), &
                 'nc_write_model_atts', 'lev comment '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  levVarID, 'axis', 'Z'),   &
                 'nc_write_model_atts', 'lev axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  levVarID, 'units', 'meters'),  &
                 'nc_write_model_atts', 'lev units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  levVarID, 'positive', 'up'),  &
                 'nc_write_model_atts', 'lev positive '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(ivar, ncFileID, MemberDimID, unlimitedDimID, myndims, mydimids)

      ! define the variable and set the attributes

      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), &
                    xtype=progvar(ivar)%xtype, &
                    dimids = mydimids(1:myndims), varid=VarID),&
                    'nc_write_model_atts', trim(string1)//' def_var' )

      call nc_check(nf90_put_att(ncFileID, VarID, &
              'long_name', trim(progvar(ivar)%long_name)), &
              'nc_write_model_atts', trim(string1)//' put_att long_name' )
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'DART_kind', trim(progvar(ivar)%kind_string)), &
              'nc_write_model_atts', trim(string1)//' put_att dart_kind' )
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'units', trim(progvar(ivar)%units)), &
              'nc_write_model_atts', trim(string1)//' put_att units' )
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'original_coordinates', trim(progvar(ivar)%coordinates)), &
              'nc_write_model_atts', trim(string1)//' put_att coordinates' )

      ! Preserve the original missing_value/_FillValue code.

      if (  progvar(ivar)%xtype == NF90_INT ) then
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 'missing_value', progvar(ivar)%spvalINT), &
                 'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 '_FillValue',  progvar(ivar)%spvalINT), &
                 'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

      elseif (  progvar(ivar)%xtype == NF90_FLOAT ) then
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 'missing_value', progvar(ivar)%spvalR4), &
                 'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 '_FillValue',  progvar(ivar)%spvalR4), &
                 'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

      elseif (  progvar(ivar)%xtype == NF90_DOUBLE ) then
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 'missing_value', progvar(ivar)%spvalR8), &
                 'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 '_FillValue',  progvar(ivar)%spvalR8), &
                 'nc_write_model_atts', trim(string1)//' put_att _FillValue' )
      endif

   enddo

   !----------------------------------------------------------------------------
   ! Finished with dimension/variable definitions, must end 'define' mode to fill.
   ! Fill the coordinate variables
   !----------------------------------------------------------------------------

   call nc_check(nf90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

   call nc_check(nf90_put_var(ncFileID, ulonVarID, ULON ), &
                'nc_write_model_atts', 'ulon put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ulatVarID, ULAT ), &
                'nc_write_model_atts', 'ulat put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ulevVarID, ULEV ), &
                'nc_write_model_atts', 'ulev put_var '//trim(filename))

   call nc_check(nf90_put_var(ncFileID, vlonVarID, VLON ), &
                'nc_write_model_atts', 'vlon put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, vlatVarID, VLAT ), &
                'nc_write_model_atts', 'vlat put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, vlevVarID, VLEV ), &
                'nc_write_model_atts', 'vlev put_var '//trim(filename))

   call nc_check(nf90_put_var(ncFileID, wlonVarID, WLON ), &
                'nc_write_model_atts', 'wlon put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, wlatVarID, WLAT ), &
                'nc_write_model_atts', 'wlat put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, wlevVarID, WLEV ), &
                'nc_write_model_atts', 'wlev put_var '//trim(filename))

   call nc_check(nf90_put_var(ncFileID, lonVarID, LON ), &
                'nc_write_model_atts', 'lon put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, latVarID, LAT ), &
                'nc_write_model_atts', 'lat put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, levVarID, LEV ), &
                'nc_write_model_atts', 'lev put_var '//trim(filename))

endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

if (has_gcom_namelist) then
   call file_to_text('param.dat', textblock)
   call nc_check(nf90_put_var(ncFileID, nmlVarID, textblock ), &
                 'nc_write_model_atts', 'put_var param.dat')
   deallocate(textblock)
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts


!-----------------------------------------------------------------------
!>
!> Writes the model variables to a netCDF file.
!>
!> All errors are fatal, so the return code is always '0 == normal',
!> since the fatal errors stop execution.
!>
!> assim_model_mod:init_diag_output uses information from the location_mod
!>     to define the location dimension and variable ID. All we need to do
!>     is query, verify, and fill ...
!>
!> This interface is required for all applications.

function nc_write_model_vars( ncFileID, state_vector, copyindex, timeindex ) result (ierr)

integer,  intent(in) :: ncFileID      ! netCDF file identifier
real(r8), intent(in) :: state_vector(:)
integer,  intent(in) :: copyindex
integer,  intent(in) :: timeindex
integer              :: ierr          ! return value of function

! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, ncstart, nccount
character(len=NF90_MAX_NAME)          :: varname
integer :: i, ivar, VarID, ncNdims, dimlen
integer :: TimeDimID, CopyDimID
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

real(r8), allocatable :: data_1d(:)
real(r8), allocatable :: data_2d(:,:)
real(r8), allocatable :: data_3d(:,:,:)

character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
              'nc_write_model_vars', 'inquire '//trim(filename))

call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
            'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
            'nc_write_model_vars', 'inq_dimid time '//trim(filename))

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, 'state', VarID), &
                 'nc_write_model_vars', 'state inq_varid '//trim(filename))
   call nc_check(NF90_put_var(ncFileID,VarID,state_vector,start=(/1,copyindex,timeindex/)),&
                 'nc_write_model_vars', 'state put_var '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! Reshape the DART vector into the 'natural' shape of the variables.
   ! Replace DART missing values with netcdf missing value.
   !----------------------------------------------------------------------------

   do ivar = 1,nfields

      varname = trim(progvar(ivar)%varname)
      string2 = trim(filename)//' '//trim(varname)

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried
      ! by looping over the dimensions stored in the progvar type.

      call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'nc_write_model_vars', 'inq_varid '//trim(string2))

      call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'nc_write_model_vars', 'inquire '//trim(string2))

      ncstart = 1   ! These are arrays, actually
      nccount = 1
      DimCheck : do i = 1,progvar(ivar)%numdims

         write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
               'nc_write_model_vars', trim(string1))

         if (progvar(ivar)%dimnames(i) == 'time') cycle DimCheck

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*)trim(string2),' dim/dimlen ',i,dimlen, &
                            ' not ',progvar(ivar)%dimlens(i)
            write(string2,*)' but it should be.'
            call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
                            source, revision, revdate, text2=trim(string2))
         endif

         nccount(i) = dimlen

      enddo DimCheck

      where(dimIDs == CopyDimID) ncstart = copyindex
      where(dimIDs == CopyDimID) nccount = 1
      where(dimIDs == TimeDimID) ncstart = timeindex
      where(dimIDs == TimeDimID) nccount = 1

      if (do_output() .and. (debug > 99)) then
         write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',ncstart(1:ncNdims)
         write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',nccount(1:ncNdims)
      endif

      ! Since 'time' is a singleton dimension, we can use the same logic
      ! as if it the variable had one less dimension.

      if (     progvar(ivar)%rank == 1 ) then

         if ( ncNdims /= 3 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'DART netcdf variable should have 3 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_1d( progvar(ivar)%dimlens(1) ))
         call vector_to_variable(state_vector, ivar, data_1d)
         call nc_check(nf90_put_var(ncFileID, VarID, data_1d, &
             start = ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_1d)

      elseif ( progvar(ivar)%rank == 2 ) then

         if ( ncNdims /= 4 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'DART netcdf variable should have 4 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_2d( progvar(ivar)%dimlens(1),  &
                           progvar(ivar)%dimlens(2) ))
         call vector_to_variable(state_vector, ivar, data_2d)
         call nc_check(nf90_put_var(ncFileID, VarID, data_2d, &
             start = ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_2d)

      elseif ( progvar(ivar)%rank == 3 ) then

         if ( ncNdims /= 5 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'DART netcdf variable should have 5 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_3d( progvar(ivar)%dimlens(1),  &
                           progvar(ivar)%dimlens(2),  &
                           progvar(ivar)%dimlens(3) ))
         call vector_to_variable(state_vector, ivar, data_3d)
         call nc_check(nf90_put_var(ncFileID, VarID, data_3d, &
             start = ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_3d)

      else

         write(string1,*)'Do not know how to handle GCOM variables higher than rank 3'
         write(string2,*)trim(progvar(ivar)%varname),'has shape', &
                              progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
         call error_handler(E_ERR,'nc_write_model_vars',string1,source,revision,revdate)

      endif

   enddo

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync '//trim(filename))

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars


!-----------------------------------------------------------------------
!>
!> Companion interface to init_conditions. Returns a time that is
!> appropriate for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.
!>
!> This interface is required for all applications.

subroutine init_time(time)

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

string2 = "cannot run perfect_model_obs with 'start_from_restart = .false.' "
string3 = 'use gcom_to_dart to generate an initial state which contains a timestamp'
call error_handler(E_ERR,'init_time', &
                  'GCOM model_mod has no built-in default time', &
                  source, revision, revdate, &
                  text2=string2, text3=string3)

! this code never reached - just here to avoid compiler warnings
! about an intent(out) variable not being set to a value.
time = set_time(0,0)

end subroutine init_time


!-----------------------------------------------------------------------
!>
!> Returns a model state vector, x, that is some sort of appropriate
!> initial condition for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.
!>
!> This interface is required for all applications.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

string2 = "cannot run perfect_model_obs with 'start_from_restart = .false.' "
string3 = 'use gcom_to_dart to generate an initial state'
call error_handler(E_ERR,'init_conditions', &
                  'GCOM model_mod has no built-in default state', &
                  source, revision, revdate, &
                  text2=string2, text3=string3)

! this code never reached - just here to avoid compiler warnings
! about an intent(out) variable not being set to a value.
x = 0.0_r8

end subroutine init_conditions


!-----------------------------------------------------------------------
!>
!> If the model could be called as a subroutine, does a single
!> timestep advance.  GCOM cannot be called this way, so fatal error
!> if this routine is called.
!>
!> This interface is required for all applications.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'adv_1step', &
                  'GCOM model cannot be called as a subroutine; async cannot = 0', &
                  source, revision, revdate)

x = MISSING_R8 ! unreachable, just silences compiler about unused variable.
! should do something to silence about time ...

end subroutine adv_1step


!-----------------------------------------------------------------------
!>
!> Perturbs a model state for generating initial ensembles.
!> The perturbed state is returned in pert_state.
!> A model may choose to provide a NULL INTERFACE by returning
!> .false. for the interf_provided argument. This indicates to
!> the filter that if it needs to generate perturbed states, it
!> may do so by adding a perturbation to each model state
!> variable independently. The interf_provided argument
!> should be returned as .true. if the model wants to do its own
!> perturbing of states.
!>
!> This interface is required for all applications.

subroutine pert_model_state(state, pert_state, interf_provided)

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

integer :: i, ivar
logical, save :: random_seq_init = .false.

if ( .not. module_initialized ) call static_init_model

call error_handler(E_MSG,'pert_model_state','FIXME Angie UNTESTED')

interf_provided = .true.

! Initialize my random number sequence
if(.not. random_seq_init) then
   call init_random_seq(random_seq, my_task_id())
   random_seq_init = .true.
endif

! only perturb the actual ocean cells; leave the land and
! ocean floor values alone.
VARLOOP : do ivar=1,nfields
   if (do_output() .and. (debug > 3)) then
      write(string1,*)'Perturbing ',ivar,trim(progvar(ivar)%varname)
      call error_handler(E_MSG,'pert_model_state',string1)
   endif

   STATELOOP : do i = progvar(ivar)%index1, progvar(ivar)%indexN

      pert_state(i) = state(i) + random_gaussian(random_seq, 0.0_r8, &
                                      model_perturbation_amplitude)

   enddo STATELOOP
enddo VARLOOP

end subroutine pert_model_state


!-----------------------------------------------------------------------
!>
!> Given a DART location (referred to as "base") and a set of candidate
!> locations & kinds (obs, obs_kind), returns the subset close to the
!> "base", their indices, and their distances to the "base" ...
!>
!> For vertical distance computations, general philosophy is to convert all
!> vertical coordinates to a common coordinate. This coordinate type is defined
!> in the namelist with the variable "vert_localization_coord".
!>
!> This interface is required for all applications.

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         obs, obs_kind, num_close, close_ind, dist)

type(get_close_type), intent(in) :: gc
type(location_type),  intent(in) :: base_obs_loc
integer,              intent(in) :: base_obs_kind
type(location_type),  intent(in) :: obs(:)
integer,              intent(in) :: obs_kind(:)
integer,              intent(out):: num_close
integer,              intent(out):: close_ind(:)
real(r8),  OPTIONAL,  intent(out):: dist(:)

integer :: t_ind, k

if ( .not. module_initialized ) call static_init_model

! Initialize variables to missing status

num_close = 0
close_ind = -99
if (present(dist)) dist = 1.0e9   !something big and positive (far away)

! Get all the potentially close obs but no dist (optional argument dist(:)
! is not present) This way, we are decreasing the number of distance
! computations that will follow.  This is a horizontal-distance operation and
! we don't need to have the relevant vertical coordinate information yet
! (for obs).
! FIXME - confirm that the close_ind() array does not benefit from having 
! all the dry_land locations pruned out.

call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs, obs_kind, &
                       num_close, close_ind)

! Loop over potentially close subset of obs priors or state variables
if (present(dist)) then
do k = 1, num_close

   t_ind = close_ind(k)

   ! if dry land, leave original LARGE value.  otherwise, compute real dist.
   if (obs_kind(t_ind) /= KIND_DRY_LAND) then
      dist(k) = get_dist(base_obs_loc,       obs(t_ind), &
                         base_obs_kind, obs_kind(t_ind))
   endif

enddo
endif

end subroutine get_close_obs


!-----------------------------------------------------------------------
!>
!> If needed by the model interface, this is the current mean
!> for all state vector items across all ensembles. It is up to this
!> code to allocate space and save a copy if it is going to be used
!> later on.  For now, we are ignoring it.
!>
!> This interface is required for all applications.

subroutine ens_mean_for_model(ens_mean)

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model

call error_handler(E_MSG,'ens_mean_for_model','FIXME TJH UNTESTED')

end subroutine ens_mean_for_model


!-----------------------------------------------------------------------
!>
!> Shutdown and clean-up.
!>
!> This interface is required for all applications.

subroutine end_model()

if ( .not. module_initialized ) call static_init_model

! assume if one is allocated, they all were.  if no one ever
! called the init routine, don't try to dealloc something that
! was never alloc'd.
if (allocated(ULAT)) deallocate(ULAT, ULON, ULEV)
if (allocated(VLAT)) deallocate(VLAT, VLON, VLEV)
if (allocated(WLAT)) deallocate(WLAT, WLON, WLEV)
if (allocated( LAT)) deallocate( LAT,  LON,  LEV)

if (allocated(state_locations)) deallocate(state_locations)
if (allocated(state_kinds))     deallocate(state_kinds)

call get_close_obs_destroy(gc_state)

end subroutine end_model


!=======================================================================
! End of the REQUIRED public interfaces.
!=======================================================================
! The remaining (optional) PUBLIC interfaces come next.
!=======================================================================


!-----------------------------------------------------------------------
!>
!> Read the grid size from the grid netcdf file.
!>
!> The file name comes from module storage ... namelist.
!>
!> nxp1, nx
!> nyp1, nx
!> nzp1, nz

subroutine get_grid_dims()

integer :: ncid, dimid

! get the ball rolling ...

call nc_check(nf90_open(trim(gcom_geometry_file), nf90_nowrite, ncid), &
         'get_grid_dims','open ['//trim(gcom_geometry_file)//']')

! Get the number of longitude gridcell centers

call nc_check(nf90_inq_dimid(ncid, 'nx', dimid), &
         'get_grid_dims', 'no "nx" in file ['//trim(gcom_geometry_file)//']')
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nx), &
         'get_grid_dims','inquire_dimension nx '//trim(gcom_geometry_file))

! Get the number of longitude gridcell edges

call nc_check(nf90_inq_dimid(ncid, 'nxp1', dimid), &
         'get_grid_dims', 'no "nxp1" in file ['//trim(gcom_geometry_file)//']')
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nxp1), &
         'get_grid_dims','inquire_dimension nxp1 '//trim(gcom_geometry_file))


! Latitudes

call nc_check(nf90_inq_dimid(ncid, 'ny', dimid), &
         'get_grid_dims', 'no "ny" in file ['//trim(gcom_geometry_file)//']')
call nc_check(nf90_inquire_dimension(ncid, dimid, len=ny), &
         'get_grid_dims','inquire_dimension ny '//trim(gcom_geometry_file))

call nc_check(nf90_inq_dimid(ncid, 'nyp1', dimid), &
         'get_grid_dims', 'no "nyp1" in file ['//trim(gcom_geometry_file)//']')
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nyp1), &
         'get_grid_dims','inquire_dimension nyp1 '//trim(gcom_geometry_file))

! Depth

call nc_check(nf90_inq_dimid(ncid, 'nz', dimid), &
         'get_grid_dims', 'no "nz" in file ['//trim(gcom_geometry_file)//']')
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nz), &
         'get_grid_dims','inquire_dimension nz '//trim(gcom_geometry_file))
call nc_check(nf90_inq_dimid(ncid, 'nzp1', dimid), &
         'get_grid_dims', 'no "nzp1" in file ['//trim(gcom_geometry_file)//']')
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nzp1), &
         'get_grid_dims','inquire_dimension nzp1 '//trim(gcom_geometry_file))

! tidy up

call nc_check(nf90_close(ncid), &
         'get_grid_dims','close '//trim(gcom_geometry_file) )

if (do_output() .and. (debug > 99)) then
   write(*,*)'get_grid_dims: filename  ['//trim(gcom_geometry_file)//']'
   write(*,*)'get_grid_dims: num longitude gridcell centers ',nx
   write(*,*)'get_grid_dims: num latitide  gridcell centers ',ny
   write(*,*)'get_grid_dims: num levels    gridcell centers ',nz
   write(*,*)'get_grid_dims: num longitude gridcell   edges ',nxp1
   write(*,*)'get_grid_dims: num latitide  gridcell   edges ',nyp1
   write(*,*)'get_grid_dims: num levels    gridcell   edges ',nzp1
endif

end subroutine get_grid_dims


!-----------------------------------------------------------------------
!>
!> Open and read the grid from a netCDF file.
!> All latitudes and longitudes in this module are stored in degrees.
!> All longitudes are [0,360)
!> All latitudes are [-90,90]
!> All levels are depths (-Inf,0]
!>
!> FIXME ... should add a fatal exit if domain encompasses the pole.
!>
!> FIXME Angie ... fill in the description of the module variables that get modified.
!> FIXME Angie ... fill in the description of the module variables that get modified.
!> FIXME Angie ... fill in the description of the module variables that get modified.

subroutine read_grid()

integer :: ncid

call nc_check(nf90_open(trim(gcom_geometry_file), nf90_nowrite, ncid), &
      'read_grid','open ['//trim(gcom_geometry_file)//']')

call get_grid_variable(ncid, 'ulon')
call get_grid_variable(ncid, 'ulat')
call get_grid_variable(ncid, 'ulev')

call get_grid_variable(ncid, 'vlon')
call get_grid_variable(ncid, 'vlat')
call get_grid_variable(ncid, 'vlev')

call get_grid_variable(ncid, 'wlon')
call get_grid_variable(ncid, 'wlat')
call get_grid_variable(ncid, 'wlev')

call get_grid_variable(ncid, 'lon')
call get_grid_variable(ncid, 'lat')
call get_grid_variable(ncid, 'lev')

call nc_check(nf90_close(ncid), &
         'read_grid','close '//trim(gcom_geometry_file) )

! FIXME ... OPTIONAL
! if the other grids can be calculated and cannot
!           be read from the geometry file, there must
!           be a subroutine call to derive the other grids.
! call calc_upoints(LAT, LON, ULAT, ULON)
! call calc_vpoints(LAT, LON, VLAT, VLON)
! call calc_wpoints(LAT, LON, WLAT, WLON)

! Ensure the loongitudes and latitudes are within expected range.

where (ULON <   0.0_r8) ULON = ULON + 360.0_r8
where (ULON > 360.0_r8) ULON = ULON - 360.0_r8
where (ULAT < -90.0_r8) ULAT = -90.0_r8
where (ULAT >  90.0_r8) ULAT =  90.0_r8

where (VLON <   0.0_r8) VLON = VLON + 360.0_r8
where (VLON > 360.0_r8) VLON = VLON - 360.0_r8
where (VLAT < -90.0_r8) VLAT = -90.0_r8
where (VLAT >  90.0_r8) VLAT =  90.0_r8

where (WLON <   0.0_r8) WLON = WLON + 360.0_r8
where (WLON > 360.0_r8) WLON = WLON - 360.0_r8
where (WLAT < -90.0_r8) WLAT = -90.0_r8
where (WLAT >  90.0_r8) WLAT =  90.0_r8

where ( LON <   0.0_r8)  LON =  LON + 360.0_r8
where ( LON > 360.0_r8)  LON =  LON - 360.0_r8
where ( LAT < -90.0_r8)  LAT = -90.0_r8
where ( LAT >  90.0_r8)  LAT =  90.0_r8

if (do_output() .and. (debug > 3)) then

   write(*,*)'Summary of the grid variables.'

   write(*,*)'range of all T,S,p longitudes ',minval(LON), maxval(LON)
   write(*,*)'range of all T,S,p latitude   ',minval(LAT), maxval(LAT)
   write(*,*)'range of all T,S,p levels     ',minval(LEV), maxval(LEV)

   write(*,*)'range of all U     longitudes ',minval(ULON), maxval(ULON)
   write(*,*)'range of all U     latitude   ',minval(ULAT), maxval(ULAT)
   write(*,*)'range of all U     levels     ',minval(ULEV), maxval(ULEV)

   write(*,*)'range of all V     longitudes ',minval(VLON), maxval(VLON)
   write(*,*)'range of all V     latitude   ',minval(VLAT), maxval(VLAT)
   write(*,*)'range of all V     levels     ',minval(VLEV), maxval(VLEV)

   write(*,*)'range of all W     longitudes ',minval(WLON), maxval(WLON)
   write(*,*)'range of all W     latitude   ',minval(WLAT), maxval(WLAT)
   write(*,*)'range of all W     levels     ',minval(WLEV), maxval(WLEV)

   write(*,*)

endif

end subroutine read_grid


!------------------------------------------------------------------
!>  This routine checks the user input against the variables available in the
!>  input netcdf file to see if it is possible to construct the DART state vector
!>  specified by the input.nml:model_nml:gcom_variables  variable.
!>  Each variable must have 5 entries.
!>  1: variable name
!>  2: DART KIND
!>  3: minimum value - as a character string - if none, use 'NA'
!>  4: maximum value - as a character string - if none, use 'NA'
!>  5: does the variable get updated in the restart file or not ...
!>     only variables from restart files may be updated.
!>     'UPDATE'       => update the variable in the restart file
!>     'NO_COPY_BACK' => do not copy the variable back to the restart file
!>     all these variables will be updated INTERNALLY IN DART
!>
!>  The calling code should check to see if the variable exists.

subroutine parse_variable_table( variable_list, ngood, table )

character(len=*), dimension(:),   intent(in)  :: variable_list
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer :: nrows, ncols, i
character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 6

nrows = size(table,1)
ncols = size(table,2)

! This loop just repackages the 1D array of values into a 2D array.
! We can do some miniminal checking along the way.
! Determining which file to check is going to be more complicated.

ngood = 0
MyLoop : do i = 1, nrows

   varname      = trim(variable_list(ncols*i - 4))
   dartstr      = trim(variable_list(ncols*i - 3))
   minvalstring = trim(variable_list(ncols*i - 2))
   maxvalstring = trim(variable_list(ncols*i - 1))
   state_or_aux = trim(variable_list(ncols*i    ))

   call to_upper(state_or_aux)

   table(i,VT_VARNAMEINDX) = trim(varname)
   table(i,VT_KINDINDX)    = trim(dartstr)
   table(i,VT_MINVALINDX)  = trim(minvalstring)
   table(i,VT_MAXVALINDX)  = trim(maxvalstring)
   table(i,VT_STATEINDX)   = trim(state_or_aux)

   ! If the first element is empty, we have found the end of the list.
   if ( table(i,1) == ' ' ) exit MyLoop

   ! Any other condition is an error.
   if ( any(table(i,:) == ' ') ) then
      string1 = 'input.nml &model_nml:gcom_variables not fully specified'
      string2 = 'must be 5 entries per variable. Last known variable name is'
      string3 = '['//trim(table(i,1))//'] ... (without the [], naturally)'
      call error_handler(E_ERR, 'parse_variable_table', string1, &
         source, revision, revdate, text2=string2, text3=string3)
   endif

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'parse_variable_table',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if (do_output() .and. (debug > 3)) then
      write(logfileunit,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2)),' ', &
                                               trim(table(i,3)), ' ', trim(table(i,4)),' ', &
                                               trim(table(i,5))
      write(     *     ,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2)),' ', &
                                               trim(table(i,3)), ' ', trim(table(i,4)),' ', &
                                               trim(table(i,5))
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'parse_variable_table',string1,text2=string2)
endif

end subroutine parse_variable_table


!------------------------------------------------------------------
!> Compute the offsets into the state vector for the start of each
!> variable type. Requires reading variable shapes from the restart file.
!> Record the extent of the data type in the state vector.

subroutine fill_progvar()

integer  :: ncid, io
integer  :: i, ivar, index1, indexN
integer  :: TimeDimID, VarID
integer  :: dimIDs(NF90_MAX_VAR_DIMS)
integer  :: varsize, dimlen
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: spvalR8, minvalue, maxvalue
logical  :: has_singleton_dimension
character(len=NF90_MAX_NAME) :: dimname

! Open the file to get dimensions, metadata.

call nc_check(nf90_open(trim(gcom_restart_file), NF90_NOWRITE, ncid), &
              'fill_progvar','open ['//trim(gcom_restart_file)//']')

! File is not required to have a time dimension
io = nf90_inq_dimid(ncid, 'time', TimeDimID)
if (io /= NF90_NOERR) TimeDimID = MISSING_I

! Loop over all the variables used to create the DART vector

index1  = 1;
indexN  = 0;
do ivar = 1, nfields

   ! Copying the information in the variable_table to each progvar.
   ! Setting the default values.

   progvar(ivar)%varname     = trim(variable_table(ivar,VT_VARNAMEINDX))
   progvar(ivar)%kind_string = trim(variable_table(ivar,VT_KINDINDX))
   progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
   progvar(ivar)%origin      = trim(gcom_restart_file)
   progvar(ivar)%numdims     = 0
   progvar(ivar)%rank        = 0
   progvar(ivar)%maxlevels   = 0
   progvar(ivar)%dimlens     = 0
   progvar(ivar)%dimnames    = ' '
   progvar(ivar)%coordinates = ' '
   progvar(ivar)%spvalINT    = -9999        ! should come from gcom
   progvar(ivar)%spvalR4     = 1.e36_r4     ! should come from gcom
   progvar(ivar)%spvalR8     = 1.e36_r8     ! should come from gcom
   progvar(ivar)%missingINT  = MISSING_I
   progvar(ivar)%missingR4   = MISSING_R4
   progvar(ivar)%missingR8   = MISSING_R8
   progvar(ivar)%rangeRestricted   = BOUNDED_NONE
   progvar(ivar)%minvalue          = MISSING_R8
   progvar(ivar)%maxvalue          = MISSING_R8
   progvar(ivar)%has_fill_value    = .false.
   progvar(ivar)%has_missing_value = .false.
   progvar(ivar)%update            = .false.

   if (variable_table(ivar,VT_STATEINDX) == 'UPDATE') progvar(ivar)%update = .true.

   ! If the character string can be interpreted as an r8, great.
   ! If not, there is no value to be used.

   read(variable_table(ivar,VT_MINVALINDX),*,iostat=io) minvalue
   if (io == 0) progvar(ivar)%minvalue = minvalue

   read(variable_table(ivar,VT_MAXVALINDX),*,iostat=io) maxvalue
   if (io == 0) progvar(ivar)%maxvalue = maxvalue

   ! rangeRestricted == BOUNDED_NONE  == 0 ... unlimited range
   ! rangeRestricted == BOUNDED_BELOW == 1 ... minimum, but no maximum
   ! rangeRestricted == BOUNDED_ABOVE == 2 ... maximum, but no minimum
   ! rangeRestricted == BOUNDED_BOTH  == 3 ... minimum and maximum

   if (   (progvar(ivar)%minvalue       /= MISSING_R8) .and. &
          (progvar(ivar)%maxvalue       /= MISSING_R8) ) then
           progvar(ivar)%rangeRestricted = BOUNDED_BOTH
   elseif (progvar(ivar)%maxvalue       /= MISSING_R8) then
           progvar(ivar)%rangeRestricted = BOUNDED_ABOVE
   elseif (progvar(ivar)%minvalue       /= MISSING_R8) then
           progvar(ivar)%rangeRestricted = BOUNDED_BELOW
   else
           progvar(ivar)%rangeRestricted = BOUNDED_NONE
   endif

   ! Check to make sure min is less than max if both are specified.

   if ( progvar(ivar)%rangeRestricted == BOUNDED_BOTH ) then
      if (maxvalue < minvalue) then
         write(string1,*)'&model_nml state_variable input error for ', &
                          trim(progvar(ivar)%varname)
         write(string2,*)'minimum value (',minvalue,') must be less than '
         write(string3,*)'maximum value (',maxvalue,')'
         call error_handler(E_ERR,'fill_progvar',string1, &
            source,revision,revdate,text2=string2,text3=string3)
      endif
   endif

   string2 = trim(gcom_restart_file)//' '//trim(progvar(ivar)%varname)

   call nc_check(nf90_inq_varid(ncid, trim(progvar(ivar)%varname), VarID), &
            'fill_progvar', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, &
                 ndims=progvar(ivar)%numdims, xtype=progvar(ivar)%xtype), &
            'fill_progvar', 'inquire '//trim(string2))

   ! If the long_name and/or units attributes are set, get them.
   ! They are not REQUIRED to exist but are nice to use if they are present.

   if( nf90_inquire_attribute(    ncid, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(ncid, VarID, 'long_name' , progvar(ivar)%long_name), &
                  'fill_progvar', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = progvar(ivar)%varname
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(ncid, VarID, 'units' , progvar(ivar)%units), &
                  'fill_progvar', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = 'unknown'
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'coordinates') == NF90_NOERR )  then
      call nc_check( nf90_get_att(ncid, VarID, 'coordinates' , progvar(ivar)%coordinates), &
                  'fill_progvar', 'get_att coordinates '//trim(string2))
   else
      progvar(ivar)%coordinates = 'unknown'
   endif

   ! Saving any FillValue, missing_value attributes to be used later.

   if (progvar(ivar)%xtype == NF90_INT) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalINT) == NF90_NOERR) then
          progvar(ivar)%spvalINT       = spvalINT
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalINT) == NF90_NOERR) then
          progvar(ivar)%missingINT        = spvalINT
          progvar(ivar)%has_missing_value = .true.
       endif

   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalR4) == NF90_NOERR) then
          progvar(ivar)%spvalR4        = spvalR4
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR4) == NF90_NOERR) then
          progvar(ivar)%missingR4         = spvalR4
          progvar(ivar)%has_missing_value = .true.
       endif

   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalR8) == NF90_NOERR) then
          progvar(ivar)%spvalR8        = spvalR8
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR8) == NF90_NOERR) then
          progvar(ivar)%missingR8         = spvalR8
          progvar(ivar)%has_missing_value = .true.
       endif
   endif

   ! This block captures the natural shape of the variable

   varsize = 1
   dimlen  = 1
   has_singleton_dimension = .FALSE.

   DimensionLoop : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimlen), &
                                          'fill_progvar', string1)

      ! Only reserve space for a single time slice
      if (dimIDs(i) == TimeDimID) then
         has_singleton_dimension = .TRUE.
         dimlen = 1
      endif

      progvar(ivar)%dimlens( i) = dimlen
      progvar(ivar)%dimnames(i) = dimname
      varsize = varsize * dimlen

   enddo DimensionLoop

   ! Record the portion of the DART vector reserved for this variable.
   if ( has_singleton_dimension ) then
      progvar(ivar)%rank = progvar(ivar)%numdims - 1
   else
      progvar(ivar)%rank = progvar(ivar)%numdims
   endif
   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

enddo

call nc_check(nf90_close(ncid),'fill_progvar','close '//trim(string2))
ncid = 0

model_size = progvar(nfields)%indexN

if (do_output() .and. (debug > 3)) call progvar_summary()

end subroutine fill_progvar


!-----------------------------------------------------------------------
!>
!> public utility routine to return the grid sizes.

subroutine progvar_summary()

integer :: i, ivar

do ivar = 1,nfields

   write(logfileunit,*)
   write(logfileunit,*) 'variable number',ivar, 'is ['//trim(progvar(ivar)%varname)//']'
   write(logfileunit,*) '  filename    ',trim(progvar(ivar)%origin)
   write(logfileunit,*) '  update      ',progvar(ivar)%update
   write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
   write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
   write(logfileunit,*) '  xtype       ',progvar(ivar)%xtype
   write(logfileunit,*) '  numdims     ',progvar(ivar)%numdims
   write(logfileunit,*) '  rank        ',progvar(ivar)%rank
   write(logfileunit,*) '  dimlens     ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
   write(logfileunit,*) '  dimnames    ',( trim(progvar(ivar)%dimnames(i))//' ', i=1,progvar(ivar)%numdims )
   write(logfileunit,*) '  coordinates ',trim(progvar(ivar)%coordinates)
   write(logfileunit,*) '  varsize     ',progvar(ivar)%varsize
   write(logfileunit,*) '  index1      ',progvar(ivar)%index1
   write(logfileunit,*) '  indexN      ',progvar(ivar)%indexN
   write(logfileunit,*) '  dart_kind   ',progvar(ivar)%dart_kind
   write(logfileunit,*) '  kind_string ',progvar(ivar)%kind_string
   write(logfileunit,*) '  spvalINT    ',progvar(ivar)%spvalINT
   write(logfileunit,*) '  spvalR4     ',progvar(ivar)%spvalR4
   write(logfileunit,*) '  spvalR8     ',progvar(ivar)%spvalR8
   write(logfileunit,*) '  missingINT  ',progvar(ivar)%missingINT
   write(logfileunit,*) '  missingR4   ',progvar(ivar)%missingR4
   write(logfileunit,*) '  missingR8   ',progvar(ivar)%missingR8
   write(logfileunit,*) '  has_fill_value    ',progvar(ivar)%has_fill_value
   write(logfileunit,*) '  has_missing_value ',progvar(ivar)%has_missing_value
   write(logfileunit,*) '  rangeRestricted   ',progvar(ivar)%rangeRestricted
   write(logfileunit,*) '  minvalue          ',progvar(ivar)%minvalue
   write(logfileunit,*) '  maxvalue          ',progvar(ivar)%maxvalue
   write(logfileunit,*)

   write(     *     ,*)
   write(     *     ,*) 'variable number',ivar, 'is ['//trim(progvar(ivar)%varname)//']'
   write(     *     ,*) '  filename    ',trim(progvar(ivar)%origin)
   write(     *     ,*) '  update      ',progvar(ivar)%update
   write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
   write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
   write(     *     ,*) '  xtype       ',progvar(ivar)%xtype
   write(     *     ,*) '  numdims     ',progvar(ivar)%numdims
   write(     *     ,*) '  rank        ',progvar(ivar)%rank
   write(     *     ,*) '  dimlens     ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
   write(     *     ,*) '  dimnames    ',( trim(progvar(ivar)%dimnames(i))//' ', i=1,progvar(ivar)%numdims )
   write(     *     ,*) '  coordinates ',trim(progvar(ivar)%coordinates)
   write(     *     ,*) '  varsize     ',progvar(ivar)%varsize
   write(     *     ,*) '  index1      ',progvar(ivar)%index1
   write(     *     ,*) '  indexN      ',progvar(ivar)%indexN
   write(     *     ,*) '  dart_kind   ',progvar(ivar)%dart_kind
   write(     *     ,*) '  kind_string ',progvar(ivar)%kind_string
   write(     *     ,*) '  spvalINT    ',progvar(ivar)%spvalINT
   write(     *     ,*) '  spvalR4     ',progvar(ivar)%spvalR4
   write(     *     ,*) '  spvalR8     ',progvar(ivar)%spvalR8
   write(     *     ,*) '  missingINT  ',progvar(ivar)%missingINT
   write(     *     ,*) '  missingR4   ',progvar(ivar)%missingR4
   write(     *     ,*) '  missingR8   ',progvar(ivar)%missingR8
   write(     *     ,*) '  has_fill_value    ',progvar(ivar)%has_fill_value
   write(     *     ,*) '  has_missing_value ',progvar(ivar)%has_missing_value
   write(     *     ,*) '  rangeRestricted   ',progvar(ivar)%rangeRestricted
   write(     *     ,*) '  minvalue          ',progvar(ivar)%minvalue
   write(     *     ,*) '  maxvalue          ',progvar(ivar)%maxvalue
   write(     *     ,*)

enddo

write(string1,*) 'total model size is ',model_size
call error_handler(E_MSG,'progvar_summary',string1)

end subroutine progvar_summary


!-----------------------------------------------------------------------
!>
!> public utility routine to make the model_nml:gcom_restart_file accessible

subroutine get_gcom_restart_filename(filename)

character(len=*), intent(out) :: filename

if ( .not. module_initialized ) call static_init_model

filename = gcom_restart_file

end subroutine get_gcom_restart_filename


!-----------------------------------------------------------------------
!>  define_var_dims() takes the N-dimensional variable and appends the DART
!>  dimensions of 'copy' and 'time'. If the variable initially had a 'time'
!>  dimension, it is ignored because (by construction) it is a singleton
!>  dimension.

subroutine define_var_dims(ivar, ncid, memberdimid, unlimiteddimid, ndims, dimids)

integer,               intent(in)  :: ivar
integer,               intent(in)  :: ncid
integer,               intent(in)  :: memberdimid
integer,               intent(in)  :: unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames

integer :: i, mydimid

ndims = 0

DIMLOOP : do i = 1,progvar(ivar)%numdims

   if (progvar(ivar)%dimnames(i) == 'time') cycle DIMLOOP

   call nc_check(nf90_inq_dimid(ncid=ncid, name=progvar(ivar)%dimnames(i), dimid=mydimid), &
           'define_var_dims','inq_dimid ['//trim(progvar(ivar)%dimnames(i))//']')

   ndims         = ndims + 1
   dimids(ndims) = mydimid
   dimnames(ndims) = progvar(ivar)%dimnames(i)

enddo DIMLOOP

! The last two dimensions are always 'copy' and 'time'
ndims           = ndims + 1
dimids(ndims)   = memberdimid
dimnames(ndims) = 'copy'
ndims           = ndims + 1
dimids(ndims)   = unlimitedDimid
dimnames(ndims) = 'time'

if (do_output() .and. (debug > 99)) then

   write(logfileunit,*)
   write(logfileunit,*)'define_var_dims knowledge'

   write(logfileunit,*)trim(progvar(ivar)%varname),' has original     dimnames: ', &
                      (trim(progvar(ivar)%dimnames(i))//' ',i=1,progvar(ivar)%numdims)
   write(logfileunit,*)trim(progvar(ivar)%varname),' repackaging into dimnames: ', &
                      (trim(dimnames(i))//' ',i=1,ndims)
   write(logfileunit,*)'thus dimids ',dimids(1:ndims)

   write(     *     ,*)
   write(     *     ,*)'define_var_dims knowledge'
   write(     *     ,*)trim(progvar(ivar)%varname),' has original     dimnames: ', &
                      (trim(progvar(ivar)%dimnames(i))//' ',i=1,progvar(ivar)%numdims)
   write(     *     ,*)trim(progvar(ivar)%varname),' repackaging into dimnames: ', &
                      (trim(dimnames(i))//' ',i=1,ndims)
   write(     *     ,*)'thus dimids ',dimids(1:ndims)

endif

return
end subroutine define_var_dims


!-----------------------------------------------------------------------
!>
!> Reads the current time and state variables from a GCOM restart
!> file and packs them into a dart state vector.

subroutine gcom_file_to_dart_vector(filename, state_vector, state_time)

character(len=*), intent(in)    :: filename
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: state_time

integer :: i, j, k, ivar, indx, iunit

! temp space to hold data while we are reading it
real(r8), allocatable :: mytimes(:)
real(r8), allocatable :: data_1d_array(:)
real(r8), allocatable :: data_2d_array(:,:)
real(r8), allocatable :: data_3d_array(:,:,:)

integer :: dimIDs(NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: attvalue

integer :: ncid, TimeDimID, VarID, numdims, dimlen, ntimes
integer :: iday, isecond, origin_days, origin_seconds
type(time_type) :: origin_time

logical :: write_metadata
real(r8) :: mylon, mylat, mylev

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

if ( .not. file_exist(filename) ) then
   write(string1,'(A)') '['//trim(filename)//'] does not exist.'
   call error_handler(E_ERR,'gcom_file_to_dart_vector',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
        'gcom_file_to_dart_vector', 'open '//trim(filename))

! Read the time data.

call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
        'gcom_file_to_dart_vector', 'inq_dimid time')

call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=ntimes), &
        'gcom_file_to_dart_vector', 'inquire_dimension time')

allocate(mytimes(ntimes))

call nc_check(nf90_inq_varid(ncid, 'time', VarID), &
        'gcom_file_to_dart_vector', 'inq_varid time')

call nc_check(nf90_get_var(  ncid, VarID, mytimes), &
        'gcom_file_to_dart_vector', 'get_var   time')

call nc_check(nf90_get_att(ncid, VarID, 'units', attvalue), &
        'gcom_file_to_dart_vector', 'time get_att units')

origin_time = convert_times(attvalue, mytimes)

call get_time(origin_time, origin_seconds, origin_days)

! FIXME ... assuming we are using the last time ...
! The module 'model_time' is set and then the output argument
! is set to the same value.

iday       = int(mytimes(ntimes))
isecond    = (mytimes(ntimes) - iday)*86400
model_time = set_time(origin_seconds+isecond, origin_days+iday)
state_time = model_time

if (do_output()) then
    call print_time(model_time,'time for restart file '//trim(filename))
    call print_date(model_time,'date for restart file '//trim(filename))
endif

! If desired, output a file that has the full i,j,k,lat,lon,lev for every element of the DART vector
if (do_output() .and. (debug > 7)) then
   iunit = open_file('metadata_table.txt',form='formatted',action='write')
   write(*,*)'varname   indx          i          j          k'
   write_metadata = .true.
else
   write_metadata = .false.
endif

 100 format(A,2(1x,i10))              ! varname, DART index, i
 200 format(A,3(1x,i10))              ! varname, DART index, i, j
 300 format(A,4(1x,i10),3(1x,f18.10)) ! varname, DART index, i, j, k, and the 3 locations

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.

do ivar=1, nfields

   varname = progvar(ivar)%varname
   string1 = trim(filename)//' '//trim(varname)

   ! Is the netCDF variable the right shape?

   call nc_check(nf90_inq_varid(ncid, varname, VarID), &
            'gcom_file_to_dart_vector', 'inq_varid '//trim(string1))

   call nc_check(nf90_inquire_variable(ncid,VarId,dimids=dimIDs,ndims=numdims), &
            'gcom_file_to_dart_vector', 'inquire '//trim(string1))

   if (numdims /= progvar(ivar)%numdims) then
      write(string2,*) trim(string1),' has rank ', numdims
      write(string3,*) 'expected rank ', progvar(ivar)%numdims, ' ... stopping.'
      call error_handler(E_ERR,'gcom_file_to_dart_vector',string1, &
                 source,revision,revdate,text2=string2,text3=string3)
   endif

   do i = 1,numdims
      write(string1,'(''inquire dimension'',i2,A)') i,trim(string1)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'gcom_file_to_dart_vector', string1)

      ! DART only keeps 1 timestep around
      if (dimIDs(i) == TimeDimID) dimlen = 1

      if (dimlen /= progvar(ivar)%dimlens(i)) then
         write(string2,*) 'dim/dimlen',i,dimlen, &
                          'not expected value of ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'gcom_file_to_dart_vector',string1, &
                    source,revision,revdate,text2=string2)
      endif
   enddo

   indx = progvar(ivar)%index1

   ! The singleton dimension (time) is always the last dimension,
   ! so we check the rank of the variable not the netcdf number of dimensions.
   ! Some of the variables may or may not have a 'time' dimension,
   ! I haven't seen them yet, so I don't know.

   if (progvar(ivar)%rank == 1) then

      allocate(data_1d_array(progvar(ivar)%dimlens(1)))
      call DART_get_var(ncid, varname, data_1d_array)

      if (    progvar(ivar)%units == 'm/s') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'psu') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'degrees celsius') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'degrees C') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'bar') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'km/s') then
         write(string1,*) 'converting variable ', trim(progvar(ivar)%varname)
         write(string2,*) 'from ', trim(progvar(ivar)%units),' to m/s'
         call error_handler(E_MSG,'gcom_file_to_dart_vector', string1, text2=string2)

         data_1d_array = data_1d_array/1000.0_r8 ! convert km/s to meters/s
      else
         write(string1,*) 'no support for units of ', trim(progvar(ivar)%units)
         write(string2,*) 'on variable ', trim(progvar(ivar)%varname)
         call error_handler(E_ERR,'gcom_file_to_dart_vector', string1, &
                        source, revision, revdate, text2=string2)
      endif

      do i = 1, progvar(ivar)%dimlens(1)
         state_vector(indx) = data_1d_array(i)
         if (write_metadata) write(iunit,100) trim(progvar(ivar)%varname), indx, i
         indx = indx + 1
      enddo
      deallocate(data_1d_array)

   elseif (progvar(ivar)%rank == 2) then

      allocate(data_2d_array(progvar(ivar)%dimlens(1), &
                             progvar(ivar)%dimlens(2)))
      call DART_get_var(ncid, varname, data_2d_array)

      if (    progvar(ivar)%units == 'm/s') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'psu') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'degrees celsius') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'degrees C') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'bar') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'km/s') then
         write(string1,*) 'converting variable ', trim(progvar(ivar)%varname)
         write(string2,*) 'from ', trim(progvar(ivar)%units),' to m/s'
         call error_handler(E_MSG,'gcom_file_to_dart_vector', string1, text2=string2)

         data_2d_array = data_2d_array/1000.0_r8 ! convert km/s to meters/s
      else
         write(string1,*) 'no support for units of ', trim(progvar(ivar)%units)
         write(string2,*) 'on variable ', trim(progvar(ivar)%varname)
         call error_handler(E_ERR,'gcom_file_to_dart_vector', string1, &
                        source, revision, revdate, text2=string2)
      endif

      do j = 1, progvar(ivar)%dimlens(2)
      do i = 1, progvar(ivar)%dimlens(1)
         state_vector(indx) = data_2d_array(i, j)
         if (write_metadata) write(iunit,200) trim(progvar(ivar)%varname), indx, i, j
         indx = indx + 1
      enddo
      enddo
      deallocate(data_2d_array)

   elseif (progvar(ivar)%rank == 3) then

      allocate(data_3d_array(progvar(ivar)%dimlens(1), &
                             progvar(ivar)%dimlens(2), &
                             progvar(ivar)%dimlens(3)))
      call DART_get_var(ncid, varname, data_3d_array)

      if (    progvar(ivar)%units == 'm/s') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'psu') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'degrees celsius') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'degrees C') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'bar') then
         ! Units are correct, nothing to do.
      elseif (progvar(ivar)%units == 'km/s') then
         write(string1,*) 'converting variable ', trim(progvar(ivar)%varname)
         write(string2,*) 'from ', trim(progvar(ivar)%units),' to m/s'
         call error_handler(E_MSG,'gcom_file_to_dart_vector', string1, text2=string2)

         data_3d_array = data_3d_array/1000.0_r8 ! convert km/s to meters/s
      else
         write(string1,*) 'no support for units of ', trim(progvar(ivar)%units)
         write(string2,*) 'on variable ', trim(progvar(ivar)%varname)
         call error_handler(E_ERR,'gcom_file_to_dart_vector', string1, &
                        source, revision, revdate, text2=string2)
      endif

      do k = 1, progvar(ivar)%dimlens(3)
      do j = 1, progvar(ivar)%dimlens(2)
      do i = 1, progvar(ivar)%dimlens(1)
         state_vector(indx) = data_3d_array(i, j, k)
         indx = indx + 1
      enddo
      enddo
      enddo

      deallocate(data_3d_array)

      ! This is a debug block. The ugly part is to figure out which set of
      ! lats/lons go with each variable.
      if (write_metadata) then

         indx = progvar(ivar)%index1
         do k = 1, progvar(ivar)%dimlens(3)
         do j = 1, progvar(ivar)%dimlens(2)
         do i = 1, progvar(ivar)%dimlens(1)

            call get_state_lonlatlev(ivar,i,j,k,mylon,mylat,mylev)

            ! TJH DEBUG - if you want to write [-180,180]
            ! if (mylon > 180.0_r8) mylon = mylon - 360.0_r8

            write(iunit,300) trim(progvar(ivar)%varname), indx, i, j, k, mylon, mylat, mylev
            indx = indx + 1
         enddo
         enddo
         enddo
      endif

   else

      write(string1,*) 'no support for data array of dimension ', numdims
      write(string2,*) 'variable [',trim(progvar(ivar)%varname),']'
      write(string3,*) 'file [',trim(progvar(ivar)%origin),']'
      call error_handler(E_ERR,'gcom_file_to_dart_vector', string1, &
                        source, revision, revdate, text2=string2, text3=string3)
   endif

enddo

if (write_metadata) call close_file(iunit)

if (do_output() .and. (debug > 0)) call print_ranges(state_vector)

end subroutine gcom_file_to_dart_vector


!-----------------------------------------------------------------------
!>
!> Writes the current time and state variables from a dart state
!> vector (1d fortran array) into a GCOM netcdf restart file.
!> FIXME The more I think about it, the more that 'statedate'
!> is always the same as 'model_time' ... must check
!> Since the progvar(ivar)%dims is based on the GCOM file, all
!> the shape error checking may be superfluous.

subroutine dart_vector_to_gcom_file(state_vector, filename, statedate)

real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: statedate

! temp space to hold data while we are writing it
real(r8), allocatable :: data_1d_array(:)
real(r8), allocatable :: data_2d_array(:,:)
real(r8), allocatable :: data_3d_array(:,:,:)

integer  ::  dimIDs(NF90_MAX_VAR_DIMS)
integer  :: ncstart(NF90_MAX_VAR_DIMS)
integer  :: nccount(NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: varname

integer :: i, ivar, ncid, TimeDimID, VarID, dimlen, numdims
integer :: timeindex

!----------------------------------------------------------------------
! Get the show underway
!----------------------------------------------------------------------

if ( .not. module_initialized ) call static_init_model

! Check that the input file exists.

if ( .not. file_exist(filename)) then
   write(string1,*)trim(filename),' does not exist. FATAL error.'
   call error_handler(E_ERR,'dart_vector_to_gcom_file',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_WRITE, ncid), &
        'dart_vector_to_gcom_file', 'open '//trim(filename))

! make sure the time tag in the restart file matches
! the current time of the DART state ...
! find_desired_time_index will abort if it cannot find a match.

call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
        'dart_vector_to_gcom_file', 'inq_dimid time')
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=dimlen ), &
        'dart_vector_to_gcom_file', 'inquire_dimension time '//trim(varname))

timeindex = find_desired_time_index(ncid, dimlen, statedate)

if (do_output()) then
    call print_time(statedate,'time of restart file '//trim(filename))
    call print_date(statedate,'date of restart file '//trim(filename))
endif

UPDATE : do ivar=1,nfields

   varname = trim(progvar(ivar)%varname)
   string2 = trim(filename)//' '//trim(varname)

   if ( .not. progvar(ivar)%update ) then
      write(string1,*)'intentionally not updating '//trim(string2)
      write(string3,*)'as per namelist control in model_nml:gcom_variables'
      call error_handler(E_MSG, 'dart_vector_to_gcom_file', string1, text2=string3)
      cycle UPDATE
   endif

   ! Ensure netCDF variable is conformable with progvar quantity.

   call nc_check(nf90_inq_varid(ncid, varname, VarID), &
           'dart_vector_to_gcom_file','inq_varid '//trim(string2))
   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=numdims), &
           'dart_vector_to_gcom_file','inquire_variable '//trim(string2))

   if (numdims /= progvar(ivar)%numdims) then
      write(string1,*) trim(filename)//' '//trim(varname)
      write(string2,*) 'has rank ',numdims, ' /= expected value of ', &
                       progvar(ivar)%numdims
      call error_handler(E_ERR,'dart_vector_to_gcom_file',string1, &
                 source,revision,revdate,text2=string2)
   endif

   ncstart = 1
   nccount = 1

   DimCheck : do i = 1,numdims

      write(string1,'(''inquire dimension'',1x,i2,1x,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'dart_vector_to_gcom_file', string1)

      if (dimIDs(i) == TimeDimID) then
         ncstart(i) = timeindex
         nccount(i) = 1
         dimlen     = 1
      else
         ncstart(i) = 1
         nccount(i) = dimlen
      endif

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' dim/dimlen ',i,dimlen
         write(string2,*) ' not ',progvar(ivar)%dimlens(i),' as it should be.'
         call error_handler(E_ERR, 'dart_vector_to_gcom_file', string1, &
                         source, revision, revdate, text2=string2)
      endif

   enddo DimCheck

   if (do_output() .and. (debug > 99)) then
      write(*,*)'dart_vector_to_gcom_file: variable ['//trim(varname)//']'
      write(*,*)'dart_vector_to_gcom_file: start ',ncstart(1:numdims)
      write(*,*)'dart_vector_to_gcom_file: count ',nccount(1:numdims)
   endif

   ! When called with a 4th argument, vector_to_variable() replaces the DART
   ! missing code with the value in the corresponding variable in the netCDF file.
   ! Any clamping to physically meaningful values occurrs in vector_to_variable.

   if (progvar(ivar)%rank == 1) then

      allocate(data_1d_array(progvar(ivar)%dimlens(1)))
      call vector_to_variable(state_vector, ivar, data_1d_array, ncid)

      call nc_check(nf90_put_var(ncid, VarID, data_1d_array, &
              start=ncstart(1:numdims), &
              count=nccount(1:numdims)), &
              'dart_vector_to_gcom_file', 'put_var '//trim(varname))
      deallocate(data_1d_array)

   elseif (progvar(ivar)%rank == 2) then

      allocate(data_2d_array(progvar(ivar)%dimlens(1), &
                             progvar(ivar)%dimlens(2)))
      call vector_to_variable(state_vector, ivar, data_2d_array, ncid)

      call nc_check(nf90_put_var(ncid, VarID, data_2d_array, &
              start=ncstart(1:numdims), &
              count=nccount(1:numdims)), &
              'dart_vector_to_gcom_file', 'put_var '//trim(varname))
      deallocate(data_2d_array)

   elseif (progvar(ivar)%rank == 3) then

      allocate(data_3d_array(progvar(ivar)%dimlens(1), &
                             progvar(ivar)%dimlens(2), &
                             progvar(ivar)%dimlens(3)))
      call vector_to_variable(state_vector, ivar, data_3d_array, ncid)

      call nc_check(nf90_put_var(ncid, VarID, data_3d_array, &
              start=ncstart(1:numdims), &
              count=nccount(1:numdims)), &
              'dart_vector_to_gcom_file', 'put_var '//trim(varname))
      deallocate(data_3d_array)

   else
      write(string1, *) 'no support for data array of rank ', dimlen
      call error_handler(E_ERR,'dart_vector_to_gcom_file', string1, &
                        source,revision,revdate)
   endif

   ! TJH FIXME ... this works perfectly if it were not for a bug in netCDF.
   ! When they fix the bug, this will be a useful thing to restore.
   ! Make note that the variable has been updated by DART
!  call nc_check(nf90_Redef(ncid),'dart_vector_to_gcom_file', 'redef '//trim(filename))
!  call nc_check(nf90_put_att(ncid, VarID,'DART','variable modified by DART'),&
!                'dart_vector_to_gcom_file', 'modified '//trim(varname))
!  call nc_check(nf90_enddef(ncid),'dart_vector_to_gcom_file','state enddef '//trim(filename))

enddo UPDATE

call nc_check(nf90_close(ncid),'dart_vector_to_gcom_file','close '//trim(filename))

end subroutine dart_vector_to_gcom_file


!-----------------------------------------------------------------------
! overloaded routines defining DART_get_var
!     MODULE PROCEDURE get_var_1d
!     MODULE PROCEDURE get_var_2d
!     MODULE PROCEDURE get_var_3d
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!> This function will return a R8 array with the netCDF attributes applied.
!> scale_factor, offset will be applied,
!> missing_value, _FillValue will be replaced by the DART missing value ...
!>
!> If _FillValue is defined then it should be scalar and of the same type as the
!> variable. If the variable is packed using scale_factor and add_offset attributes
!> (see below), the _FillValue attribute should have the data type of the packed data.
!>
!> missing_value
!> When scale_factor and add_offset are used for packing, the value(s) of the
!> missing_value attribute should be specified in the domain of the data in the
!> file (the packed data), so that missing values can be detected before the
!> scale_factor and add_offset are applied.
!>
!> scale_factor
!> If present for a variable, the data are to be multiplied by this factor after the
!> data are read by the application that accesses the data.  If valid values are
!> specified using the valid_min, valid_max, valid_range, or _FillValue attributes,
!> those values should be specified in the domain of the data in the file
!> (the packed data), so that they can be interpreted before the scale_factor and
!> add_offset are applied.
!>
!> add_offset
!> If present for a variable, this number is to be added to the data after it is read
!> by the application that accesses the data. If both scale_factor and add_offset
!> attributes are present, the data are first scaled before the offset is added.

subroutine get_var_1d(ncid, varname, var1d)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r8),         intent(out) :: var1d(:)

integer  ::  dimIDs(NF90_MAX_VAR_DIMS)
integer  :: dimlens(NF90_MAX_VAR_DIMS)
integer  :: ncstart(NF90_MAX_VAR_DIMS)
integer  :: nccount(NF90_MAX_VAR_DIMS)
integer  :: VarID, numdims, xtype, io1, io2
integer  :: TimeDimID, time_dimlen, timeindex
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable :: intarray(:)
real(r4), allocatable ::  r4array(:)

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. (debug > 3)) then
   write(*,*)
   write(logfileunit,*)
endif

! Need to know the Time Dimension ID and length

call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
        'get_var_1d', 'inq_dimid time '//trim(varname))
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
        'get_var_1d', 'inquire_dimension time '//trim(varname))

call nc_check(nf90_inq_varid(ncid, varname, VarID), &
                            'get_var_1d', 'inq_varid '//trim(varname))
call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_1d', 'inquire_variable '//trim(varname))
call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlens(1)), &
                            'get_var_1d', 'inquire_dimension '//trim(varname))

if ((numdims /= 1) .or. (size(var1d) /= dimlens(1)) ) then
   write(string1,*) trim(varname)//' is not the expected shape/length of ', size(var1d)
   call error_handler(E_ERR,'get_var_1d',string1,source,revision,revdate)
endif

ncstart = 1
nccount = dimlens(1)

if (dimIDs(1) == TimeDimID) then
   timeindex  = find_desired_time_index(ncid, time_dimlen)
   ncstart(1) = timeindex
   nccount(1) = 1
endif

if (do_output() .and. (debug > 3)) then
   write(*,*)'get_var_1d: variable ['//trim(varname)//']'
   write(*,*)'get_var_1d: start ',ncstart(1:numdims)
   write(*,*)'get_var_1d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//trim(varname))
   var1d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1)))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//trim(varname))
   var1d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var1d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//trim(varname))

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_1d',string1,source,revision,revdate)
endif

end subroutine get_var_1d


!-----------------------------------------------------------------------
!
!> This function will return a R8 array with the netCDF attributes applied.
!> scale_factor, offset will be applied,
!> missing_value, _FillValue will be replaced by the DART missing value ...
!>
!> If _FillValue is defined then it should be scalar and of the same type as the
!> variable. If the variable is packed using scale_factor and add_offset attributes
!> (see below), the _FillValue attribute should have the data type of the packed data.
!>
!> missing_value
!> When scale_factor and add_offset are used for packing, the value(s) of the
!> missing_value attribute should be specified in the domain of the data in the
!> file (the packed data), so that missing values can be detected before the
!> scale_factor and add_offset are applied.
!>
!> scale_factor
!> If present for a variable, the data are to be multiplied by this factor after the
!> data are read by the application that accesses the data.  If valid values are
!> specified using the valid_min, valid_max, valid_range, or _FillValue attributes,
!> those values should be specified in the domain of the data in the file
!> (the packed data), so that they can be interpreted before the scale_factor and
!> add_offset are applied.
!>
!> add_offset
!> If present for a variable, this number is to be added to the data after it is read
!> by the application that accesses the data. If both scale_factor and add_offset
!> attributes are present, the data are first scaled before the offset is added.

subroutine get_var_2d(ncid, varname, var2d)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r8),         intent(out) :: var2d(:,:)

integer  ::  dimIDs(NF90_MAX_VAR_DIMS)
integer  :: dimlens(NF90_MAX_VAR_DIMS)
integer  :: ncstart(NF90_MAX_VAR_DIMS)
integer  :: nccount(NF90_MAX_VAR_DIMS)
integer  :: VarID, numdims, xtype, io1, io2, i
integer  :: TimeDimID, time_dimlen, timeindex
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable :: intarray(:,:)
real(r4), allocatable ::  r4array(:,:)

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. (debug > 3)) then
   write(*,*)
   write(logfileunit,*)
endif

! Need to know the Time Dimension ID and length

call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
        'get_var_2d', 'inq_dimid time '//trim(varname))
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen), &
        'get_var_2d', 'inquire_dimension time '//trim(varname))

call nc_check(nf90_inq_varid(ncid, varname, VarID), 'get_var_2d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_2d', 'inquire_variable '//trim(varname))

ncstart(:) = 1
nccount(:) = 1

DimCheck : do i = 1,numdims

   write(string1,'(a,i2,1x,A)') 'inquire dimension ',i,trim(varname)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
           'get_var_2d', string1)

   if ( dimIDs(i) == TimeDimID ) then

      timeindex  = find_desired_time_index(ncid, time_dimlen)
      ncstart(i) = timeindex
      dimlens(i) = 1

   elseif ( size(var2d,i) /= dimlens(i) ) then
      write(string1,*) trim(varname)//' has shape ', dimlens(1:numdims)
      write(string2,*) 'which is not the expected shape of ', size(var2d,1),size(var2d,2)
      call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate,text2=string2)
   endif

   nccount(i) = dimlens(i)

enddo DimCheck

if (do_output() .and. (debug > 3)) then
   write(*,*)'get_var_2d: variable ['//trim(varname)//']'
   write(*,*)'get_var_2d: start ',ncstart(1:numdims)
   write(*,*)'get_var_2d: count ',nccount(1:numdims)
endif

! if ( (numdims /= 2)  ) then
!    write(string1,*) trim(varname)//' is not a 2D variable as expected.'
!    call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate)
! endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//trim(varname))
   var2d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1),dimlens(2)))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//trim(varname))
   var2d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var2d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//trim(varname))

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate)
endif

end subroutine get_var_2d


!-----------------------------------------------------------------------
!
!> This function will return a 3D R8 array with the netCDF attributes applied.
!> scale_factor, offset will be applied,
!> missing_value, _FillValue will be replaced by the DART missing value ...
!>
!> If _FillValue is defined then it should be scalar and of the same type as the
!> variable. If the variable is packed using scale_factor and add_offset attributes
!> (see below), the _FillValue attribute should have the data type of the packed data.
!>
!> missing_value
!> When scale_factor and add_offset are used for packing, the value(s) of the
!> missing_value attribute should be specified in the domain of the data in the
!> file (the packed data), so that missing values can be detected before the
!> scale_factor and add_offset are applied.
!>
!> scale_factor
!> If present for a variable, the data are to be multiplied by this factor after the
!> data are read by the application that accesses the data.  If valid values are
!> specified using the valid_min, valid_max, valid_range, or _FillValue attributes,
!> those values should be specified in the domain of the data in the file
!> (the packed data), so that they can be interpreted before the scale_factor and
!> add_offset are applied.
!>
!> add_offset
!> If present for a variable, this number is to be added to the data after it is read
!> by the application that accesses the data. If both scale_factor and add_offset
!> attributes are present, the data are first scaled before the offset is added.

subroutine get_var_3d(ncid, varname, var3d)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r8),         intent(out) :: var3d(:,:,:)

integer  ::  dimIDs(NF90_MAX_VAR_DIMS)
integer  :: dimlens(NF90_MAX_VAR_DIMS)
integer  :: ncstart(NF90_MAX_VAR_DIMS)
integer  :: nccount(NF90_MAX_VAR_DIMS)
integer  :: i, TimeDimID, time_dimlen, timeindex
integer  :: VarID, numdims, xtype, io1, io2
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable :: intarray(:,:,:)
real(r4), allocatable ::  r4array(:,:,:)

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. (debug > 3)) then
   write(*,*)
   write(logfileunit,*)
endif

! Need to know the Time Dimension ID and length

call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
        'get_var_3d', 'inq_dimid time '//trim(varname))
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen), &
        'get_var_3d', 'inquire_dimension time '//trim(varname))

call nc_check(nf90_inq_varid(ncid, varname, VarID), 'get_var_3d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_3d', 'inquire_variable '//trim(varname))

ncstart(:) = 1
nccount(:) = 1

DimCheck : do i = 1,numdims

   write(string1,'(a,i2,1x,A)') 'inquire dimension ',i,trim(varname)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
           'get_var_3d', string1)

   if ( dimIDs(i) == TimeDimID ) then

      timeindex  = find_desired_time_index(ncid, time_dimlen)
      ncstart(i) = timeindex
      dimlens(i) = 1

   elseif ( size(var3d,i) /= dimlens(i) ) then
      write(string1,*) trim(varname)//' has shape ', dimlens(1:numdims)
      write(string2,*) 'which is not the expected shape of ', size(var3d,i)
      call error_handler(E_ERR,'get_var_3d',string1,source,revision,revdate,text2=string2)
   endif

   nccount(i) = dimlens(i)

enddo DimCheck

if (do_output() .and. (debug > 3)) then
   write(*,*)'get_var_3d: variable ['//trim(varname)//']'
   write(*,*)'get_var_3d: start ',ncstart(1:numdims)
   write(*,*)'get_var_3d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2),dimlens(3)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//trim(varname))
   var3d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1),dimlens(2),dimlens(3)))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//trim(varname))
   var3d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var3d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//trim(varname))

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var3d == spvalR8) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var3d == spvalR8) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_3d',string1,source,revision,revdate)
endif

end subroutine get_var_3d


!-----------------------------------------------------------------------
!>
!> rigorous test of the interpolation routine.

subroutine test_interpolation(test_casenum)

integer, intent(in) :: test_casenum

! ! This is storage just used for temporary test driver
! integer :: imain, jmain, index, istatus, nxp1_temp, nyp1_temp
! integer :: dnxp1_temp, dnyp1_temp, height
! real(r8) :: ti, tj
!
! ! Storage for testing interpolation to a different grid
! integer :: dnxp1, dnyp1
! real(r8), allocatable :: dulon(:, :), dulat(:, :), dtlon(:, :), dtlat(:, :)
! real(r8), allocatable :: reg_u_data(:, :), reg_t_data(:, :)
! real(r8), allocatable :: reg_u_x(:), reg_t_x(:), dipole_u(:, :), dipole_t(:, :)
!
! ! Test program reads in two grids; one is the interpolation grid
! ! The second is a grid to which to interpolate.
!
! ! Read of first grid
! ! Set of block options for now
! ! Lon lat for u on channel 12, v on 13, u data on 14, t data on 15
!
! ! Case 1: regular grid to dipole
! ! Case 2: dipole to regular grid
! ! Case 3: regular grid to regular grid with same grid as dipole in SH
! ! Case 4: regular grid with same grid as dipole in SH to regular grid
! ! Case 5: regular grid with same grid as dipole in SH to dipole
! ! Case 6: dipole to regular grid with same grid as dipole in SH
!
! if ( .not. module_initialized ) call static_init_model
!
! call error_handler(E_MSG,'test_interpolate','FIXME TJH UNTESTED')
!
! if(test_casenum == 1 .or. test_casenum == 3) then
!    ! Case 1 or 3: read in from regular grid
!    open(unit=12, position='rewind', action='read', file='regular_grid_u')
!    open(unit=13, position='rewind', action='read', file='regular_grid_t')
!    open(unit=14, position='rewind', action='read', file='regular_grid_u_data')
!    open(unit=15, position='rewind', action='read', file='regular_grid_t_data')
!
! else if(test_casenum == 4 .or. test_casenum == 5) then
!    ! Case 3 or 4: read regular with same grid as dipole in SH
!    open(unit=12, position='rewind', action='read', file='regular_griddi_u')
!    open(unit=13, position='rewind', action='read', file='regular_griddi_t')
!    open(unit=14, position='rewind', action='read', file='regular_griddi_u_data')
!    open(unit=15, position='rewind', action='read', file='regular_griddi_t_data')
!
! else if(test_casenum == 2 .or. test_casenum == 6) then
!    ! Case 5: read in from dipole grid
!    open(unit=12, position='rewind', action='read', file='dipole_grid_u')
!    open(unit=13, position='rewind', action='read', file='dipole_grid_t')
!    open(unit=14, position='rewind', action='read', file='dipole_grid_u_data')
!    open(unit=15, position='rewind', action='read', file='dipole_grid_t_data')
! endif
!
! ! Get the size of the grid from the input u and t files
! read(12, *) nxp1, nyp1
! read(13, *) nxp1_temp, nyp1_temp
! if(nxp1 /= nxp1_temp .or. nyp1 /= nyp1_temp) then
!    write(string1,*)'mismatch nxp1,nxp1_temp ',nxp1,nxp1_temp,' or nyp1,nyp1_temp',nyp1,nyp1_temp
!    call error_handler(E_ERR,'test_interpolation',string1,source,revision,revdate)
! endif
!
! ! Allocate stuff for the first grid (the one being interpolated from)
! allocate(ulon(nxp1, nyp1), ulat(nxp1, nyp1), tlon(nxp1, nyp1), tlat(nxp1, nyp1))
! allocate(reg_u_data(nxp1, nyp1), reg_t_data(nxp1, nyp1))
! ! The Dart format 1d data arrays
! allocate(reg_u_x(nxp1*nyp1), reg_t_x(nxp1*nyp1))
!
! do imain = 1, nxp1
!    do jmain = 1, nyp1
!       read(12, *) ti, tj, ulon(imain, jmain), ulat(imain, jmain)
!       read(13, *) ti, tj, tlon(imain, jmain), tlat(imain, jmain)
!       read(14, *) ti, tj, reg_u_data(imain, jmain)
!       read(15, *) ti, tj, reg_t_data(imain, jmain)
!    enddo
! enddo
!
! ! Load into 1D dart data array
! index = 0
! do jmain = 1, nyp1
!    do imain = 1, nxp1
!       index = index + 1
!       reg_u_x(index) = reg_u_data(imain, jmain)
!       reg_t_x(index) = reg_t_data(imain, jmain)
!    enddo
! enddo
!
! ! dummy out vertical; let height always = 1 and allow
! ! all grid points to be processed.
! kmt = 2
! kmu = 2
! height = 1
!
! ! Initialize the interpolation data structure for this grid.
! call init_interp()
!
! ! Now read in the points for the output grid
! ! Case 1: regular grid to dipole
! ! Case 2: dipole to regular grid
! ! Case 3: regular grid to regular grid with same grid as dipole in SH
! ! Case 4: regular grid with same grid as dipole in SH to regular grid
! ! Case 5: regular grid with same grid as dipole in SH to dipole
! ! Case 6: dipole to regular grid with same grid as dipole in SH
! if(test_casenum == 1 .or. test_casenum == 5) then
!    ! Output to dipole grid
!    open(unit=22, position='rewind', action='read',  file='dipole_grid_u')
!    open(unit=23, position='rewind', action='read',  file='dipole_grid_t')
!    open(unit=24, position='rewind', action='write', file='dipole_grid_u_data.out')
!    open(unit=25, position='rewind', action='write', file='dipole_grid_t_data.out')
!
! else if(test_casenum == 2 .or. test_casenum == 4) then
!    ! Output to regular grid
!    open(unit=22, position='rewind', action='read',  file='regular_grid_u')
!    open(unit=23, position='rewind', action='read',  file='regular_grid_t')
!    open(unit=24, position='rewind', action='write', file='regular_grid_u_data.out')
!    open(unit=25, position='rewind', action='write', file='regular_grid_t_data.out')
!
! else if(test_casenum == 3 .or. test_casenum == 6) then
!    ! Output to regular grid with same grid as dipole in SH
!    open(unit=22, position='rewind', action='read',  file='regular_griddi_u')
!    open(unit=23, position='rewind', action='read',  file='regular_griddi_t')
!    open(unit=24, position='rewind', action='write', file='regular_griddi_u_data.out')
!    open(unit=25, position='rewind', action='write', file='regular_griddi_t_data.out')
! endif
!
! read(22, *) dnxp1, dnyp1
! read(23, *) dnxp1_temp, dnyp1_temp
! if(dnxp1 /= dnxp1_temp .or. dnyp1 /= dnyp1_temp) then
!    write(string1,*)'mismatch dnxp1,dnxp1_temp ',dnxp1,dnxp1_temp,' or dnyp1,dnyp1_temp',dnyp1,dnyp1_temp
!    call error_handler(E_ERR,'test_interpolation',string1,source,revision,revdate)
! endif
!
! allocate(dulon(dnxp1, dnyp1), dulat(dnxp1, dnyp1), dtlon(dnxp1, dnyp1), dtlat(dnxp1, dnyp1))
! allocate(dipole_u(dnxp1, dnyp1), dipole_t(dnxp1, dnyp1))
!
! dipole_u = 0.0_r8   ! just put some dummy values in to make sure they get changed.
! dipole_t = 0.0_r8   ! just put some dummy values in to make sure they get changed.
!
! do imain = 1, dnxp1
! do jmain = 1, dnyp1
!    read(22, *) ti, tj, dulon(imain, jmain), dulat(imain, jmain)
!    read(23, *) ti, tj, dtlon(imain, jmain), dtlat(imain, jmain)
! enddo
! enddo
!
! do imain = 1, dnxp1
!    do jmain = 1, dnyp1
!       ! Interpolate U from the first grid to the second grid
!
!       call lon_lat_interpolate(reg_u_x, dulon(imain, jmain), dulat(imain, jmain), &
!          KIND_U_CURRENT_COMPONENT, height, dipole_u(imain, jmain), istatus)
!
!       if ( istatus /= 0 ) then
!          write(string1,'(''cell'',i4,i4,1x,f12.8,1x,f12.8,'' U interp failed - code '',i4)') &
!               imain, jmain, dulon(imain, jmain), dulat(imain, jmain), istatus
!          call error_handler(E_MSG,'test_interpolation',string1)
!       endif
!
!       write(24, *) dulon(imain, jmain), dulat(imain, jmain), dipole_u(imain, jmain)
!
!       ! Interpolate U from the first grid to the second grid
!
!       call lon_lat_interpolate(reg_t_x, dtlon(imain, jmain), dtlat(imain, jmain), &
!          KIND_POTENTIAL_TEMPERATURE, height, dipole_t(imain, jmain), istatus)
!
!       if ( istatus /= 0 ) then
!          write(string1,'(''cell'',i4,i4,1x,f12.8,1x,f12.8,'' T interp failed - code '',i4)') &
!               imain,jmain, dtlon(imain, jmain), dtlat(imain, jmain), istatus
!          call error_handler(E_MSG,'test_interpolation',string1)
!       endif
!
!       write(25, *) dtlon(imain, jmain), dtlat(imain, jmain), dipole_t(imain, jmain)
!
!    enddo
! enddo

call error_handler(E_ERR,'test_interpolation','not written yet', source, revision, revdate)

! just to silence compiler
write(*,*)'string1 ',test_casenum

end subroutine test_interpolation


!-----------------------------------------------------------------------
!>
!> write the time information to a cookie file so we can use it to
!> control the length of the forecast from GCOM

subroutine write_gcom_timeinfo(current_time, forecast_stop_time)
type(time_type), intent(in) :: current_time
type(time_type), intent(in) :: forecast_stop_time

integer :: iunit
integer :: forecast_days, forecast_seconds
type(time_type) :: forecast_duration
integer :: iyear, imonth, iday, ihour, imin, iseconds, seconds

if ( .not. module_initialized ) call static_init_model

iunit = open_file('dart_gcom_timeinfo.txt',form='formatted',action='write')

call get_date(current_time, iyear, imonth, iday, ihour, imin, iseconds)
write(iunit,100)"Start_Time = 'seconds since ",iyear,imonth,iday,ihour,imin,iseconds
!write(iunit,100)"Start_Time = ", iyear,imonth,iday,ihour,imin,iseconds
100 format (A,i4.4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2,"'")

forecast_duration = forecast_stop_time - current_time
call get_time(forecast_duration,  forecast_seconds, forecast_days)
write(iunit,'(''Stop_Time_sec = '',f16.1)') forecast_seconds+forecast_days*86400.0_r8

seconds = ihour*3600 + imin*60 + iseconds
write(iunit,'(''currenttimetag = '',i4.4,''-''i2.2,''-''i2.2,''-''i5.5)') &
            iyear,imonth,iday,seconds

call get_date(forecast_stop_time, iyear, imonth, iday, ihour, imin, iseconds)
seconds = ihour*3600 + imin*60 + iseconds
write(iunit,'(''forecasttimetag = '',i4.4,''-''i2.2,''-''i2.2,''-''i5.5)') &
            iyear,imonth,iday,seconds

write(iunit,*)
call print_date(      current_time,'current date',iunit=iunit)
call print_date(forecast_stop_time,'desired date',iunit=iunit)
call print_time(      current_time,'current time',iunit=iunit)
call print_time(forecast_stop_time,'desired time',iunit=iunit)
call print_time( forecast_duration,'forecast duration (days,seconds)',iunit=iunit)
call close_file(iunit)

end subroutine write_gcom_timeinfo


!=======================================================================
! The remaining (private) interfaces come last.
!=======================================================================


!-----------------------------------------------------------------------
! overloaded routines defining vector_to_variable
!          MODULE PROCEDURE vector_to_1d_variable
!          MODULE PROCEDURE vector_to_2d_variable
!          MODULE PROCEDURE vector_to_3d_variable
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!>
!> Stuff the values from a 1d array, starting at an offset,
!> into a 1d array.  not particularly challenging

subroutine vector_to_1d_variable(x, varindex, data_1d_array, ncid)

real(r8),          intent(in)  :: x(:)
integer,           intent(in)  :: varindex
real(r8),          intent(out) :: data_1d_array(:)
integer, OPTIONAL, intent(in)  :: ncid

integer :: array_length

array_length = size(data_1d_array,1)

if ( array_length /= progvar(varindex)%varsize ) then
   write(string1,*)trim(progvar(varindex)%varname),' length /= array_length ', &
                        progvar(varindex)%varsize,' /= ',array_length
   call error_handler(E_ERR,'vector_to_1d_variable',string1,source,revision,revdate)
endif

data_1d_array = x(progvar(varindex)%index1:progvar(varindex)%indexN)

if (present(ncid)) then

   write(string1,'(A)') '... FIXME clamping '//trim(progvar(varindex)%varname)//' to a static value.'
   string2 = 'FIXME should be / could be adding some variability to it.'
   call error_handler(E_MSG,'vector_to_1d_variable',string1,text2=string2)

   if ((progvar(varindex)%rangeRestricted == BOUNDED_ABOVE ) .or. &
       (progvar(varindex)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_1d_array /= MISSING_R8) .and. &
             (data_1d_array > progvar(varindex)%maxvalue)) &
              data_1d_array = progvar(varindex)%maxvalue ! FIXME ... add variability
   endif

   if ((progvar(varindex)%rangeRestricted == BOUNDED_BELOW ) .or. &
       (progvar(varindex)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_1d_array /= MISSING_R8) .and. &
             (data_1d_array < progvar(varindex)%minvalue)) &
              data_1d_array = progvar(varindex)%minvalue ! FIXME ... add variability
   endif

endif

end subroutine vector_to_1d_variable


!-----------------------------------------------------------------------
!>
!> Stuff the values from a 1d fortran array, starting at an offset,
!> into a 2d fortran array.  the 2 dims are taken from the array size.

subroutine vector_to_2d_variable(x, varindex, data_2d_array, ncid)

real(r8),          intent(in)  :: x(:)
integer,           intent(in)  :: varindex
real(r8),          intent(out) :: data_2d_array(:,:)
integer, OPTIONAL, intent(in)  :: ncid

integer :: i,j,ii
integer :: dim1,dim2

call error_handler(E_MSG,'vector_to_2d_variable','routine not tested', &
                      source, revision, revdate)

dim1 = size(data_2d_array,1)
dim2 = size(data_2d_array,2)

if (dim1 /= progvar(varindex)%dimlens(1)) then
   write(string1,*)trim(progvar(varindex)%varname),' 2d array dim 1 ',dim1, &
                ' /= ', progvar(varindex)%dimlens(1)
   call error_handler(E_ERR,'vector_to_2d_variable',string1,source,revision,revdate)
endif
if (dim2 /= progvar(varindex)%dimlens(2)) then
   write(string1,*)trim(progvar(varindex)%varname),' 2d array dim 2 ',dim2, &
                ' /= ', progvar(varindex)%dimlens(2)
   call error_handler(E_ERR,'vector_to_2d_variable',string1,source,revision,revdate)
endif

ii = progvar(varindex)%index1

do j = 1,dim2
do i = 1,dim1
   data_2d_array(i,j) = x(ii)
   ii = ii + 1
enddo
enddo

if (present(ncid)) then

   write(string1,'(A)') '... FIXME clamping '//trim(progvar(varindex)%varname)//' to a static value.'
   string2 = 'FIXME should be / could be adding some variability to it.'
   call error_handler(E_MSG,'vector_to_2d_variable',string1,text2=string2)

   if ((progvar(varindex)%rangeRestricted == BOUNDED_ABOVE ) .or. &
       (progvar(varindex)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_2d_array /= MISSING_R8) .and. &
             (data_2d_array > progvar(varindex)%maxvalue)) &
              data_2d_array = progvar(varindex)%maxvalue ! FIXME ... add variability
   endif

   if ((progvar(varindex)%rangeRestricted == BOUNDED_BELOW ) .or. &
       (progvar(varindex)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_2d_array /= MISSING_R8) .and. &
             (data_2d_array < progvar(varindex)%minvalue)) &
              data_2d_array = progvar(varindex)%minvalue ! FIXME ... add variability
   endif

endif

end subroutine vector_to_2d_variable


!-----------------------------------------------------------------------
!>
!> Stuff the values from a 1d fortran array, starting at an offset,
!> into a 3d fortran array.  the 3 dims are taken from the array size.

subroutine vector_to_3d_variable(x, varindex, data_3d_array, ncid)

real(r8),          intent(in)  :: x(:)
integer,           intent(in)  :: varindex
real(r8),          intent(out) :: data_3d_array(:,:,:)
integer, OPTIONAL, intent(in)  :: ncid

integer :: i,j,k,ii
integer :: dim1,dim2,dim3

dim1 = size(data_3d_array,1)
dim2 = size(data_3d_array,2)
dim3 = size(data_3d_array,3)

if (dim1 /= progvar(varindex)%dimlens(1)) then
   write(string1,*)trim(progvar(varindex)%varname),' 3d array dim 1 ',dim1, &
                ' /= ', progvar(varindex)%dimlens(1)
   call error_handler(E_ERR,'vector_to_3d_variable',string1,source,revision,revdate)
endif
if (dim2 /= progvar(varindex)%dimlens(2)) then
   write(string1,*)trim(progvar(varindex)%varname),' 3d array dim 2 ',dim2, &
                ' /= ', progvar(varindex)%dimlens(2)
   call error_handler(E_ERR,'vector_to_3d_variable',string1,source,revision,revdate)
endif
if (dim3 /= progvar(varindex)%dimlens(3)) then
   write(string1,*)trim(progvar(varindex)%varname),' 3d array dim 3 ',dim3, &
                ' /= ', progvar(varindex)%dimlens(3)
   call error_handler(E_ERR,'vector_to_3d_variable',string1,source,revision,revdate)
endif

ii = progvar(varindex)%index1

do k = 1,dim3
do j = 1,dim2
do i = 1,dim1
   data_3d_array(i,j,k) = x(ii)
   ii = ii + 1
enddo
enddo
enddo

if (present(ncid)) then

   write(string1,'(A)') '... FIXME clamping '//trim(progvar(varindex)%varname)//' to a static value.'
   string2 = 'FIXME should be / could be adding some variability to it.'
   call error_handler(E_MSG,'vector_to_3d_variable',string1,text2=string2)

   if ((progvar(varindex)%rangeRestricted == BOUNDED_ABOVE ) .or. &
       (progvar(varindex)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_3d_array /= MISSING_R8) .and. &
             (data_3d_array > progvar(varindex)%maxvalue)) &
              data_3d_array = progvar(varindex)%maxvalue ! FIXME ... add variability
   endif

   if ((progvar(varindex)%rangeRestricted == BOUNDED_BELOW ) .or. &
       (progvar(varindex)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_3d_array /= MISSING_R8) .and. &
             (data_3d_array < progvar(varindex)%minvalue)) &
              data_3d_array = progvar(varindex)%minvalue ! FIXME ... add variability
   endif

endif

end subroutine vector_to_3d_variable


!-----------------------------------------------------------------------
!>
!> Initializes data structures needed for GCOM interpolation for
!> either dipole or irregular grid.
!> This should be called at static_init_model time to avoid
!> having all this temporary storage in the middle of a run.




!-----------------------------------------------------------------------
!>
!>  Build the data structure for interpolation for a dipole grid.

!-----------------------------------------------------------------------
!>
!> Given a longitude and latitude in degrees returns the index of the
!> regular lon-lat box that contains the point.

subroutine get_reg_box_indices(lon, lat, x_ind, y_ind)

real(r8), intent(in)  :: lon, lat
integer,  intent(out) :: x_ind, y_ind

call error_handler(E_ERR,'get_reg_box_indices','routine not written yet', &
                      source, revision, revdate)

call get_reg_lon_box(lon, x_ind)
call get_reg_lat_box(lat, y_ind)

end subroutine get_reg_box_indices


!-----------------------------------------------------------------------
!>
!> Determine which regular longitude box a longitude is in.

subroutine get_reg_lon_box(lon, x_ind)

real(r8), intent(in)  :: lon
integer,  intent(out) :: x_ind

call error_handler(E_ERR,'get_reg_lon_box','routine not written yet', &
                      source, revision, revdate)

x_ind = int(num_reg_x * lon / 360.0_r8) + 1

! Watch out for exactly at top; assume all lats and lons in legal range
if(lon == 360.0_r8) x_ind = num_reg_x

end subroutine get_reg_lon_box


!-----------------------------------------------------------------------
!>
!> Determine which regular latitude box a latitude is in.

subroutine get_reg_lat_box(lat, y_ind)

real(r8), intent(in)  :: lat
integer,  intent(out) :: y_ind

call error_handler(E_ERR,'get_reg_lat_box','routine not written yet', &
                      source, revision, revdate)

y_ind = int(num_reg_y * (lat + 90.0_r8) / 180.0_r8) + 1

! Watch out for exactly at top; assume all lats and lons in legal range
if(lat == 90.0_r8)  y_ind = num_reg_y

end subroutine get_reg_lat_box


!-----------------------------------------------------------------------
!>
!> Find a set of regular lat lon boxes that covers all of the area covered by
!> a dipole grid qaud whose corners are given by the dimension four x_corners
!> and y_corners arrays.  The two dimensional arrays reg_lon_ind and reg_lat_ind
!> return the first and last indices of the regular boxes in latitude and
!> longitude respectively. These indices may wraparound for reg_lon_ind.
!> A special computation is needed for a dipole quad that has the true north
!> pole in its interior. The logical is_pole is set to true if this is the case.
!> This can only happen for the t grid.  If the longitude boxes overlap 0
!> degrees, the indices returned are adjusted by adding the total number of
!> boxes to the second index (e.g. the indices might be 88 and 93 for a case
!> with 90 longitude boxes).

subroutine reg_box_overlap(x_corners, y_corners, is_pole, reg_lon_ind, reg_lat_ind)

real(r8), intent(in)  :: x_corners(4), y_corners(4)
logical,  intent(in)  :: is_pole
integer,  intent(out) :: reg_lon_ind(2), reg_lat_ind(2)

real(r8) :: lat_min, lat_max, lon_min, lon_max
integer  :: i

call error_handler(E_ERR,'reg_box_overlap','routine not written yet', &
                      source, revision, revdate)

!  A quad containing the pole is fundamentally different
if(is_pole) then
   ! Need all longitude boxes
   reg_lon_ind(1) = 1
   reg_lon_ind(2) = num_reg_x
   ! Need to cover from lowest latitude to top box
   lat_min = minval(y_corners)
   reg_lat_ind(1) = int(num_reg_y * (lat_min + 90.0_r8) / 180.0_r8) + 1
   call get_reg_lat_box(lat_min, reg_lat_ind(1))
   reg_lat_ind(2) = num_reg_y
else
   ! All other quads do not contain pole (pole could be on edge but no problem)
   ! This is specific to the dipole GCOM grids that do not go to the south pole
   ! Finding the range of latitudes is cake
   lat_min = minval(y_corners)
   lat_max = maxval(y_corners)

   ! Figure out the indices of the regular boxes for min and max lats
   call get_reg_lat_box(lat_min, reg_lat_ind(1))
   call get_reg_lat_box(lat_max, reg_lat_ind(2))

   ! Lons are much trickier. Need to make sure to wraparound the
   ! right way. There is no guarantee on direction of lons in the
   ! high latitude dipole rows.
   ! All longitudes for non-pole rows have to be within 180 degrees
   ! of one another.
   lon_min = minval(x_corners)
   lon_max = maxval(x_corners)
   if((lon_max - lon_min) > 180.0_r8) then
      ! If the max longitude value is more than 180
      ! degrees larger than the min, then there must be wraparound.
      ! Then, find the smallest value > 180 and the largest < 180 to get range.
      lon_min = 360.0_r8
      lon_max = 0.0_r8
      do i=1, 4
         if(x_corners(i) > 180.0_r8 .and. x_corners(i) < lon_min) lon_min = x_corners(i)
         if(x_corners(i) < 180.0_r8 .and. x_corners(i) > lon_max) lon_max = x_corners(i)
      enddo
   endif

   ! Get the indices for the extreme longitudes
   call get_reg_lon_box(lon_min, reg_lon_ind(1))
   call get_reg_lon_box(lon_max, reg_lon_ind(2))

   ! Watch for wraparound again; make sure that second index is greater than first
   if(reg_lon_ind(2) < reg_lon_ind(1)) reg_lon_ind(2) = reg_lon_ind(2) + num_reg_x
endif

end subroutine reg_box_overlap


!-----------------------------------------------------------------------
!>
!> Grabs the corners for a given quadrilateral from the global array of lower
!> right corners. Note that corners go counterclockwise around the quad.

subroutine get_quad_corners(x, i, j, corners)
real(r8), intent(in)  :: x(:, :, :)
integer,  intent(in)  :: i
integer,  intent(in)  :: j
real(r8), intent(out) :: corners(4)

call error_handler(E_ERR,'get_quad_corners','routine not written yet', &
                      source, revision, revdate)

! just to silence compiler
corners(1) = i
corners(2) = j
corners(3) = x(1,1,1)

end subroutine get_quad_corners


!-----------------------------------------------------------------------
!>
!> Updates the data structure listing dipole quads that are in a given regular box

subroutine update_reg_list(reg_list_num, reg_list_lon, reg_list_lat, &
                           reg_lon_ind, reg_lat_ind, dipole_lon_index, dipole_lat_index)

integer, intent(inout) :: reg_list_num(:, :)
integer, intent(inout) :: reg_list_lon(:, :, :)
integer, intent(inout) :: reg_list_lat(:, :, :)
integer, intent(inout) :: reg_lon_ind(2)
integer, intent(inout) :: reg_lat_ind(2)
integer, intent(in)    :: dipole_lon_index
integer, intent(in)    :: dipole_lat_index

integer :: ind_x, index_x, ind_y

call error_handler(E_ERR,'update_reg_list','routine not written yet', &
                      source, revision, revdate)

! Loop through indices for each possible regular cell
! Have to watch for wraparound in longitude
if(reg_lon_ind(2) < reg_lon_ind(1)) reg_lon_ind(2) = reg_lon_ind(2) + num_reg_x

do ind_x = reg_lon_ind(1), reg_lon_ind(2)
   ! Inside loop, need to go back to wraparound indices to find right box
   index_x = ind_x
   if(index_x > num_reg_x) index_x = index_x - num_reg_x

   do ind_y = reg_lat_ind(1), reg_lat_ind(2)
      ! Make sure the list storage isn't full
      if(reg_list_num(index_x, ind_y) >= max_reg_list_num) then
         write(string1,*) 'max_reg_list_num (',max_reg_list_num,') is too small ... increase'
         call error_handler(E_ERR, 'update_reg_list', string1, source, revision, revdate)
      endif

      ! Increment the count
      reg_list_num(index_x, ind_y) = reg_list_num(index_x, ind_y) + 1
      ! Store this quad in the list for this regular box
      reg_list_lon(index_x, ind_y, reg_list_num(index_x, ind_y)) = dipole_lon_index
      reg_list_lat(index_x, ind_y, reg_list_num(index_x, ind_y)) = dipole_lat_index
   enddo
enddo

end subroutine update_reg_list


!-----------------------------------------------------------------------
!>
!> Three different types of grids are used here. The GCOM dipole
!> grid is referred to as a dipole grid and each region is
!> referred to as a quad, short for quadrilateral.
!> The longitude latitude rectangular grid with possibly irregular
!> spacing in latitude used for some GCOM applications and testing
!> is referred to as the irregular grid and each region is
!> called a box.
!> Finally, a regularly spaced longitude latitude grid is used
!> as a computational tool for interpolating from the dipole
!> grid. This is referred to as the regular grid and each region
!> is called a box.
!> All grids are referenced by the index of the lower left corner
!> of the quad or box.
!>
!> The dipole grid is assumed to be global for all applications.
!> The irregular grid is also assumed to be global east
!> west for all applications.
!>
!------------------------------------------------------------------


!------------------------------------------------------------------
!>
!> Subroutine to interpolate to a lon lat location given the state vector
!> for that level, x. This works just on one horizontal slice.
!> NOTE: Using array sections to pass in the x array may be inefficient on some
!> compiler/platform setups. Might want to pass in the entire array with a base
!> offset value instead of the section if this is an issue.
!> This routine works for either the dipole or a regular lat-lon grid.
!> Successful interpolation returns istatus=0.

subroutine lon_lat_interpolate(x, lon, lat, var_type, height, interp_val, istatus)

real(r8), intent(in)  :: x(:)
real(r8), intent(in)  :: lon
real(r8), intent(in)  :: lat
integer,  intent(in)  :: var_type
integer,  intent(in)  :: height
real(r8), intent(out) :: interp_val
integer,  intent(out) :: istatus


integer  :: lat_bot, lat_top, lon_bot, lon_top, num_inds, start_ind
integer  :: x_ind, y_ind
real(r8) :: p(4), x_corners(4), y_corners(4), xbot, xtop
real(r8) :: lon_fract, lat_fract
logical  :: masked

call error_handler(E_ERR,'lon_lat_interpolate','routine not written yet', &
                      source, revision, revdate)

! Succesful return has istatus of 0
istatus = 0

! Get the lower left corner for either grid type
if(dipole_grid) then
   ! Figure out which of the regular grid boxes this is in
   call get_reg_box_indices(lon, lat, x_ind, y_ind)

   ! Is this on the U or T grid?
   if(is_on_ugrid(var_type)) then
      ! On U grid
      num_inds =  u_dipole_num  (x_ind, y_ind)
      start_ind = u_dipole_start(x_ind, y_ind)

      ! If there are no quads overlapping, can't do interpolation
      if(num_inds == 0) then
         istatus = 1
         return
      endif

      ! Search the list of quads to see if (lon, lat) is in one
      !call get_dipole_quad(lon, lat, ulon, ulat, num_inds, start_ind, &
        ! u_dipole_lon_list, u_dipole_lat_list, lon_bot, lat_bot, istatus)
      ! Fail on bad istatus return
      if(istatus /= 0) return

      ! Getting corners for accurate interpolation
      !call get_quad_corners(ulon, lon_bot, lat_bot, x_corners)
      !call get_quad_corners(ulat, lon_bot, lat_bot, y_corners)

      ! Fail if point is in one of the U boxes that go through the
      ! pole (this could be fixed up if necessary)
      if(lat_bot == u_pole_y .and. (lon_bot == pole_x -1 .or. &
         lon_bot == pole_x)) then
         istatus = 4
         return
      endif

   else
      ! On T grid
      num_inds =  t_dipole_num  (x_ind, y_ind)
      start_ind = t_dipole_start(x_ind, y_ind)
      !call get_dipole_quad(lon, lat, tlon, tlat, num_inds, start_ind, &
       !  t_dipole_lon_list, t_dipole_lat_list, lon_bot, lat_bot, istatus)
      ! Fail on bad istatus return
      if(istatus /= 0) return

      ! Fail if point is in T box that covers pole
      if(lon_bot == pole_x .and. lat_bot == t_pole_y) then
         istatus = 5
         return
      endif

      ! Getting corners for accurate interpolation
      !call get_quad_corners(tlon, lon_bot, lat_bot, x_corners)
      !call get_quad_corners(tlat, lon_bot, lat_bot, y_corners)
   endif

else
   ! This is an irregular grid
   ! U and V are on velocity grid
   if (is_on_ugrid(var_type)) then
      ! Get the corner indices and the fraction of the distance between
! FIXME      call get_irreg_box(lon, lat, ulon, ulat, &
! FIXME         lon_bot, lat_bot, lon_fract, lat_fract, istatus)
   else
      ! Eta, T and S are on the T grid
      ! Get the corner indices
! FIXME      call get_irreg_box(lon, lat, tlon, tlat, &
! FIXME         lon_bot, lat_bot, lon_fract, lat_fract, istatus)
   endif

   ! Return passing through error status
   if(istatus /= 0) return

endif

! Find the indices to get the values for interpolating
lat_top = lat_bot + 1
if(lat_top > nyp1) then
   istatus = 2
   return
endif

! Watch for wraparound in longitude
lon_top = lon_bot + 1
if(lon_top > nxp1) lon_top = 1

! Get the values at the four corners of the box or quad
! Corners go around counterclockwise from lower left
p(1) = get_val(lon_bot, lat_bot, nxp1, x, var_type, height, masked)
if(masked) then
   istatus = 3
   return
endif

p(2) = get_val(lon_top, lat_bot, nxp1, x, var_type, height, masked)
if(masked) then
   istatus = 3
   return
endif

p(3) = get_val(lon_top, lat_top, nxp1, x, var_type, height, masked)
if(masked) then
   istatus = 3
   return
endif

p(4) = get_val(lon_bot, lat_top, nxp1, x, var_type, height, masked)
if(masked) then
   istatus = 3
   return
endif

! Full bilinear interpolation for quads
if(dipole_grid) then
   call quad_bilinear_interp(lon, lat, x_corners, y_corners, p, interp_val)
else
   ! Rectangular biliear interpolation
   xbot = p(1) + lon_fract * (p(2) - p(1))
   xtop = p(4) + lon_fract * (p(3) - p(4))
   ! Now interpolate in latitude
   interp_val = xbot + lat_fract * (xtop - xbot)
endif

end subroutine lon_lat_interpolate


!-----------------------------------------------------------------------
!>
!> Returns the value from a single level array given the lat and lon indices
!> 'masked' returns true if this is NOT a valid grid location (e.g. land, or
!> below the ocean floor in shallower areas).

function get_val(lon_index, lat_index, nlon, x, var_type, height, masked)

integer, intent(in)  :: lon_index
integer, intent(in)  :: lat_index
integer, intent(in)  :: nlon
integer, intent(in)  :: var_type
integer, intent(in)  :: height
real(r8),intent(in)  :: x(:)
logical, intent(out) :: masked
real(r8)             :: get_val

call error_handler(E_ERR,'get_val','routine not written yet', &
                      source, revision, revdate)

! check the land/ocean bottom map and return if not valid water cell.
if(is_dry_land(var_type, lon_index, lat_index, height)) then
   masked = .true.
   get_val = MISSING_R8
   return
endif

! Layout has lons varying most rapidly
get_val = x((lat_index - 1) * nlon + lon_index)

! this is a valid ocean water cell, not land or below ocean floor
masked = .false.

end function get_val


!-----------------------------------------------------------------------
!>
!> Given a longitude and latitude of a point (lon and lat) and the
!> longitudes and latitudes of the lower left corner of the regular grid
!> boxes, gets the indices of the grid box that contains the point and
!> the fractions along each directrion for interpolation.

!TJHsubroutine get_irreg_box(lon, lat, lon_array, lat_array, &
!TJH                         found_x, found_y, lon_fract, lat_fract, istatus)
!TJH
!TJHreal(r8), intent(in)  :: lon
!TJHreal(r8), intent(in)  :: lat
!TJHreal(r8), intent(in)  :: lon_array(nxp1, nyp1)
!TJHreal(r8), intent(in)  :: lat_array(nxp1, nyp1)
!TJHreal(r8), intent(out) :: lon_fract
!TJHreal(r8), intent(out) :: lat_fract
!TJHinteger,  intent(out) :: found_x
!TJHinteger,  intent(out) :: found_y
!TJHinteger,  intent(out) :: istatus
!TJH
!TJH! Local storage
!TJHinteger  :: lat_status, lon_top, lat_top
!TJH
!TJHcall error_handler(E_ERR,'get_irreg_box','routine not written yet', &
!TJH                      source, revision, revdate)
!TJH
!TJH! Succesful return has istatus of 0
!TJHistatus = 0
!TJH
!TJH! Get latitude box boundaries
!TJHcall lat_bounds(lat, nyp1, lat_array, found_y, lat_top, lat_fract, lat_status)
!TJH
!TJH! Check for error on the latitude interpolation
!TJHif(lat_status /= 0) then
!TJH   istatus = 1
!TJH   return
!TJHendif
!TJH
!TJH! Find out what longitude box and fraction
!TJHcall lon_bounds(lon, nxp1, lon_array, found_x, lon_top, lon_fract)
!TJH
!TJHend subroutine get_irreg_box


!-----------------------------------------------------------------------
!>
!> Given a longitude lon, the array of longitudes for grid boundaries,
!> and the number of longitudes in the grid, returns the indices of the longitude
!> below and above the location longitude and the fraction of the distance
!> between. It is assumed that the longitude wraps around for a global grid.
!> Since longitude grids are going to be regularly spaced, this could be made
!> more efficient.
!> Algorithm fails for a silly grid that has only 2 longitudes separated by 180 degrees.

subroutine lon_bounds(lon, nlons, lon_array, bot, top, fract)

real(r8), intent( in) :: lon
integer,  intent( in) :: nlons
real(r8), intent( in) :: lon_array(:, :)
integer,  intent(out) :: bot
integer,  intent(out) :: top
real(r8), intent(out) :: fract

! Local storage
integer  :: i
real(r8) :: dist_bot, dist_top

call error_handler(E_ERR,'lon_bounds','routine not written yet', &
                      source, revision, revdate)

! This is inefficient, someone could clean it up since longitudes are regularly spaced
! But note that they don't have to start at 0
do i = 2, nlons
   dist_bot = lon_dist(lon, lon_array(i - 1, 1))
   dist_top = lon_dist(lon, lon_array(i, 1))
   if(dist_bot <= 0 .and. dist_top > 0) then
      bot = i - 1
      top = i
      fract = abs(dist_bot) / (abs(dist_bot) + dist_top)
      return
   endif
enddo

! Falling off the end means it's in between; wraparound
bot = nlons
top = 1
dist_bot = lon_dist(lon, lon_array(bot, 1))
dist_top = lon_dist(lon, lon_array(top, 1))
fract = abs(dist_bot) / (abs(dist_bot) + dist_top)

end subroutine lon_bounds


!-----------------------------------------------------------------------
!>
!> Given a latitude lat, the array of latitudes for grid boundaries, and the
!> number of latitudes in the grid, returns the indices of the latitude
!> below and above the location latitude and the fraction of the distance
!> between. istatus is returned as 0 unless the location latitude is
!> south of the southernmost grid point (1 returned) or north of the
!> northernmost (2 returned). If one really had lots of polar obs would
!> want to worry about interpolating around poles.

subroutine lat_bounds(lat, nlats, lat_array, bot, top, fract, istatus)

real(r8), intent( in) :: lat
integer,  intent( in) :: nlats
real(r8), intent( in) :: lat_array(:, :)
integer,  intent(out) :: bot
integer,  intent(out) :: top
real(r8), intent(out) :: fract
integer,  intent(out) :: istatus

! Local storage
integer :: i

call error_handler(E_ERR,'lat_bounds','routine not written yet', &
                      source, revision, revdate)

! Success should return 0, failure a positive number.
istatus = 0

! Check for too far south or north
if(lat < lat_array(1, 1)) then
   istatus = 1
   return
else if(lat > lat_array(1, nlats)) then
   istatus = 2
   return
endif

! In the middle, search through
do i = 2, nlats
   if(lat <= lat_array(1, i)) then
      bot = i - 1
      top = i
      fract = (lat - lat_array(1, bot)) / (lat_array(1, top) - lat_array(1, bot))
      return
   endif
enddo

! Shouldn't get here. Might want to fail really hard through error handler
istatus = 40

end subroutine lat_bounds


!-----------------------------------------------------------------------
!>
!> Returns the smallest signed distance between lon1 and lon2 on the sphere
!> If lon1 is less than 180 degrees east of lon2 the distance is negative
!> If lon1 is less than 180 degrees west of lon2 the distance is positive

function lon_dist(lon1, lon2)

real(r8), intent(in) :: lon1
real(r8), intent(in) :: lon2
real(r8)             :: lon_dist

call error_handler(E_ERR,'lon_dist','routine not written yet', &
                      source, revision, revdate)

lon_dist = lon2 - lon1
if(lon_dist >= -180.0_r8 .and. lon_dist <= 180.0_r8) then
   return
else if(lon_dist < -180.0_r8) then
   lon_dist = lon_dist + 360.0_r8
else
   lon_dist = lon_dist - 360.0_r8
endif

end function lon_dist


!-----------------------------------------------------------------------
!>
!> Given the lon and lat of a point, and a list of the
!> indices of the quads that might contain a point at (lon, lat), determines
!> which quad contains the point.  istatus is returned as 0 if all went
!> well and 1 if the point was not found to be in any of the quads.




!-----------------------------------------------------------------------
!>
!> Return in_quad true if the point (lon, lat) is in the quad with
!> the given corners.
!?
!> Do this by line tracing in latitude for now. For non-pole point, want a vertical
!> line from the lon, lat point to intersect a side of the quad both above
!> and below the point.

function in_quad(lon, lat, x_corners, y_corners)

real(r8), intent(in) :: lon
real(r8), intent(in) :: lat
real(r8), intent(in) :: x_corners(4)
real(r8), intent(in) :: y_corners(4)
logical              :: in_quad

real(r8) :: x(2), y(2)
logical  :: cant_be_in_box, in_box
integer  :: intercepts_above(4), intercepts_below(4), i
integer  :: num_above, num_below

call error_handler(E_ERR,'in_quad','routine not written yet', &
                      source, revision, revdate)

! Default answer is point is not in quad
in_quad = .false.

! Loop through the sides and compute intercept (if any) with a vertical line
! from the point. This line has equation x=lon.
do i = 1, 4
   ! Load up the sides endpoints
   if(i <= 3) then
      x(1:2) = x_corners(i:i+1)
      y(1:2) = y_corners(i:i+1)
   else
      x(1) = x_corners(4)
      x(2) = x_corners(1)
      y(1) = y_corners(4)
      y(2) = y_corners(1)
   endif

   ! Check to see how a vertical line from the point is related to this side
   call line_intercept(x, y, lon, lat, cant_be_in_box, in_box, intercepts_above(i), &
      intercepts_below(i))

   ! If cant_be_in_box is true, can return right away
   if(cant_be_in_box) then
      in_quad = .false.
      return
   ! Return true if on a side
   else if(in_box) then
      in_quad = .true.
      return
   endif

enddo

! See if the line intercepted a side of the quad both above and below
num_above = sum(intercepts_above)
num_below = sum(intercepts_below)

if(num_above > 0 .and. num_below > 0) then
   in_quad = .true.
endif

end function in_quad


!-----------------------------------------------------------------------
!>
!> Find the intercept of a vertical line from point (x_point, y_point) and
!> a line segment with endpoints side_x and side_y.
!> For a given side have endpoints (side_x1, side_y1) and (side_x2, side_y2)
!> so equation of segment is y = side_y1 + m(x-side_x1) for y
!> between side_y1 and side_y2.
!> Intersection of vertical line and line containing side
!> occurs at y = side_y1 + m(x_point - side_x1); need this
!> y to be between side_y1 and side_y2.
!> If the vertical line is colinear with the side but the point is not on the side, return
!> cant_be_in_box as true. If the point is on the side, return in_box true.
!> If the intersection of the vertical line and the side occurs at a point above
!> the given point, return 1 for intercept_above. If the intersection occurs
!> below, return 1 for intercept_below. If the vertical line does not intersect
!> the segment, return false and 0 for all intent out arguments.
!>
!> WARNING: CERTAINLY A PROBLEM FOR THE POLE BOX!!! POLE BOX COULD
!> HAVE SIDES THAT ARE LONGER THAN 180. For now pole boxes are excluded.
!>
!> This can probably be made much cleaner and more efficient.

subroutine line_intercept(side_x_in, side_y, x_point_in, y_point, &
                          cant_be_in_box, in_box, intercept_above, intercept_below)

real(r8), intent(in)  :: side_x_in(2)
real(r8), intent(in)  :: side_y(2)
real(r8), intent(in)  :: x_point_in
real(r8), intent(in)  :: y_point
logical,  intent(out) :: cant_be_in_box
logical,  intent(out) :: in_box
integer,  intent(out) :: intercept_above
integer,  intent(out) :: intercept_below

real(r8) :: slope, y_intercept, side_x(2), x_point

call error_handler(E_ERR,'line_intercept','routine not written yet', &
                      source, revision, revdate)

! May have to adjust the longitude intent in values, so copy
side_x = side_x_in
x_point = x_point_in

! See if the side wraps around in longitude
if(maxval(side_x) - minval(side_x) > 180.0_r8) then
   if(side_x(1) < 180.0_r8)  side_x(1) =  side_x(1) + 360.0_r8
   if(side_x(2) < 180.0_r8)  side_x(2) =  side_x(2) + 360.0_r8
   if(x_point < 180.0_r8) x_point = x_point + 360.0_r8
endif

! Initialize the default returns
cant_be_in_box   = .false.
in_box           = .false.
intercept_above = 0
intercept_below = 0

! First easy check, if x_point is not between endpoints of segment doesn't intersect
if(x_point < minval(side_x) .or. x_point > maxval(side_x)) return

! Otherwise line must intersect the segment

! First subblock, slope is undefined
if(side_x(2) == side_x(1)) then
   ! The line is colinear with the side
   ! If y_point is between endpoints then point is on this side
   if(y_point <= maxval(side_y) .and. y_point >= minval(side_y)) then
      in_box = .true.
      return
   ! If not on side but colinear with side, point cant be in quad
   else
      cant_be_in_box = .true.
      return
   endif

else

   ! Second possibility; slope is defined
   ! FIXME: watch out for numerical instability.
   ! near-zero x's and large side_y's may cause overflow
   slope = (side_y(2) - side_y(1)) / (side_x(2) - side_x(1))

   ! Intercept of vertical line through is at x_point and...
   y_intercept = side_y(1) + slope * (x_point - side_x(1))

   ! Intersects the segment, is it above, below, or at the point
   if(y_intercept == y_point) then
      in_box = .true.
      return
   else if(y_intercept > y_point) then
      intercept_above = 1
      return
   else
      intercept_below = 1
      return
   endif
endif

end subroutine line_intercept


!-----------------------------------------------------------------------
!>
!> Given a longitude and latitude (lon_in, lat), the longitude and
!> latitude of the 4 corners of a quadrilateral and the values at the
!> four corners, interpolates to (lon_in, lat) which is assumed to
!> be in the quad. This is done by bilinear interpolation, fitting
!> a function of the form a + bx + cy + dxy to the four points and
!> then evaluating this function at (lon, lat). The fit is done by
!> solving the 4x4 system of equations for a, b, c, and d. The system
!> is reduced to a 3x3 by eliminating a from the first three equations
!> and then solving the 3x3 before back substituting. There is concern
!> about the numerical stability of this implementation. Implementation
!> checks showed accuracy to seven decimal places on all tests.

subroutine quad_bilinear_interp(lon_in, lat, x_corners_in, y_corners, &
                                p, interp_val)

real(r8),  intent(in) :: lon_in
real(r8),  intent(in) :: lat
real(r8),  intent(in) :: x_corners_in(4)
real(r8),  intent(in) :: y_corners(4)
real(r8),  intent(in) :: p(4)
real(r8), intent(out) :: interp_val

integer :: i
real(r8) :: m(3, 3), v(3), r(3), a, x_corners(4), lon
! real(r8) :: lon_mean

call error_handler(E_ERR,'quad_bilinear_interp','routine not written yet', &
                      source, revision, revdate)

! Watch out for wraparound on x_corners.
lon = lon_in
x_corners = x_corners_in

! See if the side wraps around in longitude. If the corners longitudes
! wrap around 360, then the corners and the point to interpolate to
! must be adjusted to be in the range from 180 to 540 degrees.
if(maxval(x_corners) - minval(x_corners) > 180.0_r8) then
   if(lon < 180.0_r8) lon = lon + 360.0_r8
   do i = 1, 4
      if(x_corners(i) < 180.0_r8) x_corners(i) = x_corners(i) + 360.0_r8
   enddo
endif


!*******
! Problems with extremes in polar cell interpolation can be reduced
! by this block, but it is not clear that it is needed for actual
! ocean grid data
! Find the mean longitude of corners and remove
!!!lon_mean = sum(x_corners) / 4.0_r8
!!!x_corners = x_corners - lon_mean
!!!lon = lon - lon_mean
! Multiply everybody by the cos of the latitude
!!!do i = 1, 4
   !!!x_corners(i) = x_corners(i) * cos(y_corners(i) * deg2rad)
!!!enddo
!!!lon = lon * cos(lat * deg2rad)

!*******


! Fit a surface and interpolate; solve for 3x3 matrix
do i = 1, 3
   ! Eliminate a from the first 3 equations
   m(i, 1) = x_corners(i) - x_corners(i + 1)
   m(i, 2) = y_corners(i) - y_corners(i + 1)
   m(i, 3) = x_corners(i)*y_corners(i) - x_corners(i + 1)*y_corners(i + 1)
   v(i) = p(i) - p(i + 1)
enddo

! Solve the matrix for b, c and d
call mat3x3(m, v, r)

! r contains b, c, and d; solve for a
a = p(4) - r(1) * x_corners(4) - r(2) * y_corners(4) - &
   r(3) * x_corners(4)*y_corners(4)


!----------------- Implementation test block
! When interpolating on dipole x3 never exceeded 1e-9 error in this test
!!!do i = 1, 4
   !!!interp_val = a + r(1)*x_corners(i) + r(2)*y_corners(i)+ r(3)*x_corners(i)*y_corners(i)
   !!!if(abs(interp_val - p(i)) > 1e-9) then
      !!!write(*, *) 'large interp residual ', interp_val - p(i)
   !!!endif
!!!enddo
!----------------- Implementation test block


! Now do the interpolation
interp_val = a + r(1)*lon + r(2)*lat + r(3)*lon*lat

!********
! Avoid exceeding maxima or minima as stopgap for poles problem
! When doing bilinear interpolation in quadrangle, can get interpolated
! values that are outside the range of the corner values
if(interp_val > maxval(p)) then
   interp_val = maxval(p)
else if(interp_val < minval(p)) then
   interp_val = minval(p)
endif
!********

end subroutine quad_bilinear_interp


!-----------------------------------------------------------------------
!>
!> Solves rank 3 linear system mr = v for r
!> using Cramer's rule. This isn't the best choice
!> for speed or numerical stability so might want to replace
!> this at some point.

subroutine mat3x3(m, v, r)

real(r8),  intent(in) :: m(3, 3)
real(r8),  intent(in) :: v(3)
real(r8), intent(out) :: r(3)

real(r8) :: m_sub(3, 3), numer, denom
integer  :: i

! Compute the denominator, det(m)
denom = deter3(m)

! Loop to compute the numerator for each component of r
do i = 1, 3
   m_sub = m
   m_sub(:, i) = v
   numer = deter3(m_sub)
   r(i) = numer / denom
enddo

end subroutine mat3x3


!-----------------------------------------------------------------------
!>
!> Computes determinant of 3x3 matrix m

function deter3(m)

real(r8), intent(in) :: m(3, 3)
real(r8)             :: deter3

deter3 = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + &
         m(1,3)*m(2,1)*m(3,2) - m(3,1)*m(2,2)*m(1,3) - &
         m(1,1)*m(2,3)*m(3,2) - m(3,3)*m(2,1)*m(1,2)

end function deter3


!-----------------------------------------------------------------------
!>

subroutine height_bounds(lheight, nheights, hgt_array, bot, top, fract, istatus)

real(r8), intent(in)  :: lheight
integer,  intent(in)  :: nheights
real(r8), intent(in)  :: hgt_array(nheights)
integer,  intent(out) :: bot
integer,  intent(out) :: top
real(r8), intent(out) :: fract
integer,  intent(out) :: istatus

! Local variables
integer   :: i

call error_handler(E_ERR,'height_bounds','routine not written yet', &
                      source, revision, revdate)

! Succesful istatus is 0
! Make any failure here return istatus in the 20s
istatus = 0

! The zc array contains the depths of the center of the vertical grid boxes
! In this case (unlike how we handle the MIT depths), positive is really down.
! FIXME: in the MIT model, we're given box widths and we compute the centers,
! and we computed them with larger negative numbers being deeper.  Here,
! larger positive numbers are deeper.

! It is assumed that the top box is shallow and any observations shallower
! than the depth of this box's center are just given the value of the
! top box.
if(lheight <= hgt_array(1)) then
   top = 1
   bot = 2
   ! NOTE: the fract definition is the relative distance from bottom to top
   fract = 1.0_r8
if (do_output() .and. (debug > 5)) print *, 'above first level in height'
if (do_output() .and. (debug > 5)) print *, 'hgt_array, top, bot, fract=', hgt_array(1), top, bot, fract
   return
endif

! Search through the boxes
do i = 2, nheights
   ! If the location is shallower than this entry, it must be in this box
   if(lheight < hgt_array(i)) then
      top = i -1
      bot = i
      fract = (hgt_array(bot) - lheight) / (hgt_array(bot) - hgt_array(top))
if (do_output() .and. (debug > 5)) print *, 'i, hgt_array, top, bot, fract=', i, hgt_array(i), top, bot, fract
      return
   endif
enddo

! Falling off the end means the location is lower than the deepest height
istatus = 20

end subroutine height_bounds



!-----------------------------------------------------------------------
!>
!> Given an integer index into the state vector structure, returns the
!> associated array indices for lat, lon, and depth, as well as the type.

subroutine get_state_indices(index_in, lon_index, lat_index, lev_index, varindex)

integer, intent(in)    :: index_in
integer, intent(out)   :: lon_index
integer, intent(out)   :: lat_index
integer, intent(out)   :: lev_index
integer, intent(inout) :: varindex

integer :: n, offset, ndim1, ndim2

if ((index_in < progvar(1)%index1) .or. &
    (index_in > progvar(nfields)%indexN) ) then
   write(string1,*) 'desired index ',index_in
   write(string2,*) 'is not within the bounds of the DART address space.'
   call error_handler(E_ERR,'get_state_indices',string1, &
         source,revision,revdate,text2=string2)
endif

if (varindex < 0) then
   ! Find the DART variable that spans the index of interest.
   FindIndex : do n = 1,nfields
      if( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) then
         varindex = n
         exit FindIndex
      endif
   enddo FindIndex
endif

offset = index_in - progvar(varindex)%index1 + 1

! So now we know the variable and its shape. Remember that DART only
! stores a single timestep, so any time dimension is a singleton.
! Furthermore, netCDF requires the time/unlimited dimension be the LAST
! (in Fortran) dimension, so we can just focus on the first N dimensions.
! Relying on integer arithmetic.
!
! We know the storage order is lon,lat,lev ...
! FIXME ... may want to ensure at some point.

if     ( progvar(varindex)%rank == 1) then

   lon_index = offset

elseif ( progvar(varindex)%rank == 2) then

   ndim1 = progvar(varindex)%dimlens(1)

   lat_index = 1 + (offset - 1)/ndim1
   lon_index = offset - (lat_index - 1)*ndim1

elseif ( progvar(varindex)%rank == 3) then

   ndim1 = progvar(varindex)%dimlens(1)  ! num_fastest_dimension (Fortran leftmost)
   ndim2 = progvar(varindex)%dimlens(2)  ! num_next_fastest

   lev_index = 1 + (offset - 1)/(ndim1*ndim2)
   lat_index = 1 + (offset - (lev_index-1)*ndim1*ndim2 -1)/ndim1
   lon_index = offset - (lev_index-1)*ndim1*ndim2 - (lat_index-1)*ndim1

else
   write(string1,*) 'Does not support variables with rank ',progvar(varindex)%rank
   write(string2,*) 'variable is ',trim(progvar(varindex)%varname)
   call error_handler(E_ERR,'get_state_meta_data',string1, &
         source,revision,revdate,text2=string2)
endif

if (do_output() .and. (debug > 99)) then
   print *, 'get_state_indices: asking for meta data about index ', index_in
   print *, 'get_state_indices: lon, lat, lev index = ', lon_index, lat_index, lev_index
endif

end subroutine get_state_indices


!-----------------------------------------------------------------------
!>
!> Given the indices for the lat, lon, level (i,j,k)s
!>

function ijk_to_state_index(varindex, i_index, j_index, k_index)

integer, intent(in) :: varindex
integer, intent(in) :: i_index
integer, intent(in) :: j_index
integer, intent(in) :: k_index
integer             :: ijk_to_state_index

integer :: offset, ndim1, ndim2

offset = progvar(varindex)%index1 - 1

if     ( progvar(varindex)%rank == 1) then

   ijk_to_state_index = offset + i_index

elseif ( progvar(varindex)%rank == 2) then

   ndim1 = progvar(varindex)%dimlens(1)

   ijk_to_state_index = offset + j_index*ndim1 + i_index

elseif ( progvar(varindex)%rank == 3) then

   ndim1 = progvar(varindex)%dimlens(1)  ! num_fastest_dimension (Fortran leftmost)
   ndim2 = progvar(varindex)%dimlens(2)  ! num_next_fastest

   ijk_to_state_index = offset + (j_index-1)*ndim1 + (k_index-1)*ndim1*ndim2 + i_index

else
   write(string1,*) 'Does not support variables with rank ',progvar(varindex)%rank
   write(string2,*) 'variable is ',trim(progvar(varindex)%varname)
   call error_handler(E_ERR,'get_state_meta_data',string1, &
         source,revision,revdate,text2=string2)
endif

if (do_output() .and. (debug > 4)) then
   print *, 'ijk_to_state_index: lon, lat, lev index = ', i_index, j_index, k_index
   print *, 'ijk_to_state_index: returns       index = ', ijk_to_state_index
endif

end function ijk_to_state_index


!-----------------------------------------------------------------------
!>
!> Given an integer index into the state vector structure,
!> return the longitude, latitude and level

subroutine get_state_lonlatlev(varindex, lon_index, lat_index, lev_index, &
                               mylon, mylat, mylev)

integer,  intent(in)  :: varindex
integer,  intent(in)  :: lon_index
integer,  intent(in)  :: lat_index
integer,  intent(in)  :: lev_index
real(r8), intent(out) :: mylon
real(r8), intent(out) :: mylat
real(r8), intent(out) :: mylev

if     ((trim(progvar(varindex)%coordinates) == 'lon lat lev') .or. &
        (trim(progvar(varindex)%coordinates) == 'lon lat lev time')) then
   mylon =  LON(lon_index, lat_index, lev_index)
   mylat =  LAT(lon_index, lat_index, lev_index)
   mylev =  LEV(lon_index, lat_index, lev_index)
elseif ((trim(progvar(varindex)%coordinates) == 'ulon ulat ulev') .or. &
        (trim(progvar(varindex)%coordinates) == 'ulon ulat ulev time')) then
   mylon = ULON(lon_index, lat_index, lev_index)
   mylat = ULAT(lon_index, lat_index, lev_index)
   mylev = ULEV(lon_index, lat_index, lev_index)
elseif ((trim(progvar(varindex)%coordinates) == 'vlon vlat vlev') .or. &
        (trim(progvar(varindex)%coordinates) == 'vlon vlat vlev time')) then
   mylon = VLON(lon_index, lat_index, lev_index)
   mylat = VLAT(lon_index, lat_index, lev_index)
   mylev = VLEV(lon_index, lat_index, lev_index)
elseif ((trim(progvar(varindex)%coordinates) == 'wlon wlat wlev') .or. &
        (trim(progvar(varindex)%coordinates) == 'wlon wlat wlev time')) then
   mylon = WLON(lon_index, lat_index, lev_index)
   mylat = WLAT(lon_index, lat_index, lev_index)
   mylev = WLEV(lon_index, lat_index, lev_index)
else
   write(string1,*) 'unknown coordinate variables of ['//trim(progvar(varindex)%coordinates)//']'
   write(string2,*) 'for variable ',trim(progvar(varindex)%varname)
   call error_handler(E_ERR,'get_state_lonlatlev',string1, &
         source,revision,revdate,text2=string2)
endif

end subroutine get_state_lonlatlev


!-----------------------------------------------------------------------
!>
!> Given an integer index into the state vector structure, returns the kind,
!> and both the starting offset for this kind, as well as the offset into
!> the block of this kind.

subroutine get_state_kind(index_in, var_type, startind, offset)

integer, intent(in)  :: index_in
integer, intent(out) :: var_type
integer, intent(out) :: startind
integer, intent(out) :: offset

integer :: ivar, varindex

VARLOOP : do ivar = 1,nfields
   if ((progvar(ivar)%index1 <= index_in) .and. &
       (progvar(ivar)%indexN >= index_in)) then
      varindex = ivar
      exit VARLOOP
   endif
enddo VARLOOP

var_type  = progvar(varindex)%dart_kind
startind  = progvar(varindex)%index1
offset    = index_in - startind

if (do_output() .and. (debug > 5)) then
    print *, 'get_state_kind: index    = ', index_in
    print *, 'get_state_kind: var_type = ', var_type
    print *, 'get_state_kind: startind = ', startind
    print *, 'get_state_kind: offset   = ', offset
endif

end subroutine get_state_kind


!-----------------------------------------------------------------------
!>
!> intended for debugging use = print out the data min/max for each
!> field in the state vector

subroutine print_ranges(x,varindex)

real(r8),           intent(in) :: x(:)
integer,  optional, intent(in) :: varindex

real(r8) :: minimum, maximum
integer  :: ivar

if (present(varindex)) then

   ivar = varindex

   minimum = minval(x(progvar(ivar)%index1:progvar(ivar)%indexN))
   maximum = maxval(x(progvar(ivar)%index1:progvar(ivar)%indexN))

   write(*,*)trim(progvar(ivar)%varname),' range: ', &
             minimum, maximum, trim(progvar(ivar)%units)

else

   do ivar = 1,nfields

      minimum = minval(x(progvar(ivar)%index1:progvar(ivar)%indexN))
      maximum = maxval(x(progvar(ivar)%index1:progvar(ivar)%indexN))

      write(*,*)trim(progvar(ivar)%varname),' range: ', &
                minimum, maximum, trim(progvar(ivar)%units)
   enddo

endif

end subroutine print_ranges


!-----------------------------------------------------------------------
!>
!> returns TRUE if this point is below the ocean floor or if it is
!> on land.

function is_dry_land(obs_type, lon_index, lat_index, hgt_index)
integer, intent(in) :: obs_type
integer, intent(in) :: lon_index
integer, intent(in) :: lat_index
integer, intent(in) :: hgt_index
logical             :: is_dry_land

logical :: is_ugrid

call error_handler(E_ERR,'is_dry_land','routine not written yet', &
                      source, revision, revdate)

is_dry_land = .FALSE.    ! start out thinking everything is wet.

is_ugrid = is_on_ugrid(obs_type)
if ((      is_ugrid .and. hgt_index > KMU(lon_index, lat_index)) .or. &
    (.not. is_ugrid .and. hgt_index > KMT(lon_index, lat_index))) then
   is_dry_land = .TRUE.
   return
endif

end function is_dry_land


!-----------------------------------------------------------------------
!>
!>  returns true if U, V -- everything else is on T grid

function is_on_ugrid(obs_type)

integer, intent(in) :: obs_type
logical             :: is_on_ugrid

call error_handler(E_ERR,'is_on_ugrid','routine not written yet', &
                      source, revision, revdate)

is_on_ugrid = .FALSE.

if ((obs_type == KIND_U_CURRENT_COMPONENT)  .or.  &
    (obs_type == KIND_V_CURRENT_COMPONENT)) is_on_ugrid = .TRUE.

end function is_on_ugrid


!-----------------------------------------------------------------------
!>
!> Write the grid to an ascii file - in a format suitable for
!> subsequent use in the 'test_interpolation()' code.
!> write_grid_interptest is only possible after reading a real GCOM grid,
!> so static_init_model() must be called to gather the real GCOM grid.

subroutine write_grid_interptest()

integer  :: i, j
real(r8) :: rowmat(nxp1,1), colmat(1,nyp1), dmat(nxp1,nyp1)
real(r8) :: rtlon, rulon, rtlat, rulat, u_val, t_val

!----------------------------------------------------------------------
! Generate a 'Regular' grid with the same rough 'shape' as the dipole grid
!----------------------------------------------------------------------

call error_handler(E_ERR,'write_grid_interptest','routine not written yet', &
                      source, revision, revdate)

open(unit=12, position='rewind', action='write', file='regular_grid_u')
open(unit=13, position='rewind', action='write', file='regular_grid_t')
open(unit=14, position='rewind', action='write', file='regular_grid_u_data')
open(unit=15, position='rewind', action='write', file='regular_grid_t_data')

write(12, *) nxp1, nyp1
write(13, *) nxp1, nyp1

! Have T-grid starting at 0 and U grid offset by half
do i = 1, nxp1
   rtlon = (i - 1.0_r8) * 360.0_r8 / nxp1
   rulon = (i - 0.5_r8) * 360.0_r8 / nxp1
   do j = 1, nyp1
      rtlat = -90.0_r8 + (j - 1.0_r8) * 180.0_r8 / nyp1
      rulat = -90.0_r8 + (j - 0.5_r8) * 180.0_r8 / nyp1
      write(12, *) i, j, rulon, rulat
      write(13, *) i, j, rtlon, rtlat
      ! Now add some wave pattern data
      u_val = sin(3.0_r8*(rulat + 11.0_r8)*2.0_r8*PI/360.0_r8) * &
              sin(4.0_r8*(rulon + 17.0_r8)*2.0_r8*PI/360.0_r8)
      t_val = sin(3.0_r8*(rtlat + 11.0_r8)*2.0_r8*PI/360.0_r8) * &
              sin(4.0_r8*(rtlon + 17.0_r8)*2.0_r8*PI/360.0_r8)
      write(14, *) rulon, rulat, u_val
      write(15, *) rtlon, rtlat, t_val
   enddo
enddo

close(unit=12)
close(unit=13)
close(unit=14)
close(unit=15)

!----------------------------------------------------------------------
! GCOM grid (dipole) next
!----------------------------------------------------------------------

open(unit=12, position='rewind', action='write', file='dipole_grid_u')
open(unit=13, position='rewind', action='write', file='dipole_grid_t')
open(unit=14, position='rewind', action='write', file='dipole_grid_u_data')
open(unit=15, position='rewind', action='write', file='dipole_grid_t_data')

write(12, *) nxp1, nyp1
write(13, *) nxp1, nyp1

rowmat(:,1) = cos(PI * real((/ (i,i=0,nxp1-1) /),r8) / nxp1);
colmat(1,:) = sin(PI * real((/ (i,i=0,nyp1-1) /),r8) / nyp1);
dmat        = matmul(rowmat,colmat)

!do i = 1, nxp1
!   do j = 1, nyp1
!      write(12, *) i, j, ULON(i,j), ULAT(i,j)
!      write(13, *) i, j, TLON(i,j), TLAT(i,j)
!      write(14, *)       ULON(i,j), ULAT(i,j), dmat(i, j)
!      write(15, *)       TLON(i,j), TLAT(i,j), dmat(i, j)
!   enddo
!enddo

close(unit=12)
close(unit=13)
close(unit=14)
close(unit=15)

end subroutine write_grid_interptest


!-----------------------------------------------------------------------
!>
!> use potential temp, depth, and salinity to compute a sensible (in-situ)
!> temperature

subroutine compute_temperature(x, llon, llat, lheight, interp_val, istatus)

real(r8), intent(in)  :: x(:)
real(r8), intent(in)  :: llon
real(r8), intent(in)  :: llat
real(r8), intent(in)  :: lheight
real(r8), intent(out) :: interp_val
integer,  intent(out) :: istatus

integer  :: hstatus, hgt_bot, hgt_top
real(r8) :: hgt_fract, salinity_val, potential_temp, pres_val

call error_handler(E_ERR,'compute_temperature','routine not written yet', &
                      source, revision, revdate)

interp_val = MISSING_R8
istatus = 99

! Get the bounding vertical levels and the fraction between bottom and top
call height_bounds(lheight, nzp1, ZC, hgt_bot, hgt_top, hgt_fract, hstatus)
if(hstatus /= 0) then
   istatus = 12
   return
endif

! salinity - in msu (kg/kg).  converter will want psu (g/kg).
call do_interp(x, start_index(S_index), hgt_bot, hgt_top, hgt_fract, llon, llat, &
               KIND_SALINITY, salinity_val, istatus)
if(istatus /= 0) return
if (do_output() .and. (debug > 8)) print *, 'salinity: ', salinity_val

! potential temperature - degrees C.
call do_interp(x, start_index(T_index), hgt_bot, hgt_top, hgt_fract, llon, llat, &
               KIND_POTENTIAL_TEMPERATURE, potential_temp, istatus)
if(istatus /= 0) return
if (do_output() .and. (debug > 8)) print *, 'potential temp: ', potential_temp

! compute pressure at location between given levels.  these values are in bars;
! the convert routine wants decibars as pressure input.  hgt_fract is 0 at bottom, 1 at top
pres_val = pressure(hgt_bot) + hgt_fract * (pressure(hgt_top) - pressure(hgt_bot))
if (do_output() .and. (debug > 8)) then
   print *, 'local pressure: ', pres_val
   print *, 'bot, top, press: ', hgt_bot, pressure(hgt_bot), &
                                 hgt_top, pressure(hgt_top), pres_val
endif

! and finally, convert to sensible (in-situ) temperature.
! potential temp in degrees C, pressure in decibars, salinity in psu or pss (g/kg).
call insitu_temp(potential_temp, salinity_val*1000.0_r8, pres_val*10.0_r8, interp_val)
if (do_output() .and. (debug > 2)) print *, 's,pt,pres,t: ', salinity_val, potential_temp, pres_val, interp_val

end subroutine compute_temperature


!-----------------------------------------------------------------------
!>
!> do a 2d horizontal interpolation for the value at the bottom level,
!> then again for the top level, then do a linear interpolation in the
!> vertical to get the final value.

subroutine do_interp(x, base_offset, hgt_bot, hgt_top, hgt_fract, &
                     llon, llat, obs_type, interp_val, istatus)

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: base_offset
integer,  intent(in)  :: hgt_bot
integer,  intent(in)  :: hgt_top
real(r8), intent(in)  :: hgt_fract
real(r8), intent(in)  :: llon
real(r8), intent(in)  :: llat
integer,  intent(in)  :: obs_type
real(r8), intent(out) :: interp_val
integer,  intent(out) :: istatus

integer  :: offset
real(r8) :: bot_val, top_val

call error_handler(E_ERR,'do_interp','routine not written yet', &
                      source, revision, revdate)

! Find the base location for the bottom height and interpolate horizontally
!  on this level.  Do bottom first in case it is below the ocean floor; can
!  avoid the second horizontal interpolation.
offset = base_offset + (hgt_bot - 1) * nxp1 * nyp1
if (do_output() .and. (debug > 6)) &
   print *, 'bot, field, abs offset: ', hgt_bot, base_offset, offset

call lon_lat_interpolate(x(offset:), llon, llat, obs_type, hgt_bot, bot_val, istatus)
! Failed istatus from interpolate means give up
if(istatus /= 0) return
if (do_output() .and. (debug > 6)) &
   print *, 'bot_val = ', bot_val

! Find the base location for the top height and interpolate horizontally
!  on this level.
offset = base_offset + (hgt_top - 1) * nxp1 * nyp1
if (do_output() .and. (debug > 6)) &
   print *, 'top, field, abs offset: ', hgt_top, base_offset, offset

call lon_lat_interpolate(x(offset:), llon, llat, obs_type, hgt_top, top_val, istatus)
! Failed istatus from interpolate means give up
if(istatus /= 0) return
if (do_output() .and. (debug > 6)) &
   print *, 'top_val = ', top_val

! Then weight them by the vertical fraction and return
interp_val = bot_val + hgt_fract * (top_val - bot_val)
if (do_output() .and. (debug > 2)) print *, 'do_interp: interp val = ',interp_val

end subroutine do_interp


!-----------------------------------------------------------------------
!>
!> CODE FROM POP MODEL -
!> nsc 1 nov 2012:  i have taken the original subroutine with call:
!>  subroutine dpotmp(press,temp,s,rp,potemp)
!> and removed the original 'press' argument (setting it to 0.0 below)
!> and renamed temp -> potemp, and potemp -> insitu_t
!> i also reordered the args to be a bit more logical.  now you specify:
!> potential temp, salinity, local pressure in decibars, and you get
!> back in-situ temperature (called sensible temperature in the atmosphere;
!> what a thermometer would measure).  the original (F77 fixed format) code
!> had a computed goto which is deprecated/obsolete.  i replaced it with
!> a set of 'if() then else if()' lines.  i did try to not alter the original
!> code so much it wasn't recognizable anymore.
!>
!>  aliciak note: rp = 0 and press = local pressure as function of depth
!>  will return potemp given temp.
!>  the trick here that if you make rp = local pressure and press = 0.0,
!>  and put potemp in the "temp" variable , it will return insitu temp in the
!>  potemp variable.
!>
!> an example figure of the relationship of potential temp and in-situ temp at
!> depth: http://oceanworld.tamu.edu/resources/ocng_textbook/chapter06/chapter06_05.htm
!> see the 'potential temperature' section (note graph starts at -1000m)
!>
!>     title:
!>     *****
!>
!>       insitu_temp  -- calculate sensible (in-situ) temperature from
!>                       local pressure, salinity, and potential temperature
!>
!>     purpose:
!>     *******
!>
!>       to calculate sensible temperature, taken from a converter that
!>       went from sensible/insitu temperature to potential temperature
!>
!>       ref: N.P. Fofonoff
!>            Deep Sea Research
!>            in press Nov 1976
!>
!>     arguments:
!>     **********
!>
!>       potemp     -> potential temperature in celsius degrees
!>       s          -> salinity pss 78
!>       lpres      -> local pressure in decibars
!>       insitu_t   <- in-situ (sensible) temperature (deg c)

subroutine insitu_temp(potemp, s, lpres, insitu_t)

real(r8), intent(in)  :: potemp
real(r8), intent(in)  :: s
real(r8), intent(in)  :: lpres
real(r8), intent(out) :: insitu_t


!     local variables:
!     *********

integer  :: i,j,n
real(r8) :: dp,p,q,r1,r2,r3,r4,r5,s1,t,x

!     code:
!     ****

      s1 = s - 35.0_r8
      p  = 0.0_r8
      t  = potemp

      dp = lpres - p
      n  = int (abs(dp)/1000.0_r8) + 1
      dp = dp/n

      do i=1,n
         do j=1,4

            r1 = ((-2.1687e-16_r8 * t + 1.8676e-14_r8) * t - 4.6206e-13_r8) * p
            r2 = (2.7759e-12_r8*t - 1.1351e-10_r8) * s1
            r3 = ((-5.4481e-14_r8 * t + 8.733e-12_r8) * t - 6.7795e-10_r8) * t
            r4 = (r1 + (r2 + r3 + 1.8741e-8_r8)) * p + (-4.2393e-8_r8 * t+1.8932e-6_r8) * s1
            r5 = r4 + ((6.6228e-10_r8 * t-6.836e-8_r8) * t + 8.5258e-6_r8) * t + 3.5803e-5_r8

            x  = dp*r5

            if (j == 1) then
               t = t + 0.5_r8 * x
               q = x
               p = p + 0.5_r8 * dp

            else if (j == 2) then
               t = t + 0.29298322_r8 * (x-q)
               q = 0.58578644_r8 * x + 0.121320344_r8 * q

            else if (j == 3) then
               t = t + 1.707106781_r8 * (x-q)
               q = 3.414213562_r8*x - 4.121320344_r8*q
               p = p + 0.5_r8*dp

            else ! j must == 4
               t = t + (x - 2.0_r8 * q) / 6.0_r8

            endif

         enddo ! j loop
      enddo ! i loop

      insitu_t = t

if (do_output() .and. (debug > 2)) print *, 'potential temp, salinity, local pressure -> sensible temp'
if (do_output() .and. (debug > 2)) print *, potemp, s, lpres, insitu_t

!       potemp     -> potential temperature in celsius degrees
!       s          -> salinity pss 78
!       lpres      -> local pressure in decibars
!       insitu_t   <- in-situ (sensible) temperature (deg c)

end subroutine insitu_temp


!-----------------------------------------------------------------------
!>
!>  description:
!>  this function computes pressure in bars from depth in meters
!>  using a mean density derived from depth-dependent global
!>  average temperatures and salinities from levitus 1994, and
!>  integrating using hydrostatic balance.
!>
!>  references:
!>
!>  levitus, s., r. burgett, and t.p. boyer, world ocean atlas
!>  volume 3: salinity, noaa atlas nesdis 3, us dept. of commerce, 1994.
!>
!>  levitus, s. and t.p. boyer, world ocean atlas 1994, volume 4:
!>  temperature, noaa atlas nesdis 4, us dept. of commerce, 1994.
!>
!>  dukowicz, j. k., 2000: reduction of pressure and pressure
!>  gradient errors in ocean simulations, j. phys. oceanogr., submitted.
!>
!>  input parameters:
!>  nd     - size of arrays
!>  depth  - depth in meters. no units check is made
!>
!>  output parameters:
!>  pressure - pressure in bars

subroutine depth2pressure(nd, depth, pressure)

integer,  intent(in)  :: nd
real(r8), intent(in)  :: depth(nd)
real(r8), intent(out) :: pressure(nd)
!  local variables & parameters:
integer :: n
real(r8), parameter :: c1 = 1.0_r8

call error_handler(E_ERR,'depth2pressure','routine not written yet', &
                      source, revision, revdate)

! -----------------------------------------------------------------------
!  convert depth in meters to pressure in bars
! -------POP

      do n=1,nd
         pressure(n) = 0.059808_r8*(exp(-0.025_r8*depth(n)) - c1)  &
                     + 0.100766_r8*depth(n) + 2.28405e-7_r8*depth(n)**2
      end do

if (do_output() .and. (debug > 2)) then
   print *, 'depth->pressure conversion table.  cols are: N, depth(m), pressure(bars)'
   do n=1,nd
      print *, n, depth(n), pressure(n)
   enddo
endif

end subroutine depth2pressure


!-----------------------------------------------------------------------
!>
!> find_desired_time_index() returns the index into the time array that
!> matches the target_time.
!> If no target_time is supplied, the model_time is the target.

function find_desired_time_index(ncid, ntimes, time_wanted)
integer,                   intent(in) :: ncid
integer,                   intent(in) :: ntimes
type(time_type), OPTIONAL, intent(in) :: time_wanted
integer                               :: find_desired_time_index

integer  :: VarID
real(r8) :: mytimes(ntimes)
character(len=NF90_MAX_NAME) :: attvalue

type(time_type) :: this_time, target_time
integer :: io, itime, origin_days, origin_seconds
integer :: iday, isecond

find_desired_time_index = MISSING_I   ! initialize to failure setting

if (present(time_wanted)) then
   target_time = time_wanted
else
   target_time = model_time
endif

call nc_check(nf90_inq_varid(ncid, 'time', VarID), &
        'find_desired_time_index:', 'inq_varid time')
call nc_check(nf90_get_var(  ncid, VarID, mytimes), &
        'find_desired_time_index:', 'get_var   time')
call nc_check(nf90_get_att(ncid, VarID, 'units', attvalue), &
        'find_desired_time_index:', 'time get_att units')

this_time = convert_times(attvalue, mytimes)

call get_time(this_time, origin_seconds, origin_days)

! convert each time to a DART time and compare to desired

TIMELOOP : do itime = 1,ntimes

   iday      = int(mytimes(itime))
   isecond   = (mytimes(itime) - iday)*86400
   this_time = set_time(origin_seconds+isecond, origin_days+iday)

   if (this_time == target_time) then
      find_desired_time_index = itime
      exit TIMELOOP
   endif

enddo TIMELOOP

! FIXME ... do we actually need a perfect match ...
! or do we just use the last one if close enough.

! If we did not find one, list all the ones we did find.
if ( find_desired_time_index == MISSING_I ) then
   call print_time(target_time,str='find_desired_time_index:target time is ',iunit=logfileunit)
   call print_time(target_time,str='find_desired_time_index:target time is ')
   call print_date(target_time,str='find_desired_time_index:target date is ',iunit=logfileunit)
   call print_date(target_time,str='find_desired_time_index:target date is ')

   do itime = 1,ntimes
      iday      = int(mytimes(itime))
      isecond   = (mytimes(itime) - iday)*86400
      this_time = set_time(origin_seconds+isecond, origin_days+iday)
      write(string1,'(A,i4,A,i4,A)') 'find_desired_time_index:GCOM time index ',itime,' of ',ntimes,' is'
      call print_time(this_time,str=string1)
      call print_time(this_time,str=string1,iunit=logfileunit)
      call print_date(this_time,str=string1)
      call print_date(this_time,str=string1,iunit=logfileunit)
   enddo

   write(string1,*)'No matching time found'
   call error_handler(E_ERR, 'find_desired_time_index:', string1, &
          source, revision, revdate )
endif

if (do_output() .and. (debug > 0)) then
   call print_time(target_time,str='target time is ',iunit=logfileunit)
   call print_time(target_time,str='target time is ')
   call print_date(target_time,str='target date is ',iunit=logfileunit)
   call print_date(target_time,str='target date is ')
   write(string1,*)'matching time index is ',find_desired_time_index
   call error_handler(E_MSG, 'find_desired_time_index:', string1)
endif

end function find_desired_time_index


!-----------------------------------------------------------------------
!>
!> get_grid_variable(ncid, varname)

subroutine get_grid_variable(ncid, varname)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname

integer  ::  dimIDs(NF90_MAX_VAR_DIMS)
integer  :: dimlens(NF90_MAX_VAR_DIMS)
integer  :: varid, numdims, i

! FIXME check units of everything - particularly the levels.

call nc_check(nf90_inq_varid(ncid, varname, varid), 'get_grid_variable', &
                   'inq_varid ['//trim(varname)//']')
call nc_check(nf90_inquire_variable( ncid, varid, dimids=dimIDs, ndims=numdims), &
                   'get_grid_variable', 'inquire_variable '//trim(varname))

if ( (numdims /= 3)  ) then
   write(string1,*) trim(varname)//' is not a 3D variable as expected.'
   call error_handler(E_ERR,'get_grid_variable',string1,source,revision,revdate)
endif

DimensionLoop : do i = 1,numdims
   write(string1,*) trim(gcom_geometry_file), varname, 'dimension ',i
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
                      'get_grid_variable', string1)
enddo DimensionLoop

if ( varname == 'ulon' ) then
  allocate( ULON(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=ULON), &
               'get_grid_variable', 'get_var '//trim(varname))

elseif ( varname == 'ulat' ) then
  allocate( ULAT(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=ULAT), &
               'get_grid_variable', 'get_var '//trim(varname))

elseif ( varname == 'ulev' ) then
  allocate( ULEV(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=ULEV), &
               'get_grid_variable', 'get_var '//trim(varname))

elseif ( varname == 'vlon' ) then
  allocate( VLON(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=VLON), &
               'get_grid_variable', 'get_var '//trim(varname))

elseif ( varname == 'vlat' ) then
  allocate( VLAT(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=VLAT), &
               'get_grid_variable', 'get_var '//trim(varname))

elseif ( varname == 'vlev' ) then
  allocate( VLEV(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=VLEV), &
               'get_grid_variable', 'get_var '//trim(varname))

elseif ( varname == 'wlon' ) then
  allocate( WLON(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=WLON), &
               'get_grid_variable', 'get_var '//trim(varname))

elseif ( varname == 'wlat' ) then
  allocate( WLAT(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=WLAT), &
               'get_grid_variable', 'get_var '//trim(varname))

elseif ( varname == 'wlev' ) then
  allocate( WLEV(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=WLEV), &
               'get_grid_variable', 'get_var '//trim(varname))

elseif ( varname == 'lon' ) then
  allocate( LON(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=LON), &
               'get_grid_variable', 'get_var '//trim(varname))

elseif ( varname == 'lat' ) then
  allocate( LAT(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=LAT), &
               'get_grid_variable', 'get_var '//trim(varname))

elseif ( varname == 'lev' ) then
  allocate( LEV(dimlens(1), dimlens(2), dimlens(3)) )
  call nc_check(nf90_get_var(ncid, varid, values=LEV), &
               'get_grid_variable', 'get_var '//trim(varname))
endif

if (do_output() .and. (debug > 99)) then
   write(*,*)'get_grid_variable: variable ['//trim(varname)//']'
   write(*,*)'get_grid_variable: dimids ', dimids(1:numdims)
   write(*,*)'get_grid_variable: length ', dimlens(1:numdims)
endif

end subroutine get_grid_variable


!-----------------------------------------------------------------------
!>
!> init_interp() Initializes data structures needed for interpolation.
!>
!> General philosophy is to precompute a 'get_close' structure with
!> a list of what state vector elements are within a certain distance
!> of each other.

subroutine init_interp()

! This should be called at static_init_model time to avoid
! having all this temporary storage in the middle of a run.
!
! DART has a 'get_close' type that divides the domain into a set of boxes
! and tracks which locations are in which box. The get_close_obs() routine
! then finds the box and the subsequent distances. The trick is to set up
! the set of boxes parsimoniously.

integer :: i, j, k, num_neighbors

real(r8) :: max_lon, max_lat, max_lev, maxdist
real(r8) :: lon_dist, lat_dist, lev_dist
real(r8) :: meters_to_radians

! Need to determine the likely maximum size of any of the gridcells.
! the idea is to find a distance that allows us to find the surrounding
! gridcell locations without finding 'too many' neighboring gridcells.
!
! Since this is a regional model, I am going to use the
! lat & lon values as proxies for the physical distance.
! The maximum grid cell size from one grid should suffice.

! FIXME ... this is going to be difficult at the poles because the grid
! could span Greenwich and adjacent longitudes could be HUGELY different.
! solution might be to use xyz location mod a'la mpas_atm.

max_lon = 0.0_r8
max_lat = 0.0_r8
max_lev = 0.0_r8

do k = 1, nz
do j = 1, ny
do i = 2, nx
   lon_dist = abs(LON(i-1,j,k) - LON(i,j,k))   ! degrees
   lat_dist = abs(LAT(i-1,j,k) - LAT(i,j,k))   ! degrees
   lev_dist = abs(LEV(i-1,j,k) - LEV(i,j,k))   ! meters
   if(lon_dist > max_lon) max_lon = lon_dist
   if(lat_dist > max_lat) max_lat = lat_dist
   if(lev_dist > max_lev) max_lev = lev_dist
enddo
enddo
enddo

do k = 1, nz
do j = 2, ny
do i = 1, nx
   lon_dist = abs(LON(i,j-1,k) - LON(i,j,k))
   lat_dist = abs(LAT(i,j-1,k) - LAT(i,j,k))
   lev_dist = abs(LEV(i,j-1,k) - LEV(i,j,k))
   if(lon_dist > max_lon) max_lon = lon_dist
   if(lat_dist > max_lat) max_lat = lat_dist
   if(lev_dist > max_lev) max_lev = lev_dist
enddo
enddo
enddo

do k = 2, nz
do j = 1, ny
do i = 1, nx
   lon_dist = abs(LON(i,j,k-1) - LON(i,j,k))
   lat_dist = abs(LAT(i,j,k-1) - LAT(i,j,k))
   lev_dist = abs(LEV(i,j,k-1) - LEV(i,j,k))
   if(lon_dist > max_lon) max_lon = lon_dist
   if(lat_dist > max_lat) max_lat = lat_dist
   if(lev_dist > max_lev) max_lev = lev_dist
enddo
enddo
enddo

! 2PI radians for the circumference of the earth (PI*earth_diameter_in_m).
! meters_to_radians = 2.0_r8 * PI / (PI * 2.0_r8 * earth_radius * 1000.0_r8)

meters_to_radians = 1.0_r8 / (earth_radius * 1000.0_r8)

max_lon = max_lon / rad2deg
max_lat = max_lat / rad2deg
max_lev = max_lev * meters_to_radians
maxdist = sqrt(max_lon*max_lon + max_lat*max_lat + max_lev*max_lev)

if (do_output() .and. (debug > 1)) then
   write(*,*)
   write(*,*)'init_interp: summary'
   write(*,*)'init_interp: Maximum distance ',maxdist, ' (radians) based on '
   write(*,*)'init_interp: max_lon = ', max_lon, 'radians',max_lon*rad2deg,'degrees'
   write(*,*)'init_interp: max_lat = ', max_lat, 'radians',max_lat*rad2deg,'degrees'
   write(*,*)'init_interp: max_lev = ', max_lev, 'radians',max_lev/meters_to_radians,'meters'
   write(*,*)
   write(*,*)'There are ',1.0_r8/meters_to_radians,' meters in 1 radian (on earth).'
   write(*,*)'This could be the suggested &location_nml:vert_normalization_height.'
   write(*,*)
endif

! maxdist unscaled resulted in 928 candidates for a location in the middle of the domain
num_neighbors = 8   ! as opposed to model_size ! NANCY
call get_close_maxdist_init(gc_state, maxdist*0.8_r8 )
call get_close_obs_init(gc_state, model_size, state_locations)

interpolation_initialized = .true.

end subroutine init_interp


!-----------------------------------------------------------------------
!>
!> set_state_locations_kinds() creates an array the size of the state with
!>                       the location of each of the corresponding elements.

subroutine set_state_locations_kinds()

integer  :: ivar, indx
integer  :: lon_index, lat_index, lev_index
real(r8) :: mylon, mylat, mylev

! FIXME at some point, having the state_locations() array should obviate
! the need to keep the ULON,ULAT,ULEV ...etc. objects around. At present,
! they are needed by nc_write_model_atts() - but that could be about it.

allocate(state_locations(model_size)) ! module data that persists because it
                                      ! is used by model_interpolate
allocate(state_kinds(model_size))     ! module data that persists because it
                                      ! is used by model_interpolate

do ivar = 1,nfields
   do indx = progvar(ivar)%index1,progvar(ivar)%indexN

      call get_state_indices(indx, lon_index, lat_index, lev_index, ivar)

      call get_state_lonlatlev(ivar, lon_index, lat_index, lev_index, &
                               mylon, mylat, mylev)

      state_locations(indx) = set_location(mylon, mylat, mylev, VERTISHEIGHT)
      state_kinds(indx)     = progvar(ivar)%dart_kind
   enddo
enddo

end subroutine set_state_locations_kinds


!-----------------------------------------------------------------------
!>
!> inverse_distance_interpolation()
!>
!> This does an (inverse distance)^2 weighting over the N nearest neighbors.

subroutine inverse_distance_interpolation(x, distances, num_neighbors, obs_val, istatus)

real(r8), intent(in)  :: x(:)
real(r8), intent(in)  :: distances(:)
integer,  intent(in)  :: num_neighbors
real(r8), intent(out) :: obs_val
integer,  intent(out) :: istatus

real(r8) :: inverse_distances(num_neighbors)
real(r8) :: weights(num_neighbors)

obs_val = MISSING_R8

if ( count(distances(2:num_neighbors) <= tiny(0.0_r8)) > 0 ) then
   string1 = 'More than 1 distance is zero. Should not happen'
   call error_handler(E_MSG,'inverse_distance_interpolation', string1)
   istatus = 9
   return
   !stop
endif

! Check if one and only one has a zero distance, just use it.
! If it had exactly a zero distance it would have infinite weight.
! Since the distances are sorted, zero distance would have to be the first one.
if (distances(1) <= tiny(0.0_r8)) then
   obs_val = x(1)
else
   inverse_distances = (1.0_r8 / distances(1:num_neighbors))**2
   weights = inverse_distances / sum(inverse_distances)
   obs_val = sum(x(1:num_neighbors) * weights)
endif

istatus = 0

end subroutine inverse_distance_interpolation


!-----------------------------------------------------------------------
!>
!> horizontal_interpolation
!>

subroutine horizontal_interpolation(location, closest_index, num_wanted, indices, &
   distances, varindex, state, which_layer, interp_value, interp_dist, istatus)

type(location_type), intent(in)    :: location
integer,             intent(in)    :: closest_index
integer,             intent(in)    :: num_wanted
integer,             intent(in)    :: indices(:)
real(r8),            intent(in)    :: distances(:)
integer,             intent(inout) :: varindex
real(r8),            intent(in)    :: state(:)
character(len=*),    intent(in)    :: which_layer
real(r8),            intent(out)   :: interp_value
real(r8),            intent(out)   :: interp_dist
integer,             intent(out)   :: istatus

real(r8) :: lon_lat_lev(3)
real(r8) :: level_of_interest
integer  :: lon_index, lat_index, lev_index, dart_state_index

! Find the layer of interest.
! Find the horizontally adjacent gridpoints.

interp_value = MISSING_R8
istatus = 1

! If the distances are calculated using horizontal information only, we still need
! to find the layer of interest. For our purpose, if the observation falls exactly
! on a layer, it is considered to be on the layer 'below'.

lon_lat_lev = get_location(location)
level_of_interest = lon_lat_lev(3)

call get_state_indices(closest_index, lon_index, lat_index, lev_index, varindex)

dart_state_index = ijk_to_state_index(varindex, lon_index, lat_index, lev_index)

call error_handler(E_ERR,'horizontal_interpolation','FIXME. Not even close to being finished.')

if (do_output() .and. (debug > 99)) &
   write(*,*)'DEBUG ',closest_index, lon_index, lat_index, lev_index, dart_state_index

end subroutine horizontal_interpolation



!-----------------------------------------------------------------------
!>
!> convert_times
!>

function convert_times(timeunits, times)

character(len=*), intent(in)    :: timeunits
real(r8),         intent(inout) :: times(:)
type(time_type)                 :: convert_times

integer :: io, iyear, imonth, iday, ihour, iminute, isecond

! Expecting units of days, but sometimes they are in seconds since
! time:units = "seconds since 2004-01-01 00:00:00" ;
! time:units = "days since 2004-01-01 00:00:00" ;
!               12345678901234

if      (timeunits(1:10) == 'days since') then

   read(timeunits,'(11x,i4,5(1x,i2))',iostat=io)iyear,imonth,iday,ihour,iminute,isecond

else if (timeunits(1:10) == 'seconds si') then ! convert to days

   read(timeunits,'(14x,i4,5(1x,i2))',iostat=io)iyear,imonth,iday,ihour,iminute,isecond
   times = times / 86400.0_r8

else
   write(string1,*)'expecting time units of [days since ... ]'
   write(string2,*)'                or   [seconds since ... ]'
   write(string3,*)'read time units of ['//trim(timeunits)//']'
   call error_handler(E_ERR, 'convert_times:', string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

if (io /= 0) then
   write(string1,*)'Unable to read time units. Error status was ',io
   write(string2,*)'expected "... since YYYY-MM-DD HH:MM:SS"'
   write(string3,*)'was      "'//trim(timeunits)//'"'
   call error_handler(E_ERR, 'convert_times:', string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

convert_times = set_date(iyear, imonth, iday, ihour, iminute, isecond)

end function convert_times


!------------------------------------------------------------------
! End of model_mod
!------------------------------------------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

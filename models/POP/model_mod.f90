! Data Assimilation Research Testbed -- DART
! Copyright 2004-2009, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! This is the interface between the POP ocean model and DART.

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, SECPERDAY, MISSING_R8, rad2deg, PI
use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time, &
                             set_calendar_type, print_time, print_date, &
                             operator(*),  operator(+), operator(-), &
                             operator(>),  operator(<), operator(/), &
                             operator(/=), operator(<=)
use     location_mod, only : location_type,      get_close_maxdist_init, &
                             get_close_obs_init, get_close_obs, set_location, &
                             VERTISHEIGHT, get_location, vert_is_height, &
                             vert_is_level, vert_is_surface
use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, &
                             logfileunit, get_unit, nc_check, do_output, to_upper, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, file_exist, find_textfile_dims, file_to_text
use     obs_kind_mod, only : KIND_TEMPERATURE, KIND_SALINITY, KIND_U_CURRENT_COMPONENT, &
                             KIND_V_CURRENT_COMPONENT, KIND_SEA_SURFACE_HEIGHT
use mpi_utilities_mod, only: my_task_id
use random_seq_mod,   only : random_seq_type, init_random_seq, random_gaussian

use typesizes
use netcdf 

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,         &
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

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.
public :: POP_meta_type, get_gridsize, set_model_end_time, &
          restart_file_to_sv, sv_to_restart_file

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

character(len=256) :: msgstring
logical, save :: module_initialized = .false.

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

!------------------------------------------------------------------
!
! POP namelist section:
!  grid information comes in several files:
!   horizontal grid lat/lons in one, topography (lowest valid vert
!   level) in another, and the vertical grid spacing in a third.
!  if the grid is global, or at least wraps around the sphere in
!   longitude, set longitude_wrap to be true.
!
!------------------------------------------------------------------
! The time manager namelist variables
!------------------------------------------------------------------

character(len=128) :: horiz_grid_input_file = 'no_horiz_grid_input_file'
character(len=128) :: topography_input_file = 'no_topography_input_file'
character(len=128) :: vert_grid_input_file  = 'no_vert_grid_input_file'
logical            :: longitude_wrap        = .true.

!integer :: nblocks   - model_mod doesn't need this
!                       but advance_model will.  not sure this
!                       helps in this code, since the script
!                       needs the info, not the filter executable. 

! here are what we can get from the horiz grid file:
!   real (r8), dimension(:,:), allocatable :: &
!      ULAT,            &! latitude  (radians) of U points
!      ULON,            &! longitude (radians) of U points
!      HTN ,            &! length (cm) of north edge of T box
!      HTE ,            &! length (cm) of east  edge of T box
!      HUS ,            &! length (cm) of south edge of U box
!      HUW ,            &! length (cm) of west  edge of U box
!      ANGLE             ! angle
!
! here is what we can get from the topog file:
!   integer, dimension(:,:), allocatable :: &
!      KMT               ! k index of deepest grid cell on T grid
! maybe
!      KMU               ! k index of deepest grid cell on U grid
!      HT                ! real(r8) value of deepest valid T depth (in cm)
!      HU                ! real(r8) value of deepest valid U depth (in cm)
!
!
! the vert grid file is ascii, with 3 columns/line:
!    cell thickness(in cm)   cell center(in m)   cell bottom(in m)


! other things which can/should be in the model_nml
logical  :: output_state_vector = .true.
integer  :: assimilation_period_days = 1
integer  :: assimilation_period_seconds = 0
real(r8) :: model_perturbation_amplitude = 0.2
integer  :: debug = 0   ! turn up for more and more debug messages

character(len=32):: calendar = 'noleap'

namelist /model_nml/  &
   horiz_grid_input_file,       &
   topography_input_file,       &
   vert_grid_input_file,        &
   longitude_wrap,              &
   output_state_vector,         &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug

!------------------------------------------------------------------
!
! The DART state vector (control vector) will consist of:  S, T, U, V, SHGT
! (Salinity, Temperature, U velocity, V velocity, Sea Surface Height).
! S, T are 3D arrays, located at cell centers.  U,V are at grid cell corners.
! SHGT is a 2D field (X,Y only).  The Z direction is downward.
!
! FIXME: proposed change 1: we put SSH first, then T,U,V, then S, then
!                           any optional tracers, since SSH is the only 2D
!                           field; all tracers are 3D.  this simplifies the
!                           mapping to and from the vars to state vector.
!
! FIXME: proposed change 2: we make this completely namelist driven,
!                           both contents and order of vars.  this should
!                           wait until restart files are in netcdf format,
!                           to avoid problems with incompatible namelist
!                           and IC files.  it also complicates the mapping
!                           to and from the vars to state vector.
!------------------------------------------------------------------

integer, parameter :: n3dfields = 4
integer, parameter :: n2dfields = 1
integer, parameter :: nfields   = n3dfields + n2dfields

! (the absoft compiler likes them to all be the same length during declaration)
! we trim the blanks off before use anyway, so ...
character(len=128) :: progvarnames(nfields) = (/'SALT','TEMP','UVEL','VVEL','SHGT'/)

integer, parameter :: S_index    = 1
integer, parameter :: T_index    = 2
integer, parameter :: U_index    = 3
integer, parameter :: V_index    = 4
integer, parameter :: SHGT_index = 5

integer :: start_index(nfields)

! Grid parameters - the values will be read from a
! standard POP namelist and filled in here.

integer :: Nx=-1, Ny=-1, Nz=-1    ! grid counts for each field

! locations of cell centers (C) and edges (G) for each axis.
real(r8), allocatable :: XC(:), XG(:), YC(:), YG(:), ZC(:), ZG(:)
real(r8), allocatable :: ULAT(:,:), ULON(:,:), TLAT(:,:), TLON(:,:)


! integer, lowest valid cell number in the vertical
integer, allocatable  :: KMT(:, :), KMU(:, :)
! real, depth of lowest valid cell (0 = land).  use only if KMT/KMU not avail.
real(r8), allocatable :: HT(:,:), HU(:,:)

real(r8)        :: endTime
real(r8)        :: ocean_dynamics_timestep = 900.0_r4
integer         :: timestepcount = 0
type(time_type) :: model_time, model_timestep

integer :: model_size    ! the state vector length

! /pkg/mdsio/mdsio_write_meta.F writes the .meta files 
type POP_meta_type
!  private
   integer :: nDims
   integer :: dimList(3)
   character(len=32) :: dataprec
   integer :: reclen
   integer :: nrecords
   integer :: timeStepNumber    ! optional
end type POP_meta_type

INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
END INTERFACE


contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model. In this case,
! it reads in the grid information.

integer :: iunit, io
integer :: ss, dd

! The Plan:
!
!   read in the grid sizes from the horiz grid file and the vert grid file
!   horiz is netcdf, vert is ascii
!  
!   allocate space, and read in actual grid values
!
!   figure out model timestep.  FIXME: from where?
!
!   Compute the model size.
!
!   set the index numbers where the field types change
!

if ( module_initialized ) return ! only need to do this once.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output()) write(     *     , nml=model_nml)

! POP calendar information

call set_calendar_type(calendar)

!---------------------------------------------------------------
! get data dimensions, then allocate space, then open the files
! and actually fill in the arrays.

call get_horiz_grid_dims(Nx, Ny)
call get_vert_grid_dims(Nz)

! Allocate space for grid variables. 
allocate(XC(Nx), YC(Ny), ZC(Nz))
allocate(XG(Nx), YG(Ny), ZG(Nz))

! These will eventually replace the single dim X/Y arrays above
! once we put in interpolation in the dipole grid code.
allocate(ULAT(Nx, Ny), ULON(Nx, Ny), TLAT(Nx, Ny), TLON(Nx, Ny))
allocate(KMT(Nx, Ny), KMU(Nx, Ny))
allocate(HT(Nx, Ny), HU(Nx, Ny))

! Fill them in.
! horiz grid initializes ULAT/LON, TLAT/LON as well.
! kmt initializes HT/HU if present in input file.
call read_horiz_grid(Nx, Ny, XC, XG, YC, YG)
call read_vert_grid(Nz, ZC, ZG)
call read_kmt(Nx, Ny, KMT, KMU)

if (debug > 0) call write_grid_netcdf() ! DEBUG only

!---------------------------------------------------------------
! set the time step from the namelist for now.

!! Define the assimilation period as the model_timestep
!! Ensure model_timestep is multiple of ocean_dynamics_timestep
!
model_timestep = set_model_time_step(assimilation_period_seconds, &
                                     assimilation_period_days,    &
                                     ocean_dynamics_timestep)

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(msgstring,*)"assimilation period is ",dd," days ",ss," seconds"
call error_handler(E_MSG,'static_init_model',msgstring,source,revision,revdate)
if (do_output()) write(logfileunit,*)msgstring

!---------------------------------------------------------------
! compute the offsets into the state vector for the start of each
! different variable type.

! record where in the state vector the data type changes
! from one type to another, by computing the starting
! index for each block of data.
start_index(S_index)    = 1
start_index(T_index)    = start_index(S_index) + (Nx * Ny * Nz)
start_index(U_index)    = start_index(T_index) + (Nx * Ny * Nz)
start_index(V_index)    = start_index(U_index) + (Nx * Ny * Nz)
start_index(SHGT_index) = start_index(V_index) + (Nx * Ny * Nz)

! in spite of the staggering, all grids are the same size
! and offset by half a grid cell.  4 are 3D and 1 is 2D.
!  e.g. S,T,U,V = 256 x 225 x 70
!  e.g. SHGT = 256 x 225

if (do_output()) write(logfileunit, *) 'Using grid size : '
if (do_output()) write(logfileunit, *) '  Nx, Ny, Nz = ', Nx, Ny, Nz
if (do_output()) write(     *     , *) 'Using grid size : '
if (do_output()) write(     *     , *) '  Nx, Ny, Nz = ', Nx, Ny, Nz

model_size = (n3dfields * (Nx * Ny * Nz)) + (n2dfields * (Nx * Ny))
if (do_output()) write(*,*) 'model_size = ', model_size

end subroutine static_init_model




subroutine init_conditions(x)
!------------------------------------------------------------------
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model
 
x = 0.0_r8

end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
!
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to 
! compute a timestep, for instance for radiation computations.
! This interface is only called IF the namelist parameter
! async is set to 0 in perfect_model_obs or filter -OR- if the 
! program integrate_model is to be used to advance the model
! state as a separate executable. If none of these options
! are used (the model will only be advanced as a separate 
! model-specific executable), this can be a NULL INTERFACE.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

if (do_output()) then
   call print_time(time,'NULL interface adv_1step (no advance) DART time is')
   call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
endif

end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
!
! Returns the size of the model as an integer. Required for all
! applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
!
! Companion interface to init_conditions. Returns a time that is somehow 
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time



subroutine model_interpolate(x, location, obs_type, interp_val, istatus)
!=======================================================================
!

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: interp_val
integer,            intent(out) :: istatus

! Model interpolate will interpolate any state variable (S, T, U, V, SHGT) to
! the given location given a state vector. The type of the variable being
! interpolated is obs_type since normally this is used to find the expected
! value of an observation at some location. The interpolated value is 
! returned in interp_val and istatus is 0 for success.

! Local storage
real(r8)       :: loc_array(3), llon, llat, lheight
integer        :: base_offset, offset, ind
integer        :: hgt_bot, hgt_top
real(r8)       :: hgt_fract
real(r8)       :: top_val, bot_val
integer        :: hstatus

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the 
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

interp_val = MISSING_R8     ! the DART bad value flag
istatus = 99                ! unknown error

! Get the individual locations values
loc_array = get_location(location)
llon    = loc_array(1)
llat    = loc_array(2)
lheight = loc_array(3)

if (debug > 1) print *, 'requesting interpolation at ', llon, llat, lheight

if( vert_is_height(location) ) then
   ! Nothing to do 
elseif ( vert_is_surface(location) ) then
   ! Nothing to do 
elseif (vert_is_level(location)) then
   ! convert the level index to an actual depth 
   ind = nint(loc_array(3))
   if ( (ind < 1) .or. (ind > size(zc)) ) then 
      lheight = zc(ind)
   else
      istatus = 11
      return
   endif
else   ! if pressure or undefined, we don't know what to do
   istatus = 17
   return
endif

! Do horizontal interpolations for the appropriate levels
! Find the basic offset of this field
if(obs_type == KIND_SALINITY) then
   base_offset = start_index(1)
else if(obs_type == KIND_TEMPERATURE) then
   base_offset = start_index(2)
else if(obs_type == KIND_U_CURRENT_COMPONENT) then
   base_offset = start_index(3)
else if(obs_type == KIND_V_CURRENT_COMPONENT) then
   base_offset = start_index(4)
else if(obs_type == KIND_SEA_SURFACE_HEIGHT) then
   base_offset = start_index(5)
else
   ! Not a legal type for interpolation, return istatus error
   istatus = 15
   return
endif

if (debug > 1) print *, 'base offset now ', base_offset

! For Sea Surface Height don't need the vertical coordinate
if( vert_is_surface(location) ) then
   call lat_lon_interpolate(x(base_offset:), llon, llat, obs_type, 1, interp_val, istatus)
   return
endif

! Get the bounding vertical levels and the fraction between bottom and top
call height_bounds(lheight, nz, zc, hgt_bot, hgt_top, hgt_fract, hstatus)
if(hstatus /= 0) then
   istatus = 12
   return
endif

! Find the base location for the bottom height and interpolate horizontally 
!  on this level.  Do bottom first in case it is below the ocean floor; can
!  avoid the second horizontal interpolation.
offset = base_offset + (hgt_bot - 1) * nx * ny
if (debug > 1) print *, 'relative bot height offset = ', offset - base_offset
if (debug > 1) print *, 'absolute bot height offset = ', offset
call lat_lon_interpolate(x(offset:), llon, llat, obs_type, hgt_bot, bot_val, istatus)
! Failed istatus from interpolate means give up
if(istatus /= 0) return

! Find the base location for the top height and interpolate horizontally 
!  on this level.
offset = base_offset + (hgt_top - 1) * nx * ny
if (debug > 1) print *, 'relative top height offset = ', offset - base_offset
if (debug > 1) print *, 'absolute top height offset = ', offset
call lat_lon_interpolate(x(offset:), llon, llat, obs_type, hgt_top, top_val, istatus)
! Failed istatus from interpolate means give up
if(istatus /= 0) return


! Then weight them by the vertical fraction and return
interp_val = bot_val + hgt_fract * (top_val - bot_val)
if (debug > 1) print *, 'model_interp: interp val = ',interp_val

! All good.
istatus = 0

end subroutine model_interpolate



subroutine height_bounds(lheight, nheights, hgt_array, bot, top, fract, istatus)
!=======================================================================
!
real(r8),             intent(in) :: lheight
integer,              intent(in) :: nheights
real(r8),             intent(in) :: hgt_array(nheights)
integer,             intent(out) :: bot, top
real(r8),            intent(out) :: fract
integer,             intent(out) :: istatus

! Local variables
integer   :: i

if ( .not. module_initialized ) call static_init_model

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
if (debug > 1) print *, 'above first level in height'
if (debug > 1) print *, 'hgt_array, top, bot, fract=', hgt_array(1), top, bot, fract
endif

! Search through the boxes
do i = 2, nheights
   ! If the location is shallower than this entry, it must be in this box
   if(lheight < hgt_array(i)) then
      top = i -1
      bot = i
      fract = (hgt_array(bot) - lheight) / (hgt_array(bot) - hgt_array(top))
if (debug > 1) print *, 'i, hgt_array, top, bot, fract=', i, hgt_array(i), top, bot, fract
      return
   endif
end do

! Falling off the end means the location is lower than the deepest height
istatus = 20

end subroutine height_bounds


subroutine lat_lon_interpolate(x, llon, llat, var_type, hgt_bot, &
                               interp_val, istatus)
!=======================================================================
!

! Subroutine to interpolate to a lat lon location given that portion of state
! vector. 
! NOTE: Using array sections to pass in the x array may be inefficient on some
! compiler/platform setups. Might want to pass in the entire array with a base
! offset value instead of the section if this is an issue.

real(r8),            intent(in) :: x(:)
real(r8),            intent(in) :: llon, llat
integer,             intent(in) :: var_type
integer,             intent(in) :: hgt_bot
integer,            intent(out) :: istatus
real(r8),           intent(out) :: interp_val

! Local storage
real(r8) :: lat_array(ny), lon_array(nx)
integer  :: lat_bot, lat_top, lon_bot, lon_top
real(r8) :: lat_fract, lon_fract
real(r8) :: pa, pb, pc, pd, xbot, xtop
integer  :: lat_status, lon_status

if ( .not. module_initialized ) call static_init_model

! Succesful return has istatus of 0
! Make any failure here return istatus in the 30s
istatus = 30     ! unknown error

! Find out what latitude box and fraction
! The latitude grid being used depends on the variable type
if(is_on_ugrid(var_type)) then
   ! U and V are on the YG latitude grid
   lat_array = yg
   call lat_bounds(llat, ny, lat_array, lat_bot, lat_top, lat_fract, lat_status)
else 
   ! SHGT, T and S are on the YC latitude grid
   lat_array = yc
   call lat_bounds(llat, ny, lat_array, lat_bot, lat_top, lat_fract, lat_status)
endif

if (debug > 1) print *, 'lat bot, top, fract = ', lat_bot, lat_top, lat_fract
! Check for error on the latitude interpolation
if(lat_status /= 0) then 
   istatus = 31
   return
endif

! Find out what longitude box and fraction
if(is_on_ugrid(var_type)) then
   ! U and V velocity is on the XG grid
   lon_array = xg
   call lon_bounds(llon, nx, lon_array, lon_bot, lon_top, lon_fract, lon_status)
else
   ! SHGT, T, and S are on the XC grid
   lon_array = xc
   call lon_bounds(llon, nx, lon_array, lon_bot, lon_top, lon_fract, lon_status)
endif

if (debug > 1) print *, 'lon bot, top, fract = ', lon_bot, lon_top, lon_fract
! Check for error on the longitude interpolation
if(lat_status /= 0) then 
   istatus = 32
   return
endif

! FIXME: is this the right place to test?  make sure bottom of cell is valid
! DO WE NEED TO TEST ALL 4 lon_bot, lon_top, lat_bot, lat_tops?
if (is_not_ocean(var_type, lon_bot, lat_bot, hgt_bot) .or. &
    is_not_ocean(var_type, lon_top, lat_bot, hgt_bot) .or. &
    is_not_ocean(var_type, lon_bot, lat_top, hgt_bot) .or. &
    is_not_ocean(var_type, lon_top, lat_top, hgt_bot)) then
   istatus = 33
if (debug > 1) print *, 'one corner not ocean'
   return
endif

! Vector is laid out with lat outermost loop, lon innermost loop
! Find the bounding points for the lat lon box
if (debug > 1) print *, 'lon_bot, lon_top = ', lon_bot, lon_top
if (debug > 1) print *, 'lat_bot, lat_top = ', lat_bot, lat_top

pa = get_val(lon_bot, lat_bot, nx, x)
if (debug > 1) print *, 'pa = ', pa
pb = get_val(lon_top, lat_bot, nx, x)
if (debug > 1) print *, 'pb = ', pb
pc = get_val(lon_bot, lat_top, nx, x)
if (debug > 1) print *, 'pc = ', pc
pd = get_val(lon_top, lat_top, nx, x)
if (debug > 1) print *, 'pd = ', pd

if (debug > 1) print *, 'pa,b,c,d = ', pa, pb, pc, pd

! Finish bi-linear interpolation 
! First interpolate in longitude
if (debug > 1) print *, 'bot lon_fract, delta = ', lon_fract, (pb-pa)
xbot = pa + lon_fract * (pb - pa)
if (debug > 1) print *, 'xbot = ', xbot
if (debug > 1) print *, 'top lon_fract, delta = ', lon_fract, (pd-pc)
xtop = pc + lon_fract * (pd - pc)
if (debug > 1) print *, 'xtop = ', xtop
! Now interpolate in latitude
if (debug > 1) print *, 'lat_fract, delta = ', lat_fract, (xtop - xbot)
interp_val = xbot + lat_fract * (xtop - xbot)
if (debug > 1) print *, 'lat_lon_interp: interp_val = ', interp_val

! good return
istatus = 0

end subroutine lat_lon_interpolate




subroutine lat_bounds(llat, nlats, lat_array, bot, top, fract, istatus)

!=======================================================================
!

! Given a latitude llat, the array of latitudes for grid boundaries, and the
! number of latitudes in the grid, returns the indices of the latitude
! below and above the location latitude and the fraction of the distance
! between. istatus is returned as 0 unless the location latitude is 
! south of the southernmost grid point (1 returned) or north of the 
! northernmost (2 returned). If one really had lots of polar obs would 
! want to worry about interpolating around poles.

real(r8),          intent(in) :: llat
integer,           intent(in) :: nlats
real(r8),          intent(in) :: lat_array(nlats)
integer,          intent(out) :: bot, top
real(r8),         intent(out) :: fract
integer,          intent(out) :: istatus

! Local storage
integer    :: i

if ( .not. module_initialized ) call static_init_model

! Success should return 0, failure a positive number.
! Make any failure here return istatus in the 40s
istatus = 40

! Check for too far south or north
if(llat < lat_array(1)) then
   istatus = 41
   return
else if(llat > lat_array(nlats)) then
   istatus = 42
   return
endif

! In the middle, search through
do i = 2, nlats
   if(llat <= lat_array(i)) then
      bot = i - 1
      top = i
      fract = (llat - lat_array(bot)) / (lat_array(top) - lat_array(bot))
if (debug > 1) print *, 'bot, top vals = ', lat_array(i-1), lat_array(i)
if (debug > 1) print *, 'bot, top index = ', bot, top
if (debug > 1) print *, 'fract = ', fract
      istatus = 0
      return
   endif
end do

! shouldn't get here.

end subroutine lat_bounds



subroutine lon_bounds(llon, nlons, lon_array, bot, top, fract, istatus)

!=======================================================================
!

! Given a longitude llon, the array of longitudes for grid boundaries, and the
! number of longitudes in the grid, returns the indices of the longitude
! below and above the location longitude and the fraction of the distance
! between. istatus is returned as 0 unless the location longitude is 
! not between any of the longitude box boundaries. A module storage flag
! indicates whether the grid wraps around in longitude. 
! Algorithm fails for a silly grid that
! has only two longitudes separated by 180 degrees.

real(r8),          intent(in) :: llon
integer,           intent(in) :: nlons
real(r8),          intent(in) :: lon_array(nlons)
integer,          intent(out) :: bot, top
real(r8),         intent(out) :: fract
integer,          intent(out) :: istatus

! Local storage
integer  :: i
real(r8) :: dist_bot, dist_top

if ( .not. module_initialized ) call static_init_model

! Success should return 0, failure a positive number.
! Make any failure here return istatus in the 50s
istatus = 50

if (debug > 1) print *, 'computing bounds for = ', llon
! This is inefficient, someone could clean it up for regularly spaced longitudes
do i = 2, nlons
   dist_bot = lon_dist(llon, lon_array(i - 1))
   dist_top = lon_dist(llon, lon_array(i))
if (debug > 3) print *, 'dist top, bot = ', dist_top, dist_bot
   if(dist_bot <= 0 .and. dist_top > 0) then
      bot = i - 1
      top = i
if (debug > 1) print *, 'bot, top vals = ', lon_array(i-1), lon_array(i)
if (debug > 1) print *, 'bot, top index = ', bot, top
if (debug > 1) print *, 'numerator = ',  dist_bot
if (debug > 1) print *, 'denominator = ', dist_bot + abs(dist_top)
      fract = abs(dist_bot) / (abs(dist_bot) + dist_top)
if (debug > 1) print *, 'fract = ', fract
      istatus = 0
      return
   endif
end do

! Falling off the end means it's in between. Check for wraparound.
if(longitude_wrap) then
   bot = nlons
   top = 1
   dist_bot = lon_dist(llon, lon_array(bot))
   dist_top = lon_dist(llon, lon_array(top)) 
   fract = abs(dist_bot) / (abs(dist_bot) + dist_top)
   istatus = 0
   return
else
   ! Not in the localized longitude grid; return istatus 1
   istatus = 51
endif

end subroutine lon_bounds



function lon_dist(lon1, lon2)
!=======================================================================
!

! Returns the smallest signed distance between lon1 and lon2 on the sphere
! If lon1 is less than 180 degrees east of lon2 the distance is negative
! If lon1 is less than 180 degrees west of lon2 the distance is positive

real(r8), intent(in) :: lon1, lon2
real(r8)             :: lon_dist

if ( .not. module_initialized ) call static_init_model

lon_dist = lon2 - lon1
if(lon_dist >= -180.0_r8 .and. lon_dist <= 180.0_r8) then
   return
else if(lon_dist < -180.0_r8) then
   lon_dist = lon_dist + 360.0_r8
else
   lon_dist = lon_dist - 360.0_r8
endif

end function lon_dist


function get_val(lon_index, lat_index, nlon, x)
!=======================================================================
!

! Returns the value from a single level array given the lat and lon indices
integer,     intent(in) :: lon_index, lat_index, nlon
real(r8),    intent(in) :: x(:)
real(r8)                :: get_val

if ( .not. module_initialized ) call static_init_model

! Layout has lons varying most rapidly
if (debug > 1) print *, 'lat_index, lon_index, nlon', lat_index, lon_index, nlon
if (debug > 1) print *, 'computing index val ', (lat_index - 1) * nlon + lon_index
get_val = x((lat_index - 1) * nlon + lon_index)
if (debug > 1) print *, 'get_val = ', get_val

end function get_val



function set_model_time_step(ss, dd, dt)
!------------------------------------------------------------------
!
! Sets the model 'timestep' AKA 'assimilation period'.
! Must make sure the assimilation period is a multiple of the 
! model's dynamical timestepping requirement.

integer,  intent(in) :: ss    ! assimilation_period_seconds
integer,  intent(in) :: dd    ! assimilation_period_days
real(r8), intent(in) :: dt    ! ocean_dynamics_timestep

type(time_type) :: set_model_time_step

integer :: assim_period, ndt

if ( .not. module_initialized ) call static_init_model

assim_period = ss + dd*SECPERDAY   ! in seconds 
ndt = max(nint(assim_period / dt),1)
assim_period = nint(ndt * dt)

set_model_time_step = set_time(assim_period, 0) ! works seconds > 86400

end function set_model_time_step



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step



subroutine set_model_end_time(adv_to_offset)
!------------------------------------------------------------------
!
! sets PARM03:endTime to reflect the time when the model will stop. 
! endTime is in module storage 

type(time_type), intent(in) :: adv_to_offset

integer :: secs, days

if ( .not. module_initialized ) call static_init_model

call get_time(adv_to_offset, secs, days)

endTime = (secs + days*SECPERDAY)

end subroutine set_model_end_time



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

real(R8) :: lat, lon, depth
integer :: var_num, offset, lon_index, lat_index, depth_index, local_var

if ( .not. module_initialized ) call static_init_model

if (debug > 5) print *, 'asking for meta data about index ', index_in

if (index_in < start_index(S_index+1)) then
   local_var = KIND_SALINITY  
   var_num = S_index
else if (index_in < start_index(T_index+1)) then
   local_var = KIND_TEMPERATURE  
   var_num = T_index
else if (index_in < start_index(U_index+1)) then
   local_var = KIND_U_CURRENT_COMPONENT
   var_num = U_index
else if (index_in < start_index(V_index+1)) then
   local_var = KIND_V_CURRENT_COMPONENT
   var_num = V_index
else 
   local_var = KIND_SEA_SURFACE_HEIGHT
   var_num = SHGT_index
endif
if (present(var_type)) var_type = local_var  

if (debug > 5) print *, 'var num = ', var_num
if (debug > 5) print *, 'var type = ', local_var

! local offset into this var array
offset = index_in - start_index(var_num)

if (debug > 5) print *, 'offset = ', offset

if (var_num == SHGT_index) then
  depth = 0.0
  depth_index = 1
else
  depth_index = (offset / (Nx * Ny)) + 1
  depth = ZC(depth_index)
endif

lat_index = (offset - ((depth_index-1)*Nx*Ny)) / Nx + 1
lon_index = offset - ((depth_index-1)*Nx*Ny) - ((lat_index-1)*Nx) + 1

if (debug > 5) print *, 'lon, lat, depth index = ', lon_index, lat_index, depth_index
if (is_on_ugrid(local_var)) then
   lon = XG(lon_index)
   lat = YG(lat_index)
else
   lon = XC(lon_index)
   lat = YC(lat_index)
endif

if (debug > 5) print *, 'lon, lat, depth = ', lon, lat, depth

location = set_location(lon, lat, depth, VERTISHEIGHT)

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

! if ( .not. module_initialized ) call static_init_model

deallocate(XC, YC, ZC, XG, YG, ZG, KMT, KMU, HT, HU)

end subroutine end_model



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! TJH -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector.
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode 
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

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
integer :: NlonDimID, NlatDimID, NzDimID
integer :: ulonVarID, ulatVarID, tlonVarID, tlatVarID, ZGVarID, ZCVarID
integer :: KMTVarID, KMUVarID

! for the prognostic variables
integer :: SVarID, TVarID, UVarID, VVarID, SHGTVarID 

!----------------------------------------------------------------------
! variables for the namelist output
!----------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_pop_namelist

!----------------------------------------------------------------------
! local variables 
!----------------------------------------------------------------------

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer :: i
character(len=128)  :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename, '(a, i3)') 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                   "nc_write_model_atts", "inquire "//trim(filename))
call nc_check(nf90_Redef(ncFileID),"nc_write_model_atts",   "redef "//trim(filename))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="NMLlinelen", dimid=LineLenDimID), &
                           "nc_write_model_atts","inq_dimid NMLlinelen")
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID), &
                           "nc_write_model_atts", "copy dimid "//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID), &
                           "nc_write_model_atts", "time dimid "//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(msgstring,*)"Time Dimension ID ",TimeDimID, &
             " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", msgstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable", len=model_size, &
        dimid = StateVarDimID),"nc_write_model_atts", "state def_dim "//trim(filename))

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date" ,str1    ), &
           "nc_write_model_atts", "creation put "//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source"  ,source  ), &
           "nc_write_model_atts", "source put "//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision), &
           "nc_write_model_atts", "revision put "//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate" ,revdate ), &
           "nc_write_model_atts", "revdate put "//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model",  "POP" ), &
           "nc_write_model_atts", "model put "//trim(filename))

!-------------------------------------------------------------------------------
! Determine shape of most important namelist
!-------------------------------------------------------------------------------

call find_textfile_dims("pop_in", nlines, linelen)
if (nlines > 0) then
  has_pop_namelist = .true.
else
  has_pop_namelist = .false.
endif

if (debug > 0)    print *, 'pop namelist: nlines, linelen = ', nlines, linelen
  
if (has_pop_namelist) then 
   allocate(textblock(nlines))
   textblock = ''
   
   call nc_check(nf90_def_dim(ncid=ncFileID, name="nlines", &
                 len = nlines, dimid = nlinesDimID), &
                 'nc_write_model_atts', 'def_dim nlines ')
   
   call nc_check(nf90_def_var(ncFileID,name="pop_in", xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'nc_write_model_atts', 'def_var pop_in')
   call nc_check(nf90_put_att(ncFileID, nmlVarID, "long_name",       &
                 "contents of pop_in namelist"), 'nc_write_model_atts', 'put_att pop_in')

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
   call nc_check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
                 dimids=StateVarDimID, varid=StateVarVarID), "nc_write_model_atts", &
                 "statevariable def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,"long_name","State Variable ID"),&
                 "nc_write_model_atts","statevariable long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "units","indexical"), &
                 "nc_write_model_atts", "statevariable units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,"valid_range",(/ 1,model_size /)),&
                 "nc_write_model_atts", "statevariable valid_range "//trim(filename))

   ! Define the actual (3D) state vector, which gets filled as time goes on ... 
   call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
                 dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                 "nc_write_model_atts","state def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,"long_name","model state or fcopy"),&
                 "nc_write_model_atts", "state long_name "//trim(filename))

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID),"nc_write_model_atts","state enddef "//trim(filename))

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                 "nc_write_model_atts", "state put_var "//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to output the prognostic variables.
   !----------------------------------------------------------------------------
   ! Define the new dimensions IDs
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_def_dim(ncid=ncFileID, name="i", &
          len = Nx, dimid = NlonDimID),"nc_write_model_atts", "i def_dim "//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name="j", &
          len = Ny, dimid = NlatDimID),"nc_write_model_atts", "j def_dim "//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name="k", &
          len = Nz, dimid =   NzDimID),"nc_write_model_atts", "k def_dim "//trim(filename))
   
   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncFileID,name="POPnml", xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'nc_write_model_atts', 'def_var POPnml')
   call nc_check(nf90_put_att(ncFileID, nmlVarID, "long_name",       &
                 "namelist.input contents"), 'nc_write_model_atts', 'put_att POPnml')

   ! U,V Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name="ULON", xtype=nf90_real, &
                 dimids=(/ NlonDimID, NlatDimID /), varid=ulonVarID),&
                 "nc_write_model_atts", "ULON def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, "long_name", "longitudes of U,V grid"), &
                 "nc_write_model_atts", "ULON long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, "cartesian_axis", "X"),  &
                 "nc_write_model_atts", "ULON cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, "units", "degrees_east"), &
                 "nc_write_model_atts", "ULON units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)), &
                 "nc_write_model_atts", "ULON valid_range "//trim(filename))

   ! U,V Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name="ULAT", xtype=nf90_real, &
                 dimids=(/ NlonDimID, NlatDimID /), varid=ulatVarID),&
                 "nc_write_model_atts", "ULAT def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, "long_name", "latitudes of U,V grid"), &
                 "nc_write_model_atts", "ULAT long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, "cartesian_axis", "Y"),   &
                 "nc_write_model_atts", "ULAT cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, "units", "degrees_north"),  &
                 "nc_write_model_atts", "ULAT units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID,"valid_range",(/ -90.0_r8, 90.0_r8 /)), &
                 "nc_write_model_atts", "ULAT valid_range "//trim(filename))

   ! S,T,SHGT Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name="TLON", xtype=nf90_real, &
                 dimids=(/ NlonDimID, NlatDimID /), varid=tlonVarID),&
                 "nc_write_model_atts", "TLON def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlonVarID, "long_name", "longitudes of S,T,... grid"), &
                 "nc_write_model_atts", "TLON long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlonVarID, "cartesian_axis", "X"),   &
                 "nc_write_model_atts", "TLON cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlonVarID, "units", "degrees_east"),  &
                 "nc_write_model_atts", "TLON units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlonVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)), &
                 "nc_write_model_atts", "TLON valid_range "//trim(filename))


   ! S,T,SHGT Grid (center) Latitudes
   call nc_check(nf90_def_var(ncFileID,name="TLAT", xtype=nf90_real, &
                 dimids= (/ NlonDimID, NlatDimID /), varid=tlatVarID), &
                 "nc_write_model_atts", "TLAT def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlatVarID, "long_name", "latitudes of S,T, ... grid"), &
                 "nc_write_model_atts", "TLAT long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlatVarID, "cartesian_axis", "Y"),   &
                 "nc_write_model_atts", "TLAT cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlatVarID, "units", "degrees_north"),  &
                 "nc_write_model_atts", "TLAT units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlatVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)), &
                 "nc_write_model_atts", "TLAT valid_range "//trim(filename))

   ! Depths
   call nc_check(nf90_def_var(ncFileID,name="ZG", xtype=nf90_real, &
                 dimids=NzDimID, varid= ZGVarID), &
                 "nc_write_model_atts", "ZG def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, "long_name", "depth at grid edges"), &
                 "nc_write_model_atts", "ZG long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, "cartesian_axis", "Z"),   &
                 "nc_write_model_atts", "ZG cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, "units", "meters"),  &
                 "nc_write_model_atts", "ZG units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, "positive", "down"),  &
                 "nc_write_model_atts", "ZG units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, "comment", &
                  "more positive is closer to the center of the earth"),  &
                 "nc_write_model_atts", "ZG comment "//trim(filename))

   ! Depths
   call nc_check(nf90_def_var(ncFileID,name="ZC",xtype=nf90_real,dimids=NzDimID,varid=ZCVarID), &
                 "nc_write_model_atts", "ZC def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, "long_name", "depth at grid centroids"), &
                 "nc_write_model_atts", "ZC long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, "cartesian_axis", "Z"),   &
                 "nc_write_model_atts", "ZC cartesian_axis "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, "units", "meters"),  &
                 "nc_write_model_atts", "ZC units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, "positive", "down"),  &
                 "nc_write_model_atts", "ZC units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, "comment", &
                  "more positive is closer to the center of the earth"),  &
                 "nc_write_model_atts", "ZC comment "//trim(filename))

   ! Depth mask
   call nc_check(nf90_def_var(ncFileID,name="KMT",xtype=nf90_int, &
                 dimids= (/ NlonDimID, NlatDimID /), varid=KMTVarID), &
                 "nc_write_model_atts", "KMT def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMTVarID, "long_name", "lowest valid depth index at grid centroids"), &
                 "nc_write_model_atts", "KMT long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMTVarID, "units", "levels"),  &
                 "nc_write_model_atts", "KMT units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMTVarID, "positive", "down"),  &
                 "nc_write_model_atts", "KMT units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMTVarID, "comment", &
                  "more positive is closer to the center of the earth"),  &
                 "nc_write_model_atts", "KMT comment "//trim(filename))

   ! Depth mask
   call nc_check(nf90_def_var(ncFileID,name="KMU",xtype=nf90_int, &
                 dimids= (/ NlonDimID, NlatDimID /), varid=KMUVarID), &
                 "nc_write_model_atts", "KMU def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMUVarID, "long_name", "lowest valid depth index at grid corners"), &
                 "nc_write_model_atts", "KMU long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMUVarID, "units", "levels"),  &
                 "nc_write_model_atts", "KMU units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMUVarID, "positive", "down"),  &
                 "nc_write_model_atts", "KMU units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMUVarID, "comment", &
                  "more positive is closer to the center of the earth"),  &
                 "nc_write_model_atts", "KMU comment "//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------


   call nc_check(nf90_def_var(ncid=ncFileID, name="SALT", xtype=nf90_real, &
         dimids = (/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=SVarID),&
         "nc_write_model_atts", "S def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SVarID, "long_name", "salinity"), &
         "nc_write_model_atts", "S long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SVarID, "units", "g/g"), &
         "nc_write_model_atts", "S units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SVarID, "missing_value", NF90_FILL_REAL), &
         "nc_write_model_atts", "S missing "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SVarID, "_FillValue", NF90_FILL_REAL), &
         "nc_write_model_atts", "S fill "//trim(filename))


   call nc_check(nf90_def_var(ncid=ncFileID, name="TEMP", xtype=nf90_real, &
         dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=TVarID),&
         "nc_write_model_atts", "T def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "long_name", "Potential Temperature"), &
         "nc_write_model_atts", "T long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "units", "deg C"), &
         "nc_write_model_atts", "T units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "units_long_name", "degrees celsius"), &
         "nc_write_model_atts", "T units_long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "missing_value", NF90_FILL_REAL), &
         "nc_write_model_atts", "T missing "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, TVarID, "_FillValue", NF90_FILL_REAL), &
         "nc_write_model_atts", "T fill "//trim(filename))


   call nc_check(nf90_def_var(ncid=ncFileID, name="UVEL", xtype=nf90_real, &
         dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=UVarID),&
         "nc_write_model_atts", "U def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "long_name", "U velocity"), &
         "nc_write_model_atts", "U long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "units", "cm/s"), &
         "nc_write_model_atts", "U units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "units_long_name", "centimeters per second"), &
         "nc_write_model_atts", "U units_long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "missing_value", NF90_FILL_REAL), &
         "nc_write_model_atts", "U missing "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, UVarID, "_FillValue", NF90_FILL_REAL), &
         "nc_write_model_atts", "U fill "//trim(filename))


   call nc_check(nf90_def_var(ncid=ncFileID, name="VVEL", xtype=nf90_real, &
         dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=VVarID),&
         "nc_write_model_atts", "V def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "long_name", "V Velocity"), &
         "nc_write_model_atts", "V long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "units", "cm/s"), &
         "nc_write_model_atts", "V units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "units_long_name", "centimeters per second"), &
         "nc_write_model_atts", "V units_long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "missing_value", NF90_FILL_REAL), &
         "nc_write_model_atts", "V missing "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, VVarID, "_FillValue", NF90_FILL_REAL), &
         "nc_write_model_atts", "V fill "//trim(filename))


   call nc_check(nf90_def_var(ncid=ncFileID, name="SHGT", xtype=nf90_real, &
         dimids=(/NlonDimID,NlatDimID,MemberDimID,unlimitedDimID/),varid=SHGTVarID), &
         "nc_write_model_atts", "SHGT def_var "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SHGTVarID, "long_name", "Sea surface height"), &
         "nc_write_model_atts", "SHGT long_name "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SHGTVarID, "units", "cm"), &
         "nc_write_model_atts", "SHGT units "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SHGTVarID, "missing_value", NF90_FILL_REAL), &
         "nc_write_model_atts", "SHGT missing "//trim(filename))
   call nc_check(nf90_put_att(ncFileID, SHGTVarID, "_FillValue", NF90_FILL_REAL), &
         "nc_write_model_atts", "SHGT fill "//trim(filename))

   ! Finished with dimension/variable definitions, must end 'define' mode to fill.

   call nc_check(nf90_enddef(ncfileID), "prognostic enddef "//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables
   !----------------------------------------------------------------------------

   call nc_check(nf90_put_var(ncFileID, ulonVarID, ULON ), &
                "nc_write_model_atts", "ULON put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ulatVarID, ULAT ), &
                "nc_write_model_atts", "ULAT put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, tlonVarID, TLON ), &
                "nc_write_model_atts", "TLON put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, tlatVarID, TLAT ), &
                "nc_write_model_atts", "TLAT put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ZGVarID, ZG ), &
                "nc_write_model_atts", "ZG put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ZCVarID, ZC ), &
                "nc_write_model_atts", "ZC put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, KMTVarID, KMT ), &
                "nc_write_model_atts", "KMT put_var "//trim(filename))
   call nc_check(nf90_put_var(ncFileID, KMUVarID, KMU ), &
                "nc_write_model_atts", "KMU put_var "//trim(filename))

endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

if (has_pop_namelist) then
   call file_to_text("pop_in", textblock)
   call nc_check(nf90_put_var(ncFileID, nmlVarID, textblock ), &
                 'nc_write_model_atts', 'put_var nmlVarID')
   deallocate(textblock)
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), "nc_write_model_atts", "atts sync")

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! TJH 24 Oct 2006 -- Writes the model variables to a netCDF file.
!
! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: VarID

real(r4), dimension(Nx,Ny,Nz) :: data_3d
real(r4), dimension(Nx,Ny)    :: data_2d
character(len=128)  :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename, '(a, i3)') 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
              "nc_write_model_vars", "inquire "//trim(filename))

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, "state", VarID), &
                 "nc_write_model_vars", "state inq_varid "//trim(filename))
   call nc_check(NF90_put_var(ncFileID,VarID,statevec,start=(/1,copyindex,timeindex/)),&
                 "nc_write_model_vars", "state put_var "//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   ! Replace missing values (0.0) with netcdf missing value.
   ! Staggered grid causes some logistical problems.
   ! Hopefully, the conversion between r8 and r4 still preserves 'hard' zeros.
   !----------------------------------------------------------------------------

   call vector_to_prog_var(statevec,S_index,data_3d)
   where (data_3d == 0.0_r4) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, "SALT", VarID), &
                "nc_write_model_vars", "S inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "S put_var "//trim(filename))

   call vector_to_prog_var(statevec,T_index,data_3d)
   where (data_3d == 0.0_r4) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, "TEMP", VarID), &
                "nc_write_model_vars", "T inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "T put_var "//trim(filename))

   call vector_to_prog_var(statevec,U_index,data_3d)
   where (data_3d == 0.0_r4) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, "UVEL", VarID), &
                "nc_write_model_vars", "U inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "U put_var "//trim(filename))

   call vector_to_prog_var(statevec,V_index,data_3d)
   where (data_3d == 0.0_r4) data_3d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, "VVEL", VarID), &
                "nc_write_model_vars", "V inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "V put_var "//trim(filename))

   call vector_to_prog_var(statevec,SHGT_index,data_2d)
   where (data_2d == 0.0_r4) data_2d = NF90_FILL_REAL
   call nc_check(NF90_inq_varid(ncFileID, "SHGT", VarID), &
                "nc_write_model_vars", "SHGT inq_varid "//trim(filename))
   call nc_check(nf90_put_var(ncFileID,VarID,data_2d,start=(/1,1,copyindex,timeindex/)),&
                "nc_write_model_vars", "SHGT put_var "//trim(filename))

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), "nc_write_model_vars", "sync "//trim(filename))

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars



subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
!
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding a perturbation to each model state 
! variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

integer :: i
logical, save :: random_seq_init = .false.

if ( .not. module_initialized ) call static_init_model

interf_provided = .true.

! Initialize my random number sequence
if(.not. random_seq_init) then
   call init_random_seq(random_seq, my_task_id())
   random_seq_init = .true.
endif

! only perturb the non-zero values.  0 is a flag for missing
! ocean cells (e.g. land or under the sea floor)
do i=1,size(state)
   if (state(i) /= 0.0_r8) &
      pert_state(i) = random_gaussian(random_seq, state(i), &
                                      model_perturbation_amplitude)
enddo


end subroutine pert_model_state




subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! If needed by the model interface, this is the current mean
! for all state vector items across all ensembles. It is up to this
! code to allocate space and save a copy if it is going to be used
! later on.  For now, we are ignoring it.

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model

end subroutine ens_mean_for_model




subroutine restart_file_to_sv(filename, state_vector, model_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a POP restart
! file and packs them into a dart state vector.

character(len=*), intent(in)    :: filename 
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: model_time

! temp space to hold data while we are reading it
real(r8) :: data_2d_array(Nx,Ny), data_3d_array(Nx,Ny,Nz)
integer  :: i, j, k, l, indx

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname 
integer :: VarID, numdims, dimlen
integer :: ncid, iyear, imonth, iday, ihour, iminute, isecond, nc_rc
character(len=256) :: myerrorstring 
logical :: convert_to_ssh

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! Check that the input file exists ... 
! Read the time data. 

if ( .not. file_exist(filename) ) then
   write(msgstring,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  "restart_file_to_sv", "open "//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iyear'  , iyear), &
                  "restart_file_to_sv", "get_att iyear")
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'imonth' , imonth), &
                  "restart_file_to_sv", "get_att imonth")
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iday'   , iday), &
                  "restart_file_to_sv", "get_att iday")
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'ihour'  , ihour), &
                  "restart_file_to_sv", "get_att ihour")
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iminute', iminute), &
                  "restart_file_to_sv", "get_att iminute")
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'isecond', isecond), &
                  "restart_file_to_sv", "get_att isecond")

! FIXME: we don't allow a real year of 0 - add one for now, but
! THIS MUST BE FIXED IN ANOTHER WAY!
if (iyear == 0) then
  call error_handler(E_MSG, 'restart_file_to_sv', &
                     'WARNING!!!   year 0 not supported; setting to year 1')
  iyear = 1
endif
model_time = set_date(iyear, imonth, iday, ihour, iminute, isecond)

if (do_output()) &
    call print_time(model_time,'time for restart file '//trim(filename))
if (do_output()) &
    call print_date(model_time,'date for restart file '//trim(filename))

! Start counting and filling the state vector one item at a time,
! repacking the 3d arrays into a single 1d list of numbers.
! These must be a fixed number and in a fixed order.

indx = 1

! fill SALT, TEMP, UVEL, VVEL in that order
! The POP restart files have two time steps for each variable,
! the variables are named SALT_CUR and SALT_OLD ... for example.
! We are only interested in the CURrent time step.

do l=1, n3dfields

   varname = trim(progvarnames(l))//'_CUR'
   myerrorstring = trim(filename)//' '//trim(varname)

   ! Is the netCDF variable the right shape?

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            "restart_file_to_sv", "inq_varid "//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarId,dimids=dimIDs,ndims=numdims), &
            "restart_file_to_sv", "inquire "//trim(myerrorstring))

   if (numdims /= 3) then
      write(msgstring,*) trim(myerrorstring),' does not have exactly 3 dimensions'
      call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
   endif

   do i = 1,numdims
      write(msgstring,'(''inquire dimension'',i,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            "restart_file_to_sv", msgstring)

      if (dimlen /= size(data_3d_array,i)) then
         write(msgstring,*) trim(myerrorstring),'dim/dimlen',i,dimlen,'not',size(data_3d_array,i)
         call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
      endif
   enddo   

   ! Actually get the variable and stuff it into the array

   call nc_check(nf90_get_var(ncid, VarID, data_3d_array), "restart_file_to_sv", &
                "get_var "//trim(varname))

   do k = 1, Nz   ! size(data_3d_array,3)
   do j = 1, Ny   ! size(data_3d_array,2)
   do i = 1, Nx   ! size(data_3d_array,1)
      state_vector(indx) = data_3d_array(i, j, k)
      indx = indx + 1
   enddo
   enddo
   enddo

enddo

! and finally, SHGT (and any other 2d fields)
do l=(n3dfields+1), (n3dfields+n2dfields)

   select case ( trim(progvarnames(l)) )
   case ("SHGT")
      ! The restart files do not drag around SHGT, but they do drag
      ! around PSURF ... which can be used to calculate SHGT.
      varname = 'PSURF_CUR'
      convert_to_ssh = .true.
   case default
      varname = trim(progvarnames(l))//'_CUR'
      convert_to_ssh = .false.
   end select

   myerrorstring = trim(varname)//' '//trim(filename)

   ! Is the netCDF variable the right shape?

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            "restart_file_to_sv", "inq_varid "//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarId,dimids=dimIDs,ndims=numdims), &
            "restart_file_to_sv", "inquire "//trim(myerrorstring))

   if (numdims /= 2) then
      write(msgstring,*) trim(myerrorstring),' does not have exactly 2 dimensions'
      call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
   endif

   do i = 1,numdims
      write(msgstring,'(''inquire dimension'',i,A)')i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            "restart_file_to_sv", msgstring)

      if (dimlen /= size(data_2d_array,i)) then
         write(msgstring,*) trim(myerrorstring),'dim/dimlen',i,dimlen,'not',size(data_2d_array,i)
         call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
      endif
   enddo   

   ! Actually get the variable and stuff it into the array

   call nc_check(nf90_get_var(ncid, VarID, data_2d_array), "restart_file_to_sv", &
                "get_var "//trim(varname))

   if ( convert_to_ssh ) then  ! SSH=psurf/980.6    POP uses CGS 
      data_2d_array = data_2d_array/980.6_r8
   endif

   do j = 1, Ny   ! size(data_3d_array,2)
   do i = 1, Nx   ! size(data_3d_array,1)
      state_vector(indx) = data_2d_array(i, j)
      indx = indx + 1
   enddo
   enddo

enddo

end subroutine restart_file_to_sv



subroutine sv_to_restart_file(state_vector, filename, date1, date2)
!------------------------------------------------------------------
!
real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename 
type(time_type),  intent(in) :: date1, date2 

integer :: iyear, imonth, iday, ihour, iminute, isecond
type(time_type) :: pop_time

! temp space to hold data while we are writing it
real(r4) :: data_2d_array(Nx,Ny), data_3d_array(Nx,Ny,Nz)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname 
character(len=256)                    :: myerrorstring 

integer :: i, l, ncid, VarID, numdims, dimlen
logical :: convert_from_ssh

!----------------------------------------------------------------------
! Get the show underway
!----------------------------------------------------------------------

if ( .not. module_initialized ) call static_init_model

! Check that the input file exists. 
! make sure the time tag in the restart file matches 
! the current time of the DART state ...

if ( .not. file_exist(filename)) then
   write(msgstring,*)trim(filename),' does not exist. FATAL error.'
   call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate) 
endif

call nc_check( nf90_open(trim(filename), NF90_WRITE, ncid), &
                  "sv_to_restart_file", "open "//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iyear'  , iyear), &
                  "sv_to_restart_file", "get_att iyear")
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'imonth' , imonth), &
                  "sv_to_restart_file", "get_att imonth")
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iday'   , iday), &
                  "sv_to_restart_file", "get_att iday")
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'ihour'  , ihour), &
                  "sv_to_restart_file", "get_att ihour")
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iminute', iminute), &
                  "sv_to_restart_file", "get_att iminute")
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'isecond', isecond), &
                  "sv_to_restart_file", "get_att isecond")

pop_time = set_date(iyear, imonth, iday, ihour, iminute, isecond)

if ( pop_time /= date1 ) then
   call print_time(   date1,'DART current time',logfileunit) 
   call print_time(pop_time,'POP  current time',logfileunit) 
   write(msgstring,*)trim(filename),' current time /= model time. FATAL error.'
   call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate) 
endif

if (do_output()) &
    call print_time(pop_time,'time of restart file '//trim(filename))
if (do_output()) &
    call print_date(pop_time,'date of restart file '//trim(filename))

! fill S, T, U, V in that order
do l=1, n3dfields

   varname = trim(progvarnames(l))//'_CUR'
   myerrorstring = trim(filename)//' '//trim(varname)

   ! Is the netCDF variable the right shape?

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            "sv_to_restart_file", "inq_varid "//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarId,dimids=dimIDs,ndims=numdims), &
            "sv_to_restart_file", "inquire "//trim(myerrorstring))

   if (numdims /= 3) then
      write(msgstring,*) trim(myerrorstring),' does not have exactly 3 dimensions'
      call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate)
   endif

   do i = 1,numdims
      write(msgstring,'(''inquire dimension'',i,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            "sv_to_restart_file", msgstring)

      if (dimlen /= size(data_3d_array,i)) then
         write(msgstring,*) trim(myerrorstring),'dim/dimlen',i,dimlen,'not',size(data_3d_array,i)
         call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate)
      endif
   enddo

   call vector_to_prog_var(state_vector,start_index(l),data_3d_array)

   ! Actually stuff it into the netcdf file
   call nc_check(nf90_put_var(ncid, VarID, data_3d_array), &
            'sv_to_restart_file', 'put_var '//trim(myerrorstring))

enddo

! and finally, SHGT (and any other 2d fields)
do l=(n3dfields+1), (n3dfields+n2dfields)

   select case ( trim(progvarnames(l)) )
   case ("SHGT")
      ! The restart files do not drag around SHGT, but they do drag
      ! around PSURF ... which can be used to calculate SHGT.
      varname = 'PSURF_CUR'
      convert_from_ssh = .true.
   case default
      varname = trim(progvarnames(l))//'_CUR'
      convert_from_ssh = .false.
   end select

   myerrorstring = trim(varname)//' '//trim(filename)

   ! Is the netCDF variable the right shape?

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            "sv_to_restart_file", "inq_varid "//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarId,dimids=dimIDs,ndims=numdims), &
            "sv_to_restart_file", "inquire "//trim(myerrorstring))

   if (numdims /= 2) then
      write(msgstring,*) trim(myerrorstring),' does not have exactly 2 dimensions'
      call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate)
   endif

   do i = 1,numdims
      write(msgstring,'(''inquire dimension'',i,A)')i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            "sv_to_restart_file", msgstring)

      if (dimlen /= size(data_2d_array,i)) then
         write(msgstring,*) trim(myerrorstring),'dim/dimlen',i,dimlen,'not',size(data_2d_array,i)
         call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate)
      endif
   enddo

   call vector_to_prog_var(state_vector,start_index(l),data_2d_array)

   if ( convert_from_ssh ) then  ! SSH = psurf*980.6    POP uses CGS 
      data_2d_array = data_2d_array * 980.6_r8
   endif

   call nc_check(nf90_put_var(ncid, VarID, data_2d_array), &
            'sv_to_restart_file', 'put_var '//trim(myerrorstring))

enddo

call nc_check(nf90_close(ncid), 'sv_to_restart_file', 'close '//trim(filename))

end subroutine sv_to_restart_file



subroutine vector_to_2d_prog_var(x, varindex, data_2d_array)
!------------------------------------------------------------------
!
real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: varindex
real(r4), dimension(:,:), intent(out) :: data_2d_array

integer :: i,j,ii
integer :: dim1,dim2
character(len=128) :: varname

if ( .not. module_initialized ) call static_init_model

dim1 = size(data_2d_array,1)
dim2 = size(data_2d_array,2)

varname = progvarnames(varindex)

if (dim1 /= Nx) then
   write(msgstring,*)trim(varname),' 2d array dim 1 ',dim1,' /= ',Nx
   call error_handler(E_ERR,'vector_to_2d_prog_var',msgstring,source,revision,revdate) 
endif
if (dim2 /= Ny) then
   write(msgstring,*)trim(varname),' 2d array dim 2 ',dim2,' /= ',Ny
   call error_handler(E_ERR,'vector_to_2d_prog_var',msgstring,source,revision,revdate) 
endif

ii = start_index(varindex)

do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   data_2d_array(i,j) = x(ii)
   ii = ii + 1
enddo
enddo

end subroutine vector_to_2d_prog_var



subroutine vector_to_3d_prog_var(x, varindex, data_3d_array)
!------------------------------------------------------------------
!
real(r8), dimension(:),     intent(in)  :: x
integer,                    intent(in)  :: varindex
real(r4), dimension(:,:,:), intent(out) :: data_3d_array

integer :: i,j,k,ii
integer :: dim1,dim2,dim3
character(len=128) :: varname

if ( .not. module_initialized ) call static_init_model

dim1 = size(data_3d_array,1)
dim2 = size(data_3d_array,2)
dim3 = size(data_3d_array,3)

varname = progvarnames(varindex)

if (dim1 /= Nx) then
   write(msgstring,*)trim(varname),' 3d array dim 1 ',dim1,' /= ',Nx
   call error_handler(E_ERR,'vector_to_3d_prog_var',msgstring,source,revision,revdate) 
endif
if (dim2 /= Ny) then
   write(msgstring,*)trim(varname),' 3d array dim 2 ',dim2,' /= ',Ny
   call error_handler(E_ERR,'vector_to_3d_prog_var',msgstring,source,revision,revdate) 
endif
if (dim3 /= Nz) then
   write(msgstring,*)trim(varname),' 3d array dim 3 ',dim3,' /= ',Nz
   call error_handler(E_ERR,'vector_to_3d_prog_var',msgstring,source,revision,revdate) 
endif

ii = start_index(varindex)

do k = 1,Nz   ! vertical
do j = 1,Ny   ! latitudes
do i = 1,Nx   ! longitudes
   data_3d_array(i,j,k) = x(ii)
   ii = ii + 1
enddo
enddo
enddo

end subroutine vector_to_3d_prog_var



  subroutine get_horiz_grid_dims(Nx, Ny)
!------------------------------------------------------------------
! subroutine get_horiz_grid_dims(Nx, Ny)
!
! Read in the lon, lat grid size from the netcdf file.
! Get the actual values later.
!
! The file name comes from module storage ... namelist.

 integer, intent(out) :: Nx   ! Number of Longitudes
 integer, intent(out) :: Ny   ! Number of Latitudes

 integer :: grid_id, dimid, nc_rc

 ! get the ball rolling ...

 call nc_check(nf90_open(trim(horiz_grid_input_file), nf90_nowrite, grid_id), &
         'get_horiz_grid_dims','open '//trim(horiz_grid_input_file))

 ! Longitudes : get dimid for 'i' or 'nlon', and then get value
 nc_rc = nf90_inq_dimid(grid_id, 'i', dimid)
 if (nc_rc /= nf90_noerr) then
   nc_rc = nf90_inq_dimid(grid_id, 'nlon', dimid)
   if (nc_rc /= nf90_noerr) then
      msgstring = 'unable to find either "i" or "nlon" in file'
      call error_handler(E_ERR, 'get_horiz_grid_dims', msgstring, &
                         source,revision,revdate) 
   endif
 endif 

 call nc_check(nf90_inquire_dimension(grid_id, dimid, len=Nx), &
         'get_horiz_grid_dims','inquire_dimension i '//trim(horiz_grid_input_file))

 ! Latitudes : get dimid for 'j' or 'nlat', and then get value
 nc_rc = nf90_inq_dimid(grid_id, 'j', dimid)
 if (nc_rc /= nf90_noerr) then
   nc_rc = nf90_inq_dimid(grid_id, 'nlat', dimid)
   if (nc_rc /= nf90_noerr) then
      msgstring = 'unable to find either "j" or "nlat" in file'
      call error_handler(E_ERR, 'get_horiz_grid_dims', msgstring, &
                         source,revision,revdate) 
   endif
 endif 

 call nc_check(nf90_inquire_dimension(grid_id, dimid, len=Ny), &
         'get_horiz_grid_dims','inquire_dimension i '//trim(horiz_grid_input_file))

 ! tidy up

 call nc_check(nf90_close(grid_id), &
         'get_horiz_grid_dims','close '//trim(horiz_grid_input_file) )

end subroutine get_horiz_grid_dims



  subroutine get_vert_grid_dims(Nz)
!------------------------------------------------------------------
! subroutine get_vert_grid_dims(Nz)
!
! count the number of lines in the ascii file to figure out max
! number of vert blocks.

 integer, intent(out) :: Nz

 integer :: linelen ! disposable

 call find_textfile_dims(vert_grid_input_file, Nz, linelen) 
 
end subroutine get_vert_grid_dims



subroutine get_gridsize(num_x, num_y, num_z)
 integer, intent(out) :: num_x, num_y, num_z
!------------------------------------------------------------------
! public utility routine.

 num_x = Nx
 num_y = Ny
 num_z = Nz

end subroutine get_gridsize



  subroutine read_horiz_grid(nx, ny, XC, XG, YC, YG)
!------------------------------------------------------------------
! subroutine read_horiz_grid(nx, ny, XC, XG, YC, YG)
!
! Open the netcdf file, read in the cell corners, and compute
! the cell centers.
!
! Here is the typically spartan dump of the grid netcdf file.
! netcdf horiz_grid.x1 {
! dimensions:
!         i = 320 ;
!         j = 384 ;
! variables:
!         double ULAT(j, i) ;
!         double ULON(j, i) ;
!         double HTN(j, i) ;
!         double HTE(j, i) ;
!         double HUS(j, i) ;
!         double HUW(j, i) ;
!         double ANGLE(j, i) ;
! }
! FIXME:
! For the LANL global case, the ULAT,ULON were in radians.
! For the CCSM global case, they were in degrees.  The CCSM files do
! have units that include the word 'degrees' (degrees_north, degrees_east)
! so perhaps we can look for them.  The LANL file had no attributes
! to query.

 integer,  intent(in)  :: nx, ny
 real(r8), intent(out) :: XC(:), XG(:)
 real(r8), intent(out) :: YC(:), YG(:)

 integer  :: grid_id, i, j, nc_rc
 integer  :: ulat_id, ulon_id, tlat_id, tlon_id
 logical  :: read_tgrid, in_radians
 character(len=128) :: grid_units
 real(r8) :: wrapval

 call nc_check(nf90_open(trim(horiz_grid_input_file), nf90_nowrite, grid_id), &
         'read_horiz_grid','open '//trim(horiz_grid_input_file) )

 call nc_check(nf90_inq_varid(grid_id, 'ULAT', ulat_id), &
         'read_horiz_grid','varid ULAT '//trim(horiz_grid_input_file) )

 ! FIXME:
 ! some files have ULON, some ULONG.  try to read either.  same with T
 nc_rc = nf90_inq_varid(grid_id, 'ULON', ulon_id)
 if (nc_rc /= nf90_noerr) then
   nc_rc = nf90_inq_varid(grid_id, 'ULONG', ulon_id)
   if (nc_rc /= nf90_noerr) then
      msgstring = 'unable to find either ULON or ULONG in file'
      call error_handler(E_ERR, 'read_horiz_grid', msgstring, &
                         source,revision,revdate) 
   endif
 endif 

 call nc_check(nf90_get_var(grid_id, ulat_id, ULAT), &
         'read_horiz_grid','get ULAT '//trim(horiz_grid_input_file) )

 call nc_check(nf90_get_var(grid_id, ulon_id, ULON), &
         'read_horiz_grid','get ULON '//trim(horiz_grid_input_file) )

 ! FIXME:
 ! some files have the grid info in radians, some in degrees.
 ! see if there is a units attribute, and if it contains the string 'degrees'
 ! (e.g. 'degrees_north').  if not found, try to be tricky and look
 ! at the min/max -- THIS ASSUMES A GLOBAL RUN.
 
 nc_rc = nf90_get_att(grid_id, ulon_id, 'units', grid_units)
 if (nc_rc == nf90_noerr) then
   ! found attribute; see if it contains the substring degrees
   i = index(grid_units, 'degrees')
   if (i /= 0) then
      in_radians = .FALSE.
   else
      in_radians = .TRUE.
   endif
   if (debug > 1) print *, 'found attribute, i, in_radians = ', i, in_radians
 else
   ! FIXME: assumes global grid
   ! try the max/min and for now assume it is a global grid
   if (maxval(ULAT) < 7.0) then   ! 2 * PI plus slop.
      in_radians = .TRUE.
   else
      in_radians = .FALSE.
   endif
   if (debug > 1) print *, 'no attribute, max, in_radians = ', maxval(ULAT), in_radians
 endif

 ! check to see if grid is wrapped at 0 and make all positive.
 if (minval(ULON) < 0.0) then
    if (in_radians) then
       where (ULON < 0.0) ULON = ULON + 2.0 * PI
    else
       where (ULON < 0.0) ULON = ULON + 360.0_r8
    endif
 endif
 

 ! if the file contains the TLON, TLAT grids directly, read them in.
 ! otherwise, compute them after you have the U grids.
 nc_rc = nf90_inq_varid(grid_id, 'TLAT', tlat_id)
 if (nc_rc == nf90_noerr) then
   read_tgrid = .true.
   call nc_check(nf90_inq_varid(grid_id, 'TLAT', tlat_id), &
           'read_horiz_grid','varid TLAT '//trim(horiz_grid_input_file) )
  
   ! look for TLON or TLONG
   nc_rc = nf90_inq_varid(grid_id, 'TLON', tlon_id)
   if (nc_rc /= nf90_noerr) then
     nc_rc = nf90_inq_varid(grid_id, 'TLONG', tlon_id)
     if (nc_rc /= nf90_noerr) then
        msgstring = 'unable to find either TLON or TLONG in file'
        call error_handler(E_ERR, 'read_horiz_grid', msgstring, &
                           source,revision,revdate) 
     endif
   endif 

   call nc_check(nf90_get_var(grid_id, tlat_id, TLAT), &
           'read_horiz_grid','get TLAT '//trim(horiz_grid_input_file) )
  
   call nc_check(nf90_get_var(grid_id, tlon_id, TLON), &
           'read_horiz_grid','get TLON '//trim(horiz_grid_input_file) )

   ! check to see if grid is wrapped at 0 and make all positive.
   if (minval(TLON) < 0.0) then
      if (in_radians) then
         where (TLON < 0.0) TLON = TLON + 2.0 * PI
      else
         where (TLON < 0.0) TLON = TLON + 360.0_r8
      endif
   endif
 
 else
   read_tgrid = .false.
 endif

 ! TJH check sizes ...

 call nc_check(nf90_close(grid_id), &
         'read_horiz_grid','close '//trim(horiz_grid_input_file) )

 ! FIXME:
 !  even though the real grid is a dipole grid, we are going to treat it
 !  as reg lat/lon for initial testing.   use the first row of longs as
 !  the universal lons, and for now, take the first row of lats, because
 !  they are monotonic (i checked).

 ! lons: first row
 if (in_radians) ULON = ULON * rad2deg
 do i=1, nx
   XG(i) = ULON(i, 1)
   if ( i < 5) print *, ULON(i, 1), XG(i)
 enddo

 ! lats: first row is monotonic in dipole grid
 if (in_radians) ULAT = ULAT * rad2deg
 do j=1, ny
   YG(j) = ULAT(1, j) 
   if ( j < 5) print *, ULAT(1, j), YG(j)
 enddo

 ! FIXME: i am assuming that if the U grid was in radians or degrees, and 
 !   if the T grid is also in the file, that it is consistent units.

 if (read_tgrid) then
   ! we read in the centers directly from netcdf file, just fill in arrays

   ! lons: first row
   if (in_radians) TLON = TLON * rad2deg
   do i=1, nx
     XC(i) = TLON(i, 1)
     if ( i < 5) print *, TLON(i, 1), XC(i)
   enddo
  
   ! lats: first row is monotonic in dipole grid
   if (in_radians) TLAT = TLAT * rad2deg
   do j=1, ny
     YC(j) = TLAT(1, j) 
     if ( j < 5) print *, TLAT(1, j), YC(j)
   enddo

 else
   ! T grid not in input file, so compute centers for T grid

   ! longitudes
   do i=2, nx
     XC(i) = XG(i-1) + (XG(i) - XG(i-1)) * 0.5_r8
   enddo
  
   if (longitude_wrap) then
     if (in_radians) then
       wrapval = 2 * PI
     else
       wrapval = 360.0_r8
     endif
     XC(1) = XG(nx) + (XG(1)+wrapval - XG(nx)) * 0.5_r8 
   else
     ! FIXME: what to do here? extrapolate back .5 of box 1 for now.
     XC(1) = XG(1) - ((XG(2) - XG(1)) * 0.5_r8)
   endif
  
   ! latitudes
   do j=2, ny
     YC(j) = YG(j-1) + (YG(j) - YG(j-1)) * 0.5_r8
   enddo
  
   ! FIXME: what to do here? extrapolate back .5 of box 1 for now.
   YC(1) = YG(1) - ((YG(2) - YG(1)) * 0.5_r8)
 endif

if (debug > 5) then
  print *, 'XG, XC, YG, YC:'
  print *, '1: ' ,XG(1), XC(1), YG(1), YC(1)
  print *, '2: ' ,XG(2), XC(2), YG(2), YC(2)
  print *, '3: ' ,XG(3), XC(3), YG(3), YC(3)
  print *, '4: ' ,XG(4), XC(4), YG(4), YC(4)
endif

end subroutine read_horiz_grid



subroutine read_vert_grid(nz, ZC, ZG)
!------------------------------------------------------------------
 integer,  intent(in)  :: nz
 real(r8), intent(out) :: ZC(:), ZG(:)

 integer  :: iunit, i, ios
 real(r8) :: depth(nz), zcent(nz), zbot(nz)

 iunit = open_file(trim(vert_grid_input_file), action = 'read')

 do i=1, nz

   read(iunit,*,iostat=ios) depth(i), zcent(i), zbot(i)

   ! error
   if ( ios /= 0 ) then
      write(*,msgstring)'error reading depths, line ',i
      call error_handler(E_ERR,'read_vert_grid',msgstring,source,revision,revdate) 
   endif

 enddo 

 do i=1, nz
   ZC(i) = zcent(i)
   ZG(i) = zbot(i)   ! FIXME: this might actually be off by one, but i cannot
                     !  see that we ever use this in any interpolation.
 enddo
end subroutine read_vert_grid



  subroutine read_kmt(Nx, Ny, KMT, KMU)
!------------------------------------------------------------------
! subroutine read_kmt(Nx, Ny, KMT, KMU)
! 
! open topo file
! make sure i, j match Nx, Ny
! KMT array already allocated - fill from file
! close topo file

 integer, intent(in)  :: Nx, Ny
 integer, intent(out) :: KMT(:, :), KMU(:, :)

 integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
 integer :: ncid, varid, numdims
 integer :: dim1len, dim2len, nc_rc, i, j, k
 character(len=4) :: depthvar
  

 call nc_check(nf90_open(trim(topography_input_file), nf90_nowrite, ncid), &
         'read_kmt','open '//trim(topography_input_file) )

 ! Try to read in KMT.  If there, great.  If not, try KT.  If still not, fail.
 nc_rc = nf90_inq_varid(ncid, 'KMT', varid)
 if (nc_rc == nf90_noerr) then
    depthvar = 'KMT'
 else
    nc_rc = nf90_inq_varid(ncid, 'HT', varid)
    if (nc_rc /= nf90_noerr) then
       write(msgstring,*)trim(topography_input_file)//' cannot find KMT or HT'
       call error_handler(E_ERR,'read_kmt',msgstring,source,revision,revdate) 
    endif
    depthvar = 'HT'
 endif

 call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimIDs,ndims=numdims), &
         'read_kmt', 'inquire '//depthvar//' '//trim(topography_input_file))

 if ( numdims /= 2 ) then
    write(msgstring,*)trim(topography_input_file)//' '//depthvar//'should have 2 dims - has ',numdims
    call error_handler(E_ERR,'read_kmt',msgstring,source,revision,revdate) 
 endif

 call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dim1len), &
         'read_kmt', 'inquire dim1 of '//depthvar//' '//trim(topography_input_file))

 call nc_check(nf90_inquire_dimension(ncid, dimIDs(2), len=dim2len), &
         'read_kmt', 'inquire dim2 of '//depthvar//' '//trim(topography_input_file))

 if ( (dim1len /= Nx) .or. (dim2len /= Ny) ) then
    write(msgstring,*)trim(topography_input_file)//' '//depthvar//'shape wrong ',dim1len,dim2len
    call error_handler(E_ERR,'read_kmt',msgstring,source,revision,revdate) 
 endif

 ! If we made it this far, we can read the variable
 if (depthvar == 'KMT') then
    call nc_check(nf90_get_var(ncid, varid, KMT), &
            'read_kmt','get_var KMT '//trim(topography_input_file) )
 else
    call nc_check(nf90_get_var(ncid, varid, HT), &
            'read_kmt','get_var HT '//trim(topography_input_file) )
    
    TLATLOOP: do j=1, dim2len
       TLONLOOP: do i=1, dim1len
          KMT(i, j) = 0
          if (HT(i, j) /= 0.0_r8) then
             TDEPTHLOOP: do k=Nz, 1, -1
                if (HT(i, j) >= ZC(k)*1000.0_r8) then   ! HT in cm, ZC in m
                   KMT(i,j) = k
                   exit TDEPTHLOOP
                endif
             enddo TDEPTHLOOP
          endif
       enddo TLONLOOP
    enddo TLATLOOP

 endif

 ! Try to read in KMU.  If there, great.  If not, try HU.  If not, compute it.
 nc_rc = nf90_inq_varid(ncid, 'KMU', varid)
 if (nc_rc == nf90_noerr) then
    call nc_check(nf90_get_var(ncid, varid, KMU), &
            'read_kmt','get_var KMU '//trim(topography_input_file) )
 else 
    nc_rc = nf90_inq_varid(ncid, 'HU', varid)
    if (nc_rc == nf90_noerr) then
       call nc_check(nf90_get_var(ncid, varid, HU), &
            'read_kmt','get_var HU '//trim(topography_input_file) )

       ULATLOOP: do j=1, dim2len
          ULONLOOP: do i=1, dim1len
             KMU(i, j) = 0
             if (HU(i, j) /= 0.0_r8) then
                UDEPTHLOOP: do k=Nz, 1, -1
                   if (HU(i, j) >= ZC(k)*1000.0_r8) then   ! HU in cm, ZC in m
                      KMU(i,j) = k
                      exit UDEPTHLOOP
                   endif
                enddo UDEPTHLOOP
             endif
          enddo ULONLOOP
       enddo ULATLOOP
    else  ! neither KMU or HU found; compute from KMT
       KMU(1, 1) = 0
       do j=2, dim2len
          do i=2, dim1len
             KMU(i,j) = min(KMT(i, j), KMT(i-1, j), KMT(i, j-1), KMT(i-1, j-1)) 
          enddo
       enddo
    endif
 endif

 call nc_check(nf90_close(ncid), &
         'read_kmt','close '//trim(topography_input_file) )

end subroutine read_kmt


  function is_not_ocean(obs_type, lon_index, lat_index, hgt_index)
!------------------------------------------------------------------
! returns true if this point is below the ocean floor or if it is
! on land.
integer, intent(in)  :: obs_type
integer, intent(in)  :: lon_index, lat_index, hgt_index
logical              :: is_not_ocean

logical :: is_ugrid

is_ugrid = is_on_ugrid(obs_type)
if ((      is_ugrid .and. hgt_index > KMU(lon_index, lat_index)) .or. &
    (.not. is_ugrid .and. hgt_index > KMT(lon_index, lat_index))) then
   is_not_ocean = .TRUE.
   return
endif

end function


  function is_on_ugrid(obs_type)
!------------------------------------------------------------------
!  returns true if U, V -- everything else is on T grid
integer, intent(in) :: obs_type
logical             :: is_on_ugrid

is_on_ugrid = .FALSE.

if ((obs_type == KIND_U_CURRENT_COMPONENT)  .or.  &
    (obs_type == KIND_V_CURRENT_COMPONENT)) is_on_ugrid = .TRUE.

end function


  subroutine write_grid_netcdf()
!------------------------------------------------------------------
!
! Write the grid to a netcdf file for checking.

 integer :: ncid, NlonDimID, NlatDimID, NzDimID
 integer :: nlon, nlat, nz
 integer :: ulatVarID, ulonVarID, TLATvarid, TLONvarid
 integer :: ZGvarid, ZCvarid, KMTvarid, KMUvarid

 integer :: dimids(2);

 nlon = size(ULAT,1)
 nlat = size(ULAT,2)
 nz   = size(ZG)
 
 call nc_check(nf90_create('dart_grid.nc', NF90_CLOBBER, ncid),'write_grid_netcdf')

 ! define dimensions

 call nc_check(nf90_def_dim(ncid, 'i', nlon, NlonDimID),'write_grid_netcdf')
 call nc_check(nf90_def_dim(ncid, 'j', nlat, NlatDimID),'write_grid_netcdf')
 call nc_check(nf90_def_dim(ncid, 'k',   nz,   NzDimID),'write_grid_netcdf')

 dimids(1) = NlonDimID 
 dimids(2) = NlatDimID 

 ! define variables

 ! FIXME: we should add attributes to say what units the grids are in (degrees).
 call nc_check(nf90_def_var(ncid,  'KMT', nf90_int,     dimids,  KMTvarid),'write_grid_netcdf')
 call nc_check(nf90_def_var(ncid,  'KMU', nf90_int,     dimids,  KMUvarid),'write_grid_netcdf')
 call nc_check(nf90_def_var(ncid, 'ULON', nf90_double,  dimids, ulonVarID),'write_grid_netcdf')
 call nc_check(nf90_def_var(ncid, 'ULAT', nf90_double,  dimids, ulatVarID),'write_grid_netcdf')
 call nc_check(nf90_def_var(ncid, 'TLON', nf90_double,  dimids, TLONvarid),'write_grid_netcdf')
 call nc_check(nf90_def_var(ncid, 'TLAT', nf90_double,  dimids, TLATvarid),'write_grid_netcdf')
 call nc_check(nf90_def_var(ncid,   'ZG', nf90_double, NzDimID,   ZGvarid),'write_grid_netcdf')
 call nc_check(nf90_def_var(ncid,   'ZC', nf90_double, NzDimID,   ZCvarid),'write_grid_netcdf')

 call nc_check(nf90_put_att(ncid,ulonVarID,"long_name","U,V grid lons"), &
                                                       'write_grid_netcdf')
 call nc_check(nf90_put_att(ncid,ulatVarID,"long_name","U,V grid lats"), &
                                                       'write_grid_netcdf')
 call nc_check(nf90_put_att(ncid,tlonVarID,"long_name","S,T grid lons"), &
                                                       'write_grid_netcdf')
 call nc_check(nf90_put_att(ncid,tlatVarID,"long_name","S,T grid lats"), &
                                                      'write_grid_netcdf')

 call nc_check(nf90_enddef(ncid),'write_grid_netcdf')

 ! fill variables

 call nc_check(nf90_put_var(ncid,  KMTvarid,  KMT),'write_grid_netcdf')
 call nc_check(nf90_put_var(ncid,  KMUvarid,  KMU),'write_grid_netcdf')
 call nc_check(nf90_put_var(ncid, ulatVarID, ULAT),'write_grid_netcdf')
 call nc_check(nf90_put_var(ncid, ulonVarID, ULON),'write_grid_netcdf')
 call nc_check(nf90_put_var(ncid, TLATvarid, TLAT),'write_grid_netcdf')
 call nc_check(nf90_put_var(ncid, TLONvarid, TLON),'write_grid_netcdf')
 call nc_check(nf90_put_var(ncid,   ZGvarid,   ZG),'write_grid_netcdf')
 call nc_check(nf90_put_var(ncid,   ZCvarid,   ZC),'write_grid_netcdf')

 call nc_check(nf90_close(ncid),'write_grid_netcdf')

end subroutine write_grid_netcdf

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

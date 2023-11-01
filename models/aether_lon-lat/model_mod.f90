! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! This module was copied from models/tiegcm 
! but has restart reading and writing routines from ../gitm,
! because the lon-lat grid layout, with halos, and subdomain ("block") file structure
! seems to be the same in GITM and Aether,
! Those subroutines need to be adapted to the infrastructure in this model_mod
! and to the Aether restart files' format and contents.
! Later they will be exported to a model_mod Ben is developing from scratch.

! The model_mod.nml initially has the namelists from both tiegcm and gitm.
! Parts of both may be useful and will be merged into a new aether_lon-lat nml.


module model_mod

!-------------------------------------------------------------------------------
!
! Interface for HAO-TIEGCM 2.0
!
!-------------------------------------------------------------------------------

use        types_mod, only : r4, r8, i8, MISSING_R8, MISSING_R4, PI, RAD2DEG, &
                             earth_radius, gravity, obstypelength, MISSING_I

use time_manager_mod, only : time_type, set_calendar_type, set_time_missing,        &
                             set_time, get_time, print_time,                        &
                             set_date, get_date, print_date,                        &
                             operator(*),  operator(+), operator(-),                &
                             operator(>),  operator(<), operator(/),                &
                             operator(/=), operator(<=)

use     location_mod, only : location_type,                                         &
                             get_close_obs,                                         &
! TODO: need this from Ben's model_mod
                             loc_get_close_state => get_close_state,                &
                             set_location, get_location,                            &
                             get_dist, query_location,                              &
                             get_close_type, VERTISUNDEF,                           &
                             VERTISPRESSURE, VERTISHEIGHT, VERTISLEVEL,             &
                             vertical_localization_on, set_vertical

use    utilities_mod, only : open_file, file_exist, close_file, logfileunit,        &
                             error_handler, E_ERR, E_MSG, E_WARN, nmlfileunit,      &
                             do_output, find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, register_module,             &
                             file_to_text, find_textfile_dims, to_upper

! TODO: will need many more kinds, and maybe new kinds (6 velocity components, ...?)
use     obs_kind_mod, only : QTY_U_WIND_COMPONENT,           &
                             QTY_V_WIND_COMPONENT,           &
                             QTY_TEMPERATURE,                &! neutral temperature obs
                             QTY_PRESSURE,                   &! neutral pressure obs
                             QTY_MOLEC_OXYGEN_MIXING_RATIO,  &! neutral composition obs
                             QTY_1D_PARAMETER,               &
                             QTY_GEOPOTENTIAL_HEIGHT,        &
                             QTY_GEOMETRIC_HEIGHT,           &
                             QTY_VERTICAL_TEC,               &! total electron content
                             QTY_DENSITY_ION_OP,             &! O+
                             get_index_for_quantity

use     quad_utils_mod,  only : quad_interp_handle, init_quad_interp, &
                                set_quad_coords, finalize_quad_interp, &
                                quad_lon_lat_locate, quad_lon_lat_evaluate, &
                                GRID_QUAD_IRREG_SPACED_REGULAR,  &
                                QUAD_LOCATED_CELL_CENTERS

use mpi_utilities_mod,only : my_task_id

use default_model_mod, only : adv_1step,                                &
                              init_conditions => fail_init_conditions,  &
                              init_time => fail_init_time,              &
                              nc_write_model_vars,                      &
                              pert_model_copies

use state_structure_mod, only : add_domain, get_dart_vector_index, add_dimension_to_variable, &
                                finished_adding_domain, state_structure_info, &
                                get_domain_size, get_model_variable_indices, &
                                get_num_dims, get_dim_name, get_variable_name, &
                                get_varid_from_varname, get_num_varids_from_kind, &
                                get_varid_from_kind, get_varids_from_kind, &
                                hyperslice_domain, get_num_domains

use distributed_state_mod, only : get_state, get_state_array

use ensemble_manager_mod, only : ensemble_type

use netcdf_utilities_mod, only : nc_synchronize_file, nc_add_global_attribute,           &
                                 nc_add_global_creation_time, nc_begin_define_mode,      &
                                 nc_define_dimension, nc_end_define_mode,                &
                                 nc_put_variable,nc_add_attribute_to_variable,           &
                                 nc_define_real_variable,                                &
                                 nc_check, nc_open_file_readonly, nc_get_dimension_size, &
                                 nc_close_file, nc_get_variable,                         &
                                 nc_get_dimension_size, nc_create_file,                  &
                                 nc_define_double_variable, nc_define_double_scalar
                                 

use dart_time_io_mod,     only : write_model_time

use netcdf

implicit none
private

!DART mandatory public interfaces
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          end_model,              &
          static_init_model,      &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          get_close_obs,          &
          get_close_state,        &
          shortest_time_between_assimilations, &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          write_model_time

!DART pass through interfaces
public :: adv_1step,              &
          init_conditions,        &
          init_time,              &
          pert_model_copies

! Interfaces needed by other programs, e.g.  aether_to_dart and dart_to_aether
! block_file_name creates an Aether restart file name,
! which is useful for read_model_time calls, and others.
public :: restart_files_to_netcdf, &
          block_file_name

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = 'aether_lon-lat/model_mod.f90'
character(len=32 ), parameter :: revision = ''
character(len=128), parameter :: revdate  = ''

!-------------------------------------------------------------------------------
! namelist with default values

! TODO: Define a derived type to handle the file types which need to be read andor written?
! PRobably not; variable_table probably handles it all.
! file_root {'neutrals','ions', 'time', f10_7? ...?)
! file_ext  {'nc',  'nc',       'json', 'nc',    ...)
! num_fields {nfields_neutral, nfields_ion, 2, ?, ...)
! character(len=8), dimension(2) :: file_root = /('neutrals','ions'/)

! TODO; does it actually need filter_io_dir, or will the scripts make these programs
!       run where the restart files are?
character(len=256) :: filter_io_dir = '.'
! TODO; remove GITM namelist vars
! TODO: if filter_io_filename is in global storage it doesn't need to be in (some?) arg lists.
character(len=256) :: filter_io_filename = 'no_file_specified.nc'
integer            :: debug = 0
logical            :: estimate_f10_7 = .false.
character(len=256) :: f10_7_file_name = 'f10_7.nc'
real(r8)           :: model_res = 5.0_r8

! TODO: confirm that the units are days.
!       Better to get the actual start day of Aether's calender.
integer            :: aeth_ref_day = 2451545.0  ! cJULIAN2000 in Aether = day of date 2000/01/01.
character(len=32)  :: calendar = 'Gregorian'
! Day 0 in this calendar is (+/1 a day) -4710/11/24 0 UTC
! But what we care about is the ref time for the times in the files, which is 1964-12-31 23:30
! (from echo 2011032000 -1458345600s | ./advance_time).

integer, dimension(:) :: aeth_ref_date(5) = (/1965,1,1,0,0/)  ! y,mo,d,h,m (secs assumed 0)
type(time_type)       :: aeth_ref_time
integer               :: aeth_ref_ndays, aeth_ref_nsecs

integer               :: assimilation_period_seconds = 3600

! TODO: Aether restart files have 81 fields in them, 
!       mostly the 6 components of velocities for each ion.
!       Aaron will provide files with a few more fields; e-, f10_7, ...?
integer, parameter :: MAX_NUM_VARIABLES = 100
integer, parameter :: MAX_NUM_COLUMNS = 6
character(len=NF90_MAX_NAME) :: variables(MAX_NUM_VARIABLES * MAX_NUM_COLUMNS) = ' '

namelist /model_nml/ filter_io_dir, &
                     variables, debug, estimate_f10_7, &
                     f10_7_file_name, calendar, assimilation_period_seconds, &
                     model_res
                     
!-------------------------------------------------------------------------------
! define model parameters for creating the state NetCDF file 
! and handling interpolation, get_close, ...

! nalt is number of midpoint levels
! TODO: Replace plevs with hlevs?  Maybe not; pressure levels may be needed for interp.
!       nilev -> an aether dimension size, that's not interface levels.
!         Is not used by aether_to_dart or dart_to_aether (?).
integer                               :: nalt, nlon, nlat, nilev
real(r8),dimension(:),    allocatable :: lons, lats, alts, plevs, ilevs
! HK are plevs, pilevs per ensemble member?
real(r8)                              :: TIEGCM_reference_pressure
integer                               :: time_step_seconds
integer                               :: time_step_days
type(time_type)                       :: time_step

type(quad_interp_handle) :: quad_interp

! Codes for interpreting the columns of the variable_table
integer, parameter :: VT_VARNAMEINDX  = 1 ! variable name
integer, parameter :: VT_KINDINDX     = 2 ! DART quantity
integer, parameter :: VT_MINVALINDX   = 3 ! minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! maximum value if any
integer, parameter :: VT_ORIGININDX   = 5 ! file of origin
integer, parameter :: VT_STATEINDX    = 6 ! update (state) or not

character(len=obstypelength) :: variable_table(MAX_NUM_VARIABLES, MAX_NUM_COLUMNS)

type(time_type) :: state_time ! module-storage declaration of current model time

integer(i8)           :: model_size ! the state vector length
integer :: nfields, nfields_neutral, nfields_ion  ! numbers of aether variables in DART state

! lon and lat grid specs. 
real(r8)  :: bot_lon        = MISSING_R8
real(r8)  :: top_lon        = MISSING_R8
real(r8)  :: delta_lon      = MISSING_R8
real(r8)  :: bot_lat        = MISSING_R8
real(r8)  :: top_lat        = MISSING_R8
real(r8)  :: delta_lat      = MISSING_R8
integer   :: zero_lon_index = MISSING_I


! Obs locations are expected to be given in height [m] or level,
! and so vertical localization coordinate is *always* height.

character(len=512)    :: string1, string2, string3
logical, save         :: module_initialized = .false.

!===============================================================================
! Define Aether whole-grid and block grid dimension variables.

character(len=*), parameter :: LON_DIM_NAME = 'lon'
character(len=*), parameter :: LAT_DIM_NAME = 'lat'
character(len=*), parameter :: ALT_DIM_NAME = 'alt'

character(len=*), parameter :: LON_VAR_NAME = 'lon'
character(len=*), parameter :: LAT_VAR_NAME = 'lat'
character(len=*), parameter :: ALT_VAR_NAME = 'alt'

! {nxPerBlock,nyPerBlock} are the number of non-halo {lons,lats} PER block
!                  the number of blocks comes from UAM.in
! nzPerBlock  is the number of altitudes, which does not depend on block
! nGhost is the halo region width in the block(subdomain) files.
! TODO: n[xyz]PerBlock should probably come from a namelist (aether_to_dart.nml;
!       can that be used for dart_to_aether?)

integer :: nxPerBlock, nyPerBlock, nzPerBlock
integer, parameter :: nGhost = 2   ! number of ghost cells on all edges

! "... keep in mind that if the model resolution is 5 deg latitude,
!  the model will actually go from -87.5 to 87.5 latitude
! (even though you specify -90 to 90 in the UAM.in file),
! since the latitudes/longitudes are at cell centers,
! while the edges are at the boundaries." -- Aaron Ridley

! number of blocks along each dim
integer  :: nBlocksLon=MISSING_I, nBlocksLat=MISSING_I, nBlocksAlt=MISSING_I 
real(r8) :: LatStart=MISSING_R8, LatEnd=MISSING_R8, LonStart=MISSING_R8

contains
!===============================================================================

subroutine static_init_model()

integer :: iunit, io
character(len=*), parameter :: routine = 'static_init_model'

character(len=128) :: aether_filename

if (module_initialized) return ! only need to do this once
! Print module information to log file and stdout.
call register_module(source, revision, revdate)

module_initialized = .true.

! Read the namelist entry for model_mod from input.nml
call read_model_namelist()

if (do_output()) then
   write(     *     ,*)'static_init_model: debug level is ',debug
   write(logfileunit,*)'static_init_model: debug level is ',debug
endif

!---------------------------------------------------------------
! get whole grid dimensions and values

write(string1,'(3A)') "Now reading filter_io file ",trim(filter_io_filename),&
   " for grid information"
call error_handler(E_MSG,routine,string1,source,revision,revdate)

allocate(lons(nlon))
allocate(lats(nlat))
allocate(alts(nalt))

!---------------------------------------------------------------
! get grid dimensions and values
call get_grid_from_netcdf(filter_io_filename, lons, lats, alts)

!---------------------------------------------------------------

! mass points at cell centers
call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, nlon, nlat, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.false., spans_lon_zero=.false., pole_wrap=.false., &
                      interp_handle=quad_interp)

call set_quad_coords(quad_interp, lons, lats)

if ( debug > 0 ) then
   write(string1,'("grid: nlon, nlat, nalt =",3(1x,i5))') nlon, nlat, nalt
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
endif

if ( estimate_f10_7 ) then
   call error_handler(E_MSG, 'f10_7 part of DART state', source)
endif

! error-check, convert namelist input to variable_table, and build the
! state structure
call make_variable_table()

call set_calendar_type(calendar)

! Read and convert the time (seconds from the aether_ref_date) to a dart time.
aether_filename = block_file_name(variable_table(1,VT_ORIGININDX), 0, 0)
state_time = read_model_time(trim(aether_filename))

! Initialized in namelist.
time_step = set_time(assimilation_period_seconds, 0)

end subroutine static_init_model

!==================================================================

! Read the lon, lat, and alt arrays from the ncid

subroutine get_grid_from_netcdf(filter_io_filename, lons, lats, alts )

character(len=*), intent(in)    :: filter_io_filename
real(r8),         intent(inout) :: lons(:)
real(r8),         intent(inout) :: lats(:)
real(r8),         intent(inout) :: alts(:)

character(len=*), parameter :: routine = 'get_grid_from_netcdf'

integer :: ncid

ncid = nc_open_file_readonly(filter_io_filename, routine)

call nc_get_variable(ncid, LON_VAR_NAME, lons, routine)
call nc_get_variable(ncid, LAT_VAR_NAME, lats, routine)
call nc_get_variable(ncid, ALT_VAR_NAME, alts, routine)

call nc_close_file(ncid)

end subroutine get_grid_from_netcdf

!=================================================================

subroutine static_init_blocks(restart_dirname)

character(len=*), intent(in)  :: restart_dirname
character(len=128) :: aether_filename

character(len=*), parameter :: routine = 'static_init_blocks'

character(len=NF90_MAX_NAME)    :: varname
integer :: iunit, io, ivar
!logical :: has_gitm_namelist

if (module_initialized) return ! only need to do this once

! This prevents subroutines called from here from calling static_init_mod.
module_initialized = .true.

! Read the namelist entry for model_mod from input.nml
call read_model_namelist()

! error-check, convert namelist input to variable_table, and build the state structure
call make_variable_table()

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! TODO: Reading aether_to_dart_nml is done only in aether_to_dart?
!       filter_io_dir from here instead of redundant entry in model_mod_nml?
! ! Read the DART namelist for this model
! call find_namelist_in_file('input.nml', 'aether_to_dart_nml', iunit)
! read(iunit, nml = aether_to_dart_nml, iostat = io)
! call check_namelist_read(iunit, io, 'aether_to_dart_nml')
! 
! ! Record the namelist values used for the run
! if (do_nml_file()) write(nmlfileunit, nml=aether_to_dart_nml)
! if (do_nml_term()) write(     *     , nml=aether_to_dart_nml)

!---------------------------------------------------------------
! Set the time step ... causes gitm namelists to be read.
! Ensures model_advance_time is multiple of 'dynamics_timestep'

!TODO: Aether uses Julian time internally
!      andor a Julian calendar (days from the start of the calendar), depending on the context)
call set_calendar_type( calendar )   ! comes from model_mod_nml

!---------------------------------------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids
! 3) read them from the block restart files, could be stretched ...

call get_grid_info_from_blocks(restart_dirname, nlon, nlat, nalt, nBlocksLon, &
               nBlocksLat, nBlocksAlt, LatStart, LatEnd, LonStart)
print*,'static_init_blocks: post-get_grid_info_from_blocks; nfields_neutral = ', nfields_neutral

if( debug  > 0 ) then
    write(string1,*) 'grid dims are ',nlon,nlat,nalt
    call error_handler(E_MSG,routine,string1,source,revision,revdate)
endif

! Opens and closes the grid block file, but not the filter netcdf file.
call get_grid_from_blocks(restart_dirname, nBlocksLon, nBlocksLat, nBlocksAlt, &
   nxPerBlock, nyPerBlock, nzPerBlock, lons, lats, alts )

! Convert the Aether reference date (not calendar day = 0 date)
! to the days and seconds of the calendar set in model_mod_nml.
aeth_ref_time = set_date(aeth_ref_date(1), aeth_ref_date(2), aeth_ref_date(3), &
                     aeth_ref_date(4), aeth_ref_date(5))
call get_time(aeth_ref_time,aeth_ref_nsecs,aeth_ref_ndays)

! Get the model time from a restart file.
aether_filename = block_file_name(variable_table(1,VT_ORIGININDX), 0, 0)
state_time = read_model_time(trim(restart_dirname)//'/'//trim(aether_filename))

! TODO: Replace with aether variables check? (OR is that done when trying to read them?) 
! call verify_block_variables( gitm_block_variables, nfields )
! 
! do ivar = 1, nfields
! 
!    varname                   = trim(gitm_block_variables(ivar))
!    gitmvar(ivar)%varname     = varname
! 
!    ! This routine also checks to make sure user specified accurate GITM variables
!    call decode_gitm_indices( varname,                    &
!                              gitmvar(ivar)%gitm_varname, &
!                              gitmvar(ivar)%gitm_dim,     &
!                              gitmvar(ivar)%gitm_index,   &
!                              gitmvar(ivar)%long_name,    &
!                              gitmvar(ivar)%units)
!    if ( debug > 0 ) then
!       call print_gitmvar_info(ivar,routine)
!    endif
! enddo

if ( debug > 0 ) then
  write(string1,'("grid: nlon, nlat, nalt =",3(1x,i5))') nlon, nlat, nalt
  call error_handler(E_MSG,routine,string1,source,revision,revdate)
endif

end subroutine static_init_blocks

!==================================================================

subroutine read_model_namelist()

integer :: iunit, io

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

end subroutine read_model_namelist

!==================================================================

!> Read the grid dimensions from a restart netcdf file.
!>
!> The file name comes from module storage ... namelist.

subroutine get_grid_info_from_blocks(restart_dirname, nlon, nlat, &
                nalt, nBlocksLon, nBlocksLat, nBlocksAlt, LatStart, LatEnd, LonStart)

character(len=*), intent(in) :: restart_dirname
integer,  intent(out) :: nlon   ! Number of Longitude centers
integer,  intent(out) :: nlat   ! Number of Latitude  centers
integer,  intent(out) :: nalt   ! Number of Vertical grid centers
integer,  intent(out) :: nBlocksLon, nBlocksLat, nBlocksAlt
real(r8), intent(out) :: LatStart, LatEnd, LonStart

! TODO: get the grid info from a namelists (98 variables), instead of GITM's UAM.in.  
!       Then remove functions read_in_*.
!       The rest of the UAM.in contents are for running GITM.
!       Can wait until aether_to_dart push is done.
character(len=*), parameter :: filename = 'UAM.in'

character(len=100) :: cLine  ! iCharLen_ == 100
character(len=256) :: fileloc

integer :: i, iunit, ios

character(len=*), parameter :: routine = 'get_grid_info_from_blocks'

! get the ball rolling ...

nBlocksLon = 0
nBlocksLat = 0
nBlocksAlt = 0
LatStart   = 0.0_r8
LatEnd     = 0.0_r8
LonStart   = 0.0_r8

write(fileloc,'(a,''/'',a)') trim(restart_dirname),trim(filename)

if (debug > 4) then
   write(string1,*) 'Now opening Aether UAM file: ',trim(fileloc)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
end if


iunit = open_file(trim(fileloc), action='read')

UAMREAD : do i = 1, 1000000

   read(iunit,'(a)',iostat=ios) cLine

   if (ios /= 0) then
      ! If we get to the end of the file or hit a read error without
      ! finding what we need, die.
      write(string1,*) 'cannot find #GRID in ',trim(fileloc)
      call error_handler(E_ERR,'get_grid_info_from_blocks',string1,source,revision,revdate)
   endif

   if (cLine(1:5) .ne. "#GRID") cycle UAMREAD

   nBlocksLon = read_in_int( iunit,'NBlocksLon',trim(fileloc))
   nBlocksLat = read_in_int( iunit,'NBlocksLat',trim(fileloc))
   nBlocksAlt = read_in_int( iunit,'NBlocksAlt',trim(fileloc))
   LatStart   = read_in_real(iunit,'LatStart',  trim(fileloc))
   LatEnd     = read_in_real(iunit,'LatEnd',    trim(fileloc))
   LonStart   = read_in_real(iunit,'LonStart',  trim(fileloc))

   exit UAMREAD

enddo UAMREAD

if (debug > 4) then
   write(string1,*) 'Successfully read Aether UAM grid file:',trim(fileloc)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   nBlocksLon:',nBlocksLon
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   nBlocksLat:',nBlocksLat
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   nBlocksAlt:',nBlocksAlt
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   LatStart:',LatStart
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   LatEnd:',LatEnd
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   LonStart:',LonStart
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
end if

call close_file(iunit)

end subroutine get_grid_info_from_blocks

!==================================================================

function read_in_int(iunit,varname,filename)

integer,          intent(in) :: iunit
character(len=*), intent(in) :: varname,filename
integer                      :: read_in_int

character(len=100) :: cLine
integer :: i, ios

! Read a line 
read(iunit,'(a)',iostat=ios) cLine
if (ios /= 0) then
   write(string1,*) 'cannot find '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'get_grid_dims',string1,source,revision,revdate)
endif

! Remove anything after a space or TAB
i=index(cLine,' ');     if( i > 0 ) cLine(i:len(cLine))=' '
i=index(cLine,char(9)); if( i > 0 ) cLine(i:len(cLine))=' '

read(cLine,*,iostat=ios)read_in_int

if(ios /= 0) then
   write(string1,*)'unable to read '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'read_in_int',string1,source,revision,revdate,&
             text2=cLine)
endif

end function read_in_int

!=================================================================

function read_in_real(iunit,varname,filename)

integer,          intent(in) :: iunit
character(len=*), intent(in) :: varname,filename
real(r8)                     :: read_in_real

character(len=100) :: cLine
integer :: i, ios

! Read a line 
read(iunit,'(a)',iostat=ios) cLine
if (ios /= 0) then
   write(string1,*) 'cannot find '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'get_grid_dims',string1,source,revision,revdate)
endif

! Remove anything after a space or TAB
i=index(cLine,' ');     if( i > 0 ) cLine(i:len(cLine))=' '
i=index(cLine,char(9)); if( i > 0 ) cLine(i:len(cLine))=' '

! Now that we have a line with nothing else ... parse it
read(cLine,*,iostat=ios)read_in_real

if(ios /= 0) then
   write(string1,*)'unable to read '//trim(varname)//' in '//trim(filename)
   call error_handler(E_ERR,'read_in_real',string1,source,revision,revdate)
endif

end function read_in_real

!=================================================================

! open enough of the restart files to read in the lon, lat, alt arrays

subroutine get_grid_from_blocks(dirname, nBlocksLon, nBlocksLat, nBlocksAlt, &
                  nxPerBlock, nyPerBlock, nzPerBlock,   &
                  lons, lats, alts )

character(len=*), intent(in) :: dirname
integer, intent(in)  :: nBlocksLon ! Number of Longitude blocks
integer, intent(in)  :: nBlocksLat ! Number of Latitude  blocks
integer, intent(in)  :: nBlocksAlt ! Number of Altitude  blocks
integer, intent(out) :: nxPerBlock ! Number of non-halo Longitude centers per block
integer, intent(out) :: nyPerBlock ! Number of non-halo Latitude  centers per block
integer, intent(out) :: nzPerBlock ! Number of Vertical grid centers

real(r8), allocatable , dimension( : ), intent(inout) :: lons, lats, alts

integer :: ios, nb, offset, ncid, nboff
character(len=128) :: filename
real(r4), allocatable :: temp(:,:,:)
integer :: starts(3),ends(3), xcount, ycount, zcount

character(len=*), parameter :: routine = 'get_grid_from_blocks'

! TODO: Here it needs to read the x,y,z  from a NetCDF block file(s),
!       in order to calculate the n[xyz]PerBlock dimensions. 
!       grid_g0000.nc looks like a worthy candidate, but a restart could be used.
write (filename,'(2A)')  trim(dirname),'/grid_g0000.nc'
ncid = nc_open_file_readonly(filename, routine)

! The grid (and restart) file variables have halos, so strip them off
! to get the number of actual data values in each dimension of the block.
nxPerBlock = nc_get_dimension_size(ncid, 'x', routine) - 2*nGhost
nyPerBlock = nc_get_dimension_size(ncid, 'y', routine) - 2*nGhost
nzPerBlock = nc_get_dimension_size(ncid, 'z', routine)

nlon = nBlocksLon * nxPerBlock
nlat = nBlocksLat * nyPerBlock
nalt = nBlocksAlt * nzPerBlock     

write(string1,*)  'nlon = ', nlon
call error_handler(E_MSG,routine,string1,source,revision,revdate)
write(string1,*)  'nlat = ', nlat
call error_handler(E_MSG,routine,string1,source,revision,revdate)
write(string1,*)  'nalt = ', nalt
call error_handler(E_MSG,routine,string1,source,revision,revdate)

! This is also done in gitm's static_init_model, which is not called by aether_to_dart,
! so it's not redundant.
allocate( lons( nlon ))
allocate( lats( nlat ))
allocate( alts( nalt ))

if (debug > 4) then
   write(string1,*) 'Successfully read GITM grid file:',trim(filename)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   nxPerBlock:',nxPerBlock
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   nyPerBlock:',nyPerBlock
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*) '   nzPerBlock:',nzPerBlock
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
endif

! A temp array large enough to hold any of the 3D
! Lon,Lat or Alt arrays from a block plus ghost cells.
! The restart files have C-indexing (fastest changing dim is the last).
allocate(temp( 1:nzPerBlock, &
               1-nGhost:nyPerBlock+nGhost, &
               1-nGhost:nxPerBlock+nGhost))
temp = -888888.

print*,'shape of temp = ',shape(temp)

starts(1) = 1-nGhost
starts(2) = 1-nGhost
starts(3) = 1
ends(1)   = nxPerBlock+nGhost
ends(2)   = nyPerBlock+nGhost
ends(3)   = nzPerBlock
xcount = nxPerBlock + 2*nGhost
ycount = nyPerBlock + 2*nGhost
zcount = nzPerBlock
print*,'starts = ',starts
print*,'ends = ',ends
print*,'counts = ',xcount,ycount,zcount

! go across the south-most block row picking up all longitudes
do nb = 1, nBlocksLon

   filename = block_file_name('grid', -1, nb-1)
   ncid = open_block_file(trim(filename), 'read')

! Read 3D array and extract the longitudes of the non-halo data of this block.
!  This gets nc_get_double_3d, even though the fields are float.
!? Is there some environment setting that says float = double?
! ERROR This yields Start+count exceeds dimension bound
!     call nc_get_variable(ncid, 'Longitude', temp, routine)
! ERROR: this yields Index exceeds dimension bound
! The restart files have C-indexing (fastest changing dim is the last),
! So invert the dimension bounds.
     call nc_get_variable(ncid, 'Longitude', &
          temp(starts(3):ends(3),starts(2):ends(2),starts(1):ends(1)), &
          routine, &
        nc_count=(/zcount,ycount,xcount/))
! Shouldn't need to specify default values         nc_start=(/1,1,1/), &

!           temp(1:zcount,1:ycount,1:xcount), &
!        nc_start=(/starts(1),starts(2),starts(3)/), &
! TODO: nc_get_variable stops on error conditions, does not pass back ios.
!    if ( ios /= 0 ) then
!       print *,'size:',size(temp(1-nGhost:nxPerBlock+nGhost))
!       print *,'IO error code:',ios
!       write(string1,*)'ERROR reading file ', trim(filename)
!       write(string2,*)'longitude block ',nb,' of ',nBlocksLon
!       call error_handler(E_ERR,'get_grid',string1, &
!                  source,revision,revdate,text2=string2)
!    endif

   offset = (nxPerBlock * (nb - 1))
   lons(offset+1:offset+nxPerBlock) = temp(1,1,1:nxPerBlock)

   call nc_close_file(ncid)
enddo

! go up west-most block row picking up all latitudes
do nb = 1, nBlocksLat

   ! TODO; Aether block name counters start with 0, but the lat values can come from 
   !       any lon=const column. 
   nboff = ((nb - 1) * nBlocksLon)
   filename = block_file_name('grid', -1, nboff)
   ncid = open_block_file(trim(filename), 'read')

     call nc_get_variable(ncid, 'Latitude', &
          temp(starts(3):ends(3),starts(2):ends(2),starts(1):ends(1)), &
          routine, nc_count=(/zcount,ycount,xcount/))
        
!    if ( ios /= 0 ) then
!       write(string1,*)'ERROR reading file ', trim(filename)
!       write(string2,*)'latitude block ',nb,' of ',nBlocksLat
!       call error_handler(E_ERR,'get_grid',string1, &
!                  source,revision,revdate,text2=string2)
!    endif

   offset = (nyPerBlock * (nb - 1))
   lats(offset+1:offset+nyPerBlock) = temp(1,1:nyPerBlock,1)

   call nc_close_file(ncid)
enddo


! this code assumes UseTopography is false - that all columns share
! the same altitude array, so we can read it from the first block.
! if this is not the case, this code has to change.

filename = block_file_name('grid', -1, 0)
ncid = open_block_file(trim(filename), 'read')

temp = MISSING_R8
call nc_get_variable(ncid, 'Altitude', &
     temp(starts(3):ends(3),starts(2):ends(2),starts(1):ends(1)), &
     routine, nc_count=(/zcount,ycount,xcount/))

alts(1:nzPerBlock) = temp(1:nzPerBlock,1,1)
! print*,'temp = ',temp(:,1,1)
! print*,'alts = ',alts

call nc_close_file(ncid)

deallocate(temp)

! convert from radians into degrees
lons = lons * RAD2DEG
lats = lats * RAD2DEG

if (debug > 4) then
   print *, 'All lons ', lons
   print *, 'All lats ', lats
   print *, 'All alts ', alts
endif

if ( debug > 1 ) then ! Check dimension limits
   write(string1,*)'LON range ',minval(lons),maxval(lons)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'LAT range ',minval(lats),maxval(lats)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ALT range ',minval(alts),maxval(alts)
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
endif

end subroutine get_grid_from_blocks

!==================================================================

!> Create a filename from input file characteristics: 
!     filetype, member number, block number.
!  filetype = {'grid','neutrals','ions', [...?]}.  
!             The first part of the name of the aether file to read.
!  memnum or blocknum < 0 means don't include that part of the name.

function block_file_name(filetype, memnum, blocknum)

character(len=*), intent(in)  :: filetype  ! one of {grid,ions,neutrals}
! TODO: ? Will this need to open the grid_{below,corners,down,left} filetypes?
!       This code can handle it; a longer filetype passed in, and no member
!       ? output files?
integer,          intent(in)  :: blocknum
integer,          intent(in)  :: memnum
character(len=128) :: block_file_name

block_file_name = trim(filetype)
if (memnum   >= 0) write(block_file_name, '(A,A2,I0.4)') trim(block_file_name), '_m', memnum
if (blocknum >= 0) write(block_file_name, '(A,A2,I0.4)') trim(block_file_name), '_g', blocknum
block_file_name = trim(block_file_name)//'.nc'
! TODO: Convert print to the error handler
print*,'filename, memnum, blocknum = ' ,trim(block_file_name), memnum, blocknum 

end function block_file_name

!==================================================================

!> open the requested restart file and return the ncid

function open_block_file(filename,rw)

character(len=*), intent(in) :: filename
character(len=*), intent(in)  :: rw   ! 'read' or 'readwrite'
integer :: open_block_file

character(len=*), parameter :: routine = 'open_block_file'

if ( rw == 'read' .and. .not. file_exist(trim(filename)) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'open_block_file',string1,source,revision,revdate)
endif

if (debug > 0) then
   write(string1,*) 'Opening file ', trim(filename), ' for ', trim(rw)
   call error_handler(E_MSG,'open_block_file',string1,source,revision,revdate)
end if

open_block_file = nc_open_file_readonly(trim(filename), routine)

if (debug > 80) then
   write(string1,*) 'Returned file descriptor is ', open_block_file
   call error_handler(E_MSG,'open_block_file',string1,source,revision,revdate)
end if

end function open_block_file

!=================================================================

subroutine verify_block_variables( variable_array, ngood)

character(len=*), dimension(:),   intent(in)  :: variable_array
integer,                          intent(out) :: ngood

integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname

character(len=*), parameter :: routine = 'verify_state_variables'

nrows = size(variable_array,1)

ngood = 0
MyLoop : do i = 1, nrows

   varname   = variable_array(i)

   if ( varname  == ' ') exit MyLoop ! Found end of list.

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)
endif

end subroutine verify_block_variables

!==================================================================
!> Converts Aether restart files to a netCDF file
!> Modified from models/gitm/model_mod.f90
!>
!> This routine needs:
!>
!> 1.  A base dirname for the restart files (restart_dirname).
!> they will have the format 'dirname/bNNNN.rst'  where NNNN has
!> leading 0s and is the block number.   Blocks start in the
!> southwest corner of the lat/lon grid and go east first, 
!> then to the west end of the next row north and end in the northeast corner. 
!> The other info is in 'dirname/header.rst'
!> 
!> 2.  The name of the output file to store the netCDF variables 
!> (netcdf_output_file)
!>
!> In the process, the routine will find:
!>
!> 1. The overall grid size, {nlon,nlat,nalt} when you've read in all the blocks. 
!>    (nBlocksLon, nBlocksLat, 1)
!>
!> 2. The number of blocks in Lon and Lat (nBlocksLon, nBlocksLat)
!>
!> 3. The number of lon/lats in a single grid block  (nxPerBlock, 
!>    nyPerBlock, nzPerBlock)
!>
!> 4. The number of neutral species (and probably a mapping between
!>    the species number and the variable name)  (nSpeciesTotal, nSpecies)
!>
!> 5. The number of ion species (ditto - numbers <-> names) (nIons)
!>
!> We assume that the 'UseTopography' flag is false - that all columns
!> have the same altitude arrays.  This is true on earth but not on
!> other planets.
!>
!> In addition to reading in the state data, it fills Longitude,
!> Latitude, and Altitude arrays with the grid spacing.  This grid
!> is orthogonal and rectangular but can have irregular spacing along
!> any or all of the three dimensions.

subroutine restart_files_to_netcdf(restart_dirname, member, netcdf_output_file)

character(len=*), intent(in)  :: restart_dirname
character(len=*), intent(in)  :: netcdf_output_file
integer, intent(in) :: member

integer :: ncid

character(len=*), parameter :: routine = 'restart_files_to_netcdf'

if (module_initialized ) then
    write(string1,*)'The aether static_init_model was already initialized but ',trim(routine),&
      ' uses a separate initialization procedure'
    call error_handler(E_ERR,routine,string1,source,revision,revdate)
end if

call static_init_blocks(restart_dirname)

ncid = nc_create_file(netcdf_output_file)

! DONE: This should probably be replaced by  nc_write_model_atts(ncid).
!       That may require renaming some dimension variables.
! call add_nc_definitions(ncid)
! Enters and exits define mode;
call nc_write_model_atts(ncid, 0)

call get_data(restart_dirname, ncid, member, define=.true.)

! TODO: add_nc_dimvars has not been activated because the functionality is in nc_write_model_atts
!       but maybe it shouldn't be.  Also, we haven't settled on the mechanism for identifying
!       the state vector field names and source.
! call add_nc_dimvars(ncid)

call get_data(restart_dirname, ncid, member, define=.false.)

! TODO: this needs to be updated to write to which file?
! call write_model_time(ncid, state_time)

call nc_close_file(ncid)

end subroutine restart_files_to_netcdf

!==================================================================

subroutine add_nc_definitions(ncid)

integer, intent(in) :: ncid 

call nc_add_global_attribute(ncid, 'model', 'aether')

!-------------------------------------------------------------------------------
! Determine shape of most important namelist
!-------------------------------------------------------------------------------
!
!call find_textfile_dims('gitm_vars.nml', nlines, linelen)
!if (nlines > 0) then
!   has_gitm_namelist = .true.
!
!   allocate(textblock(nlines))
!   textblock = ''
!
!   call nc_define_dimension(ncid, 'nlines',  nlines)
!   call nc_define_dimension(ncid, 'linelen', linelen)
!   call nc_define_character_variable(ncid, 'gitm_in', (/ 'nlines ', 'linelen' /))
!   call nc_add_attribute_to_variable(ncid, 'gitm_in', 'long_name', 'contents of gitm_in namelist')
!
!else
!  has_gitm_namelist = .false.
!endif
!
!----------------------------------------------------------------------------
! output only grid info - state vars will be written by other non-model_mod code
!----------------------------------------------------------------------------

call nc_define_dimension(ncid, LON_DIM_NAME, nlon)
call nc_define_dimension(ncid, LAT_DIM_NAME, nlat)
call nc_define_dimension(ncid, ALT_DIM_NAME, nalt)
! TODO: is WL in Aether?  No; remove from model_mod.
call nc_define_dimension(ncid, 'WL',  1)  ! wavelengths - currently only 1?

!----------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and the Attributes
!----------------------------------------------------------------------------

! TODO: This defines more attributes than TIEGCM.  Prefer?  Are these accurate for Aether?
! Grid Longitudes
call nc_define_double_variable(ncid, LON_VAR_NAME, (/ LON_DIM_NAME /) )
call nc_add_attribute_to_variable(ncid, LON_VAR_NAME, 'type',           'x1d')
call nc_add_attribute_to_variable(ncid, LON_VAR_NAME, 'long_name',      'grid longitudes')
call nc_add_attribute_to_variable(ncid, LON_VAR_NAME, 'cartesian_axis', 'X')
call nc_add_attribute_to_variable(ncid, LON_VAR_NAME, 'units',          'degrees_east')
call nc_add_attribute_to_variable(ncid, LON_VAR_NAME, 'valid_range',     (/ 0.0_r8, 360.0_r8 /) )

! Grid Latitudes
call nc_define_double_variable(ncid, LAT_VAR_NAME, (/ LAT_DIM_NAME /) )
call nc_add_attribute_to_variable(ncid, LAT_VAR_NAME, 'type',           'y1d')
call nc_add_attribute_to_variable(ncid, LAT_VAR_NAME, 'long_name',      'grid latitudes')
call nc_add_attribute_to_variable(ncid, LAT_VAR_NAME, 'cartesian_axis', 'Y')
call nc_add_attribute_to_variable(ncid, LAT_VAR_NAME, 'units',          'degrees_north')
call nc_add_attribute_to_variable(ncid, LAT_VAR_NAME, 'valid_range',     (/ -90.0_r8, 90.0_r8 /) )

! Grid Altitudes
call nc_define_double_variable(ncid, ALT_VAR_NAME, (/ ALT_DIM_NAME /) )
call nc_add_attribute_to_variable(ncid, ALT_VAR_NAME, 'type',           'z1d')
call nc_add_attribute_to_variable(ncid, ALT_VAR_NAME, 'long_name',      'grid altitudes')
call nc_add_attribute_to_variable(ncid, ALT_VAR_NAME, 'cartesian_axis', 'Z')
call nc_add_attribute_to_variable(ncid, ALT_VAR_NAME, 'units',          'meters')
call nc_add_attribute_to_variable(ncid, ALT_VAR_NAME, 'positive',       'up')

! Grid wavelengths
call nc_define_double_variable(ncid, 'WL', (/ 'WL' /) )
call nc_add_attribute_to_variable(ncid, 'WL', 'type',           'x1d')
call nc_add_attribute_to_variable(ncid, 'WL', 'long_name',      'grid wavelengths')
call nc_add_attribute_to_variable(ncid, 'WL', 'cartesian_axis', 'X')
call nc_add_attribute_to_variable(ncid, 'WL', 'units',          'wavelength_index')
call nc_add_attribute_to_variable(ncid, 'WL', 'valid_range',     (/ 0.9_r8, 38.1_r8 /) )

end subroutine add_nc_definitions

!=================================================================
! open all restart files and read in the requested data item

subroutine get_data(dirname, ncid_output, member, define)

character(len=*), intent(in)  :: dirname
integer,          intent(in)  :: ncid_output, member
logical,          intent(in)  :: define

integer :: ibLoop, jbLoop
integer :: ib, jb, nb, iunit

character(len=256) :: filename


if (define) then
   ! if define, run one block.
   ! the read_data_from_block call defines the variables in the whole domain netCDF file.
   ibLoop = 1
   jbLoop = 1
   call nc_begin_define_mode(ncid_output)
else
   ! if not define, run all blocks.
   ! the read_data_from_block call adds the (ib,jb) block to a netCDF variable 
   ! in order to make a file containing the data for all the blocks.
   ibLoop = nBlocksLon
   jbLoop = nBlocksLat
end if

print*,'get_data: define = ',define
do jb = 1, jbLoop
   do ib = 1, ibLoop

      call read_data_from_block(ncid_output, dirname, ib, jb, member, define)

   enddo
enddo

if (define) call nc_end_define_mode(ncid_output)

end subroutine get_data

!==================================================================

!> Open all restart files and read in the requested data items.
!> The unpack* calls will write the data to the filter_input.nc.
!>
!> This is a two-pass method: first run through to define the NC variables
!> in the filter_input.nc (define = .true.),
!> then run again to write the data to the NC file(define = .false.)

subroutine read_data_from_block(ncid_output, dirname, ib, jb, member, define)

integer,  intent(in) :: ncid_output
character(len=*), intent(in)  :: dirname
integer,  intent(in) :: ib, jb
integer,  intent(in) :: member
logical,  intent(in) :: define

real(r4), allocatable :: temp1d(:), temp2d(:,:), temp3d(:,:,:)
real(r4), allocatable :: alt1d(:), density_ion_e(:,:,:)
real(r4) :: temp0d !Alex: single parameter has "zero dimensions"
integer :: i, j, maxsize, ivar, nb, ncid_input
integer :: block(2) = 0

logical :: no_idensity

character(len=*), parameter :: routine = 'read_data_from_block'
character(len=128) :: file_root 
character(len=256) :: filename
character(len=NF90_MAX_NAME) :: varname

block(1) = ib
block(2) = jb
! The block number, as counted in Aether.
! Lower left is 0, increase to the East, then 1 row farther north, West to East.
nb = (jb-1) * nBlocksLon + ib - 1

! a temp array large enough to hold any of the
! Lon,Lat or Alt array from a block plus ghost cells
allocate(temp1d(1-nGhost:max(nxPerBlock,nyPerBlock,nzPerBlock)+nGhost))

! treat alt specially since we want to derive TEC here
! TODO: See density_ion_e too.
allocate( alt1d(1-nGhost:max(nxPerBlock,nyPerBlock,nzPerBlock)+nGhost))

! temp array large enough to hold any 2D field 
allocate(temp2d(1-nGhost:nyPerBlock+nGhost, &
                1-nGhost:nxPerBlock+nGhost))

! TODO: We need all altitudes, but there might be vertical blocks in the future.
!       But there would be no vertical halos.
!       Make nzcount adapt to whether there are blocks.
!       And temp needs to have C-ordering, which is what the restart files have.
! temp array large enough to hold 1 species, temperature, etc
allocate(temp3d(1:nzPerBlock, &
                1-nGhost:nyPerBlock+nGhost, &
                1-nGhost:nxPerBlock+nGhost))

! save density_ion_e to compute TEC
allocate(density_ion_e(1:nzPerBlock, &
                       1-nGhost:nyPerBlock+nGhost, &
                       1-nGhost:nxPerBlock+nGhost))

! Aether gives a unique name to each (of 6) velocity components
! ! temp array large enough to hold velocity vect, etc
! maxsize = max(3, nSpecies)
! allocate(temp4d(1-nGhost:nxPerBlock+nGhost, &
!                 1-nGhost:nyPerBlock+nGhost, &
!                 1-nGhost:nzPerBlock+nGhost, maxsize))


! TODO; Does Aether need a replacement for these Density fields?  Yes.
!       But they are probably read by the loops below.
!       Don't need to fetch index because Aether has NetCDF restarts,
!       so just loop over the field names to read.
! Read the index from the first species
! call get_index_from_gitm_varname('NDensityS', inum, ivals)

! if (inum > 0) then
!    ! if i equals ival, use the data from the state vect
!    ! otherwise read/write what's in the input file
!    j = 1
!    do i = 1, nSpeciesTotal
!       if (debug > 80) then
!          write(string1,'(A,I0,A,I0,A,I0,A,I0,A)') 'Now reading species ',i,' of ',nSpeciesTotal, &
!             ' for block (',ib,',',jb,')' 
!          call error_handler(E_MSG,routine,string1,source,revision,revdate)
!       end if
!       read(iunit)  temp3d
!       if (j <= inum) then
!          if (i == gitmvar(ivals(j))%gitm_index) then
!             call unpack_data(temp3d, ivals(j), block, ncid, define)
!             j = j + 1
!          endif
!       endif
!    enddo
! else
!    if (debug > 80) then
!       write(string1,'(A)') 'Not writing the NDensityS variables to file'
!       call error_handler(E_MSG,routine,string1,source,revision,revdate)
!    end if
!    ! nothing at all from this variable in the state vector.
!    ! copy all data over from the input file to output file
!    do i = 1, nSpeciesTotal
!       read(iunit)  temp3d
!    enddo
! endif
! 
! call get_index_from_gitm_varname('IDensityS', inum, ivals)
! 
! ! assume we could not find the electron density for VTEC calculations
! no_idensity = .true.
! 
! if (inum > 0) then
!    ! one or more items in the state vector need to replace the
!    ! data in the output file.  loop over the index list in order.
!    j = 1
! ! TODO:   electron density is not in the restart files, but it's needed for TEC
!           In Aether they will be from an ions file, but now only from an output file (2023-10-30).
!    do i = 1, nIons
!       if (debug > 80) then
!          write(string1,'(A,I0,A,I0,A,I0,A,I0,A)') 'Now reading ion ',i,' of ',nIons, &
!             ' for block (',ib,',',jb,')' 
!          call error_handler(E_MSG,routine,string1,source,revision,revdate)
!       end if
!       read(iunit)  temp3d
!       if (j <= inum) then
!          if (i == gitmvar(ivals(j))%gitm_index) then
!             ! ie_, the gitm index for electron density, comes from ModEarth 
!             if (gitmvar(ivals(j))%gitm_index == ie_) then
!                ! save the electron density for TEC computation
!                density_ion_e(:,:,:) = temp3d(:,:,:)
!                no_idensity = .false.
!             end if
!             ! read from input but write from state vector
!             call unpack_data(temp3d, ivals(j), block, ncid, define)
!             j = j + 1
!          endif
!       endif
!    enddo
! else
!    ! nothing at all from this variable in the state vector.
!    ! read past this variable
!    if (debug > 80) then
!       write(string1,'(A)') 'Not writing the IDensityS variables to file'
!       call error_handler(E_MSG,routine,string1,source,revision,revdate)
!    end if
!    do i = 1, nIons
!       read(iunit)  temp3d
!    enddo
! endif

! Handle the 2 restart file types (ions and neutrals).
! Each field has a file type associated with it: variable_table(f_index,VT_ORIGININDX)
! TODO: for now require that all neutrals are listed in variable_table before the ions.

file_root = variable_table(1,VT_ORIGININDX)
filename = block_file_name(file_root, member, nb)
ncid_input = open_block_file(trim(filename), 'read')

print*,'read_data_from_block: nfields_neutral = ',nfields_neutral
do ivar = 1, nfields_neutral
   write(varname,'(A)') trim(variable_table(ivar,VT_VARNAMEINDX))

   ! TODO: Given the subroutine name, perhaps these definition sections should be 
   !       one call higher up, with the same loop around it.
   if (define) then
   ! Define the variable in the filter_input.nc file (the output from this program).
   ! The calling routine entered define mode.

      if (debug > 10) then 
         write(string1,'(A,I0,2A)') 'Defining ivar = ', ivar,':',varname
         call error_handler(E_MSG,routine,string1,source,revision,revdate)
      end if
   
      call nc_define_real_variable(ncid_output, varname, &
           (/ ALT_DIM_NAME, LAT_DIM_NAME, LON_DIM_NAME /) )
      print*,routine,': defined ivar, varname = ', ivar, varname 
! TODO: does the filter_input.nc file need all these attributes?  TIEGCM doesn't add them.
      !    They are not available from the restart files.
      !    Add them to the ions section too.
      ! call nc_add_attribute_to_variable(ncid, varname, 'long_name',    gitmvar(ivar)%long_name)
      ! call nc_add_attribute_to_variable(ncid, varname, 'units',        gitmvar(ivar)%units)
      ! !call nc_add_attribute_to_variable(ncid, varname, 'storder',     gitmvar(ivar)%storder)
      ! call nc_add_attribute_to_variable(ncid, varname, 'gitm_varname', gitmvar(ivar)%gitm_varname)
      ! call nc_add_attribute_to_variable(ncid, varname, 'gitm_dim',     gitmvar(ivar)%gitm_dim)
      ! call nc_add_attribute_to_variable(ncid, varname, 'gitm_index',   gitmvar(ivar)%gitm_index)


   else if (file_root == 'neutrals') then
   ! Read 3D array and extract the non-halo data of this block.
! TODO: There are no 2D or 1D fields in ions or neutrals, but there could be; different temp array.
      call nc_get_variable(ncid_input, varname, temp3d, routine)
      print*,'read_data_from_block: temp3d = ',temp3d(1,1,1),temp3d(15,15,15),variable_table(ivar,VT_VARNAMEINDX)
      print*,'read_data_from_block: define = ',define
      call unpack_data(temp3d, ivar, block, ncid_output, define)
   else
      write(string1,*) 'Trying to read neutrals, but variable_table(',ivar,VT_ORIGININDX, &
                       ') /= "neutrals"'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

enddo
call nc_close_file(ncid_input)

file_root = variable_table(nfields_neutral+1,VT_ORIGININDX)
filename = block_file_name(file_root, member, nb)
ncid_input = open_block_file(trim(filename), 'read')

print*,'read_data_from_block: nfields_ion = ',nfields_ion
do ivar = nfields_neutral +1,nfields_neutral + nfields_ion
   write(varname,'(A)') trim(variable_table(ivar,VT_VARNAMEINDX))

   if (define) then

      if (debug > 10) then 
         write(string1,'(A,I0,2A)') 'Defining ivar = ', ivar,':',varname
         call error_handler(E_MSG,routine,string1,source,revision,revdate)
      end if
   
      call nc_define_real_variable(ncid_output, varname, &
           (/ ALT_DIM_NAME, LAT_DIM_NAME, LON_DIM_NAME /) )
      print*,routine,': defined ivar, varname = ', ivar, varname 

   else if (file_root == 'ions') then
      call nc_get_variable(ncid_input, varname, temp3d, routine)
      call unpack_data(temp3d, ivar, block, ncid_output, define)
   else
      write(string1,*) 'Trying to read ions, but variable_table(',ivar,VT_ORIGININDX, &
                       ') /= "ions"'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

enddo
call nc_close_file(ncid_input)

! TODO: Does Aether need TEC to be calculated? Yes
! ! add the VTEC as an extended-state variable
! ! NOTE: This variable will *not* be written out to the GITM blocks to netCDF program
! call get_index_from_gitm_varname('TEC', inum, ivals)
! 
! if (inum > 0 .and. no_idensity) then
!    write(string1,*) 'Cannot compute the VTEC without the electron density'
!    call error_handler(E_ERR,routine,string1,source,revision,revdate)
! end if
! 
! if (inum > 0) then
!    if (.not. define) then
!       temp2d = 0._r8
!       ! comptue the TEC integral
!       do i =1,nzPerBlock-1 ! approximate the integral over the altitude as a sum of trapezoids
!          ! area of a trapezoid: A = (h2-h1) * (f2+f1)/2
!          temp2d(:,:) = temp2d(:,:) + ( alt1d(i+1)-alt1d(i) )  * ( density_ion_e(:,:,i+1)+density_ion_e(:,:,i) ) /2.0_r8
!       end do  
!       ! convert temp2d to TEC units
!       temp2d = temp2d/1e16_r8
!    end if
!    call unpack_data2d(temp2d, ivals(1), block, ncid, define) 
! end if

! TODO: Does Aether need f10_7 to be calculated or processed? Yes
! read(iunit)  temp0d
! !gitm_index = get_index_start(domain_id, 'VerticalVelocity')
! call get_index_from_gitm_varname('f107', inum, ivals)
! if (inum > 0) then
!   call unpack_data0d(temp0d, ivals(1), ncid, define) !see comments in the body of the subroutine
! endif
! 
! read(iunit)  temp3d
! call get_index_from_gitm_varname('Rho', inum, ivals)
! if (inum > 0) then
!    call unpack_data(temp3d, ivals(1), block, ncid, define)
! endif

!print *, 'calling dealloc'
deallocate(temp1d, temp2d, temp3d)
deallocate(alt1d, density_ion_e)

end subroutine read_data_from_block

!==================================================================

!> TODO: Activate f10_7 code?
! !> put the f107 estimate (a scalar, hence 0d) into the state vector.
! !> Written specifically
! !> for f107 since f107 is the same for all blocks. So what it does
! !> is take f107 from the first block (block = 0) and disregard
! !> f107 values from all other blocks (hopefully they are the same).
! !> written by alex
! 
! subroutine unpack_data0d(data0d, ivar, ncid, define)
! 
! real(r8), intent(in)    :: data0d
! integer,  intent(in)    :: ivar         ! index into state structure
! integer,  intent(in)    :: ncid
! logical,  intent(in)    :: define
! 
! 
! character(len=*), parameter :: routine = 'unpack_data0d'
! 
! if (define) then
!   
!    if (debug > 10) then 
!       write(string1,'(A,I0,2A)') 'Defining ivar = ', ivar,':',trim(gitmvar(ivar)%varname)
!       call error_handler(E_MSG,routine,string1,source,revision,revdate)
!    end if
! 
!    call nc_define_double_scalar(ncid,   gitmvar(ivar)%varname)
!    call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'long_name',      gitmvar(ivar)%long_name)
!    call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'units',          gitmvar(ivar)%units)
!    call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_varname',   gitmvar(ivar)%gitm_varname)
!    call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_dim',       gitmvar(ivar)%gitm_dim)
!    call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_index',     gitmvar(ivar)%gitm_index)
! 
! else
! 
!    call nc_put_variable(ncid, gitmvar(ivar)%varname, data0d, context=routine)
! 
! end if
! 
! end subroutine unpack_data0d
! 
! !==================================================================
! 
! ! put the requested data into a netcdf variable
! 
! subroutine unpack_data2d(data2d, ivar, block, ncid, define)
! 
! real(r8), intent(in)    :: data2d(1-nGhost:nxPerBlock+nGhost, &
!                                   1-nGhost:nyPerBlock+nGhost)
! 
! integer,  intent(in)    :: ivar         ! variable index
! integer,  intent(in)    :: block(2)
! integer,  intent(in)    :: ncid
! logical,  intent(in)    :: define
! 
! integer :: ib, jb
! integer :: starts(2)
! character(len=*), parameter :: routine = 'unpack_data2d'
! 
! if (define) then
!   
!    if (debug > 10) then 
!       write(string1,'(A,I0,2A)') 'Defining ivar = ', ivar,':',trim(gitmvar(ivar)%varname)
!       call error_handler(E_MSG,routine,string1,source,revision,revdate)
!    end if
! 
!    call nc_define_double_variable(ncid, gitmvar(ivar)%varname, (/ LON_DIM_NAME, LAT_DIM_NAME /) )
!    call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'long_name',      gitmvar(ivar)%long_name)
!    call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'units',          gitmvar(ivar)%units)
!    !call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'storder',        gitmvar(ivar)%storder)
!    call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_varname',   gitmvar(ivar)%gitm_varname)
!    call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_dim',       gitmvar(ivar)%gitm_dim)
!    call nc_add_attribute_to_variable(ncid, gitmvar(ivar)%varname, 'gitm_index',     gitmvar(ivar)%gitm_index)
! 
! else
!    ib = block(1)
!    jb = block(2)
! 
!    ! to compute the start, consider (ib-1)*nxPerBlock+1
!    starts(1) = (ib-1)*nxPerBlock+1
!    starts(2) = (jb-1)*nyPerBlock+1
! 
!    call nc_put_variable(ncid, gitmvar(ivar)%varname, &
!       data2d(1:nxPerBlock,1:nyPerBlock), &
!       context=routine, nc_start=starts, &
!       nc_count=(/nxPerBlock,nyPerBlock/))
! end if
! 
! end subroutine unpack_data2d

!==================================================================

! put the requested data into a netcdf variable

subroutine unpack_data(data3d, ivar, block, ncid, define)

real(r4), intent(in)    :: data3d(1:nzPerBlock, &
                                  1-nGhost:nyPerBlock+nGhost, &
                                  1-nGhost:nxPerBlock+nGhost)

integer,  intent(in)    :: ivar         ! variable index
integer,  intent(in)    :: block(2)
integer,  intent(in)    :: ncid

integer :: ib, jb
integer :: starts(3)
character(len=*), parameter :: routine = 'unpack_data'
character(len=NF90_MAX_NAME) :: varname

print*,'unpack_data: data3d = ',data3d(1,1,1),data3d(15,15,15)
print*,'unpack_data: define = ',define

write(varname,'(A)') trim(variable_table(ivar,VT_VARNAMEINDX))

ib = block(1)
jb = block(2)

! to compute the start, consider (ib-1)*nxPerBlock+1
starts(1) = 1
starts(2) = (jb-1)*nyPerBlock+1
starts(3) = (ib-1)*nxPerBlock+1

call nc_put_variable(ncid, varname, &
   data3d(1:nzPerBlock,1:nyPerBlock,1:nxPerBlock), &
   context=routine, nc_start=starts, &
   nc_count=(/nzPerBlock,nyPerBlock,nxPerBlock/))
print*,'unpack_data: filled varname = ', varname 

end subroutine unpack_data


!=================================================================
!> sort list x into order based on values in list.
!> should only be called on short ( < hundreds) of values or will be slow
!> @todo FIXME this should be using the sort module routine instead.

subroutine sortindexlist(list, x, inum)

integer, intent(inout) :: list(:)
integer, intent(inout) :: x(:)
integer, intent(in)    :: inum

integer :: tmp
integer :: j, k

!  DO A N^2 SORT - only use for short lists
do j = 1, inum - 1
   do k = j + 1, inum
      ! if list() is in wrong order, exchange both list items and
      ! items in x array.
      if(list(j) .gt. list(k)) then
         tmp = list(k)
         list(k) = list(j)
         list(j) = tmp
         tmp = x(k)
         x(k) = x(j)
         x(j) = tmp
      end if
   end do
end do
end subroutine sortindexlist


!-------------------------------------------------------------------------------

function get_model_size()
! Returns the size of the model as an integer.

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size

!==================================================================
! TODO; will be provided by Ben's model_mod.
! 
 subroutine model_interpolate(state_handle, ens_size, location, iqty, obs_val, istatus)
 ! Given a location, and a model state variable qty,
 ! interpolates the state variable field to that location.
 ! obs_val is the interpolated value for each ensemble member
 ! istatus is the success (0) or failure of the interpolation
 
 type(ensemble_type), intent(in) :: state_handle
 integer,             intent(in) :: ens_size
 type(location_type), intent(in) :: location
 integer,             intent(in) :: iqty
 real(r8),           intent(out) :: obs_val(ens_size) !< array of interpolated values
 integer,            intent(out) :: istatus(ens_size)
 
 integer  :: which_vert
 integer  :: lat_below, lat_above, lon_below, lon_above ! these are indices
 real(r8) :: lon_fract, lat_fract
 real(r8) :: lon, lat, lon_lat_lev(3)
 real(r8), dimension(ens_size) :: val11, val12, val21, val22
 real(r8) :: height
 integer  :: level, bogus_level
 integer  :: dom_id, var_id
! 
! if ( .not. module_initialized ) call static_init_model
! 
! ! Default for failure return
! istatus(:) = 1
! obs_val(:) = MISSING_R8
! 
! ! Failure codes
! ! 11 QTY_GEOPOTENTIAL_HEIGHT is unsupported
! ! 22 unsupported veritcal coordinate
! ! 33 level given < or > model levels
! ! 44 quantity not part of the state
! ! 55 outside state (can not extrapolate above or below)
! ! 66 unknown vertical stagger
! 
! ! GITM uses a vtec routine in obs_def_upper_atm_mod:get_expected_gnd_gps_vtec()
! ! TIEGCM has its own vtec routine, so we should use it. This next block ensures that.
! ! The get_expected_gnd_gps_vtec() tries to interpolate QTY_GEOPOTENTIAL_HEIGHT
! ! when it does, this will kill it. 
! 
! if ( iqty == QTY_GEOPOTENTIAL_HEIGHT ) then
!    istatus(:) = 11
!    write(string1,*)'QTY_GEOPOTENTIAL_HEIGHT currently unsupported'
!    call error_handler(E_ERR,'model_interpolate',string1,source, revision, revdate)
! endif
! 
! 
! ! Get the position
! lon_lat_lev = get_location(location)
! lon         = lon_lat_lev(1) ! degree
! lat         = lon_lat_lev(2) ! degree
! height      = lon_lat_lev(3) ! level (int) or height (real)
! level       = int(lon_lat_lev(3))
! 
! 
! which_vert = nint(query_location(location))
! 
! call compute_bracketing_lat_indices(lat, lat_below, lat_above, lat_fract)
! call compute_bracketing_lon_indices(lon, lon_below, lon_above, lon_fract)
! 
! ! Pressure is not part of the state vector
! ! pressure is static data on plevs/pilevs
! if ( iqty == QTY_PRESSURE) then
!    if (which_vert == VERTISLEVEL) then
!       ! @todo from Lanai code:
!       !   Some variables need plevs, some need pilevs
!       !   We only need the height (aka level)
!       !   the obs_def_upper_atm_mod.f90:get_expected_O_N2_ratio routines queries
!       !   for the pressure at the model levels - EXACTLY - so ...
!       !   FIXME ... at present ... the only time model_interpolate
!       !   gets called with QTY_PRESSURE is to calculate density, which
!       !   requires other variables that only live on the midpoints.
!       !   I cannot figure out how to generically decide when to
!       !   use plevs vs. pilevs
! 
!       ! Check to make sure vertical level is possible.
!       if ((level < 1) .or. (level > nalt)) then
!          istatus(:) = 33
!          return
!       else
!          obs_val(:) = plevs(level)
!          istatus(:) = 0
!       endif
!    elseif (which_vert == VERTISHEIGHT) then
! 
!       ! @todo from Lanai code:
!       !   FIXME ... is it possible to try to get a pressure with which_vert == undefined
!       !   At present, vert_interp will simply fail because height is a negative number.
!       !   @todo HK what are you supposed to do for pressure with VERTISUNDEF? level 1?
! 
!       call vert_interp(state_handle, ens_size, dom_id, var_id, lon_below, lat_below, height, iqty, val11, istatus)
!       if (any(istatus /= 0)) return  ! bail at the first failure
!       call vert_interp(state_handle, ens_size, dom_id, var_id, lon_below, lat_above, height, iqty, val12, istatus)
!       if (any(istatus /= 0)) return
!       call vert_interp(state_handle, ens_size, dom_id, var_id, lon_above, lat_below, height, iqty, val21, istatus)
!       if (any(istatus /= 0)) return
!       call vert_interp(state_handle, ens_size, dom_id, var_id, lon_above, lat_above, height, iqty, val22, istatus)
!       obs_val(:) = interpolate(ens_size, lon_fract, lat_fract, val11, val12, val21, val22)
!    else
! 
!       write(string1,*) 'vertical coordinate type:',which_vert,' cannot be handled'
!       call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
! 
!    endif ! which vert
! 
!    return
! 
! endif ! end of QTY_PRESSURE
! 
! 
! if ( iqty == QTY_VERTICAL_TEC ) then ! extrapolate vtec
! 
!    call extrapolate_vtec(state_handle, ens_size, lon_below, lat_below, val11)
!    call extrapolate_vtec(state_handle, ens_size, lon_below, lat_above, val11)
!    call extrapolate_vtec(state_handle, ens_size, lon_above, lat_below, val11)
!    call extrapolate_vtec(state_handle, ens_size, lon_above, lat_above, val11)
!    obs_val(:) = interpolate(ens_size, lon_fract, lat_fract, val11, val12, val21, val22)
!    istatus(:) = 0
! 
!    return
! endif
! 
! ! check if qty is in the state vector
! call find_qty_in_state(iqty, dom_id, var_id)
! if (dom_id < 0 ) then
!    istatus(:) = 44
!    return
! endif
! 
! if( which_vert == VERTISHEIGHT ) then
! 
!    call vert_interp(state_handle, ens_size, dom_id, var_id, lon_below, lat_below, height, iqty, val11, istatus)
!    if (any(istatus /= 0)) return  ! bail at the first failure
!    call vert_interp(state_handle, ens_size, dom_id, var_id, lon_below, lat_above, height, iqty, val12, istatus)
!    if (any(istatus /= 0)) return
!    call vert_interp(state_handle, ens_size, dom_id, var_id, lon_above, lat_below, height, iqty, val21, istatus)
!    if (any(istatus /= 0)) return
!    call vert_interp(state_handle, ens_size, dom_id, var_id, lon_above, lat_above, height, iqty, val22, istatus)
!    obs_val(:) = interpolate(ens_size, lon_fract, lat_fract, val11, val12, val21, val22)
!    istatus = 0
! elseif( which_vert == VERTISLEVEL) then
!    ! Check to make sure vertical level is possible.
!    if ((level < 1) .or. (level > nilev)) then
!      istatus(:) = 33
!      return
!    endif
! 
!    ! one use of model_interpolate is to allow other modules/routines
!    ! the ability to 'count' the model levels. To do this, create observations
!    ! with locations on model levels and 'interpolate' for QTY_GEOMETRIC_HEIGHT.
!    ! When the interpolation fails, you've gone one level too far.
!    ! HK why does it have to be QTY_GEOMETRIC_HEIGHT?
! 
!    val11(:) = get_state(get_dart_vector_index(lon_below, lat_below, level, domain_id(dom_id), var_id ), state_handle)
!    val12(:) = get_state(get_dart_vector_index(lon_below, lat_above, level, domain_id(dom_id), var_id ), state_handle)
!    val21(:) = get_state(get_dart_vector_index(lon_above, lat_below, level, domain_id(dom_id), var_id ), state_handle)
!    val22(:) = get_state(get_dart_vector_index(lon_above, lat_above, level, domain_id(dom_id), var_id ), state_handle)
!    obs_val(:) = interpolate(ens_size, lon_fract, lat_fract, val11, val12, val21, val22)
!    istatus = 0
! 
! elseif( which_vert == VERTISUNDEF) then
!    bogus_level  = 1  !HK what should this be?  Do only 2D fields have VERTISUNDEF?
!    val11(:) = get_state(get_dart_vector_index(lon_below, lat_below, bogus_level, domain_id(dom_id), var_id), state_handle)
!    val12(:) = get_state(get_dart_vector_index(lon_below, lat_above, bogus_level, domain_id(dom_id), var_id), state_handle)
!    val21(:) = get_state(get_dart_vector_index(lon_above, lat_below, bogus_level, domain_id(dom_id), var_id), state_handle)
!    val22(:) = get_state(get_dart_vector_index(lon_above, lat_above, bogus_level, domain_id(dom_id), var_id), state_handle)
!    obs_val(:) = interpolate(ens_size, lon_fract, lat_fract, val11, val12, val21, val22)
!    istatus(:) = 0
! 
! else
! 
!    write(string1,*) 'vertical coordinate type:',which_vert,' cannot be handled'
!    call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
! 
! endif
! 
 end subroutine model_interpolate

!-------------------------------------------------------------------------------
function shortest_time_between_assimilations()
type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations

!==================================================================
! 
 subroutine get_state_meta_data(index_in, location, var_qty)
 ! Given an integer index into the state vector, returns the
 ! associated location and optionally the variable quantity.
 
 integer(i8),         intent(in)  :: index_in
 type(location_type), intent(out) :: location
 integer, optional,   intent(out) :: var_qty
 
 integer  :: lon_index, lat_index, lev_index
 integer  :: local_qty, var_id, dom_id
 integer  :: seconds, days ! for f10.7 location
 real(r8) :: longitude     ! for f10.7 location
 character(len=NF90_MAX_NAME) :: dim_name
 
! if ( .not. module_initialized ) call static_init_model
! 
! call get_model_variable_indices(index_in, lon_index, lat_index, lev_index, var_id=var_id, dom_id=dom_id, kind_index=local_qty)
! 
! if(present(var_qty)) var_qty = local_qty
! 
! if (get_variable_name(dom_id, var_id) == 'f10_7') then
!    ! f10_7 is most accurately located at local noon at equator.
!    ! 360.0 degrees in 86400 seconds, 43200 secs == 12:00 UTC == longitude 0.0
! 
!    call get_time(state_time, seconds, days)
!    longitude = 360.0_r8 * real(seconds,r8) / 86400.0_r8 - 180.0_r8
!    if (longitude < 0.0_r8) longitude = longitude + 360.0_r8
!    location = set_location(longitude, 0.0_r8,  400000.0_r8, VERTISUNDEF)
!    return
! end if
! 
! ! search for either ilev or lev
! dim_name = ilev_or_lev(dom_id, var_id)
! 
! select case (trim(dim_name))
!    case ('ilev')
!       location  = set_location(lons(lon_index), lats(lat_index), ilevs(lev_index), VERTISLEVEL)
!    case (ALT_DIM_NAME)
!       location  = set_location(lons(lon_index), lats(lat_index), alts(lev_index), VERTISLEVEL)
!    case default
!     call error_handler(E_ERR, 'get_state_meta_data', 'expecting ilev or ilat dimension')
!     ! HK @todo 2D variables.
! end select
! 
 end subroutine get_state_meta_data
 
!==================================================================

subroutine end_model()
! Does any shutdown and clean-up needed for model.

end subroutine end_model

!==================================================================

! Writes the model-specific attributes to a netCDF file.
subroutine nc_write_model_atts( ncid, dom_id)

integer, intent(in)  :: ncid      ! netCDF file identifier
integer, intent(in)  :: dom_id

real(r8), allocatable :: temp_lons(:)
character(len=*), parameter :: routine = 'nc_write_model_atts'

if ( .not. module_initialized ) call static_init_model

! Write Global Attributes

call nc_add_global_creation_time(ncid, routine)

call nc_add_global_attribute(ncid, "model_source", source, routine)
call nc_add_global_attribute(ncid, "model", "Aether", routine)


! define grid dimensions
call nc_define_dimension(ncid, LON_DIM_NAME,  nlon,  routine)
call nc_define_dimension(ncid, LAT_DIM_NAME,  nlat,  routine)
call nc_define_dimension(ncid, ALT_DIM_NAME,  nalt,  routine)
call nc_define_dimension(ncid, 'ilev', nilev, routine)

! define grid variables
! longitude
call nc_define_real_variable(     ncid, LON_DIM_NAME, (/ LON_DIM_NAME /), routine)
call nc_add_attribute_to_variable(ncid, LON_DIM_NAME, 'long_name', 'geographic longitude (-west, +east)',  routine)
call nc_add_attribute_to_variable(ncid, LON_DIM_NAME, 'units', 'degrees_east', routine)

! latitude
call nc_define_real_variable(     ncid, LAT_DIM_NAME, (/ LAT_DIM_NAME /),  routine)
call nc_add_attribute_to_variable(ncid, LAT_DIM_NAME, 'long_name', 'geographic latitude (-south, +north)', routine)
call nc_add_attribute_to_variable(ncid, LAT_DIM_NAME, 'units',     'degrees_north', routine)

! alts
call nc_define_real_variable(     ncid, ALT_DIM_NAME, (/ ALT_DIM_NAME /), routine)
call nc_add_attribute_to_variable(ncid, ALT_DIM_NAME, 'long_name',      'midpoint altitudes', routine)
! DONE: vert coord is altitude, not ...
call nc_add_attribute_to_variable(ncid, ALT_DIM_NAME, 'short name',     'altitude', routine)
call nc_add_attribute_to_variable(ncid, ALT_DIM_NAME, 'positive',       'up', routine)
call nc_add_attribute_to_variable(ncid, ALT_DIM_NAME, 'standard_name',  'unknown', routine)
! call nc_add_attribute_to_variable(ncid, ALT_DIM_NAME, 'formula_terms',  'p0: p0 lev: lev', routine)
! call nc_add_attribute_to_variable(ncid, ALT_DIM_NAME, 'formula',  'p(k) = p0 * exp(-lev(k))', routine)


! ilevs
! call nc_define_real_variable(     ncid, 'ilev', (/ 'ilev' /), routine)
! call nc_add_attribute_to_variable(ncid, 'ilev', 'long_name',      'interface levels', routine)
! call nc_add_attribute_to_variable(ncid, 'ilev', 'short name',     'ln(p0/p)', routine)
! call nc_add_attribute_to_variable(ncid, 'ilev', 'positive',       'up', routine)
! call nc_add_attribute_to_variable(ncid, 'ilev', 'standard_name',  'atmosphere_ln_pressure_coordinate', routine)
! call nc_add_attribute_to_variable(ncid, 'ilev', 'formula_terms',  'p0: p0 lev: ilev', routine)
! ! TODO: Is there an interface alt?
! call nc_add_attribute_to_variable(ncid, ALT_DIM_NAME,  'formula',         'p(k) = p0 * exp(-ilev(k))', routine)


call nc_end_define_mode(ncid, routine)

!-------------------------------------------------------------------------------
! Write variables
!-------------------------------------------------------------------------------

! TODO: Should nc_write_model_atts write dimension contents, not just atts?
! Gitm had a separate routine for filling the dimensions:
! - - - - - - - - - - -
! subroutine add_nc_dimvars(ncid)
! 
! integer, intent(in) :: ncid 
! 
! !----------------------------------------------------------------------------
! ! Fill the coordinate variables
! !----------------------------------------------------------------------------
! 
! call nc_put_variable(ncid, LON_VAR_NAME, lons)
! call nc_put_variable(ncid, LAT_VAR_NAME, lats)
! call nc_put_variable(ncid, ALT_VAR_NAME, alts)
! ! what about WL?
! 
! !if (has_gitm_namelist) then
! !   call file_to_text('gitm_vars.nml', textblock)
! !   call nc_put_variable(ncid, 'gitm_in', textblock)
! !   deallocate(textblock)
! !endif
! 
! !-------------------------------------------------------------------------------
! ! Flush the buffer and leave netCDF file open
! !-------------------------------------------------------------------------------
! call nc_synchronize_file(ncid)
! 
! end subroutine add_nc_dimvars
! - - - - - - - - - - -


! Fill in the coordinate variables

! longitude - Aether uses values +/- pi, but lons has been converted already.
!             DART uses values [0,360]
allocate(temp_lons(nlon))
temp_lons = lons
where (temp_lons < 0.0_r8) temp_lons = temp_lons + 360.0_r8
! where (temp_lons >= 180.0_r8) temp_lons = temp_lons - 360.0_r8
call nc_put_variable(ncid, LON_VAR_NAME,  temp_lons,  routine)
call nc_put_variable(ncid, LAT_VAR_NAME,  lats,  routine)
call nc_put_variable(ncid, ALT_VAR_NAME,  alts,   routine)
! call nc_put_variable(ncid, 'ilev', ilevs,  routine)
deallocate(temp_lons)

! flush any pending i/o to disk
call nc_synchronize_file(ncid, routine)

end subroutine nc_write_model_atts

!==================================================================

! TODO: this will be replaced by Ben.
! Vertical localization is done only in height (ZG).
! obs vertical location is given in height (model_interpolate).
! state vertical location is given in height.
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, state_handle)

type(get_close_type),          intent(in)     :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)     :: base_type, loc_qtys(:)
integer(i8),                   intent(in)     :: loc_indx(:)
integer,                       intent(out)    :: num_close, close_ind(:)
real(r8),            optional, intent(out)    :: dist(:)
type(ensemble_type), optional, intent(in)     :: state_handle

integer :: k, q_ind
integer :: n
integer :: istatus

! n = size(locs)
! 
! if (vertical_localization_on()) then ! need to get height
!   call convert_vertical_state(state_handle, n, locs, loc_qtys, loc_indx, VERTISHEIGHT, istatus)  ! HK Do we care about istatus?
! endif
! 
! call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
!                        num_close, close_ind, dist)
! 
! ! Make the ZG part of the state vector far from everything so it does not get updated.
! ! HK Note if you have inflation on ZG has been inflated.
! ! Scroll through all the obs_loc(:) and obs_kind(:) elements
! 
! do k = 1,num_close
!    q_ind  = close_ind(k)
!    if (loc_qtys(q_ind) == QTY_GEOMETRIC_HEIGHT) then
!       if (do_output() .and. (debug > 99)) then
!          write(     *     ,*)'get_close_state ZG distance is ', &
!                      dist(k),' changing to ',10.0_r8 * PI
!          write(logfileunit,*)'get_close_state ZG distance is ', &
!                      dist(k),' changing to ',10.0_r8 * PI
!       endif
!       dist(k) = 10.0_r8 * PI
!    endif
! enddo
! 
! 
! if (estimate_f10_7) then
! ! f10_7 is given a location of latitude 0.0 and the longitude
! ! of local noon. By decreasing the distance from the observation
! ! to the dynamic f10_7 location we are allowing the already close
! ! observations to have a larger impact in the parameter estimation.
! ! 0.25 is heuristic. The 'close' observations have already been
! ! determined by the cutoff. Changing the distance here does not
! ! allow more observations to impact anything.
!    do k = 1, num_close
!       q_ind  = close_ind(k)
!       if  (loc_qtys(q_ind) == QTY_1D_PARAMETER) then
!          dist(k) = dist(k)*0.25_r8
!       endif
!    enddo
! endif
! 
! 
end subroutine get_close_state

!==================================================================

subroutine convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer,             intent(in)    :: loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus(:)

integer  :: current_vert_type, i
real(r8) :: height(1)
integer  :: local_status(1)

character(len=*), parameter :: routine = 'convert_vertical_obs'

! if ( which_vert == VERTISHEIGHT .or. which_vert == VERTISUNDEF) then
!   istatus(:) = 0
!   return
! endif
! 
! do i = 1, num
!    current_vert_type = nint(query_location(locs(i)))
!    if (( current_vert_type == which_vert ) .or. &
!        ( current_vert_type == VERTISUNDEF)) then
!       istatus(i) = 0
!       cycle
!    endif
! 
!   call model_interpolate(state_handle, 1, locs(i), QTY_GEOMETRIC_HEIGHT, height, local_status )
!   
!   if (local_status(1) == 0) call set_vertical(locs(i), height(1), VERTISHEIGHT)
!   istatus(i) = local_status(1)
! 
! enddo
! 
end subroutine convert_vertical_obs

!==================================================================
 subroutine convert_vertical_state(state_handle, num, locs, loc_qtys, loc_indx, &
                                   which_vert, istatus)
 
 type(ensemble_type), intent(in)    :: state_handle
 integer,             intent(in)    :: num
 type(location_type), intent(inout) :: locs(:)
 integer,             intent(in)    :: loc_qtys(:)
 integer(i8),         intent(in)    :: loc_indx(:)
 integer,             intent(in)    :: which_vert
 integer,             intent(out)   :: istatus
 
 integer :: var_id, dom_id, lon_index, lat_index, lev_index
 integer :: i
 real(r8) :: height(1), height1(1), height2(1)
 character(len=NF90_MAX_NAME) :: dim_name
 integer(i8) :: height_idx
 
 
! if  ( which_vert /= VERTISHEIGHT ) then
!    call error_handler(E_ERR,'convert_vertical_state', 'only supports VERTISHEIGHT')
! endif
! 
! istatus = 0 !HK what are you doing with this?
! 
! do i = 1, num
! 
!    call get_model_variable_indices(loc_indx(i), lon_index, lat_index, lev_index, var_id=var_id, dom_id=dom_id)
!    
!    ! search for either ilev or lev
!    dim_name = ilev_or_lev(dom_id, var_id)
!    
!    select case (trim(dim_name))
!       case ('ilev')
!          height_idx = get_dart_vector_index(lon_index, lat_index, lev_index, &
!                                             domain_id(SECONDARY_DOM), ivarZG)
!          height = get_state(height_idx, state_handle)/100.0_r8
!    
!       case (ALT_DIM_NAME) ! height on midpoint
!         height_idx = get_dart_vector_index(lon_index, lat_index, lev_index, &
!                                    domain_id(SECONDARY_DOM), ivarZG)
!         height1 = get_state(height_idx, state_handle)/100.0_r8
!         height_idx = get_dart_vector_index(lon_index, lat_index, lev_index+1, &
!                            domain_id(SECONDARY_DOM), ivarZG)
!         height2 = get_state(height_idx, state_handle)/100.0_r8
!         height = (height1 + height2) / 2.0_r8
!    
!       case default
!        call error_handler(E_ERR, 'convert_vertical_state', 'expecting ilev or ilat dimension')
!    end select
!    
!    locs(i) = set_location(lons(lon_index), lats(lat_index), height(1), VERTISHEIGHT)
! 
! end do
! 
end subroutine convert_vertical_state

!==================================================================

function read_model_time(filename)
type(time_type)              :: read_model_time
character(len=*), intent(in) :: filename

integer  :: ncid, i, ios
integer  :: tsimulation   ! the time read from a restart file; seconds from aeth_ref_date.
integer  :: ndays,nsecs

character(len=*), parameter :: routine = 'read_model_time'

tsimulation = MISSING_I

ncid = open_block_file(filename, 'read')
call nc_get_variable(ncid, 'time', tsimulation, routine)
call nc_close_file(ncid, routine, filename)

! Calculate the DART time of the file time.
! TODO: review calculation of ndays in read_model_time
ndays     = tsimulation/86400
nsecs     = tsimulation - ndays*86400 
! Need to subtract 1 because the ref day is not finished.
ndays     = aeth_ref_ndays -1 + ndays
read_model_time = set_time(nsecs,ndays)

if (do_output()) &
    call print_time(read_model_time,'read_model_time: time in restart file '//filename)
if (do_output()) &
    call print_date(read_model_time,'read_model_time: date in restart file '//filename)

if (debug > 8) then
   write(string1,*)'tsimulation ',tsimulation
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'ndays       ',ndays
   call error_handler(E_MSG,routine,string1,source,revision,revdate)
   write(string1,*)'nsecs       ',nsecs
   call error_handler(E_MSG,routine,string1,source,revision,revdate)

   call print_date(     aeth_ref_time, 'read_model_time:model base date')
   call print_time(     aeth_ref_time, 'read_model_time:model base time')
endif

end function read_model_time


!===============================================================================
! Routines below here are private to the module
!===============================================================================

! Fill up the variable_table from the namelist item 'variables'
! The namelist item variables is where a user specifies
! which variables they want in the DART state:
! variable name, dart qty, clamping min, clamping max, origin file, update or not

subroutine make_variable_table()

integer :: nfields_constructed   ! number of constructed state variables

integer  :: i, nrows, ncols

character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
character(len=NF90_MAX_NAME) :: minvalstring
character(len=NF90_MAX_NAME) :: maxvalstring
character(len=NF90_MAX_NAME) :: filename
character(len=NF90_MAX_NAME) :: state_or_aux

nrows = size(variable_table,1) ! these are MAX_NUM_VARIABLES, MAX_NUM_COLUMNS
ncols = size(variable_table,2)

! Convert the (input) 1D array "variables" into a table with six columns.
! The number of rows in the table correspond to the number of variables in the
! DART state vector.
! Column 1 is the netCDF variable name.
! Column 2 is the corresponding DART kind.
! Column 3 is the minimum value ("NA" if there is none) Not Applicable
! Column 4 is the maximum value ("NA" if there is none) Not Applicable
! Column 5 is the file of origin aether restart 'neutrals' or 'ions'
! Column 6 is whether or not the variable should be updated in the restart file.

nfields = 0
! TODO: TIEGCM uses 3 domains.  Aether may need only 1:
!       Do we need the 3rd category for derived fields; TEC, ...?
nfields_neutral = 0
nfields_ion = 0
nfields_constructed = 0

ROWLOOP : do i = 1, nrows

   varname      = trim(variables(ncols*i - 5))
   dartstr      = trim(variables(ncols*i - 4))
   minvalstring = trim(variables(ncols*i - 3))
   maxvalstring = trim(variables(ncols*i - 2))
   filename     = trim(variables(ncols*i - 1))
   state_or_aux = trim(variables(ncols*i    ))

! TODO: should Aether use the 6th column of namelist variable input to handle TEC, ...?
   call to_upper(state_or_aux) ! update or not

   variable_table(i,VT_VARNAMEINDX) = trim(varname)
   variable_table(i,VT_KINDINDX)    = trim(dartstr)
   variable_table(i,VT_MINVALINDX)  = trim(minvalstring)
   variable_table(i,VT_MAXVALINDX)  = trim(maxvalstring)
   variable_table(i,VT_ORIGININDX)  = trim(filename)
   variable_table(i,VT_STATEINDX)   = trim(state_or_aux)

   ! If the first element is empty, we have found the end of the list.
   if ((variable_table(i,1) == ' ') ) exit ROWLOOP

   ! Any other condition is an error.
   if ( any(variable_table(i,:) == ' ') ) then
      string1 = 'input.nml &model_nml:variables not fully specified.'
      string2 = 'Must be 6 entries per variable, last known variable name is'
      string3 = trim(variable_table(i,1))
      call error_handler(E_ERR,'get_variables_in_domain',string1, &
          source,revision,revdate,text2=string2,text3=string3)
   endif
! TODO; Modify this gitm error check for this routine?
!  ! Make sure DART kind is valid
!
!   if( get_index_for_quantity(dartstr) < 0 ) then
!      write(string1,'(3A)') 'there is no obs_kind "', trim(dartstr), '" in obs_kind_mod.f90'
!      call error_handler(E_ERR,routine,string1,source,revision,revdate)
!   endif

   nfields=nfields+1
   if (trim(variable_table(i,VT_ORIGININDX)) == 'neutrals')  then
      nfields_neutral     = nfields_neutral+1
   else if (trim(variable_table(i,VT_ORIGININDX)) == 'ions')      then
      nfields_ion         = nfields_ion+1
   else if (trim(variable_table(i,VT_ORIGININDX)) == 'CALCULATE')  then
      nfields_constructed = nfields_constructed + 1
   else
      print*,'variable_table(',i, VT_ORIGININDX,') = ', trim(variable_table(i,VT_ORIGININDX))
   endif
   print*,'make_variable_table: nfields = ',nfields, nfields_neutral, nfields_ion

enddo ROWLOOP

! Record the contents of the DART state vector
if (do_output() .and. (debug > 99)) then
   do i = 1,nfields
      write(*,'(''variable'',i4,'' is '',a12,1x,a32,4(1x,a20))') i, &
             trim(variable_table(i,1)), &
             trim(variable_table(i,2)), &
             trim(variable_table(i,3)), &
             trim(variable_table(i,4)), &
             trim(variable_table(i,5)), &
             trim(variable_table(i,6))
      write(logfileunit,'(''variable'',i4,'' is '',a12,1x,a32,4(1x,a20))') i, &
             trim(variable_table(i,1)), &
             trim(variable_table(i,2)), &
             trim(variable_table(i,3)), &
             trim(variable_table(i,4)), &
             trim(variable_table(i,5)), &
             trim(variable_table(i,6))
   enddo
endif

! TODO: Aether may need something like this.
! if (estimate_f10_7) then
!    if (nfields_constructed == 0) then
!       call error_handler(E_ERR, 'expecting f10.7 in &model_nml::variables', source)
!    endif
!    call load_up_state_structure_from_file(f10_7_file_name, nfields_constructed, 'CALCULATE', CONSTRUCT_DOM)
!    model_size = get_domain_size(RESTART_DOM) + get_domain_size(SECONDARY_DOM) &
!                           + get_domain_size(CONSTRUCT_DOM)
! else
!    model_size = get_domain_size(RESTART_DOM) + get_domain_size(SECONDARY_DOM)
! endif
! 
end subroutine make_variable_table

!==================================================================
! 
! ! Adds a domain to the state structure from a netcdf file
! ! Called from make_variable_table
! subroutine load_up_state_structure_from_file(filename, nvar, domain_name, domain_num)
! 
! character(len=*), intent(in) :: filename ! filename to read from
! integer,          intent(in) :: nvar ! number of variables in domain
! character(len=*), intent(in) :: domain_name ! restart, secondary
! integer,          intent(in) :: domain_num
! 
! integer :: i,j
! 
! character(len=NF90_MAX_NAME), allocatable :: var_names(:)
! real(r8), allocatable :: clamp_vals(:,:)
! integer, allocatable :: kind_list(:)
! logical, allocatable :: update_list(:)
! 
! 
! allocate(var_names(nvar), kind_list(nvar), &
!      clamp_vals(nvar,2), update_list(nvar))
! 
! update_list(:) = .true. ! default to update state variable
! clamp_vals(:,:) = MISSING_R8 ! default to no clamping
! 
! j = 0
! do i = 1, nfields
!    if (variable_table(i,VT_ORIGININDX) == trim(domain_name)) then
!       j = j+1
!       var_names(j) = variable_table(i, VT_VARNAMEINDX)
!       kind_list(j) = get_index_for_quantity(variable_table(i, VT_KINDINDX))
!       if (variable_table(i, VT_MINVALINDX) /= 'NA') then
!          read(variable_table(i, VT_MINVALINDX), '(d16.8)') clamp_vals(j,1)
!       endif
!       if (variable_table(i, VT_MAXVALINDX) /= 'NA') then
!         read(variable_table(i, VT_MAXVALINDX), '(d16.8)') clamp_vals(j,2)
!       endif
!       if (variable_table(i, VT_STATEINDX) == 'NO_COPY_BACK') then
!          update_list(j) = .false.
!       endif
!    endif
! enddo
! 
! domain_id(domain_num) = add_domain(filename, nvar, &
!                           var_names, kind_list, clamp_vals, update_list)
! 
! ! remove top level from all lev variables - this is the boundary condition
! call hyperslice_domain(domain_id(domain_num), ALT_DIM_NAME, nalt)
! 
! deallocate(var_names, kind_list, clamp_vals, update_list)
! 
! end subroutine load_up_state_structure_from_file
! 
!==================================================================
! 
! subroutine extrapolate_vtec(state_handle, ens_size, lon_index, lat_index, vTEC)
! !
! ! Create the vTEC from constituents in state.
! !
! 
! type(ensemble_type), intent(in)  :: state_handle
! integer,             intent(in)  :: ens_size
! integer,             intent(in)  :: lon_index, lat_index
! real(r8),            intent(out) :: vTEC(ens_size)
! 
! ! n(i)levs x ensmeble size
! real(r8), allocatable, dimension(:,:) :: NE
! real(r8), allocatable, dimension(:,:) :: TI, TE
! real(r8), allocatable, dimension(:,:) :: NEm_extended
! real(r8), allocatable, dimension(:,:)     :: NE_middle
! real(r8), dimension(ens_size)   :: GRAVITYtop, Tplasma, Hplasma
! 
! real(r8), PARAMETER :: k_constant = 1.381e-23_r8 ! m^2 * kg / s^2 / K
! real(r8), PARAMETER :: omass      = 2.678e-26_r8 ! mass of atomic oxgen kg
! 
! real(r8) :: earth_radiusm
! integer  :: naltX, nilevX, j, i, var_id
! integer(i8) :: idx
! 
! ! NE are extrapolated
! !  20 more layers for 2.5 degree resolution
! !  10 more layers for 5 degree resolution
! if (model_res == 2.5) then
!   naltX  = nalt + 20
!   nilevX  = nilev + 20
! else
!   naltX = nalt + 10
!   nilevX = nilev + 10
! endif
! 
! 
! allocate( NE(nilev, ens_size), NEm_extended(nilevX, ens_size))
! allocate( TI(nalt, ens_size), TE(nalt, ens_size) )
! allocate( NE_middle(naltX-1, ens_size) )
! 
! ! NE (interfaces)
! var_id = get_varid_from_varname(domain_id(RESTART_DOM), 'NE')
! do i = 1, nilev
!    idx = get_dart_vector_index(lon_index,lat_index, i, &
!                             domain_id(RESTART_DOM), var_id)
!    NE(i, :) = get_state(idx, state_handle)
! enddo
! 
! ! TI (midpoints)
! var_id = get_varid_from_varname(domain_id(RESTART_DOM), 'TI')
! do i = 1, nalt
!    idx = get_dart_vector_index(lon_index,lat_index, i, &
!                           domain_id(RESTART_DOM), var_id)
!    TI(i, :) = get_state(idx, state_handle)
! enddo
! 
! ! TE (midpoints)
! var_id = get_varid_from_varname(domain_id(RESTART_DOM), 'TE')
! do i = 1, nalt
!    idx = get_dart_vector_index(lon_index,lat_index, i, &
!                           domain_id(RESTART_DOM), var_id)
!    TE(i, :) = get_state(idx, state_handle)
! enddo
! 
! ! Construct vTEC given the parts
! 
! earth_radiusm = earth_radius * 1000.0_r8 ! Convert earth_radius in km to m
! NE            = NE * 1.0e+6_r8           ! Convert NE in #/cm^3 to #/m^3
! 
! ! Gravity at the top layer
! ! GRAVITYtop(:) = gravity * (earth_radiusm / (earth_radiusm + ZG(nilev,:))) ** 2
! 
! ! Plasma Temperature
! Tplasma(:) = (TI(nalt-1,:) + TE(nalt-1,:)) / 2.0_r8
! 
! ! Compute plasma scale height
! Hplasma(:) = (2.0_r8 * k_constant / omass ) * Tplasma(:) / GRAVITYtop(:)
! 
! NEm_extended(1:nilev,:) = NE
! 
! do j = nalt, naltX
!    NEm_extended(j,:) = NEm_extended(j-1,:) * exp(-0.5_r8)
! enddo
! 
! NE_middle(1:(naltX-1),:) = (NEm_extended(2:naltX,:) + NEm_extended(1:(naltX-1),:)) / 2.0_r8
! 
! do i = 1, ens_size
!    ! vTEC(i) = sum(NE_middle(:,i) * delta_ZG(:,i)) * 1.0e-16_r8 ! Convert to TECU (1.0e+16 #/m^2)
!    vTEC(i) = sum(NE_middle(:,i) ) * 1.0e-16_r8 ! Convert to TECU (1.0e+16 #/m^2)
! enddo
! 
! deallocate( NE, NEm_extended)
! deallocate( TI, TE )
! deallocate( NE_middle )
! 
! end subroutine extrapolate_vtec
! 
! !==================================================================
! 
! subroutine vert_interp(state_handle, n, dom_id, var_id, lon_index, lat_index, height, iqty, &
!                        val, istatus)
! ! returns the value at an arbitrary height on an existing horizontal grid location.
! ! istatus == 0 is success.
! 
! type(ensemble_type), intent(in) :: state_handle
! integer,          intent(in)  :: n ! ensemble_size
! integer,          intent(in)  :: dom_id
! integer,          intent(in)  :: var_id
! integer,          intent(in)  :: lon_index
! integer,          intent(in)  :: lat_index
! real(r8),         intent(in)  :: height
! integer,          intent(in)  :: iqty
! real(r8),         intent(out) :: val(n)
! integer,          intent(out) :: istatus(n)
! 
! logical :: is_pressure
! character(len=NF90_MAX_NAME) :: vertstagger
! 
! ! Presume the worst. Failure.
! istatus    = 1
! val        = MISSING_R8
! 
! is_pressure = (iqty == QTY_PRESSURE)
! if (is_pressure) then
!    vertstagger = 'ilev'
! else
!    vertstagger = ilev_or_lev(dom_id, var_id)
! endif
! 
! if (vertstagger == 'ilev') then
!   call vert_interp_ilev(state_handle, height, n, lon_index, lat_index, is_pressure, &
!                           dom_id, var_id, val, istatus)
! elseif (vertstagger == ALT_DIM_NAME) then
!   call vert_interp_lev(state_handle, height, n, lon_index, lat_index, is_pressure, &
!                           dom_id, var_id, val, istatus)
! endif
! 
! end subroutine vert_interp
! 
!==================================================================
! 
! subroutine find_qty_in_state(iqty, which_dom, var_id)
! ! Returns the variable id for a given DART qty
! ! Will return X rather than X_MN variable.
! 
! integer, intent(in)  :: iqty
! integer, intent(out) :: which_dom
! integer, intent(out) :: var_id
! 
! integer :: num_same_kind, id, k
! integer, allocatable :: multiple_kinds(:), n
! character(NF90_MAX_NAME) :: varname
! 
! which_dom = -1
! var_id = -1
! 
! do id = 1, get_num_domains() ! RESTART_DOM, SECONDARY_DOM, CONSTRUCT_DOM
! 
!    num_same_kind = get_num_varids_from_kind(domain_id(id), iqty)
!    if (num_same_kind == 0 ) cycle
!    if (num_same_kind  > 1 ) then ! need to pick which one you want
!      which_dom = id
!      allocate(multiple_kinds(num_same_kind))
!      call get_varids_from_kind(domain_id(id), iqty, multiple_kinds)
!      do k = 1, num_same_kind
!        varname = adjustl(get_variable_name(domain_id(id), multiple_kinds(k)))
!        n = len(trim(varname))
!        if (n <= 2) then ! variable name can not be X_MN
!           var_id = multiple_kinds(k)
!           exit
!        elseif (trim(varname(n-2:n)) == '_NM') then ! variable name is _MN
!           cycle ! assuming we want the X, not the X_MN
!        else
!          var_id = multiple_kinds(k)
!          exit
!        endif
!      enddo
!      deallocate(multiple_kinds)
!    else !
!       which_dom = id
!       var_id = get_varid_from_kind(domain_id(id), iqty)
!    endif
! enddo
! 
! end subroutine find_qty_in_state
! 
!==================================================================

! find enclosing lon indices
! Compute bracketing lon indices:
! TIEGCM [-180 175]  DART [180, 185, ..., 355, 0, 5, ..., 175]
subroutine compute_bracketing_lon_indices(lon, idx_below, idx_above, fraction)

real(r8), intent(in)  :: lon ! longitude
integer,  intent(out) :: idx_below, idx_above ! index in lons()
real(r8), intent(out) :: fraction ! fraction to use for interpolation

if(lon >= top_lon .and. lon < bot_lon) then     ! at wraparound point [175 <= lon < 180]
   idx_below = nlon
   idx_above = 1
   fraction = (lon - top_lon) / delta_lon
elseif (lon >= bot_lon) then                  ! [180 <= lon <= 360]
   idx_below = int((lon - bot_lon) / delta_lon) + 1
   idx_above = idx_below + 1
   fraction = (lon - lons(idx_below)) / delta_lon
else                                           ! [0 <= lon <= 175 ]
   idx_below = int((lon - 0.0_r8) / delta_lon) + zero_lon_index
   idx_above = idx_below + 1
   fraction = (lon - lons(idx_below)) / delta_lon
endif


end subroutine compute_bracketing_lon_indices

!==================================================================
! 
! ! on lev
! subroutine vert_interp_lev(state_handle, height, n, lon_index, lat_index, is_pressure, &
!                                      dom_id, var_id, val, istatus)
! 
! type(ensemble_type), intent(in) :: state_handle
! real(r8), intent(in)  :: height
! integer,  intent(in)  :: n ! ensemble size
! integer,  intent(in)  :: lon_index
! integer,  intent(in)  :: lat_index
! logical,  intent(in)  :: is_pressure
! integer,  intent(in)  :: dom_id, var_id
! real(r8), intent(out) :: val(n)  ! interpolated value
! integer,  intent(out) :: istatus(n)
! 
! integer :: lev(n), lev_minus_one(n), lev_plus_one(n)
! real(r8) :: frac_lev(n)
! 
! integer  :: k, i
! real(r8) :: delta_z(n)
! real(r8) :: zgrid_upper(n), zgrid_lower(n) ! ZG on midpoints
! real(r8) :: z_k(n), z_k_minus_one(n), z_k_plus_one(n)  ! ZG on ilves
! integer(i8)  :: indx_top(n), indx_bottom(n) ! state vector indices for qty
! integer(i8)  :: indx(n), indx_minus_one(n), indx_plus_one(n) ! state vector indices for ZG
! logical  :: found(n) ! track which ensemble members have been located
! real(r8) :: val_top(n), val_bottom(n)
! 
! istatus    = 1
! found = .false.
! 
!    ! Variable is on level midpoints, not ilevels.
!    ! Get height as the average of the ilevels.
! 
!    ! ilev index    1      2      3      4    ...  27    28    29
!    ! ilev value  -7.00, -6.50, -6.00, -5.50, ... 6.00, 6.50, 7.00 ;
!    !  lev value     -6.75, -6.25, -5.75, -5.25, ... 6.25, 6.75
!    !  lev index        1      2      3      4    ...  27    28
! 
!    !mid_level 1
!    zgrid_lower(:) = ( (get_state(get_dart_vector_index(lon_index,lat_index,1, &
!                           domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8)     +  &
!                     (get_state(get_dart_vector_index(lon_index,lat_index,2, &
!                           domain_id(SECONDARY_DOM), ivarZG), state_handle) /100.0_r8)  ) / 2.0_r8
! 
!    !mid_level nalt
!    zgrid_upper(:) = ( (get_state(get_dart_vector_index(lon_index,lat_index,nilev-1, &
!                              domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8)     +  &
!                      (get_state(get_dart_vector_index(lon_index,lat_index,nilev, &
!                           domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8) ) / 2.0_r8
! 
!    ! cannot extrapolate below bottom or beyond top
!    do i = 1, n
!       if ((zgrid_lower(i) > height) .or. (zgrid_upper(i) < height)) then
!         istatus(i) = 55
!       endif
!    enddo
!    if (any(istatus == 55)) return ! ! fail if any ensemble member fails
! 
!    ! Figure out what level is above/below, and by how much
!    h_loop_midpoint: do k = 2, nilev-1
! 
!     zgrid_upper(:) = ( (get_state(get_dart_vector_index(lon_index,lat_index,k, &
!                           domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8 )     +  &
!                        (get_state(get_dart_vector_index(lon_index,lat_index,k+1, &
!                           domain_id(SECONDARY_DOM), ivarZG), state_handle)/100.0_r8) ) / 2.0_r8
! 
!       ! per ensemble member
!       do i = 1, n
!          if (found(i)) cycle
!          if (height <= zgrid_upper(i)) then
!             found(i) = .true.
!             lev(i) = k
!             lev_minus_one(i) = lev(i) - 1
!             lev_plus_one(i) = lev(i) + 1
!             if (all(found)) exit h_loop_midpoint
!          endif
!       enddo
! 
!    enddo h_loop_midpoint
! 
!    do i = 1, n
!      indx(i) = get_dart_vector_index(lon_index,lat_index,lev(i), domain_id(SECONDARY_DOM), ivarZG)
!      indx_minus_one(i) = get_dart_vector_index(lon_index,lat_index,lev_minus_one(i), domain_id(SECONDARY_DOM), ivarZG)
!      indx_plus_one(i) = get_dart_vector_index(lon_index,lat_index,lev_plus_one(i), domain_id(SECONDARY_DOM), ivarZG)
!    enddo
! 
!    call get_state_array(z_k(:),indx(:), state_handle)
!    call get_state_array(z_k_minus_one, indx_minus_one(:), state_handle)  
!    call get_state_array(z_k_plus_one, indx_plus_one(:), state_handle)
! 
! 
!    !lower midpoint   
!    zgrid_lower(:) = ( z_k(:) + z_k_minus_one ) / 2.0_r8 / 100.0_r8
!    
!    ! upper midpoint
!    zgrid_upper(:) = ( z_k(:) + z_k_plus_one ) / 2.0_r8 / 100.0_r8
! 
!    where (zgrid_upper == zgrid_lower)  ! avoid divide by zero
!       frac_lev = 0.0_r8
!       delta_z = 0.0_r8
!    elsewhere
!       delta_z = zgrid_upper - zgrid_lower
!       frac_lev = (zgrid_upper - height)/delta_z
!    endwhere
! 
! if (is_pressure) then ! get fom plevs 
! 
!    val_top(:)    = plevs(lev(:))     !pressure at midpoint [Pa]
!    val_bottom(:) = plevs(lev_minus_one(:))  !pressure at midpoint [Pa]
!    val(:)        = exp(frac_lev(:) * log(val_bottom(:)) + (1.0 - frac_lev(:)) * log(val_top(:)))
! 
! else ! get from state vector
! 
!    do i = 1, n
!      indx_top(i) = get_dart_vector_index(lon_index,lat_index,lev(i), dom_id, var_id)
!      indx_bottom(i) = get_dart_vector_index(lon_index,lat_index,lev_minus_one(i), dom_id, var_id)
!    enddo
! 
!    call get_state_array(val_top, indx_top(:), state_handle)
!    call get_state_array(val_bottom, indx_bottom(:), state_handle)
! 
!    val(:) = frac_lev(:) * val_bottom(:)  + (1.0 - frac_lev(:)) * val_top(:)
! 
! endif
! 
! istatus(:) = 0
! 
! end subroutine vert_interp_lev
! 
! !==================================================================
! 
! ! Compute neighboring lat rows: TIEGCM [-87.5, 87.5] DART [-90, 90]
! ! Poles >|87.5| set to |87.5| 
! subroutine compute_bracketing_lat_indices(lat, idx_below, idx_above, fraction)
! 
! real(r8), intent(in)  :: lat ! latitude
! integer,  intent(out) :: idx_below, idx_above ! index in lats()
! real(r8), intent(out) :: fraction ! fraction to use for interpolation
! 
! if(lat >= bot_lat .and. lat < top_lat) then ! -87.5 <= lat < 87.5
!    idx_below = int((lat - bot_lat) / delta_lat) + 1
!    idx_above = idx_below + 1
!    fraction = (lat - lats(idx_below) ) / delta_lat
! else if(lat < bot_lat) then ! South of bottom lat
!    idx_below = 1
!    idx_above = 1
!    fraction = 1.0_r8
! else                        ! On or North of top lat
!    idx_below = nlat
!    idx_above = nlat
!    fraction = 1.0_r8
! endif
! 
! end subroutine compute_bracketing_lat_indices
! 
! !-------------------------------------------------------------------------------
! function interpolate(n, lon_fract, lat_fract, val11, val12, val21, val22) result(obs_val)
! 
! integer,  intent(in) :: n ! number of ensemble members
! real(r8), intent(in) :: lon_fract, lat_fract
! real(r8), dimension(n), intent(in) :: val11, val12, val21, val22
! real(r8), dimension(n) :: obs_val
! 
! real(r8) :: a(n, 2)
! 
! a(:, 1) = lon_fract * val21(:) + (1.0_r8 - lon_fract) * val11(:)
! a(:, 2) = lon_fract * val22(:) + (1.0_r8 - lon_fract) * val12(:)
! 
! obs_val(:) = lat_fract * a(:,2) + (1.0_r8 - lat_fract) * a(:,1)
! 
! end function interpolate
! 
! !-------------------------------------------------------------------------------
! function ilev_or_lev(dom_id, var_id) result(dim_name)
! 
! integer, intent(in) :: dom_id
! integer, intent(in) :: var_id
! character(len=NF90_MAX_NAME) :: dim_name
! 
! integer :: d
! ! search for either ilev or lev
! dim_name = 'null'
! do d = 1, get_num_dims(dom_id, var_id)
!    dim_name = get_dim_name(dom_id, var_id, d)
!    if (dim_name == 'ilev' .or. dim_name == ALT_DIM_NAME) exit
! enddo
! 
! end function ilev_or_lev
!===============================================================================
! End of model_mod
!===============================================================================
end module model_mod

! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the interface between the CICE ocean model and DART.
! author: C Bitz June 2016

! Modules that are absolutely required for use are listed
use              types_mod, only : r4, r8, i4, i8, SECPERDAY, MISSING_R8, rad2deg, &
                                   PI, metadatalength

use      time_manager_mod, only : time_type, set_time, set_date, get_date, get_time, &
                                  print_time, print_date, set_calendar_type, &
                                  operator(*),  operator(+), operator(-), &
                                  operator(>),  operator(<), operator(/), &
                                  operator(/=), operator(<=)

! location_mod is at rma-cice/location/xxx/location_mod.f90
! with the specific subdir set in path_names_xxx
! in our case, threed_sphere is the most appropriate choice
use          location_mod, only : location_type, get_dist, get_close_type,      &
                                  set_location, get_location,                   &
                                  get_close_obs,                                &
                                  loc_get_close_state => get_close_state,       &
                                  convert_vertical_obs, convert_vertical_state, &
                                  VERTISLEVEL  ! treat cat as vert level

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, nc_check, &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_end_define_mode

use location_io_mod,      only :  nc_write_location_atts, nc_get_location_varids, &
                                  nc_write_location

use default_model_mod,     only : init_time, init_conditions, adv_1step, &
                                  nc_write_model_vars

use         utilities_mod, only : register_module, error_handler,               &
                                  E_ERR, E_MSG, nmlfileunit, get_unit,  &
                                  do_output, to_upper, logfileunit,   &
                                  find_namelist_in_file, check_namelist_read,   &
                                  file_exist, find_textfile_dims, file_to_text, &
                                  do_nml_file, do_nml_term

! to add more kinds, edit ../../obs_kind/DEFAULT_obs_kind_mod.F90
use          obs_kind_mod, only : QTY_SEAICE_AGREG_CONCENTR  , &  
                                  QTY_SEAICE_AGREG_VOLUME    , &
                                  QTY_SEAICE_AGREG_SNOWVOLUME, &
                                  QTY_SEAICE_AGREG_THICKNESS , &
                                  QTY_SEAICE_AGREG_SNOWDEPTH , &
                                  QTY_SEAICE_CATEGORY        , &
                                  QTY_U_SEAICE_COMPONENT     , &
                                  QTY_V_SEAICE_COMPONENT     , &
                                  QTY_SEAICE_ALBEDODIRVIZ    , &
                                  QTY_SEAICE_ALBEDODIRNIR    , &
                                  QTY_SEAICE_ALBEDOINDVIZ    , &
                                  QTY_SEAICE_ALBEDOINDNIR    , &
                                  QTY_SEAICE_CONCENTR        , &
                                  QTY_SEAICE_VOLUME          , &
                                  QTY_SEAICE_SNOWVOLUME      , &
                                  QTY_SEAICE_SURFACETEMP     , &
                                  QTY_SEAICE_FIRSTYEARAREA   , &
                                  QTY_SEAICE_ICEAGE          , &
                                  QTY_SEAICE_LEVELAREA       , &
                                  QTY_SEAICE_LEVELVOLUME     , &
                                  QTY_SEAICE_MELTPONDAREA    , &
                                  QTY_SEAICE_MELTPONDDEPTH   , &
                                  QTY_SEAICE_MELTPONDLID     , &
                                  QTY_SEAICE_MELTPONDSNOW    , &
                                  QTY_SEAICE_SALINITY001     , &
                                  QTY_SEAICE_SALINITY002     , &
                                  QTY_SEAICE_SALINITY003     , &
                                  QTY_SEAICE_SALINITY004     , &
                                  QTY_SEAICE_SALINITY005     , &
                                  QTY_SEAICE_SALINITY006     , &
                                  QTY_SEAICE_SALINITY007     , &
                                  QTY_SEAICE_SALINITY008     , &
                                  QTY_SEAICE_ICEENTHALPY001  , &
                                  QTY_SEAICE_ICEENTHALPY002  , &
                                  QTY_SEAICE_ICEENTHALPY003  , &
                                  QTY_SEAICE_ICEENTHALPY004  , &
                                  QTY_SEAICE_ICEENTHALPY005  , &
                                  QTY_SEAICE_ICEENTHALPY006  , &
                                  QTY_SEAICE_ICEENTHALPY007  , &
                                  QTY_SEAICE_ICEENTHALPY008  , &
                                  QTY_SEAICE_SNOWENTHALPY001 , &
                                  QTY_SEAICE_SNOWENTHALPY002 , &
                                  QTY_SEAICE_SNOWENTHALPY003 , &
                                  QTY_DRY_LAND               , &
                                  QTY_SOM_TEMPERATURE        , &
                                  QTY_SEAICE_FY              , &
                                  QTY_SEAICE_AGREG_FY        , &
                                  QTY_SEAICE_AGREG_SURFACETEMP,&
                                  get_index_for_quantity     , &
                                  get_name_for_quantity

use     mpi_utilities_mod, only : my_task_id, task_count

use        random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use         dart_cice_mod, only : set_model_time_step,               &
                                  get_horiz_grid_dims,               &
                                  get_ncat_dim,                      &
                                  read_horiz_grid, read_topography,  &
                                  get_cice_restart_filename,         &
                                  set_binary_file_conversion

use  ensemble_manager_mod, only : ensemble_type

use distributed_state_mod, only : get_state

use   state_structure_mod, only : add_domain, get_model_variable_indices, &
                                  get_num_variables, get_index_start, &
                                  get_num_dims, get_domain_size, state_structure_info

use typesizes
use netcdf 

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
! these routines have code in this module
public :: get_model_size,                &
          adv_1step,                     &
          get_state_meta_data,           &
          model_interpolate,             &
          shortest_time_between_assimilations, &
          static_init_model,             &
          read_model_time,               &
          nc_write_model_atts,           &
          write_model_time,              &
          get_close_state,               &
          end_model

! the code for these routines are in other modules
public :: init_time,                     &
          init_conditions,               &
          nc_write_model_vars,           &
          pert_model_copies,             &
          get_close_obs,                 &
          convert_vertical_obs,          &
          convert_vertical_state

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: get_cice_restart_filename, test_interpolation

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! message strings
character(len=512) :: string1
character(len=512) :: string2
character(len=512) :: string3

logical, save :: module_initialized = .false.

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

! DART state vector contents are specified in the input.nml:&model_nml namelist.
integer, parameter :: max_state_variables = 10 
integer, parameter :: num_state_table_columns = 3
character(len=NF90_MAX_NAME) :: variable_table( max_state_variables, num_state_table_columns )
integer :: state_kinds_list( max_state_variables )
logical :: update_var_list( max_state_variables )

! identifiers for variable_table
integer, parameter :: VAR_NAME_INDEX = 1
integer, parameter :: VAR_QTY_INDEX = 2
integer, parameter :: VAR_UPDATE_INDEX = 3

! things which can/should be in the model_nml
integer  :: assimilation_period_days = 1
integer  :: assimilation_period_seconds = 0
real(r8) :: model_perturbation_amplitude = 0.2
logical  :: update_dry_cell_walls = .false.
character(len=metadatalength) :: model_state_variables(max_state_variables * num_state_table_columns ) = ' '
integer  :: debug = 0   ! turn up for more and more debug messages

! valid values:  native, big_endian, little_endian
character(len=64) :: binary_grid_file_format = 'big_endian'


! FIXME: currently the update_dry_cell_walls namelist value DOES
! NOTHING.  it needs additional code to detect the cells which are
! wet, but within 1 cell of the bottom/sides/etc.  

namelist /model_nml/  &
   assimilation_period_days,     &  ! for now, this is the timestep
   assimilation_period_seconds,  &
   model_perturbation_amplitude, &
   update_dry_cell_walls,        &
   model_state_variables,        &
   binary_grid_file_format,      &
   debug

!------------------------------------------------------------------
!
! The DART state vector (control vector) will consist of (from cice restart):  
! uvel(i,j), vvel(i,j)     ice velocity
! aicen(i,j,n)             area of ice (unitless, also called as conc. or frac.)
! vicen(i,j,n)             volume of ice (volume per unit area)
! vsnon(i,j,n)             volume of snow (volume per unit area)
! Tsfcn(i,j,n)             ice/snow surface temperature
! iage(i,j,n)              ice age (purely diagnostic)
! FY(i,j,n)                firstyear ice fraction of each cat area (purely diag.)
! alvl(i,j,n)    area of level ice (purely diagnostic)
! vlvl(i,j,n)    volume of level ice (purely diagnostic)
! swvdr(i,j)
! swvdf(i,j)
! swidr(i,j)
! swidf(i,j)
! apnd(i,j,n)    melt pond area fraction
! hpnd(i,j,n)    melt pond depth
! ipnd(i,j,n)    melt pond refrozen lid thickness
! dhs(i,j,n)     depth difference for snow on sea ice and pond ice
! sice001(i,j,n)           sea ice salinity per volume in each cat, layers 001-008
! qice001(i,j,n)           sea ice enthalpy per volume in each cat, layers 001-008
! qsno001(i,j,n)           snow enthalpy per snovol, layers 001-003

! The 4 main thermo variable names have 'ice' or 'sno' or 'sfc'
! with an 'n' at the end to indicate they are per category
! All are located at cell centers, except uvel, vvel are at grid cell corners.
!
! consider the 3 dimensions lat, lon, and cat
! so starting from POP code had to change vert dimension to cat throughout 
! 
! Additional variables can be read into the state vector using the 
! model_state_variables namelist by specifying the netcdf variable name
! dart kind string and an update string.  Currently the update string
! is not being used.
!------------------------------------------------------------------

! Number of fields in the state vector
integer :: nfields

! Grid parameters - the values will be read from a
! standard cice namelist and filled in here.

! grid counts for each field
! We should also make sure we never
! hardwire the value '5' anywhere - always refer to Ncat in loops.
integer :: Nx=-1, Ny=-1   ! size of the dipole (or irregular) grids. 
integer :: Ncat=-1        ! number of categories in ice-thickness dist

! These arrays store the longitude and latitude of the lower left corner of
! each of the dipole u quadrilaterals and t quadrilaterals.
real(r8), allocatable :: ULAT(:,:), ULON(:,:), TLAT(:,:), TLON(:,:)

! integer, lowest valid cell number in the vertical
integer, allocatable  :: KMT(:, :), KMU(:, :)

! real, depth of lowest valid cell (0 = land).  use only if KMT/KMU not avail.
real(r8), allocatable :: HT(:,:), HU(:,:)

real(r8)        :: endTime
real(r8)        :: ice_thermo_timestep = 3600.0_r4   ! hardwired in cice 
integer         :: timestepcount = 0
type(time_type) :: model_time, model_timestep

integer(i8) :: model_size    ! the state vector length

! NOTE (dipole/tripole grids): since both of the dipole and tripole
! grids are logically rectangular we can use the same interpolation
! scheme originally implemented for the dipole grid. Here we can 
! interchange dipole and tripole when reading the code.

! The regular grid used for dipole interpolation divides the sphere into
! a set of regularly spaced lon-lat boxes. The number of boxes in
! longitude and latitude are set by num_reg_x and num_reg_y. Making the
! number of regular boxes smaller decreases the computation required for
! doing each interpolation but increases the static storage requirements
! and the initialization computation (which seems to be pretty small).
! FIX ME: to account for various grid sizes we should dynamically 
! allocate these numbers.  To keep max_reg_list_num < 100 we can use:
!    tx0.1v2 num_reg_x = num_reg_y = 900
!    tx0.5v1 num_reg_x = num_reg_y = 180
!    gx1v6   num_reg_x = num_reg_y = 90
! Larger num_reg_(x,y) values require more temporary storage in 
! ureg_list_lon, ureg_list_lat, treg_list_lon, treg_list_lat. For now
! we can use num_reg_(x,y) = 180 and max_reg_list_num = 800 to account
! for all of the currently implemented grid types.
integer, parameter :: num_reg_x = 180, num_reg_y = 180

! The max_reg_list_num controls the size of temporary storage used for
! initializing the regular grid. Four arrays
! of size num_reg_x*num_reg_y*max_reg_list_num are needed. The initialization
! fails and returns an error if max_reg_list_num is too small. With 180 regular
! lat lon boxes a value of 30 is sufficient for the gx3 POP grid, 80 for the 
! gx1 grid, 180 for the tx0.5 grid and 800 for the tx0.1 grid.
! FIX ME: we should declare this at runtime depending on the grid size.
integer, parameter :: max_reg_list_num = 800

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

! global domain id to be used by routines in state_structure_mod
integer :: domain_id


!------------------------------------------------------------------
contains
!------------------------------------------------------------------


subroutine static_init_model()

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
!   figure out model timestep (not really used in cice)
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
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
if (do_nml_file())  write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type('Gregorian')

! Set the time step ... causes cice namelists to be read.
! Ensures model_timestep is multiple of 'ice_thermo_timestep'

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

! BEFORE calling grid routines, set the endian-ness of the binary files if needed.
call set_binary_file_conversion(binary_grid_file_format)

! get data dimensions, then allocate space, then open the files
! and actually fill in the arrays.

call get_horiz_grid_dims(Nx, Ny)
call get_ncat_dim(Ncat)

! Allocate space for grid variables. 
allocate(ULAT(Nx,Ny), ULON(Nx,Ny), TLAT(Nx,Ny), TLON(Nx,Ny))
allocate( KMT(Nx,Ny),  KMU(Nx,Ny))
allocate(  HT(Nx,Ny),   HU(Nx,Ny))

! Fill them in.
! horiz grid initializes ULAT/LON, TLAT/LON as well.
! kmt initializes HT/HU if present in input file.

call read_horiz_grid(Nx, Ny, ULAT, ULON, TLAT, TLON)
call read_topography(Nx, Ny,  KMT,  KMU)

if (debug > 2) call write_grid_netcdf()     ! DEBUG only
if (debug > 2) call write_grid_interptest() ! DEBUG only

! verify that the model_state_variables namelist was filled in correctly.  
! returns variable_table which has variable names, kinds and update strings.
call verify_state_variables(model_state_variables, nfields, variable_table, &
                            state_kinds_list, update_var_list)

! in spite of the staggering, all grids are the same size
! and offset by half a grid cell.  4 are 3D and 2 are 2D.
!  e.g. aicen,vicen,vsnon = 256 x 225 x 5
!  e.g. uvel,vvel = 256 x 225 

if (do_output()) write(logfileunit, *) 'Using grid : Nx, Ny, Ncat = ', &
                                                     Nx, Ny, Ncat
if (do_output()) write(     *     , *) 'Using grid : Nx, Ny, Ncat = ', &
                                                     Nx, Ny, Ncat

! Initialize the interpolation routines
call init_interp()

! Determine the shape of the variables from "cice.r.nc"
! The assimilate.csh, perfect_model.csh must ensure the cice restart file
! is linked to this filename.
domain_id = add_domain('cice.r.nc', nfields, &
                       var_names = variable_table(1:nfields, VAR_NAME_INDEX), &
                       kind_list = state_kinds_list(1:nfields), &
                       update_list = update_var_list(1:nfields))

if (debug > 2) call state_structure_info(domain_id)

model_size = get_domain_size(domain_id)
if (do_output()) write(*,*) 'model_size = ', model_size

end subroutine static_init_model

!------------------------------------------------------------

subroutine init_interp()

! Initializes data structures needed for cice interpolation for
! either dipole or irregular grid.
! This should be called at static_init_model time to avoid 
! having all this temporary storage in the middle of a run.

integer :: i

! Determine whether this is a irregular lon-lat grid or a dipole.
! Do this by seeing if the lons have the same values at both
! the first and last latitude row; this is not the case for dipole. 

dipole_grid = .false.
do i = 1, nx
   if(ulon(i, 1) /= ulon(i, ny)) then
      dipole_grid = .true.
      call init_dipole_interp()
      return
   endif
enddo

end subroutine init_interp

!------------------------------------------------------------

subroutine init_dipole_interp()

! Build the data structure for interpolation for a dipole grid.

! Need a temporary data structure to build this.
! These arrays keep a list of the x and y indices of dipole quads 
! that potentially overlap the regular boxes. Need one for the u 
! and one for the t grid.
integer, allocatable :: ureg_list_lon(:,:,:)
integer, allocatable :: ureg_list_lat(:,:,:)
integer, allocatable :: treg_list_lon(:,:,:)
integer, allocatable :: treg_list_lat(:,:,:)

real(r8) :: u_c_lons(4), u_c_lats(4), t_c_lons(4), t_c_lats(4), pole_row_lon
integer  :: i, j, k, pindex
integer  :: reg_lon_ind(2), reg_lat_ind(2), u_total, t_total, u_index, t_index
logical  :: is_pole
integer  :: surf_index

allocate(ureg_list_lon(num_reg_x, num_reg_y, max_reg_list_num))
allocate(ureg_list_lat(num_reg_x, num_reg_y, max_reg_list_num))
allocate(treg_list_lon(num_reg_x, num_reg_y, max_reg_list_num))
allocate(treg_list_lat(num_reg_x, num_reg_y, max_reg_list_num))

! this is the level threshold for deciding whether we are over land
! or water.  to be valid all 4 corners of the quad must have a level
! number greater than this index.  (so 0 excludes all land points.)
! if you wanted to assimilate only in regions where the water depth is
! deeper than some threshold, set this index to N and only quads where
! all the level numbers are N+1 or deeper will be used.
surf_index = 1

! Begin by finding the quad that contains the pole for the dipole t_grid. 
! To do this locate the u quad with the pole on its right boundary. This is on
! the row that is opposite the shifted pole and exactly follows a lon circle.
pole_x = nx / 2;
! Search for the row at which the longitude flips over the pole
pole_row_lon = ulon(pole_x, 1);
do i = 1, ny
   pindex = i
   if(ulon(pole_x, i) /= pole_row_lon) exit
enddo

! Pole boxes for u have indices pole_x or pole_x-1 and index - 1;
! (it's right before the flip).
u_pole_y = pindex - 1;

! Locate the T dipole quad that contains the pole.
! We know it is in either the same lat quad as the u pole edge or one higher.
! Figure out if the pole is more or less than halfway along
! the u quad edge to find the right one.
if(ulat(pole_x, u_pole_y) > ulat(pole_x, u_pole_y + 1)) then
   t_pole_y = u_pole_y;
else
   t_pole_y = u_pole_y + 1;
endif

! Loop through each of the dipole grid quads 
do i = 1, nx
   ! There's no wraparound in y, one box less than grid boundaries
   do j = 1, ny - 1
      
      ! Only update regular boxes that contain all wet corners
      if( all_corners_wet(QTY_U_SEAICE_COMPONENT,i,j) ) then
         ! Set up array of lons and lats for the corners of these u quads
         call get_quad_corners(ulon, i, j, u_c_lons)
         call get_quad_corners(ulat, i, j, u_c_lats)

         ! Get list of regular boxes that cover this u dipole quad
         ! false indicates that for the u grid there's nothing special about pole
         call reg_box_overlap(u_c_lons, u_c_lats, .false., reg_lon_ind, reg_lat_ind)         
         ! Update the temporary data structures for the u quad 
         call update_reg_list(u_dipole_num, ureg_list_lon, &
            ureg_list_lat, reg_lon_ind, reg_lat_ind, i, j)
      endif 

      ! Repeat for t dipole quads.
      ! Only update regular boxes that contain all wet corners
      if( all_corners_wet(QTY_SEAICE_AGREG_CONCENTR,i,j) ) then
         ! Set up array of lons and lats for the corners of these t quads
         call get_quad_corners(tlon, i, j, t_c_lons)
         call get_quad_corners(tlat, i, j, t_c_lats)

         ! Is this the pole quad for the T grid?
         is_pole = (i == pole_x .and. j == t_pole_y)
         
         call reg_box_overlap(t_c_lons, t_c_lats, is_pole, reg_lon_ind, reg_lat_ind)         
         call update_reg_list(t_dipole_num, treg_list_lon, &
            treg_list_lat, reg_lon_ind, reg_lat_ind, i, j)
      endif
   enddo
enddo

if (debug > 1 .and. do_output()) then
   write(*,*)'to determine (minimum) max_reg_list_num values for new grids ...'
   write(*,*)'u_dipole_num is ',maxval(u_dipole_num)
   write(*,*)'t_dipole_num is ',maxval(t_dipole_num)
endif

! Invert the temporary data structure. The total number of entries will be 
! the sum of the number of dipole cells for each regular cell. 
u_total = sum(u_dipole_num)
t_total = sum(t_dipole_num)

! Allocate storage for the final structures in module storage
allocate(u_dipole_lon_list(u_total), u_dipole_lat_list(u_total))
allocate(t_dipole_lon_list(t_total), t_dipole_lat_list(t_total))

! Fill up the long list by traversing the temporary structure. Need indices 
! to keep track of where to put the next entry.
u_index = 1
t_index = 1

! Loop through each regular grid box
do i = 1, num_reg_x
   do j = 1, num_reg_y

      ! The list for this regular box starts at the current indices.
      u_dipole_start(i, j) = u_index
      t_dipole_start(i, j) = t_index

      ! Copy all the close dipole quads for regular u box(i, j)
      do k = 1, u_dipole_num(i, j)
         u_dipole_lon_list(u_index) = ureg_list_lon(i, j, k) 
         u_dipole_lat_list(u_index) = ureg_list_lat(i, j, k) 
         u_index = u_index + 1
      enddo
      
      ! Copy all the close dipoles for regular t box (i, j)
      do k = 1, t_dipole_num(i, j)
         t_dipole_lon_list(t_index) = treg_list_lon(i, j, k) 
         t_dipole_lat_list(t_index) = treg_list_lat(i, j, k) 
         t_index = t_index + 1
      enddo

   enddo
enddo

! Confirm that the indices come out okay as debug
if(u_index /= u_total + 1) then
   string1 = 'Storage indices did not balance for U grid: : contact DART developers'
   call error_handler(E_ERR, 'init_dipole_interp', string1, source, revision, revdate)
endif
if(t_index /= t_total + 1) then
   string1 = 'Storage indices did not balance for T grid: : contact DART developers'
   call error_handler(E_ERR, 'init_dipole_interp', string1, source, revision, revdate)
endif

end subroutine init_dipole_interp

!------------------------------------------------------------

subroutine get_reg_box_indices(lon, lat, x_ind, y_ind)
 real(r8), intent(in)  :: lon, lat
 integer,  intent(out) :: x_ind, y_ind

! Given a longitude and latitude in degrees returns the index of the regular
! lon-lat box that contains the point.

call get_reg_lon_box(lon, x_ind)
call get_reg_lat_box(lat, y_ind)

end subroutine get_reg_box_indices

!------------------------------------------------------------

subroutine get_reg_lon_box(lon, x_ind)
 real(r8), intent(in)  :: lon
 integer,  intent(out) :: x_ind

! Determine which regular longitude box a longitude is in.

x_ind = int(num_reg_x * lon / 360.0_r8) + 1

! Watch out for exactly at top; assume all lats and lons in legal range
if(lon == 360.0_r8) x_ind = num_reg_x

end subroutine get_reg_lon_box

!------------------------------------------------------------

subroutine get_reg_lat_box(lat, y_ind)
 real(r8), intent(in)  :: lat
 integer,  intent(out) :: y_ind

! Determine which regular latitude box a latitude is in.

y_ind = int(num_reg_y * (lat + 90.0_r8) / 180.0_r8) + 1

! Watch out for exactly at top; assume all lats and lons in legal range
if(lat == 90.0_r8)  y_ind = num_reg_y

end subroutine get_reg_lat_box

!------------------------------------------------------------

subroutine reg_box_overlap(x_corners, y_corners, is_pole, reg_lon_ind, reg_lat_ind)
 real(r8), intent(in)  :: x_corners(4), y_corners(4)
 logical,  intent(in)  :: is_pole
 integer,  intent(out) :: reg_lon_ind(2), reg_lat_ind(2)

! Find a set of regular lat lon boxes that covers all of the area covered by 
! a dipole grid qaud whose corners are given by the dimension four x_corners 
! and y_corners arrays.  The two dimensional arrays reg_lon_ind and reg_lat_ind
! return the first and last indices of the regular boxes in latitude and
! longitude respectively. These indices may wraparound for reg_lon_ind.  
! A special computation is needed for a dipole quad that has the true north 
! pole in its interior. The logical is_pole is set to true if this is the case.
! This can only happen for the t grid.  If the longitude boxes overlap 0
! degrees, the indices returned are adjusted by adding the total number of
! boxes to the second index (e.g. the indices might be 88 and 93 for a case
! with 90 longitude boxes).

real(r8) :: lat_min, lat_max, lon_min, lon_max
integer  :: i

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
   ! This is specific to the dipole POP grids that do not go to the south pole
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

!------------------------------------------------------------

subroutine get_quad_corners(x, i, j, corners)
 real(r8), intent(in)  :: x(:, :)
 integer,  intent(in)  :: i, j
 real(r8), intent(out) :: corners(4)

! Grabs the corners for a given quadrilateral from the global array of lower
! right corners. Note that corners go counterclockwise around the quad.

integer :: ip1

! Have to worry about wrapping in longitude but not in latitude
ip1 = i + 1
if(ip1 > nx) ip1 = 1

corners(1) = x(i,   j  ) 
corners(2) = x(ip1, j  )
corners(3) = x(ip1, j+1)
corners(4) = x(i,   j+1)

end subroutine get_quad_corners

!------------------------------------------------------------

subroutine update_reg_list(reg_list_num, reg_list_lon, reg_list_lat, &
                           reg_lon_ind, reg_lat_ind, dipole_lon_index, dipole_lat_index)

 integer, intent(inout) :: reg_list_num(:, :), reg_list_lon(:, :, :), reg_list_lat(:, :, :)
 integer, intent(inout) :: reg_lon_ind(2), reg_lat_ind(2)
 integer, intent(in)    :: dipole_lon_index, dipole_lat_index
 
! Updates the data structure listing dipole quads that are in a given regular box
integer :: ind_x, index_x, ind_y

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

!------------------------------------------------------------------
!> Returns the number of items in the state vector

function get_model_size()
 integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size

!------------------------------------------------------------------

subroutine model_interpolate(state_handle, ens_size, location, obs_type, expected_obs, istatus)

 type(ensemble_type), intent(in) :: state_handle
 integer,             intent(in) :: ens_size
 type(location_type), intent(in) :: location
 integer,             intent(in) :: obs_type
 integer,            intent(out) :: istatus(ens_size)
 real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values

! Model interpolate will interpolate any state variable (all on tracer grid for now) to
! the given location given a state vector. The type of the variable being
! interpolated is obs_type since normally this is used to find the expected
! value of an observation at some location. The interpolated value is 
! returned in interp_val and istatus is 0 for success.

! Local storage
real(r8)       :: loc_array(3), llon, llat
integer(i8)    :: base_offset
integer        :: cat_index, cat_signal, icat, cat_signal_interm
real(r8)       :: expected_aggr_conc(ens_size) ! array of interpolated values

!Fei---need aicen*fyn to calculate the aggregate FY concentration------------
real(r8)       :: expected_conc(ens_size)
real(r8)       :: expected_fy(ens_size)
real(r8)       :: expected_tsfc(ens_size)
real(r8)       :: temp(ens_size)
real(r8)       :: temp1(ens_size)
!----------------------------------------------------------------------------

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the 
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

expected_obs(:) = MISSING_R8     ! the DART bad value flag
istatus(:) = 99                ! unknown error

! Get the individual locations values
loc_array = get_location(location)
llon    = loc_array(1)
llat    = loc_array(2)
cat_index = int(loc_array(3))

! Special case. Only used when trying to determine the number of ice categories.
! Note the early return. The actual lat/lon is unimportant.

if (obs_type == QTY_SEAICE_CATEGORY) then
   if (cat_index <= Ncat) then
      istatus      = 0
      expected_obs = cat_index
      RETURN
   endif
endif

if (debug > 1) then
   print *, 'requesting interpolation of ', obs_type, ' at ', llon, llat, cat_index
endif

! The base_offset is the index of state vector that corresponds to 
! a variable block without regard to level or location, so it can be 
! viewed as pointing to first cat for cice 3d vars

! set the base_offset. for aggregate quantities set it to the start of the first
! category of the appropriate kind. set the cat_signal that will be passed into 
! the lon_lat_interpolate() routine
! cat_signal says whether to aggregate multiple categories into a sum, interpolate a 3d var 
! on a particular level, or interpolate a simple 2d variable.
SELECT CASE (obs_type)
   CASE (QTY_SEAICE_AGREG_THICKNESS )  ! these kinds require aggregating 3D vars to make a 2D var
      cat_signal = -1 ! for extra special procedure to aggregate
      base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_VOLUME))  
   CASE (QTY_SEAICE_AGREG_SNOWDEPTH )  ! these kinds require aggregating 3D vars to make a 2D var
      cat_signal = -1 ! for extra special procedure to aggregate
      base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_SNOWVOLUME))  
   CASE (QTY_SEAICE_AGREG_CONCENTR )   ! these kinds require aggregating a 3D var to make a 2D var
      cat_signal = 0 ! for aggregate variable, send signal to lon_lat_interp
      base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_CONCENTR))  
   CASE (QTY_SEAICE_AGREG_VOLUME   )   ! these kinds require aggregating a 3D var to make a 2D var
      cat_signal = 0 ! for aggregate variable, send signal to lon_lat_interp
      base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_VOLUME))  
   CASE (QTY_SEAICE_AGREG_SNOWVOLUME ) ! these kinds require aggregating a 3D var to make a 2D var
      cat_signal = 0 ! for aggregate variable, send signal to lon_lat_interp
      base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_SNOWVOLUME))  
   CASE (QTY_SEAICE_AGREG_SURFACETEMP) ! FEI need aicen to average the temp, have not considered open water temp yet
      cat_signal = -3
      base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_SURFACETEMP))
   CASE (QTY_SOM_TEMPERATURE) ! these kinds are 1d variables
      cat_signal = 3
      base_offset = get_index_start(domain_id,get_varid_from_kind(QTY_SOM_TEMPERATURE))
   CASE (QTY_SEAICE_CONCENTR       , &  ! these kinds have an additional dim for category
         QTY_SEAICE_FY       , & 
         QTY_SEAICE_VOLUME         , &
         QTY_SEAICE_SNOWVOLUME     , &
         QTY_SEAICE_SURFACETEMP    , &
         QTY_SEAICE_FIRSTYEARAREA  , &
         QTY_SEAICE_ICEAGE         , &
         QTY_SEAICE_LEVELAREA      , &
         QTY_SEAICE_LEVELVOLUME    , &
         QTY_SEAICE_MELTPONDAREA   , &
         QTY_SEAICE_MELTPONDDEPTH  , &
         QTY_SEAICE_MELTPONDLID    , &
         QTY_SEAICE_MELTPONDSNOW   , &
         QTY_SEAICE_SALINITY001    , &
         QTY_SEAICE_SALINITY002    , &
         QTY_SEAICE_SALINITY003    , &
         QTY_SEAICE_SALINITY004    , &
         QTY_SEAICE_SALINITY005    , &
         QTY_SEAICE_SALINITY006    , &
         QTY_SEAICE_SALINITY007    , &
         QTY_SEAICE_SALINITY008    , &
         QTY_SEAICE_ICEENTHALPY001 , &
         QTY_SEAICE_ICEENTHALPY002 , &
         QTY_SEAICE_ICEENTHALPY003 , &
         QTY_SEAICE_ICEENTHALPY004 , &
         QTY_SEAICE_ICEENTHALPY005 , &
         QTY_SEAICE_ICEENTHALPY006 , &
         QTY_SEAICE_ICEENTHALPY007 , &
         QTY_SEAICE_ICEENTHALPY008 , &
         QTY_SEAICE_SNOWENTHALPY001, &
         QTY_SEAICE_SNOWENTHALPY002, &
         QTY_SEAICE_SNOWENTHALPY003  )
      ! move pointer to the particular category
      ! then treat as 2d field in lon_lat_interp
      base_offset = get_index_start(domain_id, get_varid_from_kind(obs_type))
      base_offset = base_offset + (cat_index-1) * Nx * Ny 
      cat_signal = 1 ! now same as boring 2d field
   CASE ( QTY_U_SEAICE_COMPONENT    , &   ! these kinds are just 2D vars
          QTY_V_SEAICE_COMPONENT    , &
          QTY_SEAICE_ALBEDODIRVIZ   , &
          QTY_SEAICE_ALBEDODIRNIR   , &
          QTY_SEAICE_ALBEDOINDVIZ   , &
          QTY_SEAICE_ALBEDOINDNIR   )
      base_offset = get_index_start(domain_id, get_varid_from_kind(obs_type))
      cat_signal = 2 ! also boring 2d field (treat same as cat_signal 1)
   CASE DEFAULT
      ! Not a legal type for interpolation, return istatus error
      istatus = 15
      return
END SELECT

! finally the interpolation for all sea ice kinds

!Fei--------need aicen and fyn-----------------------------------
!special treatment to calculate the aggregate fy = sum(aicen*fyn)/aice
if (cat_signal == -2) then
   temp = 0.0_r8
   temp1= 0.0_r8
   do icat = 1,Ncat
      !reads in aicen 
      cat_signal_interm = 1
      base_offset = get_index_start(domain_id,get_varid_from_kind(QTY_SEAICE_CONCENTR))
      base_offset = base_offset + (icat-1) * Nx * Ny
      call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, obs_type, cat_signal_interm, expected_conc, istatus)
      !reads in fyn
      cat_signal_interm = 1
      base_offset = get_index_start(domain_id,get_varid_from_kind(QTY_SEAICE_FY))
      base_offset = base_offset + (icat-1) * Nx * Ny
      call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, obs_type, cat_signal_interm, expected_fy, istatus)
   temp = temp + expected_conc * expected_fy  !sum(aicen*fyn) = FY % over ice
   temp1= temp1+ expected_conc                        !sum(aicen) = aice
  
      if (any(expected_conc<0.0) .or. any(expected_conc>1.0))then
      print*,'obstype FY expected sicn:',expected_conc
      print*,'FY sicn lat lon:',llat,llon
      endif
      if (any(expected_fy>1.0) .or. any(expected_fy<0.0)) then
      print*,'obstype FY expected fyn:',expected_fy,llat,llon
      print*,'FY fyn lat lon:',llat,llon
      endif

   end do 
   expected_obs = temp/max(temp1,1.0e-8)  !sum(aicen*fyn)/aice = FY % in the gridcell
else if (cat_signal == -3 ) then
   temp = 0.0_r8
   temp1= 0.0_r8
   do icat = 1,Ncat
      !reads in aicen 
      cat_signal_interm = 1
      base_offset = get_index_start(domain_id,get_varid_from_kind(QTY_SEAICE_CONCENTR))
      base_offset = base_offset + (icat-1) * Nx * Ny
      call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, obs_type, cat_signal_interm, expected_conc, istatus)
      !reads in Tsfcn
      cat_signal_interm = 1
      base_offset = get_index_start(domain_id,get_varid_from_kind(QTY_SEAICE_SURFACETEMP))
      base_offset = base_offset + (icat-1) * Nx * Ny
      call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, obs_type, cat_signal_interm, expected_tsfc, istatus)
      if (any(expected_conc<0.0) .or. any(expected_conc>1.0))then
      print*,'obstype TSFC expected sicn:',expected_conc
      print*,'TSFC sicn lat lon:',llat,llon
      endif
      if (any(expected_tsfc>50.0) .or. any(expected_tsfc<-100.0)) then
      print*,'obstype TSFC expected tsfcn:',expected_tsfc
      print*,'TSFC tsfcn lat lon:',llat,llon
      endif
      temp = temp + expected_conc * expected_tsfc  !sum(aicen*Tsfcn)
      temp1= temp1+ expected_conc                  !sum(aicen) = aice
   end do
   expected_obs = temp/max(temp1,1.0e-8)  !sum(aicen*Tsfcn)/aice = Tsfc ;averaged temperature over sea-ice covered portion
   if (any(expected_obs>50.0) .or. any(expected_obs<-100.0)) then
      print*,'obstype TSFC expected obs:',expected_obs
      print*,'TSFC tsfc lat lon:' ,llat,llon
      print*,'temp:',temp
      print*,'temp1:',temp1
   endif

else
    call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, obs_type, cat_signal, expected_obs, istatus)
      if (any(expected_obs<0.0))then
      print*,'obstype SIC expected concs:',expected_obs
      print*,'SIC sic negative lat lon:',llat,llon 
      endif
      if (any(expected_obs>1.0))then
      print*,'obstype SIC expected concs:',expected_obs
      print*,'SIC sic positive lat lon:',llat,llon 
      endif
endif

   

if (cat_signal == -1) then
      ! we need to know the aggregate sea ice concentration for these special cases
      base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_CONCENTR))  
      call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, obs_type, cat_signal, expected_aggr_conc, istatus)
      expected_obs = expected_obs/max(expected_aggr_conc,1.0e-8)  ! hope this is allowed so we never divide by zero
      
      if (any(expected_aggr_conc<0.0) .or. any(expected_aggr_conc>1.0))then
      print*,'obstype SIT expected conc:',expected_aggr_conc
      print*,'SIT sic lat lon:',llat,llon
      endif
    
endif

if (debug > 1) print *, 'interp val, istatus = ', expected_obs, istatus

end subroutine model_interpolate

!------------------------------------------------------------------

! Four different types of grids are used here.  
! The dipole grid is referred to as a dipole grid and each 
! region is referred to as a quad, short for quadrilateral.  
! The tripole grid uses the same interpolation scheme as the
! dipole grid, and so for most of the code we refer to
! the tripole grid as a dipole grid.
! The longitude/latitude rectangular grid with possibly irregular
! spacing in latitude used for some applications and testing
! is referred to as the irregular grid and each region is 
! called a box.
! Finally, a regularly spaced longitude/latitude grid is used
! as a computational tool for interpolating from the dipole
! grid. This is referred to as the regular grid and each region
! is called a box. 
! All grids are referenced by the index of the lower left corner
! of the quad or box.

! The dipole grid is assumed to be global for all applications.
! The irregular grid is also assumed to be global east
! west for all applications.

!------------------------------------------------------------------
! CMB notes: changed to deal with need to aggregate over categories
! offset is index of variable block w/o regard to lon, lat, or cat
subroutine lon_lat_interpolate(state_handle, ens_size, offset, lon, lat, var_type, cat_signal, expected_obs, istatus)
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer(i8),         intent(in)  :: offset 
real(r8),            intent(in)  :: lon, lat
integer,             intent(in)  :: var_type
integer,             intent(in)  :: cat_signal
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

! Subroutine to interpolate to a lon lat location given the state vector
! via state_handle. This works just on one horizontal slice.

! Local storage
integer  :: lat_bot, lat_top, lon_bot, lon_top, num_inds, start_ind
integer  :: x_ind, y_ind
real(r8) :: x_corners(4), y_corners(4)
real(r8) :: p(4,ens_size), xbot(ens_size), xtop(ens_size)
real(r8) :: work_expected_obs(ens_size)
real(r8) :: lon_fract, lat_fract
logical  :: masked
integer  :: quad_status
integer  :: e, iterations, Niterations
integer(i8) :: next_offset 

if ( .not. module_initialized ) call static_init_model

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
      call get_dipole_quad(lon, lat, ulon, ulat, num_inds, start_ind, &
         u_dipole_lon_list, u_dipole_lat_list, lon_bot, lat_bot, quad_status)
      ! Fail on bad istatus return
      if(quad_status /= 0) then
         istatus = quad_status
         return
      endif

      ! Getting corners for accurate interpolation
      call get_quad_corners(ulon, lon_bot, lat_bot, x_corners)
      call get_quad_corners(ulat, lon_bot, lat_bot, y_corners)

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

      ! If there are no quads overlapping, can't do interpolation
      if(num_inds == 0) then
         istatus = 1
         return
      endif

      call get_dipole_quad(lon, lat, tlon, tlat, num_inds, start_ind, &
         t_dipole_lon_list, t_dipole_lat_list, lon_bot, lat_bot, quad_status)

      ! Fail on bad istatus return
      if(quad_status /= 0) then
         istatus = quad_status
         return
      endif

      ! Fail if point is in T box that covers pole
      if(lon_bot == pole_x .and. lat_bot == t_pole_y) then
         istatus = 5
         return
      endif

      ! Getting corners for accurate interpolation
      call get_quad_corners(tlon, lon_bot, lat_bot, x_corners)
      call get_quad_corners(tlat, lon_bot, lat_bot, y_corners)

   endif

else
   ! This is an irregular grid
   ! U and V are on velocity grid
   if (is_on_ugrid(var_type)) then
      ! Get the corner indices and the fraction of the distance between
      call get_irreg_box(lon, lat, ulon, ulat, &
         lon_bot, lat_bot, lon_fract, lat_fract, quad_status)
   else
      ! Eta, T and S are on the T grid
      ! Get the corner indices
      call get_irreg_box(lon, lat, tlon, tlat, &
         lon_bot, lat_bot, lon_fract, lat_fract, quad_status)
   endif

   ! Return passing through error status
   if(quad_status /= 0) then
      istatus = quad_status
      return
   endif

endif

! Find the indices to get the values for interpolating
lat_top = lat_bot + 1
if(lat_top > ny) then
   istatus = 2
   return
endif

! Watch for wraparound in longitude
lon_top = lon_bot + 1
if(lon_top > nx) lon_top = 1

if ( cat_signal < 1 )  then
   Niterations = Ncat ! only iterate if aggregating over all types
else
   Niterations = 1 ! no need to iterate
endif

work_expected_obs = 0.0_r8
expected_obs = 0.0_r8

do iterations = 1, Niterations

   ! FIXME: this should use the state structure routine 'get_dart_vector_index'
   ! to get the start of the next category layer.  this code assumes it knows
   ! exactly how the state vector is laid out (reasonable, but might not be true
   ! in future versions of dart.)
   next_offset = offset + (iterations-1) * Nx * Ny 

   ! Get the values at the four corners of the box or quad
   ! Corners go around counterclockwise from lower left
   p(1, :) = get_val(lon_bot, lat_bot, nx, state_handle, next_offset, ens_size, var_type, masked)
   if(masked) then
      istatus = 3
      return
   endif
   
   p(2, :) = get_val(lon_top, lat_bot, nx, state_handle, next_offset, ens_size, var_type, masked)
   if(masked) then
      istatus = 3
      return
   endif
   
   p(3, :) = get_val(lon_top, lat_top, nx, state_handle, next_offset, ens_size, var_type, masked)
   if(masked) then
      istatus = 3
      return
   endif
   
   p(4, :) = get_val(lon_bot, lat_top, nx, state_handle, next_offset, ens_size, var_type, masked)
   if(masked) then
      istatus = 3
      return
   endif
   
   ! Full bilinear interpolation for quads
   if(dipole_grid) then
      do e = 1, ens_size
         call quad_bilinear_interp(lon, lat, x_corners, y_corners, p(:,e), ens_size, work_expected_obs(e))
      enddo
   else
      ! Rectangular bilinear interpolation - horizontal plane only
      xbot = p(1, :) + lon_fract * (p(2, :) - p(1, :))  ! bot is really south
      xtop = p(4, :) + lon_fract * (p(3, :) - p(4, :))  ! top is really north
      ! Now interpolate in latitude
      work_expected_obs = xbot + lat_fract * (xtop - xbot)
   endif

   expected_obs = expected_obs+work_expected_obs

enddo

end subroutine lon_lat_interpolate

!------------------------------------------------------------
! CMB changed so bathymetry check is land check only, height was only used 
! for masking below bathymetry, not needed in sea ice so removed

function get_val(lon_index, lat_index, nlon, state_handle, offset, &
                 ens_size, var_type, masked)

! Returns the index in state vecture structure of a variable at a single level 
! given the lat and lon indices. Used by lon_lat_interpolate only
! 'masked' returns true if this is NOT a valid grid location (e.g. land)
 integer,             intent(in)  :: lon_index, lat_index, nlon, var_type
 type(ensemble_type), intent(in)  :: state_handle
 integer(i8),         intent(in)  :: offset  ! index of the variable block
 integer,             intent(in)  :: ens_size
 logical,             intent(out) :: masked  ! T = land, F = ocean

 real(r8)    :: get_val(ens_size)
 integer(i8) :: state_index

if ( .not. module_initialized ) call static_init_model

! check for land and return if not valid water cell.
if(is_dry_land(var_type, lon_index, lat_index)) then
   masked = .true. ! cell is land
   get_val = MISSING_R8
   return
endif

masked = .false.  ! cell is water

!>@todo FIXME: this should call get_dart_vector_index() to convert from lat,lon,cat
! to offset in the state vector.
! this code assumes it knows the layout (lons varying most rapidly)
! state index must be 8byte integer
state_index = int(lat_index - 1,i8)*int(nlon,i8) + int(lon_index,i8) + int(offset-1,i8)

! this function gets the state vector value for each ensemble member
get_val = get_state(state_index, state_handle) 

end function get_val

!------------------------------------------------------------

subroutine get_irreg_box(lon, lat, lon_array, lat_array, &
                         found_x, found_y, lon_fract, lat_fract, istatus)

 real(r8),   intent(in) :: lon, lat
 real(r8),   intent(in) :: lon_array(nx, ny), lat_array(nx, ny)
 real(r8),  intent(out) :: lon_fract, lat_fract
 integer,   intent(out) :: found_x, found_y, istatus

! Given a longitude and latitude of a point (lon and lat) and the
! longitudes and latitudes of the lower left corner of the regular grid
! boxes, gets the indices of the grid box that contains the point and
! the fractions along each direction for interpolation.

! Local storage
integer  :: lat_status, lon_top, lat_top

! Succesful return has istatus of 0
istatus = 0

! Get latitude box boundaries
call lat_bounds(lat, ny, lat_array, found_y, lat_top, lat_fract, lat_status)

! Check for error on the latitude interpolation
if(lat_status /= 0) then
   istatus = 1
   return
endif

! Find out what longitude box and fraction
call lon_bounds(lon, nx, lon_array, found_x, lon_top, lon_fract)

end subroutine get_irreg_box

!------------------------------------------------------------

subroutine lon_bounds(lon, nlons, lon_array, bot, top, fract)
 real(r8),    intent(in) :: lon
 integer,     intent(in) :: nlons
 real(r8),    intent(in) :: lon_array(:, :)
 integer,    intent(out) :: bot, top
 real(r8),   intent(out) :: fract

! Given a longitude lon, the array of longitudes for grid boundaries, and the
! number of longitudes in the grid, returns the indices of the longitude
! below and above the location longitude and the fraction of the distance
! between. It is assumed that the longitude wraps around for a global grid. 
! Since longitude grids are going to be regularly spaced, this could be made more efficient.
! Algorithm fails for a silly grid that has only two longitudes separated by 180 degrees.

! Local storage
integer  :: i
real(r8) :: dist_bot, dist_top

if ( .not. module_initialized ) call static_init_model

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

!-------------------------------------------------------------

subroutine lat_bounds(lat, nlats, lat_array, bot, top, fract, istatus)
 real(r8),   intent(in) :: lat
 integer,    intent(in) :: nlats
 real(r8),   intent(in) :: lat_array(:, :)
 integer,   intent(out) :: bot, top
 real(r8),  intent(out) :: fract
 integer,   intent(out) :: istatus

! Given a latitude lat, the array of latitudes for grid boundaries, and the
! number of latitudes in the grid, returns the indices of the latitude
! below and above the location latitude and the fraction of the distance
! between. istatus is returned as 0 unless the location latitude is 
! south of the southernmost grid point (1 returned) or north of the 
! northernmost (2 returned). If one really had lots of polar obs would 
! want to worry about interpolating around poles.

! Local storage
integer :: i

if ( .not. module_initialized ) call static_init_model

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

!------------------------------------------------------------------

function lon_dist(lon1, lon2)
 real(r8), intent(in) :: lon1, lon2
 real(r8)             :: lon_dist

! Returns the smallest signed distance between lon1 and lon2 on the sphere
! If lon1 is less than 180 degrees east of lon2 the distance is negative
! If lon1 is less than 180 degrees west of lon2 the distance is positive

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

!------------------------------------------------------------

subroutine get_dipole_quad(lon, lat, qlons, qlats, num_inds, start_ind, &
                           x_inds, y_inds, found_x, found_y, istatus)

 real(r8), intent(in)  :: lon, lat, qlons(:, :), qlats(:, :)
 integer,  intent(in)  :: num_inds, start_ind, x_inds(:), y_inds(:)
 integer,  intent(out) :: found_x, found_y, istatus

! Given the lon and lat of a point, and a list of the
! indices of the quads that might contain a point at (lon, lat), determines
! which quad contains the point.  istatus is returned as 0 if all went 
! well and 1 if the point was not found to be in any of the quads.

integer :: i, my_index
real(r8) :: x_corners(4), y_corners(4)

istatus = 0

! Loop through all the quads and see if the point is inside
do i = 1, num_inds
   my_index = start_ind + i - 1
   call get_quad_corners(qlons, x_inds(my_index), y_inds(my_index), x_corners)
   call get_quad_corners(qlats, x_inds(my_index), y_inds(my_index), y_corners)

   ! Ssearch in this individual quad
   if(in_quad(lon, lat, x_corners, y_corners)) then
      found_x = x_inds(my_index)
      found_y = y_inds(my_index)
      return
   endif
enddo

! Falling off the end means search failed, return istatus 1
istatus = 1

end subroutine get_dipole_quad

!------------------------------------------------------------

function in_quad(lon, lat, x_corners, y_corners)
 real(r8), intent(in)  :: lon, lat, x_corners(4), y_corners(4)
 logical               :: in_quad

! Return in_quad true if the point (lon, lat) is in the quad with 
! the given corners.

! Do this by line tracing in latitude for now. For non-pole point, want a vertical
! line from the lon, lat point to intersect a side of the quad both above
! and below the point.

real(r8) :: x(2), y(2)
logical  :: cant_be_in_box, in_box
integer  :: intercepts_above(4), intercepts_below(4), i
integer  :: num_above, num_below

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

!------------------------------------------------------------

subroutine line_intercept(side_x_in, side_y, x_point_in, y_point, &
                          cant_be_in_box, in_box, intercept_above, intercept_below)

 real(r8), intent(in)  :: side_x_in(2), side_y(2), x_point_in, y_point
 logical,  intent(out) :: cant_be_in_box, in_box
 integer,  intent(out) :: intercept_above, intercept_below

! Find the intercept of a vertical line from point (x_point, y_point) and
! a line segment with endpoints side_x and side_y.
! For a given side have endpoints (side_x1, side_y1) and (side_x2, side_y2)
! so equation of segment is y = side_y1 + m(x-side_x1) for y 
! between side_y1 and side_y2.
! Intersection of vertical line and line containing side 
! occurs at y = side_y1 + m(x_point - side_x1); need this
! y to be between side_y1 and side_y2.
! If the vertical line is colinear with the side but the point is not on the side, return
! cant_be_in_box as true. If the point is on the side, return in_box true.
! If the intersection of the vertical line and the side occurs at a point above
! the given point, return 1 for intercept_above. If the intersection occurs 
! below, return 1 for intercept_below. If the vertical line does not intersect
! the segment, return false and 0 for all intent out arguments.

! WARNING: CERTAINLY A PROBLEM FOR THE POLE BOX!!! POLE BOX COULD
! HAVE SIDES THAT ARE LONGER THAN 180. For now pole boxes are excluded.

! This can probably be made much cleaner and more efficient.

real(r8) :: slope, y_intercept, side_x(2), x_point

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

!------------------------------------------------------------

subroutine quad_bilinear_interp(lon_in, lat, x_corners_in, y_corners, &
                                p, ens_size, expected_obs)

 real(r8),  intent(in) :: lon_in, lat, x_corners_in(4), y_corners(4), p(4)
 integer,   intent(in) :: ens_size
 real(r8), intent(out) :: expected_obs

! Given a longitude and latitude (lon_in, lat), the longitude and
! latitude of the 4 corners of a quadrilateral and the values at the
! four corners, interpolates to (lon_in, lat) which is assumed to
! be in the quad. This is done by bilinear interpolation, fitting
! a function of the form a + bx + cy + dxy to the four points and 
! then evaluating this function at (lon, lat). The fit is done by
! solving the 4x4 system of equations for a, b, c, and d. The system
! is reduced to a 3x3 by eliminating a from the first three equations
! and then solving the 3x3 before back substituting. There is concern
! about the numerical stability of this implementation. Implementation
! checks showed accuracy to seven decimal places on all tests.

integer :: i
real(r8) :: m(3, 3), v(3), r(3), a, x_corners(4), lon
! real(r8) :: lon_mean

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
expected_obs = a + r(1)*lon + r(2)*lat + r(3)*lon*lat

!********
! Avoid exceeding maxima or minima as stopgap for poles problem
! When doing bilinear interpolation in quadrangle, can get interpolated
! values that are outside the range of the corner values
if(expected_obs > maxval(p)) then
   expected_obs = maxval(p)
else if(expected_obs < minval(p)) then
   expected_obs = minval(p)
endif
!********

end subroutine quad_bilinear_interp

!------------------------------------------------------------

subroutine mat3x3(m, v, r)
 real(r8),  intent(in) :: m(3, 3), v(3)
 real(r8), intent(out) :: r(3)

! Solves rank 3 linear system mr = v for r
! using Cramer's rule. This isn't the best choice
! for speed or numerical stability so might want to replace
! this at some point.

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

!------------------------------------------------------------

function deter3(m)
 real(r8), intent(in) :: m(3, 3)
 real(r8)             :: deter3

! Computes determinant of 3x3 matrix m

deter3 = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + &
         m(1,3)*m(2,1)*m(3,2) - m(3,1)*m(2,2)*m(1,3) - &
         m(1,1)*m(2,3)*m(3,2) - m(3,3)*m(2,1)*m(1,2)

end function deter3

!------------------------------------------------------------------
!> Returns the the time step of the model; the smallest increment
!> in time that the model is capable of advancing the state in a given
!> implementation. This interface is required for all applications.

function shortest_time_between_assimilations()
 type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = model_timestep

end function shortest_time_between_assimilations

!------------------------------------------------------------------
!> Given an integer index (index_in) to point in state vector structure, returns the
!> associated location, which is an array of indices for lat, lon, and cat. 
!> A second intent(out) optional argument kind
!> can be returned if the model has more than one type of field (for
!> instance temperature and zonal wind component). This interface is
!> required for all filter applications as it is required for computing
!> the distance between observations and state variables.

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

real(r8) :: lat, lon, rcat
integer  :: lon_index, lat_index, cat_index, local_var, var_id

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, lon_index, lat_index, cat_index, var_id=var_id)
call get_state_kind(var_id, local_var)

if (is_on_ugrid(local_var)) then
   lon = ULON(lon_index, lat_index)
   lat = ULAT(lon_index, lat_index)
else
   lon = TLON(lon_index, lat_index)
   lat = TLAT(lon_index, lat_index)
endif

if (debug > 5) print *, 'lon, lat, cat_index = ', lon, lat, cat_index
rcat     = cat_index*1.0_r8
location = set_location(lon, lat, rcat, VERTISLEVEL)

if (present(var_type)) then
   var_type = local_var
   if(is_dry_land(var_type, lon_index, lat_index)) then
      var_type = QTY_DRY_LAND  
   endif
endif

end subroutine get_state_meta_data

!--------------------------------------------------------------------

function get_varid_from_kind(dart_kind)

integer, intent(in) :: dart_kind
integer             :: get_varid_from_kind

! given a kind, return what variable number it is

integer :: i

do i = 1, get_num_variables(domain_id)
   if (dart_kind == state_kinds_list(i)) then
      get_varid_from_kind = i
      return
   endif
end do

if (debug > 1) then
   write(string1, *) 'Kind ', dart_kind, ' not found in state vector'
   write(string2, *) 'AKA ', get_name_for_quantity(dart_kind), ' not found in state vector'
   call error_handler(E_MSG,'get_varid_from_kind', string1, &
                      source, revision, revdate, text2=string2)
endif

get_varid_from_kind = -1

end function get_varid_from_kind


!------------------------------------------------------------------

subroutine get_state_kind(var_ind, var_type)
 integer, intent(in)  :: var_ind
 integer, intent(out) :: var_type

! Given an integer index into the state vector structure, returns the kind,
! and both the starting offset for this kind, as well as the offset into
! the block of this kind.

if ( .not. module_initialized ) call static_init_model

var_type = state_kinds_list(var_ind)

end subroutine get_state_kind

!------------------------------------------------------------------

subroutine get_state_kind_inc_dry(index_in, var_type)
 integer(i8), intent(in)  :: index_in
 integer,     intent(out) :: var_type

! Given an integer index into the state vector structure, returns the
! type, taking into account the ocean bottom and dry land.

integer :: lon_index, lat_index, depth_index, var_id

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, lon_index, lat_index, depth_index, var_id=var_id)
call get_state_kind(var_id, var_type)

! if on land, replace type with dry land.
if(is_dry_land(var_type, lon_index, lat_index)) then
   var_type = QTY_DRY_LAND
endif

end subroutine get_state_kind_inc_dry

!------------------------------------------------------------------

subroutine end_model()

! Shutdown and clean-up.

! assume if one is allocated, they all were.  if no one ever
! called the init routine, don't try to dealloc something that
! was never alloc'd.
if (allocated(ULAT)) deallocate(ULAT, ULON, TLAT, TLON, KMT, KMU, HT, HU)

end subroutine end_model

!------------------------------------------------------------------
!> Writes the model-specific attributes to a netCDF file.
!> This includes coordinate variables and some metadata, but NOT
!> the model state vector.

subroutine nc_write_model_atts( ncid, domain_id ) 

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

! for the dimensions and coordinate variables
integer :: NlonDimID, NlatDimID, NcatDimID
integer :: ulonVarID, ulatVarID, tlonVarID, tlatVarID
integer :: KMTVarID, KMUVarID

! variables for the namelist output
integer, parameter :: MAXLINELEN = 128
character(len=8), parameter :: cice_namelist_file = 'cice_in'
character(len=MAXLINELEN), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_cice_namelist

integer :: i
character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncid', ncid

!-------------------------------------------------------------------------------
! make sure ncid refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

call nc_begin_define_mode(ncid)

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call nc_add_global_creation_time(ncid)

!-------------------------------------------------------------------------------
! Determine shape of auxiliary namelist
!-------------------------------------------------------------------------------

call find_textfile_dims(cice_namelist_file, nlines, linelen)
if (nlines > 0) then
  has_cice_namelist = .true.
else
  has_cice_namelist = .false.
endif

if (debug > 0) print *, 'cice namelist: nlines, linelen = ', nlines, linelen
  
if (has_cice_namelist) then 
   allocate(textblock(nlines))
   textblock = ''

   call nc_check(nf90_def_dim(ncid=ncid, name='linelength', &
                 len = MAXLINELEN, dimid = LineLenDimID), &
                 'nc_write_model_atts', 'def_dim linelength ')

   call nc_check(nf90_def_dim(ncid=ncid, name='nlines', &
                 len = nlines, dimid = nlinesDimID), &
                 'nc_write_model_atts', 'def_dim nlines ')

   call nc_check(nf90_def_var(ncid,name=cice_namelist_file, xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'nc_write_model_atts', 'def_var cice_in')

   call nc_check(nf90_put_att(ncid, nmlVarID, 'long_name',       &
                 'contents of cice_in namelist'), 'nc_write_model_atts', 'put_att cice_in long_name')
   call nc_check(nf90_put_att(ncid, nmlVarID, 'original file name',       &
                 cice_namelist_file), 'nc_write_model_atts', 'put_att cice_in filename')

endif

call nc_add_global_attribute(ncid, "model_source"  , source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate" , revdate )
call nc_add_global_attribute(ncid, "model"         , "CICE")

!----------------------------------------------------------------------------
! Define the new dimensions IDs
!----------------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid, name='i', &
       len = Nx, dimid = NlonDimID),'nc_write_model_atts', 'i def_dim '//trim(filename))
call nc_check(nf90_def_dim(ncid, name='j', &
       len = Ny, dimid = NlatDimID),'nc_write_model_atts', 'j def_dim '//trim(filename))
call nc_check(nf90_def_dim(ncid, name='k', &
       len = Ncat, dimid = NcatDimID),'nc_write_model_atts', 'k def_dim '//trim(filename))

!----------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and the Attributes
!----------------------------------------------------------------------------

! U,V Grid Longitudes
call nc_check(nf90_def_var(ncid,name='ULON', xtype=nf90_real, &
              dimids=(/ NlonDimID, NlatDimID /), varid=ulonVarID),&
              'nc_write_model_atts', 'ULON def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  ulonVarID, 'long_name', 'longitudes of U,V grid'), &
              'nc_write_model_atts', 'ULON long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  ulonVarID, 'cartesian_axis', 'X'),  &
              'nc_write_model_atts', 'ULON cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncid,  ulonVarID, 'units', 'degrees_east'), &
              'nc_write_model_atts', 'ULON units '//trim(filename))
call nc_check(nf90_put_att(ncid,  ulonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
              'nc_write_model_atts', 'ULON valid_range '//trim(filename))

! U,V Grid Latitudes
call nc_check(nf90_def_var(ncid,name='ULAT', xtype=nf90_real, &
              dimids=(/ NlonDimID, NlatDimID /), varid=ulatVarID),&
              'nc_write_model_atts', 'ULAT def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  ulatVarID, 'long_name', 'latitudes of U,V grid'), &
              'nc_write_model_atts', 'ULAT long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  ulatVarID, 'cartesian_axis', 'Y'),   &
              'nc_write_model_atts', 'ULAT cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncid,  ulatVarID, 'units', 'degrees_north'),  &
              'nc_write_model_atts', 'ULAT units '//trim(filename))
call nc_check(nf90_put_att(ncid,  ulatVarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
              'nc_write_model_atts', 'ULAT valid_range '//trim(filename))

! Tracer Grid Longitudes
call nc_check(nf90_def_var(ncid,name='TLON', xtype=nf90_real, &
              dimids=(/ NlonDimID, NlatDimID /), varid=tlonVarID),&
              'nc_write_model_atts', 'TLON def_var '//trim(filename))
call nc_check(nf90_put_att(ncid, tlonVarID, 'long_name', 'longitudes of tracer grid'), &
              'nc_write_model_atts', 'TLON long_name '//trim(filename))
call nc_check(nf90_put_att(ncid, tlonVarID, 'cartesian_axis', 'X'),   &
              'nc_write_model_atts', 'TLON cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncid, tlonVarID, 'units', 'degrees_east'),  &
              'nc_write_model_atts', 'TLON units '//trim(filename))
call nc_check(nf90_put_att(ncid, tlonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
              'nc_write_model_atts', 'TLON valid_range '//trim(filename))

! Tracer Grid (center) Latitudes
call nc_check(nf90_def_var(ncid,name='TLAT', xtype=nf90_real, &
              dimids= (/ NlonDimID, NlatDimID /), varid=tlatVarID), &
              'nc_write_model_atts', 'TLAT def_var '//trim(filename))
call nc_check(nf90_put_att(ncid, tlatVarID, 'long_name', 'latitudes of tracer grid'), &
              'nc_write_model_atts', 'TLAT long_name '//trim(filename))
call nc_check(nf90_put_att(ncid, tlatVarID, 'cartesian_axis', 'Y'),   &
              'nc_write_model_atts', 'TLAT cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncid, tlatVarID, 'units', 'degrees_north'),  &
              'nc_write_model_atts', 'TLAT units '//trim(filename))
call nc_check(nf90_put_att(ncid, tlatVarID, 'valid_range', (/ -90.0_r8, 90.0_r8 /)), &
              'nc_write_model_atts', 'TLAT valid_range '//trim(filename))

! Depth mask
call nc_check(nf90_def_var(ncid,name='KMT',xtype=nf90_int, &
              dimids= (/ NlonDimID, NlatDimID /), varid=KMTVarID), &
              'nc_write_model_atts', 'KMT def_var '//trim(filename))
call nc_check(nf90_put_att(ncid, KMTVarID, 'long_name', 'lowest valid depth index at grid centroids'), &
              'nc_write_model_atts', 'KMT long_name '//trim(filename))
call nc_check(nf90_put_att(ncid, KMTVarID, 'units', 'levels'),  &
              'nc_write_model_atts', 'KMT units '//trim(filename))
call nc_check(nf90_put_att(ncid, KMTVarID, 'positive', 'down'),  &
              'nc_write_model_atts', 'KMT units '//trim(filename))
call nc_check(nf90_put_att(ncid, KMTVarID, 'comment', &
               'more positive is closer to the center of the earth'),  &
              'nc_write_model_atts', 'KMT comment '//trim(filename))

! Depth mask
call nc_check(nf90_def_var(ncid,name='KMU',xtype=nf90_int, &
              dimids= (/ NlonDimID, NlatDimID /), varid=KMUVarID), &
              'nc_write_model_atts', 'KMU def_var '//trim(filename))
call nc_check(nf90_put_att(ncid, KMUVarID, 'long_name', 'lowest valid depth index at grid corners'), &
              'nc_write_model_atts', 'KMU long_name '//trim(filename))
call nc_check(nf90_put_att(ncid, KMUVarID, 'units', 'levels'),  &
              'nc_write_model_atts', 'KMU units '//trim(filename))
call nc_check(nf90_put_att(ncid, KMUVarID, 'positive', 'down'),  &
              'nc_write_model_atts', 'KMU units '//trim(filename))
call nc_check(nf90_put_att(ncid, KMUVarID, 'comment', &
               'more positive is closer to the center of the earth'),  &
              'nc_write_model_atts', 'KMU comment '//trim(filename))

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call nc_end_define_mode(ncid)

!----------------------------------------------------------------------------
! Fill the coordinate variables
!----------------------------------------------------------------------------

call nc_check(nf90_put_var(ncid, ulonVarID, ULON ), &
             'nc_write_model_atts', 'ULON put_var '//trim(filename))
call nc_check(nf90_put_var(ncid, ulatVarID, ULAT ), &
             'nc_write_model_atts', 'ULAT put_var '//trim(filename))
call nc_check(nf90_put_var(ncid, tlonVarID, TLON ), &
             'nc_write_model_atts', 'TLON put_var '//trim(filename))
call nc_check(nf90_put_var(ncid, tlatVarID, TLAT ), &
             'nc_write_model_atts', 'TLAT put_var '//trim(filename))
call nc_check(nf90_put_var(ncid, KMTVarID, KMT ), &
             'nc_write_model_atts', 'KMT put_var '//trim(filename))
call nc_check(nf90_put_var(ncid, KMUVarID, KMU ), &
             'nc_write_model_atts', 'KMU put_var '//trim(filename))

!-------------------------------------------------------------------------------
! Record the cice_in namelist
!-------------------------------------------------------------------------------

if (has_cice_namelist) then
   call file_to_text(cice_namelist_file, textblock)
   call nc_check(nf90_put_var(ncid, nmlVarID, textblock ), &
                 'nc_write_model_atts', 'put_var nmlVarID')
   deallocate(textblock)
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts


!-----------------------------------------------------------------------
!> writes the time of the current state  and (optionally) the time
!> that dictates the length of the forecast.
!>
!> @param ncfile_out name of the file
!> @param model_time the current time of the model state
!> @param adv_to_time the time in the future of the next assimilation.
!>

subroutine write_model_time(ncid, model_time, adv_to_time)

integer,         intent(in)           :: ncid
type(time_type), intent(in)           :: model_time
type(time_type), intent(in), optional :: adv_to_time

character(len=16), parameter :: routine = 'write_model_time'

integer :: io, varid, iyear, imonth, iday, ihour, imin, isec
integer :: seconds

if ( .not. module_initialized ) call static_init_model

if (present(adv_to_time)) then
   call get_date(adv_to_time, iyear, imonth, iday, ihour, imin, isec)
   write(string1,*)'CICE/DART not configured to advance CICE.'
   write(string2,*)'called with optional advance_to_time of'
   write(string3,'(i4.4,5(1x,i2.2))')iyear,imonth,iday,ihour,imin, isec
   call error_handler(E_ERR, routine, string1, &
              source, revision, revdate, text2=string2,text3=string3)
endif

call get_date(model_time, iyear, imonth, iday, ihour, imin, isec)

seconds = (ihour*60 + imin)*60 + isec

call nc_begin_define_mode(ncid)
call nc_add_global_attribute(ncid, 'nyr'   , iyear)
call nc_add_global_attribute(ncid, 'month' , imonth)
call nc_add_global_attribute(ncid, 'mday'  , iday)
call nc_add_global_attribute(ncid, 'sec'   , seconds)
call nc_end_define_mode(ncid)

end subroutine write_model_time


!------------------------------------------------------------------

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

 type(ensemble_type), intent(inout) :: state_ens_handle
 integer,             intent(in)    :: ens_size
 real(r8),            intent(in)    :: pert_amp
 logical,             intent(out)   :: interf_provided

! Perturbs state copies for generating initial ensembles.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding a perturbation to each model state 
! variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

integer     :: var_type
integer     :: j,i 
integer(i8) :: dart_index
real(r8)    :: random_number

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

if ( .not. module_initialized ) call static_init_model

interf_provided = .true.

! Initialize random number sequence
call init_random_seq(random_seq, my_task_id())

! only perturb the actual ocean cells; leave the land and
! ocean floor values alone.
do i=1,state_ens_handle%my_num_vars
   dart_index = state_ens_handle%my_vars(i)
   call get_state_kind_inc_dry(dart_index, var_type)
   do j=1, ens_size
      if (var_type /= QTY_DRY_LAND) then
         state_ens_handle%copies(j,i) = random_gaussian(random_seq, &
            state_ens_handle%copies(j,i), &
            model_perturbation_amplitude)

      endif
   enddo
enddo

end subroutine pert_model_copies

!------------------------------------------------------------------

! CMB took out bathymetry assuming KMU>0 and KMT>0 for ocn

function is_dry_land(obs_type, lon_index, lat_index)  
 integer, intent(in)  :: obs_type
 integer, intent(in)  :: lon_index, lat_index
 logical              :: is_dry_land

! returns true if this point is land

logical :: is_ugrid

if ( .not. module_initialized ) call static_init_model

is_dry_land = .FALSE.    ! start out thinking everything is wet.

is_ugrid = is_on_ugrid(obs_type)
if ((      is_ugrid .and. 0 == KMU(lon_index, lat_index)) .or. &
    (.not. is_ugrid .and. 0 == KMT(lon_index, lat_index))) then
   is_dry_land = .TRUE.
   return
endif

end function is_dry_land

!------------------------------------------------------------------

function is_on_ugrid(obs_type)

integer, intent(in) :: obs_type
logical             :: is_on_ugrid

!  returns true if U, V -- everything else is on T grid

if ( .not. module_initialized ) call static_init_model

is_on_ugrid = .FALSE.

if ((obs_type == QTY_U_SEAICE_COMPONENT)  .or.  &
    (obs_type == QTY_V_SEAICE_COMPONENT)) is_on_ugrid = .TRUE.

end function is_on_ugrid

!------------------------------------------------------------------

function all_corners_wet(obs_kind, lon_ind, lat_ind)  

integer, intent(in)  :: obs_kind, lon_ind, lat_ind
logical :: all_corners_wet

integer :: lon_ind_p1

! returns true only if all of the corners are land
 
! set to fail so we can return early.
all_corners_wet = .false. 

! Have to worry about wrapping in longitude but not in latitude
lon_ind_p1 = lon_ind + 1
if(lon_ind_p1 > nx) lon_ind_p1 = 1

if (is_dry_land(obs_kind, lon_ind,    lat_ind)) return
if (is_dry_land(obs_kind, lon_ind_p1, lat_ind)) return
if (is_dry_land(obs_kind, lon_ind_p1, lat_ind+1)) return
if (is_dry_land(obs_kind, lon_ind,    lat_ind+1)) return 

all_corners_wet = .true.

end function all_corners_wet

!------------------------------------------------------------------
!>  Write the grid to a netcdf file for checking.

subroutine write_grid_netcdf()

integer :: ncid, NlonDimID, NlatDimID, NcatDimID
integer :: nlon, nlat
integer :: ulatVarID, ulonVarID, TLATvarid, TLONvarid
integer :: KMTvarid, KMUvarid

integer :: dimids(2);

if ( .not. module_initialized ) call static_init_model

nlon = size(ULAT,1)
nlat = size(ULAT,2)

call nc_check(nf90_create('dart_grid.nc', NF90_CLOBBER, ncid),'write_grid_netcdf')

! define dimensions

call nc_check(nf90_def_dim(ncid, 'i', nlon, NlonDimID),'write_grid_netcdf')
call nc_check(nf90_def_dim(ncid, 'j', nlat, NlatDimID),'write_grid_netcdf')
call nc_check(nf90_def_dim(ncid, 'k', Ncat, NcatDimID),'write_grid_netcdf')

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

call nc_check(nf90_put_att(ncid,ulonVarID,'long_name','U,V grid lons'), &
                                                     'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,ulatVarID,'long_name','U,V grid lats'), &
                                                     'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,tlonVarID,'long_name','T grid lons'), &
                                                     'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,tlatVarID,'long_name','T grid lats'), &
                                                    'write_grid_netcdf')

call nc_check(nf90_enddef(ncid),'write_grid_netcdf')

! fill variables

call nc_check(nf90_put_var(ncid,  KMTvarid,  KMT),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid,  KMUvarid,  KMU),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid, ulatVarID, ULAT),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid, ulonVarID, ULON),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid, TLATvarid, TLAT),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid, TLONvarid, TLON),'write_grid_netcdf')

call nc_check(nf90_close(ncid),'write_grid_netcdf')

end subroutine write_grid_netcdf

!------------------------------------------------------------------

subroutine get_close_state(filt_gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                         num_close, close_indices, distances, state_handle)

type(get_close_type), intent(in)    :: filt_gc
type(location_type),  intent(inout) :: base_loc
integer,              intent(in)    :: base_type
type(location_type),  intent(inout) :: locs(:)
integer,              intent(in)    :: loc_qtys(:)
integer(i8),          intent(in)    :: loc_indx(:)
integer,              intent(out)   :: num_close
integer,              intent(out)   :: close_indices(:)
real(r8),             intent(out), optional :: distances(:)
type(ensemble_type),  intent(in),  optional :: state_handle

! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (locs, loc_qtys/indx), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

integer :: t_ind, k

! Initialize variables to missing status

num_close = 0
close_indices = -99
if (present(distances)) distances(:) = 1.0e9   !something big and positive (far away)

! Get all the potentially close obs but no dist (optional argument dist(:)
! is not present) This way, we are decreasing the number of distance
! computations that will follow.  This is a horizontal-distance operation and
! we don't need to have the relevant vertical coordinate information yet 
! (for obs).

call loc_get_close_state(filt_gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                         num_close, close_indices)

! Loop over potentially close subset of obs priors or state variables
if (present(distances)) then
   do k = 1, num_close

      t_ind = close_indices(k)

      ! if dry land, leave original 1e9 value.  otherwise, compute real dist.
      if (loc_qtys(t_ind) /= QTY_DRY_LAND) then
         distances(k) = get_dist(base_loc,      locs(t_ind), &
                                 base_type, loc_qtys(t_ind))
      endif
   
   enddo
endif

end subroutine get_close_state

!------------------------------------------------------------------

subroutine write_grid_interptest()

! Write the grid to an ascii file - in a format suitable for
! subsequent use in the 'test_interpolation()' code.
! write_grid_interptest is only possible after reading a real CICE grid,
! so static_init_model() must be called to gather the real CICE grid.

integer  :: i, j
real(r8) :: rowmat(Nx,1), colmat(1,Ny), dmat(Nx,Ny)
real(r8) :: rtlon, rulon, rtlat, rulat, u_val, t_val

if ( .not. module_initialized ) call static_init_model

!----------------------------------------------------------------------
! Generate a 'Regular' grid with the same rough 'shape' as the dipole grid
!----------------------------------------------------------------------

open(unit=12, position='rewind', action='write', file='regular_grid_u')
open(unit=13, position='rewind', action='write', file='regular_grid_t')
open(unit=14, position='rewind', action='write', file='regular_grid_u_data')
open(unit=15, position='rewind', action='write', file='regular_grid_t_data')

write(12, *) Nx, Ny
write(13, *) Nx, Ny

! Have T-grid starting at 0 and U grid offset by half
do i = 1, Nx
   rtlon = (i - 1.0_r8) * 360.0_r8 / Nx
   rulon = (i - 0.5_r8) * 360.0_r8 / Nx
   do j = 1, Ny
      rtlat = -90.0_r8 + (j - 1.0_r8) * 180.0_r8 / Ny
      rulat = -90.0_r8 + (j - 0.5_r8) * 180.0_r8 / Ny
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
! CICE grid (dipole) next
!----------------------------------------------------------------------

open(unit=12, position='rewind', action='write', file='dipole_grid_u')
open(unit=13, position='rewind', action='write', file='dipole_grid_t')
open(unit=14, position='rewind', action='write', file='dipole_grid_u_data')
open(unit=15, position='rewind', action='write', file='dipole_grid_t_data')

write(12, *) Nx, Ny
write(13, *) Nx, Ny

rowmat(:,1) = cos(PI * real((/ (i,i=0,Nx-1) /),r8) / Nx);
colmat(1,:) = sin(PI * real((/ (i,i=0,Ny-1) /),r8) / Ny);
dmat        = matmul(rowmat,colmat)

do i = 1, Nx
   do j = 1, Ny
      write(12, *) i, j, ULON(i,j), ULAT(i,j)
      write(13, *) i, j, TLON(i,j), TLAT(i,j)
      write(14, *)       ULON(i,j), ULAT(i,j), dmat(i, j)
      write(15, *)       TLON(i,j), TLAT(i,j), dmat(i, j)
   enddo
enddo

close(unit=12)
close(unit=13)
close(unit=14)
close(unit=15)

end subroutine write_grid_interptest

!------------------------------------------------------------------

subroutine test_interpolation(test_casenum)
 integer, intent(in) :: test_casenum

! Helen has destroyed this for now.

end subroutine test_interpolation

!--------------------------------------------------------------------

function read_model_time(filename)

character(len=256) :: filename
type(time_type) :: read_model_time

integer :: ncid         !< netcdf file id
integer :: nyr      , & ! year number, in cice restart
           month    , & ! month number, 1 to 12, in cice restart
           mday     , & ! day of the month, in cice restart
           sec          ! elapsed seconds into date, in cice restart
integer :: hour     , & ! hour of the day, needed for dart set_date
           minute   , & ! minute of the hour, needed for dart set_date
           secthismin 

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'read_model_time',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'read_model_time', 'open '//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'nyr'  , nyr), &
                  'read_model_time', 'get_att nyr')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'month' , month), &
                  'read_model_time', 'get_att month')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'mday'   , mday), &
                  'read_model_time', 'get_att mday')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'sec', sec), &
                  'read_model_time', 'get_att sec')

! FIXME: we don't allow a real year of 0 - add one for now, but
! THIS MUST BE FIXED IN ANOTHER WAY!
if (nyr == 0) then
  call error_handler(E_MSG, 'read_model_time', &
                     'WARNING!!!   year 0 not supported; setting to year 1')
  nyr = 1
endif

hour       = int(sec/3600) 
minute     = int((sec-hour*3600)/60)
secthismin = int(sec-hour*3600-minute*60) 

read_model_time = set_date(nyr, month, mday, hour, minute, secthismin)

end function read_model_time

!------------------------------------------------------------------
!> Verify that the namelist was filled in correctly, and check
!> that there are valid entries for the dart_kind. 
!> Returns a table with columns:  
!>
!>    netcdf_variable_name ; dart_kind_string ; update_string
!>

subroutine verify_state_variables( state_variables, ngood, table, kind_list, update_var )

character(len=*),  intent(inout) :: state_variables(:)
integer,           intent(out) :: ngood
character(len=*),  intent(out) :: table(:,:)
integer,           intent(out) :: kind_list(:)   ! kind number
logical, optional, intent(out) :: update_var(:) ! logical update

integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname, dartstr, update

if ( .not. module_initialized ) call static_init_model

nrows = size(table,1)

ngood = 0

!>@todo deprecate. Remove a hidden 'default' set of variables.
!>@     The default is provided in the input namelist.

if ( state_variables(1) == ' ' ) then ! no model_state_variables namelist provided
   call use_default_state_variables( state_variables )
   string1 = 'model_nml:model_state_variables not specified using default variables'
   call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate)
endif

MyLoop : do i = 1, nrows

   varname = trim(state_variables(3*i -2))
   dartstr = trim(state_variables(3*i -1))
   update  = trim(state_variables(3*i   ))
   
   call to_upper(update)

   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)
   table(i,3) = trim(update)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' .and. table(i,3) == ' ') exit MyLoop

   if ( table(i,1) == ' ' .or. table(i,2) == ' ' .or. table(i,3) == ' ' ) then
      string1 = 'model_nml:model_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   kind_list(i) = get_index_for_quantity(dartstr)
   if( kind_list(i)  < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif
   
   ! Make sure the update variable has a valid name

   if ( present(update_var) )then
      SELECT CASE (update)
         CASE ('UPDATE')
            update_var(i) = .true.
         CASE ('NO_COPY_BACK')
            update_var(i) = .false.
         CASE DEFAULT
            write(string1,'(A)')  'only UPDATE or NO_COPY_BACK supported in model_state_variable namelist'
            write(string2,'(6A)') 'you provided : ', trim(varname), ', ', trim(dartstr), ', ', trim(update)
            call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate, text2=string2)
      END SELECT
   endif

   ! Record the contents of the DART state vector

   if (do_output()) then
      write(string1,'(A,I2,6A)') 'variable ',i,' is ',trim(varname), ', ', trim(dartstr), ', ', trim(update)
      call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate)
   endif

   ngood = ngood + 1
enddo MyLoop

end subroutine verify_state_variables

!------------------------------------------------------------------
!> Default state_variables from model_mod.
!>@todo DEPRECATE

subroutine use_default_state_variables( state_variables )

character(len=*),  intent(inout) :: state_variables(:)

! strings must all be the same length for the gnu compiler
state_variables( 1:5*num_state_table_columns ) = &
   (/ 'CONCENTRATION             ', 'QTY_SEAICE_CONCENTR       ', 'UPDATE                    ', &
      'ICEVOLUME                 ', 'QTY_SEAICE_VOLUME         ', 'UPDATE                    ', &
      'SNOWVOLUME                ', 'QTY_SEAICE_SNOWVOLUME     ', 'UPDATE                    ', &
      'UICE                      ', 'QTY_U_SEAICE_COMPONENT    ', 'UPDATE                    ', &
      'VICE                      ', 'QTY_V_SEAICE_COMPONENT    ', 'UPDATE                    '/)

end subroutine use_default_state_variables

!------------------------------------------------------------------
! End of model_mod
!------------------------------------------------------------------

end module model_mod


! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
!----------------------------------------------------------------
!>
!> This is the interface between the ROMS ocean model and DART.
!> The required public interfaces arguments CANNOT be changed.
!>
!> Written in collaboration between Hernan Arango, Andy Moore, Chris Edwards
!> and the DART team. Uses/requires precomputed forward operator
!> values output directly from ROMS (i.e. the output of s4dvar.in:MODname).
!> The obs values are computed at the exact time of the obs as the model
!> is advancing (FGAT - First Guess at Appropriate Time). As a result,
!> there is currently NO model_interpolate() routine, since all
!> the observations already have expected values. What is required is
!> the ability to convert -on-demand- a DART observation sequence file
!> from the ROMS observation format. This is done with the convert_roms_obs
!> program in observations/ROMS.
!>
!> If required for other obs, the model_interpolate code will
!> need to be written and tested.
!>
!> The ROMS model uses a mask array to indicate dry land, so inside
!> DART and this model mod we can ignore dry land completely.  there
!> is no need to read the mask or account for it in get_state_meta_data
!> or get_close.
!>
!----------------------------------------------------------------

module model_mod

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, digits12, SECPERDAY, DEG2RAD, rad2deg, PI, &
                             MISSING_I, MISSING_R4, MISSING_R8, i4, i8, &
                             vtablenamelength

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time, &
                             print_time, print_date,                            &
                             set_calendar_type, get_calendar_type,              &
                             operator(*),  operator(+), operator(-),            &
                             operator(>),  operator(<), operator(/),            &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, set_location, get_location,         &
                             write_location, set_location_missing,              &
                             get_close_obs, get_close_state,                    &
                             convert_vertical_obs, convert_vertical_state,      &
                             VERTISHEIGHT, VERTISSURFACE, is_vertical

use    utilities_mod, only : register_module, error_handler, do_nml_term,       &
                             E_ERR, E_WARN, E_MSG, logfileunit, nmlfileunit,    &
                             get_unit, do_output, to_upper, do_nml_file,        &
                             find_namelist_in_file, check_namelist_read,        &
                             open_file, file_exist, find_textfile_dims,         &
                             file_to_text, do_output, close_file,               &
                             string_to_real, string_to_logical

use     obs_kind_mod, only : QTY_TEMPERATURE,           &
                             QTY_SALINITY,              &
                             QTY_U_CURRENT_COMPONENT,   &
                             QTY_V_CURRENT_COMPONENT,   &
                             QTY_SEA_SURFACE_HEIGHT,    &
                             QTY_SEA_SURFACE_PRESSURE,  &
                             QTY_POTENTIAL_TEMPERATURE, &
                             get_index_for_quantity,    &
                             get_name_for_quantity

use     mpi_utilities_mod, only : my_task_id

use        random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use  ensemble_manager_mod, only : ensemble_type

use distributed_state_mod, only : get_state

use   state_structure_mod, only : add_domain, get_model_variable_indices, &
                                  get_num_variables, get_index_start, &
                                  get_num_dims, get_domain_size, get_varid_from_kind, &
                                  get_dart_vector_index, state_structure_info, &
                                  get_index_start, get_index_end, get_variable_name, &
                                  get_kind_index, get_kind_string, get_dim_length, &
                                  get_dim_name, get_missing_value, get_units, &
                                  get_long_name, get_xtype, get_has_missing_value, &
                                  get_dim_lengths

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_end_define_mode, nc_check

use location_io_mod,      only : nc_write_location_atts, nc_get_location_varids, &
                                 nc_write_location

use default_model_mod,    only : pert_model_copies, nc_write_model_vars, init_conditions, &
                                 init_time, adv_1step

use quad_utils_mod,       only : quad_interp_handle, print_quad_handle, set_quad_coords, &
                                 init_quad_interp, finalize_quad_interp, &
                                 quad_lon_lat_locate, quad_lon_lat_evaluate, &
                                 GRID_QUAD_FULLY_IRREGULAR, QUAD_LOCATED_CELL_CENTERS, &
                                 QUAD_LOCATED_LON_EDGES, QUAD_LOCATED_LAT_EDGES

use typesizes
use netcdf

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.

! routines in this list have code in this module
public :: get_model_size,                &
          get_state_meta_data,           &
          model_interpolate,             &
          shortest_time_between_assimilations, &
          static_init_model,             &
          end_model,                     &
          nc_write_model_atts,           &
          write_model_time,              &
          read_model_time

! code for these routines are in other modules
public :: nc_write_model_vars,           &
          pert_model_copies,             &
          adv_1step,                     &
          init_time,                     &
          init_conditions,               &
          convert_vertical_obs,          &
          convert_vertical_state,        &
          get_close_obs,                 &
          get_close_state

! not required interfaces but useful for utility programs
public :: get_time_information!,          &
!          get_location_from_ijk

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2, string3
logical, save :: module_initialized = .false.

type(quad_interp_handle) :: ugrid_handle, vgrid_handle, tgrid_handle

! things which can/should be in the model_nml
!>@todo FIXME ... replace remaining references to VERTISHEIGHT with vert_localization_coord
integer  :: assimilation_period_days     = 1
integer  :: assimilation_period_seconds  = 0
integer  :: vert_localization_coord      = VERTISHEIGHT
integer  :: debug = 0   ! turn up for more and more debug messages
character(len=256) :: roms_grid_file     = 'roms_grid.nc'
character(len=256) :: roms_filename      = 'roms_input.nc'
integer  :: vert_transform               = 4

namelist /model_nml/  &
   assimilation_period_days,    &
   assimilation_period_seconds, &
   roms_grid_file,              &
   roms_filename,               &
   vert_transform,              &
   vert_localization_coord,     &
   debug,                       &
   variables

! DART contents are specified in the input.nml:&model_nml namelist.
!>@todo  NF90_MAX_NAME is 256 ... this makes the namelist output unreadable
integer, parameter :: MAX_STATE_VARIABLES = 8
integer, parameter :: num_state_table_columns = 5
character(len=vtablenamelength) :: variables(MAX_STATE_VARIABLES * num_state_table_columns ) = ' '
character(len=vtablenamelength) :: var_names(MAX_STATE_VARIABLES) = ' '
logical  ::                   update_list(MAX_STATE_VARIABLES) = .FALSE.
integer  ::                     kind_list(MAX_STATE_VARIABLES) = MISSING_I
real(r8) ::                    clamp_vals(MAX_STATE_VARIABLES,2) = MISSING_R8

integer :: nfields   ! This is the number of variables in the DART state vector.

integer :: domain_id ! global variable for state_structure_mod routines

! Grid parameters
! nx, ny and nz are the size of the rho grids.
integer :: Nx = -1, Ny = -1, Nz = -1

integer :: Nxi_rho
integer :: Nxi_u
integer :: Nxi_v
integer :: Neta_rho
integer :: Neta_u
integer :: Neta_v
integer :: Ns_rho
integer :: Ns_w

! Vertical grid parameters
real(r8) :: theta_s, theta_b
real(r8) :: Tcline,  hc
 
!>@todo FIXME ... nancy suggested creating pointers for each of these so
!    we could simply use the myvarid as the index in the pointer ...

!>@todo FIXME ... technically, there should be separate BATHY variables for each
!grid. Right now, there is a function in get_state_meta_data() that estimates
!the bathymetry on the U and V grids given the rho bathymetry. Similar problems
!with ssh (for vertical interpolation)

real(r8), allocatable, target :: ULAT(:,:), ULON(:,:), &
                                 TLAT(:,:), TLON(:,:), &
                                 VLAT(:,:), VLON(:,:), &
                                 BATHY(:,:)
logical, allocatable, target :: UMASK(:,:), VMASK(:,:), TMASK(:,:)

type(time_type) :: model_timestep

integer :: model_size    ! the state vector length

contains


!-----------------------------------------------------------------------
! All the REQUIRED interfaces come first - by convention.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>
!> Returns the size of the DART state vector (i.e. model) as an integer.
!> Required for all applications.
!>

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



!-----------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location. A second intent(out) optional argument kind
!> can be returned if the model has more than one type of field (for
!> instance temperature and zonal wind component). This interface is
!> required for all filter applications as it is required for computing
!> the distance between observations and state variables.
!>
!> @param index_in the index into the DART state vector
!> @param location the location at that index
!> @param var_type the DART KIND at that index
!>

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

! Local variables

integer  :: iloc, vloc, jloc
integer  :: myvarid, myqty
real(r8) :: mybathy
real(r8) :: depths(Ns_rho)

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, iloc, jloc, vloc, var_id=myvarid)

myqty = get_kind_index(domain_id, myvarid)


!>@todo Add the SSH, for now put at 0
if (myqty == QTY_U_CURRENT_COMPONENT) then
   mybathy = 0.5*(BATHY(iloc,jloc)+BATHY(iloc+1,jloc))
   call get_depths(mybathy, 0.0_r8, depths, 1, Ns_rho)
   location = set_location(ULON(iloc,jloc), ULAT(iloc,jloc), depths(vloc), VERTISHEIGHT)

elseif (myqty == QTY_V_CURRENT_COMPONENT) then
   mybathy = 0.5*(BATHY(iloc,jloc)+BATHY(iloc,jloc+1))
   call get_depths(mybathy, 0.0_r8, depths, 1, Ns_rho)
   location = set_location(VLON(iloc,jloc), VLAT(iloc,jloc), depths(vloc), VERTISHEIGHT)

elseif (myqty == QTY_SEA_SURFACE_HEIGHT) then
   location = set_location(TLON(iloc,jloc), TLAT(iloc,jloc), 0.0_r8, VERTISSURFACE)

else  ! Everything else is assumed to be on the rho points
   call get_depths(BATHY(iloc,jloc), 0.0_r8, depths, 1, Ns_rho)
   location = set_location(TLON(iloc,jloc), TLAT(iloc,jloc), depths(vloc), VERTISHEIGHT)

endif

! return state quantity for this index if requested
if (present(var_type)) var_type = myqty

end subroutine get_state_meta_data


!-----------------------------------------------------------------------
!>
!> Model interpolate will interpolate any DART state variable
!> (i.e. S, T, U, V, Eta) to the given location given a state vector.
!> The type of the variable being interpolated is obs_type since
!> normally this is used to find the expected value of an observation
!> at some location. The interpolated value is returned in expected_obs
!> and istatus is 0 for success. NOTE: This is a workhorse routine and is
!> the basis for all the forward observation operator code.
!>
!> @param state_handle DART ensemble handle
!> @param ens_size DART ensemble size
!> @param location the location of interest
!> @param obs_type the DART KIND of interest
!> @param expected_obs the estimated value of the DART state at the location
!>          of interest (the interpolated value).
!> @param istatus interpolation status ... 0 == success, /=0 is a failure
!>

subroutine model_interpolate(state_handle, ens_size, location, obs_type, &
                             expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_type
real(r8),            intent(out) :: expected_obs(:)
integer,             intent(out) :: istatus(:)

! Local storage
integer       :: icorn, imem, ilev, N_lev_un
integer       :: lstatus, hstatus
integer       :: Ns_var
integer       :: var_id, ssh_id
integer       :: hgt_bot(ens_size), hgt_top(ens_size), hgt_bot_un(2*ens_size)
integer       :: lon_corner(4), lat_corner(4)
integer(i8)   :: dart_idx
real(r8)      :: loc_array(3), llon, llat, lheight
real(r8)      :: top_val, bot_val
real(r8)      :: val_corners(4,ens_size)
real(r8)      :: val_alldepths(Ns_rho,ens_size)
real(r8)      :: ssh_corners(ens_size), hgt_fract(ens_size)
real(r8)      :: depths_tmp(Ns_rho)

if ( .not. module_initialized ) call static_init_model

! Successful istatus is 0
expected_obs = MISSING_R8
istatus = 99

var_id = get_varid_from_kind(domain_id, obs_type)

! Get the individual locations values
loc_array = get_location(location)
llon    = loc_array(1)
llat    = loc_array(2)
lheight = loc_array(3)

if (debug > 1) print *, 'requesting interpolation of ', obs_type, ' at ', &
                         llon, llat, lheight

! kind (in-situ) temperature is a combination of potential temp,
! salinity, and pressure based on depth.  call a routine that
! interpolates all three, does the conversion, and returns the
! sensible/in-situ temperature.
if(obs_type == QTY_TEMPERATURE) then
   ! we know how to interpolate this from potential temp,
   ! salinity, and pressure based on depth.
   call compute_temperature(state_handle, ens_size, llon, llat, lheight, &
                            expected_obs, istatus)
   if (debug > 1) print *, 'expected_obs, istatus = ', expected_obs, istatus
   return
endif

! Find horizontal corners
if(obs_type == QTY_U_CURRENT_COMPONENT) then
   call quad_lon_lat_locate(ugrid_handle, llon, llat, lon_corner, lat_corner, lstatus)
elseif (obs_type == QTY_V_CURRENT_COMPONENT) then
   call quad_lon_lat_locate(vgrid_handle, llon, llat, lon_corner, lat_corner, lstatus)
else
   call quad_lon_lat_locate(tgrid_handle, llon, llat, lon_corner, lat_corner, lstatus)
endif

!PRINT *, "RE: (model_interpolate) lstatus", lstatus

if (lstatus /= 0) return


!>@todo check to see if any of the corners are on dry land - this requires
!knowledge of which land mask to use 

! For Sea Surface values don't need the vertical coordinate
if( is_vertical(location, "SURFACE") ) then
   !>@todo HK CHECK surface observations
   ! Get dart vector indices
   if (get_num_dims(domain_id, var_id)>2) then
      Ns_var = get_dim_length(domain_id, var_id, 3)
   else
      Ns_var = 1
   endif
   do icorn = 1,4
      dart_idx = get_dart_vector_index(lon_corner(icorn), lat_corner(icorn), &
                                       Ns_var, domain_id, var_id)
      val_corners(icorn,:) = get_state(dart_idx, state_handle)
   enddo
else
   ssh_id = get_varid_from_kind(domain_id, QTY_SEA_SURFACE_HEIGHT)
   do icorn = 1,4
      ! get SSH at corners location for all members
      dart_idx = get_dart_vector_index(lon_corner(icorn), lat_corner(icorn), 1, &
                                       domain_id, ssh_id)
      ssh_corners = get_state(dart_idx, state_handle)
      ! Loop on ensemble members to get all indices of depth
      do imem = 1,ens_size
         call get_depths(BATHY(lon_corner(icorn), lat_corner(icorn)), &
                         ssh_corners(imem), depths_tmp, 1, Ns_rho)
         call height_bounds(lheight, Ns_rho, depths_tmp, &
                            hgt_bot(imem), hgt_top(imem), hgt_fract(imem), hstatus)
         if (hstatus > 0) then
            istatus = hstatus
            return
         endif
      enddo

      ! Get the levels needed for the vertical interpolation
      call unique_vec_from_two(hgt_bot, hgt_top, ens_size, hgt_bot_un, N_lev_un)
      do ilev = 1, N_lev_un
         dart_idx = get_dart_vector_index(lon_corner(icorn), lat_corner(icorn), &
                                          hgt_bot_un(ilev), domain_id, var_id)
         val_alldepths(hgt_bot_un(ilev),:) = get_state(dart_idx, state_handle)
      enddo
      ! Loop on ensemble members to do the vertical interpolation and get values at
      ! corners
      do imem = 1,ens_size
         top_val = val_alldepths(hgt_top(imem),imem)
         bot_val = val_alldepths(hgt_bot(imem),imem)
         val_corners(icorn,imem) = bot_val + hgt_fract(imem) * (top_val - bot_val)
      enddo
   enddo
endif



! Do the horizontal interpolation
if(obs_type == QTY_U_CURRENT_COMPONENT) then
   call quad_lon_lat_evaluate(ugrid_handle, llon, llat, lon_corner, lat_corner,&
                              ens_size, val_corners, expected_obs, lstatus)
elseif(obs_type == QTY_V_CURRENT_COMPONENT) then
   call quad_lon_lat_evaluate(vgrid_handle, llon, llat, lon_corner, lat_corner,&
                              ens_size, val_corners, expected_obs, lstatus)
else
   call quad_lon_lat_evaluate(tgrid_handle, llon, llat, lon_corner, lat_corner,&
                              ens_size, val_corners, expected_obs, lstatus)
endif

istatus = lstatus

end subroutine model_interpolate


!-----------------------------------------------------------------------
!>
!> Returns the sigma level depths for a given bathymetry (H) and 
!> sea surface elevation (zeta).
!>


subroutine get_depths(H, zeta, depths, vtype, N_vert)

real(r8), intent(in)  :: H              ! Bathymetry 
real(r8), intent(in)  :: zeta           ! Sea surface elevation
integer,  intent(in)  :: N_vert         ! Number of vertical levels
integer,  intent(in)  :: vtype          ! type of data 
real(r8), intent(out) :: depths(N_vert) ! output

! Local storage
integer  :: ik                      ! index for loops
real(r8) :: sc(N_vert),Cs(N_vert)   ! s-curves in domain [-1 < sc < 0] at vertical RHO-points.
real(r8) :: ds 
real(r8) :: h2,z0                   ! Useful temporary variables


! Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
if (vert_transform == 2) then
    if (vtype == 2) then
        ds=1.0_r8/(N_vert-1)
        sc(1) = -1.0_r8;
        sc(N_vert) = 0.0_r8;
        Cs(1) = -1.0_r8;
        Cs(N_vert) = 0.0_r8;
        do ik = 2,N_vert-1
            sc(ik) = ds*(ik-2-N_vert)
        enddo
        Cs=csf(sc, N_vert)
   else
        ds=1.0_r8/N_vert
        do ik = 1,N_vert
            sc(ik)= ds*(ik-N_vert-0.5_r8)
        enddo
        Cs=csf(sc,N_vert);
    endif
else
    if (vtype == 2) then
        do ik = 1,N_vert
            sc(ik) = (ik-2-N_vert)/(N_vert-1)
        enddo
    else
        do ik = 1,N_vert
            sc(ik)=(ik-N_vert-0.1_r8)/N_vert
        enddo
    endif
    Cs = (1.0_r8 - theta_b) * (1.0_r8/sinh(theta_s)) * sinh(theta_s*sc)+theta_b * &
         ((0.5_r8/tanh(0.5_r8*theta_s)) * tanh(theta_s*(sc+0.5_r8))-0.5_r8)
endif

! Create S-coordinate system: based on model topography h(i,j),
! fast-time-averaged free-surface field and vertical coordinate
! transformation metrics compute evolving depths of of the three-
! dimensional model grid.

if (vert_transform == 2) then
    h2=(H+hc)
    do ik = 1,N_vert
        z0=hc*sc(ik)+Cs(ik)*H
        depths(ik) = z0*H/h2 + zeta*(1.0_r8+z0*1.0_r8/h2)
    enddo
else
    do ik = 1,N_vert
        z0=hc*(sc(ik)-Cs(ik))+Cs(ik)*H;
        depths(ik) = z0 + zeta*(1.0_r8+z0*1.0_r8/H)
    enddo
endif

! for DART, depth is positive:
depths = - depths


end subroutine get_depths


!-----------------------------------------------------------------------
!>
!> Computation of C function for vertical levels
!>

function csf (sc, N_vert)

real(r8), intent(in) :: sc(N_vert)   !s-curves in domain [-1 < sc < 0]
integer,  intent(in) :: N_vert       !size of the input variable

! variables declaration
real(r8), dimension(N_vert) :: csf                 !returned value
integer(r8)                 :: ik                  !index for loops
real(r8), dimension(N_vert) :: csurf               ! Csurface

IF (theta_s > 0.0_r8 ) THEN
    DO ik = 1,N_vert
        csurf(ik)=(1.0_r8-COSH(sc(ik)*theta_s))/(COSH(theta_s)-1.0_r8)
    ENDDO
ELSE
    DO ik = 1,N_vert
        csurf(ik)=-(sc(ik)**2)
    ENDDO
ENDIF
IF (theta_b > 0.0_r8) THEN
    DO ik = 1,N_vert
        csf(ik) = (EXP(theta_b*csurf(ik))-1.0_r8)/(1.0_r8-EXP(-theta_b))
    ENDDO
ELSE
    csf  = csurf
ENDIF

end function csf



!-----------------------------------------------------------------------
!>
!> Get the top and bottom levels for a given depth.
!>

subroutine height_bounds(lheight, nheights, hgt_array, bot, top, fract, istatus)

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
istatus = 0

! The zc array contains the depths of the center of the vertical grid boxes

! It is assumed that the top box is shallow and any observations shallower
! than the depth of this boxes center are just given the value of the
! top box.
if(lheight <= hgt_array(nheights)) then
   top = nheights
   bot = nheights - 1
   ! NOTE: the fract definition is the relative distance from bottom to top
   ! ??? Make sure this is consistent with the interpolation
   fract = 1.0_r8
   return
endif

! Search through the boxes
do i = nheights-1, 1, -1
   ! If the location is shallower than this entry, it must be in this box
   if(lheight <= hgt_array(i)) then
      top = i +1
      bot = i
      fract = (lheight - hgt_array(bot)) / (hgt_array(top) - hgt_array(bot))
      return
   endif
end do

! Falling off the end means the location is lower than the deepest height
! Fail with istatus 2 in this case
istatus = 2

end subroutine height_bounds




!-----------------------------------------------------------------------
!>
!> Interpolate the rho centered variable var_rho to a field at ctype 
!> points (U or V)
!>


function  rho2u(var_rho, ctype, size_x, size_y)
!!---------------------------------------------------------------------
!!                  ***  FUNCTION  rho2u  ***
!!
!! ** Purpose : interpole the rho centered variable var_rho
!!                          to a field at ctype points (U or V)
!!              size_x,size_y are the size of the 2D variable
!!
!! ** Method  : straight forward function
!!
!!---------------------------------------------------------------------
INTEGER, INTENT(in)                      :: size_x,size_y  ! size of the input 2D variable
REAL(r8), DIMENSION(size_x,size_y)   :: var_rho    ! 3D REAL 4 with the variable
CHARACTER(LEN=1)                         :: ctype        ! position wanted of the variable
REAL(r8), DIMENSION(:,:),  ALLOCATABLE   :: rho2u   ! 3D REAL 4 with the modified variable

INTEGER(r8)                              :: ij,ii   ! indices for loops

!!----------------------------------------------------------

SELECT CASE (ctype )
   CASE ('U','u')
      ALLOCATE ( rho2u(size_x-1,size_y))
      DO ij=1,size_y
         DO ii=1,size_x-1
            rho2u(ii,ij)=0.5_r8*(var_rho(ii,ij)+var_rho(ii+1,ij))
         ENDDO
      ENDDO
   CASE ('V','v')
      ALLOCATE ( rho2u(size_x,size_y-1))
      DO ij=1,size_y-1
         DO ii=1,size_x
            rho2u(ii,ij)=0.5_r8*(var_rho(ii,ij)+var_rho(ii,ij+1))
         ENDDO
      ENDDO
END SELECT

end function rho2u


!-----------------------------------------------------------------------
!>
!> Return vector with unique values of vect
!>
!>@todo this function appears to be unused ... remove

subroutine unique_vec(val, N, val_un, N_un)

integer, intent(in)  :: N
integer, intent(in)  :: val(N)
integer, intent(out) :: val_un(N)
integer, intent(out) :: N_un

! Local variables
integer :: min_val, max_val, i

val_un = -1000

min_val = minval(val)-1
max_val = maxval(val)
do while (min_val<max_val)
   i = i+1
   min_val = minval(val, mask=val>min_val)
   val_un(i) = min_val
enddo
N_un = i

end subroutine unique_vec


!-----------------------------------------------------------------------
!>
!> Return vector with unique values of vect
!>

subroutine unique_vec_from_two(val1, val2, N, val_un, N_un)

integer, intent(in)  :: N
integer, intent(in)  :: val1(N), val2(N)
integer, intent(out) :: val_un(:)
integer, intent(out) :: N_un

! Local variables
integer :: min_val, max_val, i
integer :: val(2*N)

val(1:N) = val1
val((N+1):2*N) = val2

val_un = -1000

min_val = minval(val)-1
max_val = maxval(val)
i = 0
do while (min_val<max_val)
   i = i+1
   min_val = minval(val, mask=val>min_val)
   val_un(i) = min_val
enddo
N_un = i

end subroutine unique_vec_from_two


!------------------------------------------------------------------
!> use potential temp, depth, and salinity to compute a sensible (in-situ)
!> temperature

subroutine compute_temperature(state_handle, ens_size, llon, llat, lheight, &
                               expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
real(r8),            intent(in)  :: llon, llat, lheight
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

! local variables
integer     :: icorn, imem, ilev, N_lev_un
integer     :: ssh_id, temp_id, salt_id
integer     :: lstatus, hstatus
integer     :: lon_corner(4), lat_corner(4)
integer     :: hgt_bot(ens_size), hgt_top(ens_size), hgt_bot_un(2*ens_size)
integer(i8) :: dart_idx
real(r8)    :: top_val, bot_val, pres_val
real(r8)    :: salt_val(ens_size), potential_temp(ens_size)
real(r8)    :: ssh_corners(ens_size), hgt_fract(ens_size)
real(r8)    :: depths_tmp(Ns_rho)
real(r8)    :: temp_alldepths(Ns_rho,ens_size)
real(r8)    :: salt_alldepths(Ns_rho,ens_size)
real(r8)    :: temp_corners(4,ens_size), salt_corners(4,ens_size)

expected_obs(:) = MISSING_R8
istatus = 98

! Get variable indices in state vector
ssh_id = get_varid_from_kind(domain_id, QTY_SEA_SURFACE_HEIGHT)
temp_id = get_varid_from_kind(domain_id, QTY_TEMPERATURE)
salt_id = get_varid_from_kind(domain_id, QTY_SALINITY)

! Find pressure equivalent to depth
pres_val = dpth2pres(lheight)

! Find horizontal corners
call quad_lon_lat_locate(tgrid_handle, llon, llat, lon_corner, lat_corner, lstatus)

if (lstatus /= 0) return

! Get values of T and S at corners
do icorn = 1,4
   ! get SSH at corners location for all members
   dart_idx = get_dart_vector_index(lon_corner(icorn), lat_corner(icorn), 1, &
                                    domain_id, ssh_id)
   ssh_corners = get_state(dart_idx, state_handle)
   do imem = 1,ens_size
      call get_depths(BATHY(lon_corner(icorn), lat_corner(icorn)), &
                      ssh_corners(imem), depths_tmp, 1, Ns_rho)
      call height_bounds(lheight, Ns_rho, depths_tmp, hgt_bot(imem), &
                         hgt_top(imem), hgt_fract(imem), hstatus)
      if (hstatus > 0) then
         istatus = hstatus
         return
      endif
   enddo
   ! Get the levels needed for the vertical interpolation
   call unique_vec_from_two(hgt_bot, hgt_top, ens_size, hgt_bot_un, N_lev_un)
   do ilev = 1, N_lev_un
      dart_idx = get_dart_vector_index(lon_corner(icorn), lat_corner(icorn), &
                                       hgt_bot_un(ilev),domain_id, temp_id)
      temp_alldepths(hgt_bot_un(ilev),:) = get_state(dart_idx, state_handle)
      dart_idx = get_dart_vector_index(lon_corner(icorn), lat_corner(icorn), &
                                       hgt_bot_un(ilev),domain_id, salt_id)
      salt_alldepths(hgt_bot_un(ilev),:) = get_state(dart_idx, state_handle)
   enddo
   ! Loop on ensemble members to do the vertical interpolation and get values at
   ! corners
   do imem = 1,ens_size
      top_val = temp_alldepths(hgt_top(imem),imem)
      bot_val = temp_alldepths(hgt_bot(imem),imem)
      temp_corners(icorn,imem) = bot_val + hgt_fract(imem) * (top_val - bot_val)
      top_val = salt_alldepths(hgt_top(imem),imem)
      bot_val = salt_alldepths(hgt_bot(imem),imem)
      salt_corners(icorn,imem) = bot_val + hgt_fract(imem) * (top_val - bot_val)
   enddo
enddo

! Do the horizontal interpolation: get T and S at position
call quad_lon_lat_evaluate(tgrid_handle, llon, llat, lon_corner, lat_corner, &
                           ens_size, temp_corners, potential_temp, lstatus)
call quad_lon_lat_evaluate(tgrid_handle, llon, llat, lon_corner, lat_corner, &
                           ens_size, salt_corners, salt_val, lstatus)

if (lstatus > 0) then
   istatus = lstatus
   return
endif

! and finally, convert to sensible (in-situ) temperature.
! potential temp in degrees C, pressure in decibars, salinity in psu or pss
! (g/kg).
do imem = 1, ens_size !> @todo should this vectorize inside insitu_temp?
   call insitu_temp(potential_temp(imem), salt_val(imem), pres_val*10.0_r8, &
                      expected_obs(imem))
enddo

istatus = 0

end subroutine compute_temperature


!-----------------------------------------------------------------------
!>
!> This function computes pressure in bars from depth in meters
!> using a mean density derived from depth-dependent global
!> average temperatures and salinities from levitus 1994, and
!> integrating using hydrostatic balance.
!>
!> references:
!>
!> levitus, s., r. burgett, and t.p. boyer, world ocean atlas
!> volume 3: salinity, noaa atlas nesdis 3, us dept. of commerce, 1994.
!>
!> levitus, s. and t.p. boyer, world ocean atlas 1994, volume 4:
!> temperature, noaa atlas nesdis 4, us dept. of commerce, 1994.
!>
!> dukowicz, j. k., 2000: reduction of pressure and pressure
!> gradient errors in ocean simulations, j. phys. oceanogr., submitted.
!>

function dpth2pres(depth)

real(r8), intent(in)  :: depth
real(r8)              :: dpth2pres

!  input parameters:
!  nd     - size of arrays
!  depth  - depth in meters. no units check is made

!  output parameters:
!  pressure - pressure in bars

!  local variables & parameters:
real(r8), parameter :: c1 = 1.0_r8

! -----------------------------------------------------------------------
!  convert depth in meters to pressure in bars
! -----------------------------------------------------------------------

dpth2pres = 0.059808_r8 * (exp(-0.025_r8*depth) - c1)  &
              + 0.100766_r8*depth + 2.28405e-7_r8*depth**2

if (debug > 2 .and. do_output()) then
   print *, 'depth->pressure conversion table.  cols are: depth(m),pressure(bars)'
   print *, depth, dpth2pres
endif

end function dpth2pres



!-----------------------------------------------------------------------
!>
!> Get the in-situ temperature from the potential temperature, salinity
!> and pressure
!>

subroutine insitu_temp(potemp, s, lpres, insitu_t)

real(r8), intent(in)  :: potemp, s, lpres
real(r8), intent(out) :: insitu_t

! CODE FROM POP MODEL -
! nsc 1 nov 2012:  i have taken the original subroutine with call:
!  subroutine dpotmp(press,temp,s,rp,potemp)
! and removed the original 'press' argument (setting it to 0.0 below)
! and renamed temp -> potemp, and potemp -> insitu_t
! i also reordered the args to be a bit more logical.  now you specify:
! potential temp, salinity, local pressure in decibars, and you get
! back in-situ temperature (called sensible temperature in the atmosphere;
! what a thermometer would measure).  the original (F77 fixed format) code
! had a computed goto which is deprecated/obsolete.  i replaced it with
! a set of 'if() then else if()' lines.  i did try to not alter the original
! code so much it wasn't recognizable anymore.
!
!  aliciak note: rp = 0 and press = local pressure as function of depth
!  will return potemp given temp.
!  the trick here that if you make rp = local pressure and press = 0.0,
!  and put potemp in the "temp" variable , it will return insitu temp in the
!  potemp variable.

! an example figure of the relationship of potential temp and in-situ temp
! at depth:
! http://oceanworld.tamu.edu/resources/ocng_textbook/chapter06/chapter06_05.htm
! see the 'potential temperature' section (note graph starts at -1000m)

!     title:
!     *****

!       insitu_temp  -- calculate sensible (in-situ) temperature from
!                       local pressure, salinity, and potential temperature

!     purpose:
!     *******

!       to calculate sensible temperature, taken from a converter that
!       went from sensible/insitu temperature to potential temperature
!
!       ref: N.P. Fofonoff
!            Deep Sea Research
!            in press Nov 1976

!     arguments:
!     **********

!       potemp     -> potential temperature in celsius degrees
!       s          -> salinity pss 78
!       lpres      -> local pressure in decibars
!       insitu_t   <- in-situ (sensible) temperature (deg c)


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
            r4 = (r1 + (r2 + r3 + 1.8741e-8_r8)) * p + &
                 (-4.2393e-8_r8 * t+1.8932e-6_r8) * s1
            r5 = r4 + ((6.6228e-10_r8 * t-6.836e-8_r8) * t + 8.5258e-6_r8) * t + &
                 3.5803e-5_r8

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

if (debug > 2) print *, 'potential temp, salinity, local pressure -> sensible temp'
if (debug > 2) print *, potemp, s, lpres, insitu_t

!       potemp     -> potential temperature in celsius degrees
!       s          -> salinity pss 78
!       lpres      -> local pressure in decibars
!       insitu_t   <- in-situ (sensible) temperature (deg c)

end subroutine insitu_temp


!-----------------------------------------------------------------------
!>
!> Returns the the time step of the model; the smallest increment in
!> time that the model is capable of advancing the ROMS state.
!>

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = model_timestep

end function shortest_time_between_assimilations


!-----------------------------------------------------------------------
!>
!> Called to do one time initialization of the model.
!> In this case, it reads in the grid information, the namelist
!> containing the variables of interest, where to get them, their size,
!> their associated DART KIND, etc.
!>
!> In addition to harvesting the model metadata (grid,
!> desired model advance step, etc.), it also fills a structure
!> containing information about what variables are where in the DART
!> framework.

subroutine static_init_model()

integer :: iunit, io
integer :: ss, dd
integer :: ncid

character(len=32) :: calendar

type(time_type) :: model_time

if ( module_initialized ) return

! The Plan:
!
! * read in the grid sizes from grid file
! * allocate space, and read in actual grid values
! * figure out model timestep
! * Compute the model size.
! * set the index numbers where the field types change

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model:',string1,source,revision,revdate)

call nc_check( nf90_open(trim(roms_filename), NF90_NOWRITE, ncid), &
                  'static_init_model', 'open '//trim(roms_filename))

call get_time_information(roms_filename, ncid, 'ocean_time', 'ocean_time', &
                          calendar=calendar, last_time=model_time)

call set_calendar_type( trim(calendar) )

! Get the ROMS grid -- sizes and variables.
call get_grid_dimensions()
call get_grid()

! initialize the quad interp code
call init_quad_interp(grid_type = GRID_QUAD_FULLY_IRREGULAR, &
                      num_lons = Nxi_u, num_lats = Neta_u, &
                      cell_relative = QUAD_LOCATED_LON_EDGES, &
                      global = .false., spans_lon_zero = .false., &
                      pole_wrap = .false., interp_handle = ugrid_handle)

call init_quad_interp(grid_type = GRID_QUAD_FULLY_IRREGULAR, &
                      num_lons = Nxi_v, num_lats = Neta_v, &
                      cell_relative = QUAD_LOCATED_LAT_EDGES, &
                      global = .false., spans_lon_zero = .false., &
                      pole_wrap = .false., interp_handle = vgrid_handle)

call init_quad_interp(grid_type = GRID_QUAD_FULLY_IRREGULAR, &
                      num_lons = Nxi_rho, num_lats = Neta_rho, &
                      cell_relative = QUAD_LOCATED_CELL_CENTERS, &
                      global = .false., spans_lon_zero = .false., &
                      pole_wrap = .false., interp_handle = tgrid_handle)

call set_quad_coords(ugrid_handle, ULON, ULAT, UMASK)
call set_quad_coords(vgrid_handle, VLON, VLAT, VMASK)
call set_quad_coords(tgrid_handle, TLON, TLAT, TMASK)

! parse_variable_input() fills var_names, kind_list, clamp_vals, update_list
call parse_variable_input(variables, nfields)

domain_id = add_domain(roms_filename, nfields, &
                    var_names, kind_list, clamp_vals, update_list )

if (debug > 2) call state_structure_info(domain_id)

call nc_check( nf90_close(ncid), &
                  'static_init_model', 'close '//trim(roms_filename))

model_size = get_domain_size(domain_id)

call write_roms_time_information(model_time)

end subroutine static_init_model


!-----------------------------------------------------------------------
!>
!> Does any shutdown and clean-up needed for model.
!>

subroutine end_model()

! good style ... perhaps you could deallocate stuff (from static_init_model?).
! deallocate(state_loc)

call finalize_quad_interp(ugrid_handle)
call finalize_quad_interp(vgrid_handle)
call finalize_quad_interp(tgrid_handle)

if (allocated(ULAT))  deallocate(ULAT)
if (allocated(ULON))  deallocate(ULON)
if (allocated(UMASK)) deallocate(UMASK)

if (allocated(VLAT))  deallocate(VLAT)
if (allocated(VLON))  deallocate(VLON)
if (allocated(VMASK)) deallocate(VMASK)

if (allocated(TLAT))  deallocate(TLAT)
if (allocated(TLON))  deallocate(TLON)
if (allocated(TMASK)) deallocate(TMASK)

if (allocated(BATHY)) deallocate(BATHY)

end subroutine end_model


!-----------------------------------------------------------------------
!>
!> Writes the model-specific attributes to a DART 'diagnostic' netCDF file.
!> This includes coordinate variables and some metadata, but NOT the
!> actual DART state.
!>
!> @param ncid the netCDF handle of the DART diagnostic file opened by
!>                 assim_model_mod:init_diag_output
!> @param model_writes_state have the state structure write out all of the
!>                 state variables

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

! for the dimensions and coordinate variables
integer :: nxirhoDimID, nxiuDimID, nxivDimID
integer :: netarhoDimID, netauDimID, netavDimID
integer :: nsrhoDimID, nswDimID
integer :: VarID

! local variables

character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.

write(filename,*) 'ncid', ncid

! Write Global Attributes

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model", "ROMS")

! We need to output the grid information
! Define the new dimensions IDs

call nc_check(nf90_def_dim(ncid, name='xi_rho',  len = Nxi_rho, &
     dimid = nxirhoDimID),'nc_write_model_atts', 'xi_rho def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='eta_rho', len = Neta_rho,&
     dimid = netarhoDimID),'nc_write_model_atts', 'eta_rho def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='s_rho',   len = Ns_rho,&
     dimid = nsrhoDimID),'nc_write_model_atts', 's_rho def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='s_w',   len = Ns_w,&
     dimid = nswDimID),'nc_write_model_atts', 's_w def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='xi_u',    len = Nxi_u,&
     dimid = nxiuDimID),'nc_write_model_atts', 'xi_u def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='xi_v',    len = Nxi_v,&
     dimid = nxivDimID),'nc_write_model_atts', 'xi_v def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='eta_u',   len = Neta_u,&
     dimid = netauDimID),'nc_write_model_atts', 'eta_u def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='eta_v',   len = Neta_v,&
     dimid = netavDimID),'nc_write_model_atts', 'eta_v def_dim '//trim(filename))

! Create the Coordinate Variables and give them Attributes
! The values will be added in a later block of code.

call nc_check(nf90_def_var(ncid,name='lon_rho', xtype=nf90_double, &
              dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
              'nc_write_model_atts', 'lon_rho def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'rho longitudes'), &
              'nc_write_model_atts', 'lon_rho long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_east'), &
              'nc_write_model_atts', 'lon_rho units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='lat_rho', xtype=nf90_double, &
              dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
              'nc_write_model_atts', 'lat_rho def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'rho latitudes'), &
              'nc_write_model_atts', 'lat_rho long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_north'), &
              'nc_write_model_atts', 'lat_rho units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='lon_u', xtype=nf90_double, &
              dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
              'nc_write_model_atts', 'lon_u def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'u longitudes'), &
              'nc_write_model_atts', 'lon_u long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_east'), &
              'nc_write_model_atts', 'lon_u units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='lat_u', xtype=nf90_double, &
              dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
              'nc_write_model_atts', 'lat_u def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'u latitudes'), &
              'nc_write_model_atts', 'lat_u long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_north'), &
              'nc_write_model_atts', 'lat_u units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='lon_v', xtype=nf90_double, &
              dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
              'nc_write_model_atts', 'lon_v def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'v longitudes'), &
              'nc_write_model_atts', 'lon_v long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_east'), &
              'nc_write_model_atts', 'lon_v units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='lat_v', xtype=nf90_double, &
              dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
              'nc_write_model_atts', 'lat_v def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'v latitudes'), &
              'nc_write_model_atts', 'lat_v long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_north'), &
              'nc_write_model_atts', 'lat_v units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='z_rho', xtype=nf90_double, &
              dimids=(/ nxirhoDimID, netarhoDimID, nsrhoDimID /), varid=VarID),&
              'nc_write_model_atts', 'z_rho def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
              'nc_write_model_atts', 'z_rho long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
              'nc_write_model_atts', 'z_rho units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='z_u', xtype=nf90_double, &
              dimids=(/ nxiuDimID, netauDimID, nsrhoDimID /), varid=VarID),&
              'nc_write_model_atts', 'z_u def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
              'nc_write_model_atts', 'z_u long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
              'nc_write_model_atts', 'z_u units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='z_v', xtype=nf90_double, &
              dimids=(/ nxivDimID, netavDimID, nsrhoDimID /), varid=VarID),&
              'nc_write_model_atts', 'z_v def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
              'nc_write_model_atts', 'z_v long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
              'nc_write_model_atts', 'z_v units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='z_w', xtype=nf90_double, &
              dimids=(/ nxirhoDimID, netarhoDimID, nswDimID /), varid=VarID),&
              'nc_write_model_atts', 'z_w def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
              'nc_write_model_atts', 'z_w long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
              'nc_write_model_atts', 'z_w units '//trim(filename))

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call nc_end_define_mode(ncid)

! Fill the coordinate variable values
! the RHO grid

call nc_check(NF90_inq_varid(ncid, 'lon_rho', VarID), &
              'nc_write_model_atts', 'lon_rho inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, TLON ), &
             'nc_write_model_atts', 'lon_rho put_var '//trim(filename))

call nc_check(NF90_inq_varid(ncid, 'lat_rho', VarID), &
              'nc_write_model_atts', 'lat_rho inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, TLAT ), &
             'nc_write_model_atts', 'lat_rho put_var '//trim(filename))

! the U grid

call nc_check(NF90_inq_varid(ncid, 'lon_u', VarID), &
              'nc_write_model_atts', 'lon_u inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, ULON ), &
             'nc_write_model_atts', 'lon_u put_var '//trim(filename))

call nc_check(NF90_inq_varid(ncid, 'lat_u', VarID), &
              'nc_write_model_atts', 'lat_u inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, ULAT ), &
             'nc_write_model_atts', 'lat_u put_var '//trim(filename))

! the V grid

call nc_check(NF90_inq_varid(ncid, 'lon_v', VarID), &
              'nc_write_model_atts', 'lon_v inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, VLON ), &
             'nc_write_model_atts', 'lon_v put_var '//trim(filename))

call nc_check(NF90_inq_varid(ncid, 'lat_v', VarID), &
              'nc_write_model_atts', 'lat_v inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, VLAT ), &
             'nc_write_model_atts', 'lat_v put_var '//trim(filename))

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)


end subroutine nc_write_model_atts

!-----------------------------------------------------------------------
!> writes the time of the current state and (optionally) the time
!> to be conveyed to ROMS to dictate the length of the forecast.
!> This file is then used by scripts to modify the ROMS run.
!> The format in the time information is totally at your discretion.
!>
!> @param ncfile_out name of the file
!> @param model_time the current time of the model state
!> @param adv_to_time the time in the future of the next assimilation.
!>

subroutine write_model_time(ncid, model_time, adv_to_time)
integer,         intent(in)           :: ncid
type(time_type), intent(in)           :: model_time
type(time_type), intent(in), optional :: adv_to_time

integer :: io, varid, seconds, days
type(time_type) :: origin_time, deltatime
real(digits12)  :: run_duration

if ( .not. module_initialized ) call static_init_model

if (present(adv_to_time)) then
   string3 = time_to_string(adv_to_time)
   write(string1,*)'ROMS/DART not configured to advance ROMS.'
   write(string2,*)'called with optional advance_to_time of'
   call error_handler(E_ERR, 'write_model_time', string1, &
              source, revision, revdate, text2=string2,text3=string3)
endif

! If the ocean_time variable exists, we are updating a ROMS file,
! if not ... must be updating a DART diagnostic file.

io = nf90_inq_varid(ncid,'ocean_time',varid)
if (io == NF90_NOERR) then
   call get_time_information('unknown', ncid, 'ocean_time', 'ocean_time', &
                myvarid=varid, origin_time=origin_time)
   deltatime = model_time - origin_time
   call get_time(deltatime, seconds, days)
   run_duration = real(days,digits12)*86400.0_digits12 + real(seconds,digits12)
   call nc_check(nf90_put_var(ncid, varid, run_duration), 'write_model_time', 'put_var')
   return
endif

io = nf90_inq_varid(ncid,'time',varid)
if (io == NF90_NOERR) then
   call get_time_information('unknown', ncid, 'time', 'time', &
                myvarid=varid, origin_time=origin_time)
   deltatime = model_time - origin_time
   call get_time(deltatime, seconds, days)
   run_duration = real(days,digits12)*86400.0_digits12 + real(seconds,digits12)
   call nc_check(nf90_put_var(ncid, varid, run_duration), 'write_model_time', 'put_var')
   return
endif

end subroutine write_model_time

!--------------------------------------------------------------------
!>
!> read the time from the input file
!>
!> @param filename name of file that contains the time
!>

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type)              :: read_model_time

integer :: ncid

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'read_model_time',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'read_model_time', 'open '//trim(filename))

call get_time_information(filename, ncid, 'ocean_time', 'ocean_time', &
                          last_time=read_model_time)

call nc_check( nf90_close(ncid), 'read_model_time', 'close '//trim(filename))

end function read_model_time



!-----------------------------------------------------------------------
! The remaining (private) interfaces come last.
! None of the private interfaces need to call static_init_model()
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>
!> Set the desired minimum model advance time. This is generally NOT the
!> dynamical timestep of the model, but rather the shortest forecast length
!> you are willing to make. This impacts how frequently the observations
!> may be assimilated.
!>

function set_model_time_step()

type(time_type) :: set_model_time_step

! assimilation_period_seconds, assimilation_period_days are from the namelist

!>@todo FIXME make sure set_model_time_step is an integer multiple of
!> the dynamical timestep or whatever strategy ROMS employs.
! TJH NHIST*DT ... from the inputfile ... can remove assim_* from DART input namelist

!>@todo FIXME : JPH we should really be getting this from the history file??
set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!-----------------------------------------------------------------------
!>
!> Read the grid dimensions from the ROMS grid netcdf file.
!> By reading the dimensions first, we can use them in variable
!> declarations later - which is faster than using allocatable arrays.
!>

subroutine get_grid_dimensions()

integer :: ncid

! Read the (static) grid dimensions from the ROMS grid file.

call nc_check(nf90_open(trim(roms_grid_file), nf90_nowrite, ncid), &
              'get_grid_dimensions', 'open '//trim(roms_grid_file))

Nxi_rho   = get_dimension_length(ncid, 'xi_rho',   roms_grid_file)
Nxi_u     = get_dimension_length(ncid, 'xi_u',     roms_grid_file)
Nxi_v     = get_dimension_length(ncid, 'xi_v',     roms_grid_file)
Neta_rho  = get_dimension_length(ncid, 'eta_rho',  roms_grid_file)
Neta_u    = get_dimension_length(ncid, 'eta_u',    roms_grid_file)
Neta_v    = get_dimension_length(ncid, 'eta_v',    roms_grid_file)

call nc_check(nf90_close(ncid), &
              'get_grid_dimensions','close '//trim(roms_grid_file))

! Read the vertical dimensions from the dedicated file.

call nc_check(nf90_open(trim(roms_grid_file), nf90_nowrite, ncid), &
               'get_grid_dimensions', 'open '//trim(roms_grid_file))

Ns_rho    = get_dimension_length(ncid, 's_rho',    roms_grid_file)
Ns_w      = get_dimension_length(ncid, 's_w'  ,    roms_grid_file)

call nc_check(nf90_close(ncid), &
              'get_grid_dimensions','close '//trim(roms_grid_file))

Nx =  Nxi_rho  ! Setting the nominal value of the 'global' variables
Ny = Neta_rho  ! Setting the nominal value of the 'global' variables
Nz =   Ns_rho  ! Setting the nominal value of the 'global' variables

end subroutine get_grid_dimensions


!-----------------------------------------------------------------------
!>
!> Read the actual grid values from the ROMS netcdf file.
!>
!>@todo FIXME:  the original implementation opened 3 different files
!> to get the grid info - the namelist was:
!>    roms_ini_filename            = '../data/wc13_ini.nc'
!>    grid_definition_filename     = '../data/wc13_grd.nc'
!>    depths_definition_filename   = '../data/wc13_depths.nc'
!>
!> these have been consolidated by hernan for the santa cruz version
!> into a single file.  check with the other rutgers folks to see if
!> they still need to open 3 different files.  if so, we might need
!> to restore the 3 namelist items and we can use the same file for
!> all 3 types of grid info in the first case, and 3 different files
!> for the second case.
!>

subroutine get_grid()

integer  :: ncid, VarID

real(r8), parameter :: all_land = 0.001_r8

real(r8), allocatable :: mask(:,:)

if (.not. allocated(ULAT)) allocate(ULAT(Nxi_u, Neta_u))
if (.not. allocated(ULON)) allocate(ULON(Nxi_u, Neta_u))
if (.not. allocated(UMASK)) allocate(UMASK(Nxi_u, Neta_u))

if (.not. allocated(VLAT)) allocate(VLAT(Nxi_v, Neta_v))
if (.not. allocated(VLON)) allocate(VLON(Nxi_v, Neta_v))
if (.not. allocated(VMASK)) allocate(VMASK(Nxi_v, Neta_v))

if (.not. allocated(TLAT)) allocate(TLAT(Nxi_rho, Neta_rho))
if (.not. allocated(TLON)) allocate(TLON(Nxi_rho, Neta_rho))
if (.not. allocated(TMASK)) allocate(TMASK(Nxi_rho, Neta_rho))

if (.not. allocated(BATHY)) allocate(BATHY(Nxi_rho, Neta_rho))

! Assume everything is ocean ... is_masked() requires land = .true..
UMASK = .false.
VMASK = .false.
TMASK = .false.

call nc_check(nf90_open(trim(roms_grid_file), nf90_nowrite, ncid), &
      'get_grid', 'open '//trim(roms_grid_file))

! Read the vertical grid parameters

call nc_check(nf90_inq_varid(ncid, 'theta_s', VarID), &
   'get_grid', 'inq_varid theta_s '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, theta_s), &
      'get_grid', 'get_var theta_s '//trim(roms_grid_file))

call nc_check(nf90_inq_varid(ncid, 'theta_b', VarID), &
   'get_grid', 'inq_varid theta_b '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, theta_b), &
      'get_grid', 'get_var theta_b '//trim(roms_grid_file))

call nc_check(nf90_inq_varid(ncid, 'Tcline', VarID), &
   'get_grid', 'inq_varid Tcline '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, Tcline), &
      'get_grid', 'get_var Tcline '//trim(roms_grid_file))

call nc_check(nf90_inq_varid(ncid, 'hc', VarID), &
   'get_grid', 'inq_varid hc '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, hc), &
      'get_grid', 'get_var hc '//trim(roms_grid_file))

! Read the rest of the grid information from the traditional grid file

call nc_check(nf90_inq_varid(ncid, 'lon_rho', VarID), &
   'get_grid', 'inq_varid lon_rho '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, TLON), &
      'get_grid', 'get_var lon_rho '//trim(roms_grid_file))

where (TLON < 0.0_r8) TLON = TLON + 360.0_r8

call nc_check(nf90_inq_varid(ncid, 'lat_rho', VarID), &
      'get_grid', 'inq_varid lat_rho '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, TLAT), &
      'get_grid', 'get_var lat_rho '//trim(roms_grid_file))

call nc_check(nf90_inq_varid(ncid, 'lon_u', VarID), &
      'get_grid', 'inq_varid lon_u '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, ULON), &
      'get_grid', 'get_var lon_u '//trim(roms_grid_file))

where (ULON < 0.0_r8) ULON = ULON + 360.0_r8

call nc_check(nf90_inq_varid(ncid, 'lat_u', VarID), &
      'get_grid', 'inq_varid lat_u '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, ULAT), &
      'get_grid', 'get_var lat_u '//trim(roms_grid_file))

call nc_check(nf90_inq_varid(ncid, 'lon_v', VarID), &
      'get_grid', 'inq_varid lon_v '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, VLON), &
      'get_grid', 'get_var lon_v '//trim(roms_grid_file))

where (VLON < 0.0_r8) VLON = VLON + 360.0_r8

call nc_check(nf90_inq_varid(ncid, 'lat_v', VarID), &
      'get_grid', 'inq_varid lat_v '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, VLAT), &
      'get_grid', 'get_var lat_v '//trim(roms_grid_file))

call nc_check(nf90_inq_varid(ncid, 'h', VarID), &
      'get_grid', 'inq_varid h '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, BATHY), &
      'get_grid', 'get_var h '//trim(roms_grid_file))

allocate(mask(Nxi_u, Neta_u))
call nc_check(nf90_inq_varid(ncid, 'mask_u', VarID), &
      'get_grid', 'inq_varid mask_u '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, mask), &
      'get_grid', 'get_var mask_u '//trim(roms_grid_file))
where(mask < 0.5_r8) UMASK = .true.  ! 1.0 is water, 0.0 is land

deallocate(mask)
allocate(mask(Nxi_v, Neta_v))
call nc_check(nf90_inq_varid(ncid, 'mask_v', VarID), &
      'get_grid', 'inq_varid mask_v '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, mask), &
      'get_grid', 'get_var mask_v '//trim(roms_grid_file))
where(mask < 0.5_r8) VMASK = .true.  ! 1.0 is water, 0.0 is land

deallocate(mask)
allocate(mask(Nxi_rho, Neta_rho))
call nc_check(nf90_inq_varid(ncid, 'mask_rho', VarID), &
      'get_grid', 'inq_varid mask_rho '//trim(roms_grid_file))
call nc_check(nf90_get_var( ncid, VarID, mask), &
      'get_grid', 'get_var mask_rho '//trim(roms_grid_file))
where(mask < 0.5_r8) TMASK = .true.  ! 1.0 is water, 0.0 is land

deallocate(mask)

! Be aware that all the depths are negative values.
! The surface of the ocean is 0.0, the deepest is a big negative value.

if (do_output() .and. debug > 0) then
    write(string1,*)'    min/max ULON ',minval(ULON), maxval(ULON)
    write(string2,*)    'min/max ULAT ',minval(ULAT), maxval(ULAT)
    call error_handler(E_MSG,'get_grid',string1, text2=string2)

    write(string1,*)'    min/max VLON ',minval(VLON), maxval(VLON)
    write(string2,*)    'min/max VLAT ',minval(VLAT), maxval(VLAT)
    call error_handler(E_MSG,'get_grid',string1, text2=string2)

    write(string1,*)'    min/max TLON ',minval(TLON), maxval(TLON)
    write(string2,*)    'min/max TLAT ',minval(TLAT), maxval(TLAT)
    call error_handler(E_MSG,'get_grid',string1, text2=string2)

    write(string1,*)    'min/max BATHY ',minval(BATHY), maxval(BATHY)
    call error_handler(E_MSG,'get_grid',string1)

endif

end subroutine get_grid


!-----------------------------------------------------------------------
!>
!> Fill the array of requested variables, dart kinds, possible min/max
!> values and whether or not to update the field in the output file.
!>
!>@param state_variables the list of variables and kinds from model_mod_nml
!>@param ngood the number of variable/KIND pairs specified

subroutine parse_variable_input( state_variables, ngood )

character(len=*), intent(in)  :: state_variables(:)
integer,          intent(out) :: ngood

character(len=*), parameter :: routine = 'parse_variable_input'
integer :: i
character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 5   change to updateable

ngood = 0
MyLoop : do i = 1, MAX_STATE_VARIABLES

   varname      = trim(state_variables(num_state_table_columns*i-4))
   dartstr      = trim(state_variables(num_state_table_columns*i-3))
   minvalstring = trim(state_variables(num_state_table_columns*i-2))
   maxvalstring = trim(state_variables(num_state_table_columns*i-1))
   state_or_aux = trim(state_variables(num_state_table_columns*i  ))

   if ( varname == ' ' .and. dartstr == ' ' ) exit MyLoop ! Found end of list.

   if ( varname == ' ' .or. dartstr == ' ' ) then
      string1 = 'model_nml:model "variables" not fully specified'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(''there is no quantity <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

   call to_upper(minvalstring)
   call to_upper(maxvalstring)
   call to_upper(state_or_aux)

   var_names(   i) = varname
   kind_list(   i) = get_index_for_quantity(dartstr)
   clamp_vals(i,1) = string_to_real(minvalstring)
   clamp_vals(i,2) = string_to_real(maxvalstring)
   update_list( i) = string_to_logical(state_or_aux, 'UPDATE')

   ngood = ngood + 1

enddo MyLoop

if (ngood == MAX_STATE_VARIABLES) then
   string1 = 'WARNING: There is a possibility you need to increase ''MAX_STATE_VARIABLES'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)
endif

end subroutine parse_variable_input


!-----------------------------------------------------------------------
!
!> Find the named variable (often 'ocean_time') in a ROMS netCDF file.
!> If it is not found, it is a fatal error.
!>
!> @param filename the name of the ROMS netCDF file
!>                 (used to generate useful error messages).
!> @param ncid the netCDF handle to the ROMS netCDF file.
!> @param variable name which contains the time
!> @param calendar the character string indicating the calendar in use
!> @param last_time_index the value of the last time dimension
!> @param last_time the time/date of the last time
!> @param origin_time the base time other times are relative to
!> @param all_times an array of all times in the variable
!>
!>@todo FIXME Make sure the calculation is correct.
!>  A 64bit real can support whole numbers that overflow a 32 bit integer.

subroutine get_time_information(filename, ncid, var_name, dim_name, myvarid, &
                    calendar, last_time_index, last_time, origin_time, all_times)

character(len=*),            intent(in)  :: filename
integer,                     intent(in)  :: ncid
character(len=*),            intent(in)  :: var_name
character(len=*),            intent(in)  :: dim_name
integer,           optional, intent(out) :: myvarid
character(len=32), optional, intent(out) :: calendar
integer,           optional, intent(out) :: last_time_index
type(time_type),   optional, intent(out) :: last_time
type(time_type),   optional, intent(out) :: origin_time
type(time_type),   optional, intent(out) :: all_times(:)


character(len=*), parameter :: routine = 'get_time_information'
integer :: ios, DimID, VarID, dimlen, i
character(len=64) :: unitstring
character(len=32) :: calendarstring

integer :: year, month, day, hour, minute, second, rc
real(digits12), allocatable :: these_times(:)
type(time_type) :: time_offset, base_time

logical :: offset_in_seconds  ! if .false., assuming offset in days

integer :: original_calendar_type

!>@todo FIXME get the variable length from the varid and remove the need for the dimension name

call nc_check(nf90_inq_dimid(ncid,dim_name,dimid=DimID), &
       routine,'cannot find "'//trim(dim_name)//'" dimension in '//trim(filename))

call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
       routine, 'inquire_dimension '//trim(dim_name)//' from '//trim(filename))
if (present(last_time_index)) last_time_index = dimlen

call nc_check(nf90_inq_varid(ncid, var_name, VarID), &
       routine, 'inq_varid '//trim(var_name)//' from '//trim(filename))
if (present(myvarid)) myvarid = VarID

! assume gregorian calendar unless there's a calendar attribute saying elsewise
rc = nf90_get_att(ncid, VarID, 'calendar', calendarstring)
if (rc /= nf90_noerr) calendarstring = 'gregorian'

if (index(calendarstring,'gregorian') == 0) then
   write(string1,*)'expecting '//trim(var_name)//' calendar of "gregorian"'
   write(string2,*)'got '//trim(calendarstring)
   write(string3,*)'from file "'//trim(filename)//'"'
   call error_handler(E_MSG,routine, string1, &
             source, revision, revdate, text2=string2, text3=string3)
else
   ! coerce all forms of gregorian to the one DART supports
   ! 'gregorian_proleptic' needs to be changed, for example.
   calendarstring = 'gregorian'
endif

if (present(calendar)) calendar = trim(calendarstring)

if (present(last_time) .or. present(origin_time) .or. present(all_times)) then

   ! May need to put the calendar back to some original value
   original_calendar_type = get_calendar_type()

   ! We need to set the calendar to interpret the time values
   ! do we need to preserve the original calendar setting if there is one?

   call set_calendar_type( trim(calendarstring) )

   ! Make sure the calendar is expected form
   ! var_name:units    = "seconds since 1999-01-01 00:00:00" ;
   !                      1234567890123
   !   OR
   ! var_name:units    = "days since 1999-01-01 00:00:00" ;
   !                      1234567890

   call nc_check(nf90_get_att(ncid, VarID, 'units', unitstring), &
          routine, 'get_att '//trim(var_name)//' units '//trim(filename))

   ! decode the start time of the time variable - expecting time to be coded
   ! as an offset to some base

   if (unitstring(1:13) == 'seconds since') then
      read(unitstring,'(14x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
      if (ios /= 0) then
         write(string1,*)'Unable to read time variable units. Error status was ',ios
         write(string2,*)'expected "seconds since YYYY-MM-DD HH:MM:SS"'
         write(string3,*)'was      "'//trim(unitstring)//'"'
         call error_handler(E_ERR, routine, string1, &
                source, revision, revdate, text2=string2, text3=string3)
      endif
      offset_in_seconds = .true.

   else if (unitstring(1:10) == 'days since') then
      read(unitstring,'(11x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
      if (ios /= 0) then
         write(string1,*)'Unable to read time variable units. Error status was ',ios
         write(string2,*)'expected "days since YYYY-MM-DD HH:MM:SS"'
         write(string3,*)'was      "'//trim(unitstring)//'"'
         call error_handler(E_ERR, routine, string1, &
                source, revision, revdate, text2=string2, text3=string3)
      endif
      offset_in_seconds = .false.

   else
      write(string1,*)'expecting time attribute units of "seconds since ..." -OR-'
      write(string2,*)'                              "days since ..."'
      write(string3,*)'got "'//trim(unitstring)//'"'
      call error_handler(E_ERR,routine, string1, &
                source, revision, revdate, text2=string2, text3=string3)
   endif

   base_time = set_date(year, month, day, hour, minute, second)

   if (present(origin_time)) origin_time = base_time

   if (present(last_time) .or. present(all_times)) then

      ! big_integer may overflow a 32bit integer, so declare it 64bit
      ! and parse it into an integer number of days and seconds, both
      ! of which can be 32bit. Our set_time, set_date routines need 32bit integers.

      allocate(these_times(dimlen))

      call nc_check(nf90_get_var( ncid, VarID, these_times), &
             routine, 'get_var '//trim(var_name)//' from '//trim(filename))

      if (present(last_time)) then
         time_offset = convert_to_time_offset(these_times(dimlen), offset_in_seconds)
         last_time = base_time + time_offset
      endif

      if (present(all_times)) then
         do i=1, dimlen
            time_offset = convert_to_time_offset(these_times(i), offset_in_seconds)
            all_times(i) = base_time + time_offset
         enddo
      endif

      if (do_output() .and. debug > 0 .and. present(last_time)) then
         call print_time(last_time, str='last roms time is ',iunit=logfileunit)
         call print_time(last_time, str='last roms time is ')
         call print_date(last_time, str='last roms date is ',iunit=logfileunit)
         call print_date(last_time, str='last roms date is ')
      endif

      deallocate(these_times)

   endif

   call set_calendar_type(original_calendar_type)

endif

end subroutine get_time_information

!-----------------------------------------------------------------------

!> convert a fractional day to a dart time type

function convert_to_time_offset(offset, offset_in_seconds)

real(digits12), intent(in) :: offset
logical,        intent(in) :: offset_in_seconds
type(time_type) :: convert_to_time_offset

integer(i8) :: big_integer
integer :: some_seconds, some_days

if (offset_in_seconds) then
   big_integer  = int(offset,i8)
   some_days    = big_integer / (24*60*60)
   big_integer  = big_integer - (some_days * (24*60*60))
   some_seconds = int(big_integer,i4)
else
   ! offset in fractional days
   some_days    = int(offset)
   some_seconds = (offset - some_days) * (24*60*60)
endif

convert_to_time_offset = set_time(some_seconds, some_days)

end function convert_to_time_offset

!-----------------------------------------------------------------------
!>
!> convert DART time type into a character string with the
!> format of YYYYMMDDhh ... or DDhh
!>
!> @param time_to_string the character string containing the time
!> @param t the time
!> @param interval logical flag describing if the time is to be
!>                 interpreted as a calendar date or a time increment.
!>                 If the flag is merely present, the time is to be
!>                 interpreted as an increment and the format is simply
!>                 DDhh. If the flag is not present, the time is a full
!>                 calendar (Gregorian) date and will be renedered with
!>                 the YYYYMMDDhh format.
!>

function time_to_string(t, interval)

character(len=19)              :: time_to_string
type(time_type),   intent(in) :: t
logical, optional, intent(in) :: interval

! local variables

integer :: iyear, imonth, iday, ihour, imin, isec
integer :: ndays, nsecs
logical :: dointerval

if (present(interval)) then
   dointerval = interval
else
   dointerval = .false.
endif

! for interval output, output the number of days, then hours, mins, secs
! for date output, use the calendar routine to get the year/month/day hour:min:sec
if (dointerval) then
   call get_time(t, nsecs, ndays)
   if (ndays > 99) then
      write(string1, *) 'interval number of days is ', ndays
      call error_handler(E_ERR,'time_to_string:', 'interval days cannot be > 99', &
                         source, revision, revdate, text2=string1)
   endif
   ihour = nsecs / 3600
   nsecs = nsecs - (ihour * 3600)
   imin  = nsecs / 60
   nsecs = nsecs - (imin * 60)
   isec  = nsecs
!   write(time_to_string, '(I2.2,3(A1,I2.2))') &
!                        ndays, '_', ihour, ':', imin, ':', isec
   write(time_to_string, '(I2.2,I2.2)') &
                        ndays, ihour
else
   call get_date(t, iyear, imonth, iday, ihour, imin, isec)
   write(time_to_string, '(I4.4,5(A1,I2.2))') &
          iyear, '-', imonth, '-', iday, ' ', ihour, ':', imin, ':', isec
endif

end function time_to_string


!-----------------------------------------------------------------------
!>
!> gets the length of a netCDF dimension given the dimension name.
!> This bundles the nf90_inq_dimid and nf90_inquire_dimension routines
!> into a slightly easier-to-use function.
!>
!> @param dimlen the length of the netCDF dimension in question
!> @param ncid the netCDF file handle
!> @param dimension_name the character string of the dimension name
!> @param filename the name of the netCDF file (for error message purposes)
!>

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

!-----------------------------------------------------------------------
!>
!> writes the time of the current state and (optionally) the time
!> to be conveyed to ROMS to dictate the length of the forecast.
!> This file is then used by scripts to modify the ROMS run.
!> The format in the time information is totally at your discretion.
!>
!> @param model_time the current time of the model state
!>

subroutine write_roms_time_information( model_time )
type(time_type), intent(in) :: model_time

integer :: iunit, day, second, ios, ncid
real(digits12) :: dstart
type(time_type) :: base_time, forecast_time

base_time = model_time

call nc_check( nf90_open(trim(roms_filename), NF90_NOWRITE, ncid), &
                  'write_roms_time_information', 'open '//trim(roms_filename))

call get_time_information(roms_filename, ncid, 'ocean_time', 'ocean_time', &
                          origin_time=base_time)

call nc_check( nf90_close(ncid), &
                 'write_roms_time_information', 'close '//trim(roms_filename))

forecast_time = model_time - base_time

call get_time(forecast_time,second,day)

dstart = real(day,digits12) + second/86400.0_digits12

iunit = open_file('new_dstart.txt', action='write')

write(iunit,'(A,1x,f25.5)', iostat=ios) 'DSTART ', dstart
if (ios /= 0) then
   write(string1,*)'Unable to write new DSTART. Error status was ',ios
   write(string2,*)'dstart = ', dstart
   call error_handler(E_ERR, 'write_roms_time_information:', string1, &
          source, revision, revdate, text2=string2)
endif

! The next records are totally optional - not checking write status
write(iunit,'(A,1x,f25.5)') 'ROMS_offset ', dstart*86400.0_digits12
call print_date(model_time, str='ROMS_date ',iunit=iunit)
call print_time(model_time, str='DART_time ',iunit=iunit)

string1 = time_to_string(model_time)
write(iunit,'(A,1x,A)') 'YYYYMMDD',trim(string1)

call close_file(iunit)

end subroutine write_roms_time_information


!-----------------------------------------------------------------------
!>
!> Returns the lat,lon,depth given a fractional i,j,k and a specified kind
!>
!> @param filoc fractional x index
!> @param fjloc fractional y index
!> @param fkloc fractional vert index
!> @param dart_kind
!> @param locatiation location at fractional i,j,k
!>
!>  Each grid cell is oriented in a counter clockwise direction
!>  for interpolating locations.  First we interpolate in latitude
!>  and longitude, then interpolate in height.  The hgt of each grid
!>  cell can very on each interpolation, so we have to be careful how
!>  we interpolate in the horizontal.  Using the 4 different heights
!>  and lat_frac, lon_frac, hgt_frac we can do a simple trilinear
!>  interpolation to find the location given fractional indicies.
!>
!>              (i ,j+1) ----- (i+1,j+1)   hgt(1) hgt(2) hgt(3) hgt(4)
!>               hgt(4)         hgt(3)       |      |      |      |
!>                 |               |         |      *      |      |
!>  lon_frac _____ |       X       |         |      |      |      *
!>                 |               |         |             *      |
!>                 |               |         *             |
!>              (i ,j)   ----- (i+1,j)       |
!>               hgt(1)   |     hgt(2)             * - hgt_frac location
!>                        |
!>                     lat_frac
!>
!>    ISTATUS : 10 - bad incoming dart_kind
!>    ISTATUS : 11 - fkloc out of range
!>    ISTATUS : 12 - filoc or fjloc out of range for u grid
!>    ISTATUS : 13 - filoc or fjloc out of range for v grid
!>    ISTATUS : 14 - filoc or fjloc out of range for rho grid
!>    ISTATUS : 99 - initalized istatus, this should not happen

!function get_location_from_ijk(filoc, fjloc, fkloc, dart_kind, location) result(istatus)
!real(r8),            intent(in)  :: filoc
!real(r8),            intent(in)  :: fjloc
!real(r8),            intent(in)  :: fkloc
!integer,             intent(in)  :: dart_kind
!type(location_type), intent(out) :: location
!integer :: istatus ! 0 good, else bad
!
!integer  :: var_id, iloc, jloc, vloc
!integer  :: vert_type, my_kind
!real(r8) :: lon_fract, lat_fract, hgt_fract
!real(r8) :: lon_val, lat_val, hgt_val, tmp_hgt(2), hgt(4)
!real(r8), pointer :: mylon(:,:), mylat(:,:), mydep(:,:,:)
!logical, save :: first_time = .true.
!
!
!write(string1,*)'Routine not finished.'
!call error_handler(E_ERR, 'get_location_from_ijk:', string1, &
!                      source, revision, revdate)
!
!! start out assuming bad istatus
!istatus  = 99
!
!var_id = get_varid_from_kind(domain_id, dart_kind)
!if (var_id < 0) then
!  istatus = 10 ! can not find variable id for dart_kind
!  return
!endif
!
!! check that we have a valid vertical location.
!! allow obs above the top rho point but below the
!! surface to be considered at the top rho point.
!! the commented out code is the test for excluding
!! obs above the top rho point in addition to obs
!! below the bottom rho point.
!!if (fkloc < 0.5 .or. fkloc > Ns_rho - 0.5) then
!if (fkloc < 0.5 .or. fkloc > Ns_rho) then
!  istatus = 11
!  location = set_location_missing()
!  return
!endif
!
!iloc = FLOOR(filoc)
!jloc = FLOOR(fjloc)
!vloc = FLOOR(fkloc)
!
!lon_fract = filoc - iloc
!lat_fract = fjloc - jloc
!hgt_fract = fkloc - vloc
!
!my_kind = get_kind_index(domain_id,var_id)
!if (my_kind==QTY_U_CURRENT_COMPONENT) then
!   write(string1,*)'Not interpolating ', get_name_for_quantity(my_kind), ' at the moment.'
!   write(string2,*)'Need to check that we are using the right grid for location interpolation'
!   call error_handler(E_ERR, 'get_location_from_ijk:', string1, &
!                      source, revision, revdate, text2=string2)
!   if (filoc < 1 .or. filoc > Nxi_u-1 .or. &
!       fjloc < 1 .or. fjloc > Neta_u-1 ) then
!     istatus = 12
!     location = set_location_missing()
!     return
!   endif
!   mylon => ULON
!   mylat => ULAT
!elseif (my_kind==QTY_V_CURRENT_COMPONENT) then
!   write(string1,*)'Not interpolating ', get_name_for_quantity(my_kind), ' at the moment.'
!   write(string2,*)'Need to check that we are using the right grid for location interpolation'
!   call error_handler(E_ERR, 'get_location_from_ijk:', string1, &
!                      source, revision, revdate, text2=string2)
!   if (filoc < 1 .or. filoc > Nxi_v-1 .or. &
!       fjloc < 1 .or. fjloc > Neta_v-1 ) then
!     istatus = 13
!     location = set_location_missing()
!     return
!   endif
!   mylon => VLON
!   mylat => VLAT
!   mydep => VDEP
!else  ! Everything else is assumed to be on the rho points
!   if (filoc < 1 .or. filoc > Nxi_rho-1 .or. &
!       fjloc < 1 .or. fjloc > Neta_rho-1 ) then
!
!     write(*,*)
!     write(*,*)'filoc, Nxi_rho-1       = ',filoc, Nxi_rho-1
!     write(*,*)'fjloc, Neta_rho-1      = ',fjloc, Neta_rho-1
!     write(*,*)'fkloc, vloc, hgt_fract = ',fkloc,vloc,hgt_fract
!
!     write(logfileunit,*)
!     write(logfileunit,*)'filoc, Nxi_rho-1     = ',filoc, Nxi_rho-1
!     write(logfileunit,*)'fjloc, Neta_rho-1    = ',fjloc, Neta_rho-1
!     write(logfileunit,*)'fkloc,vloc,hgt_fract = ',fkloc,vloc,hgt_fract
!
!     istatus = 14
!     location = set_location_missing()
!     return
!   endif
!   mylon => TLON
!   mylat => TLAT
!   mydep => TDEP
!endif
!
!lon_val = (1.0-lon_fract)*(mylon(iloc,jloc)) + (lon_fract)*(mylon(iloc+1,jloc))
!lat_val = (1.0-lat_fract)*(mylat(iloc,jloc)) + (lat_fract)*(mylat(iloc  ,jloc+1))
!
!if( get_kind_index(domain_id,var_id) == QTY_SEA_SURFACE_HEIGHT ) then
!   hgt_val    = 0.0_r8
!   vert_type  = VERTISSURFACE
!else
!   if (fkloc > Ns_rho - 0.5) then
!      ! leave something in the log to remind us we're doing this and
!      ! make sure it's ok with everyone.
!      if (first_time) then
!         call error_handler(E_MSG, 'ROMS model_mod', &
!           'NOTE: Locations above the top Rho point are moved to that point')
!         first_time = .false.
!      endif
!      ! special case for obs at or closer to surface than top rho point
!      ! make them 100% at the top point
!      vloc = vloc - 1
!      hgt_fract = 1.0_r8
!   endif
!
!   ! fractional heights in the vertical for each corner of the horizontal box
!   hgt(1) = (1.0-hgt_fract)*(mydep(iloc  ,jloc,  vloc)) + hgt_fract*(mydep(iloc  ,jloc  ,vloc+1))
!   hgt(2) = (1.0-hgt_fract)*(mydep(iloc+1,jloc  ,vloc)) + hgt_fract*(mydep(iloc+1,jloc  ,vloc+1))
!   hgt(3) = (1.0-hgt_fract)*(mydep(iloc+1,jloc+1,vloc)) + hgt_fract*(mydep(iloc+1,jloc+1,vloc+1))
!   hgt(4) = (1.0-hgt_fract)*(mydep(iloc  ,jloc+1,vloc)) + hgt_fract*(mydep(iloc  ,jloc+1,vloc+1))
!
!   tmp_hgt(1) = (1.0-lon_fract)*hgt(1) + lon_fract*hgt(2)
!   tmp_hgt(2) = (1.0-lon_fract)*hgt(4) + lon_fract*hgt(3)
!
!   hgt_val   = (1.0-lat_fract)*tmp_hgt(1)+lat_fract*tmp_hgt(2)
!   vert_type = VERTISHEIGHT
!endif
!
!istatus  = 0
!location = set_location(lon_val, lat_val, hgt_val, vert_type)
!
!if (debug > 5) then
!   print*,' i,j,k', filoc, fjloc, fkloc
!
!   print*,' lon(i  ,j  ), lat(i  ,j  ) : (', mylon(iloc  ,jloc  ),',',mylat(iloc  ,jloc  ),')'
!   print*,' lon(i+1,j  ), lat(i+1,j  ) : (', mylon(iloc+1,jloc  ),',',mylat(iloc+1,jloc  ),')'
!   print*,' lon(i+1,j+1), lat(i+1,j+1) : (', mylon(iloc+1,jloc+1),',',mylat(iloc+1,jloc+1),')'
!   print*,' lon(i  ,j+1), lat(i  ,j+1) : (', mylon(iloc  ,jloc+1),',',mylon(iloc  ,jloc+1),')'
!   print*,' '
!   print*,' tmp_hgt(1)  = ', tmp_hgt(1)
!   print*,' tmp_hgt(2)  = ', tmp_hgt(2)
!   print*, ' '
!   print*,' lon_frac    = ', lon_fract
!   print*,' lat_frac    = ', lat_fract
!   print*,' hgt_frac    = ', hgt_fract
!   print*,' '
!   print*,' lon_val     = ', lon_val
!   print*,' lat_val     = ', lat_val
!   print*,' hgt_val     = ', hgt_val
!   print*,' '
!   print*,' WDEP(i  , j  , k  )', WDEP(iloc,   jloc  , vloc)
!   print*,' WDEP(i  , j+1, k  )', WDEP(iloc,   jloc+1, vloc)
!   print*,' WDEP(i+1, j+1, k  )', WDEP(iloc+1, jloc+1, vloc)
!   print*,' WDEP(i+1, j  , k  )', WDEP(iloc+1, jloc  , vloc)
!   print*,' '
!   print*,' WDEP(i  , j  , k+1)', WDEP(iloc,   jloc  , vloc+1)
!   print*,' WDEP(i  , j+1, k+1)', WDEP(iloc,   jloc+1, vloc+1)
!   print*,' WDEP(i+1, j+1, k+1)', WDEP(iloc+1, jloc+1, vloc+1)
!   print*,' WDEP(i+1, j  , k+1)', WDEP(iloc+1, jloc  , vloc+1)
!endif
!
!end function get_location_from_ijk


!===================================================================
! End of model_mod
!===================================================================

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

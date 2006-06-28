! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next six lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Id$
! $Name$

! Assimilation interface for WRF model

!-----------------------------------------------------------------------
!
!     interface for WRF
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------
!-----------------------------------------------------------------------

use         types_mod, only : r8, deg2rad, missing_r8, ps0, earth_radius, &
                              gas_constant, gas_constant_v, gravity
use  time_manager_mod, only : time_type, set_time, set_calendar_type, GREGORIAN
use      location_mod, only : location_type, get_location, set_location, &
                              get_dist, horiz_dist_only, &
                              LocationDims, LocationName, LocationLName, &
                              query_location, vert_is_undef, vert_is_surface, &
                              vert_is_level, vert_is_pressure, vert_is_height
use     utilities_mod, only : file_exist, open_file, close_file, &
                              register_module, error_handler, E_ERR, &
                              E_MSG, logfileunit, find_namelist_in_file, check_namelist_read
use      obs_kind_mod, only : KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, &
                              KIND_SURFACE_PRESSURE, KIND_TEMPERATURE, &
                              KIND_SPECIFIC_HUMIDITY, &
                              KIND_PRESSURE, KIND_VERTICAL_VELOCITY, &
                              KIND_RAINWATER_MIXING_RATIO, KIND_DENSITY, &
                              KIND_GRAUPEL_MIXING_RATIO, KIND_SNOW_MIXING_RATIO, &
                              KIND_VORTEX_LAT, KIND_VORTEX_LON, &
                              KIND_VORTEX_PMIN, KIND_VORTEX_WMAX
use         map_utils, only : proj_info, map_init, map_set, latlon_to_ij, &
                              PROJ_LATLON, PROJ_MERC, PROJ_LC, PROJ_PS, &
                              ij_to_latlon, gridwind_to_truewind
use netcdf
use typesizes

implicit none
private

!  public routines and data for the WRF model

public ::  get_model_size,                    &
           get_state_meta_data,               &
           model_interpolate,                 &
           get_model_time_step,               &
           static_init_model,                 &
           model_get_close_states,            &
           pert_model_state,                  &
           nc_write_model_atts,               &
           nc_write_model_vars

!  public stubs 

public ::  adv_1step,       &
           end_model,       &
           init_time,       &
           init_conditions

!-----------------------------------------------------------------------
! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

!-----------------------------------------------------------------------
! Model namelist parameters with default values.
!
! center_search_size:  grid size for searching minimum pressure at grid point
! center_spline_scale: coarse grid to spline interp. fine grid ratio
!-----------------------------------------------------------------------

logical :: output_state_vector  = .false.     ! output prognostic variables
integer :: num_moist_vars       = 0
integer :: num_domains          = 1
integer :: calendar_type        = GREGORIAN
integer :: assimilation_period_seconds = 21600
logical :: surf_obs             = .false.
logical :: h_diab               = .false.
character(len = 72) :: adv_mod_command = 'wrf.exe'
integer :: center_search_size       = 25
integer :: center_spline_grid_scale = 10

namelist /model_nml/ output_state_vector, num_moist_vars, &
                     num_domains, calendar_type, surf_obs, h_diab, &
                     adv_mod_command, assimilation_period_seconds, &
                     center_search_size, center_spline_grid_scale

!-----------------------------------------------------------------------

! Private definition of domain map projection use by WRF

integer, parameter :: map_sphere = 0, map_lambert = 1, map_polar_stereo = 2, map_mercator = 3

! Private definition of model variable types

integer, parameter :: TYPE_U   = 1,   TYPE_V   = 2,  TYPE_W  = 3,  &
                      TYPE_GZ  = 4,   TYPE_T   = 5,  TYPE_MU = 6,  &
                      TYPE_QV  = 7,   TYPE_QC  = 8,  TYPE_QR = 9,  &
                      TYPE_QI  = 10,  TYPE_QS  = 11, TYPE_QG = 12, &
                      TYPE_U10 = 13,  TYPE_V10 = 14, TYPE_T2 = 15, &
                      TYPE_Q2  = 16,  TYPE_PS  = 17, TYPE_TSLB = 18, &
                      TYPE_TSK = 19,  TYPE_HDIAB = 20, TYPE_QNICE = 21


real (kind=r8), PARAMETER    :: kappa = 2.0_r8/7.0_r8 ! gas_constant / cp
real (kind=r8), PARAMETER    :: ts0 = 300.0_r8        ! Base potential temperature for all levels.

!---- private data ----

TYPE wrf_static_data_for_dart

   integer  :: bt, bts, sn, sns, we, wes, sls
   real(r8) :: dx, dy, dt, p_top
   integer  :: map_proj
   real(r8) :: cen_lat,cen_lon
   type(proj_info) :: proj

   integer  :: n_moist
   logical  :: surf_obs
   real(r8), dimension(:),     pointer :: znu, dn, dnw, zs
   real(r8), dimension(:,:),   pointer :: mub, latitude, longitude, hgt
   real(r8), dimension(:,:),   pointer :: mapfac_m, mapfac_u, mapfac_v
   real(r8), dimension(:,:,:), pointer :: phb

   integer :: number_of_wrf_variables
   integer, dimension(:,:), pointer :: var_index
   integer, dimension(:,:), pointer :: var_size
   integer, dimension(:),   pointer :: var_type
   integer, dimension(:,:), pointer :: land

   integer, dimension(:,:,:,:), pointer :: dart_ind

end type wrf_static_data_for_dart

type wrf_dom
   type(wrf_static_data_for_dart), pointer :: dom(:)
   integer :: model_size
end type wrf_dom

type(wrf_dom) :: wrf


contains

!#######################################################################

subroutine static_init_model()

! INitializes class data for WRF???

integer :: ncid
integer :: io, iunit

character (len=80)    :: name
character (len=1)     :: idom
logical, parameter    :: debug = .false.
integer               :: var_id, ind, i, j, k, id, dart_index, model_type

integer  :: proj_code
real(r8) :: stdlon,truelat1,truelat2,dt

!----------------------------------------------------------------------

! Register the module
call register_module(source, revision, revdate)

! Begin by reading the namelist input
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
write(logfileunit, nml=model_nml)
write(     *     , nml=model_nml)

allocate(wrf%dom(num_domains))

wrf%dom(:)%n_moist = num_moist_vars

if( num_moist_vars > 7) then
   write(*,'(''num_moist_vars = '',i3)')num_moist_vars
   call error_handler(E_ERR,'static_init_model', &
        'num_moist_vars is too large.', source, revision,revdate)
endif

wrf%dom(:)%surf_obs = surf_obs

if ( debug ) then
   if ( output_state_vector ) then
      write(*,*)'netcdf file in state vector format'
   else
      write(*,*)'netcdf file in prognostic vector format'
   endif
endif

call set_calendar_type(calendar_type)

dart_index = 1

call read_dt_from_wrf_nml()

do id=1,num_domains

   write( idom , '(I1)') id

   write(*,*) '******************'
   write(*,*) '**  DOMAIN # ',idom,'  **'
   write(*,*) '******************'

   if(file_exist('wrfinput_d0'//idom)) then

      call check( nf90_open('wrfinput_d0'//idom, NF90_NOWRITE, ncid) )

   else

      call error_handler(E_ERR,'static_init_model', &
           'Please put wrfinput_d0'//idom//' in the work directory.', source, revision,revdate)

   endif

   if(debug) write(*,*) ' ncid is ',ncid

! get wrf grid dimensions

   call check( nf90_inq_dimid(ncid, "bottom_top", var_id) )
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%bt) )

   call check( nf90_inq_dimid(ncid, "bottom_top_stag", var_id) ) ! reuse var_id, no harm
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%bts) )

   call check( nf90_inq_dimid(ncid, "south_north", var_id) )
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%sn) )

   call check( nf90_inq_dimid(ncid, "south_north_stag", var_id)) ! reuse var_id, no harm
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%sns) )

   call check( nf90_inq_dimid(ncid, "west_east", var_id) )
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%we) )

   call check( nf90_inq_dimid(ncid, "west_east_stag", var_id) )  ! reuse var_id, no harm
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%wes) )

   call check( nf90_inq_dimid(ncid, "soil_layers_stag", var_id) )  ! reuse var_id, no harm
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%sls) )

   if(debug) then
      write(*,*) ' dimensions bt, sn, we are ',wrf%dom(id)%bt, &
           wrf%dom(id)%sn, wrf%dom(id)%we
      write(*,*) ' staggered  bt, sn, we are ',wrf%dom(id)%bts, &
           wrf%dom(id)%sns,wrf%dom(id)%wes
   endif

! get meta data and static data we need

   call check( nf90_get_att(ncid, nf90_global, 'DX', wrf%dom(id)%dx) )
   call check( nf90_get_att(ncid, nf90_global, 'DY', wrf%dom(id)%dy) )
   call check( nf90_get_att(ncid, nf90_global, 'DT', dt) )
   print*,'dt from wrfinput is: ',dt
   print*,'Using dt from namelist.input: ',wrf%dom(id)%dt
   if(debug) write(*,*) ' dx, dy, dt are ',wrf%dom(id)%dx, &
        wrf%dom(id)%dy, wrf%dom(id)%dt

   call check( nf90_get_att(ncid, nf90_global, 'MAP_PROJ', wrf%dom(id)%map_proj) )
   if(debug) write(*,*) ' map_proj is ',wrf%dom(id)%map_proj

   call check( nf90_get_att(ncid, nf90_global, 'CEN_LAT', wrf%dom(id)%cen_lat) )
   if(debug) write(*,*) ' cen_lat is ',wrf%dom(id)%cen_lat

   call check( nf90_get_att(ncid, nf90_global, 'CEN_LON', wrf%dom(id)%cen_lon) )

   call check( nf90_get_att(ncid, nf90_global, 'TRUELAT1', truelat1) )
   if(debug) write(*,*) ' truelat1 is ',truelat1

   call check( nf90_get_att(ncid, nf90_global, 'TRUELAT2', truelat2) )
   if(debug) write(*,*) ' truelat2 is ',truelat2

   call check( nf90_get_att(ncid, nf90_global, 'STAND_LON', stdlon) )

   call check( nf90_inq_varid(ncid, "P_TOP", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%p_top) )

!  get 1D (z) static data defining grid levels

   allocate(wrf%dom(id)%dn(1:wrf%dom(id)%bt))
   call check( nf90_inq_varid(ncid, "DN", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%dn) )
   if(debug) write(*,*) ' dn ',wrf%dom(id)%dn

   allocate(wrf%dom(id)%znu(1:wrf%dom(id)%bt))
   call check( nf90_inq_varid(ncid, "ZNU", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%znu) )
   if(debug) write(*,*) ' znu is ',wrf%dom(id)%znu

   allocate(wrf%dom(id)%dnw(1:wrf%dom(id)%bt))
   call check( nf90_inq_varid(ncid, "DNW", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%dnw) )
   if(debug) write(*,*) ' dnw is ',wrf%dom(id)%dnw

   allocate(wrf%dom(id)%zs(1:wrf%dom(id)%sls))
   call check( nf90_inq_varid(ncid, "ZS", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%zs) )

!  get 2D (x,y) base state for mu, latitude, longitude

   allocate(wrf%dom(id)%mub(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "MUB", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%mub) )
   if(debug) then
      write(*,*) ' corners of mub '
      write(*,*) wrf%dom(id)%mub(1,1),wrf%dom(id)%mub(wrf%dom(id)%we,1),  &
           wrf%dom(id)%mub(1,wrf%dom(id)%sn),wrf%dom(id)%mub(wrf%dom(id)%we, &
           wrf%dom(id)%sn)
   end if

   allocate(wrf%dom(id)%longitude(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "XLONG", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%longitude) )

   allocate(wrf%dom(id)%latitude(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "XLAT", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%latitude) )

   allocate(wrf%dom(id)%land(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "XLAND", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%land) )
   if(debug) then
      write(*,*) ' corners of land '
      write(*,*) wrf%dom(id)%land(1,1),wrf%dom(id)%land(wrf%dom(id)%we,1),  &
           wrf%dom(id)%land(1,wrf%dom(id)%sn),wrf%dom(id)%land(wrf%dom(id)%we, &
           wrf%dom(id)%sn)
   end if

   if(debug) then
      write(*,*) ' corners of lat '
      write(*,*) wrf%dom(id)%latitude(1,1),wrf%dom(id)%latitude(wrf%dom(id)%we,1),  &
           wrf%dom(id)%latitude(1,wrf%dom(id)%sn), &
           wrf%dom(id)%latitude(wrf%dom(id)%we,wrf%dom(id)%sn)
      write(*,*) ' corners of long '
      write(*,*) wrf%dom(id)%longitude(1,1),wrf%dom(id)%longitude(wrf%dom(id)%we,1),  &
           wrf%dom(id)%longitude(1,wrf%dom(id)%sn), &
           wrf%dom(id)%longitude(wrf%dom(id)%we,wrf%dom(id)%sn)
   end if

   allocate(wrf%dom(id)%mapfac_m(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "MAPFAC_M", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%mapfac_m) )

   allocate(wrf%dom(id)%hgt(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "HGT", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%hgt) )

   allocate(wrf%dom(id)%mapfac_u(1:wrf%dom(id)%wes,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "MAPFAC_U", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%mapfac_u) )

   allocate(wrf%dom(id)%mapfac_v(1:wrf%dom(id)%we,1:wrf%dom(id)%sns))
   call check( nf90_inq_varid(ncid, "MAPFAC_V", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%mapfac_v) )

! get 3D base state geopotential

   allocate(wrf%dom(id)%phb(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%bts))
   call check( nf90_inq_varid(ncid, "PHB", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%phb) )
   if(debug) then
      write(*,*) ' corners of phb '
      write(*,*) wrf%dom(id)%phb(1,1,1),wrf%dom(id)%phb(wrf%dom(id)%we,1,1),  &
           wrf%dom(id)%phb(1,wrf%dom(id)%sn,1),wrf%dom(id)%phb(wrf%dom(id)%we, &
           wrf%dom(id)%sn,1)
      write(*,*) wrf%dom(id)%phb(1,1,wrf%dom(id)%bts), &
           wrf%dom(id)%phb(wrf%dom(id)%we,1,wrf%dom(id)%bts),  &
           wrf%dom(id)%phb(1,wrf%dom(id)%sn,wrf%dom(id)%bts), &
           wrf%dom(id)%phb(wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bts)
   end if

! close data file, we have all we need

   call check( nf90_close(ncid) )

! Initializes the map projection structure to missing values

   call map_init(wrf%dom(id)%proj)

! Populate the map projection structure

   if(wrf%dom(id)%map_proj == map_sphere) then
      proj_code = PROJ_LATLON
   elseif(wrf%dom(id)%map_proj == map_lambert) then
      proj_code = PROJ_LC
   elseif(wrf%dom(id)%map_proj == map_polar_stereo) then
      proj_code = PROJ_PS
   elseif(wrf%dom(id)%map_proj == map_mercator) then
      proj_code = PROJ_MERC
   else
      call error_handler(E_ERR,'static_init_model', &
        'Map projection no supported.', source, revision, revdate)
   endif

   call map_set(proj_code,wrf%dom(id)%latitude(1,1),wrf%dom(id)%longitude(1,1), &
        1.0_r8,1.0_r8,wrf%dom(id)%dx,stdlon,truelat1,truelat2,wrf%dom(id)%proj)

!  build the map into the 1D DART vector for WRF data

   wrf%dom(id)%number_of_wrf_variables = 8 + wrf%dom(id)%n_moist
   if( wrf%dom(id)%surf_obs ) then
      wrf%dom(id)%number_of_wrf_variables = wrf%dom(id)%number_of_wrf_variables + 5
   endif
   if( h_diab ) then
      wrf%dom(id)%number_of_wrf_variables = wrf%dom(id)%number_of_wrf_variables + 1
   endif
   allocate(wrf%dom(id)%var_type(wrf%dom(id)%number_of_wrf_variables))
   wrf%dom(id)%var_type(1)  = TYPE_U
   wrf%dom(id)%var_type(2)  = TYPE_V
   wrf%dom(id)%var_type(3)  = TYPE_W
   wrf%dom(id)%var_type(4)  = TYPE_GZ
   wrf%dom(id)%var_type(5)  = TYPE_T
   wrf%dom(id)%var_type(6)  = TYPE_MU
   wrf%dom(id)%var_type(7)  = TYPE_TSLB
   wrf%dom(id)%var_type(8)  = TYPE_TSK
   ind = 8
   if( wrf%dom(id)%n_moist >= 1) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind)  = TYPE_QV
   end if
   if( wrf%dom(id)%n_moist >= 2) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind)  = TYPE_QC
   end if
   if( wrf%dom(id)%n_moist >= 3) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_QR
   end if
   if( wrf%dom(id)%n_moist >= 4) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_QI
   end if
   if( wrf%dom(id)%n_moist >= 5) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_QS
   end if
   if( wrf%dom(id)%n_moist >= 6) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_QG
   end if
   if( wrf%dom(id)%n_moist == 7) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_QNICE
   end if
   if( wrf%dom(id)%surf_obs ) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_U10
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_V10
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_T2
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_Q2
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_PS
   end if
   if( h_diab ) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_HDIAB
   end if

! indices into 1D array
   allocate(wrf%dom(id)%dart_ind(wrf%dom(id)%wes,wrf%dom(id)%sns,wrf%dom(id)%bts,20))
   allocate(wrf%dom(id)%var_index(2,wrf%dom(id)%number_of_wrf_variables))
! dimension of variables
   allocate(wrf%dom(id)%var_size(3,wrf%dom(id)%number_of_wrf_variables))

   wrf%dom(id)%dart_ind = 0

   ind = 1                         ! *** u field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%wes
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_U) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** v field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sns
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_V) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** w field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bts
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_W) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** geopotential field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bts
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_GZ) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** theta field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_T) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** mu field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = 1
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_MU) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** tslb field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%sls
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_TSLB) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** tsk field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = 1
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_TSK) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   do model_type = TYPE_QV, TYPE_QV + wrf%dom(id)%n_moist - 1
      ind = ind + 1                   ! *** moisture field ***
      wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
      wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
      wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
      wrf%dom(id)%var_index(1,ind) = dart_index
      do k=1,wrf%dom(id)%var_size(3,ind)
         do j=1,wrf%dom(id)%var_size(2,ind)
            do i=1,wrf%dom(id)%var_size(1,ind)
               wrf%dom(id)%dart_ind(i,j,k,model_type) = dart_index
               dart_index = dart_index + 1
            enddo
         enddo
      enddo
      wrf%dom(id)%var_index(2,ind) = dart_index - 1
   enddo

   if(wrf%dom(id)%surf_obs ) then
      do model_type = TYPE_U10, TYPE_PS
         ind = ind + 1                   ! *** Surface variable ***
         wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
         wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
         wrf%dom(id)%var_size(3,ind) = 1
         wrf%dom(id)%var_index(1,ind) = dart_index
         do k=1,wrf%dom(id)%var_size(3,ind)
            do j=1,wrf%dom(id)%var_size(2,ind)
               do i=1,wrf%dom(id)%var_size(1,ind)
                  wrf%dom(id)%dart_ind(i,j,k,model_type) = dart_index
                  dart_index = dart_index + 1
               enddo
            enddo
         enddo
         wrf%dom(id)%var_index(2,ind) = dart_index - 1
      enddo
   end if

   if(h_diab ) then
      ind = ind + 1                   ! *** h_diabatic variable ***
      wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
      wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
      wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
      wrf%dom(id)%var_index(1,ind) = dart_index
      do k=1,wrf%dom(id)%var_size(3,ind)
         do j=1,wrf%dom(id)%var_size(2,ind)
            do i=1,wrf%dom(id)%var_size(1,ind)
               wrf%dom(id)%dart_ind(i,j,k,TYPE_HDIAB) = dart_index
               dart_index = dart_index + 1
            enddo
         enddo
      enddo
      wrf%dom(id)%var_index(2,ind) = dart_index - 1
   end if

enddo

write(*,*)

wrf%model_size = dart_index - 1
if(debug) write(*,*) ' wrf model size is ',wrf%model_size

contains

  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'static_init_model', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end subroutine static_init_model


!#######################################################################

function get_model_size()

integer :: get_model_size

get_model_size = wrf%model_size

end function get_model_size

!#######################################################################

function get_model_time_step()
!------------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the time step of the model. In the long run should be replaced
! by a more general routine that returns details of a general time-stepping
! capability.
!
! toward that end ... we are now reading a namelist variable for the
! width of the assimilation time window.

type(time_type) :: get_model_time_step
integer :: model_dt, assim_dt

! We need to coordinate the desired assimilation window to be a 
! multiple of the model time step (which has no precision past integer seconds).

model_dt = nint(wrf%dom(1)%dt)

! The integer arithmetic does its magic.
assim_dt = (assimilation_period_seconds / model_dt) * model_dt

get_model_time_step = set_time(assim_dt)

end function get_model_time_step

!#######################################################################

subroutine get_state_meta_data(index_in, location, var_type_out, id_out)
!---------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type_out, id_out

integer  :: var_type
integer  :: index, ip, jp, kp
integer  :: nz, ny, nx
logical  :: var_found
real(r8) :: lon, lat, lev

integer :: i, id
logical, parameter :: debug = .false.
character(len=129) :: errstring

if(debug) then
   write(errstring,*)' index_in = ',index_in
   call error_handler(E_MSG,'get_state_meta_data',errstring,' ',' ',' ')
endif

index = index_in
var_found = .false.

!  first find var_type

if(debug) then
   do id=1,num_domains
      do i=1, wrf%dom(id)%number_of_wrf_variables
         write(errstring,*)' domain, var, var_type(i) = ',id,i,wrf%dom(id)%var_type(i)
         call error_handler(E_MSG,'get_state_meta_data',errstring,' ',' ',' ')
      enddo
   enddo
endif

i = 0
id = 1
do while (.not. var_found)
   i = i + 1
   if(i .gt. wrf%dom(id)%number_of_wrf_variables) then
      i = 1
      if (id < num_domains) then
         id = id + 1
      else
         write(errstring,*)' size of vector ',wrf%model_size
         call error_handler(E_MSG,'get_state_meta_data', errstring, ' ', ' ', ' ')
         write(errstring,*)' dart_index ',index_in
         call error_handler(E_ERR,'get_state_meta_data', 'index out of range', &
              source, revision, revdate)
      end if
   end if
   if( (index .ge. wrf%dom(id)%var_index(1,i) ) .and.  &
       (index .le. wrf%dom(id)%var_index(2,i) )       )  then
      var_found = .true.
      var_type = wrf%dom(id)%var_type(i)
      index = index - wrf%dom(id)%var_index(1,i) + 1
   end if
end do

!  now find i,j,k location.
!  index has been normalized such that it is relative to
!  array starting at (1,1,1)

nx = wrf%dom(id)%var_size(1,i)
ny = wrf%dom(id)%var_size(2,i)
nz = wrf%dom(id)%var_size(3,i)

kp = 1 + (index-1)/(nx*ny)
jp = 1 + (index - (kp-1)*nx*ny - 1)/nx
ip = index - (kp-1)*nx*ny - (jp-1)*nx

if(debug) write(*,*) ' ip, jp, kp for index ',ip,jp,kp,index
if(debug) write(*,*) ' Var type: ',var_type

call get_wrf_horizontal_location( ip, jp, var_type, id, lon, lat )

if( (var_type == type_w ) .or. (var_type == type_gz) ) then
   lev = real(kp) - 0.5_r8 ! This is the index of the vertical
else
   lev = real(kp) ! This is the index of the vertical
endif

if(debug) write(*,*) 'lon, lat, lev: ',lon, lat, lev

! lev is an index here, so which_vert is OK to be hardwired to a 1
location = set_location(lon, lat, lev, 1) 

if(present(var_type_out)) var_type_out = var_type
if(present(id_out)) id_out = id

end subroutine get_state_meta_data

!#######################################################################


 subroutine model_interpolate(x, location, obs_kind, obs_val, istatus)
!----------------------------------------------------------------------
!subroutine model_interpolate(x, location, obs_kind, obs_val, istatus)
!
! modified 26 June 2006 to accomodate vortex attributes

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_kind
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

logical, parameter  :: debug = .false.
real(r8)            :: xloc, yloc, zloc, xloc_u, yloc_v, xyz_loc(3)
integer             :: i, i_u, j, j_v, k, k2, i1,i2
real(r8)            :: dx,dy,dz,dxm,dym,dzm,dx_u,dxm_u,dy_v,dym_v
real(r8)            :: a1,utrue,vtrue,ugrid,vgrid
integer             :: id

character(len=129) :: errstring

real(r8), dimension(2) :: fld
real(r8), allocatable, dimension(:) :: v_h, v_p

! local vars, used in finding sea-level pressure and vortex center
real(r8), allocatable, dimension(:)   :: t1d, p1d, qv1d, z1d
real(r8), allocatable, dimension(:,:) :: psea, pp, pd
real(r8), allocatable, dimension(:)   :: x1d,y1d,xx1d,yy1d
integer  :: xlen, ylen, xxlen, yylen, ii1, ii2
real(r8) :: clat, clon, cxmin, cymin

! center_track_*** used to define center search area
integer :: center_track_xmin, center_track_ymin, &
           center_track_xmax, center_track_ymax

! local vars, used in calculating density, pressure
real(r8)            :: rho1 , rho2 , rho3, rho4
real(r8)            :: pres1, pres2, pres3, pres4, pres

logical  :: dom_found

dom_found = .false.

istatus = 0

xyz_loc = get_location(location)

id = num_domains
do while (.not. dom_found)

   ! Checking for exact equality on real variable types is generally a bad idea.

   if( (wrf%dom(id)%proj%hemi ==  1.0_r8 .and. xyz_loc(2) == -90.0_r8) .or. &
       (wrf%dom(id)%proj%hemi == -1.0_r8 .and. xyz_loc(2) ==  90.0_r8) .or. &
       (wrf%dom(id)%proj%code == PROJ_MERC .and. abs(xyz_loc(2)) >= 90.0_r8) ) then

   else

      call latlon_to_ij(wrf%dom(id)%proj, xyz_loc(2), xyz_loc(1), xloc,yloc)

      if ( (xloc >= 1 .and. xloc <= wrf%dom(id)%we .and. &
           yloc >= 1 .and. yloc <= wrf%dom(id)%sn) ) then

         dom_found = .true.

      endif

   endif

   if (.not. dom_found) then
      if (id == 1) then
         obs_val = missing_r8
         istatus = 1
         return
      else
         id = id - 1
      endif
   endif

end do

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

allocate(v_h(0:wrf%dom(id)%bt), v_p(0:wrf%dom(id)%bt))

call toGrid(xloc,i,dx,dxm)
call toGrid(yloc,j,dy,dym)

if(vert_is_level(location)) then

   ! If obs is by model level
   zloc = xyz_loc(3)
   if(debug) print*,' obs is by model level and zloc =',zloc
else if(vert_is_pressure(location)) then
   if(xyz_loc(3) < 10000.0_r8) then
      obs_val = missing_r8
      istatus = 1
      deallocate(v_h, v_p)
      return
   endif
   ! get model pressure profile
   call get_model_pressure_profile(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,x,id,v_p)
   ! get pressure vertical co-ordinate
   call pres_to_zk(xyz_loc(3), v_p, wrf%dom(id)%bt,zloc)
   if(debug.and.obs_kind /= KIND_SURFACE_PRESSURE) print*,' obs is by pressure and zloc =',zloc
   if(debug) print*,'model pressure profile'
   if(debug) print*,v_p
else if(vert_is_height(location)) then
   ! get model height profile
   call get_model_height_profile(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,x,id,v_h)
   if(debug) print*,'model height profile'
   if(debug) print*,v_h
   ! get height vertical co-ordinate
   call height_to_zk(xyz_loc(3), v_h, wrf%dom(id)%bt,zloc)
   if(debug) print*,' obs is by height and zloc =',zloc
else if(vert_is_surface(location)) then
   ! get terrain height
   v_h(0) = dym*( dxm*wrf%dom(id)%hgt(i,j) + &
                   dx*wrf%dom(id)%hgt(i+1,j) ) + &
             dy*( dxm*wrf%dom(id)%hgt(i,j+1) + &
                   dx*wrf%dom(id)%hgt(i+1,j+1) )
   zloc = 0.0_r8
   if(debug) print*,' obs is at the surface = ', xyz_loc(3)
else
   write(errstring,*) 'wrong option for which_vert ',nint(query_location(location,'which_vert'))
   call error_handler(E_ERR,'model_interpolate', errstring, &
        source, revision, revdate)
end if

if(zloc == missing_r8) then
   obs_val = missing_r8
   istatus = 1
   deallocate(v_h, v_p)
   return
endif

k = max(1,int(zloc))

! Get the desired field to be interpolated
if( obs_kind == KIND_U_WIND_COMPONENT .or. obs_kind == KIND_V_WIND_COMPONENT) then        ! U, V

   if(.not. vert_is_surface(location)) then

      xloc_u = xloc + 0.5
      yloc_v = yloc + 0.5
      call toGrid(xloc_u,i_u,dx_u,dxm_u)
      call toGrid(yloc_v,j_v,dy_v,dym_v)

      if(i_u >= 1 .and. i_u < wrf%dom(id)%var_size(1,TYPE_U) .and. &
           j   >= 1 .and. j   < wrf%dom(id)%var_size(2,TYPE_U) .and. &
           i   >= 1 .and. i   < wrf%dom(id)%var_size(1,TYPE_V) .and. &
           j_v >= 1 .and. j_v < wrf%dom(id)%var_size(2,TYPE_V)) then

         do k2=1,2

            ! For memory economy over speed, dart_ind should be removed from
            ! the wrf structure and get_wrf_index reactivated.

!!$            i1 = get_wrf_index(i_u,j  ,k+k2-1,TYPE_U,id)
!!$            i2 = get_wrf_index(i_u,j+1,k+k2-1,TYPE_U,id)

            i1 = wrf%dom(id)%dart_ind(i_u,j  ,k+k2-1,TYPE_U)
            i2 = wrf%dom(id)%dart_ind(i_u,j+1,k+k2-1,TYPE_U)

            ugrid = dym*( dxm_u*x(i1) + dx_u*x(i1+1) ) + dy*( dxm_u*x(i2) + dx_u*x(i2+1) )

!!$            i1 = get_wrf_index(i,j_v  ,k+k2-1,TYPE_V,id)
!!$            i2 = get_wrf_index(i,j_v+1,k+k2-1,TYPE_V,id) 

            i1 = wrf%dom(id)%dart_ind(i,j_v  ,k+k2-1,TYPE_V)
            i2 = wrf%dom(id)%dart_ind(i,j_v+1,k+k2-1,TYPE_V) 

            vgrid = dym_v*( dxm*x(i1) + dx*x(i1+1) ) + dy_v*( dxm*x(i2) + dx*x(i2+1) )

            call gridwind_to_truewind(xyz_loc(1), wrf%dom(id)%proj, ugrid, vgrid, &
                 utrue, vtrue)

            if( obs_kind == KIND_U_WIND_COMPONENT) then

               fld(k2) = utrue

            else   ! must want v

               fld(k2) = vtrue

            endif

         enddo

      else

         fld(:) = missing_r8

      endif

   else

      if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
         j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
         wrf%dom(id)%surf_obs) then

         if( obs_kind == KIND_U_WIND_COMPONENT) then

!!$            i1 = get_wrf_index(i,j,1,TYPE_U10,id)
!!$            i2 = get_wrf_index(i,j+1,1,TYPE_U10,id)
            i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_U10)
            i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_U10)

         else
!!$            i1 = get_wrf_index(i,j,1,TYPE_V10,id)
!!$            i2 = get_wrf_index(i,j+1,1,TYPE_V10,id)
            i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_V10)
            i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_V10)

         endif

         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

      else

         fld(1) = missing_r8

      endif

   endif

else if( obs_kind == KIND_TEMPERATURE ) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)) then

      if(.not. vert_is_surface(location)) then

!!$         i1 = get_wrf_index(i,j  ,k,TYPE_T,id)
!!$         i2 = get_wrf_index(i,j+1,k,TYPE_T,id)
         i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_T)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_T)
         a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
         pres1 = model_pressure_t(i  ,j  ,k,id,x)
         pres2 = model_pressure_t(i+1,j  ,k,id,x)
         pres3 = model_pressure_t(i  ,j+1,k,id,x)
         pres4 = model_pressure_t(i+1,j+1,k,id,x)
         pres = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )
         fld(1) = (ts0 + a1)*(pres/ps0)**kappa

!!$         i1 = get_wrf_index(i,j  ,k+1,TYPE_T,id)
!!$         i2 = get_wrf_index(i,j+1,k+1,TYPE_T,id)
         i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_T)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_T)
         a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
         pres1 = model_pressure_t(i  ,j  ,k+1,id,x)
         pres2 = model_pressure_t(i+1,j  ,k+1,id,x)
         pres3 = model_pressure_t(i  ,j+1,k+1,id,x)
         pres4 = model_pressure_t(i+1,j+1,k+1,id,x)
         pres = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )
         fld(2) = (ts0 + a1)*(pres/ps0)**kappa

      else

         if(wrf%dom(id)%surf_obs) then

!!$            i1 = get_wrf_index(i,j,1,TYPE_T2,id)
!!$            i2 = get_wrf_index(i,j+1,1,TYPE_T2,id)
            i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_T2)
            i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_T2)
            fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

         else

            fld(1) = missing_r8

         endif

      endif

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_DENSITY ) then

   ! Rho calculated at mass points, and so is like "TYPE_T" 
   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)) then

      ! calculate full rho at corners of interp box
      ! and interpolate to desired horizontal location

      rho1 = model_rho_t(i  ,j  ,k,id,x)
      rho2 = model_rho_t(i+1,j  ,k,id,x)
      rho3 = model_rho_t(i  ,j+1,k,id,x)
      rho4 = model_rho_t(i+1,j+1,k,id,x)
      fld(1) = dym*( dxm*rho1 + dx*rho2 ) + dy*( dxm*rho3 + dx*rho4 )

      rho1 = model_rho_t(i  ,j  ,k+1,id,x)
      rho2 = model_rho_t(i+1,j  ,k+1,id,x)
      rho3 = model_rho_t(i  ,j+1,k+1,id,x)
      rho4 = model_rho_t(i+1,j+1,k+1,id,x)
      fld(2) = dym*( dxm*rho1 + dx*rho2 ) + dy*( dxm*rho3 + dx*rho4 )

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_VERTICAL_VELOCITY ) then

   zloc = zloc + 0.5_r8
   k = max(1,int(zloc))

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_W) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_W)) then

!!$      i1 = get_wrf_index(i,j  ,k,TYPE_W,id)
!!$      i2 = get_wrf_index(i,j+1,k,TYPE_W,id)
      i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_W)
      i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_W)
      fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

!!$      i1 = get_wrf_index(i,j  ,k+1,TYPE_W,id)
!!$      i2 = get_wrf_index(i,j+1,k+1,TYPE_W,id)
      i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_W)
      i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_W)
      fld(2) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_SPECIFIC_HUMIDITY ) then

!  Convert water vapor mixing ratio to specific humidity:

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)) then

      if ( wrf%dom(id)%n_moist >= 1) then

         if(.not. vert_is_surface(location)) then

!!$            i1 = get_wrf_index(i,j  ,k,TYPE_QV,id)
!!$            i2 = get_wrf_index(i,j+1,k,TYPE_QV,id)
            i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_QV)
            i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_QV)
            a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
            fld(1) = a1 /(1.0_r8 + a1)

!!$            i1 = get_wrf_index(i,j  ,k+1,TYPE_QV,id)
!!$            i2 = get_wrf_index(i,j+1,k+1,TYPE_QV,id)
            i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_QV)
            i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_QV)
            a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
            fld(2) = a1 /(1.0_r8 + a1)

         else

            if(wrf%dom(id)%surf_obs) then

!!$               i1 = get_wrf_index(i,j,1,TYPE_Q2,id)
!!$               i2 = get_wrf_index(i,j+1,1,TYPE_Q2,id)
               i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_Q2)
               i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_Q2)
               fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
               fld(1) = fld(1) / (1.0_r8 + fld(1))

            else

               fld(1) = missing_r8

            endif

         endif

      else

         fld(:) = 0.0_r8

      endif

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_RAINWATER_MIXING_RATIO) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)) then

      if( wrf%dom(id)%n_moist >= 3) then

!!$      i1 = get_wrf_index(i,j  ,k,TYPE_QR,id)
!!$      i2 = get_wrf_index(i,j+1,k,TYPE_QR,id)
         i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_QR)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_QR)
         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

!!$      i1 = get_wrf_index(i,j  ,k+1,TYPE_QR,id)
!!$      i2 = get_wrf_index(i,j+1,k+1,TYPE_QR,id)
         i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_QR)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_QR)
         fld(2) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

         fld = max(0.0_r8, fld)

      else

         fld(:) = 0.0_r8

      endif


   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_GRAUPEL_MIXING_RATIO) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)) then

      if( wrf%dom(id)%n_moist >= 6) then

         i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_QG)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_QG)
         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

         i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_QG)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_QG)
         fld(2) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
         fld = max(0.0_r8, fld)

      else

         fld(:) = 0.0_r8

      endif


   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_SNOW_MIXING_RATIO) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)) then

      if( wrf%dom(id)%n_moist >= 5) then

         i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_QS)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_QS)
         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

         i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_QS)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_QS)
         fld(2) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

         fld = max(0.0_r8, fld)

      else

         fld(:) = 0.0_r8

      endif

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_PRESSURE) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)) then

      pres1 = model_pressure_t(i  ,j  ,k,id,x)
      pres2 = model_pressure_t(i+1,j  ,k,id,x)
      pres3 = model_pressure_t(i  ,j+1,k,id,x)
      pres4 = model_pressure_t(i+1,j+1,k,id,x)
      fld(1) = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )

      pres1 = model_pressure_t(i  ,j  ,k+1,id,x)
      pres2 = model_pressure_t(i+1,j  ,k+1,id,x)
      pres3 = model_pressure_t(i  ,j+1,k+1,id,x)
      pres4 = model_pressure_t(i+1,j+1,k+1,id,x)
      fld(2) = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_SURFACE_PRESSURE) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
      wrf%dom(id)%surf_obs) then

!!$      i1 = get_wrf_index(i,j,1,TYPE_PS,id)
!!$      i2 = get_wrf_index(i,j+1,1,TYPE_PS,id)
      i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_PS)
      i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_PS)
      if(x(i1) /= 0.0_r8 .and. x(i1+1) /= 0.0_r8 .and. &
           x(i2) /= 0.0_r8 .and. x(i2+1) /= 0.0_r8) then

         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

      else

         fld(1) = missing_r8

      endif

   else

      fld(1) = missing_r8

   endif

!!-------------------------------------------------------------------------
!! find vortex center location and with minimum sea level pressure
!! steps: 1) define a search box centered at the obs. location (i,j)
!!        2) compute sea-level pressure (slp)
!!        3) spline-interpolation of slp within the refined search box 
!!           -- by a factor of center_spline_grid_scale (default=10)
!!        4) find the minimum slp and grid point location 
!!           (default=1/10 precision)
!!        5) return the minimum slp and lat/lon of this grid point
!!-------------------------------------------------------------------------
else if( obs_kind == KIND_VORTEX_LAT .or. &
         obs_kind == KIND_VORTEX_LON .or. &
         obs_kind == KIND_VORTEX_PMIN .or. &
         obs_kind == KIND_VORTEX_WMAX) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)) then

!!   define a search box bounded by center_track_***
     center_track_xmin = max(1,i-center_search_size/2)
     center_track_xmax = min(wrf%dom(id)%var_size(1,TYPE_MU),i+center_search_size/2)
     center_track_ymin = max(1,j-center_search_size/2)
     center_track_ymax = min(wrf%dom(id)%var_size(2,TYPE_MU),j+center_search_size/2)
     if(center_track_xmin<1 .or. center_track_xmax>wrf%dom(id)%var_size(1,TYPE_MU) .or. &
        center_track_ymin<1 .or. center_track_ymax>wrf%dom(id)%var_size(2,TYPE_MU) .or. &
        center_track_xmin >= center_track_xmax .or. &
        center_track_ymin >= center_track_ymax) then
          print*,'i,j,center_search_size,center_track_xmin(max),center_track_ymin(max)'
          print*,i,j,center_search_size,center_track_xmin,center_track_xmax,center_track_ymin,center_track_ymax
         write(errstring,*)'Wrong setup in center_track_nml'
         call error_handler(E_ERR,'model_interpolate', errstring, source, revision, revdate)
     endif 

!!   define spline interpolation box dimensions
     xlen = center_track_xmax - center_track_xmin + 1
     ylen = center_track_ymax - center_track_ymin + 1
     xxlen = (center_track_xmax - center_track_xmin)*center_spline_grid_scale + 1
     yylen = (center_track_ymax - center_track_ymin)*center_spline_grid_scale + 1
     allocate(p1d(wrf%dom(id)%bt), t1d(wrf%dom(id)%bt))
     allocate(qv1d(wrf%dom(id)%bt), z1d(wrf%dom(id)%bt))
     allocate(psea(xlen,ylen))
     allocate(pd(xlen,ylen))
     allocate(pp(xxlen,yylen))
     allocate(x1d(xlen))
     allocate(y1d(ylen))
     allocate(xx1d(xxlen))
     allocate(yy1d(yylen))

!!   compute sea-level pressure
     do i1 = center_track_xmin, center_track_xmax
     do i2 = center_track_ymin, center_track_ymax
        do k = 1,wrf%dom(id)%var_size(3,TYPE_T)
           p1d(k) = model_pressure_t(i1,i2,k,id,x)
           t1d(k) = ts0 + x(wrf%dom(id)%dart_ind(i1,i2,k,TYPE_T))
           qv1d(k)= x(wrf%dom(id)%dart_ind(i1,i2,k,TYPE_QV))
           z1d(k) = ( x(wrf%dom(id)%dart_ind(i1,i2,k,TYPE_GZ))+wrf%dom(id)%phb(i1,i2,k) + &
                      x(wrf%dom(id)%dart_ind(i1,i2,k+1,TYPE_GZ))+wrf%dom(id)%phb(i1,i2,k+1) &
                    )*0.5_r8/gravity
        enddo
        call compute_seaprs(wrf%dom(id)%bt, z1d, t1d, p1d, qv1d, &
                          psea(i1-center_track_xmin+1,i2-center_track_ymin+1),debug)
     enddo
     enddo

!!   spline-interpolation
     do i1 = 1,xlen
        x1d(i1) = (i1-1)+center_track_xmin
     enddo
     do i2 = 1,ylen
        y1d(i2) = (i2-1)+center_track_ymin
     enddo
     do ii1 = 1,xxlen
        xx1d(ii1) = center_track_xmin+real(ii1-1,r8)*1_r8/real(center_spline_grid_scale,r8)
     enddo
     do ii2 = 1,yylen
        yy1d(ii2) = center_track_ymin+real(ii2-1,r8)*1_r8/real(center_spline_grid_scale,r8)
     enddo

     call splie2(x1d,y1d,psea,xlen,ylen,pd)

     pres1 = 1.e20
     cxmin = -1
     cymin = -1
     do ii1=1,xxlen
     do ii2=1,yylen
        call splin2(x1d,y1d,psea,pd,xlen,ylen,xx1d(ii1),yy1d(ii2),pp(ii1,ii2))
        if (pres1 .gt. pp(ii1,ii2)) then
           pres1=pp(ii1,ii2)
           cxmin=xx1d(ii1)
           cymin=yy1d(ii2)
        endif
     enddo
     enddo

     call ij_to_latlon(wrf%dom(id)%proj, cxmin, cymin, clat, clon)

     if( obs_kind == KIND_VORTEX_LAT ) then
        fld(1) = clat
     else if( obs_kind == KIND_VORTEX_LON ) then
        fld(1) = clon
     else if( obs_kind == KIND_VORTEX_PMIN ) then
        fld(1) = pres1
     else if( obs_kind == KIND_VORTEX_WMAX ) then
        fld(1) = missing_r8
        ! not implemented yet
     endif

     deallocate(p1d, t1d, qv1d, z1d)
     deallocate(psea,pd,pp,x1d,y1d,xx1d,yy1d)


   else

      fld(1) = missing_r8

   endif
else
   write(errstring,*)'do not recognize obs kind ',obs_kind
   call error_handler(E_ERR,'model_interpolate', errstring, &
        source, revision, revdate)
end if

! Do vertical interpolation
if(.not. vert_is_surface(location)) then

   call toGrid(zloc, k, dz, dzm)

   if( k >= 1 ) then
      obs_val = dzm*fld(1) + dz*fld(2)
   elseif( k == 0 ) then
      obs_val = fld(1) - (fld(2)-fld(1))*dzm
   else
      obs_val = missing_r8
   endif

else

   obs_val = fld(1)

endif

if(obs_val == missing_r8) istatus = 1

if(debug) print*,' interpolated value= ',obs_val

deallocate(v_h, v_p)

end subroutine model_interpolate

!#######################################################################

subroutine model_get_close_states(o_loc, radius, num_total, indices, dist, x)

type(location_type), intent(in) :: o_loc
real(r8),            intent(in) :: radius
integer,            intent(out) :: num_total, indices(:)
real(r8),           intent(out) :: dist(:)
real(r8),            intent(in) :: x(:)

integer                         :: u_pts, v_pts, p_pts, w_pts
integer                         :: i,k,indmax, ii, jj, kk, id
integer                         :: max_size
integer,            allocatable :: lon_ind(:), lat_ind(:), vert_ind(:)
real(r8),           allocatable :: close_dist(:)

indmax = size(indices)

num_total = 0

do id=1,num_domains

! We're allocating enough space for all u, v, w, and p points.

   max_size = (wrf%dom(id)%we*wrf%dom(id)%sn + wrf%dom(id)%wes*wrf%dom(id)%sn + &
        wrf%dom(id)%we*(wrf%dom(id)%sns))*wrf%dom(id)%bt + &
        wrf%dom(id)%we*wrf%dom(id)%sn*wrf%dom(id)%bts
   allocate(lon_ind(max_size), lat_ind(max_size), vert_ind(max_size), &
        close_dist(max_size))

! Look for close grid points.
   call grid_close_states( o_loc, &
        radius, id, lon_ind, lat_ind, vert_ind, close_dist, &
        u_pts, v_pts, p_pts, w_pts, x )

! start with u_pts

   do i = 1, u_pts

      ii = lon_ind(i)
      jj = lat_ind(i)
      k = vert_ind(i)

      num_total = num_total + 1
      if(num_total <= indmax) then
!!$         indices(num_total) = get_wrf_index(ii,jj,k,type_u,id)
         indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_u)
         dist(num_total) = close_dist(i)
      end if

   enddo

! now v_pts

   do i = u_pts+1, u_pts+v_pts

      ii = lon_ind(i)
      jj = lat_ind(i)
      k = vert_ind(i)

      num_total = num_total + 1
      if(num_total <= indmax) then
!!$         indices(num_total) = get_wrf_index(ii,jj,k,type_v,id)
         indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_v)
         dist(num_total) = close_dist(i)
      end if

   enddo

! now w_pts (w and gz)

   do i = u_pts+v_pts+1, u_pts+v_pts+w_pts

      ii = lon_ind(i)
      jj = lat_ind(i)
      k = vert_ind(i)

      num_total = num_total + 1
      if(num_total <= indmax) then
!!$         indices(num_total) = get_wrf_index(ii,jj,k,type_w,id)
         indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_w)
         dist(num_total) = close_dist(i)
      end if
      num_total = num_total + 1
      if(num_total <= indmax) then
!!$         indices(num_total) = get_wrf_index(ii,jj,k,type_gz,id)
         indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_gz)
         dist(num_total) = close_dist(i)
      end if

! 2D (near-surface) fields.

      if (k == 1) then

         num_total = num_total + 1
         if(num_total <= indmax) then
!!$            indices(num_total) = get_wrf_index(ii,jj,k,type_mu,id)
            indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_mu)
            dist(num_total) = close_dist(i)
         end if

         num_total = num_total + 1
         if(num_total <= indmax) then
!!$            indices(num_total) = get_wrf_index(ii,jj,k,type_tsk,id)
            indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_tsk)
            dist(num_total) = close_dist(i)
         end if

         do kk = 1, wrf%dom(id)%sls
            num_total = num_total + 1
            if(num_total <= indmax) then
!!$               indices(num_total) = get_wrf_index(ii,jj,kk,type_tslb,id)
               indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,kk,type_tslb)
               dist(num_total) = close_dist(i)
            end if
         enddo

         if(wrf%dom(id)%surf_obs ) then

            num_total = num_total + 1
            if(num_total <= indmax) then
!!$               indices(num_total) = get_wrf_index(ii,jj,k,type_u10,id)
               indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_u10)
               dist(num_total) = close_dist(i)
            end if

            num_total = num_total + 1
            if(num_total <= indmax) then
!!$               indices(num_total) = get_wrf_index(ii,jj,k,type_v10,id)
               indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_v10)
               dist(num_total) = close_dist(i)
            end if

            num_total = num_total + 1
            if(num_total <= indmax) then
!!$               indices(num_total) = get_wrf_index(ii,jj,k,type_t2,id)
               indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_t2)
               dist(num_total) = close_dist(i)
            end if

            num_total = num_total + 1
            if(num_total <= indmax) then
!!$               indices(num_total) = get_wrf_index(ii,jj,k,type_q2,id)
               indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_q2)
               dist(num_total) = close_dist(i)
            end if

            num_total = num_total + 1
            if(num_total <= indmax) then
!!$               indices(num_total) = get_wrf_index(ii,jj,k,type_ps,id)
               indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_ps)
               dist(num_total) = close_dist(i)
            end if

         end if

      end if

   enddo

! last -> p_pts (t + wrf%dom(id)%n_moist variables)

   do i = u_pts+v_pts+w_pts+1, u_pts+v_pts+w_pts+p_pts

      ii = lon_ind(i)
      jj = lat_ind(i)
      k = vert_ind(i)

      num_total = num_total + 1
      if(num_total <= indmax) then
!!$         indices(num_total) = get_wrf_index(ii,jj,k,type_t,id)
         indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_t)
         dist(num_total) = close_dist(i)
      end if
      if( wrf%dom(id)%n_moist >= 1) then
         num_total = num_total + 1
         if(num_total <= indmax) then
!!$            indices(num_total) = get_wrf_index(ii,jj,k,type_qv,id)
            indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_qv)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( wrf%dom(id)%n_moist >= 2) then
         num_total = num_total + 1
         if(num_total <= indmax) then
!!$            indices(num_total) = get_wrf_index(ii,jj,k,type_qc,id)
            indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_qc)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( wrf%dom(id)%n_moist >= 3) then
         num_total = num_total + 1
         if(num_total <= indmax) then
!!$            indices(num_total) = get_wrf_index(ii,jj,k,type_qr,id)
            indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_qr)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( wrf%dom(id)%n_moist >= 4) then
         num_total = num_total + 1
         if(num_total <= indmax) then
!!$            indices(num_total) = get_wrf_index(ii,jj,k,type_qi,id)
            indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_qi)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( wrf%dom(id)%n_moist >= 5) then
         num_total = num_total + 1
         if(num_total <= indmax) then
!!$            indices(num_total) = get_wrf_index(ii,jj,k,type_qs,id)
            indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_qs)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( wrf%dom(id)%n_moist >= 6) then
         num_total = num_total + 1
         if(num_total <= indmax) then
!!$            indices(num_total) = get_wrf_index(ii,jj,k,type_qg,id)
            indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_qg)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( wrf%dom(id)%n_moist == 7) then
         num_total = num_total + 1
         if(num_total <= indmax) then
!!$            indices(num_total) = get_wrf_index(ii,jj,k,type_qnice,id)
            indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_qnice)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( h_diab) then
         num_total = num_total + 1
         if(num_total <= indmax) then
!!$            indices(num_total) = get_wrf_index(ii,jj,k,type_hdiab,id)
            indices(num_total) = wrf%dom(id)%dart_ind(ii,jj,k,type_hdiab)
            dist(num_total) = close_dist(i)
         end if
      end if

   enddo

   deallocate(lon_ind, lat_ind, vert_ind, close_dist)

enddo

end subroutine model_get_close_states

!-----------------------------------------------------------------

function get_wrf_index( i,j,k,var_type,id )

integer, intent(in) :: i,j,k,var_type,id

integer :: get_wrf_index
integer :: in

character(len=129) :: errstring

do in = 1, wrf%dom(id)%number_of_wrf_variables
   if(var_type == wrf%dom(id)%var_type(in) ) then
      exit
   endif
enddo

if(i >= 1 .and. i <= wrf%dom(id)%var_size(1,in) .and. &
   j >= 1 .and. j <= wrf%dom(id)%var_size(2,in) .and. &
   k >= 1 .and. k <= wrf%dom(id)%var_size(3,in)) then

   get_wrf_index = wrf%dom(id)%dart_ind(i,j,k,var_type)

!!$   get_wrf_index = wrf%dom(id)%var_index(1,in)-1 +   &
!!$        i + wrf%dom(id)%var_size(1,in)*((j-1) + &
!!$        wrf%dom(id)%var_size(2,in)*(k-1))

else

  write(errstring,*)'Indices ',i,j,k,' exceed grid dimensions: ', &
       wrf%dom(id)%var_size(1,in), &
       wrf%dom(id)%var_size(2,in),wrf%dom(id)%var_size(3,in)
  call error_handler(E_ERR,'get_wrf_index', errstring, &
       source, revision, revdate)

endif

end function get_wrf_index

!-----------------------------------------------------------------

subroutine grid_close_states( o_loc, radius_in, id,  &
                              close_lon_ind, close_lat_ind, close_vert_ind, &
                              close_dist, u_pts, v_pts, p_pts, w_pts, x )

! Finds close state points from a particular grid for the WRF model

type(location_type), intent(in) :: o_loc
real(r8),            intent(in) :: radius_in
integer,             intent(in) :: id
integer,            intent(out) :: close_lon_ind(:), close_lat_ind(:), close_vert_ind(:)
real(r8),           intent(out) :: close_dist(:)
integer,            intent(out) :: u_pts, v_pts, p_pts, w_pts
real(r8),            intent(in) :: x(:)

real(r8)            :: radius
real(r8)            :: rad, radn, dxr, dyr, sdx, sdy, gdist
integer             :: i_closest, j_closest, ixmin, jymin, ixmax, jymax

integer             :: i, j, k, n, m, fis, is, ie, num

logical,  parameter :: debug = .false.  

if(debug) write(*,*) ' in grid_close_states '

! use meaningful units for radius -- convert radians to meters 
radius = radius_in*earth_radius*1000.0_r8
if(debug) write(*,*) ' radius in grid_close_states is ',radius

! Get index to closest lat and lon for this observation

n = wrf%dom(id)%we
m = wrf%dom(id)%sn

if(debug) write(*,*) 'obs location: ',get_location(o_loc)

rad = get_dist_wrf(1,1,1, type_t, o_loc, id, x)
i_closest = 1
j_closest = 1

! brute force search
do j=1,m
   do i=1,n
      radn = get_dist_wrf(i,j,1, type_t, o_loc, id, x)
      if( radn .lt. rad ) then
         rad = radn
         i_closest = i
         j_closest = j
      end if
      if(debug) write(*,*) i,j,radn,rad
   enddo
enddo

if(debug) write(*,*) ' closest grid point is ',i_closest,j_closest
if(debug) write(*,*) ' at distance ',rad

! define box edges for radius check
dxr = 1.0_r8 + radius/wrf%dom(id)%dx  !  radius in multiples of dx
dyr = 1.0_r8 + radius/wrf%dom(id)%dy  !  radius in multiples of dy

if(debug) write(*,*) ' dxr, dyr in grid_close_states ',dxr,dyr


j = j_closest
ixmin = max(1,i_closest - 1)
sdx   = 1.0_r8/wrf%dom(id)%mapfac_u(i_closest, j)
do while( sdx .lt. dxr )
   ixmin = max(1,ixmin - 1)
   sdx = sdx + 1.0_r8/wrf%dom(id)%mapfac_u(ixmin + 1, j)
   if(ixmin <= 1) sdx = 1.1_r8*dxr
enddo

ixmax = min(wrf%dom(id)%we,i_closest + 1)
sdx   = 1.0_r8/wrf%dom(id)%mapfac_u(ixmax, j)
do while( sdx .lt. dxr )
   ixmax = min(wrf%dom(id)%we,ixmax + 1)
   sdx = sdx + 1.0_r8/wrf%dom(id)%mapfac_u(ixmax, j)
   if(ixmax >= wrf%dom(id)%we) sdx = 1.1_r8*dxr
enddo

i = i_closest
jymin = max(1,j_closest - 1)
sdy   = 1.0_r8/wrf%dom(id)%mapfac_u(i_closest, jymin)
do while( sdy .lt. dyr )
   jymin = max(1,jymin - 1)
   sdy = sdy + 1.0_r8/wrf%dom(id)%mapfac_u(i, jymin + 1)
   if(jymin <= 1) sdy = 1.1_r8*dyr
enddo

jymax = min(wrf%dom(id)%sn,j_closest + 1)
sdy   = 1.0_r8/wrf%dom(id)%mapfac_u(i, jymax)
do while( sdy .lt. dyr )
   jymax = min(wrf%dom(id)%sn,jymax + 1)
   sdy = sdy + 1.0_r8/wrf%dom(id)%mapfac_u(i, jymax)
   if(jymax >= wrf%dom(id)%sn) sdy = 1.1_r8*dyr
enddo

if(debug) then
   write(*,*) ' ixmin, ixmax, jymin, jymax are '
   write(*,*) ixmin, ixmax, jymin, jymax
endif

!  we have bounding box, get and check distances.
!  first, convert radius back to radians

radius = radius_in
num = 0

if(debug) write(*,*) 'radius (radian) ',radius

! check distance for u points, expand box so that
! we don't leave possible points out of check

k = 1
do j = jymin, jymax
   do i = max(1,ixmin-1), ixmax+1

      gdist = get_dist_wrf(i,j,k, type_u, o_loc, id, x)
      if ( gdist <= radius ) then
         num = num + 1
         if(debug) write(*,*) ' u pt ',num,i,j,gdist
         close_lon_ind(num) = i
         close_lat_ind(num) = j
         close_vert_ind(num) = k
         close_dist(num) = gdist
      end if

   enddo
enddo

u_pts = num

if(horiz_dist_only) then

   fis = 1
   is = num + 1
   ie = num + u_pts
   do k = 2, wrf%dom(id)%bt
      close_lon_ind(is:ie) = close_lon_ind(fis:num)
      close_lat_ind(is:ie) = close_lat_ind(fis:num)
      close_vert_ind(is:ie) = k
      close_dist(is:ie) = close_dist(fis:num)
      is = is + u_pts
      ie = ie + u_pts
   enddo

   num = num + u_pts*(wrf%dom(id)%bt - 1)

else

   do k = 2, wrf%dom(id)%bt
      do j = jymin, jymax
         do i = max(1,ixmin-1), ixmax+1

            gdist = get_dist_wrf(i,j,k, type_u, o_loc, id, x)
            if ( gdist <= radius ) then
               num = num + 1
               if(debug) write(*,*) ' u pt ',num,i,j,k, id,gdist
               close_lon_ind(num) = i
               close_lat_ind(num) = j
               close_vert_ind(num) = k
               close_dist(num) = gdist
            end if

         enddo
      enddo
   enddo

endif

u_pts = num

k = 1
do j = max(1,jymin-1), jymax+1
   do i = ixmin, ixmax

      gdist = get_dist_wrf(i,j,k, type_v, o_loc, id, x)
      if ( gdist <= radius ) then
         num = num + 1
         close_lon_ind(num) = i
         close_lat_ind(num) = j
         close_vert_ind(num) = k
         close_dist(num) = gdist
         if(debug) write(*,*) ' v pt ',num,i,j,gdist
      end if
   enddo
enddo

v_pts = num - u_pts

if(horiz_dist_only) then

   fis = u_pts + 1
   is = num + 1
   ie = num + v_pts
   do k = 2, wrf%dom(id)%bt
      close_lon_ind(is:ie) = close_lon_ind(fis:num)
      close_lat_ind(is:ie) = close_lat_ind(fis:num)
      close_vert_ind(is:ie) = k
      close_dist(is:ie) = close_dist(fis:num)
      is = is + v_pts
      ie = ie + v_pts
   enddo

   num = num + v_pts*(wrf%dom(id)%bt - 1)

else

   do k = 2, wrf%dom(id)%bt
      do j = max(1,jymin-1), jymax+1
         do i = ixmin, ixmax

            gdist = get_dist_wrf(i,j,k, type_v, o_loc, id, x)
            if ( gdist <= radius ) then
               num = num + 1
               close_lon_ind(num) = i
               close_lat_ind(num) = j
               close_vert_ind(num) = k
               close_dist(num) = gdist
               if(debug) write(*,*) ' v pt ',num,i,j,gdist
            end if
         enddo
      enddo
   enddo

endif

v_pts = num - u_pts

k = 1
do j = jymin, jymax
   do i = ixmin, ixmax

      gdist = get_dist_wrf(i,j,k, type_w, o_loc, id, x)
      if ( gdist <= radius ) then
         num = num + 1
         close_lon_ind(num) = i
         close_lat_ind(num) = j
         close_vert_ind(num) = k
         close_dist(num) = gdist
         if(debug) write(*,*) ' w pt ',num,i,j,gdist
      end if
   enddo
enddo

w_pts = num - u_pts - v_pts

if(horiz_dist_only) then

   fis = u_pts + v_pts + 1
   is = num + 1
   ie = num + w_pts
   do k = 2, wrf%dom(id)%bts
      close_lon_ind(is:ie) = close_lon_ind(fis:num)
      close_lat_ind(is:ie) = close_lat_ind(fis:num)
      close_vert_ind(is:ie) = k
      close_dist(is:ie) = close_dist(fis:num)
      is = is + w_pts
      ie = ie + w_pts
   enddo

   num = num + w_pts*(wrf%dom(id)%bts - 1)

else

   do k = 2, wrf%dom(id)%bts
      do j = jymin, jymax
         do i = ixmin, ixmax

            gdist = get_dist_wrf(i,j,k, type_w, o_loc, id, x)
            if ( gdist <= radius ) then
               num = num + 1
               close_lon_ind(num) = i
               close_lat_ind(num) = j
               close_vert_ind(num) = k
               close_dist(num) = gdist
               if(debug) write(*,*) ' w pt ',num,i,j,gdist
            end if
         enddo
      enddo
   enddo

endif

w_pts = num - u_pts - v_pts

k = 1
do j = jymin, jymax
   do i = ixmin, ixmax

      gdist = get_dist_wrf(i,j,k, type_t, o_loc, id, x)
      if ( gdist <= radius ) then
         num = num + 1
         close_lon_ind(num) = i
         close_lat_ind(num) = j
         close_vert_ind(num) = k
         close_dist(num) = gdist
         if(debug) write(*,*) ' p pt ',num,i,j,gdist
      end if

   enddo
enddo

p_pts = num - u_pts - v_pts - w_pts

if(horiz_dist_only) then

   fis = u_pts + v_pts + w_pts + 1
   is = num + 1
   ie = num + p_pts
   do k = 2, wrf%dom(id)%bt
      close_lon_ind(is:ie) = close_lon_ind(fis:num)
      close_lat_ind(is:ie) = close_lat_ind(fis:num)
      close_vert_ind(is:ie) = k
      close_dist(is:ie) = close_dist(fis:num)
      is = is + p_pts
      ie = ie + p_pts
   enddo

   num = num + p_pts*(wrf%dom(id)%bt - 1)

else

   do k = 2, wrf%dom(id)%bt
      do j = jymin, jymax
         do i = ixmin, ixmax

            gdist = get_dist_wrf(i,j,k, type_t, o_loc, id, x)
            if ( gdist <= radius ) then
               num = num + 1
               close_lon_ind(num) = i
               close_lat_ind(num) = j
               close_vert_ind(num) = k
               close_dist(num) = gdist
               if(debug) write(*,*) ' p pt ',num,i,j,gdist
            end if

         enddo
      enddo
   enddo

endif

p_pts = num - u_pts - v_pts - w_pts

end subroutine grid_close_states

!***********************************************************************

function get_dist_wrf( i, j, k, var_type, o_loc, id, x )

real(r8)                        :: get_dist_wrf 
type(location_type), intent(in) :: o_loc
integer,             intent(in) :: i,j,k,var_type, id
real(r8),            intent(in) :: x(:)

type(location_type) :: s_loc, obs_loc
real(r8)            :: long, lat, vloc, xyz_loc(3)
integer             :: which_vert

character(len=129) :: errstring

! get distance for input var_type

call get_wrf_horizontal_location( i, j, var_type, id, long, lat )

! We will set which_vert based on the incoming location type
! This will ensure compatibility for measuring 'horizontal'
! and cannot think of anything better for vertical.
!
! For (near) surface observations,
! use height as the vertical coordinate for localization.

xyz_loc = get_location(o_loc)
which_vert = nint(query_location(o_loc,'which_vert'))

if(horiz_dist_only .or. vert_is_undef(o_loc)) then

   vloc = missing_r8

else

   if(vert_is_level(o_loc)) then

      if( (var_type == type_w ) .or. (var_type == type_gz) ) then
         vloc = real(k) - 0.5_r8
      else
         vloc = real(k)
      endif

   elseif(vert_is_pressure(o_loc)) then

      vloc = model_pressure(i,j,k,id,var_type,x)

   elseif(vert_is_height(o_loc) .or. vert_is_surface(o_loc)) then

      vloc = model_height(i,j,k,id,var_type,x)
      which_vert = 3

   else
      write(errstring, *) 'Which_vert = ',which_vert, ' is not allowed.'
      call error_handler(E_ERR,'get_dist_wrf', errstring, source, revision, revdate)
   endif

   if (vloc == missing_r8 ) then
      print*, i, j, k, id, var_type, which_vert
      write(errstring, *) 'Unable to define vloc.'
      call error_handler(E_ERR,'get_dist_wrf', errstring, source, revision, revdate)
   endif

endif

s_loc        = set_location(long, lat, vloc, which_vert)
obs_loc      = set_location(xyz_loc(1), xyz_loc(2), xyz_loc(3), which_vert)
get_dist_wrf = get_dist(s_loc,obs_loc)

end function get_dist_wrf



subroutine get_wrf_horizontal_location( i, j, var_type, id, long, lat )

integer,  intent(in)  :: i,j,var_type, id
real(r8), intent(out) :: long, lat

! find lat and long, must
! correct for possible u or v staggering in x, y

if(var_type == type_u) then

   if (i == 1) then
      long = wrf%dom(id)%longitude(1,j) - &
           0.5_r8*(wrf%dom(id)%longitude(2,j)-wrf%dom(id)%longitude(1,j))
      if ( abs(wrf%dom(id)%longitude(2,j) - wrf%dom(id)%longitude(1,j)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat = wrf%dom(id)%latitude(1,j) - &
           0.5_r8*(wrf%dom(id)%latitude(2,j)-wrf%dom(id)%latitude(1,j))
   else if (i == wrf%dom(id)%wes) then
      long = wrf%dom(id)%longitude(i-1,j) + &
           0.5_r8*(wrf%dom(id)%longitude(i-1,j)-wrf%dom(id)%longitude(i-2,j))
      if ( abs(wrf%dom(id)%longitude(i-1,j) - wrf%dom(id)%longitude(i-2,j)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat = wrf%dom(id)%latitude(i-1,j) + &
           0.5_r8*(wrf%dom(id)%latitude(i-1,j)-wrf%dom(id)%latitude(i-2,j))
   else
      long = 0.5_r8*(wrf%dom(id)%longitude(i,j)+wrf%dom(id)%longitude(i-1,j))
      if ( abs(wrf%dom(id)%longitude(i,j) - wrf%dom(id)%longitude(i-1,j)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat = 0.5_r8*(wrf%dom(id)%latitude(i,j) +wrf%dom(id)%latitude(i-1,j))
   end if

else if( var_type == type_v) then

   if (j == 1) then
      long = wrf%dom(id)%longitude(i,1) - &
           0.5_r8*(wrf%dom(id)%longitude(i,2)-wrf%dom(id)%longitude(i,1))
      if ( abs(wrf%dom(id)%longitude(i,2) - wrf%dom(id)%longitude(i,1)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat = wrf%dom(id)%latitude(i,1) - &
           0.5_r8*(wrf%dom(id)%latitude(i,2)-wrf%dom(id)%latitude(i,1))
   else if (j == wrf%dom(id)%sns) then
      long = wrf%dom(id)%longitude(i,j-1) + &
           0.5_r8*(wrf%dom(id)%longitude(i,j-1)-wrf%dom(id)%longitude(i,j-2))
      if ( abs(wrf%dom(id)%longitude(i,j-1) - wrf%dom(id)%longitude(i,j-2)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat = wrf%dom(id)%latitude(i,j-1) + &
           0.5_r8*(wrf%dom(id)%latitude(i,j-1)-wrf%dom(id)%latitude(i,j-2))
   else
      long = 0.5_r8*(wrf%dom(id)%longitude(i,j)+wrf%dom(id)%longitude(i,j-1))
      if ( abs(wrf%dom(id)%longitude(i,j) - wrf%dom(id)%longitude(i,j-1)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat  = 0.5_r8*(wrf%dom(id)%latitude(i,j) +wrf%dom(id)%latitude(i,j-1))

   end if

else

   long = wrf%dom(id)%longitude(i,j)
   lat  = wrf%dom(id)%latitude(i,j)

end if

do while (long <   0.0_r8)
   long = long + 360.0_r8
end do
do while (long > 360.0_r8)
   long = long - 360.0_r8
end do

end subroutine get_wrf_horizontal_location



function nc_write_model_atts( ncFileID ) result (ierr)
!-----------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! A. Caya May 7 2003
! T. Hoar Mar 8 2004 writes prognostic flavor

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!-----------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarDimID, StateVarVarID, StateVarID, TimeDimID

integer, dimension(num_domains) :: weDimID, weStagDimID, snDimID, snStagDimID, &
     btDimID, btStagDimID, slSDimID

integer :: MemberDimID, DomDimID
integer :: DXVarID, DYVarID, TRUELAT1VarID, TRUELAT2VarID
integer :: CEN_LATVarID, CEN_LONVarID, MAP_PROJVarID

integer, dimension(num_domains) :: DNVarID, ZNUVarID, DNWVarID, phbVarID, &
     MubVarID, LonVarID, LatVarID, ilevVarID, XlandVarID, &
     MapFacMVarID, MapFacUVarID, MapFacVVarID, hgtVarID

integer :: var_id
integer :: i, id

character(len=129) :: errstring

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

character (len=1)     :: idom

!-----------------------------------------------------------------

ierr = 0     ! assume normal termination

!-----------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-----------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_Redef(ncFileID))

!-----------------------------------------------------------------
! We need the dimension ID for the number of copies 
!-----------------------------------------------------------------

call check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID))
call check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID))

if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)'Time Dimension ID ',TimeDimID, &
        ' must match Unlimited Dimension ID ',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', errstring, source, revision, revdate)
endif

!-----------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-----------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=wrf%model_size, dimid = StateVarDimID))

!-----------------------------------------------------------------
! Write Global Attributes 
!-----------------------------------------------------------------
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","WRF"))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source",source))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate))

! how about namelist input? might be nice to save ...

!-----------------------------------------------------------------
! Define the dimensions IDs
!-----------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="domain",           &
          len = num_domains,  dimid = DomDimID))

do id=1,num_domains
   write( idom , '(I1)') id
   call check(nf90_def_dim(ncid=ncFileID, name="west_east_d0"//idom,        &
        len = wrf%dom(id)%we,  dimid = weDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="west_east_stag_d0"//idom,   &
        len = wrf%dom(id)%wes, dimid = weStagDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="south_north_d0"//idom,      &
        len = wrf%dom(id)%sn,  dimid = snDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="south_north_stag_d0"//idom, &
        len = wrf%dom(id)%sns, dimid = snStagDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="bottom_top_d0"//idom,       &
        len = wrf%dom(id)%bt,  dimid = btDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="bottom_top_stag_d0"//idom,  &
        len = wrf%dom(id)%bts, dimid = btStagDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="soil_layers_stag_d0"//idom,  &
        len = wrf%dom(id)%sls, dimid = slSDimID(id)))
enddo

!-----------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Create the (empty) static variables and their attributes
! Commented block is from wrfinput
!-----------------------------------------------------------------

call check(nf90_def_var(ncFileID, name="DX", xtype=nf90_real, &
     dimids= DomDimID, varid=DXVarID) )
call check(nf90_put_att(ncFileID, DXVarID, "long_name", &
     "X HORIZONTAL RESOLUTION"))
call check(nf90_put_att(ncFileID, DXVarID, "units", "m"))

call check(nf90_def_var(ncFileID, name="DY", xtype=nf90_real, &
     dimids= DomDimID, varid=DYVarID) )
call check(nf90_put_att(ncFileID, DYVarID, "long_name", &
     "Y HORIZONTAL RESOLUTION"))
call check(nf90_put_att(ncFileID, DYVarID, "units", "m"))

call check(nf90_def_var(ncFileID, name="TRUELAT1", xtype=nf90_real, &
     dimids= DomDimID, varid=TRUELAT1VarID) )
call check(nf90_put_att(ncFileID, TRUELAT1VarID, "long_name", &
     "first standard parallel"))
call check(nf90_put_att(ncFileID, TRUELAT1VarID, "units", &
     "degrees, negative is south"))

call check(nf90_def_var(ncFileID, name="TRUELAT2", xtype=nf90_real, &
     dimids= DomDimID, varid=TRUELAT2VarID) )
call check(nf90_put_att(ncFileID, TRUELAT2VarID, "long_name", &
     "second standard parallel"))
call check(nf90_put_att(ncFileID, TRUELAT2VarID, "units", &
     "degrees, negative is south"))

call check(nf90_def_var(ncFileID, name="CEN_LAT", xtype=nf90_real, &
     dimids= DomDimID, varid=CEN_LATVarID) )
call check(nf90_put_att(ncFileID, CEN_LATVarID, "long_name", &
     "center latitude"))
call check(nf90_put_att(ncFileID, CEN_LATVarID, "units", &
     "degrees, negative is south"))

call check(nf90_def_var(ncFileID, name="CEN_LON", xtype=nf90_real, &
     dimids= DomDimID, varid=CEN_LONVarID) )
call check(nf90_put_att(ncFileID, CEN_LONVarID, "long_name", &
     "central longitude"))
call check(nf90_put_att(ncFileID, CEN_LONVarID, "units", &
     "degrees, negative is west"))

call check(nf90_def_var(ncFileID, name="MAP_PROJ", xtype=nf90_real, &
     dimids= DomDimID, varid=MAP_PROJVarID) )
call check(nf90_put_att(ncFileID, MAP_PROJVarID, "long_name", &
     "domain map projection"))
call check(nf90_put_att(ncFileID, MAP_PROJVarID, "units", &
     "0=none, 1=Lambert, 2=polar, 3=Mercator"))

do id=1,num_domains
   write( idom , '(I1)') id

   call check(nf90_def_var(ncFileID, name="DN_d0"//idom, xtype=nf90_real, &
        dimids= btDimID(id), varid=DNVarID(id)) )
   call check(nf90_put_att(ncFileID, DNVarID(id), "long_name", &
        "dn values on half (mass) levels"))
   call check(nf90_put_att(ncFileID, DNVarID(id), "units", "dimensionless"))

   call check(nf90_def_var(ncFileID, name="ZNU_d0"//idom, xtype=nf90_real, &
        dimids= btDimID(id), varid=ZNUVarID(id)) )
   call check(nf90_put_att(ncFileID, ZNUVarID(id), "long_name", &
        "eta values on half (mass) levels"))
   call check(nf90_put_att(ncFileID, ZNUVarID(id), "units", "dimensionless"))

   call check(nf90_def_var(ncFileID, name="DNW_d0"//idom, xtype=nf90_real, &
        dimids= btDimID(id), varid=DNWVarID(id)) )
   call check(nf90_put_att(ncFileID, DNWVarID(id), "long_name", &
        "dn values on full (w) levels"))
   call check(nf90_put_att(ncFileID, DNWVarID(id), "units", "dimensionless"))

!
!    float MUB(Time, south_north, west_east) ;
!            MUB:FieldType = 104 ;
!            MUB:MemoryOrder = "XY " ;
!            MUB:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="MUB_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snDimID(id) /), varid=MubVarID(id)) )
   call check(nf90_put_att(ncFileID, MubVarID(id), "long_name", &
        "base state dry air mass in column"))
   call check(nf90_put_att(ncFileID, MubVarID(id), "units", "Pa"))

! Longitudes
!      float XLONG(Time, south_north, west_east) ;
!         XLONG:FieldType = 104 ;
!         XLONG:MemoryOrder = "XY " ;
!         XLONG:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="XLON_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snDimID(id) /), varid=LonVarID(id)) )
   call check(nf90_put_att(ncFileID, LonVarID(id), "long_name", "longitude"))
   call check(nf90_put_att(ncFileID, LonVarID(id), "units", "degrees_east"))
   call check(nf90_put_att(ncFileID, LonVarID(id), "valid_range", &
        (/ -180.0_r8, 180.0_r8 /)))
   call check(nf90_put_att(ncFileID, LonVarID(id), "description", &
        "LONGITUDE, WEST IS NEGATIVE"))

! Latitudes
!      float XLAT(Time, south_north, west_east) ;
!         XLAT:FieldType = 104 ;
!         XLAT:MemoryOrder = "XY " ;
!         XLAT:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="XLAT_d0"//idom, xtype=nf90_real, &
        dimids=(/ weDimID(id), snDimID(id) /), varid=LatVarID(id)) ) 
   call check(nf90_put_att(ncFileID, LatVarID(id), "long_name", "latitude"))
   call check(nf90_put_att(ncFileID, LatVarID(id), "units", "degrees_north"))
   call check(nf90_put_att(ncFileID, LatVarID(id), "valid_range", &
        (/ -90.0_r8, 90.0_r8 /)))
   call check(nf90_put_att(ncFileID, LatVarID(id), "description", &
        "LATITUDE, SOUTH IS NEGATIVE"))

! grid levels
   call check(nf90_def_var(ncFileID, name="level_d0"//idom, xtype=nf90_short, &
        dimids=btDimID(id), varid=ilevVarID(id)) )
   call check(nf90_put_att(ncFileID, ilevVarID(id), "long_name", &
        "placeholder for level"))
   call check(nf90_put_att(ncFileID, ilevVarID(id), "units", &
        "at this point, indexical"))

! Land Mask
!    float XLAND(Time, south_north, west_east) ;
!            XLAND:FieldType = 104 ;
!            XLAND:MemoryOrder = "XY " ;
!            XLAND:units = "NA" ;
!            XLAND:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="XLAND_d0"//idom, xtype=nf90_short, &
        dimids= (/ weDimID(id), snDimID(id) /), varid=XlandVarID(id)) )
   call check(nf90_put_att(ncFileID, XlandVarID(id), "long_name", "land mask"))
   call check(nf90_put_att(ncFileID, XlandVarID(id), "units", "NA"))
   call check(nf90_put_att(ncFileID, XlandVarID(id), "valid_range", (/ 1, 2 /)))
   call check(nf90_put_att(ncFileID, XlandVarID(id), "description", &
        "1 = LAND, 2 = WATER"))

! Map Scale Factor on m-grid
!    float MAPFAC_M(Time, south_north, west_east) ;
!            MAPFAC_M:FieldType = 104 ;
!            MAPFAC_M:MemoryOrder = "XY " ;
!            MAPFAC_M:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="MAPFAC_M_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snDimID(id) /), varid=MapFacMVarID(id)) )
   call check(nf90_put_att(ncFileID, MapFacMVarID(id), "long_name", &
        "Map scale factor on mass grid"))
   call check(nf90_put_att(ncFileID, MapFacMVarID(id), "units", "dimensionless"))

! Map Scale Factor on u-grid
!    float MAPFAC_U(Time, south_north, west_east_stag) ;
!            MAPFAC_U:FieldType = 104 ;
!            MAPFAC_U:MemoryOrder = "XY " ;
!            MAPFAC_U:stagger = "X" ;
   call check(nf90_def_var(ncFileID, name="MAPFAC_U_d0"//idom, xtype=nf90_real, &
        dimids= (/ weStagDimID(id), snDimID(id) /), varid=MapFacUVarID(id)) )
   call check(nf90_put_att(ncFileID, MapFacUVarID(id), "long_name", &
        "Map scale factor on u-grid"))
   call check(nf90_put_att(ncFileID, MapFacUVarID(id), "units", "dimensionless"))

! Map Scale Factor on v-grid
!    float MAPFAC_V(Time, south_north_stag, west_east) ;
!            MAPFAC_V:FieldType = 104 ;
!            MAPFAC_V:MemoryOrder = "XY " ;
!            MAPFAC_V:stagger = "Y" ;
   call check(nf90_def_var(ncFileID, name="MAPFAC_V_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snStagDimID(id) /), varid=MapFacVVarID(id)) )
   call check(nf90_put_att(ncFileID, MapFacVVarID(id), "long_name", &
        "Map scale factor on v-grid"))
   call check(nf90_put_att(ncFileID, MapFacVVarID(id), "units", "dimensionless"))

! PHB
!    float PHB(Time, bottom_top_stag, south_north, west_east) ;
!            PHB:FieldType = 104 ;
!            PHB:MemoryOrder = "XYZ" ;
!            PHB:stagger = "Z" ;
   call check(nf90_def_var(ncFileID, name="PHB_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snDimID(id), btStagDimID(id) /), varid=phbVarId(id)) )
   call check(nf90_put_att(ncFileID, phbVarId(id), "long_name", &
        "base-state geopotential"))
   call check(nf90_put_att(ncFileID, phbVarId(id), "units", "m^2/s^2"))
   call check(nf90_put_att(ncFileID, phbVarId(id), "units_long_name", "m{2} s{-2}"))

   call check(nf90_def_var(ncFileID, name="HGT_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snDimID(id) /), varid=hgtVarId(id)) )
   call check(nf90_put_att(ncFileID, hgtVarId(id), "long_name", "Terrain Height"))
   call check(nf90_put_att(ncFileID, hgtVarId(id), "units", "m"))
   call check(nf90_put_att(ncFileID, hgtVarId(id), "units_long_name", "meters"))

enddo

if ( output_state_vector ) then

   !-----------------------------------------------------------------
   ! Create attributes for the state vector 
   !-----------------------------------------------------------------

   ! Define the state vector coordinate variable

   call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID))

   call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", &
        "State Variable ID"))
   call check(nf90_put_att(ncFileID, StateVarVarID, "units", &
        "indexical") )
   call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", &
        (/ 1, wrf%model_size /)))

   ! Define the actual state vector

   call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
              varid=StateVarID))
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", &
        "model state or fcopy"))
   call check(nf90_put_att(ncFileID, StateVarId, "U_units","m/s"))
   call check(nf90_put_att(ncFileID, StateVarId, "V_units","m/s"))
   call check(nf90_put_att(ncFileID, StateVarId, "W_units","m/s"))
   call check(nf90_put_att(ncFileID, StateVarId, "GZ_units","m2/s2"))
   call check(nf90_put_att(ncFileID, StateVarId, "T_units","K"))
   call check(nf90_put_att(ncFileID, StateVarId, "MU_units","Pa"))
   call check(nf90_put_att(ncFileID, StateVarId, "TSLB_units","K"))
   call check(nf90_put_att(ncFileID, StateVarId, "TSK_units","K"))
   if( wrf%dom(num_domains)%n_moist >= 1) then
      call check(nf90_put_att(ncFileID, StateVarId, "QV_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist >= 2) then
      call check(nf90_put_att(ncFileID, StateVarId, "QC_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist >= 3) then
      call check(nf90_put_att(ncFileID, StateVarId, "QR_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist >= 4) then
      call check(nf90_put_att(ncFileID, StateVarId, "QI_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist >= 5) then
      call check(nf90_put_att(ncFileID, StateVarId, "QS_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist >= 6) then
      call check(nf90_put_att(ncFileID, StateVarId, "QG_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist == 7) then
      call check(nf90_put_att(ncFileID, StateVarId, "QNICE_units","kg-1"))
   endif
   if(wrf%dom(num_domains)%surf_obs ) then
      call check(nf90_put_att(ncFileID, StateVarId, "U10_units","m/s"))
      call check(nf90_put_att(ncFileID, StateVarId, "V10_units","m/s"))
      call check(nf90_put_att(ncFileID, StateVarId, "T2_units","K"))
      call check(nf90_put_att(ncFileID, StateVarId, "Q2_units","kg/kg"))
      call check(nf90_put_att(ncFileID, StateVarId, "PS_units","Pa"))
   endif

   ! Leave define mode so we can actually fill the variables.

   call check(nf90_enddef(ncfileID))

   call check(nf90_put_var(ncFileID, StateVarVarID, &
        (/ (i,i=1,wrf%model_size) /) ))

else

do id=1,num_domains
   write( idom , '(I1)') id

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and their attributes
   !----------------------------------------------------------------------------

   !      float U(Time, bottom_top, south_north, west_east_stag) ;
   !         U:FieldType = 104 ;
   !         U:MemoryOrder = "XYZ" ;
   !         U:stagger = "X" ;
   call check(nf90_def_var(ncid=ncFileID, name="U_d0"//idom, xtype=nf90_real, &
         dimids = (/ weStagDimID(id), snDimId(id), btDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "x-wind component"))
   call check(nf90_put_att(ncFileID, var_id, "units", "m/s"))
   call check(nf90_put_att(ncFileID, var_id, "units_long_name", "m s{-1}"))


   !      float V(Time, bottom_top, south_north_stag, west_east) ;
   !         V:FieldType = 104 ;
   !         V:MemoryOrder = "XYZ" ;
   !         V:stagger = "Y" ;
   call check(nf90_def_var(ncid=ncFileID, name="V_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snStagDimID(id), btDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "y-wind component"))
   call check(nf90_put_att(ncFileID, var_id, "units", "m/s"))
   call check(nf90_put_att(ncFileID, var_id, "units_long_name", "m s{-1}"))


   !      float W(Time, bottom_top_stag, south_north, west_east) ;
   !         W:FieldType = 104 ;
   !         W:MemoryOrder = "XYZ" ;
   !         W:stagger = "Z" ;
   call check(nf90_def_var(ncid=ncFileID, name="W_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snDimID(id), btStagDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "z-wind component"))
   call check(nf90_put_att(ncFileID, var_id, "units", "m/s"))
   call check(nf90_put_att(ncFileID, var_id, "units_long_name", "m s{-1}"))


   !      float PH(Time, bottom_top_stag, south_north, west_east) ;
   !         PH:FieldType = 104 ;
   !         PH:MemoryOrder = "XYZ" ;
   !         PH:stagger = "Z" ;
   call check(nf90_def_var(ncid=ncFileID, name="PH_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snDimID(id), btStagDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", &
        "perturbation geopotential"))
   call check(nf90_put_att(ncFileID, var_id, "units", "m^2/s^2"))
   call check(nf90_put_att(ncFileID, var_id, "units_long_name", "m{2} s{-2}"))


   !      float T(Time, bottom_top, south_north, west_east) ;
   !         T:FieldType = 104 ;
   !         T:MemoryOrder = "XYZ" ;
   !         T:units = "K" ;
   !         T:stagger = "" ;
   call check(nf90_def_var(ncid=ncFileID, name="T_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "temperature"))
   call check(nf90_put_att(ncFileID, var_id, "units", "K"))
   call check(nf90_put_att(ncFileId, var_id, "description", &
        "perturbation potential temperature (theta-t0)"))


   !      float MU(Time, south_north, west_east) ;
   !         MU:FieldType = 104 ;
   !         MU:MemoryOrder = "XY " ;
   !         MU:description = "perturbation dry air mass in column" ;
   !         MU:units = "pascals" ;
   !         MU:stagger = "" ;
   call check(nf90_def_var(ncid=ncFileID, name="MU_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "mu field"))
   call check(nf90_put_att(ncFileID, var_id, "units", "pascals"))
   call check(nf90_put_att(ncFileId, var_id, "description", &
        "perturbation dry air mass in column"))


   !      float TSLB(Time, soil_layers_stag, south_north, west_east) ;
   !         TSLB:FieldType = 104 ;
   !         TSLB:MemoryOrder = "XYZ" ;
   !         TSLB:description = "SOIL TEMPERATURE" ;
   !         TSLB:units = "K" ;
   !         TSLB:stagger = "Z" ;
   call check(nf90_def_var(ncid=ncFileID, name="TSLB_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snDimID(id), slSDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "soil temperature"))
   call check(nf90_put_att(ncFileID, var_id, "units", "K"))
   call check(nf90_put_att(ncFileId, var_id, "description", &
        "SOIL TEMPERATURE"))


   !      float TSK(Time, south_north, west_east) ;
   !         TSK:FieldType = 104 ;
   !         TSK:MemoryOrder = "XY " ;
   !         TSK:description = "SURFACE SKIN TEMPERATURE" ;
   !         TSK:units = "K" ;
   !         TSK:stagger = "" ;
   call check(nf90_def_var(ncid=ncFileID, name="TSK_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "tsk field"))
   call check(nf90_put_att(ncFileID, var_id, "units", "K"))
   call check(nf90_put_att(ncFileId, var_id, "description", &
        "SURFACE SKIN TEMPERATURE"))


   !      float QVAPOR(Time, bottom_top, south_north, west_east) ;
   !         QVAPOR:FieldType = 104 ;
   !         QVAPOR:MemoryOrder = "XYZ" ;
   !         QVAPOR:stagger = "" ;
   if( wrf%dom(id)%n_moist >= 1) then
      call check(nf90_def_var(ncid=ncFileID, name="QVAPOR_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileId, var_id, "description", &
           "Water vapor mixing ratio"))
   endif


   !      float QCLOUD(Time, bottom_top, south_north, west_east) ;
   !         QCLOUD:FieldType = 104 ;
   !         QCLOUD:MemoryOrder = "XYZ" ;
   !         QCLOUD:stagger = "" ;
   if( wrf%dom(id)%n_moist >= 2) then
      call check(nf90_def_var(ncid=ncFileID, name="QCLOUD_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, "description", &
           "Cloud water mixing ratio"))
   endif


   !      float QRAIN(Time, bottom_top, south_north, west_east) ;
   !         QRAIN:FieldType = 104 ;
   !         QRAIN:MemoryOrder = "XYZ" ;
   !         QRAIN:stagger = "" ;
   if( wrf%dom(id)%n_moist >= 3) then
      call check(nf90_def_var(ncid=ncFileID, name="QRAIN_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, "description", &
           "Rain water mixing ratio"))
   endif

   if( wrf%dom(id)%n_moist >= 4) then
      call check(nf90_def_var(ncid=ncFileID, name="QICE_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, "description", &
           "Ice mixing ratio"))
   endif

   if( wrf%dom(id)%n_moist >= 5) then
      call check(nf90_def_var(ncid=ncFileID, name="QSNOW_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, &
           "description", "Snow mixing ratio"))
   endif

   if( wrf%dom(id)%n_moist >= 6) then
      call check(nf90_def_var(ncid=ncFileID, name="QGRAUP_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, "description", &
           "Graupel mixing ratio"))
   endif

   if( wrf%dom(id)%n_moist == 7) then
      call check(nf90_def_var(ncid=ncFileID, name="QNICE_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg-1"))
      call check(nf90_put_att(ncFileID, var_id, "description", &
           "Ice Number concentration"))
   endif

   if(wrf%dom(id)%surf_obs ) then

      call check(nf90_def_var(ncid=ncFileID, name="U10_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), MemberDimID, unlimitedDimID /), &
           varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "m/s"))
      call check(nf90_put_att(ncFileID, var_id, "description", "U at 10 m"))

      call check(nf90_def_var(ncid=ncFileID, name="V10_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), MemberDimID, unlimitedDimID /), &
           varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "m/s"))
      call check(nf90_put_att(ncFileID, var_id, "description", "V at 10 m"))

      call check(nf90_def_var(ncid=ncFileID, name="T2_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), MemberDimID, unlimitedDimID /), &
           varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "K"))
      call check(nf90_put_att(ncFileID, var_id, "description", "TEMP at 2 m"))

      call check(nf90_def_var(ncid=ncFileID, name="Q2_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), MemberDimID, unlimitedDimID /), &
           varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, "description", "QV at 2 m"))

      call check(nf90_def_var(ncid=ncFileID, name="PS_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), MemberDimID, unlimitedDimID /), &
           varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "Pa"))
      call check(nf90_put_att(ncFileID, var_id, "description", "Total surface pressure"))

   end if

enddo

endif

!-----------------------------------------------------------------
! Fill the variables we can
!-----------------------------------------------------------------
call check(nf90_enddef(ncfileID))

call check(nf90_put_var(ncFileID,       DXVarID, wrf%dom(1:num_domains)%dx        ))
call check(nf90_put_var(ncFileID,       DYVarID, wrf%dom(1:num_domains)%dy        ))
call check(nf90_put_var(ncFileID, TRUELAT1VarID, wrf%dom(1:num_domains)%proj%truelat1  ))
call check(nf90_put_var(ncFileID, TRUELAT2VarID, wrf%dom(1:num_domains)%proj%truelat2  ))
call check(nf90_put_var(ncFileID,  CEN_LATVarID, wrf%dom(1:num_domains)%cen_lat   ))
call check(nf90_put_var(ncFileID,  CEN_LONVarID, wrf%dom(1:num_domains)%cen_lon   ))
call check(nf90_put_var(ncFileID, MAP_PROJVarID, wrf%dom(1:num_domains)%map_proj  ))

do id=1,num_domains

! defining grid levels
   call check(nf90_put_var(ncFileID,       DNVarID(id), wrf%dom(id)%dn        ))
   call check(nf90_put_var(ncFileID,      ZNUVarID(id), wrf%dom(id)%znu       ))
   call check(nf90_put_var(ncFileID,      DNWVarID(id), wrf%dom(id)%dnw       ))

! defining horizontal
   call check(nf90_put_var(ncFileID,      mubVarID(id), wrf%dom(id)%mub       ))
   call check(nf90_put_var(ncFileID,      LonVarID(id), wrf%dom(id)%longitude ))
   call check(nf90_put_var(ncFileID,      LatVarID(id), wrf%dom(id)%latitude  ))
   call check(nf90_put_var(ncFileID,     ilevVarID(id), (/ (i,i=1,wrf%dom(id)%bt) /) ))
   call check(nf90_put_var(ncFileID,    XlandVarID(id), wrf%dom(id)%land      ))
   call check(nf90_put_var(ncFileID,  MapFacMVarID(id), wrf%dom(id)%mapfac_m  ))
   call check(nf90_put_var(ncFileID,  MapFacUVarID(id), wrf%dom(id)%mapfac_u  ))
   call check(nf90_put_var(ncFileID,  MapFacVVarID(id), wrf%dom(id)%mapfac_v  ))
   call check(nf90_put_var(ncFileID,      phbVarID(id), wrf%dom(id)%phb       ))
   call check(nf90_put_var(ncFileID,      hgtVarID(id), wrf%dom(id)%hgt       ))

enddo

!-----------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-----------------------------------------------------------------

call check(nf90_sync(ncFileID))

write (*,*)'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...'

contains

  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(istatus)
    integer, intent ( in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_atts', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!-----------------------------------------------------------------
! Writes the model-specific variables to a netCDF file
! TJH 25 June 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.


integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!-----------------------------------------------------------------

logical, parameter :: debug = .false.  
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, VarID, id
integer :: i,j
real(r8), allocatable, dimension(:,:)   :: temp2d
real(r8), allocatable, dimension(:,:,:) :: temp3d
character(len=10) :: varname
character(len=1) :: idom

ierr = 0     ! assume normal termination

!-----------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! then get all the Variable ID's we need.
!-----------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

if ( output_state_vector ) then

   call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)))

else

j = 0

do id=1,num_domains

   write( idom , '(I1)') id

   !----------------------------------------------------------------------------
   ! Fill the variables, the order is CRITICAL  ...   U,V,W,GZ,T,MU,TSLB,TSK,QV,QC,QR,...
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   varname = 'U_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%wes * wrf%dom(id)%sn * wrf%dom(id)%bt - 1 
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%wes,wrf%dom(id)%sn,wrf%dom(id)%bt 
   allocate ( temp3d(wrf%dom(id)%wes, wrf%dom(id)%sn, wrf%dom(id)%bt) )
   temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%wes, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'V_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sns * wrf%dom(id)%bt - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sns,wrf%dom(id)%bt
   allocate ( temp3d(wrf%dom(id)%we, wrf%dom(id)%sns, wrf%dom(id)%bt) )
   temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sns, wrf%dom(id)%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'W_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bts - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bts
   allocate ( temp3d(wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bts) )
   temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bts /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))


   !----------------------------------------------------------------------------
   varname = 'PH_d0'//idom       ! AKA "GZ"
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bts - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bts
   temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bts /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'T_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
   allocate ( temp3d(wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt) )
   temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'MU_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
   allocate ( temp2d(wrf%dom(id)%we, wrf%dom(id)%sn) )
   temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp2d, &
                            start=(/ 1, 1, copyindex, timeindex /) ))


   !----------------------------------------------------------------------------
   varname = 'TSLB_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%sls - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%sls
   allocate ( temp3d(wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%sls) )
   temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%sls /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'TSK_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
   temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp2d, &
                            start=(/ 1, 1, copyindex, timeindex /) ))


   if( wrf%dom(id)%n_moist >= 1) then
      !----------------------------------------------------------------------------
      varname = 'QVAPOR_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      allocate ( temp3d(wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt) )
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif


   if( wrf%dom(id)%n_moist >= 2) then
      !----------------------------------------------------------------------------
      varname = 'QCLOUD_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif


   if( wrf%dom(id)%n_moist >= 3) then
      !----------------------------------------------------------------------------
      varname = 'QRAIN_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif


   if( wrf%dom(id)%n_moist >= 4) then
      !----------------------------------------------------------------------------
      varname = 'QICE_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif
   if( wrf%dom(id)%n_moist >= 5) then
      !----------------------------------------------------------------------------
      varname = 'QSNOW_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif
   if( wrf%dom(id)%n_moist >= 6) then
      !----------------------------------------------------------------------------
      varname = 'QGRAUP_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif
   if( wrf%dom(id)%n_moist == 7) then
      !----------------------------------------------------------------------------
      varname = 'QNICE_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif

   deallocate(temp3d)

   if ( wrf%dom(id)%n_moist > 7 ) then
      write(*,'(''wrf%dom(id)%n_moist = '',i3)')wrf%dom(id)%n_moist
      call error_handler(E_ERR,'nc_write_model_vars', &
               'num_moist_vars is too large.', source, revision, revdate)
   endif

   if(wrf%dom(id)%surf_obs ) then

      !----------------------------------------------------------------------------
      varname = 'U10_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
      temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp2d, &
           start=(/ 1, 1, copyindex, timeindex /) ))

      !----------------------------------------------------------------------------
      varname = 'V10_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
      temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp2d, &
           start=(/ 1, 1, copyindex, timeindex /) ))

      !----------------------------------------------------------------------------
      varname = 'T2_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
      temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp2d, &
           start=(/ 1, 1, copyindex, timeindex /) ))

      !----------------------------------------------------------------------------
      varname = 'Q2_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
      temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp2d, &
           start=(/ 1, 1, copyindex, timeindex /) ))

      !----------------------------------------------------------------------------
      varname = 'PS_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
      temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp2d, &
           start=(/ 1, 1, copyindex, timeindex /) ))

   endif

   deallocate(temp2d)

   if( h_diab) then
      j  = j + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt
   end if

enddo

endif

!-----------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-----------------------------------------------------------------

write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
write (*,*)'netCDF file is synched ...'

contains

  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(istatus)
    integer, intent ( in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_vars', &
         trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end function nc_write_model_vars

!-------------------------------

!  public stubs

!**********************************************

subroutine adv_1step(x, Time)

! Does single time-step advance with vector state as
! input and output.

  real(r8), intent(inout) :: x(:)

! Time is needed for more general models like this; need to add in to 
! low-order models
  type(time_type), intent(in) :: Time

end subroutine adv_1step

!**********************************************

subroutine end_model()
end subroutine end_model

!**********************************************

subroutine init_time(i_time)
! For now returns value of Time_init which is set in initialization routines.

  type(time_type), intent(out) :: i_time

!Where should initial time come from here?
! WARNING: CURRENTLY SET TO 0
  i_time = set_time(0, 0)

end subroutine init_time

!**********************************************

subroutine init_conditions(x)
! Reads in restart initial conditions and converts to vector

! Following changed to intent(inout) for ifc compiler;should be like this
  real(r8), intent(inout) :: x(:)

end subroutine init_conditions



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

subroutine pres_to_zk(pres, mdl_v, n3, zk)

! Calculate the model level "zk" on half (mass) levels,
! corresponding to pressure "pres".

  integer,  intent(in)  :: n3
  real(r8), intent(in)  :: pres
  real(r8), intent(in)  :: mdl_v(0:n3)
  real(r8), intent(out) :: zk

  integer  :: k

  zk = missing_r8

  if (pres > mdl_v(0) .or. pres < mdl_v(n3)) return

  do k = 0,n3-1
     if(pres <= mdl_v(k) .and. pres >= mdl_v(k+1)) then
        zk = real(k) + (mdl_v(k) - pres)/(mdl_v(k) - mdl_v(k+1))
        exit
     endif
  enddo

end subroutine pres_to_zk

!#######################################################################

subroutine height_to_zk(obs_v, mdl_v, n3, zk)

! Calculate the model level "zk" on half (mass) levels,
! corresponding to height "obs_v".

  real(r8), intent(in)  :: obs_v
  integer,  intent(in)  :: n3
  real(r8), intent(in)  :: mdl_v(0:n3)
  real(r8), intent(out) :: zk

  integer   :: k

  zk = missing_r8

  if (obs_v < mdl_v(0) .or. obs_v > mdl_v(n3)) return

  do k = 0,n3-1
     if(obs_v >= mdl_v(k) .and. obs_v <= mdl_v(k+1)) then
        zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
        exit
     endif
  enddo

end subroutine height_to_zk

!#######################################################

subroutine get_model_pressure_profile(i,j,dx,dy,dxm,dym,n,x,id,v_p)

! Calculate the full model pressure profile on half (mass) levels,
! horizontally interpolated at the observation location.

integer,  intent(in)  :: i,j,n,id
real(r8), intent(in)  :: dx,dy,dxm,dym
real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: v_p(0:n)

integer  :: i1,i2,k
real(r8) :: pres1, pres2, pres3, pres4

if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
   j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)) then

   do k=1,n
      pres1 = model_pressure_t(i  ,j  ,k,id,x)
      pres2 = model_pressure_t(i+1,j  ,k,id,x)
      pres3 = model_pressure_t(i  ,j+1,k,id,x)
      pres4 = model_pressure_t(i+1,j+1,k,id,x)
      v_p(k) = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )
   enddo

   if(wrf%dom(id)%surf_obs ) then

!!$      i1 = get_wrf_index(i,j,1,TYPE_PS,id)
!!$      i2 = get_wrf_index(i,j+1,1,TYPE_PS,id)
      i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_PS)
      i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_PS)
      if(x(i1) /= 0.0_r8 .and. x(i1+1) /= 0.0_r8 .and. &
           x(i2) /= 0.0_r8 .and. x(i2+1) /= 0.0_r8) then

         v_p(0) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

      else

         pres1 = model_pressure_t(i  ,j  ,2,id,x)
         pres2 = model_pressure_t(i+1,j  ,2,id,x)
         pres3 = model_pressure_t(i  ,j+1,2,id,x)
         pres4 = model_pressure_t(i+1,j+1,2,id,x)
         v_p(0) = (3.0_r8*v_p(1) - &
              dym*( dxm*pres1 + dx*pres2 ) - dy*( dxm*pres3 + dx*pres4 ))/2.0_r8

      endif

   else

      pres1 = model_pressure_t(i  ,j  ,2,id,x)
      pres2 = model_pressure_t(i+1,j  ,2,id,x)
      pres3 = model_pressure_t(i  ,j+1,2,id,x)
      pres4 = model_pressure_t(i+1,j+1,2,id,x)
      v_p(0) = (3.0_r8*v_p(1) - &
           dym*( dxm*pres1 + dx*pres2 ) - dy*( dxm*pres3 + dx*pres4 ))/2.0_r8

   endif

else

   v_p(:) = missing_r8

endif

end subroutine get_model_pressure_profile

!#######################################################

function model_pressure(i,j,k,id,var_type,x)

! Calculate the pressure at grid point (i,j,k), domain id.
! The grid is defined according to var_type.

integer,  intent(in)  :: i,j,k,id,var_type
real(r8), intent(in)  :: x(:)
real(r8)              :: model_pressure

integer  :: ips, imu
real(r8) :: pres1, pres2

model_pressure = missing_r8

if( (var_type == type_w) .or. (var_type == type_gz) ) then

   if( k == 1 ) then

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i,j,k+1,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   elseif( k == wrf%dom(id)%var_size(3,TYPE_W) ) then

      pres1 = model_pressure_t(i,j,k-1,id,x)
      pres2 = model_pressure_t(i,j,k-2,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   else

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i,j,k-1,id,x)
      model_pressure = (pres1 + pres2)/2.0_r8

   endif

elseif( var_type == type_u ) then

   if( i == wrf%dom(id)%var_size(1,TYPE_U) ) then

      pres1 = model_pressure_t(i-1,j,k,id,x)
      pres2 = model_pressure_t(i-2,j,k,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   elseif( i == 1 ) then

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i+1,j,k,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   else

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i-1,j,k,id,x)
      model_pressure = (pres1 + pres2)/2.0_r8

   endif

elseif( var_type == type_v ) then

   if( j == wrf%dom(id)%var_size(2,TYPE_V) ) then

      pres1 = model_pressure_t(i,j-1,k,id,x)
      pres2 = model_pressure_t(i,j-2,k,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   elseif( j == 1 ) then

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i,j+1,k,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   else

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i,j-1,k,id,x)
      model_pressure = (pres1 + pres2)/2.0_r8

   endif

elseif( var_type == type_mu  .or. var_type == type_tslb .or. &
        var_type == type_ps  .or. var_type == type_u10 .or. &
        var_type == type_v10 .or. var_type == type_t2 .or. &
        var_type == type_q2  .or. var_type == type_tsk) then

   if(wrf%dom(id)%surf_obs ) then

!!$      ips = get_wrf_index(i,j,1,TYPE_PS,id)
      ips = wrf%dom(id)%dart_ind(i,j,1,TYPE_PS)
      model_pressure = x(ips)

   else

!!$      imu = get_wrf_index(i,j,1,TYPE_MU,id)
      imu = wrf%dom(id)%dart_ind(i,j,1,TYPE_MU)
      model_pressure = wrf%dom(id)%p_top + wrf%dom(id)%mub(i,j) + x(imu)

   endif

else

   model_pressure = model_pressure_t(i,j,k,id,x)

endif

end function model_pressure

!#######################################################

function model_pressure_t(i,j,k,id,x)

! Calculate total pressure on mass point (half (mass) levels, T-point).

integer,  intent(in)  :: i,j,k,id
real(r8), intent(in)  :: x(:)
real(r8)              :: model_pressure_t

real (kind=r8), PARAMETER    :: rd_over_rv = gas_constant / gas_constant_v
real (kind=r8), PARAMETER    :: cpovcv = 1.4_r8        ! cp / (cp - gas_constant)

integer  :: iqv,it
real(r8) :: qvf1,rho

model_pressure_t = missing_r8

! Adapted the code from WRF module_big_step_utilities_em.F ----
!         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)

! Simplification: alb*mub = (phb(i,j,k+1) - phb(i,j,k))/dnw(k)

!!$iqv = get_wrf_index(i,j,k,TYPE_QV,id)
!!$it  = get_wrf_index(i,j,k,TYPE_T,id)
iqv = wrf%dom(id)%dart_ind(i,j,k,TYPE_QV)
it  = wrf%dom(id)%dart_ind(i,j,k,TYPE_T)

qvf1 = 1.0_r8 + x(iqv) / rd_over_rv

rho = model_rho_t(i,j,k,id,x)

! .. total pressure:
model_pressure_t = ps0 * ( (gas_constant*(ts0+x(it))*qvf1) / &
     (ps0/rho) )**cpovcv

end function model_pressure_t

!#######################################################

function model_rho_t(i,j,k,id,x)

! Calculate the total density on mass point (half (mass) levels, T-point).

integer,  intent(in)  :: i,j,k,id
real(r8), intent(in)  :: x(:)
real(r8)              :: model_rho_t

integer  :: imu,iph,iphp1
real(r8) :: ph_e

model_rho_t = missing_r8

! Adapted the code from WRF module_big_step_utilities_em.F ----
!         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)

! Simplification: alb*mub = (phb(i,j,k+1) - phb(i,j,k))/dnw(k)

!!$imu = get_wrf_index(i,j,1,TYPE_MU,id)
!!$iph = get_wrf_index(i,j,k,TYPE_GZ,id)
!!$iphp1 = get_wrf_index(i,j,k+1,TYPE_GZ,id)
imu = wrf%dom(id)%dart_ind(i,j,1,TYPE_MU)
iph = wrf%dom(id)%dart_ind(i,j,k,TYPE_GZ)
iphp1 = wrf%dom(id)%dart_ind(i,j,k+1,TYPE_GZ)

ph_e = ( (x(iphp1) + wrf%dom(id)%phb(i,j,k+1)) - (x(iph) + wrf%dom(id)%phb(i,j,k)) ) &
     /wrf%dom(id)%dnw(k)

! now calculate rho = - mu / dphi/deta

model_rho_t = - (wrf%dom(id)%mub(i,j)+x(imu)) / ph_e

end function model_rho_t

!#######################################################

subroutine get_model_height_profile(i,j,dx,dy,dxm,dym,n,x,id,v_h)

! Calculate the model height profile on half (mass) levels,
! horizontally interpolated at the observation location.

integer,  intent(in)  :: i,j,n,id
real(r8), intent(in)  :: dx,dy,dxm,dym
real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: v_h(0:n)

real(r8)  :: fll(n+1)
integer   :: i1,i2,k

if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_GZ) .and. &
   j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_GZ)) then

   do k = 1, wrf%dom(id)%var_size(3,TYPE_GZ)
!!$      i1 = get_wrf_index(i,j,k,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j+1,k,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j,k,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_GZ)
      fll(k) = (dym*( dxm*(wrf%dom(id)%phb(i,j,k)+x(i1)) + &
           dx*(wrf%dom(id)%phb(i+1,j,k)+x(i1+1))) + &
           dy*(dxm*(wrf%dom(id)%phb(i,j+1,k)+x(i2)) + &
           dx*(wrf%dom(id)%phb(i+1,j+1,k)+x(i2+1)) ))/gravity
   end do

   do k=1,n
      v_h(k) = 0.5_r8*(fll(k) + fll(k+1) )
   end do

   v_h(0) = dym*( dxm*wrf%dom(id)%hgt(i,j) + &
                   dx*wrf%dom(id)%hgt(i+1,j) ) + &
             dy*( dxm*wrf%dom(id)%hgt(i,j+1) + &
                   dx*wrf%dom(id)%hgt(i+1,j+1) )

else

   print*,'Not able the get height_profile'
   print*,i,j,dx,dy,dxm,dym,n,id,wrf%dom(id)%var_size(1,TYPE_GZ), &
        wrf%dom(id)%var_size(2,TYPE_GZ)

   v_h(:) =  missing_r8

endif

end subroutine get_model_height_profile

!#######################################################

function model_height(i,j,k,id,var_type,x)

integer,  intent(in)  :: i,j,k,id,var_type
real(r8), intent(in)  :: x(:)
real(r8)              :: model_height

integer   :: i1, i2, i3, i4

model_height = missing_r8

if( (var_type == type_w) .or. (var_type == type_gz) ) then

!!$   i1 = get_wrf_index(i,j,k,TYPE_GZ,id)
   i1 = wrf%dom(id)%dart_ind(i,j,k,TYPE_GZ)
   model_height = (wrf%dom(id)%phb(i,j,k)+x(i1))/gravity

elseif( var_type == type_u ) then

   if( i == wrf%dom(id)%var_size(1,TYPE_U) ) then

!!$      i1 = get_wrf_index(i-1,j,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i-1,j,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i-1,j,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i-1,j,k+1,TYPE_GZ)

      model_height = ( 3.0_r8*(wrf%dom(id)%phb(i-1,j,k  )+x(i1)) &
                      +3.0_r8*(wrf%dom(id)%phb(i-1,j,k+1)+x(i2)) &
                             -(wrf%dom(id)%phb(i-2,j,k  )+x(i1-1)) &
                             -(wrf%dom(id)%phb(i-2,j,k+1)+x(i2-1)) )/(4.0_r8*gravity)

   elseif( i == 1 ) then

!!$      i1 = get_wrf_index(i,j,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j,k+1,TYPE_GZ)

      model_height = ( 3.0_r8*(wrf%dom(id)%phb(i  ,j,k  )+x(i1)) &
                      +3.0_r8*(wrf%dom(id)%phb(i  ,j,k+1)+x(i2)) &
                             -(wrf%dom(id)%phb(i+1,j,k  )+x(i1+1)) &
                             -(wrf%dom(id)%phb(i+1,j,k+1)+x(i2+1)) )/(4.0_r8*gravity)

   else

!!$      i1 = get_wrf_index(i,j,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j,k+1,TYPE_GZ)

      model_height = ( (wrf%dom(id)%phb(i  ,j,k  )+x(i1)) &
                      +(wrf%dom(id)%phb(i  ,j,k+1)+x(i2)) &
                      +(wrf%dom(id)%phb(i-1,j,k  )+x(i1-1)) &
                      +(wrf%dom(id)%phb(i-1,j,k+1)+x(i2-1)) )/(4.0_r8*gravity)

   endif

elseif( var_type == type_v ) then

   if( j == wrf%dom(id)%var_size(2,TYPE_V) ) then

!!$      i1 = get_wrf_index(i,j-1,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j-1,k+1,TYPE_GZ,id)
!!$      i3 = get_wrf_index(i,j-2,k  ,TYPE_GZ,id)
!!$      i4 = get_wrf_index(i,j-2,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j-1,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j-1,k+1,TYPE_GZ)
      i3 = wrf%dom(id)%dart_ind(i,j-2,k  ,TYPE_GZ)
      i4 = wrf%dom(id)%dart_ind(i,j-2,k+1,TYPE_GZ)

      model_height = ( 3.0_r8*(wrf%dom(id)%phb(i,j-1,k  )+x(i1)) &
                      +3.0_r8*(wrf%dom(id)%phb(i,j-1,k+1)+x(i2)) &
                             -(wrf%dom(id)%phb(i,j-2,k  )+x(i3)) &
                             -(wrf%dom(id)%phb(i,j-2,k+1)+x(i4)) )/(4.0_r8*gravity)

   elseif( j == 1 ) then

!!$      i1 = get_wrf_index(i,j  ,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j  ,k+1,TYPE_GZ,id)
!!$      i3 = get_wrf_index(i,j+1,k  ,TYPE_GZ,id)
!!$      i4 = get_wrf_index(i,j+1,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j  ,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_GZ)
      i3 = wrf%dom(id)%dart_ind(i,j+1,k  ,TYPE_GZ)
      i4 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_GZ)

      model_height = ( 3.0_r8*(wrf%dom(id)%phb(i,j  ,k  )+x(i1)) &
                      +3.0_r8*(wrf%dom(id)%phb(i,j  ,k+1)+x(i2)) &
                             -(wrf%dom(id)%phb(i,j+1,k  )+x(i3)) &
                             -(wrf%dom(id)%phb(i,j+1,k+1)+x(i4)) )/(4.0_r8*gravity)

   else

!!$      i1 = get_wrf_index(i,j  ,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j  ,k+1,TYPE_GZ,id)
!!$      i3 = get_wrf_index(i,j-1,k  ,TYPE_GZ,id)
!!$      i4 = get_wrf_index(i,j-1,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j  ,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_GZ)
      i3 = wrf%dom(id)%dart_ind(i,j-1,k  ,TYPE_GZ)
      i4 = wrf%dom(id)%dart_ind(i,j-1,k+1,TYPE_GZ)

      model_height = ( (wrf%dom(id)%phb(i,j  ,k  )+x(i1)) &
                      +(wrf%dom(id)%phb(i,j  ,k+1)+x(i2)) &
                      +(wrf%dom(id)%phb(i,j-1,k  )+x(i3)) &
                      +(wrf%dom(id)%phb(i,j-1,k+1)+x(i4)) )/(4.0_r8*gravity)

   endif

elseif( var_type == type_mu .or. var_type == type_ps .or. &
        var_type == type_tsk) then

   model_height = wrf%dom(id)%hgt(i,j)

elseif( var_type == type_tslb) then

   model_height = wrf%dom(id)%hgt(i,j) - wrf%dom(id)%zs(k)

elseif( var_type == type_u10 .or. var_type == type_v10 ) then

   model_height = wrf%dom(id)%hgt(i,j) + 10.0_r8

elseif( var_type == type_t2 .or. var_type == type_q2 ) then

   model_height = wrf%dom(id)%hgt(i,j) + 2.0_r8

else

!!$   i1 = get_wrf_index(i,j,k  ,TYPE_GZ,id)
!!$   i2 = get_wrf_index(i,j,k+1,TYPE_GZ,id)
   i1 = wrf%dom(id)%dart_ind(i,j,k  ,TYPE_GZ)
   i2 = wrf%dom(id)%dart_ind(i,j,k+1,TYPE_GZ)

   model_height = ( (wrf%dom(id)%phb(i,j,k  )+x(i1)) &
                   +(wrf%dom(id)%phb(i,j,k+1)+x(i2)) )/(2.0_r8*gravity)

endif

end function model_height



subroutine pert_model_state(state, pert_state, interf_provided)
!----------------------------------------------------------------------
! subroutine pert_model_state(state, pert_state, interf_provided)
!
! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

interf_provided = .false.
pert_state = state

end subroutine pert_model_state

!#######################################################

subroutine read_dt_from_wrf_nml()

real(r8) :: dt

integer :: time_step, time_step_fract_num, time_step_fract_den
integer :: max_dom, feedback, smooth_option
integer, dimension(3) :: s_we, e_we, s_sn, e_sn, s_vert, e_vert
integer, dimension(3) :: dx, dy, grid_id, parent_id
integer, dimension(3) :: i_parent_start, j_parent_start, parent_grid_ratio
integer, dimension(3) :: parent_time_step_ratio
integer :: io, iunit, id

namelist /domains/ time_step, time_step_fract_num, time_step_fract_den
namelist /domains/ max_dom
namelist /domains/ s_we, e_we, s_sn, e_sn, s_vert, e_vert
namelist /domains/ dx, dy, grid_id, parent_id
namelist /domains/ i_parent_start, j_parent_start, parent_grid_ratio
namelist /domains/ parent_time_step_ratio
namelist /domains/ feedback, smooth_option


! Begin by reading the namelist input
call find_namelist_in_file("namelist.input", "domains", iunit)
read(iunit, nml = domains, iostat = io)
call check_namelist_read(iunit, io, "domains")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'read_dt_from_wrf_nml','domains namelist values are',' ',' ',' ')
write(logfileunit, nml=domains)
write(     *     , nml=domains)

if (max_dom /= num_domains) then

   write(*,*) 'max_dom in namelist.input = ',max_dom
   write(*,*) 'num_domains in input.nml  = ',num_domains
   call error_handler(E_ERR,'read_dt_from_wrf_nml', &
        'Make them consistent.', source, revision,revdate)

endif

if (time_step_fract_den /= 0) then
   dt = real(time_step) + real(time_step_fract_num) / real(time_step_fract_den)
else
   dt = real(time_step)
endif

do id=1,num_domains
   wrf%dom(id)%dt = dt / real(parent_time_step_ratio(id))
enddo

end subroutine read_dt_from_wrf_nml



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

      character(len=129) :: errstring
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
      PARAMETER (NN=100)
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
      PARAMETER (NN=100)
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
      PARAMETER (NMAX=500)
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
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
END subroutine splint

!cys_add_end

end module model_mod

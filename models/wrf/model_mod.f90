! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next three lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Id$

! Assimilation interface for WRF model

!-----------------------------------------------------------------------
!
!     interface for WRF
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------
!-----------------------------------------------------------------------

use         types_mod, only : r8, deg2rad, missing_r8, ps0, earth_radius
use  time_manager_mod, only : time_type, set_time, set_calendar_type, GREGORIAN
use      location_mod, only : location_type, get_location, set_location, &
                              get_dist, &
                              LocationDims, LocationName, LocationLName, &
                              query_location
use     utilities_mod, only : file_exist, open_file, check_nml_error, &
                              close_file, &
                              register_module, error_handler, E_ERR, &
                              E_MSG, logfileunit
use      obs_kind_mod, only : KIND_U, KIND_V, KIND_PS, KIND_T, KIND_QV, &
                              KIND_P, KIND_W, KIND_QR, KIND_TD, KIND_VR, &
                              KIND_REF, KIND_U10, KIND_V10, KIND_T2, KIND_Q2, &
                              KIND_TD2

use netcdf
use typesizes

implicit none
private

!  public routines and data for the WRF model

public     get_model_size,                    &
           get_state_meta_data,               &
           model_interpolate,                 &
           get_model_time_step,               &
           static_init_model,                 &
           model_get_close_states,            &
           pert_model_state,                  &
           nc_write_model_atts,               &
           nc_write_model_vars

!  public stubs 

public     adv_1step,       &
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
!-----------------------------------------------------------------------

logical :: output_state_vector  = .true.     ! output prognostic variables
integer :: num_moist_vars       = 0
integer :: num_domains          = 1
integer :: calendar_type        = GREGORIAN
logical :: surf_obs             = .false.
character(len = 32) :: adv_mod_command = 'wrf.exe'

namelist /model_nml/ output_state_vector, num_moist_vars, &
                     num_domains, calendar_type, surf_obs, &
                     adv_mod_command


! Public definition of model variable types

integer, parameter :: TYPE_U   = 1,   TYPE_V   = 2,  TYPE_W  = 3,  &
                      TYPE_GZ  = 4,   TYPE_T   = 5,  TYPE_MU = 6,  &
                      TYPE_QV  = 7,   TYPE_QC  = 8,  TYPE_QR = 9,  &
                      TYPE_QI  = 10,  TYPE_QS  = 11, TYPE_QG = 12, &
                      TYPE_U10 = 13,  TYPE_V10 = 14, TYPE_T2 = 15, &
                      TYPE_Q2  = 16


!-----------------------------------------------------------------------

real (kind=r8), PARAMETER    :: gas_constant = 287.04_r8
real (kind=r8), PARAMETER    :: gas_constant_v = 461.51_r8
real (kind=r8), PARAMETER    :: cp = 1004.0_r8

real (kind=r8), PARAMETER    :: kappa = gas_constant / cp

!  Earth constants:
real (kind=r8), PARAMETER    :: gravity = 9.81_r8
INTEGER, PARAMETER           :: v_interp_p = 1, v_interp_h = 2
real (kind=r8), PARAMETER    :: ts0 = 300.0_r8
!

!---- private data ----

TYPE wrf_static_data_for_dart

   integer  :: bt, bts, sn, sns, we, wes
   real(r8) :: p_top, dx, dy, dt
   integer  :: map_proj
   real(r8) :: cen_lat,cen_lon,truelat1,truelat2,stand_lon
   real(r8) :: cone_factor,ycntr,psi1

   integer  :: n_moist
   logical  :: surf_obs
   real(r8), dimension(:),     pointer :: znu, dn, dnw
   real(r8), dimension(:,:),   pointer :: mub, latitude, longitude, hgt
   real(r8), dimension(:,:),   pointer :: mapfac_m, mapfac_u, mapfac_v
   real(r8), dimension(:,:,:), pointer :: phb

   integer :: number_of_wrf_variables
   integer, dimension(:,:), pointer :: var_index
   integer, dimension(:,:), pointer :: var_size
   integer, dimension(:),   pointer :: var_type
   integer, dimension(:,:), pointer :: land

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
integer :: io, ierr, iunit

character (len=80)    :: name
character (len=1)     :: idom
logical, parameter    :: debug = .false.
integer               :: var_id, ind, i, id, dart_index

real(r8) :: theta1,theta2,cell,cell2,psx

!----------------------------------------------------------------------

! Register the module
call register_module(source, revision, revdate)

! Reading the namelist input
if(file_exist('input.nml')) then

   iunit = open_file('input.nml', action = 'read')
   read(iunit, nml = model_nml, iostat = io )
   ierr = check_nml_error(io, 'model_nml')
   call close_file(iunit)

endif

! Record the namelist values in the logfile
write(logfileunit, nml=model_nml)

allocate(wrf%dom(num_domains))

wrf%dom(:)%n_moist = num_moist_vars

if( num_moist_vars > 6) then
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

dart_index = 0

do id=1,num_domains

   write( idom , '(I1)') id

   call check( nf90_open('wrfinput_d0'//idom, NF90_NOWRITE, ncid) )
   if(debug) write(6,*) ' ncid is ',ncid

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

   if(debug) then
      write(6,*) ' dimensions bt, sn, we are ',wrf%dom(id)%bt, &
           wrf%dom(id)%sn, wrf%dom(id)%we
      write(6,*) ' staggered  bt, sn, we are ',wrf%dom(id)%bts, &
           wrf%dom(id)%sns,wrf%dom(id)%wes
   endif

! get meta data and static data we need

   call check( nf90_get_att(ncid, nf90_global, 'DX', wrf%dom(id)%dx) )
   call check( nf90_get_att(ncid, nf90_global, 'DY', wrf%dom(id)%dy) )
   call check( nf90_get_att(ncid, nf90_global, 'DT', wrf%dom(id)%dt) )
   if(debug) write(6,*) ' dx, dy, dt are ',wrf%dom(id)%dx, &
        wrf%dom(id)%dy, wrf%dom(id)%dt

   call check( nf90_inq_varid(ncid, "P_TOP", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%p_top) )
   if(debug) write(6,*) ' p_top is ',wrf%dom(id)%p_top

   call check( nf90_get_att(ncid, nf90_global, 'MAP_PROJ', wrf%dom(id)%map_proj) )
   if(debug) write(6,*) ' map_proj is ',wrf%dom(id)%map_proj

   call check( nf90_get_att(ncid, nf90_global, 'CEN_LAT', wrf%dom(id)%cen_lat) )
   if(debug) write(6,*) ' cen_lat is ',wrf%dom(id)%cen_lat

   call check( nf90_get_att(ncid, nf90_global, 'CEN_LON', wrf%dom(id)%cen_lon) )

   call check( nf90_get_att(ncid, nf90_global, 'TRUELAT1', wrf%dom(id)%truelat1) )
   if(debug) write(6,*) ' truelat1 is ',wrf%dom(id)%truelat1

   call check( nf90_get_att(ncid, nf90_global, 'TRUELAT2', wrf%dom(id)%truelat2) )
   if(debug) write(6,*) ' truelat2 is ',wrf%dom(id)%truelat2

   call check( nf90_get_att(ncid, nf90_global, 'STAND_LON', wrf%dom(id)%stand_lon) )

   if ( abs(wrf%dom(id)%truelat1-wrf%dom(id)%truelat2) .gt. 0.1_r8 ) then
      theta1 = (90.0_r8 - wrf%dom(id)%truelat1)*deg2rad
      theta2 = (90.0_r8 - wrf%dom(id)%truelat2)*deg2rad
      wrf%dom(id)%cone_factor = (log(sin(theta1)) - log(sin(theta2))) &
           / (log(tan(theta1*0.5_r8)) - log(tan(theta2*0.5_r8)))
   else
      wrf%dom(id)%cone_factor = &
           sign(1.0_r8,wrf%dom(id)%truelat1)*sin(wrf%dom(id)%truelat1 * deg2rad)
   end if

   if(debug) write(6,*) 'cone factor = ', wrf%dom(id)%cone_factor

   IF (wrf%dom(id)%map_proj.EQ.1 .OR. wrf%dom(id)%map_proj.EQ.2) THEN
      IF(wrf%dom(id)%cen_lat.LT.0.0_r8)THEN 
         wrf%dom(id)%psi1 = -(90.0_r8+wrf%dom(id)%truelat1)
      ELSE
         wrf%dom(id)%psi1 = 90.0_r8-wrf%dom(id)%truelat1
      ENDIF
   ELSE
      wrf%dom(id)%psi1 = 0.0_r8
   ENDIF

   wrf%dom(id)%psi1 = deg2rad * wrf%dom(id)%psi1

   IF (wrf%dom(id)%map_proj.NE.3) THEN
      psx = (90.0_r8 - wrf%dom(id)%cen_lat)*deg2rad
      IF (wrf%dom(id)%map_proj.EQ.1) THEN
         cell  = earth_radius*SIN(wrf%dom(id)%psi1)/wrf%dom(id)%cone_factor
         cell2 = (TAN(psx/2.0_r8))/(TAN(wrf%dom(id)%psi1/2.0_r8))
      ENDIF
      IF (wrf%dom(id)%map_proj.EQ.2) THEN
         cell  = earth_radius*SIN(psx)/wrf%dom(id)%cone_factor
         cell2 = (1.0_r8 + COS(wrf%dom(id)%psi1))/(1.0_r8 + COS(psx))
      ENDIF
      wrf%dom(id)%ycntr = - cell*(cell2)**wrf%dom(id)%cone_factor
   ENDIF
   ! -----FOR MERCATOR PROJECTION, THE PROJECTION IS TRUE AT LAT AT PHI1
   IF (wrf%dom(id)%map_proj.EQ.3) THEN
      cell = COS(wrf%dom(id)%cen_lat*deg2rad)/(1.0_r8+SIN(wrf%dom(id)%cen_lat*deg2rad))
      wrf%dom(id)%ycntr = - earth_radius*COS(wrf%dom(id)%psi1)* log(cell)
   ENDIF

!  get 1D (z) static data defining grid levels

   allocate(wrf%dom(id)%dn(1:wrf%dom(id)%bt))
   call check( nf90_inq_varid(ncid, "DN", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%dn) )
   if(debug) write(6,*) ' dn ',wrf%dom(id)%dn

   allocate(wrf%dom(id)%znu(1:wrf%dom(id)%bt))
   call check( nf90_inq_varid(ncid, "ZNU", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%znu) )
   if(debug) write(6,*) ' znu is ',wrf%dom(id)%znu

   allocate(wrf%dom(id)%dnw(1:wrf%dom(id)%bt))
   call check( nf90_inq_varid(ncid, "DNW", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%dnw) )
   if(debug) write(6,*) ' dnw is ',wrf%dom(id)%dnw

!  get 2D (x,y) base state for mu, latitude, longitude

   allocate(wrf%dom(id)%mub(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "MUB", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%mub) )
   if(debug) then
      write(6,*) ' corners of mub '
      write(6,*) wrf%dom(id)%mub(1,1),wrf%dom(id)%mub(wrf%dom(id)%we,1),  &
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
      write(6,*) ' corners of land '
      write(6,*) wrf%dom(id)%land(1,1),wrf%dom(id)%land(wrf%dom(id)%we,1),  &
           wrf%dom(id)%land(1,wrf%dom(id)%sn),wrf%dom(id)%land(wrf%dom(id)%we, &
           wrf%dom(id)%sn)
   end if

   if(debug) then
      write(6,*) ' corners of lat '
      write(6,*) wrf%dom(id)%latitude(1,1),wrf%dom(id)%latitude(wrf%dom(id)%we,1),  &
           wrf%dom(id)%latitude(1,wrf%dom(id)%sn), &
           wrf%dom(id)%latitude(wrf%dom(id)%we,wrf%dom(id)%sn)
      write(6,*) ' corners of long '
      write(6,*) wrf%dom(id)%longitude(1,1),wrf%dom(id)%longitude(wrf%dom(id)%we,1),  &
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
      write(6,*) ' corners of phb '
      write(6,*) wrf%dom(id)%phb(1,1,1),wrf%dom(id)%phb(wrf%dom(id)%we,1,1),  &
           wrf%dom(id)%phb(1,wrf%dom(id)%sn,1),wrf%dom(id)%phb(wrf%dom(id)%we, &
           wrf%dom(id)%sn,1)
      write(6,*) wrf%dom(id)%phb(1,1,wrf%dom(id)%bts), &
           wrf%dom(id)%phb(wrf%dom(id)%we,1,wrf%dom(id)%bts),  &
           wrf%dom(id)%phb(1,wrf%dom(id)%sn,wrf%dom(id)%bts), &
           wrf%dom(id)%phb(wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bts)
   end if

! close data file, we have all we need

   call check( nf90_close(ncid) )

!  build the map into the 1D DART vector for WRF data

   wrf%dom(id)%number_of_wrf_variables = 6 + wrf%dom(id)%n_moist
   if( wrf%dom(id)%surf_obs ) then
      wrf%dom(id)%number_of_wrf_variables = wrf%dom(id)%number_of_wrf_variables + 4
   endif
   allocate(wrf%dom(id)%var_type(wrf%dom(id)%number_of_wrf_variables))
   wrf%dom(id)%var_type(1)  = TYPE_U
   wrf%dom(id)%var_type(2)  = TYPE_V
   wrf%dom(id)%var_type(3)  = TYPE_W
   wrf%dom(id)%var_type(4)  = TYPE_GZ
   wrf%dom(id)%var_type(5)  = TYPE_T
   wrf%dom(id)%var_type(6)  = TYPE_MU
   ind = 6
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
   if( wrf%dom(id)%n_moist == 6) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_QG
   end if
   if( wrf%dom(id)%surf_obs ) then
      wrf%dom(id)%var_type(ind + 1) = TYPE_U10
      wrf%dom(id)%var_type(ind + 2) = TYPE_V10
      wrf%dom(id)%var_type(ind + 3) = TYPE_T2
      wrf%dom(id)%var_type(ind + 4) = TYPE_Q2
   end if

! indices into 1D array
   allocate(wrf%dom(id)%var_index(2,wrf%dom(id)%number_of_wrf_variables))
! dimension of variables
   allocate(wrf%dom(id)%var_size(3,wrf%dom(id)%number_of_wrf_variables))

   ind = 1                         ! *** u field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%wes
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
   dart_index = dart_index + 1
   wrf%dom(id)%var_index(1,ind) = dart_index
   dart_index = dart_index - 1 +  &
        wrf%dom(id)%var_size(1,ind)*wrf%dom(id)%var_size(2,ind)*wrf%dom(id)%var_size(3,ind)
   wrf%dom(id)%var_index(2,ind) = dart_index

   ind = ind + 1                   ! *** v field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sns
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
   dart_index = dart_index + 1
   wrf%dom(id)%var_index(1,ind) = dart_index
   dart_index = dart_index - 1 +  &
        wrf%dom(id)%var_size(1,ind)*wrf%dom(id)%var_size(2,ind)*wrf%dom(id)%var_size(3,ind)
   wrf%dom(id)%var_index(2,ind) = dart_index

   ind = ind + 1                   ! *** w field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bts
   dart_index = dart_index + 1
   wrf%dom(id)%var_index(1,ind) = dart_index
   dart_index = dart_index - 1 +  &
        wrf%dom(id)%var_size(1,ind)*wrf%dom(id)%var_size(2,ind)*wrf%dom(id)%var_size(3,ind)
   wrf%dom(id)%var_index(2,ind) = dart_index

   ind = ind + 1                   ! *** geopotential field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bts
   dart_index = dart_index + 1
   wrf%dom(id)%var_index(1,ind) = dart_index
   dart_index = dart_index - 1 +  &
        wrf%dom(id)%var_size(1,ind)*wrf%dom(id)%var_size(2,ind)*wrf%dom(id)%var_size(3,ind)
   wrf%dom(id)%var_index(2,ind) = dart_index

   ind = ind + 1                   ! *** theta field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
   dart_index = dart_index + 1
   wrf%dom(id)%var_index(1,ind) = dart_index
   dart_index = dart_index - 1 +  &
        wrf%dom(id)%var_size(1,ind)*wrf%dom(id)%var_size(2,ind)*wrf%dom(id)%var_size(3,ind)
   wrf%dom(id)%var_index(2,ind) = dart_index

   ind = ind + 1                   ! *** mu field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = 1
   dart_index = dart_index + 1
   wrf%dom(id)%var_index(1,ind) = dart_index
   dart_index = dart_index - 1 +  &
        wrf%dom(id)%var_size(1,ind)*wrf%dom(id)%var_size(2,ind)*wrf%dom(id)%var_size(3,ind)
   wrf%dom(id)%var_index(2,ind) = dart_index

   do i = 1, wrf%dom(id)%n_moist
      ind = ind + 1                   ! *** moisture field ***
      wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
      wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
      wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
      dart_index = dart_index + 1
      wrf%dom(id)%var_index(1,ind) = dart_index
      dart_index = dart_index - 1 +  &
           wrf%dom(id)%var_size(1,ind)*wrf%dom(id)%var_size(2,ind)*wrf%dom(id)%var_size(3,ind)
      wrf%dom(id)%var_index(2,ind) = dart_index
   enddo

   if(wrf%dom(id)%surf_obs ) then
      do i = 1, 4
         ind = ind + 1                   ! *** Surface variable ***
         wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
         wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
         wrf%dom(id)%var_size(3,ind) = 1
         dart_index = dart_index + 1
         wrf%dom(id)%var_index(1,ind) = dart_index
         dart_index = dart_index - 1 +  &
              wrf%dom(id)%var_size(1,ind)*wrf%dom(id)%var_size(2,ind)*wrf%dom(id)%var_size(3,ind)
         wrf%dom(id)%var_index(2,ind) = dart_index
      enddo
   end if

enddo

wrf%model_size = dart_index
if(debug) write(6,*) ' wrf model size is ',wrf%model_size

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

type(time_type) :: get_model_time_step

! Need to translate from wrf model timestep (in seconds) to
! DART time increment

get_model_time_step = set_time(nint(wrf%dom(1)%dt), 0)

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

if(debug) write(6,*) ' ip, jp, kp for index ',ip,jp,kp,index
if(debug) write(6,*) ' Var type: ',var_type

call get_wrf_horizontal_location( ip, jp, var_type, id, lon, lat )

lev = float(kp) ! This is the index of the vertical

if(debug) write(6,*) 'lon, lat, lev: ',lon, lat, lev

! lev is an index here, so which_vert is OK to be hardwired to a 1
location = set_location(lon, lat, lev, 1) 

if(present(var_type_out)) var_type_out = var_type
if(present(id_out)) id_out = id

end subroutine get_state_meta_data

!#######################################################################

subroutine model_interpolate(x, location, obs_kind, obs_val, istatus)

logical, parameter              :: debug = .false.
real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_kind
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus
real(r8)                        :: xloc, yloc, zloc, xyz_loc(3)
integer                         :: which_vert
integer                         :: i, j, k, i1,i2,i3
real(r8)                        :: dx,dy,dxm,dym
real(r8)                        :: p1,p2,p3,p4,a1,alpha
integer                         :: in, ii, id

real(r8), allocatable, dimension(:) :: v_h, v_p, fld, fldu, fldv, fll

istatus = 0

xyz_loc(:) = get_location(location)
which_vert = nint(query_location(location,'which_vert'))
call llxy(xyz_loc(1),xyz_loc(2),xloc,yloc,id)

allocate(v_h(wrf%dom(id)%bt), v_p(wrf%dom(id)%bt), fld(wrf%dom(id)%bt), &
     fldu(wrf%dom(id)%bt), fldv(wrf%dom(id)%bt), fll(wrf%dom(id)%bts))

call toGrid(xloc,i,dx,dxm)
call toGrid(yloc,j,dy,dym)

!  get model pressure profile
call get_model_pressure_profile(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,x,id,v_p, &
                                p1,p2,p3,p4)

if(which_vert == 1) then

!  If obs is by level
   zloc = xyz_loc(3)
   if(debug) print*,' obs is by location and zloc =',zloc
else if(which_vert == 2) then
   ! get pressure vertical co-ordinate
   call to_zk(xyz_loc(3), v_p, wrf%dom(id)%bt, v_interp_p,zloc)
   if(debug.and.obs_kind /= KIND_PS) print*,' obs is by pressure and zloc =',zloc
   if(debug) print*,'model pressure profile'
   if(debug) print*,v_p
else if(which_vert == 3) then
!  get model height profile
   call get_model_height_profile(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,x,id,v_h)
!  get height vertical co-ordinate
   call to_zk(xyz_loc(3), v_h, wrf%dom(id)%bt, v_interp_h,zloc)
   if(debug) print*,' obs is by height and zloc =',zloc
else if(which_vert == -1) then
   ! get height vertical co-ordinate
   zloc = xyz_loc(3)
   if(debug) print*,' obs is at the surface = ', xyz_loc(3)
else
   call error_handler(E_ERR,'model_interpolate', 'wrong option for which_vert', &
        source, revision, revdate)
end if

! Get the desired field to be interpolated
if( obs_kind == KIND_U .or. obs_kind == KIND_V) then        ! U, V

   if(i >= 1 .and. i+2 <= wrf%dom(id)%var_size(1,TYPE_U) .and. &
      j >= 1 .and. j   <  wrf%dom(id)%var_size(2,TYPE_U)) then

      do k=1,wrf%dom(id)%bt
         i1 = get_wrf_index(i,j,k,TYPE_U,id)
         i2 = get_wrf_index(i,j+1,k,TYPE_U,id)

         fldu(k) = dym*(dxm*(x(i1) + x(i1+1))*0.5_r8 + dx*(x(i1+1)+x(i1+2))*0.5_r8) + &
                   dy*(dxm*(x(i2) + x(i2+1))*0.5_r8 + dx*(x(i2+1)+x(i2+2))*0.5_r8)
         if(debug) print*,k,' model u profile ',fldu(k)
      end do

   else

      fldu(:) = missing_r8

   endif

   if(i >= 1 .and. i   <  wrf%dom(id)%var_size(1,TYPE_V) .and. &
      j >= 1 .and. j+2 <= wrf%dom(id)%var_size(2,TYPE_V)) then

      do k=1,wrf%dom(id)%bt
         i1 = get_wrf_index(i,j,k,TYPE_V,id)
         i2 = get_wrf_index(i,j+1,k,TYPE_V,id) 
         i3 = get_wrf_index(i,j+2,k,TYPE_V,id)

         fldv(k) = dym*(dxm*(x(i1) + x(i2))*0.5_r8 + dx*(x(i1+1)+x(i2+1))*0.5_r8) + &
                   dy*(dxm*(x(i2) + x(i3))*0.5_r8 + dx*(x(i2+1)+x(i3+1))*0.5_r8)
         if(debug) print*,k,' model v profile ',fldv(k)
      end do

   else

      fldv(:) = missing_r8

   endif

   alpha = map_to_sphere( xyz_loc(1), xyz_loc(2), id )

   if( obs_kind == KIND_U) then

      fld(:) = fldv(:)*sin(alpha) + fldu(:)*cos(alpha)

   else   ! must want v

      fld(:) = fldv(:)*cos(alpha) - fldu(:)*sin(alpha)

   endif

else if( obs_kind == KIND_T ) then                ! T

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)) then

      do k=1,wrf%dom(id)%bt
         i1 = get_wrf_index(i,j,k,TYPE_T,id)
         i2 = get_wrf_index(i,j+1,k,TYPE_T,id)
         a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
         fld(k) = (ts0 + a1)*(v_p(k)/ps0)**kappa
         if(debug) print*,k,' model temp profile ',fld(k)
      end do

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_W ) then                ! W

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_W) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_W)) then

      do k=1,wrf%dom(id)%var_size(3,TYPE_W)
         i1 = get_wrf_index(i,j,k,TYPE_W,id)
         i2 = get_wrf_index(i,j+1,k,TYPE_W,id)
         fll(k) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
      end do

      do k=1,wrf%dom(id)%bt
         fld(k) = 0.5_r8*(fll(k) + fll(k+1) )
         if(debug) print*,k,' model w profile ',fld(k)
      end do

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_QV) then                ! QV

   in = 0
   do ii = 1, wrf%dom(id)%number_of_wrf_variables 
      if(TYPE_QV == wrf%dom(id)%var_type(ii) ) in = ii
   enddo

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,in) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,in)) then

      do k=1,wrf%dom(id)%bt

         i1 = get_wrf_index(i,j,k,TYPE_QV,id)
         i2 = get_wrf_index(i,j+1,k,TYPE_QV,id)
         a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
         fld(k) = a1 /(1.0_r8 + a1)
         if(debug) print*,k,' model qv profile ',fld(k),i1,i2, &
              x(i1),x(i1+1),x(i2),x(i2+1)
      end do

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_QR) then                ! QR

   in = 0
   do ii = 1, wrf%dom(id)%number_of_wrf_variables 
      if(TYPE_QR == wrf%dom(id)%var_type(ii) ) in = ii
   enddo

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,in) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,in)) then

      do k=1,wrf%dom(id)%bt

         i1 = get_wrf_index(i,j,k,TYPE_QR,id)
         i2 = get_wrf_index(i,j+1,k,TYPE_QR,id)
         fld(k) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
         if(debug) print*,k,' model qr profile ',fld(k),i1,i2, &
              x(i1),x(i1+1),x(i2),x(i2+1)
      end do

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_P) then                 ! Pressure

   fld(:) = v_p(:)

else if( obs_kind == KIND_PS) then                ! Surface pressure

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_MU) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_MU) .and. &
      p1 /= missing_r8 .and. p2 /= missing_r8 .and. &
      p3 /= missing_r8 .and. p4 /= missing_r8) then

      obs_val = wrf%dom(id)%p_top                              +&
           dym*(dxm*wrf%dom(id)%mub(i,j)+dx*wrf%dom(id)%mub(i+1,j))    +&
           dy *(dxm*wrf%dom(id)%mub(i,j+1)+dx*wrf%dom(id)%mub(i+1,j+1))+&
           dym*(dxm*p1 + dx*p2) + dy*(dxm*p3 + dx*p4)
      if(debug) print*,' for sfc model val =',obs_val

   else

      obs_val = missing_r8

   endif

else if( obs_kind == KIND_U10 ) then                ! 10-m U-wind

   in = 0
   do ii = 1, wrf%dom(id)%number_of_wrf_variables 
      if(TYPE_U10 == wrf%dom(id)%var_type(ii) ) in = ii
   enddo

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,in) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,in)) then

      i1 = get_wrf_index(i,j,1,TYPE_U10,id)
      i2 = get_wrf_index(i,j+1,1,TYPE_U10,id)
      obs_val = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

   else

      obs_val = missing_r8

   endif

else if( obs_kind == KIND_V10 ) then                ! 10-m V-wind

   in = 0
   do ii = 1, wrf%dom(id)%number_of_wrf_variables 
      if(TYPE_V10 == wrf%dom(id)%var_type(ii) ) in = ii
   enddo

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,in) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,in)) then

      i1 = get_wrf_index(i,j,1,TYPE_V10,id)
      i2 = get_wrf_index(i,j+1,1,TYPE_V10,id)
      obs_val = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

   else

      obs_val = missing_r8

   endif

else if( obs_kind == KIND_T2 ) then                ! 2-m Temperature

   in = 0
   do ii = 1, wrf%dom(id)%number_of_wrf_variables 
      if(TYPE_T2 == wrf%dom(id)%var_type(ii) ) in = ii
   enddo

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,in) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,in)) then

      i1 = get_wrf_index(i,j,1,TYPE_T2,id)
      i2 = get_wrf_index(i,j+1,1,TYPE_T2,id)
      obs_val = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

   else

      obs_val = missing_r8

   endif

else if( obs_kind == KIND_Q2 ) then                ! 2-m Specific humidity

   in = 0
   do ii = 1, wrf%dom(id)%number_of_wrf_variables 
      if(TYPE_Q2 == wrf%dom(id)%var_type(ii) ) in = ii
   enddo

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,in) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,in)) then

      i1 = get_wrf_index(i,j,1,TYPE_Q2,id)
      i2 = get_wrf_index(i,j+1,1,TYPE_Q2,id)
      obs_val = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

   else

      obs_val = missing_r8

   endif

else
   call error_handler(E_ERR,'model_interpolate', 'wrong obs kind', &
        source, revision, revdate)
end if

! Do 1D interpolation
if(obs_kind /= KIND_PS .and. obs_kind /= KIND_U10 .and. obs_kind /= KIND_V10 &
     .and. obs_kind /= KIND_T2 .and. obs_kind /= KIND_Q2) then
   call Interp_lin_1D(fld, wrf%dom(id)%bt, zloc, obs_val)
endif

if(obs_val == missing_r8) istatus = 1

if(debug) print*,' interpolated value= ',obs_val

deallocate(v_h, v_p, fld, fldu, fldv, fll)

end subroutine model_interpolate

!#######################################################################

subroutine model_get_close_states(o_loc, radius, num_total, indices, dist, x)

type(location_type), intent(in) :: o_loc
real(r8),            intent(in) :: radius
integer,            intent(out) :: num_total, indices(:)
real(r8),           intent(out) :: dist(:)
real(r8),            intent(in) :: x(:)

integer                         :: u_pts, v_pts, p_pts, w_pts
integer                         :: i,k,indmax, ii, jj, id
integer                         :: num, max_size
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
        radius, id, num, lon_ind, lat_ind, vert_ind, close_dist, &
        u_pts, v_pts, p_pts, w_pts, x )

! start with u_pts

   do i = 1, u_pts

      ii = lon_ind(i)
      jj = lat_ind(i)
      k = vert_ind(i)

      num_total = num_total + 1
      if(num_total <= indmax) then
         indices(num_total) = get_wrf_index(ii,jj,k,type_u,id)
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
         indices(num_total) = get_wrf_index(ii,jj,k,type_v,id)
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
         indices(num_total) = get_wrf_index(ii,jj,k,type_w,id)
         dist(num_total) = close_dist(i)
      end if
      num_total = num_total + 1
      if(num_total <= indmax) then
         indices(num_total) = get_wrf_index(ii,jj,k,type_gz,id)
         dist(num_total) = close_dist(i)
      end if

! 2D (near-surface) fields.

      if (k == 1) then

         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_mu,id)
            dist(num_total) = close_dist(i)
         end if

         if(wrf%dom(id)%surf_obs ) then

            num_total = num_total + 1
            if(num_total <= indmax) then
               indices(num_total) = get_wrf_index(ii,jj,k,type_u10,id)
               dist(num_total) = close_dist(i)
            end if

            num_total = num_total + 1
            if(num_total <= indmax) then
               indices(num_total) = get_wrf_index(ii,jj,k,type_v10,id)
               dist(num_total) = close_dist(i)
            end if

            num_total = num_total + 1
            if(num_total <= indmax) then
               indices(num_total) = get_wrf_index(ii,jj,k,type_t2,id)
               dist(num_total) = close_dist(i)
            end if

            num_total = num_total + 1
            if(num_total <= indmax) then
               indices(num_total) = get_wrf_index(ii,jj,k,type_q2,id)
               dist(num_total) = close_dist(i)
            end if

         end if

      end if

   enddo

! last -> p_pts (t + wrf%dom(id)%n_moist variables)

   do i = u_pts+v_pts+w_pts+1, num

      ii = lon_ind(i)
      jj = lat_ind(i)
      k = vert_ind(i)

      num_total = num_total + 1
      if(num_total <= indmax) then
         indices(num_total) = get_wrf_index(ii,jj,k,type_t,id)
         dist(num_total) = close_dist(i)
      end if
      if( wrf%dom(id)%n_moist >= 1) then
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qv,id)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( wrf%dom(id)%n_moist >= 2) then
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qc,id)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( wrf%dom(id)%n_moist >= 3) then
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qr,id)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( wrf%dom(id)%n_moist >= 4) then
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qi,id)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( wrf%dom(id)%n_moist >= 5) then
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qs,id)
            dist(num_total) = close_dist(i)
         end if
      end if
      if( wrf%dom(id)%n_moist == 6) then
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qg,id)
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
integer :: in, ii

character(len=129) :: errstring

in = 0
do ii = 1, wrf%dom(id)%number_of_wrf_variables 
   if(var_type == wrf%dom(id)%var_type(ii) ) in = ii
enddo

if(i >= 1 .and. i <= wrf%dom(id)%var_size(1,in) .and. &
   j >= 1 .and. j <= wrf%dom(id)%var_size(2,in) .and. &
   k >= 1 .and. k <= wrf%dom(id)%var_size(3,in)) then

   get_wrf_index = wrf%dom(id)%var_index(1,in)-1 +   &
        i + wrf%dom(id)%var_size(1,in)*(j-1) + &
        (wrf%dom(id)%var_size(1,in)*wrf%dom(id)%var_size(2,in))*(k-1)

else

  write(errstring,*)'Indices ',i,j,k,' exceed grid dimensions: ', &
       wrf%dom(id)%var_size(1,in), &
       wrf%dom(id)%var_size(2,in),wrf%dom(id)%var_size(3,in)
  call error_handler(E_ERR,'get_wrf_index', errstring, &
       source, revision, revdate)

endif

end function get_wrf_index

!-----------------------------------------------------------------

subroutine grid_close_states( o_loc, radius_in, id, num,  &
                              close_lon_ind, close_lat_ind, close_vert_ind, &
                              close_dist, u_pts, v_pts, p_pts, w_pts, x )

! Finds close state points from a particular grid for the WRF model

type(location_type), intent(in) :: o_loc
real(r8),            intent(in) :: radius_in
integer,             intent(in) :: id
integer,            intent(out) :: num
integer,            intent(out) :: close_lon_ind(:), close_lat_ind(:), close_vert_ind(:)
real(r8),           intent(out) :: close_dist(:)
integer,            intent(out) :: u_pts, v_pts, p_pts, w_pts
real(r8),            intent(in) :: x(:)

real(r8)            :: radius
real(r8)            :: rad, radn, dxr, dyr, sdx, sdy, gdist
integer             :: i_closest, j_closest, ixmin, jymin, ixmax, jymax

integer             :: i, j, k, n, m

logical,  parameter :: debug = .false.  

if(debug) write(6,*) ' in grid_close_states '

! use meaningful units for radius -- convert radians to meters 
radius = radius_in*earth_radius*1000.0_r8
if(debug) write(6,*) ' radius in grid_close_states is ',radius

! Get index to closest lat and lon for this observation

n = size( wrf%dom(id)%latitude, 1 )
m = size( wrf%dom(id)%latitude, 2 )

if(debug) write(6,*) 'obs location: ',get_location(o_loc)

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
      if(debug) write(6,*) i,j,radn,rad
   enddo
enddo

if(debug) write(6,*) ' closest wrf long and lat is ',i_closest,j_closest
if(debug) write(6,*) ' at distance ',rad

! define box edges for radius check
dxr = 1.0_r8 + radius/wrf%dom(id)%dx  !  radius in multiples of dx
dyr = 1.0_r8 + radius/wrf%dom(id)%dy  !  radius in multiples of dy

if(debug) write(6,*) ' dxr, dyr in grid_close_states ',dxr,dyr


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
   write(6,*) ' ixmin, ixmax, jymin, jymax are '
   write(6,*) ixmin, ixmax, jymin, jymax
endif

!  we have bounding box, get and check distances.
!  first, convert radius back to radians

radius = radius_in
num = 0

if(debug) write(6,*) 'radius (radian) ',radius

! check distance for u points, expand box so that
! we don't leave possible points out of check

do k = 1, wrf%dom(id)%bt
   do j = jymin, jymax
      do i = max(1,ixmin-1), ixmax+1

         gdist = get_dist_wrf(i,j,k, type_u, o_loc, id, x)
         if ( gdist <= radius ) then
            num = num + 1
            if(debug) write(6,*) ' u pt ',num,i,j,k, id,gdist
            close_lon_ind(num) = i
            close_lat_ind(num) = j
            close_vert_ind(num) = k
            close_dist(num) = gdist
         end if

      enddo
   enddo
enddo

u_pts = num

do k = 1, wrf%dom(id)%bt
   do j = max(1,jymin-1), jymax+1
      do i = ixmin, ixmax

         gdist = get_dist_wrf(i,j,k, type_v, o_loc, id, x)
         if ( gdist <= radius ) then
            num = num + 1
            close_lon_ind(num) = i
            close_lat_ind(num) = j
            close_vert_ind(num) = k
            close_dist(num) = gdist
            if(debug) write(6,*) ' v pt ',num,i,j,gdist
         end if
      enddo
   enddo
enddo

v_pts = num - u_pts

do k = 1, wrf%dom(id)%bts
   do j = jymin, jymax
      do i = ixmin, ixmax

         gdist = get_dist_wrf(i,j,k, type_w, o_loc, id, x)
         if ( gdist <= radius ) then
            num = num + 1
            close_lon_ind(num) = i
            close_lat_ind(num) = j
            close_vert_ind(num) = k
            close_dist(num) = gdist
            if(debug) write(6,*) ' w pt ',num,i,j,gdist
         end if
      enddo
   enddo
enddo

w_pts = num - u_pts - v_pts

do k = 1, wrf%dom(id)%bt
   do j = jymin, jymax
      do i = ixmin, ixmax

         gdist = get_dist_wrf(i,j,k, type_t, o_loc, id, x)
         if ( gdist <= radius ) then
            num = num + 1
            close_lon_ind(num) = i
            close_lat_ind(num) = j
            close_vert_ind(num) = k
            close_dist(num) = gdist
            if(debug) write(6,*) ' p pt ',num,i,j,gdist
         end if

      enddo
   enddo
enddo

p_pts = num - u_pts - v_pts - w_pts

end subroutine grid_close_states

!***********************************************************************

function get_dist_wrf( i, j, k, var_type, o_loc, id, x )

real(r8)                        :: get_dist_wrf 
type(location_type), intent(in) :: o_loc
integer,             intent(in) :: i,j,k,var_type, id
real(r8),            intent(in) :: x(:)

type(location_type) :: loc
real(r8)            :: long, lat, vloc
integer             :: which_vert

real(r8)                            :: p1,p2,p3,p4
real(r8), allocatable, dimension(:) :: fld, fll

character(len=129) :: errstring

allocate(fld(wrf%dom(id)%bt), fll(wrf%dom(id)%bts))

! get distance for input var_type

call get_wrf_horizontal_location( i, j, var_type, id, long, lat )

! We will set which_vert based on the incoming location type
! This will ensure compatibility for measuring 'horizontal'
! and cannot think of anything better for vertical.

which_vert = nint(query_location(o_loc,'which_vert'))

if(which_vert == 1 ) then

   if( (var_type == type_w ) .or. (var_type == type_gz) ) then
      vloc = float(k) - 0.5_r8
   else
      vloc = float(k)
   endif

elseif(which_vert == 2 ) then

   call get_model_pressure_profile(i,j,0.0_r8,0.0_r8,1.0_r8,1.0_r8, &
        wrf%dom(id)%bt,x,id,fld, p1,p2,p3,p4)
   if( (var_type == type_w ) .or. (var_type == type_gz) ) then
      if (k == 1 ) then
         vloc = (3.0_r8*fld(1) - fld(2))/2.0_r8
      elseif (k == wrf%dom(id)%bts) then
         vloc = (3.0_r8*fld(wrf%dom(id)%bt) - fld(wrf%dom(id)%bt-1))/2.0_r8
      else
         vloc = 0.5_r8*(fld(k)+fld(k-1))
      endif
   else
      vloc = fld(k)
   endif

elseif(which_vert == 3 ) then

   if( (var_type == type_w ) .or. (var_type == type_gz) ) then
      call get_model_height_profile_stag(i,j,0.0_r8,0.0_r8,1.0_r8,1.0_r8, &
           wrf%dom(id)%bt,x,id,fll)
      vloc = fll(k)
   else
      call get_model_height_profile(i,j,0.0_r8,0.0_r8,1.0_r8,1.0_r8, &
           wrf%dom(id)%bt,x,id,fld)
      vloc = fld(k)
   endif

elseif(which_vert == -1 ) then

   vloc = 0.0_r8

else
   write(errstring, *) 'Which_vert = ',which_vert, ' is not allowed.'
   call error_handler(E_ERR,'get_dist_wrf', errstring, source, revision, revdate)
endif

if (vloc == missing_r8 ) then
   print*,i,j,k, id, p1,var_type,which_vert
   print*,fld
   write(errstring, *) 'Unable to define vloc.'
   call error_handler(E_ERR,'get_dist_wrf', errstring, source, revision, revdate)
endif

loc          = set_location(long, lat, vloc, which_vert)
get_dist_wrf = get_dist(loc,o_loc)

deallocate(fld, fll)

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
     btDimID, btStagDimID

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
   if( wrf%dom(num_domains)%n_moist == 6) then
      call check(nf90_put_att(ncFileID, StateVarId, "QG_units","kg/kg"))
   endif
   if(wrf%dom(num_domains)%surf_obs ) then
      call check(nf90_put_att(ncFileID, StateVarId, "U10_units","m/s"))
      call check(nf90_put_att(ncFileID, StateVarId, "V10_units","m/s"))
      call check(nf90_put_att(ncFileID, StateVarId, "T2_units","K"))
      call check(nf90_put_att(ncFileID, StateVarId, "Q2_units","kg/kg"))
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

   if( wrf%dom(id)%n_moist == 6) then
      call check(nf90_def_var(ncid=ncFileID, name="QGRAUP_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, "description", &
           "Graupel mixing ratio"))
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

   end if

enddo

endif

!-----------------------------------------------------------------
! Fill the variables we can
!-----------------------------------------------------------------
call check(nf90_enddef(ncfileID))

do id=1,num_domains
   call check(nf90_put_var(ncFileID,       DXVarID, wrf%dom(id)%dx        ))
   call check(nf90_put_var(ncFileID,       DYVarID, wrf%dom(id)%dy        ))
   call check(nf90_put_var(ncFileID, TRUELAT1VarID, wrf%dom(id)%truelat1  ))
   call check(nf90_put_var(ncFileID, TRUELAT2VarID, wrf%dom(id)%truelat2  ))
   call check(nf90_put_var(ncFileID,  CEN_LATVarID, wrf%dom(id)%cen_lat   ))
   call check(nf90_put_var(ncFileID,  CEN_LONVarID, wrf%dom(id)%cen_lon   ))
   call check(nf90_put_var(ncFileID, MAP_PROJVarID, wrf%dom(id)%map_proj  ))

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

do id=1,num_domains

   write( idom , '(I1)') id

   !----------------------------------------------------------------------------
   ! Fill the variables, the order is CRITICAL  ...   U,V,W,GZ,T,MU,QV,QC,QR,...
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   varname = 'U_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = 1
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
   ! deallocate(temp3d)  remaining 3D variables same size ...


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


   if( wrf%dom(id)%n_moist >= 1) then
      !----------------------------------------------------------------------------
      varname = 'QVAPOR_d0'//idom
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
   if( wrf%dom(id)%n_moist == 6) then
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

   deallocate(temp3d)

   if ( wrf%dom(id)%n_moist > 6 ) then
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

   endif

   deallocate(temp2d)

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



subroutine llxy (xloni,xlatj,x,y,id)
!-----------------------------------------------------------------
!
!                 ROUTINE LLXY
!                **************
!
!
! PURPOSE:  CALCULATES THE (X,Y) LOCATION (DOT) IN THE MESOSCALE GRIDS
! -------   FROM LATITUDES AND LONGITUDES
!
!
!  INPUT:
!  -----
!   XLAT:    LATITUDES
!   XLON:    LONGITUDES
!
! OUTPUT:
! -----
!   X:        THE COORDINATE IN X (I)-DIRECTION.
!   Y:        THE COORDINATE IN Y (J)-DIRECTION.
!
!-----------------------------------------------------------------

  real(r8), intent(in)  :: xloni, xlatj
  real(r8), intent(out) :: x, y
  integer,  intent(out) :: id

  real(r8) :: dxlon
  real(r8) :: xlat, xlon
  real(r8) :: xx, yy, xc, yc
  real(r8) :: cell, psi0, psx, r, flp
  real(r8) :: centri, centrj
  real(r8) :: ds
  real(r8) :: bb,c2
  logical  :: dom_found

  dom_found = .false.

  xlon = xloni
  xlat = xlatj
  xlat = max (xlat, -89.9999_r8)
  xlat = min (xlat, +89.9999_r8)

  id = num_domains
  do while (.not. dom_found)

     ds = 0.001_r8 *wrf%dom(id)%dx
     c2 = earth_radius * COS(wrf%dom(id)%psi1)

     if (wrf%dom(id)%map_proj == 3) then
        xc = 0.0_r8
        yc = wrf%dom(id)%ycntr

        cell = cos(xlat*deg2rad)/(1.0_r8+sin(xlat*deg2rad))
        yy = -c2*log(cell)
        xx = c2*(xlon-wrf%dom(id)%cen_lon)*deg2rad

     else

        psi0 = ( 90.0_r8 - wrf%dom(id)%cen_lat)*deg2rad
        xc = 0.0_r8

!-----CALCULATE X,Y COORDS. RELATIVE TO POLE

        dxlon = xlon - wrf%dom(id)%cen_lon
        if (dxlon >  180.0_r8) dxlon = dxlon - 360.0_r8
        if (dxlon < -180.0_r8) dxlon = dxlon + 360.0_r8

        flp = wrf%dom(id)%cone_factor*dxlon*deg2rad

        psx = ( 90.0_r8 - xlat )*deg2rad

        if (wrf%dom(id)%map_proj == 2) then
! ...... Polar stereographics:
           bb = 2.0_r8*(cos(wrf%dom(id)%psi1/2.0_r8)**2)
           yc = -earth_radius*bb*tan(psi0/2.0_r8)
           r = -earth_radius*bb*tan(psx/2.0_r8)
        else
! ...... Lambert conformal:
           bb = -earth_radius/wrf%dom(id)%cone_factor*sin(wrf%dom(id)%psi1)
           yc = bb*(tan(psi0/2.0_r8)/tan(wrf%dom(id)%psi1/2.0_r8))**wrf%dom(id)%cone_factor
           r = bb*(tan(psx /2.0_r8)/tan(wrf%dom(id)%psi1/2.0_r8))**wrf%dom(id)%cone_factor
        endif

        if (wrf%dom(id)%cen_lat < 0.0_r8) then
           xx = r*sin(flp)
           yy = r*cos(flp)
        else
           xx = -r*sin(flp)
           yy =  r*cos(flp)
        endif

     endif

! TRANSFORM (1,1) TO THE ORIGIN
! the location of the center in the coarse domain

     centri = real (wrf%dom(id)%we)/2.0_r8
     centrj = real (wrf%dom(id)%sn)/2.0_r8
! the (X,Y) coordinates in the coarse domain
     x = ( xx - xc )/ds + centri
     y = ( yy - yc )/ds + centrj

!--only add 0.5 so that x/y is relative to first cross points (MM5 input):

     x = (x - 1.0_r8) + 0.5_r8
     y = (y - 1.0_r8) + 0.5_r8

     if ( (x >= 1 .and. x <= wrf%dom(id)%we .and. &
          y >= 1 .and. y <= wrf%dom(id)%sn) .or. id == 1 ) then

        dom_found = .true.

     else

        id = id - 1

     endif

  end do

end subroutine llxy

!**********************************************
subroutine Interp_lin_1D(fi1d, n1, z, fo1d)

  integer,  intent(in)  :: n1
  real(r8), intent(in)  :: fi1d(n1)
  real(r8), intent(in)  :: z
  real(r8), intent(out) :: fo1d

  integer   :: k
  real(r8)  :: dz, dzm

  call toGrid(z, k, dz, dzm)
  if(k >= 1 .and. k < n1) then
     fo1d = dzm*fi1d(k) + dz*fi1d(k+1)
  else
     fo1d = missing_r8
  endif

end subroutine Interp_lin_1D

!**********************************************
subroutine Interp_lin_2D(fi2d,n1,n2, x,y, fo2d)

  integer,  intent(in)  :: n1, n2
  real(r8), intent(in)  :: fi2d(n1,n2)
  real(r8), intent(in)  :: x, y
  real(r8), intent(out) :: fo2d

  integer   :: i, j
  real(r8)  :: dx, dxm, dy, dym

  fo2d = missing_r8

  call toGrid (x,i,dx,dxm)
  call toGrid (y,j,dy,dym)

  if(i >= 1 .and. i < n1 .and. j >= 1 .and. j < n2) then
     fo2d = dym*(dxm*fi2d(i,j  ) + dx*fi2d(i+1,j  )) &
          + dy *(dxm*fi2d(i,j+1) + dx*fi2d(i+1,j+1))
  endif

end subroutine Interp_lin_2D
   
!**********************************************

subroutine Interp_lin_3D(fi3d,n1,n2,n3, x,y,z,fo3d)

  integer,  intent(in)  :: n1,n2,n3
  real(r8), intent(in)  :: fi3d(n1,n2,n3)
  real(r8), intent(in)  :: x, y, z
  real(r8), intent(out) :: fo3d

  integer   :: i, j, k
  real(r8)  :: dx, dxm, dy, dym, dz, dzm
  real(r8)  :: fiz(n3)

  fo3d = missing_r8

  call toGrid (x,i,dx,dxm)
  call toGrid (y,j,dy,dym)
  call toGrid (z,k,dz,dzm)

  if(i >= 1 .and. i < n1 .and. j >= 1 .and. j < n2 .and. k >= 1 .and. k < n3) then
     fiz(1:n3) = dym*(dxm*fi3d(i, j,   1:n3) + dx *fi3d(i+1 ,j  ,1:n3))&
          + dy *(dxm*fi3d(i, j+1, 1:n3) + dx *fi3d(i+1, j+1,1:n3))
     fo3d = dzm*fiz(k) + dz*fiz(k+1)
  endif

end subroutine Interp_lin_3D

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

subroutine to_zk(obs_v, mdl_v, n3, v_interp_optn, zk)

   real(r8), intent(in)  :: obs_v
   integer,  intent(in)  :: n3, v_interp_optn
   real(r8), intent(in)  :: mdl_v(n3)
   real(r8), intent(out) :: zk

   integer   :: k

   zk = missing_r8

   if(v_interp_optn == v_interp_p) then

      if (obs_v > mdl_v(1) .or. obs_v < mdl_v(n3)) return

      do k = 1,n3-1
         if(obs_v <= mdl_v(k) .and. obs_v >= mdl_v(k+1)) then
            zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
            exit
         endif
      enddo
   else if(v_interp_optn == v_interp_h) then
      if (obs_v < mdl_v(1) .or. obs_v > mdl_v(n3)) return

      do k = 1,n3-1
         if(obs_v >= mdl_v(k) .and. obs_v <= mdl_v(k+1)) then
            zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
            exit
         endif
      enddo
   endif

end subroutine to_zk

!#######################################################

subroutine get_model_pressure_profile(i_in,j_in,dx,dy,dxm,dym,n,x,id,fld, &
                                      pp11,pp21,pp31,pp41)

integer,  intent(in)  :: i_in,j_in,n,id
real(r8), intent(in)  :: dx,dy,dxm,dym
real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: fld(n)
real(r8), intent(out) :: pp11,pp21,pp31,pp41

integer               :: i1,i2,k,q1,q2,q3,q4
real(r8)              :: qv1,qv2,qv3,qv4
real(r8), dimension(n):: pp1,pp2,pp3,pp4,pb,pp
integer               :: in, ii, i, j

i = i_in
j = j_in

in = 0
do ii = 1, wrf%dom(id)%number_of_wrf_variables
   if(TYPE_QV == wrf%dom(id)%var_type(ii) ) in = ii
enddo

if((dx /= 0.0_r8) .or. (dy /= 0.0_r8)) then

if(i >= 1 .and. i < wrf%dom(id)%var_size(1,in) .and. &
   j >= 1 .and. j < wrf%dom(id)%var_size(2,in)) then

   do k=1,n
      pb(k) = wrf%dom(id)%p_top + wrf%dom(id)%znu(k)*  &
           (  dym*(dxm*wrf%dom(id)%mub(i,j  ) + dx*wrf%dom(id)%mub(i+1,j  )) + &
              dy *(dxm*wrf%dom(id)%mub(i,j+1) + dx*wrf%dom(id)%mub(i+1,j+1)) )
   end do
   i1 = get_wrf_index(i,j,1,TYPE_MU,id)
   i2 = get_wrf_index(i,j+1,1,TYPE_MU,id)
   q1 = get_wrf_index(i,j,n,TYPE_QV,id)
   q2 = get_wrf_index(i,j+1,n,TYPE_QV,id)
   qv1 = x(q1)/(1.0_r8+x(q1))
   qv2 = x(q1+1)/(1.0_r8+x(q1+1))
   qv3 = x(q2)/(1.0_r8+x(q2))
   qv4 = x(q2+1)/(1.0_r8+x(q2+1))
   pp1(n) = -0.5_r8 *(x(i1)  +qv1*wrf%dom(id)%mub(i  ,j  ))*wrf%dom(id)%dnw(n)*(1.0_r8 + x(q1))
   pp2(n) = -0.5_r8 *(x(i1+1)+qv2*wrf%dom(id)%mub(i+1,j  ))*wrf%dom(id)%dnw(n)*(1.0_r8 + x(q1+1))
   pp3(n) = -0.5_r8 *(x(i2)  +qv3*wrf%dom(id)%mub(i  ,j+1))*wrf%dom(id)%dnw(n)*(1.0_r8 + x(q2))
   pp4(n) = -0.5_r8 *(x(i2+1)+qv4*wrf%dom(id)%mub(i+1,j+1))*wrf%dom(id)%dnw(n)*(1.0_r8 + x(q2+1))
   pp(n)  = dym*(dxm*pp1(n)+dx*pp2(n)) + dy*(dxm*pp3(n)+dx*pp4(n))
   fld(n) = pp(n) + pb(n)
   do k= n-1,1,-1
      q1 = get_wrf_index(i,j,k,TYPE_QV,id)
      q2 = get_wrf_index(i,j+1,k,TYPE_QV,id)
      q3 = get_wrf_index(i,j,k+1,TYPE_QV,id)
      q4 = get_wrf_index(i,j+1,k+1,TYPE_QV,id)
      qv1 = 0.5_r8*(x(q1)+x(q3))/(1.0_r8+0.5_r8*(x(q1)+x(q3)))
      qv2 = 0.5_r8*(x(q1+1)+x(q3+1))/(1.0_r8+0.5_r8*(x(q1+1)+x(q3+1)))
      qv3 = 0.5_r8*(x(q2)+x(q4))/(1.0_r8+0.5_r8*(x(q2)+x(q4)))
      qv4 = 0.5_r8*(x(q2+1)+x(q4+1))/(1.0_r8+0.5_r8*(x(q2+1)+x(q4+1)))

      pp1(k) = pp1(k+1) -(x(i1)  +qv1*wrf%dom(id)%mub(i  ,j  ))*wrf%dom(id)%dn(k+1)* &
           (1.0_r8 + 0.5_r8*(x(q1) + x(q3)))
      pp2(k) = pp2(k+1) -(x(i1+1)+qv2*wrf%dom(id)%mub(i+1,j  ))*wrf%dom(id)%dn(k+1)* &
           (1.0_r8 + 0.5_r8*(x(q1+1) + x(q3+1)))
      pp3(k) = pp3(k+1) -(x(i2)  +qv3*wrf%dom(id)%mub(i  ,j+1))*wrf%dom(id)%dn(k+1)* &
           (1.0_r8 + 0.5_r8*(x(q2) + x(q4)))
      pp4(k) = pp4(k+1) -(x(i2+1)  +qv4*wrf%dom(id)%mub(i+1,j+1))*wrf%dom(id)%dn(k+1)* &
           (1.0_r8 + 0.5_r8*(x(q2+1) + x(q4+1)))

      pp(k)  = dym*(dxm*pp1(k)+dx*pp2(k)) + dy*(dxm*pp3(k)+dx*pp4(k))
      fld(k) = pp(k) + pb(k)
   end do

   pp11 = pp1(1)
   pp21 = pp2(1)
   pp31 = pp3(1)
   pp41 = pp4(1)

else

   fld(:) = missing_r8
   pp11 = missing_r8
   pp21 = missing_r8
   pp31 = missing_r8
   pp41 = missing_r8

endif

else

! Cheating to include the last grid point of U or V.

if(i == wrf%dom(id)%var_size(1,in)+1) i = wrf%dom(id)%var_size(1,in)
if(j == wrf%dom(id)%var_size(2,in)+1) j = wrf%dom(id)%var_size(2,in)

if(i >= 1 .and. i <= wrf%dom(id)%var_size(1,in) .and. &
   j >= 1 .and. j <= wrf%dom(id)%var_size(2,in)) then

   do k=1,n
      pb(k) = wrf%dom(id)%p_top + wrf%dom(id)%znu(k)*wrf%dom(id)%mub(i,j)
   end do
   i1 = get_wrf_index(i,j,1,TYPE_MU,id)
   q1 = get_wrf_index(i,j,n,TYPE_QV,id)
   qv1 = x(q1)/(1.0_r8+x(q1))
   pp1(n) = -0.5_r8 *(x(i1)  +qv1*wrf%dom(id)%mub(i,j))*wrf%dom(id)%dnw(n)*(1.0_r8 + x(q1))
   pp(n)  = pp1(n)
   fld(n) = pp(n) + pb(n)
   do k= n-1,1,-1
      q1 = get_wrf_index(i,j,k,TYPE_QV,id)
      q3 = get_wrf_index(i,j,k+1,TYPE_QV,id)
      qv1 = 0.5_r8*(x(q1)+x(q3))/(1.0_r8+0.5_r8*(x(q1)+x(q3)))

      pp1(k) = pp1(k+1) -(x(i1)  +qv1*wrf%dom(id)%mub(i,j))*wrf%dom(id)%dn(k+1)* &
           (1.0_r8 + 0.5_r8*(x(q1) + x(q3)))

      pp(k)  = pp1(k)
      fld(k) = pp(k) + pb(k)
   end do

   pp11 = pp1(1)
   pp21 = missing_r8
   pp31 = missing_r8
   pp41 = missing_r8

else

   fld(:) = missing_r8
   pp11 = missing_r8
   pp21 = missing_r8
   pp31 = missing_r8
   pp41 = missing_r8

endif

endif

end subroutine get_model_pressure_profile

!#######################################################

subroutine get_model_height_profile(i_in,j_in,dx,dy,dxm,dym,n,x,id,fld)

integer,  intent(in)  :: i_in,j_in,n,id
real(r8), intent(in)  :: dx,dy,dxm,dym
real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: fld(n)

real(r8)  :: fll(n+1)
integer   :: i, j, k

i = i_in
j = j_in

! Cheating to include the last grid point of U or V.

if((dx == 0.0_r8) .and. (dy == 0.0_r8)) then
   if(i == wrf%dom(id)%var_size(1,TYPE_GZ)+1) i = wrf%dom(id)%var_size(1,TYPE_GZ)
   if(j == wrf%dom(id)%var_size(2,TYPE_GZ)+1) j = wrf%dom(id)%var_size(2,TYPE_GZ)
endif

if(i >= 1 .and. i <= wrf%dom(id)%var_size(1,TYPE_GZ) .and. &
     j >= 1 .and. j <= wrf%dom(id)%var_size(2,TYPE_GZ)) then

   call get_model_height_profile_stag(i,j,dx,dy,dxm,dym,n,x,id,fll)

   do k=1,n
      fld(k) = 0.5_r8*(fll(k) + fll(k+1) )
   end do

else

   print*,'Not able the get height_profile'
   print*,i,j,dx,dy,dxm,dym,n,id,wrf%dom(id)%var_size(1,TYPE_GZ), &
        wrf%dom(id)%var_size(2,TYPE_GZ)

   fld(:) =  missing_r8

endif

end subroutine get_model_height_profile



!#######################################################

subroutine get_model_height_profile_stag(i_in,j_in,dx,dy,dxm,dym,n,x,id,fll)

integer,  intent(in)  :: i_in,j_in,n,id
real(r8), intent(in)  :: dx,dy,dxm,dym
real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: fll(n+1)

integer   :: i1,i2,i,j,k

i = i_in
j = j_in

if((dx /= 0.0_r8) .or. (dy /= 0.0_r8)) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_GZ) .and. &
        j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_GZ)) then

      do k = 1, wrf%dom(id)%var_size(3,TYPE_GZ)
         i1 = get_wrf_index(i,j,k,TYPE_GZ,id)
         i2 = get_wrf_index(i,j+1,k,TYPE_GZ,id)
         fll(k) = (dym*( dxm*(wrf%dom(id)%phb(i,j,k)+x(i1)) + &
              dx*(wrf%dom(id)%phb(i+1,j,k)+x(i1+1))) + &
              dy*(dxm*(wrf%dom(id)%phb(i,j+1,k)+x(i2)) + &
              dx*(wrf%dom(id)%phb(i+1,j+1,k)+x(i2+1)) ))/gravity
      end do

   else

      fll(:) =  missing_r8

   endif

else

   if(i >= 1 .and. i <= wrf%dom(id)%var_size(1,TYPE_GZ) .and. &
        j >= 1 .and. j <= wrf%dom(id)%var_size(2,TYPE_GZ)) then

      do k = 1, wrf%dom(id)%var_size(3,TYPE_GZ)
         i1 = get_wrf_index(i,j,k,TYPE_GZ,id)
         fll(k) = (wrf%dom(id)%phb(i,j,k)+x(i1))/gravity
      end do

   else

      fll(:) =  missing_r8

   endif

endif

end subroutine get_model_height_profile_stag



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

end subroutine pert_model_state

!#######################################################

function map_to_sphere( longitude, latitude, id )
!
!   INPUT      : longitude, latitude (degrees)
!  OUTPUT      : map_to_sphere (radians), the angle by which to rotate
!              : the wrf wind vector to obtain the meridional and zonal winds.

real(r8), intent(in) :: longitude, latitude
integer,  intent(in) :: id
real(r8)             :: map_to_sphere

real(r8) :: cen_long, diff

if(wrf%dom(id)%map_proj .eq. 0) then  ! no projection
   map_to_sphere = 0.0_r8
end if

cen_long = wrf%dom(id)%stand_lon      ! Default
!!$cen_long = wrf%dom(id)%cen_lon     ! if stand_lon doesn't exist.

!!$cone = 1.0_r8
!!$if( wrf%dom(id)%map_proj .eq. 1) then    ; Lambert Conformal mapping
!!$   if( (fabs(wrf%dom(id)%truelat1 - wrf%dom(id)%truelat2) .gt. 0.1_r8) .and.  \
!!$      (fabs(wrf%dom(id)%truelat2 - 90.0_r8 )      .gt. 0.1_r8)       ) then
!!$      cone = 10^(cos(wrf%dom(id)%truelat1*deg2rad)) \
!!$      -10^(cos(wrf%dom(id)%truelat2*deg2rad))
!!$      cone = cone/(10^(tan(45.0_r8 -fabs(wrf%dom(id)%truelat1/2.0_r8)*deg2rad)) - \
!!$      10^(tan(45.0_r8 -fabs(wrf%dom(id)%truelat2/2.0_r8)*deg2rad))   )
!!$   else
!!$      cone = sin(fabs(wrf%dom(id)%truelat1)*deg2rad)
!!$   end if
!!$end if
!!$if(wrf%dom(id)%map_proj .eq. 2) then      ; polar stereographic
!!$   cone = 1.0_r8
!!$end if
!!$if(wrf%dom(id)%map_proj .eq. 3) then      ; Mercator
!!$   cone = 0.0_r8
!!$end if

diff = longitude - cen_long

if(diff .gt. 180.0_r8) then
   diff = diff - 360.0_r8
end if
if(diff .lt. -180.0_r8) then
   diff = diff + 360.0_r8
end if

!      map_to_sphere = diff * cone * deg2rad *sign(1.0_r8,latitude)
!!$map_to_sphere = diff

if(latitude .lt. 0.0_r8) then
   map_to_sphere = - diff * wrf%dom(id)%cone_factor * deg2rad
else
   map_to_sphere = diff * wrf%dom(id)%cone_factor * deg2rad
end if

!!$if(variable .eq. "umet") then
!!$   var = v*sin(map_to_sphere) + u*cos(map_to_sphere)
!!$else  ; must want vmet
!!$   var = v*cos(map_to_sphere) - u*sin(map_to_sphere)
!!$end if

end function map_to_sphere

end module model_mod

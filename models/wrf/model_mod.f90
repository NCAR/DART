module model_mod

! <next three lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$

! Assimilation interface for WRF model

!-----------------------------------------------------------------------
!
!     interface for WRF
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------
!-----------------------------------------------------------------------

use        types_mod, only : r8, deg2rad, rad2deg, missing_r
use time_manager_mod, only : time_type, set_time
use     location_mod, only : location_type, get_location, set_location, get_dist, &
                             LocationDims, LocationName, LocationLName, query_location
use    utilities_mod, only : file_exist, open_file, check_nml_error, close_file, &
                             error_handler, E_ERR

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
           TYPE_U, TYPE_V, TYPE_W, TYPE_GZ,   &
           TYPE_T, TYPE_MU,                   &
           TYPE_QV, TYPE_QC, TYPE_QR,         &
           pert_model_state

!  public stubs 

public     adv_1step,           &
           end_model,           &
           init_time,           &
           init_conditions,     &
           nc_write_model_atts, &
           nc_write_model_vars

!-----------------------------------------------------------------------
! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

!-----------------------------------------------------------------------
! model namelist parameters
!-----------------------------------------------------------------------

logical :: output_state_vector = .true.  ! output prognostic variables
integer :: num_moist_vars = 0            ! default value

namelist /model_nml/ output_state_vector, num_moist_vars


! Public definition of variable types

integer, parameter :: TYPE_U = 1, TYPE_V = 2, TYPE_W = 3, TYPE_GZ = 4, &
                      TYPE_T = 5,  TYPE_MU = 6,                        &
                      TYPE_QV = 7,  TYPE_QC = 8,  TYPE_QR = 9,         &
                      TYPE_QI = 10,  TYPE_QS = 11, TYPE_QG = 12


!-----------------------------------------------------------------------

real (kind=r8), PARAMETER    :: gas_constant = 287.04_r8
real (kind=r8), PARAMETER    :: gas_constant_v = 461.51_r8
real (kind=r8), PARAMETER    :: cp = 1004.0_r8
real (kind=r8), PARAMETER    :: t_kelvin = 273.15_r8
real (kind=r8), PARAMETER    :: gamma = 1.4_r8

real (kind=r8), PARAMETER    :: kappa = gas_constant / cp
real (kind=r8), PARAMETER    :: rd_over_rv = gas_constant / gas_constant_v
real (kind=r8), PARAMETER    :: rd_over_rv1 = 1.0 - rd_over_rv

!  Earth constants:
real (kind=r8), PARAMETER    :: gravity = 9.81_r8
real (kind=r8), PARAMETER    :: earth_radius = 6378.15_r8
INTEGER, PARAMETER           :: v_interp_p = 1, v_interp_h = 2
real (kind=r8), PARAMETER    :: ts0 = 300.0_r8
real (kind=r8), PARAMETER    :: ps0 = 100000.0_r8
!

!---- private data ----

TYPE wrf_static_data_for_dart

   integer  :: bt, bts, sn, sns, we, wes 
   real(r8) :: p_top, dx, dy, dt 
   integer  :: map_proj 
   real(r8) :: cen_lat,cen_lon,truelat1,truelat2,cone_factor,ycntr,psi1 

   integer  :: n_moist
   real(r8), dimension(:),     pointer :: znu, dn, dnw
   real(r8), dimension(:,:),   pointer :: mub, latitude, longitude
   real(r8), dimension(:,:),   pointer :: mapfac_m, mapfac_u, mapfac_v
   real(r8), dimension(:,:,:), pointer :: phb

   integer :: model_size, number_of_wrf_variables
   integer, dimension(:,:), pointer :: var_index
   integer, dimension(:,:), pointer :: var_size
   integer, dimension(:),   pointer :: var_type
   integer, dimension(:,:), pointer :: land
   character(len=8), dimension(:), pointer :: var_name

end type

type(wrf_static_data_for_dart) :: wrf

! ***** define other auxillary model data here as needed

!------------------------------------------------------------

contains

!#######################################################################

subroutine static_init_model()

! INitializes class data for WRF???

character (len = 80)      :: path
     integer              :: mode
     integer              :: ncid, bt_id, we_id, sn_id
     integer              :: io, ierr, unit

integer :: status
character (len=80) :: name
logical, parameter :: debug = .true.  
integer :: var_id, ind, i, map_proj 
integer, dimension(5) :: count, start, stride, map
real(r8)    :: zero_d(1)

real(r8) :: dx, dy, dt

real(r8) :: cen_lat, cen_lon, truelat1, truelat2
real(r8), dimension(:,:,:), pointer :: mub_test
real(r8), allocatable, dimension(:,:) :: temp

integer :: n_values

real(r8)    :: theta1,theta2,cell,cell2,psx
!----------

! Begin by reading the namelist input                                           
if(file_exist('input.nml')) then

   unit = open_file(file = 'input.nml', action = 'read')
   read(unit, nml = model_nml, iostat = io )
   ierr = check_nml_error(io, 'model_nml')
   call close_file(unit)                                                        

   wrf%n_moist = num_moist_vars

   if ( debug ) then
      write(*,'(''num_moist_vars = '',i3)')num_moist_vars
      write(*,'(''wrf%n_moist = '',i3)')wrf%n_moist
      if ( output_state_vector ) then
         write(*,*)'netcdf file in state vector format'
      else
         write(*,*)'netcdf file in prognostic vector format'
      endif
   endif
endif


mode = 0
call check( nf90_open('wrfinput', mode, ncid) )
if(debug) write(6,*) ' ncid is ',ncid

! get wrf grid dimensions

call check( nf90_inq_dimid(ncid, "bottom_top", bt_id) )
call check( nf90_inquire_dimension(ncid, bt_id, name, wrf%bt) )

call check( nf90_inq_dimid(ncid, "bottom_top_stag", bt_id) ) ! reuse bt_id, no harm
call check( nf90_inquire_dimension(ncid, bt_id, name, wrf%bts) )

call check( nf90_inq_dimid(ncid, "south_north", sn_id) )
call check( nf90_inquire_dimension(ncid, sn_id, name, wrf%sn) )

call check( nf90_inq_dimid(ncid, "south_north_stag", sn_id)) ! reuse sn_id, no harm
call check( nf90_inquire_dimension(ncid, sn_id, name, wrf%sns) )

call check( nf90_inq_dimid(ncid, "west_east", we_id) )
call check( nf90_inquire_dimension(ncid, we_id, name, wrf%we) )

call check( nf90_inq_dimid(ncid, "west_east_stag", we_id) )  ! reuse we_id, no harm
call check( nf90_inquire_dimension(ncid, we_id, name, wrf%wes) )

if(debug) then
   write(6,*) ' dimensions bt, sn, we are ',wrf%bt, wrf%sn, wrf%we
   write(6,*) ' staggered  bt, sn, we are ',wrf%bts,wrf%sns,wrf%wes
endif

! get meta data and static data we need

count  = 1
start  = 1
stride = 1
map    = 1

call check( nf90_get_att(ncid, nf90_global, 'DX', dx) )
wrf%dx = dx
if(debug) write(6,*) ' dx is ',dx

call check( nf90_get_att(ncid, nf90_global, 'DY', dy) )
wrf%dy = dy
if(debug) write(6,*) ' dy is ',dy

call check( nf90_get_att(ncid, nf90_global, 'DT', dt) )
wrf%dt = dt
if(debug) write(6,*) ' dt is ',dt

call netcdf_read_write_var( "P_TOP", ncid, var_id, zero_d,        &
                            start, count, stride, map, 'INPUT ', debug, 1  )
wrf%p_top = zero_d(1)
if(debug) write(6,*) ' p_top is ',wrf%p_top

call check( nf90_get_att(ncid, nf90_global, 'MAP_PROJ', map_proj) )


wrf%map_proj = map_proj

if(debug) write(6,*) ' map_proj is ',map_proj
call check( nf90_get_att(ncid, nf90_global, 'CEN_LAT', cen_lat) )

wrf%cen_lat = cen_lat
if(debug) write(6,*) ' cen_lat is ',wrf%cen_lat

call check( nf90_get_att(ncid, nf90_global, 'CEN_LON', cen_lon) )

wrf%cen_lon = cen_lon
call check( nf90_get_att(ncid, nf90_global, 'TRUELAT1', truelat1) )

wrf%truelat1 = truelat1
call check( nf90_get_att(ncid, nf90_global, 'TRUELAT2', truelat2) )

wrf%truelat2 = truelat2

if(debug) write(6,*) ' truelat2 is ',wrf%truelat2

if ( abs(truelat1-truelat2) .gt. 1.e-1 ) then
   theta1 = (90.0 - wrf%truelat1)*deg2rad
   theta2 = (90.0 - wrf%truelat2)*deg2rad
   wrf%cone_factor = (log(sin(theta1)) - log(sin(theta2))) &
        / (log(tan(theta1*0.5)) - log(tan(theta2*0.5)))
else
   wrf%cone_factor = sign(1.0,truelat1)*sin(truelat1 * deg2rad)
end if

write(unit=*, fmt='(2(a, e16.6))')'cone_factor  =', wrf%cone_factor  

IF (wrf%map_proj.EQ.1 .OR. wrf%map_proj.EQ.2) THEN
   IF(wrf%cen_lat.LT.0)THEN 
      wrf%psi1 = -(90.+wrf%truelat1)
   ELSE
      wrf%psi1 = 90.-wrf%truelat1            
   ENDIF
ELSE
   wrf%psi1 = 0.            
ENDIF

wrf%psi1 = deg2rad * wrf%psi1             

IF (wrf%map_proj.NE.3) THEN
   psx = (90.0 - wrf%cen_lat)*deg2rad    
   IF (wrf%map_proj.EQ.1) THEN
      cell  = earth_radius*SIN(wrf%psi1)/wrf%cone_factor
      cell2 = (TAN(psx/2.))/(TAN(wrf%psi1/2.))
   ENDIF
   IF (wrf%map_proj.EQ.2) THEN
      cell  = earth_radius*SIN(psx)/wrf%cone_factor
      cell2 = (1. + COS(wrf%psi1))/(1. + COS(psx))
   ENDIF
   wrf%ycntr = - cell*(cell2)**wrf%cone_factor
ENDIF
! -----FOR MERCATOR PROJECTION, THE PROJECTION IS TRUE AT LAT AT PHI1
IF (wrf%map_proj.EQ.3) THEN
   cell      = COS(wrf%cen_lat*deg2rad)/(1.0+SIN(wrf%cen_lat*deg2rad))
   wrf%ycntr = - earth_radius*COS(wrf%psi1)* log(cell)
ENDIF

! wrf%n_moist = 3  ! determined in namelist input 

!  get 1D (z) static data defining grid levels

count(1)  = wrf%bt
allocate(wrf%dn(1:wrf%bt))
call netcdf_read_write_var( "DN",ncid, var_id, wrf%dn,        &
                            start, count, stride, map, 'INPUT ', debug, 2  )
if(debug) write(6,*) ' dn ',wrf%dn

count(1)  = wrf%bt
allocate(wrf%znu(1:wrf%bt))
call netcdf_read_write_var( "ZNU",ncid, var_id, wrf%znu,        &
                            start, count, stride, map, 'INPUT ', debug, 2  )
if(debug) write(6,*) ' znu is ',wrf%znu
!
count(1)  = wrf%bt
allocate(wrf%dnw(1:wrf%bt))
call netcdf_read_write_var( "DNW",ncid, var_id, wrf%dnw,        &
                            start, count, stride, map, 'INPUT ', debug, 2  )
if(debug) write(6,*) ' dnw is ',wrf%dnw

!  get 2D (x,y) base state for mu, latitude, longitude

count(1)  = wrf%we
count(2)  = wrf%sn
count(3)  = 1
allocate(wrf%mub(1:wrf%we,1:wrf%sn))
call netcdf_read_write_var( "MUB",ncid, var_id, wrf%mub,        &
                            start, count, stride, map, 'INPUT ', debug, 3  )
if(debug) then
    write(6,*) ' corners of mub '
    write(6,*) wrf%mub(1,1),wrf%mub(wrf%we,1),  &
               wrf%mub(1,wrf%sn),wrf%mub(wrf%we,wrf%sn)
end if

allocate(wrf%longitude(1:wrf%we,1:wrf%sn))
call netcdf_read_write_var( "XLONG",ncid, var_id, wrf%longitude,        &
                            start, count, stride, map, 'INPUT ', debug, 3  )
allocate(wrf%latitude(1:wrf%we,1:wrf%sn))
call netcdf_read_write_var( "XLAT",ncid, var_id, wrf%latitude,        &
                            start, count, stride, map, 'INPUT ', debug, 3  )
allocate(wrf%land(1:wrf%we,1:wrf%sn), temp(1:wrf%we,1:wrf%sn))
call netcdf_read_write_var( "XLAND",ncid, var_id, temp,        &
                            start, count, stride, map, 'INPUT ', debug, 3  )
wrf%land = nint(temp)       ! coerce from float to integer ...
deallocate(temp)

if(debug) then
    write(6,*) ' corners of lat '
    write(6,*) wrf%latitude(1,1),wrf%latitude(wrf%we,1),  &
               wrf%latitude(1,wrf%sn),wrf%latitude(wrf%we,wrf%sn)
    write(6,*) ' corners of long '
    write(6,*) wrf%longitude(1,1),wrf%longitude(wrf%we,1),  &
               wrf%longitude(1,wrf%sn),wrf%longitude(wrf%we,wrf%sn)
end if

allocate(wrf%mapfac_m(1:wrf%we,1:wrf%sn))
call netcdf_read_write_var( "MAPFAC_M",ncid, var_id, wrf%mapfac_m,        &
                            start, count, stride, map, 'INPUT ', debug, 3  )

count(1)  = wrf%wes
count(2)  = wrf%sn
allocate(wrf%mapfac_u(1:wrf%wes,1:wrf%sn))
call netcdf_read_write_var( "MAPFAC_U",ncid, var_id, wrf%mapfac_u,        &
                            start, count, stride, map, 'INPUT ', debug, 4  )

count(1)  = wrf%we
count(2)  = wrf%sns
allocate(wrf%mapfac_v(1:wrf%we,1:wrf%sns))
call netcdf_read_write_var( "MAPFAC_V",ncid, var_id, wrf%mapfac_V,        &
                            start, count, stride, map, 'INPUT ', debug, 4  )

! get 3D base state geopotential

count(1)  = wrf%we
count(2)  = wrf%sn
count(3)  = wrf%bts
allocate(wrf%phb(1:wrf%we,1:wrf%sn,1:wrf%bts))
call netcdf_read_write_var( "PHB",ncid, var_id, wrf%phb,        &
                            start, count, stride, map, 'INPUT ', debug, 4  )
if(debug) then
    write(6,*) ' corners of phb '
    write(6,*) wrf%phb(1,1,1),wrf%phb(wrf%we,1,1),  &
               wrf%phb(1,wrf%sn,1),wrf%phb(wrf%we,wrf%sn,1)
    write(6,*) wrf%phb(1,1,wrf%bts),wrf%phb(wrf%we,1,wrf%bts),  &
               wrf%phb(1,wrf%sn,wrf%bts),wrf%phb(wrf%we,wrf%sn,wrf%bts)
end if

! close data file, we have all we need

call check( nf90_close(ncid) )

!  build the map into the 1D DART vector for WRF data

  wrf%number_of_wrf_variables = 6 + wrf%n_moist
  allocate(wrf%var_type(12), wrf%var_name(12))    ! always define all 12, may use less
  wrf%var_type(1)  = TYPE_U
  wrf%var_type(2)  = TYPE_V
  wrf%var_type(3)  = TYPE_W
  wrf%var_type(4)  = TYPE_GZ
  wrf%var_type(5)  = TYPE_T
  wrf%var_type(6)  = TYPE_MU
  wrf%var_type(7)  = TYPE_QV
  wrf%var_type(8)  = TYPE_QC
  wrf%var_type(9)  = TYPE_QR
  wrf%var_type(10) = TYPE_QI
  wrf%var_type(11) = TYPE_QS
  wrf%var_type(12) = TYPE_QG

  ! These are going to be used as the netCDF variable names

  wrf%var_name( 1)  = "U"            ! The first 6 are the "base variables"
  wrf%var_name( 2)  = "V"            ! We always use AT LEAST these.
  wrf%var_name( 3)  = "W"
  wrf%var_name( 4)  = "GZ"
  wrf%var_name( 5)  = "T"
  wrf%var_name( 6)  = "MU"
  wrf%var_name( 7)  = "QVAPOR"       ! wrf%nmoist = 1
  wrf%var_name( 8)  = "QCLOUD"       ! wrf%nmoist = 2
  wrf%var_name( 9)  = "QRAIN"        ! wrf%nmoist = 3
  wrf%var_name(10)  = "QICE"         ! wrf%nmoist = 4
  wrf%var_name(11)  = "QSNOW"        ! wrf%nmoist = 5
  wrf%var_name(12)  = "QGRAP"        ! wrf%nmoist = 6

  allocate(wrf%var_index(2,6 + wrf%n_moist)) ! indices into 1D array
  allocate(wrf%var_size(3,6 + wrf%n_moist)) ! dimension of variables
  
  ind = 1                         ! *** u field ***
  wrf%var_size(1,ind) = wrf%wes
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = wrf%bt
  wrf%var_index(1,ind) = 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +  &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** v field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sns
  wrf%var_size(3,ind) = wrf%bt
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +   &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** w field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = wrf%bts
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +    &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** geopotential field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = wrf%bts
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +   &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** theta field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = wrf%bt
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +   &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** mu field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = 1
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +   &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  do i = 1, wrf%n_moist
    ind = ind + 1                   ! *** moisture field ***
    wrf%var_size(1,ind) = wrf%we
    wrf%var_size(2,ind) = wrf%sn
    wrf%var_size(3,ind) = wrf%bt
    wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
    wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +   &
                         wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)
  enddo

  wrf%model_size = wrf%var_index(2,wrf%number_of_wrf_variables)
  write(6,*) ' wrf model size is ',wrf%model_size

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'static_init_model', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end subroutine static_init_model


!**********************************************************************

subroutine netcdf_read_write_var( variable, ncid, var_id, var,          &
                                  start, count, stride, map, in_or_out, debug, ndims )

integer :: ncid, var_id, ndims
real(r8), dimension(ndims) :: var
character (len=6) :: in_or_out
integer, dimension(ndims) :: start, count, stride, map
integer :: status
character (len=*) :: variable
logical :: debug
character (len=129) :: error_string

if(debug) write(6,*) ' var for io is ',variable
call check(  nf90_inq_varid(ncid, variable, var_id) )
if(debug) write(6,*) variable, ' id = ',var_id

if( in_or_out(1:5) == "INPUT" ) then
  if(debug) write(6,*) ' call netcdf read ', ncid, var_id
  call check( nf90_get_var(ncid, var_id, var, start, count, stride ) )
  if(debug) write(6,*) ' returned netcdf read '
else if( in_or_out(1:6) == "OUTPUT" ) then
  call check( nf90_put_var(ncid, var_id, var, start, count, stride, map) )
else
  write(error_string,*)' unknown IO function for var_id ',var_id, in_or_out
  call error_handler(E_ERR,'netcdf_read_write_var', &
       error_string, source, revision,revdate)
end if

contains
  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus
    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'netcdf_read_write_var', &
       trim(nf90_strerror(istatus)), source, revision, revdate)
  end subroutine check
end subroutine netcdf_read_write_var

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
! Returns the the time step of the model. In the long run should be replaced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

! Time_step_atmos is global static storage

! need to translate from wrf model timestep (in seconds) to
! DART time increment

get_model_time_step =  set_time(nint(wrf%dt), 0)

end function get_model_time_step

!#######################################################################

subroutine get_state_meta_data(index_in, location, var_type)
!---------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?
!  SHOULD THIS ALSO RETURN THE TYPE OF THIS VARIABLE???
! YES NEED TO RETURN VARIABLE TYPE HERE

integer, intent(in) :: index_in
! Temporary kluge of location type
!integer, intent(in) :: location
type(location_type), intent(out) :: location
integer, intent(out), optional :: var_type
integer :: var_type_out

integer :: index, ip, jp, kp
integer :: nz, ny, nx
logical :: var_found
real(r8)    :: lon, lat, lev

integer :: i, number_of_wrf_variables
logical, parameter :: debug = .false.  

if(debug) then
  write(6,*) ' in get_state_meta_data '
  write(6,*) ' index in is ',index_in
endif

index = index_in
var_found = .false.

!  first find var_type

  if(debug) then
    do i=1, wrf%number_of_wrf_variables
      write(6,*) ' i, var_type(i) ',i,wrf%var_type(i)
    enddo
  endif

   i = 0
   do while (.not. var_found)
     i = i + 1
     if( (index .ge. wrf%var_index(1,i) ) .and.  &
         (index .le. wrf%var_index(2,i) )       )  then
       var_found = .true.
       var_type_out = wrf%var_type(i)
       index = index - wrf%var_index(1,i) + 1
     end if
     if(i .gt. wrf%number_of_wrf_variables) then
       write(6,*) ' index out of range in get_state_meta_data '
       write(6,*) ' index is ',index
       write(6,*) ' size of vector ',wrf%model_size
       stop
     end if
   end do

!  now find i,j,k location.
!  index has been normalized such that it is relative to 
!  array starting at (1,1,1)

   nx = wrf%var_size(1,i)
   ny = wrf%var_size(2,i)
   nz = wrf%var_size(3,i)

   kp = 1 + (index-1)/(nx*ny)
   jp = 1 + (index - (kp-1)*nx*ny - 1)/nx
   ip = index - (kp-1)*nx*ny - (jp-1)*nx

  if(debug) write(6,*) ' ip, jp, kp for index ',ip,jp,kp,index
  if(debug) write(6,*) ' Var type: ',var_type_out

! find lat and long, must
! correct for possible u or v staggering in x, y

   if (var_type_out == TYPE_U) then

     if(ip == 1) then

        lon = wrf%longitude(1,jp) - 0.5*(wrf%longitude(2,jp)-wrf%longitude(1,jp))
        if (wrf%longitude(2,jp) < wrf%longitude(1,jp)) lon = lon - 180.0
        lat = wrf%latitude(1,jp)  - 0.5*(wrf%latitude(2,jp)-wrf%latitude(1,jp))

     else if(ip == nx) then

        lon = wrf%longitude(nx-1,jp) + 0.5*(wrf%longitude(nx-1,jp)-wrf%longitude(nx-2,jp))
        if (wrf%longitude(nx-1,jp) < wrf%longitude(nx-2,jp)) lon = lon + 540.0
        lat = wrf%latitude(nx-1,jp)  + 0.5*(wrf%latitude(nx-1,jp)-wrf%latitude(nx-2,jp))

     else

        lon = 0.5*(wrf%longitude(ip,jp)+wrf%longitude(ip-1,jp))
        if (wrf%longitude(ip,jp) < wrf%longitude(ip-1,jp)) lon = lon + 180.0
        lat = 0.5*(wrf%latitude(ip,jp) +wrf%latitude(ip-1,jp))

     end if

   else if (var_type_out == TYPE_V) then

     if(jp == 1)  then
       lon = wrf%longitude(ip,1) - 0.5*(wrf%longitude(ip,2)-wrf%longitude(ip,1))
       if (wrf%longitude(ip,2) < wrf%longitude(ip,1)) lon = lon - 180.0
       lat = wrf%latitude(ip,1)  - 0.5*(wrf%latitude(ip,2)-wrf%latitude(ip,1))
     else if(jp == ny) then
       lon = wrf%longitude(ip,ny-1) + 0.5*(wrf%longitude(ip,ny-1)-wrf%longitude(ip,ny-2))
       if (wrf%longitude(ip,ny-1) < wrf%longitude(ip,ny-2)) lon = lon + 540.0
       lat = wrf%latitude(ip,ny-1)  + 0.5*(wrf%latitude(ip,ny-1)-wrf%latitude(ip,ny-2))
     else 
       lon = 0.5*(wrf%longitude(ip,jp)+wrf%longitude(ip,jp-1))
       if (wrf%longitude(ip,jp) < wrf%longitude(ip,jp-1)) lon = lon + 180.0
       lat = 0.5*(wrf%latitude(ip,jp) +wrf%latitude(ip,jp-1))
     end if

   else  ! regularily staggered variable

     lon = wrf%longitude(ip,jp)
     lat = wrf%latitude(ip,jp)

   end if

   if (lon < 0.0_r8) lon = lon + 360.0

   if (lon > 360.0_r8) lon = lon - 360.0

   lev = float(kp)

   if(debug) write(6,*) 'lon, lat, lev: ',lon, lat, lev
          
   location = set_location(lon, lat, lev, 1)

   if(present(var_type)) var_type = var_type_out

end subroutine get_state_meta_data

!#######################################################################

function model_interpolate(x, location, obs_kind)
!!!function model_interpolate(x, lon, lat, level, type)

logical, parameter :: debug = .false.  
real(r8) :: model_interpolate
real(r8), intent(in) :: x(:)
type(location_type), intent(in) :: location
integer, intent(in) :: obs_kind
real (r8)           :: xloc, yloc, zloc, xyz_loc(3)
integer             :: i, j, k, i1,i2,i3
real                :: dx,dy,dxm,dym
real(r8), dimension (wrf%bt)   :: v_h, v_p
real(r8), dimension (wrf%bt  ) :: fld
real(r8)                       :: p1,p2,p3,p4,a1

xyz_loc(:) = get_location(location)
call llxy(xyz_loc(1),xyz_loc(2),xloc,yloc)
call toGrid(xloc,wrf%we,i,dx,dxm)
call toGrid(yloc,wrf%sn,j,dy,dym)
!  get model pressure profile
call get_model_pressure_profile(i,j,dx,dy,dxm,dym,wrf%bt,x,v_p,p1,p2,p3,p4)
!  get model height profile
call get_model_height_profile(i,j,dx,dy,dxm,dym,wrf%bt,x,v_h)

  if(nint(query_location(location,'which_vert')) == 1) then

!  If obs is by level
     zloc = xyz_loc(3)
     if(debug) print*,' obs is by location and zloc =',zloc
  else if(nint(query_location(location,'which_vert')) == 2) then
     ! get pressure vertical co-ordinate
     call to_zk(xyz_loc(3), v_p, wrf%bt, v_interp_p,zloc) 
     if(debug.and.obs_kind /=3) print*,' obs is by pressure and zloc =',zloc
!     print*,'model pressure profile'
!     print*,v_p
  else if(nint(query_location(location,'which_vert')) == 3) then
     ! get height vertical co-ordinate
     call to_zk(xyz_loc(3), v_h, wrf%bt, v_interp_h,zloc) 
     if(debug) print*,' obs is by height and zloc =',zloc
  else if(nint(query_location(location,'which_vert')) == -1) then
     ! get height vertical co-ordinate
     if(debug) print*,' obs is surface pressure   = ', xyz_loc(3)
  else
     print*,' wrong option for which_vert'
     stop

  end if
! Get the desired field to be interpolated
 if( obs_kind == 1 ) then                 ! U
  do k=1,wrf%bt
     i1 = get_wrf_index(i,j,k,TYPE_U)   
     i2 = get_wrf_index(i,j+1,k,TYPE_U)   
    
     fld(k) = dym*(dxm*(x(i1) + x(i1+1))*.5 + dx*(x(i1+1)+x(i1+2))*0.5) + &
               dy*(dxm*(x(i2) + x(i2+1))*.5 + dx*(x(i2+1)+x(i2+2))*0.5)     
!     if(debug) print*,k,' model u profile ',fld(k)
  end do 
 else if( obs_kind == 2) then                 ! V
  do k=1,wrf%bt
     i1 = get_wrf_index(i,j,k,TYPE_V) 
     i2 = get_wrf_index(i,j+1,k,TYPE_V) 
     i3 = get_wrf_index(i,j+2,k,TYPE_V) 
!
     fld(k) = dym*(dxm*(x(i1) + x(i2))*.5 + dx*(x(i1+1)+x(i2+1))*0.5) + &
               dy*(dxm*(x(i2) + x(i3))*.5 + dx*(x(i2+1)+x(i3+1))*0.5)     
     if(debug) print*,k,' model v profile ',fld(k)
  end do 
 else if( obs_kind == 4 ) then                ! T  
  do k=1,wrf%bt
     i1 = get_wrf_index(i,j,k,TYPE_T) 
     i2 = get_wrf_index(i,j+1,k,TYPE_T) 
     a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
     fld(k) = (ts0 + a1)*(v_p(k)/ps0)**kappa
     if(debug) print*,k,' model temp profile ',fld(k)
  end do 
 else if( obs_kind == 5) then                ! Q
  do k=1,wrf%bt

     i1 = get_wrf_index(i,j,k,TYPE_QV)   
     i2 = get_wrf_index(i,j+1,k,TYPE_QV)   
     a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
     fld(k) = a1 /(1.0 + a1)      
     if(debug) print*,k,' model q profile ',fld(k)
  end do 
! Do 1D interpolation
 else if( obs_kind == 3)then
!  surfacce pressure
  model_interpolate = wrf%p_top                                  +&
                     dym*(dxm*wrf%mub(i,j)+dx*wrf%mub(i+1,j))    +&
                     dy *(dxm*wrf%mub(i,j+1)+dx*wrf%mub(i+1,j+1))+&
                     dym*(dxm*p1 + dx*p2) + dy*(dxm*p3 + dx*p4)
  if(debug) print*,' for sfc model val =',model_interpolate
 else
  print*,'Do not know what to do for this observation kind type ',obs_kind
  stop
 end if
  if(obs_kind /= 3) call Interp_lin_1D(fld, wrf%bt, zloc, model_interpolate)
if(debug) print*,' interpolated value= ',model_interpolate

end function model_interpolate

!#######################################################################

function get_val(x, lon_index, lat_index, level, obs_kind)

real(r8) :: get_val
real(r8), intent(in) :: x(:)
integer, intent(in) :: lon_index, lat_index, level, obs_kind

get_val = x(get_wrf_index(lon_index,lat_index,level,obs_kind))

end function get_val

!#######################################################################

subroutine model_get_close_states(o_loc, radius, number, indices, dist)

type(location_type), intent(in) :: o_loc
real(r8), intent(in) :: radius
integer, intent(out) :: number, indices(:)
real(r8), intent(out) :: dist(:)

real(r8) :: loc_array(3), o_lon, o_lat
integer :: num, max_size
integer, allocatable :: lon_ind(:), lat_ind(:)
real(r8), allocatable :: close_dist(:)

integer :: u_pts, v_pts, p_pts
integer :: i,j,k,indmax, num_total, ii, jj

! Number found starts at 0
number = 0

! Num of close horizontal grid points starts at 0, too
num = 0

indmax = size(indices)

! For now, just allocate enough space for all grid points, may want
! to make this smaller at some point for big models.

! we're allocating enough space for all u, v, and p points
! on a horizontal plane

max_size = wrf%we*wrf%sn + wrf%wes*wrf%sn + wrf%we*(wrf%sns)
allocate(lon_ind(max_size), lat_ind(max_size), close_dist(max_size))

! Look for close grid points on the horizontal grid (2D)
call grid_close_states( o_loc, wrf%latitude, wrf%longitude, radius,  &
                        num, lon_ind, lat_ind, close_dist, u_pts, v_pts, p_pts )

! next check vertical spacing, find out how many points
! we have in 3D

   ! w+gz + t + moist + mu
max_size = p_pts*( 2*(wrf%bts) + (1 + wrf%n_moist)*wrf%bts)

   ! u and v
max_size = max_size + (u_pts+v_pts)*wrf%bt
!allocate( dist_3d(max_size), lat_3d(max_size),   &
!          long_3d(max_size), h_3d(max_size)     )


num_total = 0

! start with p_pts (t + wrf%n_moist variables)

do k = 1, wrf%bt
   do i = 1, p_pts

      ii = lon_ind(i)
      jj = lat_ind(i)

      num_total = num_total + 1
      if(num_total <= indmax) then
         indices(num_total) = get_wrf_index(ii,jj,k,type_t)
         dist(num_total) = close_dist(i)
      end if
      if( wrf%n_moist == 1) then
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qv)
            dist(num_total) = close_dist(i)
         end if
      else if( wrf%n_moist == 3) then
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qv)
            dist(num_total) = close_dist(i)
         end if
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qc)
            dist(num_total) = close_dist(i)
         end if
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qr)
            dist(num_total) = close_dist(i)
         end if
      else if( wrf%n_moist == 6) then
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qv)
            dist(num_total) = close_dist(i)
         end if
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qc)
            dist(num_total) = close_dist(i)
         end if
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qr)
            dist(num_total) = close_dist(i)
         end if
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qi)
            dist(num_total) = close_dist(i)
         end if
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qs)
            dist(num_total) = close_dist(i)
         end if
         num_total = num_total + 1
         if(num_total <= indmax) then
            indices(num_total) = get_wrf_index(ii,jj,k,type_qg)
            dist(num_total) = close_dist(i)
         end if
      end if

   enddo
enddo

! now w_pts (w and gz)

do k = 1, wrf%bts
   do i = 1, p_pts

      ii = lon_ind(i)
      jj = lat_ind(i)

      num_total = num_total + 1
      if(num_total <= indmax) then
         indices(num_total) = get_wrf_index(ii,jj,k,type_w)
         dist(num_total) = close_dist(i)
      end if
      num_total = num_total + 1
      if(num_total <= indmax) then
         indices(num_total) = get_wrf_index(ii,jj,k,type_gz)
         dist(num_total) = close_dist(i)
      end if

   enddo
enddo

! now mu_pts (surface pressure)
!  we're assuming that vertical location does not matter?

k = 1
do i = 1, p_pts
   ii = lon_ind(i)
   jj = lat_ind(i)

   num_total = num_total + 1
   if(num_total <= indmax) then
      indices(num_total) = get_wrf_index(ii,jj,k,type_mu)
      dist(num_total) = close_dist(i)
   end if
enddo

! now u_pts 

do k = 1, wrf%bt
   do i = p_pts+1, p_pts+u_pts

      ii = lon_ind(i)
      jj = lat_ind(i)

      num_total = num_total + 1
      if(num_total <= indmax) then
         indices(num_total) = get_wrf_index(ii,jj,k,type_u)
         dist(num_total) = close_dist(i)
      end if

   enddo
enddo

! last -> v_pts 

do k = 1, wrf%bt
   do i = p_pts+u_pts+1, num

      ii = lon_ind(i)
      jj = lat_ind(i)

      num_total = num_total + 1
      if(num_total <= indmax) then
         indices(num_total) = get_wrf_index(ii,jj,k,type_v)
         dist(num_total) = close_dist(i)
      end if

   enddo
enddo

deallocate(lon_ind, lat_ind, close_dist)

number = num_total

end subroutine model_get_close_states

!-----------------------------------------------------------------

function get_wrf_index( i,j,k,var_type )

integer, intent(in) :: i,j,k,var_type
integer :: get_wrf_index
integer :: in

integer :: ii

  in = 0
  do ii = 1, wrf%number_of_wrf_variables 
    if(var_type == wrf%var_type(ii) ) in = ii
  enddo

  get_wrf_index = wrf%var_index(1,in)-1 +   &
    i + wrf%var_size(1,in)*(j-1) + (wrf%var_size(1,in)*wrf%var_size(2,in))*(k-1)

end function get_wrf_index

!-----------------------------------------------------------------

subroutine grid_close_states( o_loc, lat, lon, radius_in, num,  &
                              close_lon_ind, close_lat_ind,     &
                              close_dist, u_pts, v_pts, p_pts     )

! Finds close state points from a particular grid for the WRF model

type(location_type), intent(in) :: o_loc
integer, intent(inout) :: num
integer, intent(out) :: close_lon_ind(:), close_lat_ind(:)
integer, intent(out) :: u_pts, v_pts, p_pts
real(r8), intent(in) :: radius_in
real(r8), intent(out) :: close_dist(:)
real(r8), dimension(3) :: loc_array
real(r8) :: o_lon, o_lat

real(r8) :: radius
real(r8), intent(in) :: lat(:,:), lon(:,:)
real(r8) :: rad, radn, dxr, dyr, sdx, sdy, gdist
integer :: i_closest, j_closest, ixmin, jymin, ixmax, jymax

integer :: i, j, n, m
real(r8), parameter :: r_earth = 6.37e+06 ! earth radius in meters
logical, parameter :: debug= .false.  

type(location_type) :: loc

if(debug) write(6,*) ' in grid_close_states '

! use meaningful units for radius -- convert radians to meters 
radius = radius_in*r_earth
if(debug) write(6,*) ' radius in grid_close_states is ',radius

! Get the lat and lon from the loc
loc_array = get_location(o_loc)

o_lon = loc_array(1)
o_lat = loc_array(2)

if(debug) write(6,*) ' observations long and lat ',o_lon,o_lat

! Get index to closest lat and lon for this observation

n = size( lat, 1 )
m = size( lat, 2 )
!<<<alain

!===alainrad = (o_lon-lon(1,1))**2 + (o_lat-lat(1,1))**2
 rad = get_dist_wrf(1,1,0, type_t, o_loc)
i_closest = 1
j_closest = 1

! brute force search
do j=1,m
do i=1,n
!===alain  radn = (o_lon-lon(i,j))**2 + (o_lat-lat(i,j))**2
   radn = get_dist_wrf(i,j,0, type_t, o_loc)
  if( radn .lt. rad ) then
    rad = radn
    i_closest = i
    j_closest = j
  end if
enddo
enddo

!>>>alain

if(debug) write(6,*) ' closest wrf long and lat is ',i_closest,j_closest,lon(i_closest,j_closest),lat(i_closest,j_closest)
if(debug) write(6,*) ' radius is ',radius

! define box edges for radius check
dxr = 1 + radius/wrf%dx  !  radius in multiples of dx
dyr = 1 + radius/wrf%dy  !  radius in multiples of dy

if(debug) write(6,*) ' dxr, dyr in grid_close_states ',dxr,dyr


  j = j_closest
  ixmin = max(1,i_closest - 1)
  sdx   = 1./wrf%mapfac_u(i_closest, j)
  do while( sdx .lt. dxr )
    ixmin = max(1,ixmin - 1)
    sdx = sdx + 1./wrf%mapfac_u(ixmin + 1, j)
    if(ixmin <= 1) sdx = 1.1*dxr
  enddo

  ixmax = min(wrf%we,i_closest + 1)
  sdx   = 1./wrf%mapfac_u(i_closest+1, j)
  do while( sdx .lt. dxr )
    ixmax = min(wrf%we,ixmax + 1)
    sdx = sdx + 1./wrf%mapfac_u(ixmax, j)
    if(ixmax >= wrf%we) sdx = 1.1*dxr
  enddo

  jymin = max(1,j_closest - 1)
  i = i_closest
  sdy   = 1./wrf%mapfac_u(i_closest, jymin)
  do while( sdy .lt. dyr )
    jymin = max(1,jymin - 1)
    sdy = sdy + 1./wrf%mapfac_u(i, jymin + 1)
    if(jymin <= 1) sdy = 1.1*dyr
  enddo

  jymax = min(wrf%sn,j_closest + 1)
  sdy   = 1./wrf%mapfac_u(i, j_closest+1)
  do while( sdy .lt. dyr )
    jymax = min(wrf%sn,jymax + 1)
    sdy = sdy + 1./wrf%mapfac_u(i, jymax)
    if(jymax >= wrf%sn) sdy = 1.1*dyr
  enddo

  if(debug) then
    write(6,*) ' ixmin, ixmax, jymin, jymax are '
    write(6,*) ixmin, ixmax, jymin, jymax
  endif

!  we have bounding box, get and check distances.
!  first, convert radius back to radians

  radius = radius_in
  num = 0

  do j = jymin, jymax
  do i = ixmin, ixmax
    gdist = get_dist_wrf(i,j,0, type_t, o_loc)
    if ( gdist <= radius ) then
      num = num + 1
      close_lon_ind(num) = i
      close_lat_ind(num) = j
      close_dist(num) = gdist
      if(debug) write(6,*) ' p pt ',num,i,j,gdist
    end if
  enddo
  enddo

  p_pts = num

! check distance for u points, expand box so that 
! we don't leave possible points out of check

  do j = jymin, jymax
  do i = max(1,ixmin-1), ixmax+1

     gdist = get_dist_wrf(i,j,0, type_u, o_loc)
     if ( gdist <= radius ) then
       num = num + 1
       close_lon_ind(num) = i
       close_lat_ind(num) = j
       close_dist(num) = gdist
      if(debug) write(6,*) ' u pt ',num,i,j,gdist
     end if
   enddo
   enddo

    u_pts = num - p_pts

   do j = max(1,jymin-1), jymax+1
   do i = ixmin, ixmax

      gdist = get_dist_wrf(i,j,0, type_v, o_loc)
      if ( gdist <= radius ) then
        num = num + 1
        close_lon_ind(num) = i
        close_lat_ind(num) = j
        close_dist(num) = gdist
      if(debug) write(6,*) ' v pt ',num,i,j,gdist
      end if

   enddo
   enddo

     v_pts = num - u_pts - p_pts

end subroutine grid_close_states

!***********************************************************************

function get_dist_wrf( i,j,k,var_type,o_loc )

type(location_type), intent(in) :: o_loc
integer, intent(in) :: i,j,k,var_type

type(location_type) :: loc

real(r8) :: get_dist_wrf

real(r8) :: long, lat, lev



!  get distance for input var_type
   if(var_type == type_u) then

     if (i == 1) then
       long = wrf%longitude(1,j) - 0.5*(wrf%longitude(2,j)-wrf%longitude(1,j))
       if (wrf%longitude(2,j) < wrf%longitude(1,j)) long = long - 180.0
       lat  = wrf%latitude(1,j)  - 0.5*(wrf%latitude(2,j)-wrf%latitude(1,j))
     else if (i == wrf%wes) then
       long = wrf%longitude(i-1,j) + 0.5*(wrf%longitude(i-1,j)-wrf%longitude(i-2,j))
       if (wrf%longitude(i-1,j) < wrf%longitude(i-2,j)) long = long + 540.0
       lat  = wrf%latitude(i-1,j)  + 0.5*(wrf%latitude(i-1,j)-wrf%latitude(i-2,j))
     else
       long = 0.5*(wrf%longitude(i,j)+wrf%longitude(i-1,j))
       if (wrf%longitude(i,j) < wrf%longitude(i-1,j)) long = long + 180.0
       lat  = 0.5*(wrf%latitude(i,j) +wrf%latitude(i-1,j))
     end if

   else if( var_type == type_v) then

     if (j == 1) then
       long = wrf%longitude(i,1) - 0.5*(wrf%longitude(i,2)-wrf%longitude(i,1))
       if (wrf%longitude(i,2) < wrf%longitude(i,1)) long = long - 180.0
       lat  = wrf%latitude(i,1)  - 0.5*(wrf%latitude(i,2)-wrf%latitude(i,1))
     else if (j == wrf%sns) then
       long = wrf%longitude(i,j-1) + 0.5*(wrf%longitude(i,j-1)-wrf%longitude(i,j-2))
       if (wrf%longitude(i,j-1) < wrf%longitude(i,j-2)) long = long + 540.0
       lat  = wrf%latitude(i,j-1)  + 0.5*(wrf%latitude(i,j-1)-wrf%latitude(i,j-2))
     else
       long = 0.5*(wrf%longitude(i,j)+wrf%longitude(i,j-1))
       if (wrf%longitude(i,j) < wrf%longitude(i,j-1)) long = long + 180.0
       lat  = 0.5*(wrf%latitude(i,j) +wrf%latitude(i,j-1))
     end if

   else

      long = wrf%longitude(i,j)
      lat = wrf%latitude(i,j)
   end if
   if( (var_type == type_w ) .or. &
        (var_type == type_gz) .or. &
        (var_type == type_mu)     ) then

      lev = float(k)-1.0

   else

      lev = float(k)-0.5

   endif

   if (long < 0.0_r8) long = long + 360.0

   if (long > 360.0_r8) long = long - 360.0

   loc = set_location(long,lat,lev, 1)
   get_dist_wrf = get_dist(loc,o_loc)

   end function get_dist_wrf

!---------------------------------

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
integer :: weDimID, weStagDimID
integer :: snDimID, snStagDimID
integer :: btDimID, btStagDimId, MemberDimID
integer :: DNVarID, ZNUVarID, DNWVarID, phbVarID
integer :: MubVarID, LonVarID, LatVarID, ilevVarID, XlandVarID 
integer :: MapFacMVarID, MapFacUVarID, MapFacVVarID
integer :: UVarID, VVarID, WVarID, PHVarID, tVarID, MuVarID
integer :: QVVarID, QCVarID, QRVarID
integer :: i
integer :: ncid            ! for wrfinput reading

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
  write(*,*)'ERROR: nc_write_model_atts: Time      dimension is ',TimeDimID
  write(*,*)'ERROR: nc_write_model_atts: unlimited dimension is ',unlimitedDimId
  write(*,*)'ERROR: they must be the same.'
  stop
endif

!-----------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-----------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=wrf%model_size, dimid = StateVarDimID))

!-----------------------------------------------------------------
! Write Global Attributes 
!-----------------------------------------------------------------

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source",source))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "DX", wrf%dx))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "DY", wrf%dy))

! how about namelist input? might be nice to save ...

!-----------------------------------------------------------------
! Define the dimensions IDs
!-----------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="west_east",        &
          len = wrf%we,  dimid = weDimID))
call check(nf90_def_dim(ncid=ncFileID, name="west_east_stag",   &
          len = wrf%wes, dimid = weStagDimID))
call check(nf90_def_dim(ncid=ncFileID, name="south_north",      &
          len = wrf%sn,  dimid = snDimID))
call check(nf90_def_dim(ncid=ncFileID, name="south_north_stag", &
          len = wrf%sns, dimid = snStagDimID))
call check(nf90_def_dim(ncid=ncFileID, name="bottom_top",       &
          len = wrf%bt,  dimid = btDimID))
call check(nf90_def_dim(ncid=ncFileID, name="bottom_top_stag",  &
          len = wrf%bts, dimid = btStagDimID))

!-----------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!-----------------------------------------------------------------

if ( output_state_vector ) then

   !-----------------------------------------------------------------
   ! Create attributes for the state vector 
   !-----------------------------------------------------------------

   call check(nf90_def_var(ncFileID, name="west_east", xtype=nf90_double, &
             dimids = weDimID, varid=LonVarID) )

   call check(nf90_put_att(ncFileID, LonVarID, "long_name", "longitude"))
   call check(nf90_put_att(ncFileID, LonVarID, "cartesian_axis", "X"))
   call check(nf90_put_att(ncFileID, LonVarID, "units", "degrees_east"))
   call check(nf90_put_att(ncFileID, LonVarID, "valid_range", (/ -180.0_r8, 180.0_r8 /)))
   call check(nf90_put_att(ncFileID, LonVarID, "description", "LONGITUDE, WEST IS NEGATIVE"))

   ! Latitudes
   call check(nf90_def_var(ncFileID, name="south_north", xtype=nf90_double, &
             dimids = snDimID, varid=LatVarID) ) 

   call check(nf90_put_att(ncFileID, LatVarID, "long_name", "latitude"))
   call check(nf90_put_att(ncFileID, LatVarID, "cartesian_axis", "Y"))
   call check(nf90_put_att(ncFileID, LatVarID, "units", "degrees_north"))
   call check(nf90_put_att(ncFileID, LatVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)))
   call check(nf90_put_att(ncFileID, LatVarID, "description", "LATITUDE, SOUTH IS NEGATIVE"))

   ! grid levels
   call check(nf90_def_var(ncFileID, name="bottom_top", xtype=nf90_double, &
             dimids = btDimID, varid=ilevVarID) )

   call check(nf90_put_att(ncFileID, ilevVarID, "long_name", "level"))
   call check(nf90_put_att(ncFileID, ilevVarID, "cartesian_axis", "Z"))
   call check(nf90_put_att(ncFileID, ilevVarID, "units", "at this point, indexical"))

   ! Define the state vector coordinate variable

   call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID))

   call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"))
   call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical") )
   call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, wrf%model_size /)))

   ! Define the actual state vector

   call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), varid=StateVarID))
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"))
   call check(nf90_put_att(ncFileID, StateVarId, "U_units","m/s"))
   call check(nf90_put_att(ncFileID, StateVarId, "V_units","m/s"))
   call check(nf90_put_att(ncFileID, StateVarId, "W_units","m/s"))
   call check(nf90_put_att(ncFileID, StateVarId, "GZ_units","m2/s2"))
   call check(nf90_put_att(ncFileID, StateVarId, "T_units","K"))
   call check(nf90_put_att(ncFileID, StateVarId, "MU_units","Pa"))
   call check(nf90_put_att(ncFileID, StateVarId, "QV_units","kg/kg"))
   call check(nf90_put_att(ncFileID, StateVarId, "QC_units","kg/kg"))
   call check(nf90_put_att(ncFileID, StateVarId, "QR_units","kg/kg"))

   ! Leave define mode so we can actually fill the variables.

   call check(nf90_enddef(ncfileID))

   call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,wrf%model_size) /) ))
   call check(nf90_put_var(ncFileID,      LonVarID, wrf%longitude(1:wrf%we,1)    ))
   call check(nf90_put_var(ncFileID,      LatVarID, wrf%latitude(1,1:wrf%sn)     ))
   call check(nf90_put_var(ncFileID,      iLevVarID, wrf%dn(1:wrf%bt)             ))

else

   !----------------------------------------------------------------------------
   ! We need to retain the prognostic variables.
   !----------------------------------------------------------------------------
   ! We could copy a lot of attributes from wrfinput to the output,
   ! call check(nf90_open('wrfinput',NF90_NOWRITE, ncid))


   !----------------------------------------------------------------------------
   ! Create the (empty) static variables and their attributes
   ! Commented block is from wrfinput
   !----------------------------------------------------------------------------
   call check(nf90_def_var(ncFileID, name="DN", xtype=nf90_double, &
              dimids= btDimID, varid=DNVarID) )
   call check(nf90_put_att(ncFileID, DNVarID, "long_name", "dn values on half (mass) levels"))
   call check(nf90_put_att(ncFileID, DNVarID, "cartesian_axis", "X"))
   call check(nf90_put_att(ncFileID, DNVarID, "units", "dimensionless"))

   call check(nf90_def_var(ncFileID, name="ZNU", xtype=nf90_double, &
              dimids= btDimID, varid=ZNUVarID) )
   call check(nf90_put_att(ncFileID, ZNUVarID, "long_name", "eta values on half (mass) levels"))
   call check(nf90_put_att(ncFileID, ZNUVarID, "cartesian_axis", "X"))
   call check(nf90_put_att(ncFileID, ZNUVarID, "units", "dimensionless"))

   call check(nf90_def_var(ncFileID, name="DNW", xtype=nf90_double, &
              dimids= btDimID, varid=DNWVarID) )
   call check(nf90_put_att(ncFileID, DNWVarID, "long_name", "dn values on full (w) levels"))
   call check(nf90_put_att(ncFileID, DNWVarID, "cartesian_axis", "X"))
   call check(nf90_put_att(ncFileID, DNWVarID, "units", "dimensionless"))

   !
   !    float MUB(Time, south_north, west_east) ;
   !            MUB:FieldType = 104 ;
   !            MUB:MemoryOrder = "XY " ;
   !            MUB:description = "base state dry air mass in column" ;
   !            MUB:units = "pascals" ;
   !            MUB:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="MUB", xtype=nf90_double, &
                 dimids= (/ weDimID, snDimID/), varid=MubVarID) )
   call check(nf90_put_att(ncFileID, MubVarID, "long_name", "base state dry air mass in column"))
   call check(nf90_put_att(ncFileID, MubVarID, "cartesian_axis", "N"))
   call check(nf90_put_att(ncFileID, MubVarID, "units", "pascals"))
   call check(nf90_put_att(ncFileID, MubVarID, "description", "base state dry air mass in column"))

   ! Longitudes
   !      float XLONG(Time, south_north, west_east) ;
   !         XLONG:FieldType = 104 ;
   !         XLONG:MemoryOrder = "XY " ;
   !         XLONG:description = "LONGITUDE, WEST IS NEGATIVE" ;
   !         XLONG:units = "degree" ;
   !         XLONG:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="XLON", xtype=nf90_double, &
                 dimids= (/ weDimID, snDimID/), varid=LonVarID) )
   call check(nf90_put_att(ncFileID, LonVarID, "long_name", "longitude"))
   call check(nf90_put_att(ncFileID, LonVarID, "cartesian_axis", "X"))
   call check(nf90_put_att(ncFileID, LonVarID, "units", "degrees_east"))
   call check(nf90_put_att(ncFileID, LonVarID, "valid_range", (/ -180.0_r8, 180.0_r8 /)))
   call check(nf90_put_att(ncFileID, LonVarID, "description", "LONGITUDE, WEST IS NEGATIVE"))

   ! Latitudes
   !      float XLAT(Time, south_north, west_east) ;
   !         XLAT:FieldType = 104 ;
   !         XLAT:MemoryOrder = "XY " ;
   !         XLAT:description = "LATITUDE, SOUTH IS NEGATIVE" ;
   !         XLAT:units = "degree" ;
   !         XLAT:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="XLAT", xtype=nf90_double, &
                 dimids=(/ weDimID, snDimID /), varid=LatVarID) ) 
   call check(nf90_put_att(ncFileID, LatVarID, "long_name", "latitude"))
   call check(nf90_put_att(ncFileID, LatVarID, "cartesian_axis", "Y"))
   call check(nf90_put_att(ncFileID, LatVarID, "units", "degrees_north"))
   call check(nf90_put_att(ncFileID, LatVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)))
   call check(nf90_put_att(ncFileID, LatVarID, "description", "LATITUDE, SOUTH IS NEGATIVE"))

   ! grid levels
   call check(nf90_def_var(ncFileID, name="level", xtype=nf90_short, &
                 dimids=btDimID, varid=ilevVarID) )
   call check(nf90_put_att(ncFileID, ilevVarID, "long_name", "placeholder for level"))
   call check(nf90_put_att(ncFileID, ilevVarID, "cartesian_axis", "Z"))
   call check(nf90_put_att(ncFileID, ilevVarID, "units", "at this point, indexical"))

   ! Land Mask
   !    float XLAND(Time, south_north, west_east) ;
   !            XLAND:FieldType = 104 ;
   !            XLAND:MemoryOrder = "XY " ;
   !            XLAND:description = "LAND MASK (1 FOR LAND, 2 FOR WATER)" ;
   !            XLAND:units = "NA" ;
   !            XLAND:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="XLAND", xtype=nf90_short, &
                 dimids= (/ weDimID, snDimID/), varid=XlandVarID) )
   call check(nf90_put_att(ncFileID, XlandVarID, "long_name", "land mask"))
   call check(nf90_put_att(ncFileID, XlandVarID, "units", "NA"))
   call check(nf90_put_att(ncFileID, XlandVarID, "valid_range", (/ 1, 2 /)))
   call check(nf90_put_att(ncFileID, XlandVarID, "description", "1 = LAND, 2 = WATER"))

   ! Map Scale Factor on m-grid
   !    float MAPFAC_M(Time, south_north, west_east) ;
   !            MAPFAC_M:FieldType = 104 ;
   !            MAPFAC_M:MemoryOrder = "XY " ;
   !            MAPFAC_M:description = "Map scale factor on mass grid" ;
   !            MAPFAC_M:units = "dimensionless" ;
   !            MAPFAC_M:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="MAPFAC_M", xtype=nf90_real, &
                 dimids= (/ weDimID, snDimID/), varid=MapFacMVarID) )
   call check(nf90_put_att(ncFileID, MapFacMVarID, "long_name", "Map scale factor on mass grid"))
   call check(nf90_put_att(ncFileID, MapFacMVarID, "units", "dimensionless"))

   ! Map Scale Factor on u-grid
   !    float MAPFAC_U(Time, south_north, west_east_stag) ;
   !            MAPFAC_U:FieldType = 104 ;
   !            MAPFAC_U:MemoryOrder = "XY " ;
   !            MAPFAC_U:description = "Map scale factor on u-grid" ;
   !            MAPFAC_U:units = "dimensionless" ;
   !            MAPFAC_U:stagger = "X" ;
   call check(nf90_def_var(ncFileID, name="MAPFAC_U", xtype=nf90_real, &
                 dimids= (/ weStagDimID, snDimID/), varid=MapFacUVarID) )
   call check(nf90_put_att(ncFileID, MapFacUVarID, "long_name", "Map scale factor on u-grid"))
   call check(nf90_put_att(ncFileID, MapFacUVarID, "units", "dimensionless"))

   ! Map Scale Factor on v-grid
   !    float MAPFAC_V(Time, south_north_stag, west_east) ;
   !            MAPFAC_V:FieldType = 104 ;
   !            MAPFAC_V:MemoryOrder = "XY " ;
   !            MAPFAC_V:description = "Map scale factor on v-grid" ;
   !            MAPFAC_V:units = "dimensionless" ;
   !            MAPFAC_V:stagger = "Y" ;
   call check(nf90_def_var(ncFileID, name="MAPFAC_V", xtype=nf90_real, &
                 dimids= (/ weDimID, snStagDimID/), varid=MapFacVVarID) )
   call check(nf90_put_att(ncFileID, MapFacVVarID, "long_name", "Map scale factor on v-grid"))
   call check(nf90_put_att(ncFileID, MapFacVVarID, "units", "dimensionless"))

   ! PHB
   !    float PHB(Time, bottom_top_stag, south_north, west_east) ;
   !            PHB:FieldType = 104 ;
   !            PHB:MemoryOrder = "XYZ" ;
   !            PHB:description = "base-state geopotential" ;
   !            PHB:units = "m{2} s{-2}" ;
   !            PHB:stagger = "Z" ;
   call check(nf90_def_var(ncFileID, name="PHB", xtype=nf90_real, &
                 dimids= (/ weDimID, snDimID, btStagDimID /), varid=phbVarId) )
   call check(nf90_put_att(ncFileID, phbVarId, "long_name", "base-state geopotential"))
   call check(nf90_put_att(ncFileID, phbVarId, "units", "m^2/s^2"))
   call check(nf90_put_att(ncFileID, phbVarId, "units_long_name", "m{2} s{-2}"))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and their attributes
   !----------------------------------------------------------------------------

   !      float U(Time, bottom_top, south_north, west_east_stag) ;
   !         U:FieldType = 104 ;
   !         U:MemoryOrder = "XYZ" ;
   !         U:description = "x-wind component" ;
   !         U:units = "m s{-1}" ;
   !         U:stagger = "X" ;
   call check(nf90_def_var(ncid=ncFileID, name="U", xtype=nf90_real, &
         dimids = (/ weStagDimID, snDimId, btDimID, MemberDimID, unlimitedDimID /), &
         varid  = UVarID))
   call check(nf90_put_att(ncFileID, UVarID, "long_name", "x-wind component"))
   call check(nf90_put_att(ncFileID, UVarID, "units", "m/s"))
   call check(nf90_put_att(ncFileID, UVarID, "units_long_name", "m s{-1}"))


   !      float V(Time, bottom_top, south_north_stag, west_east) ;
   !         V:FieldType = 104 ;
   !         V:MemoryOrder = "XYZ" ;
   !         V:description = "y-wind component" ;
   !         V:units = "m s{-1}" ;
   !         V:stagger = "Y" ;
   call check(nf90_def_var(ncid=ncFileID, name="V", xtype=nf90_real, &
         dimids = (/ weDimID, snStagDimID, btDimID, MemberDimID, unlimitedDimID /), &
         varid  = VVarID))
   call check(nf90_put_att(ncFileID, VVarID, "long_name", "y-wind component"))
   call check(nf90_put_att(ncFileID, VVarID, "units", "m/s"))
   call check(nf90_put_att(ncFileID, VVarID, "units_long_name", "m s{-1}"))


   !      float W(Time, bottom_top_stag, south_north, west_east) ;
   !         W:FieldType = 104 ;
   !         W:MemoryOrder = "XYZ" ;
   !         W:description = "z-wind component" ;
   !         W:units = "m s{-1}" ;
   !         W:stagger = "Z" ;
   call check(nf90_def_var(ncid=ncFileID, name="W", xtype=nf90_real, &
         dimids = (/ weDimID, snDimID, btStagDimID, MemberDimID, unlimitedDimID /), &
         varid  = WVarID))
   call check(nf90_put_att(ncFileID, WVarID, "long_name", "z-wind component"))
   call check(nf90_put_att(ncFileID, WVarID, "units", "m/s"))
   call check(nf90_put_att(ncFileID, WVarID, "units_long_name", "m s{-1}"))


   !      float PH(Time, bottom_top_stag, south_north, west_east) ;               
   !         PH:FieldType = 104 ;
   !         PH:MemoryOrder = "XYZ" ;
   !         PH:description = "perturbation geopotential" ;
   !         PH:units = "m{2} s{-2}" ;
   !         PH:stagger = "Z" ;
   call check(nf90_def_var(ncid=ncFileID, name="PH", xtype=nf90_real, &
         dimids = (/ weDimID, snDimID, btStagDimID, MemberDimID, unlimitedDimID /), &
         varid  = PHVarID))
   call check(nf90_put_att(ncFileID, PHVarID, "long_name", "perturbation geopotential"))
   call check(nf90_put_att(ncFileID, PHVarID, "units", "m^2/s^2"))
   call check(nf90_put_att(ncFileID, PHVarID, "units_long_name", "m{2} s{-2}"))


   !      float T(Time, bottom_top, south_north, west_east) ;
   !         T:FieldType = 104 ;
   !         T:MemoryOrder = "XYZ" ;
   !         T:description = "perturbation potential temperature (theta-t0)" ;
   !         T:units = "K" ;
   !         T:stagger = "" ;
   call check(nf90_def_var(ncid=ncFileID, name="T", xtype=nf90_real, &
         dimids = (/ weDimID, snDimID, btDimID, MemberDimID, unlimitedDimID /), &
         varid  = tVarID))
   call check(nf90_put_att(ncFileID, tVarID, "long_name", "temperature"))
   call check(nf90_put_att(ncFileID, tVarID, "units", "degrees Kelvin"))
   call check(nf90_put_att(ncFileId, tVarID, "description", "perturbation potential temperature (theta-t0)"))


   !      float MU(Time, south_north, west_east) ;
   !         MU:FieldType = 104 ;
   !         MU:MemoryOrder = "XY " ;
   !         MU:description = "perturbation dry air mass in column" ;
   !         MU:units = "pascals" ;
   !         MU:stagger = "" ;
   call check(nf90_def_var(ncid=ncFileID, name="MU", xtype=nf90_real, &
         dimids = (/ weDimID, snDimID, MemberDimID, unlimitedDimID /), &
         varid  = MuVarID))
   call check(nf90_put_att(ncFileID, MuVarID, "long_name", "mu field"))
   call check(nf90_put_att(ncFileID, MuVarID, "units", "pascals"))
   call check(nf90_put_att(ncFileId, MuVarID, "description", "perturbation dry air mass in column"))


   !      float QVAPOR(Time, bottom_top, south_north, west_east) ;
   !         QVAPOR:FieldType = 104 ;
   !         QVAPOR:MemoryOrder = "XYZ" ;
   !         QVAPOR:description = "-" ;
   !         QVAPOR:units = "-" ;
   !         QVAPOR:stagger = "" ;
   call check(nf90_def_var(ncid=ncFileID, name="QVAPOR", xtype=nf90_real, &
         dimids = (/ weDimID, snDimID, btDimID, MemberDimID, unlimitedDimID /), &
         varid  = QVVarID))
   call check(nf90_put_att(ncFileID, QVVarID, "long_name", "-"))
   call check(nf90_put_att(ncFileID, QVVarID, "units", "kg/kg"))
   call check(nf90_put_att(ncFileId, QVVarID, "description", "-"))


   !      float QCLOUD(Time, bottom_top, south_north, west_east) ;
   !         QCLOUD:FieldType = 104 ;
   !         QCLOUD:MemoryOrder = "XYZ" ;
   !         QCLOUD:description = "-" ;
   !         QCLOUD:units = "-" ;
   !         QCLOUD:stagger = "" ;
   call check(nf90_def_var(ncid=ncFileID, name="QCLOUD", xtype=nf90_real, &
         dimids = (/ weDimID, snDimID, btDimID, MemberDimID, unlimitedDimID /), &
         varid  = QCVarID))
   call check(nf90_put_att(ncFileID, QCVarID, "long_name", "-"))
   call check(nf90_put_att(ncFileID, QCVarID, "units", "kg/kg"))
   call check(nf90_put_att(ncFileId, QCVarID, "description", "-"))


   !      float QRAIN(Time, bottom_top, south_north, west_east) ;
   !         QRAIN:FieldType = 104 ;
   !         QRAIN:MemoryOrder = "XYZ" ;
   !         QRAIN:description = "-" ;
   !         QRAIN:units = "-" ;
   !         QRAIN:stagger = "" ;
   call check(nf90_def_var(ncid=ncFileID, name="QRAIN", xtype=nf90_real, &
         dimids = (/ weDimID, snDimID, btDimID, MemberDimID, unlimitedDimID /), &
         varid  = QRVarID))
   call check(nf90_put_att(ncFileID, QRVarID, "long_name", "-"))
   call check(nf90_put_att(ncFileID, QRVarID, "units", "kg/kg"))
   call check(nf90_put_att(ncFileId, QRVarID, "description", "-"))

   if ( wrf%n_moist > 3 ) then
      write(*,*)' YO DUMMY -- NEED TO INITIALIZE THE SOLID PHASE WATER VARS'
      write(*,*)' YO DUMMY -- NEED TO INITIALIZE THE SOLID PHASE WATER VARS'
      write(*,*)' YO DUMMY -- NEED TO INITIALIZE THE SOLID PHASE WATER VARS'
      write(*,*)' YO DUMMY -- NEED TO INITIALIZE THE SOLID PHASE WATER VARS'
      stop
   endif

   !-----------------------------------------------------------------
   ! Fill the variables we can
   !-----------------------------------------------------------------
   call check(nf90_enddef(ncfileID))

   ! defining grid levels
   call check(nf90_put_var(ncFileID,       DNVarID, wrf%dn        ))
   call check(nf90_put_var(ncFileID,      ZNUVarID, wrf%znu       ))
   call check(nf90_put_var(ncFileID,      DNWVarID, wrf%dnw       ))
   ! defining horizontal
   call check(nf90_put_var(ncFileID,      mubVarID, wrf%mub       )) 
   call check(nf90_put_var(ncFileID,      LonVarID, wrf%longitude ))
   call check(nf90_put_var(ncFileID,      LatVarID, wrf%latitude  )) 
   call check(nf90_put_var(ncFileID,     ilevVarID, (/ (i,i=1,wrf%bt) /) )) 
   call check(nf90_put_var(ncFileID,    XlandVarID, wrf%land      )) 
   call check(nf90_put_var(ncFileID,  MapFacMVarID, wrf%mapfac_m  )) 
   call check(nf90_put_var(ncFileID,  MapFacUVarID, wrf%mapfac_u  )) 
   call check(nf90_put_var(ncFileID,  MapFacVVarID, wrf%mapfac_v  )) 
   call check(nf90_put_var(ncFileID,      phbVarID, wrf%phb       )) 

endif

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

logical, parameter :: debug = .true.  
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, VarID
integer :: i,j
real(r8), allocatable, dimension(:,:)   :: temp2d
real(r8), allocatable, dimension(:,:,:) :: temp3d
character(len=8) :: varname

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

   !----------------------------------------------------------------------------
   ! Fill the variables, the order is CRITICAL  ...   U,V,W,GZ,T,MU,QV,QC,QR
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   varname = 'U'
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = 1
   j       = i + wrf%wes * wrf%sn * wrf%bt - 1 
   if (debug) write(*,'(a7,'' = statevec('',i6,'':'',i6,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%wes,wrf%sn,wrf%bt 
   allocate ( temp3d(wrf%wes, wrf%sn, wrf%bt) )
   temp3d  = reshape(statevec(i:j), (/ wrf%wes, wrf%sn, wrf%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'V'
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%we * wrf%sns * wrf%bt - 1
   if (debug) write(*,'(a7,'' = statevec('',i6,'':'',i6,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%we,wrf%sns,wrf%bt
   allocate ( temp3d(wrf%we, wrf%sns, wrf%bt) )
   temp3d  = reshape(statevec(i:j), (/ wrf%we, wrf%sns, wrf%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'W'
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%we * wrf%sn * wrf%bts - 1
   if (debug) write(*,'(a7,'' = statevec('',i6,'':'',i6,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%we,wrf%sn,wrf%bts
   allocate ( temp3d(wrf%we, wrf%sn, wrf%bts) )
   temp3d  = reshape(statevec(i:j), (/ wrf%we, wrf%sn, wrf%bts /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'PH'       ! AKA "GZ"
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%we * wrf%sn * wrf%bts - 1
   if (debug) write(*,'(a7,'' = statevec('',i6,'':'',i6,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%we,wrf%sn,wrf%bts
   allocate ( temp3d(wrf%we, wrf%sn, wrf%bts) )
   temp3d  = reshape(statevec(i:j), (/ wrf%we, wrf%sn, wrf%bts /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'T  '
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%we * wrf%sn * wrf%bt - 1
   if (debug) write(*,'(a7,'' = statevec('',i6,'':'',i6,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%we,wrf%sn,wrf%bt
   allocate ( temp3d(wrf%we, wrf%sn, wrf%bt) )
   temp3d  = reshape(statevec(i:j), (/ wrf%we, wrf%sn, wrf%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   ! deallocate(temp3d)  remaining 3D variables same size ...


   !----------------------------------------------------------------------------
   varname = 'MU '
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%we * wrf%sn - 1
   if (debug) write(*,'(a7,'' = statevec('',i6,'':'',i6,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%we,wrf%sn
   allocate ( temp2d(wrf%we, wrf%sn) )
   temp2d  = reshape(statevec(i:j), (/ wrf%we, wrf%sn /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp2d, &
                            start=(/ 1, 1, copyindex, timeindex /) ))
   deallocate(temp2d)


   !----------------------------------------------------------------------------
   varname = 'QVAPOR '
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%we * wrf%sn * wrf%bt - 1
   if (debug) write(*,'(a7,'' = statevec('',i6,'':'',i6,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%we,wrf%sn,wrf%bt
   temp3d  = reshape(statevec(i:j), (/ wrf%we, wrf%sn, wrf%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))


   !----------------------------------------------------------------------------
   varname = 'QCLOUD '
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%we * wrf%sn * wrf%bt - 1
   if (debug) write(*,'(a7,'' = statevec('',i6,'':'',i6,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%we,wrf%sn,wrf%bt
   temp3d  = reshape(statevec(i:j), (/ wrf%we, wrf%sn, wrf%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))


   !----------------------------------------------------------------------------
   varname = 'QRAIN '
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%we * wrf%sn * wrf%bt - 1
   if (debug) write(*,'(a7,'' = statevec('',i6,'':'',i6,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%we,wrf%sn,wrf%bt
   temp3d  = reshape(statevec(i:j), (/ wrf%we, wrf%sn, wrf%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))

   deallocate(temp3d)


   if ( wrf%n_moist > 3 ) then
      write(*,*)' YO DUMMY -- NEED TO OUTPUT THE SOLID PHASE WATER VARS'
      write(*,*)' YO DUMMY -- NEED TO OUTPUT THE SOLID PHASE WATER VARS'
      write(*,*)' YO DUMMY -- NEED TO OUTPUT THE SOLID PHASE WATER VARS'
      write(*,*)' YO DUMMY -- NEED TO OUTPUT THE SOLID PHASE WATER VARS'
      stop
   endif

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



subroutine llxy (xloni,xlatj,x,y)
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

   real(r8) :: dxlon
   real(r8) :: xlat, xlon
   real(r8) :: xx, yy, xc, yc
   real(r8) :: cell, psi0, psx, r, flp
   real(r8) :: centri, centrj
   real(r8) :: ds       
   real(r8) :: bb,c2
   
!-----------------------------------------------------------------
   ds = 0.001 *wrf%dx
   xlon = xloni
   xlat = xlatj
   xlat = max (xlat, -89.9999)
   xlat = min (xlat, +89.9999)
   
!-----------------------------------------------------------------
   c2 = earth_radius * COS(wrf%psi1)

   if (wrf%map_proj == 3) then
      xc = 0.0
      yc = wrf%ycntr 

      cell = cos(xlat*deg2rad)/(1.0+sin(xlat*deg2rad))
      yy = -c2*log(cell)
      xx = c2*(xlon-wrf%cen_lon)*deg2rad

   else

      psi0 = ( 90.0 - wrf%cen_lat)*deg2rad
      xc = 0.0

!-----CALCULATE X,Y COORDS. RELATIVE TO POLE

      dxlon = xlon - wrf%cen_lon
      if (dxlon >  180) dxlon = dxlon - 360.
      if (dxlon < -180) dxlon = dxlon + 360.
   
      flp = wrf%cone_factor*dxlon*deg2rad
   
      psx = ( 90.0 - xlat )*deg2rad
   
      if (wrf%map_proj == 2) then
! ...... Polar stereographics:
         bb = 2.0*(cos(wrf%psi1/2.0)**2)
         yc = -earth_radius*bb*tan(psi0/2.0)
          r = -earth_radius*bb*tan(psx/2.0)
      else
! ...... Lambert conformal:
         bb = -earth_radius/wrf%cone_factor*sin(wrf%psi1)
         yc = bb*(tan(psi0/2.0)/tan(wrf%psi1/2.0))**wrf%cone_factor
          r = bb*(tan(psx /2.0)/tan(wrf%psi1/2.0))**wrf%cone_factor
      endif

      if (wrf%cen_lat < 0.0) then
         xx = r*sin(flp)
         yy = r*cos(flp)
      else
         xx = -r*sin(flp)
         yy =  r*cos(flp)
      endif

   endif

! TRANSFORM (1,1) TO THE ORIGIN
! the location of the center in the coarse domain

   centri = real (wrf%we)/2.0  
   centrj = real (wrf%sn)/2.0  
! the (X,Y) coordinates in the coarse domain
   x = ( xx - xc )/ds + centri  
   y = ( yy - yc )/ds + centrj  
             

!--only add 0.5 so that x/y is relative to first cross points (MM5 input):

   x = (x - 1.0) + 0.5
   y = (y - 1.0) + 0.5

end subroutine llxy
   
!**********************************************
SUBROUTINE XYLL(XX,YY,XLAT,XLON)
!               
!   PURPOSE      : CALCULATES THE LATITUDES AND LONGITUDES FROM THE
!                  (X,Y) LOCATION (DOT) IN THE MESOSCALE GRIDS.
!   ON ENTRY     :   
!   X            : THE COORDINATE IN X (J)-DIRECTION.
!   Y            : THE COORDINATE IN Y (I)-DIRECTION.
!
!   ON EXIT      :                      
!   XLAT         : LATITUDES 
!   XLON         : LONGITUDES 
!

   REAL(R8), INTENT(IN)  :: XX, YY
   REAL(R8), INTENT(OUT) :: XLAT,XLON
        
   REAL(R8) :: flp, flpp, r, cell, cel1, cel2, c2
   REAL(R8) :: psx,Rcone_factor
   REAL(R8) :: centri, centrj, x, y, xcntr

   c2 = earth_radius * COS(wrf%psi1)
   centri = wrf%we / 2.0_r8
   centrj = wrf%sn / 2.0_r8
!   CNTRI = float(coarse_iy+1)/2.
!   CNTRJ = float(coarse_jx+1)/2. 
   
   xcntr = 0.0_r8

!-----CALCULATE X AND Y POSITIONS OF GRID
   X = ( XCNTR+(XX-1.0) +(1.0 - CENTRI) ) * wrf%dx*0.001 
   Y = ( wrf%YCNTR+(YY-1.0) +(1.0 - CENTRJ) ) * wrf%dx*0.001
!-----NOW CALCULATE LAT AND LON OF THIS POINT

   IF (wrf%map_proj.NE.3) THEN
      IF(Y.EQ.0.) THEN      
        IF(X.GE.0.0) FLP =  90.0*deg2rad 
        IF(X.LT.0.0) FLP = -90.0*deg2rad
      ELSE
        IF (wrf%cen_lat.LT.0.0)THEN
            FLP = ATAN2(X,Y)   
        ELSE
            FLP = ATAN2(X,-Y) 
        ENDIF
      ENDIF 
      FLPP = (FLP/wrf%cone_factor)*rad2deg+wrf%cen_lon
      IF (FLPP.LT.-180.) FLPP = FLPP + 360    
      IF (FLPP.GT.180.)  FLPP = FLPP - 360.  
      XLON = FLPP 
!--------NOW SOLVE FOR LATITUDE
      R = SQRT(X*X+Y*Y)  
      IF (wrf%cen_lat.LT.0.0) R = -R  
      IF (wrf%map_proj.EQ.1) THEN   
         CELL = (R*wrf%cone_factor)/(earth_radius*SIN(wrf%PSI1))    
         Rcone_factor  = 1.0/wrf%cone_factor   
         CEL1 = TAN(wrf%PSI1/2.)*(CELL)**Rcone_factor    
      ENDIF 
      IF (wrf%map_proj.EQ.2) THEN
         CELL = R/earth_radius        
         CEL1 = CELL/(1.0+COS(wrf%PSI1))  
      ENDIF 
      CEL2 = ATAN(CEL1)    
      PSX  = 2.*CEL2*rad2deg 
      XLAT = 90.0-PSX 
   ENDIF   
!-----CALCULATIONS FOR MERCATOR LAT,LON    
   IF (wrf%map_proj.EQ.3) THEN   
      XLON = wrf%cen_lon + ((X-XCNTR)/C2)*rad2deg
      IF (XLON.LT.-180.) XLON = XLON + 360
      IF (XLON.GT.180.)  XLON = XLON - 360.
      CELL = EXP(Y/C2)  
      XLAT = 2.*(rad2deg*ATAN(CELL))-90.0 
   ENDIF
end subroutine xyll
   
!**********************************************
subroutine Interp_lin_1D(fi1d, n1, z, fo1d)                        

  integer,  intent(in)  :: n1
  real(r8), intent(in)  :: fi1d(n1)
  real(r8), intent(in)  :: z
  real(r8), intent(out) :: fo1d

  integer   :: k
  real(r8)  :: dz, dzm

  fo1d = missing_r

     if(z > 0.0_r8) then
        call toGrid(z,n1, k, dz, dzm)
        fo1d = dzm*fi1d(k) + dz*fi1d(k+1)
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

  call toGrid (x,n1,i,dx,dxm)
  call toGrid (y,n2,j,dy,dym)

  fo2d   = dym*(dxm*fi2d(i,j  ) + dx*fi2d(i+1,j  )) &
         + dy *(dxm*fi2d(i,j+1) + dx*fi2d(i+1,j+1))
end subroutine Interp_lin_2D
   
!**********************************************

subroutine Interp_lin_3D(fi3d,n1,n2,n3, x,y,z,fo3d)

! real(r8),dimension(n1,n2,n3), intent(in)  :: fi3d(n1,n2,n3)
  integer,  intent(in)  :: n1,n2,n3
  real(r8), intent(in)  :: fi3d(n1,n2,n3)
  real(r8), intent(in)  :: x, y, z
  real(r8), intent(out) :: fo3d

  integer   :: i, j, k, kk
  real(r8)  :: dx, dxm, dy, dym, dz, dzm
  real(r8)  :: fiz (n3)

  call toGrid (x,n1,i,dx,dxm)
  call toGrid (y,n2,j,dy,dym)

  fiz(1:n3) = dym*(dxm*fi3d(i, j,   1:n3) + dx *fi3d(i+1 ,j  ,1:n3))&
            + dy *(dxm*fi3d(i, j+1, 1:n3) + dx *fi3d(i+1, j+1,1:n3))

  fo3d = missing_r

     if(z > 0.0) then
        call toGrid(z,n3, k, dz, dzm)
        fo3d = dzm*fiz(k) + dz*fiz(k+1)
     endif

end subroutine Interp_lin_3D

!#######################################################################
subroutine toGrid (x, jx, j, dx, dxm)
   
!  Transfer obs. x to grid j and calculate its
!  distance to grid j and j+1

   real(r8), intent(in)  :: x
   integer,  intent(in)  :: jx
   real(r8), intent(out) :: dx, dxm
   integer,  intent(out) :: j
   
   j = int (x)

   if (j <=  0) j = 1
   if (j >= jx) j = jx - 1

   dx = x - real (j)

   dxm= 1.0 - dx

end subroutine toGrid
!#######################################################################

subroutine to_zk(obs_v, mdl_v, n3, v_interp_optn, zk)

   real(r8), intent(in)  :: obs_v
   integer,  intent(in)  :: n3, v_interp_optn
   real(r8), intent(in)  :: mdl_v(n3)
   real(r8), intent(out) :: zk

   integer   :: k

   zk = missing_r

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
subroutine get_model_pressure_profile(i,j,dx,dy,dxm,dym,n,x,fld,&
                                      pp11,pp21,pp31,pp41)

integer,  intent(in)  :: i,j,n
real(r8), intent(in)  :: dx,dy,dxm,dym
real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: fld(n)      
real(r8), intent(out) :: pp11,pp21,pp31,pp41
 
integer               :: i1,i2,k,q1,q2,q3,q4
real(r8)              :: qv1,qv2,qv3,qv4
real(r8), dimension(n):: pp1,pp2,pp3,pp4,pb,pp


do k=1,n
   pb(k) = wrf%p_top + wrf%znu(k)*  &
          (  dym*(dxm*wrf%mub(i,j  ) + dx*wrf%mub(i+1,j  )) + &
             dy *(dxm*wrf%mub(i,j+1) + dx*wrf%mub(i+1,j+1)) )  
end do
  i1 = get_wrf_index(i,j,1,TYPE_MU)
  i2 = get_wrf_index(i,j+1,1,TYPE_MU)
  q1 = get_wrf_index(i,j,n,TYPE_QV)
  q2 = get_wrf_index(i,j+1,n,TYPE_QV)
  qv1 = x(q1)/(1.0+x(q1))
  qv2 = x(q1+1)/(1.0+x(q1+1))
  qv3 = x(q2)/(1.0+x(q2))
  qv4 = x(q2+1)/(1.0+x(q2+1))
  pp1(n) = -0.5 *(x(i1)  +qv1*wrf%mub(i  ,j  ))*wrf%dnw(n)*(1.0 + x(q1))
  pp2(n) = -0.5 *(x(i1+1)+qv2*wrf%mub(i+1,j  ))*wrf%dnw(n)*(1.0 + x(q1+1))
  pp3(n) = -0.5 *(x(i2)  +qv3*wrf%mub(i  ,j+1))*wrf%dnw(n)*(1.0 + x(q2))
  pp4(n) = -0.5 *(x(i2+1)+qv4*wrf%mub(i+1,j+1))*wrf%dnw(n)*(1.0 + x(q2+1))
  pp(n)  = dym*(dxm*pp1(n)+dx*pp2(n)) + dy*(dxm*pp3(n)+dx*pp4(n))
  fld(n) = pp(n) + pb(n)                                            
do k= n-1,1,-1   
   q1 = get_wrf_index(i,j,k,TYPE_QV)
   q2 = get_wrf_index(i,j+1,k,TYPE_QV)
   q3 = get_wrf_index(i,j,k+1,TYPE_QV)
   q4 = get_wrf_index(i,j+1,k+1,TYPE_QV)
   qv1 = 0.5*(x(q1)+x(q3))/(1.0+0.5*(x(q1)+x(q3)))
   qv2 = 0.5*(x(q1+1)+x(q3+1))/(1.0+0.5*(x(q1+1)+x(q3+1)))
   qv3 = 0.5*(x(q2)+x(q4))/(1.0+0.5*(x(q2)+x(q4)))
   qv4 = 0.5*(x(q2+1)+x(q4+1))/(1.0+0.5*(x(q2+1)+x(q4+1)))

   pp1(k) = pp1(k+1) -(x(i1)  +qv1*wrf%mub(i  ,j  ))*wrf%dn(k+1)* &
                (1.0 + 0.5*(x(q1) + x(q3)))
   pp2(k) = pp2(k+1) -(x(i1+1)+qv2*wrf%mub(i+1,j  ))*wrf%dn(k+1)* &
                (1.0 + 0.5*(x(q1+1) + x(q3+1)))
   pp3(k) = pp3(k+1) -(x(i2)  +qv3*wrf%mub(i  ,j+1))*wrf%dn(k+1)* &
                (1.0 + 0.5*(x(q2) + x(q4)))
   pp4(k) = pp4(k+1) -(x(i2+1)  +qv4*wrf%mub(i+1,j+1))*wrf%dn(k+1)* &
                (1.0 + 0.5*(x(q2+1) + x(q4+1)))

   pp(k)  = dym*(dxm*pp1(k)+dx*pp2(k)) + dy*(dxm*pp3(k)+dx*pp4(k))
   fld(k) = pp(k) + pb(k)                                             
end do
   pp11 = pp1(1)
   pp21 = pp2(1)
   pp31 = pp3(1)
   pp41 = pp4(1)
end subroutine get_model_pressure_profile
!#######################################################
subroutine get_model_height_profile(i,j,dx,dy,dxm,dym,n,x,fld)

integer,  intent(in)  :: i,j,n
real(r8), intent(in)  :: dx,dy,dxm,dym
real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: fld(n)

real(r8)  :: fll(n+1) 
integer   :: i1,i2,k     
            
do k = 1, wrf%var_size(3,TYPE_GZ) 
   i1 = get_wrf_index(i,j,k,TYPE_GZ)
   i2 = get_wrf_index(i,j+1,k,TYPE_GZ)
   fll(k) = (dym*( dxm*(wrf%phb(i,j,k)+x(i1))+dx*(wrf%phb(i+1,j,k)+x(i1+1)))+&
   dy*(dxm*(wrf%phb(i,j+1,k)+x(i2)) + dx*(wrf%phb(i+1,j+1,k)+x(i2+1)) ))/gravity
end do

do k=1,n        
   fld(k) = 0.5*(fll(k) + fll(k+1) )
end do
end subroutine get_model_height_profile



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



end module model_mod

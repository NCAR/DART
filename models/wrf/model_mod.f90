module model_mod

! Assimilation interface for WRF model

!-----------------------------------------------------------------------
!
!     interface for WRF
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------

!-----------------------------------------------------------------------
use time_manager_mod, only : time_type, set_time
use location_mod    , only : location_type, get_location, set_location, get_dist, &
                             LocationDims, LocationName, LocationLName
use types_mod

implicit none
private

!  public routines and data for the WRF model

public     get_model_size,                    &
           get_state_meta_data,               &
           model_interpolate,                 &
           get_model_time_step,               &
           static_init_model,                 &
           model_get_close_states,            &
           TYPE_MU, TYPE_T, TYPE_U, TYPE_V,    &
           TYPE_W, TYPE_GZ,                   &
           TYPE_QV, TYPE_QC, TYPE_QR

!  public stubs 

public     adv_1step,           &
           end_model,           &
           init_time,           &
           init_conditions,     &
           nc_write_model_atts , &
        nc_write_model_vars

!-----------------------------------------------------------------------

character(len=128) :: version = "$Id$"

character(len=128) :: tag = "$Name$"

character(len=128) :: &                                                                              
   source = "$Source$", &                    
   revision = "$Revision$", &                                                                 
   revdate  = "$Date$"

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

      logical :: output_state_vector = .true.  ! output prognostic variables

! Public definition of variable types

integer, parameter :: TYPE_MU = 0, TYPE_T = 1, TYPE_U = 2, TYPE_V = 3, &
                      TYPE_W = 4,  TYPE_GZ = 5,                        &
                      TYPE_QV = 6,  TYPE_QC = 7,  TYPE_QR = 8,         &
                      TYPE_QI = 9,  TYPE_QS = 10, TYPE_QG = 11


!-----------------------------------------------------------------------
!---- private data ----

TYPE wrf_static_data_for_dart

   integer :: bt, sn, we

   real :: p_top, dx, dy, dt
   real, dimension(:), pointer :: dn, dnw
   integer :: n_moist
   real, dimension(:,:), pointer :: mub, latitude, longitude
   real, dimension(:,:), pointer :: mapfac_m, mapfac_u, mapfac_v
   real, dimension(:,:,:), pointer :: phb

   integer :: model_size, number_of_wrf_variables
   integer, dimension(:), pointer :: var_type
   integer, dimension(:,:), pointer :: var_index
   integer, dimension(:,:), pointer :: var_size

end type

type(wrf_static_data_for_dart) :: wrf

! ***** define other auxillary model data here as needed

!------------------------------------------------------------

contains

!#######################################################################

subroutine static_init_model()

! INitializes class data for WRF???

use netcdf

implicit none

character (len = 80)      :: path
     integer              :: mode
     integer              :: ncid, bt_id, we_id, sn_id

integer :: status
character (len=80) :: name
logical, parameter :: debug = .false.
integer :: var_id, ind, i
integer, dimension(5) :: count, start, stride, map
real    :: zero_d(1)

real :: dx, dy, dt

real, dimension(:,:,:), pointer :: mub_test

integer :: n_values

!----------

mode = 0
status = nf90_open('wrfinput', mode, ncid) 
if (status /= nf90_noerr) call handle_err(1,status)
if(debug) write(6,*) ' ncid is ',ncid

! get wrf grid dimensions

status = nf90_inq_dimid(ncid, "bottom_top", bt_id)
if (status /= nf90_noerr) call handle_err(2,status)
status = nf90_inquire_dimension(ncid, bt_id, name, wrf%bt)
if (status /= nf90_noerr) call handle_err(3,status)

status = nf90_inq_dimid(ncid, "south_north", sn_id)
if (status /= nf90_noerr) call handle_err(2,status)
status = nf90_inquire_dimension(ncid, sn_id, name, wrf%sn)
if (status /= nf90_noerr) call handle_err(3,status)

status = nf90_inq_dimid(ncid, "west_east", we_id)
if (status /= nf90_noerr) call handle_err(2,status)
status = nf90_inquire_dimension(ncid, we_id, name, wrf%we)
if (status /= nf90_noerr) call handle_err(3,status)

if(debug) write(6,*) ' dimensions bt, sn, we are ',wrf%bt,wrf%sn,wrf%we

! get meta data adn static data we need

count  = 1
start  = 1
stride = 1
map    = 1

status = nf90_get_att(ncid, nf90_global, 'DX', dx)
if (status /= nf90_noerr) call handle_err(4,status)
wrf%dx = dx
if(debug) write(6,*) ' dx is ',dx

status = nf90_get_att(ncid, nf90_global, 'DY', dy)
if (status /= nf90_noerr) call handle_err(4,status)
wrf%dy = dy
if(debug) write(6,*) ' dy is ',dy

status = nf90_get_att(ncid, nf90_global, 'DT', dt)
if (status /= nf90_noerr) call handle_err(4,status)
wrf%dt = dt
if(debug) write(6,*) ' dt is ',dt

call netcdf_read_write_var( "P_TOP", ncid, var_id, zero_d,        &
                            start, count, stride, map, 'INPUT ', debug, 1  )
wrf%p_top = zero_d(1)
if(debug) write(6,*) ' p_top is ',wrf%p_top


wrf%n_moist = 3  ! hardwired for now, need a way to find this

!  get 1D (z) static data defining grid levels

count(1)  = wrf%bt
allocate(wrf%dn(1:wrf%bt))
call netcdf_read_write_var( "DN",ncid, var_id, wrf%dn,        &
                            start, count, stride, map, 'INPUT ', debug, 2  )
if(debug) write(6,*) ' dn ',wrf%dn

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

if(debug) then
    write(6,*) ' corners of lat '
    write(6,*) wrf%latitude(1,1),wrf%latitude(wrf%we,1),  &
               wrf%latitude(1,wrf%sn),wrf%latitude(wrf%we,wrf%sn)
end if

if(debug) then
    write(6,*) ' corners of long '
    write(6,*) wrf%longitude(1,1),wrf%longitude(wrf%we,1),  &
               wrf%longitude(1,wrf%sn),wrf%longitude(wrf%we,wrf%sn)
end if

allocate(wrf%mapfac_m(1:wrf%we,1:wrf%sn))
call netcdf_read_write_var( "MAPFAC_M",ncid, var_id, wrf%mapfac_m,        &
                            start, count, stride, map, 'INPUT ', debug, 3  )

count(1)  = wrf%we+1
count(2)  = wrf%sn
allocate(wrf%mapfac_u(1:wrf%we+1,1:wrf%sn))
call netcdf_read_write_var( "MAPFAC_U",ncid, var_id, wrf%mapfac_u,        &
                            start, count, stride, map, 'INPUT ', debug, 4  )

count(1)  = wrf%we
count(2)  = wrf%sn+1
allocate(wrf%mapfac_v(1:wrf%we,1:wrf%sn+1))
call netcdf_read_write_var( "MAPFAC_V",ncid, var_id, wrf%mapfac_V,        &
                            start, count, stride, map, 'INPUT ', debug, 4  )

! get 3D base state geopotential

count(1)  = wrf%we
count(2)  = wrf%sn
count(3)  = wrf%bt+1
allocate(wrf%phb(1:wrf%we,1:wrf%sn,1:wrf%bt+1))
call netcdf_read_write_var( "PHB",ncid, var_id, wrf%phb,        &
                            start, count, stride, map, 'INPUT ', debug, 4  )
if(debug) then
    write(6,*) ' corners of phb '
    write(6,*) wrf%phb(1,1,1),wrf%phb(wrf%we,1,1),  &
               wrf%phb(1,wrf%sn,1),wrf%phb(wrf%we,wrf%sn,1)
    write(6,*) wrf%phb(1,1,wrf%bt+1),wrf%phb(wrf%we,1,wrf%bt+1),  &
               wrf%phb(1,wrf%sn,wrf%bt+1),wrf%phb(wrf%we,wrf%sn,wrf%bt+1)
end if

! close data file, we have all we need

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(4,status)

!  build the map into the 1D DART vector for WRF data

  wrf%number_of_wrf_variables = 6 + wrf%n_moist
  allocate(wrf%var_type(12))  ! map for var_type, we'll always carry all
  wrf%var_type(1)  = TYPE_U   !  possible types
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

  allocate(wrf%var_index(2,6 + wrf%n_moist)) ! indices into 1D array
  allocate(wrf%var_size(3,6 + wrf%n_moist)) ! dimension of variables
  
  ind = 1                         ! *** u field ***
  wrf%var_size(1,ind) = wrf%we + 1  
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = wrf%bt
  wrf%var_index(1,ind) = 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +  &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** v field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn + 1
  wrf%var_size(3,ind) = wrf%bt
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +   &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** w field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = wrf%bt + 1
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +    &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** geopotential field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = wrf%bt + 1
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

! we're finished here

!**********************************************************************

end subroutine static_init_model

!****************************1***********************************

subroutine handle_err(ifn,status)
use netcdf
implicit none
integer :: ifn, status

write(6,*) ' error for netcdf function ',ifn
write(6,*) ' status code = ',status
write(6,'(a80)') nf90_strerror(status)

stop

end subroutine handle_err

!**********************************************************************

subroutine netcdf_read_write_var( variable, ncid, var_id, var,          &
                                  start, count, stride, map, in_or_out, debug, ndims )

use netcdf
implicit none
integer :: ncid, var_id, ndims
real, dimension(ndims) :: var
character (len=6) :: in_or_out
integer, dimension(ndims) :: start, count, stride, map
integer :: status
character (len=*) :: variable
logical :: debug

if(debug) write(6,*) ' var for io is ',variable
status = nf90_inq_varid(ncid, variable, var_id)
if(debug) write(6,*) variable, ' id = ',var_id
if (status /= nf90_noerr) call handle_err(4,status)

if( in_or_out(1:5) == "INPUT" ) then
  if(debug) write(6,*) ' call netcdf read ', ncid, var_id
  status = nf90_get_var(ncid, var_id, var, start, count, stride )
  if(debug) write(6,*) ' returned netcdf read '
else if( in_or_out(1:6) == "OUTPUT" ) then
  status = nf90_put_var(ncid, var_id, var, start, count, stride, map)
else
  write(6,*) ' unknown IO function for var_id ',var_id, in_or_out
  write(6,*) ' error stop in netcdf_read_write_var in wrf model_mod '
  stop
end if
if (status /= nf90_noerr) call handle_err(100,status)

end subroutine netcdf_read_write_var

!#######################################################################

function get_model_size()

implicit none

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

implicit none

integer, intent(in) :: index_in
! Temporary kluge of location type
!integer, intent(in) :: location
type(location_type), intent(out) :: location
integer, intent(out), optional :: var_type
integer :: var_type_out

integer :: index, ip, jp, kp
integer :: nz, ny, nx
logical :: var_found
real    :: lon, lat, lev

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
        lat = wrf%latitude(i,jp)  - 0.5*(wrf%latitude(2,jp)-wrf%latitude(1,jp))
     else if(ip == nx) then
        lon = wrf%longitude(nx-1,jp) + 0.5*(wrf%longitude(nx-1,jp)-wrf%longitude(nx-2,jp))
        lat = wrf%latitude(nx-1,jp)  + 0.5*(wrf%latitude(nx-1,jp)-wrf%latitude(nx-2,jp))
     else
        lon = 0.5*(wrf%longitude(ip,jp)+wrf%longitude(ip-1,jp))
        lat = 0.5*(wrf%latitude(ip,jp)+wrf%latitude(ip-1,jp))
     end if

   else if (var_type_out == TYPE_V) then

     if(jp == 1)  then
       lon = wrf%longitude(ip,1) - 0.5*(wrf%longitude(ip,2)-wrf%longitude(ip,1))
       lat = wrf%latitude(ip,1)  - 0.5*(wrf%latitude(ip,2)-wrf%latitude(ip,1))
     else if(jp == ny) then
       lon = wrf%longitude(ip,ny-1) + 0.5*(wrf%longitude(ip,ny-1)-wrf%longitude(ip,ny-2))
       lat = wrf%latitude(ip,ny-1)  + 0.5*(wrf%latitude(ip,ny-1)-wrf%latitude(ip,ny-2))
     else 
       lon = 0.5*(wrf%longitude(ip,jp)+wrf%longitude(ip,jp-1))
       lat = 0.5*(wrf%latitude(ip,jp) +wrf%latitude(ip,jp-1))
     end if

   else  ! regularily staggered variable

     lon = wrf%longitude(ip,jp)
     lat = wrf%latitude(ip,jp)

   end if

   if (lon < 0.0_r8) lon = lon + 360.0

   lev = float(kp)

   if(debug) write(6,*) 'lon, lat, lev: ',lon, lat, lev

   location = set_location(lon, lat, lev)

   if(present(var_type)) var_type = var_type_out

end subroutine get_state_meta_data

!#######################################################################

function model_interpolate(x, location, type)
!!!function model_interpolate(x, lon, lat, level, type)

implicit none

!!! EVENTUALLY NEED TO DO SOMETHING WITH TYPE; BUT FOR NOW THIS JUST FINDS
!!! THE HORIZONTAL PART OF THE INTPERPOLATION??? SHOULD ARGUMENT BE A HYBRID
!!! LOCATION TYPE VARIABLE HERE???

real :: model_interpolate
real, intent(in) :: x(:)
type(location_type), intent(in) :: location
integer, intent(in) :: type

write(6,*) ' need to build function interpolate '
write(6,*) ' error step '
stop

model_interpolate = 0.

end function model_interpolate

!#######################################################################

function get_val(x, lon_index, lat_index, level, type)

implicit none

real :: get_val
real, intent(in) :: x(:)
integer, intent(in) :: lon_index, lat_index, level, type

write(6,*) ' need to build function get_val '
write(6,*) ' error step '
stop

get_val = 0.

end function get_val

!#######################################################################

subroutine model_get_close_states(o_loc, radius, number, indices, dist)

implicit none

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

max_size = wrf%we*wrf%sn + (wrf%we+1)*wrf%sn + wrf%we*(wrf%sn+1)
allocate(lon_ind(max_size), lat_ind(max_size), close_dist(max_size))

! Look for close grid points on the horizontal grid (2D)

call grid_close_states( o_loc, wrf%latitude, wrf%longitude, radius,  &
                        num, lon_ind, lat_ind, close_dist, u_pts, v_pts, p_pts )
!!$write(*,*) ' back from first grid_close_states num = ', num
!!$write(*,*) ' p_pts, u_pts, v_pts ', p_pts, u_pts, v_pts
!!$
!!$write(*,*) ' '
!!$write(*,*) ' p_pt, i, j, dist '
!!$do i = 1, p_pts
!!$  write(*,*) i, lon_ind(i), lat_ind(i), close_dist(i)
!!$enddo
!!$
!!$write(*,*) ' '
!!$write(*,*) ' u_pt, i, j, dist '
!!$do i = p_pts+1, p_pts+u_pts
!!$  write(*,*) i, lon_ind(i), lat_ind(i), close_dist(i)
!!$enddo
!!$
!!$write(*,*) ' '
!!$write(*,*) ' v_pt, i, j, dist '
!!$do i = p_pts+u_pts+1, p_pts+u_pts+v_pts
!!$  write(*,*) i, lon_ind(i), lat_ind(i), close_dist(i)
!!$enddo

! next check vertical spacing, find out how many points
! we have in 3D

   ! w+gz + t + moist + mu
max_size = p_pts*( 2*(wrf%bt+1) + (1 + wrf%n_moist)*wrf%bt + 1)

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

do k = 1, wrf%bt+1
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

!      write(6,*) ' ind, u, i,j,k,index ',num_total,ii,jj,k,indices(num_total)
!      write(6,*) ' dist, radius ',close_dist(i),radius

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

implicit none
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

implicit none

type(location_type), intent(in) :: o_loc
integer, intent(inout) :: num
integer, intent(out) :: close_lon_ind(:), close_lat_ind(:)
integer, intent(out) :: u_pts, v_pts, p_pts
real(r8), intent(in) :: radius_in
real(r8), intent(out) :: close_dist(:)
real(r8), dimension(3) :: loc_array
real(r8) :: o_lon, o_lat

real(r8) :: radius
real, intent(in) :: lat(:,:), lon(:,:)
real :: rad, radn, dxr, dyr, sdx, sdy, gdist
integer :: i_closest, j_closest, ixmin, jymin, ixmax, jymax

integer :: i, j, n, m
real(r8), parameter :: r_earth = 6.37e+06 ! earth radius in meters
logical, parameter :: debug=.false.

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

implicit none
type(location_type), intent(in) :: o_loc
integer, intent(in) :: i,j,k,var_type

type(location_type) :: loc

real :: get_dist_wrf

real :: long, lat


!  get distance for input var_type

   if(var_type == type_u) then

     if (i == 1) then
       long = wrf%longitude(i,j)   &
                - 0.5*(wrf%longitude(i+1,j)-wrf%longitude(i,j))
       lat  = wrf%latitude(i,j)   &
                - 0.5*(wrf%latitude(i+1,j)-wrf%latitude(i,j))
     else if (i == wrf%we + 1) then
       long = wrf%longitude(i-1,j) &
                + 0.5*(wrf%longitude(i-1,j)-wrf%longitude(i-2,j))
       lat  = wrf%latitude(i-1,j) &
                + 0.5*(wrf%latitude(i-1,j)-wrf%latitude(i-2,j))
     else
       long = 0.5*(wrf%longitude(i,j)+wrf%longitude(i-1,j))
       lat  = 0.5*(wrf%latitude(i,j)+wrf%latitude(i-1,j))
     end if

!     write(6,*) ' calling set location, u point,i,j ',i,j
!     write(6,*) ' lon, lat, lev ',long,lat,float(k)-0.5

     if (long < 0.0_r8) long = long + 360.0

     loc = set_location(long,lat,float(k)-0.5)
     get_dist_wrf = get_dist(loc,o_loc)

   else if( var_type == type_v) then

     if (j == 1) then
       long = wrf%longitude(i,j)   &
                - 0.5*(wrf%longitude(i,j+1)-wrf%longitude(i,j))
       lat  = wrf%latitude(i,j)   &
                - 0.5*(wrf%latitude(i,j+1)-wrf%latitude(i,j))
     else if (j == wrf%sn + 1) then
       long = wrf%longitude(i,j-1) &
                + 0.5*(wrf%longitude(i,j-1)-wrf%longitude(i,j-2))
       lat  = wrf%latitude(i,j-1) &
                + 0.5*(wrf%latitude(i,j-1)-wrf%latitude(i,j-2))
     else
       long = 0.5*(wrf%longitude(i,j)+wrf%longitude(i,j-1))
       lat  = 0.5*(wrf%latitude(i,j)+wrf%latitude(i,j-1))
     end if

!     write(6,*) ' calling set location, point 3,i,j ',i,j
!     write(6,*) ' lon, lat, lev ',long,lat,float(k)-0.5

     if (long < 0.0_r8) long = long + 360.0

     loc = set_location(long,lat,float(k)-0.5)

     get_dist_wrf = get_dist(loc,o_loc)

   else if( (var_type == type_w ) .or. &
            (var_type == type_gz) .or. &
            (var_type == type_mu)     ) then

!     write(6,*) ' lon, lat, lev ',wrf%longitude(i,j),wrf%latitude(i,j),float(k)-1.0

      long = wrf%longitude(i,j)
      lat = wrf%latitude(i,j)

     if (long < 0.0_r8) long = long + 360.0

     loc = set_location(long,lat,float(k)-1.0)
     get_dist_wrf = get_dist(loc,o_loc)

   else

!     write(6,*) ' calling set location, p point, i, j ',i,j
!     write(6,*) ' lon, lat, lev ',wrf%longitude(i,j),wrf%latitude(i,j),float(k)-0.5

      long = wrf%longitude(i,j)
      lat = wrf%latitude(i,j)

     if (long < 0.0_r8) long = long + 360.0

     loc = set_location(long,lat,float(k)-0.5)
     get_dist_wrf = get_dist(loc,o_loc)

   end if

!   get_dist_wrf = get_dist_wrf*180./3.1415928

   end function get_dist_wrf

!---------------------------------

function nc_write_model_atts( ncFileID ) result (ierr)
!-----------------------------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! A. Caya May 7 2003
!

use typeSizes
use netcdf
implicit none

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!-----------------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: longDimID, latDimID, levDimID, MemberDimID
integer :: longVarID, latVarID, levVarID, StateVarID
integer :: StateVarDimID, StateVarVarID, TimeDimID
integer :: i

!-----------------------------------------------------------------------------------------

ierr = 0     ! assume normal termination

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_Redef(ncFileID))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies 
!-------------------------------------------------------------------------------

call check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID))
call check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID))

if ( TimeDimID /= unlimitedDimId ) then
  write(*,*)'ERROR: nc_write_model_atts: Time      dimension is ',TimeDimID
  write(*,*)'ERROR: nc_write_model_atts: unlimited dimension is ',unlimitedDimId
  write(*,*)'ERROR: they must be the same.'
  stop
endif

!-------------------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-------------------------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=wrf%model_size, dimid = StateVarDimID))

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source",source))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "DX", wrf%dx))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "DY", wrf%dy))

! how about namelist input? might be nice to save ...

!-------------------------------------------------------------------------------
! Define the dimensions IDs
!-------------------------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="west_east",   len = wrf%we, dimid = longDimID)) 
call check(nf90_def_dim(ncid=ncFileID, name="south_north", len = wrf%sn, dimid = latDimID)) 
call check(nf90_def_dim(ncid=ncFileID, name="bottom_top",  len = wrf%bt,  dimid =  levDimID)) 

! should implement "trajectory-like" coordinate defn ... a'la section 5.4, 5.5 of CF standard
! call check(nf90_def_dim(ncid=ncFileID, name="locationrank", &
!   len = LocationDims, dimid = LocationDimID))

!-------------------------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!-------------------------------------------------------------------------------

! Longitudes
call check(nf90_def_var(ncFileID, name="west_east", xtype=nf90_double, &
                                               dimids=LongDimID, varid=LongVarID) )
call check(nf90_put_att(ncFileID, LongVarID, "long_name", "longitude"))
call check(nf90_put_att(ncFileID, LongVarID, "cartesian_axis", "X"))
call check(nf90_put_att(ncFileID, LongVarID, "units", "degrees_east"))
call check(nf90_put_att(ncFileID, LongVarID, "valid_range", (/ -180.0_r8, 180.0_r8 /)))

! Latitudes
call check(nf90_def_var(ncFileID, name="south_north", xtype=nf90_double, &
                                               dimids=LatDimID, varid=LatVarID) )
call check(nf90_put_att(ncFileID, LatVarID, "long_name", "latitude"))
call check(nf90_put_att(ncFileID, LatVarID, "cartesian_axis", "Y"))
call check(nf90_put_att(ncFileID, LatVarID, "units", "degrees_north"))
call check(nf90_put_att(ncFileID, LatVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)))

! (Common) grid levels
call check(nf90_def_var(ncFileID, name="bottom_top", xtype=nf90_double, &
                                                dimids=levDimID, varid=levVarID) )
call check(nf90_put_att(ncFileID, levVarID, "long_name", "level"))
call check(nf90_put_att(ncFileID, levVarID, "cartesian_axis", "Z"))
call check(nf90_put_att(ncFileID, levVarID, "units", "???"))

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create attributes for the state vector 
   !----------------------------------------------------------------------------

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

!-------------------------------------------------------------------------------
! Leave define mode so we can actually fill the variables.
!-------------------------------------------------------------------------------

call check(nf90_enddef(ncfileID))

   ! Fill the state variable coordinate variable
   call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,wrf%model_size) /) ))

endif

!-------------------------------------------------------------------------------
! Fill the variables
!-------------------------------------------------------------------------------

call check(nf90_put_var(ncFileID, LongVarID, wrf%longitude(1:wrf%we,1) ))
call check(nf90_put_var(ncFileID, LatVarID, wrf%latitude(1,1:wrf%sn) ))

call check(nf90_put_var(ncFileID,  levVarID, wrf%dn(1:wrf%bt) ))
!call check(nf90_put_var(ncFileID,  levVarID, (/ (i,i=1,   wrf%bt) /) ))

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call check(nf90_sync(ncFileID))

write (*,*)'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...'

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus

    if(istatus /= nf90_noerr) then
       print *,'model_mod:nc_write_model_atts'
       print *, trim(nf90_strerror(istatus))
       ierr = istatus
       stop
    end if
  end subroutine check

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!-----------------------------------------------------------------------------------------
! Writes the model-specific variables to a netCDF file
! TJH 25 June 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!                                                                                
! There are two different (staggered) 3D grids being used simultaneously here. 
! The routines "prog_var_to_vector" and "vector_to_prog_var", 
! packs the prognostic variables into
! the requisite array for the data assimilation routines. That routine
! is the basis for the information stored in the netCDF files.
!
! TemperatureGrid : surface pressure  vars%ps(tis:tie, tjs:tje) 
!                 : temperature       vars%t (tis:tie, tjs:tje, klb:kup)
!                 : tracers           vars%r (tis:tie, tjs:tje, klb:kub, 1:vars%ntrace)
! VelocityGrid    : u                 vars%u (vis:vie, vjs:vje, klb:kub) 
!                 : v                 vars%v (vis:vie, vjs:tje, klb:kup)
!
! So there are six different dimensions and five different variables as long as
! simply lump "tracers" into one. 

use typeSizes
use netcdf
implicit none

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!-----------------------------------------------------------------------------------------
real, dimension(SIZE(statevec)) :: x
!!$type(prog_var_type) :: Var

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, psVarID, tVarID, rVarID, uVarID, vVarID
integer :: tis, tie, tjs, tje       ! temperature grid start/stop
integer :: vis, vie, vjs, vje       ! velocity    grid start/stop
integer :: kub, klb
integer :: nTmpI, nTmpJ, nVelI, nVelJ, nlev, ntracer, i

ierr = 0     ! assume normal termination

!-------------------------------------------------------------------------------
! Get the bounds for storage on Temp and Velocity grids
! ?JEFF why can't I use the components of the prog_var_type?  
! More precisely, why doesn't prog_var_type drag around the necessary
! indices instead of just the extents?
!-------------------------------------------------------------------------------

!!$tis = Dynam%Hgrid%Tmp%is; tie = Dynam%Hgrid%Tmp%ie
!!$tjs = Dynam%Hgrid%Tmp%js; tje = Dynam%Hgrid%Tmp%je
!!$vis = Dynam%Hgrid%Vel%is; vie = Dynam%Hgrid%Vel%ie
!!$vjs = Dynam%Hgrid%Vel%js; vje = Dynam%Hgrid%Vel%je
!!$kub = Var_dt%kub
!!$klb = Var_dt%klb
!!$
!!$nTmpI   = tie - tis + 1
!!$nTmpJ   = tje - tjs + 1
!!$nlev    = Var_dt%kub - Var_dt%klb + 1
!!$ntracer = Var_dt%ntrace 
!!$nVelI   = vie - vis + 1
!!$nVelJ   = vje - vjs + 1

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! then get all the Variable ID's we need.
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

!!$if ( output_state_vector ) then

   call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)))                               

!!$else
!!$   
!!$   !----------------------------------------------------------------------------
!!$   ! Fill the variables
!!$   ! TemperatureGrid : surface pressure  Var%ps(tis:tie, tjs:tje) 
!!$   !                 : temperature       Var%t (tis:tie, tjs:tje, klb:kub)
!!$   !                 : tracers           Var%r (tis:tie, tjs:tje, klb:kub, 1:vars%ntrace)
!!$   ! VelocityGrid    : u                 Var%u (vis:vie, vjs:vje, klb:kub) 
!!$   !                 : v                 Var%v (vis:vie, vjs:tje, klb:kub)
!!$   !----------------------------------------------------------------------------
!!$
!!$   x = statevec ! Unfortunately, have to explicity cast it ...
!!$                ! the filter uses a type=double,
!!$                ! the vector_to_prog_var function expects a single.
!!$   call vector_to_prog_var(x, get_model_size(), Var)
!!$   
!!$   
!!$   call check(NF90_inq_varid(ncFileID, "ps", psVarID))
!!$   call check(nf90_put_var( ncFileID, psVarID, Var%ps(tis:tie, tjs:tje), &
!!$                            start=(/ 1, 1, copyindex, timeindex /) ))
!!$
!!$   call check(NF90_inq_varid(ncFileID,  "t",  tVarID))
!!$   call check(nf90_put_var( ncFileID,  tVarID, Var%t( tis:tie, tjs:tje, klb:kub ), &
!!$                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
!!$
!!$   call check(NF90_inq_varid(ncFileID,  "u",  uVarID))
!!$   call check(nf90_put_var( ncFileID,  uVarId, Var%u( vis:vie, vjs:vje, klb:kub ), &
!!$                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
!!$
!!$   call check(NF90_inq_varid(ncFileID,  "v",  vVarID))
!!$   call check(nf90_put_var( ncFileID,  vVarId, Var%v( vis:vie, vjs:vje, klb:kub ), &
!!$                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
!!$
!!$   if ( ntracer > 0 ) then
!!$      call check(NF90_inq_varid(ncFileID,  "r",  rVarID))
!!$      call check(nf90_put_var( ncFileID,  rVarID, &
!!$                    Var%r( tis:tie, tjs:tje, klb:kub, 1:ntracer ), & 
!!$                   start=(/   1,       1,       1,     1,    copyindex, timeindex /) ))
!!$   endif
!!$endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
write (*,*)'netCDF file is synched ...'

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus

    if(istatus /= nf90_noerr) then
      print *,'model_mod:nc_write_model_vars'
      print *, trim(nf90_strerror(istatus))
      stop
    end if
  end subroutine check

end function nc_write_model_vars

!-------------------------------

!  public stubs

!**********************************************

subroutine adv_1step(x, Time)

! Does single time-step advance for B-grid model with vector state as
! input and output. This is a modified version of subroutine atmosphere
! in original bgrid_solo_atmosphere driver.

implicit none

real, intent(inout) :: x(:)
! Time is needed for more general models like this; need to add in to 
! low-order models
type(time_type), intent(in) :: Time

end subroutine adv_1step

!**********************************************

subroutine end_model()
end subroutine end_model

!**********************************************

subroutine init_time(i_time)

implicit none

! For now returns value of Time_init which is set in initialization routines.

type(time_type), intent(out) :: i_time

end subroutine init_time

!**********************************************

subroutine init_conditions(x)

implicit none

! Reads in restart initial conditions from B-grid and converts to vector

! Following changed to intent(inout) for ifc compiler;should be like this
real, intent(inout) :: x(:)

end subroutine init_conditions

!#######################################################################
end module model_mod

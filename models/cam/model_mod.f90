! netCDF filename; where will this come from in DART?
!                  it's created by CAM, and written out to a file somewhere; read it in?

! Do we need other functionality form bgrid model_mod?

!#######################################################################
module model_mod

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!
! method: Read CAM 'initial' file (netCDF format).
!         Reform fields into a state vector.
!         (Let DART modify values state vector.)
!         Reform state vector back into CAM fields.
!         Replace those fields on the CAM initial file with the new values,
!         preserving all other information on the file.
!
! author: Kevin Raeder 2/14/03
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
!
!----------------------------------------------------------------------

use netcdf
use types_mod
use time_manager_mod, only : time_type, set_time
use location_mod         , only: location_type, get_location, set_location, get_dist, &
                                 LocationDims, LocationName, LocationLName


implicit none
private

public prog_var_to_vector, vector_to_prog_var, read_cam_init, &
   read_cam_init_size, &
   write_cam_init, get_model_size, static_init_model, &
   get_state_meta_data, get_model_time_step, model_interpolate, &
   init_conditions, init_time, adv_1step, end_model, &
   model_get_close_states, nc_write_model_atts, &
   vars, x, &
   TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q, TYPE_TRACER


!-----------------------------------------------------------------------
! Public definition of variable types
integer, parameter :: TYPE_PS = 0, TYPE_T = 1, TYPE_U = 2, TYPE_V = 3, TYPE_Q = 4, TYPE_TRACER = 5


!----------------------------------------------------------------------

! Global storage for describing cam model class
integer :: model_size, num_lons, num_lats, num_levs
type(time_type) :: Time_step_atmos

integer, parameter :: n3dflds=4     ! # of 3d fields to read from file
!                                   including Q, but not tracers
integer, parameter :: pcnst  =0     ! advected tracers (don't include Q in this)
integer, parameter :: pnats  =0     ! nonadvected tracers
integer, parameter :: n2dflds=1     ! # of 2d fields to read from file
! derived parameters
integer, parameter :: n3tflds  = n3dflds+pcnst+pnats   ! # fields to read
integer, parameter :: nflds  = n3tflds+n2dflds         ! # fields to read

! Arrays to store lat and lon indices
real(r8), allocatable :: lons(:), lats(:)

! This should be local inside the conversion routines???
real(r8), allocatable :: vars(:,:,:,:)                 !3d and 2D variables read from CAM

! list variables according to the category to which they belong, 
! in the order the categories appear above (n3dflds,pcnst,pnats,n2dflds).
! from ncdump of CAM standard configuration initial file:
!   lat = 64 ; lon = 128 ; lev = 26 ;
!   hyai, hybi, hyam, hybm, gw
!   U, V, T, Q, PS, 
!   (names of advected and nonadv constituents)
!   PHIS, SGH, LANDM, 
!   TS, TSICE, TS1, TS2, TS3, TS4
!   SNOWHICE, LANDFRAC, CWAT
character (len=8),dimension(nflds) :: cflds = &
          (/'U       ','V       ','T       ','Q       ','PS      ' /)
integer, dimension(nflds) :: nlevs 
data nlevs/ n3tflds*0, n2dflds*1/

! netCDF filename; where will this come from in DART?
!                  it's created by CAM, and written out to a file somewhere; read it in?
character (len=128) :: filename = '/home/raeder/DAI/test.nc'

!----------------------------------------------------------------------
! dimension DART arrays

real(r8), allocatable :: x(:)          ! state vector to pass to DART, size (siz)

!----------------------------------------------------------------------
! netCDF parameters and variables
integer londimid, levdimid, latdimid ! Dimension ID's
integer ncfileid                     ! netCDF file ID 
integer ncfldid                      ! netCDF field ID 

! kdr; uncomment for installation into DART
!-----------------------------------------------------------------------
! let CVS fill strings ... DO NOT EDIT ...
! 
! character(len=128) :: version = "$Id$"
! character(len=128) :: tag = "$Name$"
! 
! character(len=128) :: &
!    source = "$Source$", &
!    revision = "$Revision$", &
!    revdate  = "$Date$"
! 
!-----------------------------------------------------------------------
!---- namelist (saved in file input.nml) ----
!-----------------------------------------------------------------------

contains



!#######################################################################

subroutine read_cam_init_size(file_name, num_lons, num_lats, num_levs)

! Gets the number of lons, lats and levels from a CAM init netcdf file
! in file_name

character(len = *), intent(in) :: file_name
integer, intent(out) :: num_lons, num_lats, num_levs

character (len=NF90_MAX_NAME) :: clon,clat,clev

write(*, *) 'file_name in read_cam_init is ', trim(file_name)

!----------------------------------------------------------------------
! read CAM 'initial' file domain info
call check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid))

! get dimension 'id's
call check(nf90_inq_dimid(ncfileid, 'lon', londimid))
call check(nf90_inq_dimid(ncfileid, 'lat', latdimid))
call check(nf90_inq_dimid(ncfileid, 'lev', levdimid))

! get dimension sizes
call check(nf90_inquire_dimension(ncfileid, londimid, clon , num_lons ))
call check(nf90_inquire_dimension(ncfileid, latdimid, clat , num_lats ))
call check(nf90_inquire_dimension(ncfileid, levdimid, clev , num_levs ))

end subroutine read_cam_init_size

!#######################################################################
subroutine read_cam_init

!----------------------------------------------------------------------
! Local workspace
integer :: i,j,ifld  ! grid and constituent indices
integer :: siz, plon, plat, plev

character (len=NF90_MAX_NAME) :: clon,clat,clev

!----------------------------------------------------------------------
! read CAM 'initial' file domain info
call check(nf90_open(path = trim(filename), mode = nf90_write, ncid = ncfileid))

! get dimension 'id's
call check(nf90_inq_dimid(ncfileid, 'lon', londimid))
call check(nf90_inq_dimid(ncfileid, 'lat', latdimid))
call check(nf90_inq_dimid(ncfileid, 'lev', levdimid))

! get dimension sizes
call check(nf90_inquire_dimension(ncfileid, londimid, clon , plon ))
call check(nf90_inquire_dimension(ncfileid, latdimid, clat , plat ))
call check(nf90_inquire_dimension(ncfileid, levdimid, clev , plev ))

! allocate space for arrays based on domain size
allocate ( vars(plon,plev,plat,nflds) )
siz = (plon*plat) * ((n3dflds+pcnst+pnats)*plev + n2dflds)
allocate ( x(siz))
PRINT*,'in read_cam_init siz = ',siz

! specify # vertical levels for 3D fields (2d have already been filled)
i=1
do while (nlevs(i)==0 .and. i<=nflds)
   nlevs(i) = plev
   i=i+1
end do
WRITE(*,*) 'nlevs = ',(nlevs(j),j=1,nflds)

! read CAM 'initial' file fields desired
! 3d fields
do ifld=1,n3tflds
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 4D; lon, lev, lat, TIME(=1) 
   call check(nf90_get_var(ncfileid, ncfldid, vars(:,:,:,ifld) &
             ,start=(/1,1,1,1/) ,count=(/plon,nlevs(ifld),plat,1/) ))
end do

!2d fields
do ifld=n3tflds+1,n3tflds+n2dflds
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 3D; lon, lat, TIME(=1)
   call check(nf90_get_var(ncfileid, ncfldid, vars(:,1,:,ifld) &
             ,start=(/1,1,1/) ,count=(/plon,plat,1/) ))
end do

return
end subroutine read_cam_init

!#######################################################################
subroutine prog_var_to_vector(vars, x, siz)

implicit none

integer, intent(in) :: siz

real(r8), intent(in), dimension(num_lons, num_levs, num_lats, nflds) :: vars
real(r8), intent(out),dimension(siz) :: x

integer :: i, j, k, nf, nt, n0, index

! Start copying fields to straight vector
index = 0
do i = 1, num_lons
   do j = 1, num_lats
!     u,v,t,q, and tracers at successively lower levels
      n0 = n3dflds
      do k = 1, num_levs
         do nf=1,n3dflds
            index = index + 1
            x(index) = vars(i, k, j, nf)
         end do
         do nt = n0+1, n0+pcnst+pnats
            IF (i==1 .and. j==1) PRINT*,'filling tracers'
            index = index + 1
            x(index) = vars(i, k, j, nt)
         end do
      end do
!     Surface pressure and other 2d flds are last
      n0 = n3dflds+pcnst+pnats
      do nf=n0+1,n0+n2dflds
         index = index + 1
         x(index) = vars(i, 1, j, nf)
      end do
   end do
end do

! Temporary check
if(index /= siz) then
   write(*, *) 'prog_var_to_vector bad index sum '
   write(*, *) 'index, siz ', index, siz
   stop
endif

return
end subroutine prog_var_to_vector

!#######################################################################

subroutine vector_to_prog_var(x, siz, vars) 

implicit none

 integer, intent(in) :: siz
 real(r8), intent(in),dimension(siz) :: x
 real(r8), dimension(num_lons, num_levs, num_lats,nflds), intent(out) :: vars

integer :: i, j, k, nf, nt, n0, index

! Start copying fields from straight vector
index = 0
do i = 1, num_lons
   do j = 1, num_lats
!     u,v,t,q  and tracers at successive levels
      n0 = n3dflds
      do k = 1, num_levs
         do nf = 1,n3dflds
            index = index + 1
            vars(i, k, j, nf) = x(index)
         end do 
         do nt = n0+1, n0+pcnst+pnats
            index = index + 1
            vars(i, k, j, nt) = x(index)
         end do
      end do
! Surface pressure and other 2d fields are last
      n0 = n3dflds+pcnst+pnats
      do nf=n0+1,n0+n2dflds
         index = index + 1
         vars(i, 1, j, nf) = x(index)
      end do
   end do
end do

! Temporary check
if(index /= siz) then
   write(*, *) 'prog_var_to_vector bad index sum '
   write(*, *) 'index, siz ', index, siz
   stop
endif

return
end subroutine vector_to_prog_var

!#######################################################################
subroutine write_cam_init

implicit none

integer ifld

! write CAM 'initial' file fields that have been updated
do ifld=1,n3tflds
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid, vars(:,:,:,ifld) &
             ,start=(/1,1,1,1/) ,count=(/num_lons,nlevs(ifld),num_lats,1/) ))
end do
do ifld=n3tflds+1 ,n3tflds+n2dflds
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid, vars(:,nlevs(ifld),:,ifld) &
             ,start=(/1,1,1/) ,count=(/num_lons,num_lats,1/) ))
end do

call check(nf90_close(ncfileid))

return
end subroutine write_cam_init

!#######################################################################
subroutine check(status)

integer, intent ( in) :: status

if (status /= nf90_noerr) then
   print *, trim(nf90_strerror(status))
end if

return
end subroutine check

!#######################################################################

function get_model_size()

implicit none

integer :: get_model_size

get_model_size = model_size

end function get_model_size

!#######################################################################

subroutine static_init_model()

! INitializes class data for CAM model (all the stuff that needs to
! be done once. For now, does this by reading info from a fixed
! name netcdf file. Need to make this file a namelist parameter
! at some point.

implicit none

! Get num lons, lats and levs from netcdf and put in global storage
call read_cam_init_size('H12-24icl.nc', num_lons, num_lats, num_levs)

! Compute the model size (Need to be sure where all these come from
! and if they themselves need to be initialized here)

! Need to set Time_step_atmos here to 1 hour for now 
Time_step_atmos = set_time(3600, 0)

! Compute overall model size and put in global storage
model_size = num_lons * num_lats * (n2dflds + num_levs * &
   (n3dflds + pcnst + pnats))

write(*, *) 'CAM size initialized as ', model_size

end subroutine static_init_model

!#######################################################################

subroutine adv_1step(x, Time)

! Does single time-step advance for B-grid model with vector state as
! input and output. This is a modified version of subroutine atmosphere
! in original bgrid_solo_atmosphere driver.

implicit none

real, intent(inout) :: x(:)
! Time is needed for more general models like this; need to add in to
! low-order models
type(time_type), intent(in) :: Time

! This is a no-op for CAM; only asynch integration

end subroutine adv_1step

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
! Types for this CAM model are, TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q

implicit none

integer, intent(in) :: index_in
! Temporary kluge of location type
!integer, intent(in) :: location
type(location_type), intent(out) :: location
integer, intent(out), optional :: var_type

integer :: index, num_per_col, col_num, col_elem, lon_index, lat_index
real(r8) :: lon, lat, lev
integer :: local_var_type, var_type_temp


! Easier to compute with a 0 to size - 1 index
index = index_in - 1

! Compute number of items per column
num_per_col = num_levs * (n3dflds + pcnst + pnats) + n2dflds

! What column is this index in
col_num = index_in / num_per_col 
col_elem = index - col_num * num_per_col

! What lon and lat index for this column
lon_index = col_num / num_lats
lat_index = col_num - lon_index * num_lats

! Get actual lon lat values from static_init arrays ???
lon = lons(lon_index + 1)
lat = lats(lat_index + 1)

! Now figure out which beast in column this is
! Surface pressure is the last element
lev = col_elem / (n3dflds + pcnst + pnats) + 1
! If we're out of the levels
if(lev > num_levs) then
   local_var_type = TYPE_PS
   lev = -1
else
   var_type_temp = mod(col_elem - 1, n3dflds + pcnst + pnats)
   if(var_type_temp == 0) then
      local_var_type = TYPE_U
   else if(var_type_temp == 1) then
      local_var_type = TYPE_V
   else if(var_type_temp == 2) then
      local_var_type = TYPE_T
   else
      local_var_type = TYPE_Q
   endif
endif



!write(*, *) 'lon, lat, and lev ', lon, lat, lev
location = set_location(lon, lat, lev)

! If the type is wanted, return it
if(present(var_type)) var_type = local_var_type


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

! No op for now
if(1 == 1) then
   write(*, *) 'STOP: CAM Model interpolate is not implemented yet'
   stop
endif
model_interpolate = 0.0

end function model_interpolate

!#######################################################################


function get_model_time_step()
!------------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

! Time_step_atmos is global static storage
get_model_time_step =  Time_step_atmos

end function get_model_time_step

!#######################################################################

subroutine end_model()

implicit none

! At some point, this stub should coordinate with atmosphere_end but
! that requires an instance variable.

end subroutine end_model

!#######################################################################

subroutine init_time(i_time)

implicit none

! For now returns value of Time_init which is set in initialization routines.

type(time_type), intent(out) :: i_time

!Where should initial time come from here?
! WARNING: CURRENTLY SET TO 0
i_time = set_time(0, 0)

end subroutine init_time

!#######################################################################

subroutine init_conditions(x)

implicit none

! Reads in restart initial conditions from B-grid and converts to vector

! Following changed to intent(inout) for ifc compiler;should be like this
real, intent(inout) :: x(:)

! Can this be a noop for CAM???

end subroutine init_conditions
!#######################################################################

subroutine model_get_close_states(o_loc, radius, number, indices, dist)

implicit none

type(location_type), intent(in) :: o_loc
real(r8), intent(in) :: radius
integer, intent(out) :: number, indices(:)
real(r8), intent(out) :: dist(:)

real(r8) :: loc_array(3), o_lon, o_lat
integer :: nlon, nlat, num, max_size, i, j, num1, nlev
integer :: hsize, num_per_col, grid_size, col_base_index
integer, allocatable :: lon_ind(:), lat_ind(:)
real(r8), allocatable :: close_dist(:)



! Number found starts at 0
number = 0

! Assume that grid size is known from static initialized storage

! Num of close horizontal grid points starts at 0, too
num = 0
! For now, just allocate enough space for all grid points, may want
! to make this smaller at some point for big models.
max_size = nlon * nlat
allocate(lon_ind(max_size), lat_ind(max_size), close_dist(max_size))

! Look for close grid points on the 
call grid_close_states(o_loc, lons, lats, nlon, nlat, radius, &
   num, lon_ind, lat_ind, close_dist)
write(*, *) 'back from grid_close_states num = ', num

! Compute size of grid storage for full levels
hsize = nlon * nlat
num_per_col = nlev * (n3dflds + pcnst + pnats) + n2dflds
grid_size = num_per_col * hsize

! Add all variables in this column to the close list with this distance
do i = 1, num
   col_base_index = ((lon_ind(i) - 1) * nlat + lat_ind(i) - 1) * num_per_col
   do j = 1, num_per_col
      number = number + 1
      if(number <= size(indices)) indices(number) = col_base_index + j
      if(number <= size(dist)) dist(number) = close_dist(i)
   end do
end do

deallocate(lon_ind, lat_ind, close_dist)

end subroutine model_get_close_states

!#######################################################################

subroutine grid_close_states(o_loc, lons, lats, nlon, nlat, radius, &
   num, close_lon_ind, close_lat_ind, close_dist)

! WARNING: THIS ROUTINE DOES NOT FUNCTION PROPERLY FOR LARGE RADII. It
! SACRIFICES GETTING ALL THE POINTS AT THE MOST DISTANT NEGATIVE
! LATITUDE OFFSETS TO AVOID GETTING REDUNDANT POINT (WHICH WOULD BE
! VERY BAD POTENTIALLY). AT SOME POINT SOME TALENTED PROGRAMMER
! SHOULD CLEAN THIS UP. THERE IS ALSO A DANGER OF GETTING
! REDUNDANT POINTS NEAR THE POLES IF THE LONGITUDE DISTRIBUTION
! HAS AN ODD NUMBER OF POINTS SO THAT THE REV_BLON IS NOT 180
! DEGREES DIFFERENT. A CHECK IS IN FOR THIS.

! Finds close state points from a particular grid;

implicit none

type(location_type), intent(in) :: o_loc
integer, intent(in) :: nlon, nlat
real(r8), intent(in) :: lons(nlon), lats(nlat), radius
integer, intent(inout) :: num
integer, intent(inout) :: close_lon_ind(:), close_lat_ind(:)
real(r8), intent(out) :: close_dist(:)

real(r8) :: rev_lon, glat, glon, loc_array(3), o_lon, o_lat, o_lev
real(r8) :: gdist, diff
integer :: blat_ind, blon_ind, i, j, rev_blon_ind, lat_ind, lon_ind
integer :: final_blon_ind
real(r8), parameter :: glev = 1.0
type(location_type) :: loc


! Get the lat and lon from the loc
loc_array = get_location(o_loc)
o_lon = loc_array(1)
o_lat = loc_array(2)

! Get index to closest lat and lon for this observation
blat_ind = get_closest_lat_index(o_lat, lats, nlat)
write(*, *) 'closest latitude in grid is ', blat_ind, lats(blat_ind)

blon_ind = get_closest_lon_index(o_lon, lons, nlon)
write(*, *) 'closest longitude in grid is ', blon_ind, lons(blon_ind)

! Also may need to work with the longitude on opposite side for pole wrap
rev_blon_ind = blon_ind + nlon / 2
if(rev_blon_ind > nlon) rev_blon_ind = rev_blon_ind - nlon

! If the rev_blon and blon are not 180 degrees apart this algorithm
! could produce redundant close points, don't want that to happen.
diff = dabs(lons(rev_blon_ind) - lons(blon_ind))
if(dabs(diff - 180.0) > 0.000001) then
! Actually, it can fail anyway because of the wraparound longitude search
   write(*, *) 'ERROR in subroutine grid_close_states in bgrid model_mod'
   write(*, *) 'Algorithm can fail if grids do not have even number of '
   write(*, *) 'longitudes'
   stop
endif
write(*, *) 'closest rev long in grid is ', rev_blon_ind, lons(rev_blon_ind)

! Begin a search along the latitude axis in the positive direction
final_blon_ind = blon_ind
do i = 0, nlat / 2
   lat_ind = blat_ind + i
   if(lat_ind > nlat) then
      lat_ind = nlat - (lat_ind - nlat - 1)
! This will change once and for ever after wrapping over pole
      final_blon_ind = rev_blon_ind
   endif
   glat = lats(lat_ind)
! Take care of storage round-off
   if(glat < -90.0) glat = 0.0
   if(glat > 90.0) glat = 90.0
! Search in positive longitude offset
   do j = 0, nlon / 2
      lon_ind = final_blon_ind + j
      if(lon_ind > nlon) lon_ind = lon_ind - nlon
      glon = lons(lon_ind)
      if(glon > 360.0) glon = glon - 360.0
      if(glon < 0.0) glon = glon + 360.0
      loc = set_location(glon, glat, glev)
      gdist = get_dist(loc, o_loc)
      if(gdist <= radius) then
         num = num + 1
         close_lon_ind(num) = lon_ind
         close_lat_ind(num) = lat_ind
         close_dist(num) = gdist
! If radius is too far for closest longitude, no need to search further in lat
      else if (j == 0) then
         goto 11
      else
! Look in negative longitude offset direction next
         goto 21
      endif
   end do
! Search in negative longitude offset
21 do j = 1, (nlon + 1) / 2 - 1
      lon_ind = final_blon_ind - j
      if(lon_ind < 1) lon_ind = nlon + lon_ind
      glon = lons(lon_ind)
      if(glon > 360.0) glon = glon - 360.0
      if(glon < 0.0) glon = glon + 360.0
      loc = set_location(glon, glat, glev)
      gdist = get_dist(loc, o_loc)
      if(gdist <= radius) then
         num = num + 1
         close_lon_ind(num) = lon_ind
         close_lat_ind(num) = lat_ind
         close_dist(num) = gdist
      else
         goto 31
      endif
   end do
31 continue
end do
11 continue


! Search in the negative lat direction, be sure not to be redundant
final_blon_ind = blon_ind
do i = 1, (nlat + 1) / 2 - 1
   lat_ind = blat_ind - i
   if(lat_ind < 1) then
      lat_ind = -1 * lat_ind + 1
      final_blon_ind = rev_blon_ind
   endif
   glat = lats(lat_ind)
! Take care of storage round-off
   if(glat < -90.0) glat = 0.0
   if(glat > 90.0) glat = 90.0
! Search in positive longitude offset
   do j = 0, nlon / 2
      lon_ind = final_blon_ind + j
      if(lon_ind > nlon) lon_ind = lon_ind - nlon
      glon = lons(lon_ind)
      if(glon > 360.0) glon = glon - 360.0
      if(glon < 0.0) glon = glon + 360.0
      loc = set_location(glon, glat, glev)
      gdist = get_dist(loc, o_loc)
      if(gdist <= radius) then
         num = num + 1
         close_lon_ind(num) = lon_ind
         close_lat_ind(num) = lat_ind
         close_dist(num) = gdist
! If radius is too far for closest longitude, no need to search further in lat
      else if (j == 0) then
         goto 111
      else
! Look in negative longitude offset direction next
         goto 121
      endif
   end do
! Search in negative longitude offset
121 do j = 1, (nlon + 1) / 2 - 1
      lon_ind = final_blon_ind - j
      if(lon_ind < 1) lon_ind = nlon + lon_ind
      glon = lons(lon_ind)
      if(glon > 360.0) glon = glon - 360.0
      if(glon < 0.0) glon = glon + 360.0
      loc = set_location(glon, glat, glev)
      gdist = get_dist(loc, o_loc)
      if(gdist <= radius) then
         num = num + 1
         close_lon_ind(num) = lon_ind
         close_lat_ind(num) = lat_ind
         close_dist(num) = gdist
      else
         goto 131
      endif
   end do
131 continue
end do
111 continue

end subroutine grid_close_states

!#######################################################################

function nc_write_model_atts( ncFileID ) result (ierr)

use typeSizes
use netcdf
implicit none

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

! Normal null return is 0
ierr = 0


end function nc_write_model_atts

!#######################################################################

function get_closest_lat_index(o_lat, lats, nlat)

implicit none
integer, intent(in) :: nlat
real(r8), intent(in) :: o_lat, lats(nlat)
integer :: get_closest_lat_index

real(r8) :: lat_bot, lat_top, lat_int, diff
integer :: lower_ind

! Find closest lat
lat_bot = lats(1)
lat_top = lats(nlat)
lat_int = lats(2) - lats(1)
if(o_lat <= lat_bot) then
   get_closest_lat_index = 1
else if(o_lat >= lat_top) then
   get_closest_lat_index = nlat
else
   diff = (o_lat - lat_bot) / lat_int
   lower_ind = int(diff) + 1
   if(diff - int(diff) < 0.5) then
      get_closest_lat_index = lower_ind
   else
      get_closest_lat_index = lower_ind + 1
   endif
endif

end function get_closest_lat_index

!#######################################################################

function get_closest_lon_index(o_lon, lons, nlon)

implicit none
integer, intent(in) :: nlon
real(r8), intent(in) :: o_lon, lons(nlon)
integer :: get_closest_lon_index

real(r8) :: diff, lon_bot, lon_top, lon_int
integer :: lower_ind, blon_ind

! Find closest longitude on grid to given longitude
lon_bot = lons(1)
lon_top = lons(nlon)
lon_int = lons(2) - lons(1)
if(o_lon <= lon_bot) then
   diff = (lon_bot - o_lon) / lon_int
   if(diff > 0.5) then
      get_closest_lon_index = nlon
   else
      get_closest_lon_index = 1
   end if
else if(o_lon >= lon_top) then
   diff = (o_lon - lon_top) / lon_int
   if(diff > 0.5) then
      get_closest_lon_index = 1
   else
      get_closest_lon_index = nlon
   end if
else
   diff = (o_lon - lon_bot) / lon_int
   lower_ind = int(diff) + 1
   if(diff - int(diff) < 0.5) then
      get_closest_lon_index = lower_ind
   else
      get_closest_lon_index = lower_ind + 1
   end if
end if

end function get_closest_lon_index

!#######################################################################

end module model_mod



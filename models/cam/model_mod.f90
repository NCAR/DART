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
use utilities_mod, only : file_exist, open_file, check_nml_error, close_file
use time_manager_mod, only : time_type, set_time
use location_mod         , only: location_type, get_location, set_location, get_dist, &
                                 LocationDims, LocationName, LocationLName


implicit none
private

public model_type, prog_var_to_vector, vector_to_prog_var, read_cam_init, &
   read_cam_init_size, init_model_instance, &
   write_cam_init, get_model_size, static_init_model, &
   get_state_meta_data, get_model_time_step, model_interpolate, &
   init_conditions, init_time, adv_1step, end_model, &
   model_get_close_states, nc_write_model_atts, nc_write_model_vars, &
   TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_Q, TYPE_TRACER

!-----------------------------------------------------------------------
! let CVS fill strings ... DO NOT EDIT ...
 
character(len=128) :: version = "$Id$"
character(len=128) :: tag = "$Name$"
 
character(len=128) :: &
   source = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Public definition of variable types
integer, parameter :: TYPE_PS = 0, TYPE_T = 1, TYPE_U = 2, TYPE_V = 3, TYPE_Q = 4, TYPE_TRACER = 5

!----------------------------------------------------------------------
! A type for cam model, very simple for now for conversion only
type model_type
!   private
   real(r8), pointer :: vars_2d(:, :, :)
   real(r8), pointer :: vars_3d(:, :, :, :)
   real(r8), pointer ::  tracers(:, :, :, :)
end type model_type

!-----------------------------------------------------------------------                       
! need a model namelist
! output_state_vector = .true.     results in a "state-vector" netCDF file
! output_state_vector = .false.    results in a "prognostic-var" netCDF file
                                                                                               
logical :: output_state_vector = .true.  ! output state-vector 
                                                                                               
namelist /model_nml/ output_state_vector


!----------------------------------------------------------------------
! File where basic info about model configuration can be found; should be namelist

character(len = 128) :: model_config_file = 'caminput.nc'
!character(len = 128) :: model_config_file = 'T5H0-12icl.cam2.i.0001-09-01-43200.nc'
!----------------------------------------------------------------------

!
! Global storage for describing cam model class
integer :: model_size, num_lons, num_lats, num_levs
type(time_type) :: Time_step_atmos

integer, parameter :: n3dflds=4     ! # of 3d fields to read from file
!                                   including Q, but not tracers
integer, parameter :: pcnst  =0     ! advected tracers (don't include Q in this)
integer, parameter :: pnats  =0     ! nonadvected tracers
integer, parameter :: num_tracers = pcnst + pnats
integer, parameter :: n2dflds=1     ! # of 2d fields to read from file
! derived parameters
integer, parameter :: n3tflds  = n3dflds+pcnst+pnats   ! # fields to read
integer, parameter :: nflds  = n3tflds+n2dflds         ! # fields to read

! Arrays to store lat and lon indices
real(r8), allocatable :: lons(:), lats(:)

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
          (/'PS      ','T       ','U       ','V       ','Q       ' /)
!          (/'U       ','V       ','T       ','Q       ','PS      ' /)

!  ??? Need to understand what Kevin was doing with nlevs; I may be missing something essential
integer, dimension(nflds) :: nlevs 
data nlevs/ n3tflds*0, n2dflds*1/

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
integer :: londimid, levdimid, latdimid, ncfileid

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
subroutine read_cam_init(file_name, var)

implicit none

character(len = *), intent(in) :: file_name
type(model_type), intent(out) :: var

!----------------------------------------------------------------------
! Local workspace
integer :: i,j,ifld  ! grid and constituent indices
integer :: plon, plat, plev

character (len=NF90_MAX_NAME) :: clon,clat,clev
integer :: londimid, levdimid, latdimid, ncfileid, ncfldid

!----------------------------------------------------------------------
! read CAM 'initial' file domain info
call check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid))

! Could do this for storage size error check later
!call check(nf90_inq_dimid(ncfileid, 'lon', londimid))
!call check(nf90_inq_dimid(ncfileid, 'lat', latdimid))
!call check(nf90_inq_dimid(ncfileid, 'lev', levdimid))

! get dimension sizes
!call check(nf90_inquire_dimension(ncfileid, londimid, clon , plon ))
!call check(nf90_inquire_dimension(ncfileid, latdimid, clat , plat ))
!call check(nf90_inquire_dimension(ncfileid, levdimid, clev , plev ))

!!! ??? I'm not sure what this next part is really doing, verify
! specify # vertical levels for 3D fields (2d have already been filled)
!i=1
!do while (nlevs(i)==0 .and. i<=nflds)
!   nlevs(i) = plev
!   i=i+1
!end do
!WRITE(*,*) 'nlevs = ',(nlevs(j),j=1,nflds)

! read CAM 'initial' file fields desired

!2d fields; hard coded for only 1 2d field here (generalize???)
do ifld= 1, 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 3D; lon, lat, TIME(=1)
   call check(nf90_get_var(ncfileid, ncfldid, var%vars_2d(:, :, 1) &
             ,start=(/1,1,1/) ,count=(/num_lons, num_lats, 1/) ))
end do

! 3d fields; hard coded for only 4 here now (generalize???)
do ifld=2, 5
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   PRINT*,'reading ',cflds(ifld),' using id ',ncfldid
!  fields on file are 4D; lon, lev, lat, TIME(=1) 
   call check(nf90_get_var(ncfileid, ncfldid, var%vars_3d(:, :, :, ifld - 1) &
             ,start=(/1,1,1,1/) ,count=(/num_lons, num_levs, num_lats,1/) ))

!!! ?? WARNING: DOES THE NUMBER OF VERTICAL LEVELS PER 3D FIELD VARY
! AS KEVIN"S CODE SUGGESTED. IF SO NEED TO FIX
end do

end subroutine read_cam_init

!#######################################################################
subroutine prog_var_to_vector(var, x)

implicit none

type(model_type), intent(in) :: var
real(r8), intent(out) :: x(:)

integer :: i, j, k, nf, nt, index

! Do order as ps, t, u, v, q, tracers to be consistent with b-grid

! Start copying fields to straight vector
index = 0
do i = 1, num_lons
   do j = 1, num_lats
!  Surface pressure and other 2d flds are first
      do nf = 1, n2dflds
         index = index + 1
         x(index) = var%vars_2d(i, j, nf)
      end do
!     u,v,t,q, and tracers at successively lower levels
      do k = 1, num_levs
         do nf= 1, n3dflds
            index = index + 1
            x(index) = var%vars_3d(i, k, j, nf)
         end do
         do nt = 1, num_tracers
            IF (i==1 .and. j==1) PRINT*,'filling tracers'
            index = index + 1
            x(index) = var%tracers(i, k, j, nt)
         end do
      end do
   end do
end do

! Temporary check
if(index /= model_size) then
   write(*, *) 'prog_var_to_vector bad index sum '
   write(*, *) 'index, model_size ', index, model_size
   stop
endif

end subroutine prog_var_to_vector

!#######################################################################

subroutine vector_to_prog_var(x, var) 

implicit none

real(r8), intent(in) :: x(:)
type(model_type), intent(out) :: var

integer :: i, j, k, nf, nt, n0, index

! Start copying fields from straight vector
index = 0
do i = 1, num_lons
   do j = 1, num_lats
! Surface pressure and other 2d fields are first
      do nf = 1, n2dflds
         index = index + 1
         var%vars_2d(i, j, nf) = x(index)
      end do
!     u,v,t,q  and tracers at successive levels
      do k = 1, num_levs
         do nf = 1, n3dflds
            index = index + 1
            var%vars_3d(i, k, j, nf) = x(index)
         end do 
         do nt = 1, num_tracers
            index = index + 1
            var%tracers(i, k, j, nt) = x(index)
         end do
      end do
   end do
end do

! Temporary check
if(index /= model_size) then
   write(*, *) 'prog_var_to_vector bad index sum '
   write(*, *) 'index, model_size ', index,  model_size
   stop
endif

end subroutine vector_to_prog_var

!#######################################################################
subroutine write_cam_init(file_name, var)

implicit none

character (len = *), intent(in) :: file_name
type(model_type), intent(in) :: var

integer ifld, ncfileid, ncfldid

! Read CAM 'initial' file domain info
call check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid))

! Try doing this in revised order
! 2d fields are first
do ifld = 1, 1
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid, var%vars_2d(:, :, ifld), &
      start=(/1, 1, 1/), count = (/num_lons, num_lats, 1/)))
end do 

! write CAM 'initial' file fields that have been updated
! ????? WARNING: PROBLEMS WITH TRACERS HERE
! ??? THIS WHOLE MODULE NEEDS FURTHER REVISION WITH KEVIN"S HELP
! Now do 3d fields, ignoring tracers for now
do ifld= 2, 5
   call check(nf90_inq_varid(ncfileid, trim(cflds(ifld)), ncfldid))
   call check(nf90_put_var(ncfileid, ncfldid, var%vars_3d(:,:,:,ifld - 1) &
             ,start=(/1,1,1,1/) ,count=(/num_lons, num_levs, num_lats,1/) ))
end do

call check(nf90_close(ncfileid))

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

! Initializes class data for CAM model (all the stuff that needs to
! be done once. For now, does this by reading info from a fixed
! name netcdf file. Need to make this file a namelist parameter
! at some point.

implicit none

integer :: i, j, unit, ierr, io

! Begin by reading the namelist input
if(file_exist('input.nml')) then
   unit = open_file(file = 'input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(unit, nml = model_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'model_nml')
   enddo
 11 continue
   call close_file(unit)
endif

! Temporary namelist validation
write(*, *)'namelist read: values are'
write(*, *)'output_state_vector is ', output_state_vector

write(*,*)'model attributes:'
write(*,*)'   ',source
write(*,*)'   ',revision
write(*,*)'   ',revdate

! Get num lons, lats and levs from netcdf and put in global storage
call read_cam_init_size(model_config_file, num_lons, num_lats, num_levs)

! Compute the model size (Need to be sure where all these come from
! and if they themselves need to be initialized here)

! Need to set Time_step_atmos here to 1 hour for now 
Time_step_atmos = set_time(43200, 0)

! Compute overall model size and put in global storage
model_size = num_lons * num_lats * (n2dflds + num_levs * &
   (n3dflds + pcnst + pnats))

! Allocate space for longitude and latitude global arrays
allocate(lons(num_lons), lats(num_lats))

! Need values for lons and lats, too; should come from netcdf file in read_cam_init_size
do i = 1, num_lons
   lons(i) = 360.0 * (i - 1.0) / num_lons
end do

do i = 1, num_lats
   lats(i) = -90.0 + 180.0 * (i - 0.5) / num_lats
end do

write(*, *) 'CAM size initialized as ', model_size

end subroutine static_init_model

!#######################################################################

subroutine init_model_instance(var)

! Initializes an instance of a cam model state variable

implicit none

type(model_type), intent(out) :: var

! Initialize the storage space and return
allocate(var%vars_2d(num_lons, num_lats, n2dflds), &
   var%vars_3d(num_lons, num_levs, num_lats, n3dflds), &
   var%tracers(num_lons, num_levs, num_lats, num_tracers))

end subroutine init_model_instance

!----------------------------------------------------------------------

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
! Can be used to test the assim capabilities with a null advance

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
col_num = index / num_per_col 
col_elem = index - col_num * num_per_col

! What lon and lat index for this column
lon_index = col_num / num_lats
lat_index = col_num - lon_index * num_lats

! Get actual lon lat values from static_init arrays ???
lon = lons(lon_index + 1)
lat = lats(lat_index + 1)

! Now figure out which beast in column this is
! Surface pressure is the first element
lev = (col_elem - 1) / (n3dflds + pcnst + pnats) + 1
if(col_elem == 0) then
   local_var_type = TYPE_PS
   lev = -1
else
   var_type_temp = mod(col_elem - 1, n3dflds + pcnst + pnats)
   if(var_type_temp == 0) then
      local_var_type = TYPE_T
   else if(var_type_temp == 1) then
      local_var_type = TYPE_U
   else if(var_type_temp == 2) then
      local_var_type = TYPE_V
   else
      local_var_type = TYPE_Q
   endif
endif

!write(*, 11) lon, lat, lev, local_var_type
11 format(1x, 3(f6.2, 1x), i3)
location = set_location(lon, lat, lev)

! If the type is wanted, return it
if(present(var_type)) var_type = local_var_type


end subroutine get_state_meta_data

!#######################################################################

function model_interpolate(x, location, type)

implicit none

real :: model_interpolate
real, intent(in) :: x(:)
type(location_type), intent(in) :: location
integer, intent(in) :: type

integer :: lon_below, lon_above, lat_below, lat_above, i
real :: bot_lon, top_lon, delta_lon, bot_lat, top_lat
real :: lon_fract, lat_fract, val(2, 2), temp_lon, a(2)
real :: lon, lat, level, lon_lat_lev(3)

! First version only allows (level specifies verical???)
lon_lat_lev = get_location(location)
lon = lon_lat_lev(1); lat = lon_lat_lev(2); level = lon_lat_lev(3)

! Get lon and lat grid specs, num_lons, num_lats are globally defined for cam
   bot_lon = lons(1)
   top_lon = lons(num_lons)
   delta_lon = lons(2) - lons(1)
   bot_lat = lats(1)
   top_lat = lats(num_lats)

! Compute bracketing lon indices
if(lon >= bot_lon .and. lon <= top_lon) then
   lon_below = int((lon - bot_lon) / delta_lon) + 1
   lon_above = lon_below + 1
   lon_fract = (lon - ((lon_below - 1) * delta_lon + bot_lon)) / delta_lon
else
! At wraparound point
   lon_below = num_lons
   lon_above = 1
   if(lon < bot_lon) then
      temp_lon = lon + 360.0
   else
      temp_lon = lon
   endif
   lon_fract = (temp_lon - top_lon) / delta_lon
endif


! Next, compute neighboring lat rows
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE MAY BE WRONG
! Inefficient search used for latitudes in Gaussian grid. Might want to speed up.
if(lat >= bot_lat .and. lat <= top_lat) then

   do i = 2, num_lats
      if(lat <= lats(i)) then
         lat_above = i
         lat_below = i - 1
         lat_fract = (lat - lats(lat_below)) / (lats(lat_above) - lats(lat_below))
         goto 20
      end if 
   end do

else if(lat <= bot_lat) then
! South of bottom lat NEED TO DO BETTER: NOT REALLY REGULAR
   lat_below = 1
   lat_above = 2
   lat_fract = 0.0
else
! North of top lat NEED TO DO BETTER: NOT REALLY REGULAR
   lat_below = num_lats - 1
   lat_above = num_lats
   lat_fract = 1.0
endif

! Level is obvious for now

! Now, need to find the values for the four corners
20 continue
val(1, 1) =  get_val(x, lon_below, lat_below, int(level), type)
val(1, 2) =  get_val(x, lon_below, lat_above, int(level), type)
val(2, 1) =  get_val(x, lon_above, lat_below, int(level), type)
val(2, 2) =  get_val(x, lon_above, lat_above, int(level), type)

! Do the weighted average for interpolation
!write(*, *) 'fracts ', lon_fract, lat_fract
do i = 1, 2
   a(i) = lon_fract * val(2, i) + (1.0 - lon_fract) * val(1, i)
end do

model_interpolate = lat_fract * a(2) + (1.0 - lat_fract) * a(1)

end function model_interpolate

!#######################################################################

function get_val(x, lon_index, lat_index, level, type)

implicit none

real :: get_val
real, intent(in) :: x(:)
integer, intent(in) :: lon_index, lat_index, level, type

integer :: per_col, index

! Compute size of grid storage in a column; includes tracers
! Single 2D state vector is pressure
per_col = 1 + num_levs * n3tflds

! Find the starting index for this column
index = per_col * (lat_index - 1 + (lon_index - 1) * num_lats)

! Pressure is first 
if(type == 3) then
   index = index + 1
else
! For interior fields compute the base for their level and add offset
   index = index + 1 + (level - 1) * n3tflds
! Temperature
   if(type == 4) then
      index = index + 1
! U wind component
   else if(type == 1) then
      index = index + 2
! V wind component
   else if(type == 2) then
      index = index + 3
! Tracers
   else if(type > 4) then
      index = index + type - 1
   end if
endif
   
get_val = x(index)

end function get_val

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
integer :: num, max_size, i, j, num1
integer :: hsize, num_per_col, col_base_index
integer, allocatable :: lon_ind(:), lat_ind(:)
real(r8), allocatable :: close_dist(:)

write(*, *) 'in model_get_close_states', radius
loc_array = get_location(o_loc)
write(*, *) 'oloc is ', loc_array(:)

! Number found starts at 0
number = 0

! Assume that grid size is known from static initialized storage

! Num of close horizontal grid points starts at 0, too
num = 0
! For now, just allocate enough space for all grid points, may want
! to make this smaller at some point for big models.
max_size = num_lons * num_lats
allocate(lon_ind(max_size), lat_ind(max_size), close_dist(max_size))

! Look for close grid points on the 
call grid_close_states2(o_loc, lons, lats, num_lons, num_lats, radius, &
   num, lon_ind, lat_ind, close_dist)
write(*, *) 'back from grid_close_states num = ', num

! Compute size of grid storage for full levels
hsize = num_lons * num_lats
num_per_col = num_levs * (n3dflds + pcnst + pnats) + n2dflds

! Add all variables in this column to the close list with this distance
write(*, *) 'available space is size(indices) ', size(indices)
do i = 1, num
   col_base_index = ((lon_ind(i) - 1) * num_lats + lat_ind(i) - 1) * num_per_col
   do j = 1, num_per_col
      number = number + 1
      if(number <= size(indices)) indices(number) = col_base_index + j
      if(number <= size(dist)) dist(number) = close_dist(i)
   end do
end do
write(*, *) 'number at end is ', number

deallocate(lon_ind, lat_ind, close_dist)

end subroutine model_get_close_states

!#######################################################################


subroutine grid_close_states2(o_loc, lons, lats, nlon, nlat, radius, &
   num, close_lon_ind, close_lat_ind, close_dist)

! Finds close state points from a particular grid;

implicit none

type(location_type), intent(in) :: o_loc
integer, intent(in) :: nlon, nlat
real(r8), intent(in) :: lons(nlon), lats(nlat), radius
integer, intent(inout) :: num
integer, intent(inout) :: close_lon_ind(:), close_lat_ind(:)
real(r8), intent(out) :: close_dist(:)

real(r8) :: glat, glon, loc_array(3), o_lon, o_lat, o_lev
real(r8) :: gdist, diff, row_dist(nlon)
integer :: blat_ind, blon_ind, i, j, lat_ind, lon_ind
integer :: row_lon_ind(nlon), row_num
real(r8), parameter :: glev = 1.0
type(location_type) :: loc

! Get the lat and lon from the loc
loc_array = get_location(o_loc)
o_lon = loc_array(1)
o_lat = loc_array(2)

! Get index to closest lat and lon for this observation
blat_ind = get_closest_lat_index(o_lat, lats, nlat)
!write(*, *) 'closest latitude in grid is ', blat_ind, lats(blat_ind)
blon_ind = get_closest_lon_index(o_lon, lons, nlon)
!write(*, *) 'closest longitude in grid is ', blon_ind, lons(blon_ind)

! Begin a search along the latitude axis in the positive direction
do lat_ind = blat_ind, nlat
   glat = lats(lat_ind)
! Take care of storage round-off
   if(glat < -90.0) glat = 0.0
   if(glat > 90.0) glat = 90.0

! Search all the contiguous close longitudes around the base longitude
   call lon_search(glat, glev, blon_ind, o_loc, radius, lons, &
      row_lon_ind, row_dist, row_num)
! If none are found, it's time to search in the negative latitude direction
   if(row_num == 0) goto 11
! Copy the points found in the row into summary storage
   close_lon_ind(num+1 : num+row_num) = row_lon_ind(1:row_num)
   close_lat_ind(num+1 : num+row_num) = lat_ind
   close_dist(num+1 : num+row_num) = row_dist(1:row_num)
   num = num + row_num
end do

! Search in the negative lat direction
11 continue
do lat_ind = blat_ind - 1, 1, -1
   glat = lats(lat_ind)
! Take care of storage round-off
   if(glat < -90.0) glat = 0.0
   if(glat > 90.0) glat = 90.0

! Search all the contiguous close longitudes around the base longitude
   call lon_search(glat, glev, blon_ind, o_loc, radius, lons, &
      row_lon_ind, row_dist, row_num)
! If none are found, it's time to give up
   if(row_num == 0) return
! Copy the points found in the row into summary storage
   close_lon_ind(num+1 : num+row_num) = row_lon_ind(1:row_num)
   close_lat_ind(num+1 : num+row_num) = lat_ind
   close_dist(num+1 : num+row_num) = row_dist(1:row_num)
   num = num + row_num
end do

end subroutine grid_close_states2

!------------------------------------------------------------------------

subroutine lon_search(glat, glev, blon_ind, o_loc, radius, lons, &
   close_lon_ind, close_dist, num)

! Given an observation location and radius and a latitude row from a grid,
! searches to find all longitude points in this row that are within radius
! of the observation location and returns their latitude index, longitude
! index, and the distance between them and the observation.

real(r8), intent(in) :: glat, glev, radius, lons(:)
integer, intent(in) :: blon_ind
type(location_type), intent(in) :: o_loc
integer, intent(out) :: close_lon_ind(:), num
real(r8), intent(out) :: close_dist(:)

integer :: nlon, j, max_pos, lon_ind
type(location_type) :: loc
real(r8) :: glon, gdist

! Total number found is 0 at start
num = 0
nlon = size(lons)

! Search as far as possible in the positive direction
do j = 0, nlon - 1
   max_pos = j
   lon_ind = blon_ind + j
   if(lon_ind > nlon) lon_ind = lon_ind - nlon
   glon = lons(lon_ind)
! Correct for longitude storage round-off
   if(glon > 360.0) glon = 360.0
   if(glon < 0.0) glon = 0.0
   loc = set_location(glon, glat, glev)
   gdist = get_dist(loc, o_loc)
   if(gdist <= radius) then
      num = num + 1
      close_lon_ind(num) = lon_ind
      close_dist(num) = gdist
! If radius is too far for closest longitude, no need to search further or to search other side
   else if (j == 0) then
      return
   else
! Look in negative longitude offset direction next
      goto 21
   endif
end do
! Falling off end means the whole longitude circle has been searched; move along
return

! Search around the other way
21 continue
do j = 1, nlon - 1 - max_pos
   lon_ind = blon_ind - j
   if(lon_ind < 1) lon_ind = nlon + lon_ind
   glon = lons(lon_ind)
! Correct for longitude storage round-off
   if(glon > 360.0) glon = 360.0
   if(glon < 0.0) glon = 0.0
   loc = set_location(glon, glat, glev)
   gdist = get_dist(loc, o_loc)
   if(gdist <= radius) then
      num = num + 1
      close_lon_ind(num) = lon_ind
      close_dist(num) = gdist
   else
! No more longitudes in negative direction
      return
   endif
end do

end subroutine lon_search
                                              
!#######################################################################

function nc_write_model_atts( ncFileID ) result (ierr)

use typeSizes
use netcdf
implicit none

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!-----------------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: IDimID, JDimID, levDimID, MemberDimID, StateVarDimID, TimeDimID
integer :: IVarID, JVarID, levVarID, StateVarID, StateVarVarID
integer :: i

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
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                        len=model_size, dimid = StateVarDimID))

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source",source))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate))

! how about namelist input? might be nice to save ...

!-------------------------------------------------------------------------------
! Define the new dimensions IDs
!-------------------------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="I",   len = num_lons,   dimid =   IDimID))
call check(nf90_def_dim(ncid=ncFileID, name="J",   len = num_lats,   dimid =   JDimID))
call check(nf90_def_dim(ncid=ncFileID, name="lev", len = num_levs,   dimid = levDimID))

!-------------------------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!-------------------------------------------------------------------------------

! Temperature Grid Longitudes
call check(nf90_def_var(ncFileID, name="I", xtype=nf90_double, &
                                               dimids=IDimID, varid=IVarID) )
call check(nf90_put_att(ncFileID, IVarID, "long_name", "longitude"))
call check(nf90_put_att(ncFileID, IVarID, "units", "degrees_east"))
call check(nf90_put_att(ncFileID, IVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)))

! Temperature Grid Latitudes
call check(nf90_def_var(ncFileID, name="J", xtype=nf90_double, &
                                               dimids=JDimID, varid=JVarID) )
call check(nf90_put_att(ncFileID, JVarID, "long_name", "latitude"))
call check(nf90_put_att(ncFileID, JVarID, "units", "degrees_north"))
call check(nf90_put_att(ncFileID, JVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)))

! (Common) grid levels
call check(nf90_def_var(ncFileID, name="level", xtype=nf90_int, &
                                                dimids=levDimID, varid=levVarID) )
call check(nf90_put_att(ncFileID, levVarID, "long_name", "level"))

! Number of Tracers
! call check(nf90_def_var(ncFileID, name="tracer", xtype=nf90_int, &
!                                                  dimids=tracerDimID, varid=tracerVarID) )
! call check(nf90_put_att(ncFileID, tracerVarID, "long_name", "tracer identifier"))
! write (*,*)'tracerVarID determined'

if ( output_state_vector ) then
   !----------------------------------------------------------------------------
   ! Create attributes for the state vector
   !----------------------------------------------------------------------------
                                                                                               
  ! Define the state vector coordinate variable
   call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID))
   call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"))
   call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical") )
   call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)))
                                                                                               
   ! Define the actual state vector
   call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), varid=StateVarID))
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"))
   call check(nf90_put_att(ncFileID, StateVarId, "vector_to_prog_var","CAM"))

   ! Leave define mode so we can fill                                                          
   call check(nf90_enddef(ncfileID))                                                           
                                                                                               
   ! Fill the state variable coordinate variable                                               
   call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ))                
                                                                                               
else

   write(*,*)'ERROR:CAM:model_mod:    trying to output the prognostic variables.'
   write(*,*)'      That is not implemented yet.'
   write(*,*)'      TJH 27 June 2003'
   stop

endif

!-------------------------------------------------------------------------------
! Fill the variables
!-------------------------------------------------------------------------------

call check(nf90_put_var(ncFileID,      IVarID, lons ))
call check(nf90_put_var(ncFileID,      JVarID, lats ))
call check(nf90_put_var(ncFileID,    levVarID, (/ (i,i=1,   num_levs) /) ))
! call check(nf90_put_var(ncFileID, tracerVarID, (/ (i,i=1,ntracer) /) ))

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
! type(prog_var_type) :: Var

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, psVarID, tVarID, rVarID, uVarID, vVarID
integer :: tis, tie, tjs, tje       ! temperature grid start/stop
integer :: vis, vie, vjs, vje       ! velocity    grid start/stop
integer :: kub, klb
integer :: nTmpI, nTmpJ, nVelI, nVelJ, nlev, ntracer, i

ierr = 0     ! assume normal termination


!-------------------------------------------------------------------------------
! Get the bounds for storage on Temp and Velocity grids
! Everything in here is from the bgrid ... don't have a clue for CAM ... yet.
!-------------------------------------------------------------------------------

! tis = Dynam%Hgrid%Tmp%is; tie = Dynam%Hgrid%Tmp%ie
! tjs = Dynam%Hgrid%Tmp%js; tje = Dynam%Hgrid%Tmp%je
! vis = Dynam%Hgrid%Vel%is; vie = Dynam%Hgrid%Vel%ie
! vjs = Dynam%Hgrid%Vel%js; vje = Dynam%Hgrid%Vel%je
! kub = Var_dt%kub
! klb = Var_dt%klb

nTmpI   = tie - tis + 1
nTmpJ   = tje - tjs + 1
! nlev    = Var_dt%kub - Var_dt%klb + 1
! ntracer = Var_dt%ntrace 
nVelI   = vie - vis + 1
nVelJ   = vje - vjs + 1

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! then get all the Variable ID's we need.
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

if ( output_state_vector ) then

   call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)))                               

else

   write(*,*)'FATAL ERROR: CAM:nc_write_model_vars  is not operational.'
   write(*,*)'   you should not be here ...'
   write(*,*)'model_nml: output_state_vector  MUST be .true.    for now.'
   write(*,*)'TJH ... 27 June 2003'
   stop
   
   !----------------------------------------------------------------------------
   ! Fill the variables
   !----------------------------------------------------------------------------

 ! x = statevec ! Unfortunately, have to explicity cast it ...
                ! the filter uses a type=double,
                ! the vector_to_prog_var function expects a single.
 !  call vector_to_prog_var(x, get_model_size(), Var)
   
   
 ! call check(NF90_inq_varid(ncFileID, "ps", psVarID))
 ! call check(nf90_put_var( ncFileID, psVarID, Var%ps(tis:tie, tjs:tje), &
 !                          start=(/ 1, 1, copyindex, timeindex /) ))

 ! call check(NF90_inq_varid(ncFileID,  "t",  tVarID))
 ! call check(nf90_put_var( ncFileID,  tVarID, Var%t( tis:tie, tjs:tje, klb:kub ), &
 !                          start=(/ 1, 1, 1, copyindex, timeindex /) ))

 ! call check(NF90_inq_varid(ncFileID,  "u",  uVarID))
 ! call check(nf90_put_var( ncFileID,  uVarId, Var%u( vis:vie, vjs:vje, klb:kub ), &
 !                          start=(/ 1, 1, 1, copyindex, timeindex /) ))

 ! call check(NF90_inq_varid(ncFileID,  "v",  vVarID))
 ! call check(nf90_put_var( ncFileID,  vVarId, Var%v( vis:vie, vjs:vje, klb:kub ), &
 !                          start=(/ 1, 1, 1, copyindex, timeindex /) ))

 ! if ( ntracer > 0 ) then
 !    call check(NF90_inq_varid(ncFileID,  "r",  rVarID))
 !    call check(nf90_put_var( ncFileID,  rVarID, &
 !                  Var%r( tis:tie, tjs:tje, klb:kub, 1:ntracer ), & 
 !                 start=(/   1,       1,       1,     1,    copyindex, timeindex /) ))
 ! endif
endif

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
      ierr = istatus
      stop
    end if
  end subroutine check

end function nc_write_model_vars



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



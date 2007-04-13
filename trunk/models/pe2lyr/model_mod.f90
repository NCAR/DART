! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!-----------------------------------------------------------------------
! Assimilation interface for 2 layer PE model
!-----------------------------------------------------------------------

!---------------- m o d u l e   i n f o r m a t i o n ------------------

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time
use    utilities_mod, only : file_exist, open_file, close_file, &
                             register_module, error_handler, E_ERR, E_MSG
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use     location_mod, only : location_type, get_location, set_location, get_dist, &
                             LocationDims, LocationName, LocationLName, &
                             get_close_maxdist_init, get_close_obs_init, get_close_obs


use  pe2lyr_mod, only : modelvars, modelvars_init, rk4
use spharmt_mod, only : sphere, getvrtdiv, spharm, spharmt_init, getuv

!-----------------------------------------------------------------------

implicit none
private

include "resolt31.h"

public :: get_model_size, &
          adv_1step, &
          get_state_meta_data, &
          model_interpolate, &
          get_model_time_step, &
          end_model, &
          static_init_model,  &
          init_time, &
          init_conditions, &
          TYPE_u, TYPE_v, TYPE_z, &
          nc_write_model_atts, &
          nc_write_model_vars, &
          pert_model_state, &
          get_close_maxdist_init, get_close_obs_init, get_close_obs, ens_mean_for_model


! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!-----------------------------------------------------------------------
! Public definition of variable types
integer, parameter :: TYPE_u = 1, TYPE_v = 2, TYPE_z = 901

type(time_type) :: time_step
type (modelvars):: model_dat
type (sphere) :: sphere_dat
!real(r8), allocatable :: lons(:), lats(:)

integer :: levs(2)
real :: lons(nlons), lats(nlats)

contains

!#######################################################################


!#######################################################################

subroutine adv_1step(x, Time)

!! Does single time-step advance for Pe2lyr model with vector state as
!! input and output. call rh4 in pe2lyr.mod to advance model

implicit none

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: Time
integer :: i,j,mm,k
real, dimension(nlons,nlats,2) :: ug,vg,zg,delpig


!! update model_dat from X

mm =1
do i = 1,nlats
  do j=1, nlons
    ug(j,i,1)=x(mm) 
    ug(j,i,2)=x(nlats*nlons+mm) 
    vg(j,i,1)=x(2*nlats*nlons+mm) 
    vg(j,i,2)  =x(3*nlats*nlons+mm) 
    zg(j,i,1) =x(4*nlats*nlons+mm) 
    zg(j,i,2) = x(5*nlats*nlons+mm) 
    mm = mm+1
  enddo
enddo 

!! convert zg to delpig (pseudo-height to exner function)

delpig(:,:,2) =cp- zg(:,:,2)*g/theta1 -model_dat%pitop
delpig(:,:,1) =cp- zg(:,:,1)*g/theta1 -model_dat%pitop - delpig(:,:,2)

!! get spectral data

do k =1,2
 call getvrtdiv(sphere_dat,model_dat%vrtnm(:,k),&
    model_dat%divnm(:,k),ug(:,:,k),vg(:,:,k))
 call spharm(sphere_dat,delpig(:,:,k),&
    model_dat%delpinm(:,k),1)
enddo

!! forward model_dat

call rk4(model_dat,sphere_dat)

!! get grid data

do k=1,2
  call getuv(sphere_dat,model_dat%vrtnm(:,k),model_dat%divnm(:,k),&
             ug(:,:,k),vg(:,:,k))
  call spharm(sphere_dat,delpig(:,:,k),model_dat%delpinm(:,k),-1)
enddo   

zg(:,:,2)=  (theta1/g)*(cp - (delpig(:,:,2)+model_dat%pitop))
zg(:,:,1)=  (theta1/g)*(cp - (delpig(:,:,1)+model_dat%pitop+delpig(:,:,2)))

!! updata x from model grid data

mm =1
do i = 1,nlats
  do j=1, nlons
    x(mm) = ug(j,i,1)
    x(nlats*nlons+mm) = ug(j,i,2)
    x(2*nlats*nlons+mm) =  vg(j,i,1)
    x(3*nlats*nlons+mm) =  vg(j,i,2)  
    x(4*nlats*nlons+mm) = zg(j,i,1)
    x(5*nlats*nlons+mm) = zg(j,i,2)
    mm = mm+1
  enddo
enddo 


!!===change x vector to model_dat


end subroutine adv_1step

!#######################################################################

subroutine static_init_model()

! INitializes class data for a pe2lyr model (all the stuff that needs to
! be done once.

implicit none
integer :: i
real :: pi

time_step =  set_time(int(dt),0)

pi = 4.*atan(1.0)

call spharmt_init(sphere_dat,nlons,nlats,ntrunc,rsphere)
call modelvars_init(model_dat,sphere_dat,&
   cp,g,r,p0,efold,ndiss,dt,delta_pi,ztop,rot,tdrag,tdiab,theta1,delth,hmax)

! evenly spaced lons (starting at Greenwich, going east, no wraparound).
do i = 1, nlons
   lons(i) = 360.0 * (i - 1.0) / nlons
end do


! gaussian lats (north to south).
do i = 1, nlats
   lats(i) = (180./pi)*asin(sphere_dat%gaulats(i))
end do


levs(1) =1
levs(2) =2

end subroutine static_init_model


!#######################################################################

subroutine init_conditions(x)

implicit none

real(r8), intent(inout) :: x(:)
integer :: k,i,j,mm
real, dimension(nlons,nlats,2) :: ug,vg,delpig,zg

! convert model_dat object to state vector.
!xue
!! initialize spectral transform structure (derived data type).
!call spharmt_init(sphere_dat,nlons,nlats,ntrunc,rsphere)
! initialize model structure (derived data type).
! parameters from include file.
!call modelvars_init(model_dat,sphere_dat,&
!  cp,g,r,p0,efold,ndiss,dt,delta_pi,ztop,rot,tdrag,tdiab,theta1,delth,hmax)

!xue

do k=1,2
  call getuv(sphere_dat,model_dat%vrtnm(:,k),model_dat%divnm(:,k),&
             ug(:,:,k),vg(:,:,k))
  call spharm(sphere_dat,delpig(:,:,k),model_dat%delpinm(:,k),-1)
enddo   

zg(:,:,2)=  (theta1/g)*(cp - (delpig(:,:,2)+model_dat%pitop))
zg(:,:,1)=  (theta1/g)*(cp - (delpig(:,:,2)+model_dat%pitop+delpig(:,:,1)))

mm =1
do i = 1,nlats
  do j=1, nlons
    x(mm) = ug(j,i,1)
    x(nlats*nlons+mm) = ug(j,i,2)
    x(2*nlats*nlons+mm) =  vg(j,i,1)
    x(3*nlats*nlons+mm) =  vg(j,i,2)  
    x(4*nlats*nlons+mm) = zg(j,i,1)
    x(5*nlats*nlons+mm) = zg(j,i,2)
    mm = mm+1
 enddo
enddo 

end subroutine init_conditions

!#######################################################################

function get_model_size()

implicit none

integer :: get_model_size

get_model_size = model_size

end function get_model_size

!#######################################################################

function get_model_time_step()
!------------------------------------------------------------------------

type(time_type) :: get_model_time_step

! model_dat%dt

get_model_time_step =  time_step

end function get_model_time_step

!#######################################################################

subroutine get_state_meta_data(indx, location, var_type)
!---------------------------------------------------------------------

implicit none

integer, intent(in) :: indx
type(location_type), intent(out) :: location
integer, intent(out), optional :: var_type

!integer :: indx

integer :: lat_index, lon_index,lev_index
integer :: u_indxmax,v_indxmax,z_indxmax
real(r8) :: lon, lat, lev

!print *,'get_state_meta_data'
u_indxmax = 2*nlats*nlons
v_indxmax = 4*nlats*nlons
z_indxmax = 6*nlats*nlons

!! interface height obs type = 901 
!! ( a big number so no conflict in future ))
if(indx <= u_indxmax) then
   lev_index =  (indx-1)/(nlats*nlons) +1
   lat_index = ((indx-1) - ((lev_index-1)*nlats*nlons)) / nlons +1
   lon_index =  (indx-1) - ((lev_index-1)*nlats*nlons) - ((lat_index-1)*nlons) +1
   if(present(var_type)) var_type = TYPE_u
else if(indx <= v_indxmax) then
   lev_index =  (indx-u_indxmax-1)/(nlats*nlons) +1
   lat_index = ((indx-u_indxmax-1) - ((lev_index-1)*nlats*nlons)) / nlons +1
   lon_index =  (indx-u_indxmax-1) - ((lev_index-1)*nlats*nlons) - ((lat_index-1)*nlons) +1
   if(present(var_type)) var_type = TYPE_v
else
   lev_index =  (indx-v_indxmax-1)/(nlats*nlons) +1
   lat_index = ((indx-v_indxmax-1) - ((lev_index-1)*nlats*nlons)) / nlons +1
   lon_index =  (indx-v_indxmax-1) - ((lev_index-1)*nlats*nlons) - ((lat_index-1)*nlons) +1
   if(present(var_type)) var_type = TYPE_z
endif 
 
! With the threed_sphere location module ... you specify that the 
! vertical coordinate is a 'level' by 'which_vert' == 1

lon = lons(lon_index)
lat = lats(lat_index)
lev = levs(lev_index)

location = set_location(lon, lat, lev, 1)

end subroutine get_state_meta_data

!#######################################################################

subroutine model_interpolate(x, location, type, obs_val, istatus)

implicit none


real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: type
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

integer :: lon_below, lon_above, lat_below, lat_above, i
real :: bot_lon, top_lon, delta_lon, bot_lat, top_lat
real :: lon_fract, lat_fract, val(2, 2), temp_lon, a(2)
real :: lon, lat, level, lon_lat_lev(3)

! All interpolations okay for now
istatus = 0

lon_lat_lev = get_location(location)
lon = lon_lat_lev(1); lat = lon_lat_lev(2); level = lon_lat_lev(3)

! Get lon and lat grid specs are globally defined for pe2lyr
   bot_lon = lons(1)
   top_lon = lons(nlons)
   delta_lon = lons(2) - lons(1)
   bot_lat = lats(1)
   top_lat = lats(nlats)

! Compute bracketing lon indices
if(lon >= bot_lon .and. lon <= top_lon) then
   lon_below = int((lon - bot_lon) / delta_lon) + 1
   lon_above = lon_below + 1
   lon_fract = (lon - ((lon_below - 1) * delta_lon + bot_lon)) / delta_lon
else
! At wraparound point
   lon_below = nlons
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
if(lat <= bot_lat .and. lat >= top_lat) then

   do i = 2, nlats
      if(lat >= lats(i)) then
         lat_above = i
         lat_below = i - 1
         lat_fract = (lat - lats(lat_below)) / (lats(lat_above) - lats(lat_below))
         goto 20
      end if 
   end do

else if(lat >= bot_lat) then
! North of bottom lat NEED TO DO BETTER: NOT REALLY REGULAR
   lat_below = 1
   lat_above = 2
   lat_fract = 1.0
else
! South of top lat NEED TO DO BETTER: NOT REALLY REGULAR
   lat_below = nlats - 1
   lat_above = nlats
   lat_fract = 0.0
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

obs_val = lat_fract * a(2) + (1.0 - lat_fract) * a(1)


end subroutine model_interpolate

!#######################################################################

function get_val(x, lon_index, lat_index, level, type)

implicit none

real :: get_val
real(r8), intent(in) :: x(:)
integer, intent(in) :: lon_index, lat_index, level, type

integer :: indx


! order is u,v,z 
if(type == TYPE_u) then
   indx = (level-1)*nlats*nlons+(lat_index-1)*nlons + lat_index
else if(type == TYPE_v) then
   indx = 2*nlats*nlons+(level-1)*nlats*nlons+(lat_index-1)*nlons + lon_index
else if(type == TYPE_z) then
   indx = 4*nlats*nlons+(level-1)*nlats*nlons+(lat_index-1)*nlons + lon_index
endif
   
get_val = x(indx)

end function get_val
!########################################################
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
print *,'init_time'
i_time = set_time(0, 0)

end subroutine init_time


!#############################################################
function nc_write_model_atts( ncFileID ) result (ierr)

use typeSizes
use netcdf
implicit none
integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!--------------------------------------------------------------------
! netCDF variables
!--------------------------------------------------------------------

integer :: latDimID, lonDimID, levDimID,  MemberDimID
integer :: levVarID
integer :: StateVarDimID, TimeDimID
integer :: latVarID, lonVarID, uVarID, vVarID,zVarID
integer :: i

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------
character(len=129) :: errstring
character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

ierr = 0                      ! assume normal termination
!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file 
!--------------------------------------------------------------------
call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_Redef(ncFileID))

!--------------------------------------------------------------------
! Determine ID's from stuff already in the netCDF file
!--------------------------------------------------------------------

! make sure time is unlimited dimid

call check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID))
call check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID))

if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)'Time Dimension ID ',TimeDimID,' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', errstring, source, revision, revdate)
endif
!--------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!--------------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=model_size, dimid = StateVarDimID))

!--------------------------------------------------------------------
! Write Global Attributes 
!--------------------------------------------------------------------
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source",source))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","pe2lyr"))

!--------------------------------------------------------------------
! Define the new dimensions IDs
!--------------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="lat",   len = nlats,   dimid = latDimID)) 
call check(nf90_def_dim(ncid=ncFileID, name="lon",   len = nlons,   dimid = lonDimID)) 
call check(nf90_def_dim(ncid=ncFileID, name="lev",   len = nlevs,   dimid = levDimID)) 

!--------------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!--------------------------------------------------------------------

! Temperature Grid Longitudes
call check(nf90_def_var(ncFileID, name="lon", xtype=nf90_double, &
                                               dimids=lonDimID, varid=lonVarID) )
call check(nf90_put_att(ncFileID, lonVarID, "long_name", "longitude"))
call check(nf90_put_att(ncFileID, lonVarID, "cartesian_axis", "X"))
call check(nf90_put_att(ncFileID, lonVarID, "units", "degrees_east"))
call check(nf90_put_att(ncFileID, lonVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)))

call check(nf90_def_var(ncFileID, name="lat", xtype=nf90_double, &
                                               dimids=latDimID, varid=latVarID) )
call check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude"))
call check(nf90_put_att(ncFileID, latVarID, "cartesian_axis", "Y"))
call check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north"))
call check(nf90_put_att(ncFileID, latVarID, "valid_range", (/  -90.0_r8, 90.0_r8 /)))

call check(nf90_def_var(ncFileID, name="level", xtype=nf90_int, &
                                                dimids=levDimID, varid=levVarID) )
call check(nf90_put_att(ncFileID, levVarID, "long_name", "level"))
call check(nf90_put_att(ncFileID, levVarID, "cartesian_axis", "Z"))
call check(nf90_put_att(ncFileID, levVarID, "units", " "))
call check(nf90_put_att(ncFileID, levVarID, "positive", "down"))


call check(nf90_def_var(ncid=ncFileID, name="u", xtype=nf90_real, &
         dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
         varid  = uVarID))
call check(nf90_put_att(ncFileID, uVarID, "long_name", "zonal wind component"))
call check(nf90_put_att(ncFileID, uVarID, "units", "m/s"))

call check(nf90_def_var(ncid=ncFileID, name="v", xtype=nf90_real, &
         dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
         varid  = vVarID))
call check(nf90_put_att(ncFileID, vVarID, "long_name", "meridional wind component"))
call check(nf90_put_att(ncFileID, vVarID, "units", "m/s"))

call check(nf90_def_var(ncid=ncFileID, name="z", xtype=nf90_real, &
         dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
         varid  = zVarID))
call check(nf90_put_att(ncFileID, vVarID, "long_name", "surface/interface height"))
call check(nf90_put_att(ncFileID, vVarID, "units", "meters"))

call check(nf90_enddef(ncfileID))

!--------------------------------------------------------------------
! Fill the variables
!--------------------------------------------------------------------

call check(nf90_put_var(ncFileID,   lonVarID, lons(:) ))
call check(nf90_put_var(ncFileID,   latVarID, lats(:) ))
call check(nf90_put_var(ncFileID,    levVarID, (/ (i,i=1,   nlevs) /) ))

!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------
call check(nf90_sync(ncFileID))
write(errstring,*)'netCDF file ',ncFileID,' is synched ...'
call error_handler(E_MSG,'nc_write_model_atts',errstring,source,revision, revdate)

contains
  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_atts', &
         trim(nf90_strerror(istatus)), source, revision, revdate)
  end subroutine check
end function nc_write_model_atts

!##########################################################

function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
use typeSizes
use netcdf
implicit none

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          


integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: uVarID,vVarID,zVarID

integer :: mm,i,j
real, dimension(nlons,nlats,2) :: ug,vg,zg
!integer :: i
real, dimension(SIZE(statevec)) :: x

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

x = statevec

mm =1
do i = 1,nlats
  do j=1, nlons
    ug(j,i,1)=x(mm) 
    ug(j,i,2)=x(nlats*nlons+mm) 
    vg(j,i,1)=x(2*nlats*nlons+mm) 
    vg(j,i,2)  =x(3*nlats*nlons+mm) 
    zg(j,i,1) =x(4*nlats*nlons+mm) 
    zg(j,i,2) = x(5*nlats*nlons+mm) 
    mm = mm+1
  enddo
enddo 


call check(NF90_inq_varid(ncFileID,  "u",  uVarID))
call check(nf90_put_var( ncFileID,  uVarId, ug( :, :, : ), &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))

call check(NF90_inq_varid(ncFileID,  "v",  vVarID))
call check(nf90_put_var( ncFileID,  vVarId, vg( :, :, : ), &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))

call check(NF90_inq_varid(ncFileID,  "z",  zVarID))
call check(nf90_put_var( ncFileID,  zVarId, zg( :, :, : ), &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))

!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------

call check(nf90_sync(ncFileID))

contains
  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_vars', &
            trim(nf90_strerror(istatus)), source, revision,revdate)
  end subroutine check

end function nc_write_model_vars




subroutine pert_model_state(state, pert_state, interf_provided)

! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_state




subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

end subroutine ens_mean_for_model



!#######################################################################
end module model_mod

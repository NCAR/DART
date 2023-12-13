! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

!-----------------------------------------------------------------------
! Assimilation interface for 2 layer PE model
!-----------------------------------------------------------------------

!---------------- m o d u l e   i n f o r m a t i o n ------------------

use        types_mod, only : r8, MISSING_R8
use time_manager_mod, only : time_type, set_time
use    utilities_mod, only : file_exist, open_file, close_file, &
                             register_module, error_handler, E_ERR, E_MSG
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use     location_mod, only : location_type, get_location, set_location, get_dist, &
                             LocationDims, LocationName, LocationLName, &
                             get_close_maxdist_init, get_close_obs_init, get_close_obs
use     obs_kind_mod, only : QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, &
                             QTY_GEOPOTENTIAL_HEIGHT

use    netcdf_utilities_mod, only : nc_check
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
          nc_write_model_atts, &
          nc_write_model_vars, &
          pert_model_state, &
          get_close_maxdist_init, get_close_obs_init, get_close_obs, ens_mean_for_model


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------

type(time_type) :: time_step
type(modelvars) :: model_dat
type(sphere)    :: sphere_dat
character(len=129) :: errstring

integer  :: levs(2)
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
    ug(j,i,1)=x(              mm) 
    ug(j,i,2)=x(  nlats*nlons+mm) 
    vg(j,i,1)=x(2*nlats*nlons+mm) 
    vg(j,i,2)=x(3*nlats*nlons+mm) 
    zg(j,i,1)=x(4*nlats*nlons+mm) 
    zg(j,i,2)=x(5*nlats*nlons+mm) 
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
    x(              mm) = ug(j,i,1)
    x(  nlats*nlons+mm) = ug(j,i,2)
    x(2*nlats*nlons+mm) = vg(j,i,1)
    x(3*nlats*nlons+mm) = vg(j,i,2)  
    x(4*nlats*nlons+mm) = zg(j,i,1)
    x(5*nlats*nlons+mm) = zg(j,i,2)
    mm = mm+1
  enddo
enddo 


!!===change x vector to model_dat


end subroutine adv_1step

!#######################################################################

subroutine static_init_model()

! Initializes class data for a pe2lyr model (all the stuff that needs to
! be done once.

implicit none
integer  :: i
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


levs(1) = 1
levs(2) = 2

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
    x(              mm) = ug(j,i,1)
    x(  nlats*nlons+mm) = ug(j,i,2)
    x(2*nlats*nlons+mm) = vg(j,i,1)
    x(3*nlats*nlons+mm) = vg(j,i,2)  
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

integer  :: lat_index, lon_index,lev_index
integer  :: u_indxmax,v_indxmax,z_indxmax
real(r8) :: lon, lat, lev

!print *,'get_state_meta_data'
u_indxmax = 2*nlats*nlons
v_indxmax = 4*nlats*nlons
z_indxmax = 6*nlats*nlons

! avoid out-of-range queries
if (indx > z_indxmax) then
   write(errstring,*)'indx ',indx,' must be between 1 and ', z_indxmax
   call error_handler(E_ERR,'model_mod:get_state_meta_data', errstring, source, revision, revdate)
endif

if(indx <= u_indxmax) then
   lev_index =  (indx-1)/(nlats*nlons) +1
   lat_index = ((indx-1) - ((lev_index-1)*nlats*nlons)) / nlons +1
   lon_index =  (indx-1) - ((lev_index-1)*nlats*nlons) - ((lat_index-1)*nlons) +1
   if(present(var_type)) var_type = QTY_U_WIND_COMPONENT
else if(indx <= v_indxmax) then
   lev_index =  (indx-u_indxmax-1)/(nlats*nlons) +1
   lat_index = ((indx-u_indxmax-1) - ((lev_index-1)*nlats*nlons)) / nlons +1
   lon_index =  (indx-u_indxmax-1) - ((lev_index-1)*nlats*nlons) - ((lat_index-1)*nlons) +1
   if(present(var_type)) var_type = QTY_V_WIND_COMPONENT
else
   lev_index =  (indx-v_indxmax-1)/(nlats*nlons) +1
   lat_index = ((indx-v_indxmax-1) - ((lev_index-1)*nlats*nlons)) / nlons +1
   lon_index =  (indx-v_indxmax-1) - ((lev_index-1)*nlats*nlons) - ((lat_index-1)*nlons) +1
   if(present(var_type)) var_type = QTY_GEOPOTENTIAL_HEIGHT
endif 
 
! With the threed_sphere location module ... you specify that the 
! vertical coordinate is a 'level' by 'which_vert' == 1

lon = lons(lon_index)
lat = lats(lat_index)
lev = levs(lev_index)

location = set_location(lon, lat, lev, 1)

end subroutine get_state_meta_data

!#######################################################################

subroutine model_interpolate(x, location, mytype, obs_val, istatus)

implicit none


real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: mytype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

integer  :: lon_below, lon_above, lat_below, lat_above, i
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
val(1, 1) =  get_val(x, lon_below, lat_below, int(level), mytype)
val(1, 2) =  get_val(x, lon_below, lat_above, int(level), mytype)
val(2, 1) =  get_val(x, lon_above, lat_below, int(level), mytype)
val(2, 2) =  get_val(x, lon_above, lat_above, int(level), mytype)

! Do the weighted average for interpolation
!write(*, *) 'fracts ', lon_fract, lat_fract
do i = 1, 2
   a(i) = lon_fract * val(2, i) + (1.0 - lon_fract) * val(1, i)
end do

obs_val = lat_fract * a(2) + (1.0 - lat_fract) * a(1)


end subroutine model_interpolate

!#######################################################################

function get_val(x, lon_index, lat_index, level, mytype)

implicit none

real :: get_val
real(r8), intent(in) :: x(:)
integer,  intent(in) :: lon_index, lat_index, level, mytype

integer :: indx

! the set_location() call does simple error checks for out of range vals
! for lat and lon, but cannot know what valid levels are.  make sure the
! index numbers computed here are within range.
if (level < 1 .or. level > 2) then
   write(errstring,*)'level ',level,' must be 1 or 2'
   call error_handler(E_ERR,'model_mod:get_val', errstring, source, revision, revdate)
endif

! order is u,v,z 
if(mytype == QTY_U_WIND_COMPONENT) then
   indx =               (level-1)*nlats*nlons+(lat_index-1)*nlons + lon_index
else if(mytype == QTY_V_WIND_COMPONENT) then
   indx = 2*nlats*nlons+(level-1)*nlats*nlons+(lat_index-1)*nlons + lon_index
else if(mytype == QTY_GEOPOTENTIAL_HEIGHT) then
   indx = 4*nlats*nlons+(level-1)*nlats*nlons+(lat_index-1)*nlons + lon_index
endif
   
!! should not be possible now; but this error check can be commented back in.
!! (it is out for performance reasons, but if you get any strange values, this
!! is a good first check to reenable.)
!if (indx < 1 .or. indx > size(x)) then
!   write(errstring,*)'index ',indx,' not between 1 and ', size(x), ' (should not be possible)'
!   call error_handler(E_ERR,'model_mod:get_val', errstring, source, revision, revdate)
!endif

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
character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1
character(len=128) :: filename

ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file 
!--------------------------------------------------------------------
call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
         "nc_write_model_atts",'Inquire, '//trim(filename))
call nc_check(nf90_Redef(ncFileID),"nc_write_model_atts",'Redef, '//trim(filename))

!--------------------------------------------------------------------
! Determine ID's from stuff already in the netCDF file
!--------------------------------------------------------------------

! make sure time is unlimited dimid

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID), &
                           "nc_write_model_atts",'copy inquire, '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID), &
                           "nc_write_model_atts",'time inquire, '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)'Time Dimension ID ',TimeDimID, &
  ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", errstring, source, revision, revdate)
endif
!--------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!--------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable", len=model_size, &
  dimid = StateVarDimID),"nc_write_model_atts",'StateVariable def_dim'//trim(filename))

!--------------------------------------------------------------------
! Write Global Attributes 
!--------------------------------------------------------------------
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date" ,str1    ), &
                    "nc_write_model_atts",'put_att creation_date, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source"  ,source  ), &
                    "nc_write_model_atts",'put_att model_source, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision), &
                    "nc_write_model_atts",'put_att model_revision, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate" ,revdate ), &
                    "nc_write_model_atts",'put_att model_revdate, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model"         ,"pe2lyr"), &
                    "nc_write_model_atts",'put_att model, '//trim(filename))

!--------------------------------------------------------------------
! Define the new dimensions IDs
!--------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name="lat", len = nlats, dimid = latDimID), &
                 "nc_write_model_atts",'def_dim lat, '//trim(filename)) 
call nc_check(nf90_def_dim(ncid=ncFileID, name="lon", len = nlons, dimid = lonDimID), &
                 "nc_write_model_atts",'def_dim lon, '//trim(filename)) 
call nc_check(nf90_def_dim(ncid=ncFileID, name="lev", len = nlevs, dimid = levDimID), &
                 "nc_write_model_atts",'def_dim lev, '//trim(filename)) 

!--------------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!--------------------------------------------------------------------


call nc_check(nf90_def_var(ncFileID, name="lon", xtype=nf90_double, &
     dimids=lonDimID, varid=lonVarID), &
            "nc_write_model_atts",'def_var lon, '//trim(filename) )

call nc_check(nf90_put_att(ncFileID, lonVarID, "long_name" , "longitude"), &
             "nc_write_model_atts",'put_att lon:long_name, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, lonVarID, "cartesian_axis", "X"), &
             "nc_write_model_atts",'put_att lon:cartesian_axis, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, lonVarID, "units", "degrees_east"), &
             "nc_write_model_atts",'put_att lon:units, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, lonVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)), &
             "nc_write_model_atts",'put_att lon:valid_range, '//trim(filename))


call nc_check(nf90_def_var(ncFileID, name="lat", xtype=nf90_double, &
     dimids=latDimID, varid=latVarID), &
            "nc_write_model_atts",'def_var lat, '//trim(filename) )

call nc_check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude"), &
             "nc_write_model_atts",'put_att lat:long_name, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, latVarID, "cartesian_axis", "Y"), &
             "nc_write_model_atts",'put_att lat:cartesian_axis, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north"), &
             "nc_write_model_atts",'put_att lat:units, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, latVarID, "valid_range", (/  -90.0_r8, 90.0_r8 /)), &
             "nc_write_model_atts",'put_att lat:valid_range, '//trim(filename))


call nc_check(nf90_def_var(ncFileID, name="lev", xtype=nf90_int, &
     dimids=levDimID, varid=levVarID), &
            "nc_write_model_atts",'def_var lev'//trim(filename) )

call nc_check(nf90_put_att(ncFileID, levVarID, "long_name", "level"), &
           "nc_write_model_atts",'put_att lev:long_name, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, levVarID, "cartesian_axis", "Z"), &
           "nc_write_model_atts",'put_att lev:cartesian_axis, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, levVarID, "units", " "), &
           "nc_write_model_atts",'put_att lev:units, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, levVarID, "positive", "down"), &
           "nc_write_model_atts",'put_att lev:positive, '//trim(filename))


call nc_check(nf90_def_var(ncid=ncFileID, name="u", xtype=nf90_real, &
         dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
         varid  = uVarID),"nc_write_model_atts",'def_var u, '//trim(filename))

call nc_check(nf90_put_att(ncFileID, uVarID, "long_name", "zonal wind component"), &
             "nc_write_model_atts",'put_att u:long_name, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, uVarID, "units", "m/s"), &
             "nc_write_model_atts",'put_att u:units, '//trim(filename)) 
call nc_check(nf90_put_att(ncFileID, uVarID, "_FillValue", NF90_FILL_REAL), &
             "nc_write_model_atts", 'put_att u:FillValue, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, uVarID, "missing_value", NF90_FILL_REAL), &
             "nc_write_model_atts",'put_att u:missing_value, '//trim(filename))


call nc_check(nf90_def_var(ncid=ncFileID, name="v", xtype=nf90_real, &
         dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
         varid  = vVarID),"nc_write_model_atts",'def_var v, '//trim(filename))

call nc_check(nf90_put_att(ncFileID, vVarID, "long_name", "meridional wind component"), &
             "nc_write_model_atts",'put_att v:long_name, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, vVarID, "units", "m/s"), &
             "nc_write_model_atts",'put_att v:units, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, vVarID, "_FillValue", NF90_FILL_REAL), &
             "nc_write_model_atts", 'put_att v:FillValue, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, vVarID, "missing_value", NF90_FILL_REAL), &
             "nc_write_model_atts",'put_att v:missing_value, '//trim(filename))


call nc_check(nf90_def_var(ncid=ncFileID, name="z", xtype=nf90_real, &
         dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
         varid  = zVarID),"nc_write_model_atts",'def_var z, '//trim(filename))

call nc_check(nf90_put_att(ncFileID, zVarID, "long_name", "interface height"), &
             "nc_write_model_atts",'put_att z:long_name, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, zVarID, "units", "meters"),&
             "nc_write_model_atts",'put_att z:units, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, zVarID, "_FillValue", NF90_FILL_REAL), &
             "nc_write_model_atts", 'put_att z:FillValue, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, zVarID, "missing_value", NF90_FILL_REAL), &
             "nc_write_model_atts",'put_att z:missing_value, '//trim(filename))

call nc_check(nf90_enddef(ncfileID),"nc_write_model_atts",'enddef, '//trim(filename))

!--------------------------------------------------------------------
! Fill the variables
!--------------------------------------------------------------------

call nc_check(nf90_put_var(ncFileID, lonVarID, lons(:) ), &
                "nc_write_model_atts",'put_var lons, '//trim(filename))
call nc_check(nf90_put_var(ncFileID, latVarID, lats(:) ), &
                "nc_write_model_atts",'put_var lats, '//trim(filename))
call nc_check(nf90_put_var(ncFileID, levVarID, (/ (i,i=1,nlevs) /) ), &
                "nc_write_model_atts",'put_var lev, '//trim(filename))

!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID),"nc_write_model_atts",'sync, '//trim(filename))
write(errstring,*)'netCDF file ',ncFileID,' is synched ...'
call error_handler(E_MSG,'nc_write_model_atts',errstring,source,revision, revdate)

end function nc_write_model_atts




function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
use typeSizes
use netcdf
implicit none

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: uVarID,vVarID,zVarID
integer :: mm,i,j
real(r8), dimension(nlons,nlats,2) :: ug,vg,zg
character(len=128) :: filename

ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!--------------------------------------------------------------------
! unpack the state vector into prognostic variables
!--------------------------------------------------------------------

mm =1
do i = 1,nlats
  do j=1, nlons
    ug(j,i,1) = statevec(              mm) 
    ug(j,i,2) = statevec(  nlats*nlons+mm) 
    vg(j,i,1) = statevec(2*nlats*nlons+mm) 
    vg(j,i,2) = statevec(3*nlats*nlons+mm) 
    zg(j,i,1) = statevec(4*nlats*nlons+mm) 
    zg(j,i,2) = statevec(5*nlats*nlons+mm) 
    mm = mm+1
  enddo
enddo 

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
     "nc_write_model_vars", 'Inquire, '//trim(filename))

call nc_check(NF90_inq_varid(ncFileID,  "u",  uVarID), &
       "nc_write_model_vars", 'inq_varid u, '//trim(filename))
call nc_check(nf90_put_var( ncFileID,  uVarId, ug( :, :, : ), &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ), &
               "nc_write_model_vars", 'put_var ug, '//trim(filename))

call nc_check(NF90_inq_varid(ncFileID,  "v",  vVarID), &
       "nc_write_model_vars", 'inq_varid v, '//trim(filename))
call nc_check(nf90_put_var( ncFileID,  vVarId, vg( :, :, : ), &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ), &
               "nc_write_model_vars", 'put_var vg, '//trim(filename))

call nc_check(NF90_inq_varid(ncFileID,  "z",  zVarID), &
       "nc_write_model_vars", 'inq_varid z, '//trim(filename))
call nc_check(nf90_put_var( ncFileID,  zVarId, zg( :, :, : ), &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ), &
               "nc_write_model_vars", 'put_var zg, '//trim(filename))

!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID),"nc_write_model_vars",'sync, '//trim(filename))

end function nc_write_model_vars




subroutine pert_model_state(state, pert_state, interf_provided)

! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

interf_provided = .false.
pert_state = state ! just to satisfy INTENT(OUT)

end subroutine pert_model_state




subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

end subroutine ens_mean_for_model



!#######################################################################
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

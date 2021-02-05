! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! BEGIN DART PREPROCESS TYPE DEFINITIONS
!TOWER_AIR_TEMPERATURE,          QTY_TEMPERATURE,              COMMON_CODE
!TOWER_GLOBAL_RADIATION,         QTY_RADIATION,                COMMON_CODE
!TOWER_NET_CARBON_FLUX,          QTY_NET_CARBON_FLUX,          COMMON_CODE
!TOWER_SOIL_TEMPERATURE,         QTY_TEMPERATURE,              COMMON_CODE
!TOWER_U_WIND_COMPONENT,         QTY_U_WIND_COMPONENT,         COMMON_CODE
!TOWER_V_WIND_COMPONENT,         QTY_V_WIND_COMPONENT,         COMMON_CODE
!SOIL_RESPIRATION_FLUX,          QTY_SOIL_RESPIRATION_FLUX
!TOWER_ER_FLUX,                  QTY_ER_FLUX
!TOWER_GPP_FLUX,                 QTY_GROSS_PRIMARY_PROD_FLUX
!TOWER_LATENT_HEAT_FLUX,         QTY_LATENT_HEAT_FLUX
!TOWER_NETC_ECO_EXCHANGE,        QTY_NET_CARBON_PRODUCTION
!TOWER_SENSIBLE_HEAT_FLUX,       QTY_SENSIBLE_HEAT_FLUX
! END DART PREPROCESS TYPE DEFINITIONS

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_tower_mod, only : get_scalar_from_history
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(TOWER_LATENT_HEAT_FLUX)
!     call get_scalar_from_history('EFLX_LH_TOT_R', state_handle, ens_size, &
!                    copy_indices, location, obs_time, expected_obs, istatus)
!  case(TOWER_SENSIBLE_HEAT_FLUX)
!     call get_scalar_from_history('FSH', state_handle, ens_size, &
!                    copy_indices, location, obs_time, expected_obs, istatus)
!  case(TOWER_NETC_ECO_EXCHANGE)
!     call get_scalar_from_history('NEP', state_handle, ens_size, &
!                    copy_indices, location, obs_time, expected_obs, istatus)
!  case(TOWER_GPP_FLUX)
!     call get_scalar_from_history('GPP', state_handle, ens_size, &
!                    copy_indices, location, obs_time, expected_obs, istatus)
!  case(TOWER_ER_FLUX)
!     call get_scalar_from_history('ER', state_handle, ens_size, &
!                    copy_indices, location, obs_time, expected_obs, istatus)
!  case(SOIL_RESPIRATION_FLUX)
!     call get_scalar_from_history('SR', state_handle, ens_size, &
!                    copy_indices, location, obs_time, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!    case(TOWER_LATENT_HEAT_FLUX, &
!         TOWER_SENSIBLE_HEAT_FLUX, &
!         TOWER_NETC_ECO_EXCHANGE, &
!         TOWER_GPP_FLUX, &
!         TOWER_ER_FLUX, &
!         SOIL_RESPIRATION_FLUX)
!       continue
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!    case(TOWER_LATENT_HEAT_FLUX, &
!         TOWER_SENSIBLE_HEAT_FLUX, &
!         TOWER_NETC_ECO_EXCHANGE, &
!         TOWER_GPP_FLUX, &
!         TOWER_ER_FLUX, &
!         SOIL_RESPIRATION_FLUX)
!       continue
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!    case(TOWER_LATENT_HEAT_FLUX, &
!         TOWER_SENSIBLE_HEAT_FLUX, &
!         TOWER_NETC_ECO_EXCHANGE, &
!         TOWER_GPP_FLUX, &
!         TOWER_ER_FLUX, &
!         SOIL_RESPIRATION_FLUX)
!       continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_tower_mod

! This is the forward observation operator code for CLM for flux tower observations.
! The DART model state for CLM generally does not have all the pieces necessary
! to apply the forward observation operator directly, so this code gets what it
! needs from a CLM history file. These history files have become complicated now
! that CLM is trying to support unstructured grids. Sometimes the variables of
! interest are shaped NEP(time, lat, lon), sometimes NEP(time, lndgrid).
! 'single column' runs may appear as either lat=lon=1 or lndgrid=1
!
use        types_mod, only : r4, r8, digits12, MISSING_R8, PI

use     location_mod, only : location_type, get_location, get_dist, &
                             set_location, write_location, VERTISUNDEF

use time_manager_mod, only : time_type, get_date, set_date, print_date, print_time, &
                             get_time, set_time, operator(-), operator(/=)

use    utilities_mod, only : register_module, E_ERR, E_MSG, error_handler, &
                             check_namelist_read, find_namelist_in_file,   &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             file_exist, is_longitude_between

use netcdf_utilities_mod, only : nc_check

use      assim_model_mod, only : interpolate

use ensemble_manager_mod, only : ensemble_type

use typesizes
use netcdf

implicit none
private

public :: get_scalar_from_history

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_tower_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical            :: module_initialized = .false.
logical            :: unstructured = .false.
character(len=512) :: string1, string2, string3
integer            :: nlon, nlat, ntime, ens_size
type(time_type)    :: initialization_time

character(len=256), allocatable, dimension(:) :: fname
integer,            allocatable, dimension(:) :: ncid
real(r8),           allocatable, dimension(:) :: lon, lat, area
real(digits12),     allocatable, dimension(:) :: rtime

real(r8), parameter :: RAD2KM = 40030.0_r8/(2.0_r8 * PI) ! (mean radius of earth ~6371km)

! namelist items
character(len=256) :: casename = 'clm_dart'
logical            :: debug = .false.
integer            :: hist_nhtfrq = -24
! CLM variable hist_nhtfrq ... controls how often to write out the history files.
! Negative value means the output frequency is the absolute value (in hours).

namelist /obs_def_tower_nml/ casename, debug, hist_nhtfrq

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of executable routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine initialize_module(model_time)
type(time_type), intent(in) :: model_time

! Called once to set values and allocate space, open all the CLM files
! that have the observations, etc.

integer :: iunit, io, i
integer :: varid
integer :: year, month, day, hour, minute, second, leftover
integer, allocatable, dimension(:) :: yyyymmdd,sssss
type(time_type) :: tower_time

! Prevent multiple calls from executing this code more than once.
if (module_initialized) then
   if (initialization_time /= model_time) then
      string1 = 'model time does not match initialization time'
      string2 = 'model time does not match initialization time'
      string3 = 'model time does not match initialization time'
      call error_handler(E_ERR, 'obs_def_tower.initialize_routine', string1, &
                     source, revision, revdate, text2=string2,text3=string3)
   endif
   return
else
   initialization_time = model_time
endif

module_initialized = .true.

! Log the version of this source file.
call register_module(source, revision, revdate)

! Read the namelist entry.
call find_namelist_in_file("input.nml", "obs_def_tower_nml", iunit)
read(iunit, nml = obs_def_tower_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_tower_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_def_tower_nml)
if (do_nml_term()) write(     *     , nml=obs_def_tower_nml)

! Need to know what day we are trying to assimilate.
! The CLM h1 files contain everything STARTING with the time in their filename.
! all output intervals from that run are simply appended to that file.
! Consequently, we need to know the filename from the START of the model advance
! that resulted in the current model state. To construct the filename, we need
! to know one of the namelist variables from CLM. We are mandating that this
! value gets passed to DART via the obs_def_tower_nml instead of reading
! the CESM namelist .... hist_nhtfrq must be a negative number in a DART
! application of CLM ... and the STOP_OPTION must be "HOURS".
!
!  | start of the model advance ... *.h1.* file starts getting written
!  |
!  X==X==X==X==X==X==X==X==X==X==X==[O] (CLM model advance)
!  |<------- hist_nhtfrq ------->|   |
!                                |   | END of the model advance
!                                |
!                                | END of the data in the *.h1.* file

second = abs(hist_nhtfrq)*60*60
tower_time = model_time - set_time(second,0)

call get_date(tower_time, year, month, day, hour, minute, second)
second = second + minute*60 + hour*3600

! Figure out how many files (i.e. ensemble size) and construct their names.
! The CLM h0 files are constructed such that the midnight that starts the
! day is IN the file. The last time in the file is 23:30 ...

! In a perfect model scenario (i.e. 1 instance) CESM constructs filenames
! without the instance counter in them. We check for that first. If it
! exists we know the filename and ensemble size.

100 format (A,'.clm2_',I4.4,'.h1.',I4.4,'-',I2.2,'-',I2.2,'-',I5.5,'.nc')
110 format (A,'.clm2'      ,'.h1.',I4.4,'-',I2.2,'-',I2.2,'-',I5.5,'.nc')

write(string1,110) trim(casename),year,month,day,second

if( file_exist(string1) ) then
   ! We know we are in a perfect model scenario
   ens_size = 1
   if (debug .and. do_output()) write(*,*)'Ensemble size is believed to be ',ens_size
   goto 200
endif

ens_size = 0
ENSEMBLESIZE : do i = 1,200

   write(string1,100) trim(casename),i,year,month,day,second

   if( file_exist(string1) ) then
      if(debug .and. do_output()) &
            write(*,*)'observation file "',trim(string1),'" exists.'
      ens_size = ens_size + 1
   else
      if(debug .and. do_output()) &
            write(*,*)'WARNING observation file "',trim(string1),'" does not exist.'
      exit ENSEMBLESIZE
   endif

enddo ENSEMBLESIZE

if (ens_size < 2) then

   write(string1,110) trim(casename),year,month,day,second
   write(string2,*)'cannot find files to use for observation operator.'
   write(string3,*)'trying files with names like "',trim(string1),'"'
   call error_handler(E_ERR, 'obs_def_tower.initialize_routine', string2, &
                  source, revision, revdate, text2=string3)

elseif (ens_size >= 200) then

   write(string2,*)'ensemble size (',ens_size,') is unnaturally large.'
   write(string3,*)'trying files with names like "',trim(string1),'"'
   call error_handler(E_ERR, 'obs_def_tower.initialize_routine', string2, &
                  source, revision, revdate, text2=string3)

else

   if (debug .and. do_output()) then
      write(string1,*)'Ensemble size is believed to be ',ens_size
      call error_handler(E_MSG, 'obs_def_tower.initialize_routine', string1, &
                  source, revision, revdate )
   endif

endif

200 continue

allocate(fname(ens_size),ncid(ens_size))
ncid = 0

ENSEMBLE : do i = 1,ens_size
   if (ens_size == 1) then
      write(fname(i),110) trim(casename),year,month,day,second
      call nc_check(nf90_open(trim(fname(i)), nf90_nowrite, ncid(i)), &
          'obs_def_tower.initialize_routine','open '//trim(fname(i)))
      exit ENSEMBLE
   else
      write(fname(i),100) trim(casename),i,year,month,day,second
      call nc_check(nf90_open(trim(fname(i)), nf90_nowrite, ncid(i)), &
          'obs_def_tower.initialize_routine','open '//trim(fname(i)))
   endif
enddo ENSEMBLE

i = 1

! FIXME All other files will be opened to make sure they have the same dimensions.

call GetDimensions(ncid(i), fname(i))  ! determines values for nlon, nlat, ntime

allocate(lon(nlon), lat(nlat), rtime(ntime), yyyymmdd(ntime), sssss(ntime))

call nc_check(nf90_inq_varid(ncid(i), 'lon',    varid), &
              'obs_def_tower.initialize_routine','inq_varid lon '//trim(fname(i)))
call nc_check(nf90_get_var(  ncid(i), varid,      lon), &
              'obs_def_tower.initialize_routine','get_var lon'//trim(fname(i)))

call nc_check(nf90_inq_varid(ncid(i), 'lat',    varid), &
              'obs_def_tower.initialize_routine','inq_varid lat '//trim(fname(i)))
call nc_check(nf90_get_var(  ncid(i), varid,      lat), &
              'obs_def_tower.initialize_routine','get_var lat'//trim(fname(i)))

call nc_check(nf90_inq_varid(ncid(i), 'mcdate', varid), &
              'obs_def_tower.initialize_routine','inq_varid mcdate '//trim(fname(i)))
call nc_check(nf90_get_var(  ncid(i), varid, yyyymmdd), &
              'obs_def_tower.initialize_routine','get_var yyyymmdd'//trim(fname(i)))

call nc_check(nf90_inq_varid(ncid(i), 'mcsec',  varid), &
              'obs_def_tower.initialize_routine','inq_varid mcsec '//trim(fname(i)))
call nc_check(nf90_get_var(  ncid(i), varid,    sssss), &
              'obs_def_tower.initialize_routine','get_var sssss'//trim(fname(i)))

if ( (nlon == 1) .and. (nlat ==1) ) then
   allocate(area(nlon))
   call nc_check(nf90_inq_varid(ncid(i), 'area', varid), &
        'obs_def_tower.initialize_routine','inq_varid area '//trim(fname(i)))
   call nc_check(nf90_get_var(ncid(i), varid, area), &
        'obs_def_tower.initialize_routine', 'get_var area'//trim(fname(i)))
   if (debug .and. do_output()) write(*,*)'obs_def_tower      area',area(nlon)
endif

! Convert time in file to a time compatible with the observation sequence file.
do i = 1,ntime

   year     = yyyymmdd(i)/10000
   leftover = yyyymmdd(i) - year*10000
   month    = leftover/100
   day      = leftover - month*100

   hour     = sssss(i)/3600
   leftover = sssss(i) - hour*3600
   minute   = leftover/60
   second   = leftover - minute*60

   tower_time = set_date(year, month, day, hour, minute, second)
   call get_time(tower_time, second, day)

   rtime(i) = real(day,digits12) + real(second,digits12)/86400.0_digits12

   if (debug .and. do_output()) then
      write(*,*)'timestep yyyymmdd sssss',i,yyyymmdd(i),sssss(i)
      call print_date(tower_time,'tower_mod date')
      call print_time(tower_time,'tower_mod time')
      write(*,*)'tower_mod time as a real ',rtime(i)
   endif

enddo

if (debug .and. do_output()) write(*,*)'obs_def_tower      lon',lon
if (debug .and. do_output()) write(*,*)'obs_def_tower      lat',lat

!>@todo FIXME check all other ensemble member history files to make sure metadata is the same.

deallocate(yyyymmdd, sssss)

end subroutine initialize_module


!======================================================================


subroutine GetDimensions(ncid, fname)

! Harvest information from the first observation file.
! The SingleColumMode files have
!        float lat(lndgrid) ;
!        float lon(lndgrid) ;
!        float area(lndgrid) ;
! while the 2D files have
!        float lat(lat) ;
!        float lon(lon) ;
!        float area(lat, lon) ;

integer,          intent(in) :: ncid
character(len=*), intent(in) :: fname

! integer, intent(out) :: nlon, nlat, ntime ... module variables

integer, dimension(NF90_MAX_VAR_DIMS) :: londimids, latdimids
integer :: lonvarid, lonndims
integer :: latvarid, latndims
integer :: dimid

call nc_check(nf90_inq_varid(ncid,  'lon', lonvarid), &
              'obs_def_tower.GetDimensions','inq_varid lon '//trim(fname))
call nc_check(nf90_inq_varid(ncid,  'lat', latvarid), &
              'obs_def_tower.GetDimensions','inq_varid lat '//trim(fname))

call nc_check(nf90_inquire_variable(ncid, lonvarid, ndims=lonndims, dimids=londimids),&
              'obs_def_tower.GetDimensions','inquire lon '//trim(fname))
call nc_check(nf90_inquire_variable(ncid, latvarid, ndims=latndims, dimids=latdimids),&
              'obs_def_tower.GetDimensions','inquire lat '//trim(fname))

if ( (lonndims /= 1) .or. (latndims /= 1) ) then
   write(string1,*) 'Require "lon" and "lat" variables to be 1D. They are ', &
                     lonndims, latndims
   call error_handler(E_ERR,'obs_def_tower.GetDimensions',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(ncid, londimids(1), len=nlon), &
              'obs_def_tower.GetDimensions','inquire_dimension lon '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, latdimids(1), len=nlat), &
              'obs_def_tower.GetDimensions','inquire_dimension lat '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'time', dimid), &
              'obs_def_tower.GetDimensions','inq_dimid time '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=ntime), &
              'obs_def_tower.GetDimensions','inquire_dimension time '//trim(fname))

if ( nf90_inq_dimid(ncid, 'lndgrid', dimid) == NF90_NOERR) then
   unstructured = .true.
   if ((nlon /= 1) .or. (nlat /= 1)) then
      string1 = 'unstructured grids with more than a single gridcell are not supported.'
      call error_handler(E_ERR,'obs_def_tower.GetDimensions',string1,source,revision,revdate)
   endif
endif

end subroutine GetDimensions


!======================================================================

subroutine get_scalar_from_history(varstring, state_handle, ens_size, copy_indices, location, &
                                   obs_time, obs_val, istatus)

character(len=*),    intent(in)  :: varstring
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: copy_indices(ens_size) ! global ens_indicies
type(location_type), intent(in)  :: location
type(time_type),     intent(in)  :: obs_time
real(r8),            intent(out) :: obs_val(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer :: imem

! This initialize_module call is going to open all files on all tasks
! when doing a distributed forward operator.
if ( .not. module_initialized ) call initialize_module(state_handle%current_time)

if ( unstructured ) then
   do imem = 1, ens_size
      call get_scalar_from_2Dhistory(varstring, copy_indices(imem), location, obs_time, obs_val(imem), istatus(imem))
   enddo
else
   do imem = 1, ens_size
      call get_scalar_from_3Dhistory(varstring, copy_indices(imem), location, obs_time,obs_val(imem), istatus(imem))
   enddo
endif

end subroutine get_scalar_from_history


!======================================================================


subroutine get_scalar_from_3Dhistory(varstring, ens_index, location, obs_time, &
                                     obs_val, istatus)
! the routine must return values for:
! obs_val -- the computed forward operator value
! istatus -- return code: 0=ok, > 0 is error, < 0 reserved for system use
!
! The requirement is that the history file variable is a 3D variable shaped similarly:
!
! float NEP(time, lat, lon) ;
!          NEP:long_name = "net ecosystem production, blah, blah, blah" ;
!          NEP:units = "gC/m^2/s" ;
!          NEP:cell_methods = "time: mean" ;
!          NEP:_FillValue = 1.e+36f ;
!          NEP:missing_value = 1.e+36f ;

character(len=*),    intent(in)  :: varstring
integer,             intent(in)  :: ens_index
type(location_type), intent(in)  :: location
type(time_type),     intent(in)  :: obs_time
real(r8),            intent(out) :: obs_val
integer,             intent(out) :: istatus

integer,  dimension(NF90_MAX_VAR_DIMS) :: dimids
real(r8), dimension(3) :: loc
integer,  dimension(3) :: ncstart, nccount
integer,  dimension(1) :: loninds, latinds, timeinds
integer                :: gridloni, gridlatj, timei
integer                :: varid, xtype, ndims, natts, dimlen
integer                :: io1, io2, second, day
real(r8)               :: loc_lon, loc_lat, radius, distance
real(r4), dimension(1) :: hyperslab
real(r4)               :: spvalR4
real(r8)               :: scale_factor, add_offset
real(digits12)         :: otime
character(len=NF90_MAX_NAME+20)      :: strshort

type(location_type) :: gridloc

obs_val = MISSING_R8
istatus = 1

!----------------------------------------------------------------------
! if observation is outside region encompassed in the history file - fail
loc      = get_location(location) ! loc is in DEGREES
loc_lon  = loc(1)
loc_lat  = loc(2)

if ( (nlon==1) .and. (nlat==1) ) then

   ! Defining the region if running in an unstructured grid is tricky.
   ! Have lat, lon, and the area of the gridcell which we assume to be basically square.
   ! The square root of the area defines the length of the edge of the gridcell.
   ! Half the hypotenuse defines the radius of a circle. Any ob within
   ! that radius is close enough.

   gridloc   = set_location(lon(1),lat(1), 0.0_r8, VERTISUNDEF)
   distance  = get_dist(gridloc, location, no_vert = .TRUE.) * RAD2KM ! planet earth
   radius    = sqrt(2.0_r8 * area(1))/2.0_r8

   if (debug .and. do_output()) then
      write(string1,*)'    observation lon, lat is ',loc_lon, loc_lat
      write(string2,*)'gridcell    lon, lat is ',lon(1),lat(1)
      write(string3,*)'area,radius is ',area(1),radius,' distance ',distance
      call error_handler(E_MSG, 'obs_def_tower.get_scalar_from_3Dhistory', &
                 string1, source, revision, revdate, text2=string2, text3=string3)
   endif

   if ( distance > radius ) return

else
   if ( .not. is_longitude_between(loc_lon, lon(1), lon(nlon), doradians=.FALSE.)) return
   if ((loc_lat < lat(1)) .or. (loc_lat > lat(nlat))) return
endif

!----------------------------------------------------------------------
! Now that we know the observation operator is possible, continue ...

write(strshort,'(''ens_index '',i4,1x,A)')ens_index,trim(varstring)

if (ens_index > ens_size) then
   write(string1,*)'Known to have ',ens_size,'ensemble members for observation operator.'
   write(string2,*)'asking to use operator for ensemble member ',ens_index
   call error_handler(E_ERR, 'obs_def_tower.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate, text2=string2)
endif

!----------------------------------------------------------------------
! bombproofing ... make sure the netcdf file is open.

call nc_check(nf90_inquire(ncid(ens_index)), &
              'obs_def_tower.get_scalar_from_3Dhistory', 'inquire '//trim(strshort))

! bombproofing ... make sure the variable is the shape and size we expect

call nc_check(nf90_inq_varid(ncid(ens_index), trim(varstring), varid), &
        'obs_def_tower.get_scalar_from_3Dhistory', 'inq_varid '//trim(strshort))
call nc_check(nf90_inquire_variable(ncid(ens_index), varid, xtype=xtype, ndims=ndims, &
        dimids=dimids, natts=natts), &
        'obs_def_tower.get_scalar_from_3Dhistory','inquire variable '//trim(strshort))

if (ndims /= 3) then
   write(string1,*)trim(varstring),' is supposed to have 3 dimensions, it has',ndims
   call error_handler(E_ERR, 'obs_def_tower.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! If the variable is not a NF90_FLOAT, then the assumptions for processing
! the missing_value, _FillValue, etc., may not be correct.
if (xtype /= NF90_FLOAT) then
   write(string1,*)trim(varstring),' is supposed to be a 32 bit real. xtype = ', &
                   NF90_FLOAT,' it is ',xtype
   call error_handler(E_ERR, 'obs_def_tower.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 1 is longitude
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(1), len=dimlen), &
        'obs_def_tower.get_scalar_from_3Dhistory', 'inquire_dimension 1 '//trim(strshort))
if (dimlen /= nlon) then
   write(string1,*)'LON has length',nlon,trim(varstring),' has ',dimlen,'longitudes.'
   call error_handler(E_ERR, 'obs_def_tower.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 2 is latitude
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(2), len=dimlen), &
        'obs_def_tower.get_scalar_from_3Dhistory', 'inquire_dimension 2 '//trim(strshort))
if (dimlen /= nlat) then
   write(string1,*)'LAT has length',nlat,trim(varstring),' has ',dimlen,'latitudes.'
   call error_handler(E_ERR, 'obs_def_tower.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 3 is time
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(3), len=dimlen), &
        'obs_def_tower.get_scalar_from_3Dhistory', 'inquire_dimension 3'//trim(strshort))
if (dimlen /= ntime) then
   write(string1,*)'TIME has length',ntime,trim(varstring),' has ',dimlen,'times.'
   call error_handler(E_ERR, 'obs_def_tower.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

!----------------------------------------------------------------------
! Find the grid cell and timestep of interest
! FIXME ... since the history file contents are for the previous 30m,
! perhaps the closest time is not the best approximation.
! Get the individual locations values

call get_time(obs_time, second, day)
otime    = real(day,digits12) + real(second,digits12)/86400.0_digits12

latinds  = minloc(abs(lat - loc_lat))   ! these return 'arrays' ...
loninds  = minloc(abs(lon - loc_lon))   ! these return 'arrays' ...
timeinds = minloc(abs(rtime - otime))   ! these return 'arrays' ...

gridlatj = latinds(1)
gridloni = loninds(1)
timei    = timeinds(1)

if (debug .and. do_output()) then
   write(*,*)'obs_def_tower.get_scalar_from_3Dhistory:targetlon, lon, lon index is ', &
                                           loc_lon,lon(gridloni),gridloni
   write(*,*)'obs_def_tower.get_scalar_from_3Dhistory:targetlat, lat, lat index is ', &
                                           loc_lat,lat(gridlatj),gridlatj
   write(*,*)'obs_def_tower.get_scalar_from_3Dhistory:  targetT,   T,   T index is ', &
                                           otime,rtime(timei),timei
endif

if ( abs(otime - rtime(timei)) > 30*60 ) then
   if (debug .and. do_output()) then
      write(*,*)'obs_def_tower.get_scalar_from_3Dhistory: no close time ... skipping observation'
      call print_time(obs_time,'obs_def_tower.get_scalar_from_3Dhistory:observation time')
      call print_date(obs_time,'obs_def_tower.get_scalar_from_3Dhistory:observation date')
   endif
   istatus = 2
   return
endif

!----------------------------------------------------------------------
! Grab exactly the scalar we want.

ncstart = (/ gridloni, gridlatj, timei /)
nccount = (/        1,        1,     1 /)

call nc_check(nf90_get_var(ncid(ens_index),varid,hyperslab,start=ncstart,count=nccount), &
     'obs_def_tower.get_scalar_from_3Dhistory', 'get_var')

obs_val = hyperslab(1)

!----------------------------------------------------------------------
! Apply any netCDF attributes ...

io1 = nf90_get_att(ncid(ens_index), varid, '_FillValue' , spvalR4)
if ((io1 == NF90_NOERR) .and. (hyperslab(1) == spvalR4)) obs_val = MISSING_R8

io2 = nf90_get_att(ncid(ens_index), varid, 'missing_value' , spvalR4)
if ((io2 == NF90_NOERR) .and. (hyperslab(1) == spvalR4)) obs_val = MISSING_R8

io1 = nf90_get_att(ncid(ens_index), varid, 'scale_factor', scale_factor)
io2 = nf90_get_att(ncid(ens_index), varid, 'add_offset'  , add_offset)

if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
   if (obs_val /= MISSING_R8) obs_val = obs_val * scale_factor + add_offset
elseif (io1 == NF90_NOERR) then
   if (obs_val /= MISSING_R8) obs_val = obs_val * scale_factor
elseif (io2 == NF90_NOERR) then
   if (obs_val /= MISSING_R8) obs_val = obs_val + add_offset
endif

if (obs_val /= MISSING_R8) istatus = 0

end subroutine get_scalar_from_3Dhistory


!======================================================================


subroutine get_scalar_from_2Dhistory(varstring, ens_index, location, obs_time, &
                                     obs_val, istatus)
! the routine must return values for:
! obs_val -- the computed forward operator value
! istatus -- return code: 0=ok, > 0 is error, < 0 reserved for system use
!
! The requirement is that the history file variable is a 2D variable shaped similarly:
!
! float NEP(time, lndgrid) ;
!          NEP:long_name = "net ecosystem production, blah, blah, blah" ;
!          NEP:units = "gC/m^2/s" ;
!          NEP:cell_methods = "time: mean" ;
!          NEP:_FillValue = 1.e+36f ;
!          NEP:missing_value = 1.e+36f ;
!
! Just because it is 2D does not mean it is a single column,
! although single columns are all that is really supported right now.

character(len=*),    intent(in)  :: varstring
integer,             intent(in)  :: ens_index
type(location_type), intent(in)  :: location
type(time_type),     intent(in)  :: obs_time
real(r8),            intent(out) :: obs_val
integer,             intent(out) :: istatus

integer,  dimension(NF90_MAX_VAR_DIMS) :: dimids
real(r8), dimension(3) :: loc
integer,  dimension(2) :: ncstart, nccount
integer,  dimension(1) :: timeinds
integer                :: gridij, timei
integer                :: varid, xtype, ndims, natts, dimlen
integer                :: io1, io2, second, day
real(r8)               :: loc_lon, loc_lat, radius, distance
real(r4), dimension(1) :: hyperslab
real(r4)               :: spvalR4
real(r8)               :: scale_factor, add_offset
real(digits12)         :: otime
character(len=NF90_MAX_NAME+20)      :: strshort

type(location_type) :: gridloc

obs_val = MISSING_R8
istatus = 1

!----------------------------------------------------------------------
! if observation is outside region encompassed in the history file - fail
loc      = get_location(location) ! loc is in DEGREES
loc_lon  = loc(1)
loc_lat  = loc(2)

! Defining the region if running in an unstructured grid is tricky.
! We have lat, lon, and the area of the gridcell which we assume to be basically square.
! The square root of the area defines the length of the edge of the gridcell.
! Half the hypotenuse defines the radius of a circle. Any ob within
! that radius is close enough.

! TJH FIXME This does not work with unstructured grid.
! latinds  = minloc(abs(lat - loc_lat))   ! these return 'arrays' ...
! loninds  = minloc(abs(lon - loc_lon))   ! these return 'arrays' ...
! gridij = the closest location

! The "1" in the following reflect the fact that only a single gridcell
! is currently supported in the unstructured grid configuration.
! (see GetDimensions())

gridij   = 1
gridloc  = set_location(lon(gridij),lat(gridij), 0.0_r8, VERTISUNDEF)
distance = get_dist(gridloc, location, no_vert = .TRUE.) * RAD2KM
radius   = sqrt(2.0_r8 * area(gridij))/2.0_r8

if (debug .and. do_output()) then
   write(string1,*)'    observation lon, lat is ',loc_lon, loc_lat
   write(string2,*)'gridcell    lon, lat is ',lon(gridij),lat(gridij)
   write(string3,*)'area,radius is ',area(gridij),radius,' distance ',distance
   call error_handler(E_MSG, 'obs_def_tower.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate, text2=string2, text3=string3)
endif

if ( distance > radius ) return

!----------------------------------------------------------------------
! Now that we know the observation operator is possible, continue ...

write(strshort,'(''ens_index '',i4,1x,A)')ens_index,trim(varstring)

if (ens_index > ens_size) then
   write(string1,*)'believed to have ',ens_size,'ensemble members for observation operator.'
   write(string2,*)'asking to use operator for ensemble member ',ens_index
   call error_handler(E_ERR, 'obs_def_tower.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate, text2=string2)
endif

!----------------------------------------------------------------------
! bombproofing ... make sure the netcdf file is open.

call nc_check(nf90_inquire(ncid(ens_index)), &
              'obs_def_tower.get_scalar_from_2Dhistory', 'inquire '//trim(strshort))

! bombproofing ... make sure the variable is the shape and size we expect

call nc_check(nf90_inq_varid(ncid(ens_index), trim(varstring), varid), &
        'obs_def_tower.get_scalar_from_2Dhistory', 'inq_varid '//trim(strshort))
call nc_check(nf90_inquire_variable(ncid(ens_index), varid, xtype=xtype, ndims=ndims, &
        dimids=dimids, natts=natts), &
        'obs_def_tower.get_scalar_from_2Dhistory','inquire variable '//trim(strshort))

if (ndims /= 2) then
   write(string1,*)trim(varstring),' is supposed to have 2 dimensions, it has',ndims
   call error_handler(E_ERR, 'obs_def_tower.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate)
endif

! If the variable is not a NF90_FLOAT, then the assumptions for processing
! the missing_value, _FillValue, etc., may not be correct.
if (xtype /= NF90_FLOAT) then
   write(string1,*)trim(varstring),' is supposed to be a 32 bit real. xtype = ', &
                   NF90_FLOAT,' it is ',xtype
   call error_handler(E_ERR, 'obs_def_tower.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 1 is spatial
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(1), len=dimlen), &
        'obs_def_tower.get_scalar_from_2Dhistory', 'inquire_dimension 1 '//trim(strshort))
if (dimlen /= nlon) then
   write(string1,*)'LON has length',nlon,trim(varstring),' has ',dimlen,'longitudes.'
   call error_handler(E_ERR, 'obs_def_tower.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 2 is time
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(2), len=dimlen), &
        'obs_def_tower.get_scalar_from_2Dhistory', 'inquire_dimension 2 '//trim(strshort))
if (dimlen /= ntime) then
   write(string1,*)'TIME has length',ntime,trim(varstring),' has ',dimlen,'times.'
   call error_handler(E_ERR, 'obs_def_tower.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate)
endif

!----------------------------------------------------------------------
! Find the timestep of interest
! FIXME ... since the history file contents are for the previous 30m,
! perhaps the closest time is not the best approximation.

call get_time(obs_time, second, day)
otime    = real(day,digits12) + real(second,digits12)/86400.0_digits12
timeinds = minloc(abs(rtime - otime))   ! these return 'arrays' ...
timei    = timeinds(1)

if (debug .and. do_output()) then
   write(*,*)'obs_def_tower.get_scalar_from_2Dhistory:  targetT,   T,   T index is ', &
                                           otime,rtime(timei),timei
endif

if ( abs(otime - rtime(timei)) > 30*60 ) then
   if (debug .and. do_output()) then
      write(*,*)'obs_def_tower.get_scalar_from_2Dhistory: no close time ... skipping observation'
      call print_time(obs_time,'obs_def_tower.get_scalar_from_2Dhistory:observation time')
      call print_date(obs_time,'obs_def_tower.get_scalar_from_2Dhistory:observation date')
   endif
   istatus = 2
   return
endif

!----------------------------------------------------------------------
! Grab exactly the scalar we want.

ncstart = (/ gridij, timei /)
nccount = (/      1,     1 /)

call nc_check(nf90_get_var(ncid(ens_index),varid,hyperslab,start=ncstart,count=nccount), &
     'obs_def_tower.get_scalar_from_2Dhistory', 'get_var')

obs_val = hyperslab(1)

!----------------------------------------------------------------------
! Apply any netCDF attributes ...

io1 = nf90_get_att(ncid(ens_index), varid, '_FillValue' , spvalR4)
if ((io1 == NF90_NOERR) .and. (hyperslab(1) == spvalR4)) obs_val = MISSING_R8

io2 = nf90_get_att(ncid(ens_index), varid, 'missing_value' , spvalR4)
if ((io2 == NF90_NOERR) .and. (hyperslab(1) == spvalR4)) obs_val = MISSING_R8

io1 = nf90_get_att(ncid(ens_index), varid, 'scale_factor', scale_factor)
io2 = nf90_get_att(ncid(ens_index), varid, 'add_offset'  , add_offset)

if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
   if (obs_val /= MISSING_R8) obs_val = obs_val * scale_factor + add_offset
elseif (io1 == NF90_NOERR) then
   if (obs_val /= MISSING_R8) obs_val = obs_val * scale_factor
elseif (io2 == NF90_NOERR) then
   if (obs_val /= MISSING_R8) obs_val = obs_val + add_offset
endif

if (obs_val /= MISSING_R8) istatus = 0

end subroutine get_scalar_from_2Dhistory


!======================================================================


subroutine test_block
! Defining the region if running in a single column is tricky.
! We have lat, lon, and the area of the gridcell which we assume to be basically square.
! The square root of the area defines the length of the edge of the gridcell which
! can then be interpreted as the diameter of a circle. Any observation within this
! distance is close enough.

real(r8) :: x1, x2, y1, y2, d1, d2, d3, d4, radius, distance
type(location_type) :: gridloc, testloc

gridloc  = set_location(0.0_r8, 0.0_r8, 0.0_r8, VERTISUNDEF)
testloc  = set_location(1.0_r8, 0.0_r8, 0.0_r8, VERTISUNDEF)
distance = get_dist(gridloc, testloc, no_vert = .TRUE.) * RAD2KM
write(*,*)'TJH DEBUG: 1 degree at the equator has distance ',distance,' km.'

x1 = 286.8750               ! gridcell longitude is 287.50
x2 = 288.1250
y1 = 42.40837669372559      ! gridcell latitude is 42.8795814514160
y2 = 43.35078620910645

! compute the distance along the top of the grid cell
gridloc = set_location(x1, y1, 0.0_r8, VERTISUNDEF)
testloc = set_location(x2, y1, 0.0_r8, VERTISUNDEF)
d1      = get_dist(gridloc, testloc, no_vert = .TRUE.) * RAD2KM
gridloc = set_location(x2, y2, 0.0_r8, VERTISUNDEF)
d2      = get_dist(gridloc, testloc, no_vert = .TRUE.) * RAD2KM
testloc = set_location(x1, y2, 0.0_r8, VERTISUNDEF)
d3      = get_dist(gridloc, testloc, no_vert = .TRUE.) * RAD2KM
gridloc = set_location(x1, y1, 0.0_r8, VERTISUNDEF)
d4      = get_dist(gridloc, testloc, no_vert = .TRUE.) * RAD2KM

write(*,*)
write(*,*)'lengths         are ',d1,d2,d3,d4
write(*,*)'rectangular area is ',((d1+d3)/2.0_r8)*((d2+d4)/2.0_r8)

d1 = sqrt(area(1))
d2 = sqrt(2.0_r8 * d1**2)
d3 = d2 / 2.0_r8
radius = sqrt(2.0_r8 * area(1))/2.0_r8
write(*,*)'length of one side of a square is x = sqrt(area) = ',d1
write(*,*)'diagonal is                sqrt(x**2 + x**2)     = ',d2
write(*,*)'radius   is                sqrt(x**2 + x**2)/2.0 = ',d3
write(*,*)'-or-                       sqrt(area + area)/2.0 = ',radius

write(*,*)'TJH DEBUG: Radius,area is ',distance,PI*distance**2, &
          ' gridcell area is ',area(1)

end subroutine test_block


!======================================================================


end module obs_def_tower_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------


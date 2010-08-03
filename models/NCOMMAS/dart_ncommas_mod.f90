! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module dart_ncommas_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8, rad2deg, PI, SECPERDAY
use time_manager_mod, only : time_type, get_date, set_date, get_time, set_time, &
                             print_date, print_time, &
                             operator(==), operator(-), operator(+)
use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler, nc_check, &
                             find_namelist_in_file, check_namelist_read, &
                             E_ERR, E_MSG, timestamp, find_textfile_dims, &
                             logfileunit, do_output

use typesizes
use netcdf

implicit none
private

public :: set_model_time_step, grid_type, get_grid_dims, get_grid, &
          get_base_time, get_state_time, &
          write_ncommas_namelist, get_ncommas_restart_filename

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = '$URL$', &
   revision = '$Revision$', &
   revdate  = '$Date$'

character(len=256) :: string1, string2
logical, save :: module_initialized = .false.

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.
logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

type grid_type
   private
   integer  :: nx, ny, nz         ! determines if by level, height, pressure, ...
   real(r8), pointer ::  lon(:,:) ! lon stored in radians
   real(r8), pointer ::  lat(:,:) ! lat stored in radians
   real(r8), pointer :: vloc(:,:) ! height stored in meters
end type grid_type

type(grid_type) :: Ugrid
type(grid_type) :: Vgrid
type(grid_type) :: Wgrid
type(grid_type) ::  grid

!------------------------------------------------------------------
! The ncommas restart manager namelist variables
!------------------------------------------------------------------

character(len=256) :: ic_filename      = 'ncommas.nc'
!character(len=256) :: restart_filename = 'dart_ncommas_mod_restart_filename_not_set'
character(len= 64) :: ew_boundary_type, ns_boundary_type

namelist /restart_nml/ ic_filename, ew_boundary_type, ns_boundary_type

INTERFACE get_base_time
      MODULE PROCEDURE get_base_time_ncid
      MODULE PROCEDURE get_base_time_fname
END INTERFACE

INTERFACE get_state_time
      MODULE PROCEDURE get_state_time_ncid
      MODULE PROCEDURE get_state_time_fname
END INTERFACE

!======================================================================
contains
!======================================================================


subroutine initialize_module
!------------------------------------------------------------------
integer :: iunit, io

! Make sure we have a ncommas restart file (for grid dims)
if ( .not. file_exist(ic_filename) ) then
   string1 = trim(ic_filename)//' not found'
   call error_handler(E_ERR,'initialize_module', &
          string1, source, revision, revdate)
endif

module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

end subroutine initialize_module



subroutine get_grid_dims(NXC, NXE, NYC, NYE, NZC, NZE)
!------------------------------------------------------------------
!
! Read the grid dimensions from the restart netcdf file.
!
! The file name comes from module storage ... namelist.

integer, intent(out) :: NXC   ! Number of Longitude centers
integer, intent(out) :: NXE   ! Number of Longitude edges
integer, intent(out) :: NYC   ! Number of Latitude  centers
integer, intent(out) :: NYE   ! Number of Latitude  edges
integer, intent(out) :: NZC   ! Number of Vertical grid centers
integer, intent(out) :: NZE   ! Number of Vertical grid edges

integer :: grid_id, dimid

if ( .not. module_initialized ) call initialize_module

! get the ball rolling ...

call nc_check(nf90_open(trim(ic_filename), nf90_nowrite, grid_id), &
            'get_grid_dims','open '//trim(ic_filename))

! Longitudes : get dimid for 'XC' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'XC', dimid), &
            'get_grid_dims','inq_dimid XC '//trim(ic_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NXC), &
            'get_grid_dims','inquire_dimension XC '//trim(ic_filename))

! Longitudes : get dimid for 'XE and then get value

call nc_check(nf90_inq_dimid(grid_id, 'XE', dimid), &
            'get_grid_dims','inq_dimid XE '//trim(ic_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NXE), &
            'get_grid_dims','inquire_dimension XE '//trim(ic_filename))

! Latitudes : get dimid for 'YC' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'YC', dimid), &
            'get_grid_dims','inq_dimid YC '//trim(ic_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NYC), &
            'get_grid_dims','inquire_dimension YC '//trim(ic_filename))

! Latitudes : get dimid for 'YE' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'YE', dimid), &
            'get_grid_dims','inq_dimid YE '//trim(ic_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NYE), &
            'get_grid_dims','inquire_dimension YE '//trim(ic_filename))

! Vertical Levels : get dimid for 'ZC' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'ZC', dimid), &
            'get_grid_dims','inq_dimid ZC '//trim(ic_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NZC), &
            'get_grid_dims','inquire_dimension ZC '//trim(ic_filename))

! Vertical Levels : get dimid for 'ZE' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'ZE', dimid), &
            'get_grid_dims','inq_dimid ZE '//trim(ic_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=NZE), &
            'get_grid_dims','inquire_dimension ZE '//trim(ic_filename))

! tidy up

call nc_check(nf90_close(grid_id), &
         'get_grid_dims','close '//trim(ic_filename) )

end subroutine get_grid_dims



subroutine get_grid(NXC, NXE, NYC, NYE, NZC, NZE, &
                    ULAT, ULON, VLAT, VLON, WLAT, WLON, ZC, ZE)
!------------------------------------------------------------------
!
! Read the grid dimensions from the restart netcdf file.
!
! The file name comes from module storage ... namelist.

integer, intent(in) :: NXC   ! Number of Longitude centers
integer, intent(in) :: NXE   ! Number of Longitude edges
integer, intent(in) :: NYC   ! Number of Latitude  centers
integer, intent(in) :: NYE   ! Number of Latitude  edges
integer, intent(in) :: NZC   ! Number of Vertical grid centers
integer, intent(in) :: NZE   ! Number of Vertical grid edges

real(r8), dimension(:,:), intent(out) :: ULAT, ULON, VLAT, VLON, WLAT, WLON
real(r8), dimension( : ), intent(out) :: ZC, ZE

! type(grid_type), intent(out) :: Ugrid  ! (ZC, YC, XE)
! type(grid_type), intent(out) :: Vgrid  ! (ZC, YE, XC)
! type(grid_type), intent(out) :: Wgrid  ! (ZE, YC, XC)
! type(grid_type), intent(out) ::  grid  ! (ZC, YC, XC)

real(r8), dimension(NXC) :: XC
real(r8), dimension(NXE) :: XE
real(r8), dimension(NYC) :: YC
real(r8), dimension(NYE) :: YE

! real(r8), dimension(nx,ny), intent(out) :: ULAT, ULON, TLAT, TLON

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
integer                               :: VarID, numdims, dimlen
integer                               :: ncid, dimid


if ( .not. module_initialized ) call initialize_module

! get the ball rolling ...

call nc_check(nf90_open(trim(ic_filename), nf90_nowrite, ncid), 'get_grid', 'open '//trim(ic_filename))


! fixme - in a perfect world - 
! Get the variable ID
! Check to make sure it is the right shape
! Read it
call nc_check(nf90_inq_varid(ncid, 'XC', VarID), 'get_grid', 'inq_varid XC '//trim(ic_filename))
!call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
!                        'get_grid', 'inquire_variable XC '//trim(ic_filename))
call nc_check(nf90_get_var(ncid, VarID, XC), 'get_grid', 'get_var XC '//trim(ic_filename))

call nc_check(nf90_inq_varid(ncid, 'XE', VarID), 'get_grid', 'inq_varid XE '//trim(ic_filename))
!call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
!                        'get_grid', 'inquire_variable XE '//trim(ic_filename))
call nc_check(nf90_get_var(ncid, VarID, XE), 'get_grid', 'get_var XE '//trim(ic_filename))

call nc_check(nf90_inq_varid(ncid, 'YC', VarID), 'get_grid', 'inq_varid YC '//trim(ic_filename))
!call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
!                        'get_grid', 'inquire_variable YC '//trim(ic_filename))
call nc_check(nf90_get_var(ncid, VarID, YC), 'get_grid', 'get_var YC '//trim(ic_filename))

call nc_check(nf90_inq_varid(ncid, 'YE', VarID), 'get_grid', 'inq_varid YE '//trim(ic_filename))
!call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
!                        'get_grid', 'inquire_variable YE '//trim(ic_filename))
call nc_check(nf90_get_var(ncid, VarID, YE), 'get_grid', 'get_var YE '//trim(ic_filename))

call nc_check(nf90_inq_varid(ncid, 'ZC', VarID), 'get_grid', 'inq_varid ZC '//trim(ic_filename))
!call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
!                        'get_grid', 'inquire_variable ZC '//trim(ic_filename))
call nc_check(nf90_get_var(ncid, VarID, ZC), 'get_grid', 'get_var ZC '//trim(ic_filename))

call nc_check(nf90_inq_varid(ncid, 'ZE', VarID), 'get_grid', 'inq_varid ZE '//trim(ic_filename))
!call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
!                        'get_grid', 'inquire_variable ZE '//trim(ic_filename))
call nc_check(nf90_get_var(ncid, VarID, ZE), 'get_grid', 'get_var ZE '//trim(ic_filename))

! FIXME - need to convert these things to radians to store in the grid variables.
! FIXME - need to allocate the pointers in the grid variables, etc.

! call calc_tpoints(nx, ny, ULAT, ULON, TLAT, TLON)

! convert from radians to degrees

!ULAT = ULAT * rad2deg
!ULON = ULON * rad2deg
!TLAT = TLAT * rad2deg
!TLON = TLON * rad2deg

! ensure [0,360) [-90,90]

! where (ULON <   0.0_r8) ULON = ULON + 360.0_r8
! where (ULON > 360.0_r8) ULON = ULON - 360.0_r8
! where (TLON <   0.0_r8) TLON = TLON + 360.0_r8
! where (TLON > 360.0_r8) TLON = TLON - 360.0_r8
!
! where (ULAT < -90.0_r8) ULAT = -90.0_r8
! where (ULAT >  90.0_r8) ULAT =  90.0_r8
! where (TLAT < -90.0_r8) TLAT = -90.0_r8
! where (TLAT >  90.0_r8) TLAT =  90.0_r8

! tidy up

call nc_check(nf90_close(ncid), 'get_grid','close '//trim(ic_filename) )

end subroutine get_grid



function get_base_time_ncid( ncid )
!------------------------------------------------------------------
! The restart netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

type(time_type) :: get_base_time_ncid

integer, intent(in) :: ncid

integer :: year, month, day, hour, minute, second

if ( .not. module_initialized ) call initialize_module

call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'YEAR'  , year), &
                  'get_base_time', 'get_att year')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MONTH' , month), &
                  'get_base_time', 'get_att month')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'DAY'   , day), &
                  'get_base_time', 'get_att day')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'HOUR'  , hour), &
                  'get_base_time', 'get_att hour')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MINUTE', minute), &
                  'get_base_time', 'get_att minute')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'SECOND', second), &
                  'get_base_time', 'get_att second')

get_base_time_ncid = set_date(year, month, day, hour, minute, second)

end function get_base_time_ncid



function get_base_time_fname(filename)
!------------------------------------------------------------------
! The restart netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

type(time_type) :: get_base_time_fname

character(len=*), intent(in) :: filename

integer :: ncid, year, month, day, hour, minute, second

if ( .not. module_initialized ) call initialize_module

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'get_base_time',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'get_base_time', 'open '//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'YEAR'  , year), &
                  'get_base_time', 'get_att year')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MONTH' , month), &
                  'get_base_time', 'get_att month')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'DAY'   , day), &
                  'get_base_time', 'get_att day')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'HOUR'  , hour), &
                  'get_base_time', 'get_att hour')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MINUTE', minute), &
                  'get_base_time', 'get_att minute')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'SECOND', second), &
                  'get_base_time', 'get_att second')
call nc_check(nf90_close(ncid), 'get_base_time', 'close '//trim(filename))

get_base_time_fname = set_date(year, month, day, hour, minute, second)

end function get_base_time_fname



function get_state_time_ncid( ncid, filename )
!------------------------------------------------------------------
! the initialize_module ensures that the ncommas namelists are read.
! The restart times in the ncommas_in&restart_nml are used to define
! appropriate assimilation timesteps.
!
type(time_type)              :: get_state_time_ncid
integer,          intent(in) :: ncid
character(len=*), intent(in) :: filename

integer         :: VarID, numdims, dimlen
type(time_type) :: model_offset, base_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
integer, allocatable, dimension(:)    :: mytimes

if ( .not. module_initialized ) call initialize_module

base_time = get_base_time(ncid)

call nc_check( nf90_inq_varid(ncid, 'TIME', VarID), &
                  'get_state_time', 'inq_varid TIME '//trim(filename))

call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
                  'get_state_time', 'inquire TIME '//trim(filename))

if ( numdims > 1 ) then
   write(string1,*) 'TIME is not expected to have ',numdims,' dimensions.'
   call error_handler(E_ERR,'get_state_time',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlen), &
            'get_state_time', 'inquire time dimension length '//trim(filename))

allocate(mytimes(dimlen))

call nc_check( nf90_get_var(ncid, VarID, mytimes ), &
                  'get_state_time', 'get_var TIME '//trim(filename))

write(*,*)' temporal offset is (in seconds) is ',maxval(mytimes)
model_offset = set_time(maxval(mytimes))

get_state_time_ncid = base_time + model_offset

if (do_output()) &
    call print_time(get_state_time_ncid,'time for restart file '//trim(filename))
if (do_output()) &
    call print_date(get_state_time_ncid,'date for restart file '//trim(filename))

deallocate(mytimes)

end function get_state_time_ncid



function get_state_time_fname(filename)
!------------------------------------------------------------------
! the initialize_module ensures that the ncommas namelists are read.
! The restart times in the ncommas_in&restart_nml are used to define
! appropriate assimilation timesteps.
!
type(time_type) :: get_state_time_fname
character(len=*), intent(in) :: filename

integer         :: ncid, VarID, numdims, dimlen
type(time_type) :: model_offset, base_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
integer, allocatable, dimension(:)    :: mytimes

if ( .not. module_initialized ) call initialize_module

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'get_state_time',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'get_base_time', 'open '//trim(filename))

base_time = get_base_time(ncid)

call nc_check( nf90_inq_varid(ncid, 'TIME', VarID), &
                  'get_state_time', 'inq_varid TIME '//trim(filename))

call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims), &
                  'get_state_time', 'inquire TIME '//trim(filename))

if ( numdims > 1 ) then
   write(string1,*) 'TIME is not expected to have ',numdims,' dimensions.'
   call error_handler(E_ERR,'get_state_time',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlen), &
            'get_state_time', 'inquire time dimension length '//trim(filename))

allocate(mytimes(dimlen))

call nc_check( nf90_get_var(ncid, VarID, mytimes ), &
                  'get_state_time', 'get_var TIME '//trim(filename))
call nc_check(nf90_close(ncid), 'get_state_time', 'close '//trim(filename))

write(*,*)' temporal offset is (in seconds) is ',maxval(mytimes)
model_offset = set_time(maxval(mytimes))

get_state_time_fname = base_time + model_offset

if (do_output()) &
    call print_time(get_state_time_fname,'time for restart file '//trim(filename))
if (do_output()) &
    call print_date(get_state_time_fname,'date for restart file '//trim(filename))

deallocate(mytimes)

end function get_state_time_fname



function set_model_time_step()
!------------------------------------------------------------------
! the initialize_module ensures that the ncommas namelists are read.
! The restart times in the ncommas_in&restart_nml are used to define
! appropriate assimilation timesteps.
!
type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call initialize_module

! FIXME - determine when we can stop the model

   set_model_time_step = set_time(0, 1) ! (seconds, days)

end function set_model_time_step




subroutine write_ncommas_namelist(model_time, adv_to_time)
!------------------------------------------------------------------
!
type(time_type), INTENT(IN) :: model_time, adv_to_time
type(time_type) :: offset

integer :: iunit, secs, days

if ( .not. module_initialized ) call initialize_module

offset = adv_to_time - model_time
call get_time(offset, secs, days)

if (secs /= 0 ) then
   write(string1,*)'adv_to_time has seconds == ',secs,' must be zero'
   call error_handler(E_ERR,'write_ncommas_namelist', string1, source, revision, revdate)
endif

! call print_date( model_time,'write_ncommas_namelist:dart model date')
! call print_date(adv_to_time,'write_ncommas_namelist:advance_to date')
! call print_time( model_time,'write_ncommas_namelist:dart model time')
! call print_time(adv_to_time,'write_ncommas_namelist:advance_to time')
! call print_time(     offset,'write_ncommas_namelist:a distance of')
! write( *,'(''write_ncommas_namelist:TIME_MANAGER_NML   STOP_COUNT '',i10,'' days'')') days

!Convey the information to the namelist 'stop option' and 'stop count'

!if ( trim(stop_option) == 'nday' ) then
!   stop_count = days
!else
   call error_handler(E_ERR,'write_ncommas_namelist', &
              'stop_option must be "nday"', source, revision, revdate)
!endif

iunit = open_file('ncommas_in.DART',form='formatted',action='rewind')
write(iunit, nml=restart_nml)
write(iunit, '('' '')')
close(iunit)

end subroutine write_ncommas_namelist





  subroutine calc_tpoints(nx, ny, ULAT, ULON, TLAT, TLON)
!------------------------------------------------------------------
! subroutine calc_tpoints(nx, ny, ULAT, ULON, TLAT, TLON)
!
! mimic ncommas grid.F90:calc_tpoints(), but for one big block.

integer,                    intent( in) :: nx, ny
real(r8), dimension(nx,ny), intent( in) :: ULAT, ULON
real(r8), dimension(nx,ny), intent(out) :: TLAT, TLON

integer  :: i, j
real(r8) :: xc,yc,zc,xs,ys,zs,xw,yw,zw   ! Cartesian coordinates for
real(r8) :: xsw,ysw,zsw,tx,ty,tz,da      ! nbr points

real(r8), parameter ::  c0 = 0.000_r8, c1 = 1.000_r8
real(r8), parameter ::  c2 = 2.000_r8, c4 = 4.000_r8
real(r8), parameter :: p25 = 0.250_r8, p5 = 0.500_r8
real(r8)            :: pi, pi2, pih, radian

if ( .not. module_initialized ) call initialize_module

! Define some constants as in ncommas

pi     = c4*atan(c1)
pi2    = c2*pi
pih    = p5*pi
radian = 180.0_r8/pi

do j=2,ny
do i=2,nx

   !*** convert neighbor U-cell coordinates to 3-d Cartesian coordinates 
   !*** to prevent problems with averaging near the pole

   zsw = cos(ULAT(i-1,j-1))
   xsw = cos(ULON(i-1,j-1))*zsw
   ysw = sin(ULON(i-1,j-1))*zsw
   zsw = sin(ULAT(i-1,j-1))

   zs  = cos(ULAT(i  ,j-1))
   xs  = cos(ULON(i  ,j-1))*zs
   ys  = sin(ULON(i  ,j-1))*zs
   zs  = sin(ULAT(i  ,j-1))

   zw  = cos(ULAT(i-1,j  ))
   xw  = cos(ULON(i-1,j  ))*zw
   yw  = sin(ULON(i-1,j  ))*zw
   zw  = sin(ULAT(i-1,j  ))

   zc  = cos(ULAT(i  ,j  ))
   xc  = cos(ULON(i  ,j  ))*zc
   yc  = sin(ULON(i  ,j  ))*zc
   zc  = sin(ULAT(i  ,j  ))

   !*** straight 4-point average to T-cell Cartesian coords

   tx = p25*(xc + xs + xw + xsw)
   ty = p25*(yc + ys + yw + ysw)
   tz = p25*(zc + zs + zw + zsw)

   !*** convert to lat/lon in radians

   da = sqrt(tx**2 + ty**2 + tz**2)

   TLAT(i,j) = asin(tz/da)

   if (tx /= c0 .or. ty /= c0) then
      TLON(i,j) = atan2(ty,tx)
   else
      TLON(i,j) = c0
   endif

end do
end do

!*** for bottom row of domain where sw 4pt average is not valid,
!*** extrapolate from interior
!*** NOTE: THIS ASSUMES A CLOSED SOUTH BOUNDARY - WILL NOT
!***       WORK CORRECTLY FOR CYCLIC OPTION

do i=1,nx
   TLON(i,1) =    TLON(i,1+1)
   TLAT(i,1) = c2*TLAT(i,1+1) - TLAT(i,1+2)
end do

where (TLON(:,:) > pi2) TLON(:,:) = TLON(:,:) - pi2
where (TLON(:,:) < c0 ) TLON(:,:) = TLON(:,:) + pi2

!*** this leaves the leftmost/western edge to be filled 
!*** if the longitudes wrap, this is easy.
!*** the gx3v5 grid TLON(:,2) and TLON(:,nx) are both about 2pi,
!*** so taking the average is reasonable.
!*** averaging the latitudes is always reasonable.

if ( trim(ew_boundary_type) == 'cyclic' ) then

   TLAT(1,:) = (TLAT(2,:) + TLAT(nx,:))/c2
   TLON(1,:) = (TLON(2,:) + TLON(nx,:))/c2

else
   write(string1,'(''ncommas_in&domain_nml:ew_boundary_type '',a,'' unknown.'')') &
                                    trim(ew_boundary_type)
   call error_handler(E_ERR,'calc_tpoints',string1,source,revision,revdate)
endif

end subroutine calc_tpoints


!------------------------------------------------------------------


subroutine get_ncommas_restart_filename( filename )
character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call initialize_module

filename   = trim(ic_filename)

end subroutine get_ncommas_restart_filename




!From netCDF file we need the following variables:
!
!xg_pos(1), yg_pos(1), lat (variable), lon(variable)
!
!
!xctrue(:) = xc(:) + xg_pos(1)
!yctrue(:) = yc(:) + yg_pos(1)
!
!xetrue(:) = xe(:) + xg_pos(1)
!yetrue(:) = ye(:) + yg_pos(1)
!
!
!DO j = 1,ny-1
! DO i = 1,nx-1
!    CALL XY_TO_LL(new_lat, new_lon, 0, xctrue(i), yctrue(j), lat, lon)
!    slat(i,j) = new_lat
!    slon(i,j) = new_lon 
! ENDDO
!ENDDO
!
!DO j = 1,ny-1
! DO i = 1,nx
!    CALL XY_TO_LL(new_lat, new_lon, 0, xetrue(i), yctrue(j), lat, lon)
!    ulat(i,j) = new_lat
!    ulon(i,j) = new_lon 
! ENDDO
!ENDDO
!
!DO j = 1,ny
! DO i = 1,nx-1
!    CALL XY_TO_LL(new_lat, new_lon, 0, xctrue(i), yetrue(j), lat, lon)
!    vlat(i,j) = new_lat
!    vlon(i,j) = new_lon 
! ENDDO
!ENDDO
!
!
!!############################################################################
!!
!!     ##################################################################
!!     ######                                                                                                                              ######
!!     ######                                     SUBROUTINE XY_TO_LL                                         ######
!!     ######                                                                                                                              ######
!!     ##################################################################
!!
!!
!!     PURPOSE:
!!
!!     This subroutine computes the projected (lat, lon) coordinates of the
!!     point (x, y) relative to (lat0, lon0).  Various map projections
!!     are possible.
!!
!!############################################################################
!!
!!     Author:  David Dowell
!!
!!     Creation Date:  25 February 2005
!!     Modified:  12 April 2005 (changed units of rearth, x, and y from km to m)
!!
!!############################################################################
!
! SUBROUTINE XY_TO_LL(lat, lon, map_proj, x, y, lat0, lon0)
!
! implicit none
!
!! Passed variables
!
!   integer map_proj            ! map projection:
!                               !   0 = flat earth
!                               !   1 = oblique azimuthal
!                               !   2 = Lambert conformal
!
!   real x, y                   ! distance (m)
!   real lat0, lon0             ! coordinates (rad) of origin (where x=0, y=0)
!
!! Returned variables
!
!   real lat, lon               ! coordinates (rad) of point
!
!! Local variables
!
!   real rearth; parameter(rearth=1000.0 * 6367.0)      ! radius of earth (m)
!
!   if (map_proj.eq.0) then
!     lat = lat0 + y / rearth
!     lon = lon0 + x / ( rearth * cos(0.5*(lat0+lat)) )
!   else
!     write(*,*) 'map projection unavailable:  ', map_proj
!     stop
!   endif
!
!   RETURN
!   END



end module dart_ncommas_mod

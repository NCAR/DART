! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module dart_cice_mod

use        types_mod, only : r8, rad2deg, PI, SECPERDAY, digits12
use time_manager_mod, only : time_type, get_date, set_date, get_time, set_time, &
                             set_calendar_type, get_calendar_string, &
                             print_date, print_time, operator(==), operator(-)
use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler, &
                             find_namelist_in_file, check_namelist_read, &
                             E_ERR, E_MSG, find_textfile_dims

use  netcdf_utilities_mod, only : nc_check
                             

use typesizes
use netcdf

implicit none
private

public :: get_cice_calendar, set_model_time_step, &
          get_horiz_grid_dims, get_ncat_dim,      &
          read_horiz_grid, read_topography,       &
          write_cice_namelist, get_cice_restart_filename, &
          set_binary_file_conversion

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

character(len=512) :: msgstring
logical, save :: module_initialized = .false.

character(len=256) :: ic_filename      = 'cice.r.nc'

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.
logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

! if the binary grid files are big_endian and you are running on
! a little_endian machine, set this string so the convert will
! happen correctly.  options are:  'native', 'big_endian', 'little_endian'

character(len=64) :: conversion = 'native'    ! no conversion


!------------------------------------------------------------------
! The CICE grid info namelist 
!------------------------------------------------------------------
!
! CICE grid information comes in several files:
!   horizontal grid lat/lons in one, 
!   topography (used to get land/ocean cells) in another
!
!------------------------------------------------------------------
!
! Here is what we can get from the (binary) horiz grid file:
!   real (r8), dimension(:,:), allocatable :: &
!      ULAT,            &! latitude  (radians) of U points
!      ULON,            &! longitude (radians) of U points
!      HTN ,            &! length (cm) of north edge of T box
!      HTE ,            &! length (cm) of east  edge of T box
!      HUS ,            &! length (cm) of south edge of U box
!      HUW ,            &! length (cm) of west  edge of U box
!      ANGLE             ! angle
!
! Here is what we can get from the topography file:
!   integer, dimension(:,:), allocatable :: &
!      KMT               ! k index of deepest grid cell on T grid
!
! These must be derived or come from someplace else ...
!      KMU               ! k index of deepest grid cell on U grid
!      HT                ! real(r8) value of deepest valid T depth (in cm)
!      HU                ! real(r8) value of deepest valid U depth (in cm)
!
!------------------------------------------------------------------


character(len=256) :: grid_file, grid_format,  grid_type, &     ! CMB verified types
                      gridcpl_file, kmt_file
integer :: kcatbound

! CMB stole code from ice_domain.F90
   integer, parameter :: log_kind  = kind(.true.), &
                         int_kind  = selected_int_kind(6)

   character (len=80) :: &
      ew_boundary_type, &! type of domain bndy in each logical
      ns_boundary_type   !    direction (ew is i, ns is j)

   logical (kind=log_kind) :: &
      maskhalo_dyn   , & ! if true, use masked halo updates for dynamics
      maskhalo_remap , & ! if true, use masked halo updates for transport
      maskhalo_bound     ! if true, use masked halo updates for bound_state

   character (len=80) :: &
       distribution_type,   &! method to use for distributing blocks
       distribution_wght     ! method for weighting work per block 

   integer (int_kind) :: nprocs                ! num of processors

   character (len=80) :: processor_shape


! code from cesm1_5_beta06c/cime/share/csm_share/shr/shr_kind_mod.F990
   integer,parameter :: SHR_QTY_CS = 80                     ! short char
   integer,parameter :: SHR_QTY_IN = kind(1)                ! native integer
   
! code from cesm1_5_beta06c/cime/driver_cpl/shr/seq_timemgr_mod.F90
    character(SHR_QTY_CS)  :: calendar              ! Calendar type
    character(SHR_QTY_CS)  :: stop_option           ! Stop option units
    integer(SHR_QTY_IN)    :: stop_n                ! Number until stop
    integer(SHR_QTY_IN)    :: stop_ymd              ! Stop date (YYYYMMDD)
    integer(SHR_QTY_IN)    :: stop_tod              ! Stop time-of-day
    character(SHR_QTY_CS)  :: restart_option        ! Restart option units
    integer(SHR_QTY_IN)    :: restart_n             ! Number until restart interval
    integer(SHR_QTY_IN)    :: restart_ymd           ! Restart date (YYYYMMDD)
    character(SHR_QTY_CS)  :: pause_option
    integer(SHR_QTY_IN)    :: pause_n
    character(SHR_QTY_CS)  :: pause_component_list
    character(SHR_QTY_CS)  :: history_option        ! History option units
    integer(SHR_QTY_IN)    :: history_n             ! Number until history interval
    integer(SHR_QTY_IN)    :: history_ymd           ! History date (YYYYMMDD)
    character(SHR_QTY_CS)  :: histavg_option        ! Histavg option units
    integer(SHR_QTY_IN)    :: histavg_n             ! Number until histavg interval
    integer(SHR_QTY_IN)    :: histavg_ymd           ! Histavg date (YYYYMMDD)
    character(SHR_QTY_CS)  :: barrier_option        ! Barrier option units
    integer(SHR_QTY_IN)    :: barrier_n             ! Number until barrier interval
    integer(SHR_QTY_IN)    :: barrier_ymd           ! Barrier date (YYYYMMDD)
    character(SHR_QTY_CS)  :: tprof_option          ! tprof option units
    integer(SHR_QTY_IN)    :: tprof_n               ! Number until tprof interval
    integer(SHR_QTY_IN)    :: tprof_ymd             ! tprof date (YYYYMMDD)
    integer(SHR_QTY_IN)    :: start_ymd             ! Start date (YYYYMMDD)
    integer(SHR_QTY_IN)    :: start_tod             ! Start time of day (seconds)
    integer(SHR_QTY_IN)    :: curr_ymd              ! Current ymd (YYYYMMDD)
    integer(SHR_QTY_IN)    :: curr_tod              ! Current tod (seconds)
    integer(SHR_QTY_IN)    :: ref_ymd               ! Reference date (YYYYMMDD)
    integer(SHR_QTY_IN)    :: ref_tod               ! Reference time of day (seconds)
    integer(SHR_QTY_IN)    :: atm_cpl_dt            ! Atmosphere coupling interval
    integer(SHR_QTY_IN)    :: lnd_cpl_dt            ! Land coupling interval
    integer(SHR_QTY_IN)    :: ice_cpl_dt            ! Sea-Ice coupling interval
    integer(SHR_QTY_IN)    :: ocn_cpl_dt            ! Ocean coupling interval
    integer(SHR_QTY_IN)    :: glc_cpl_dt            ! Glc coupling interval
    character(SHR_QTY_CS)  :: glc_avg_period
    integer(SHR_QTY_IN)    :: rof_cpl_dt            ! Runoff coupling interval
    integer(SHR_QTY_IN)    :: wav_cpl_dt            ! Wav coupling interval
    integer(SHR_QTY_IN)    :: esp_cpl_dt            ! Esp coupling interval
    integer(SHR_QTY_IN)    :: atm_cpl_offset        ! Atmosphere coupling interval
    integer(SHR_QTY_IN)    :: lnd_cpl_offset        ! Land coupling interval
    integer(SHR_QTY_IN)    :: ice_cpl_offset        ! Sea-Ice coupling interval
    integer(SHR_QTY_IN)    :: ocn_cpl_offset        ! Ocean coupling interval
    integer(SHR_QTY_IN)    :: glc_cpl_offset        ! Glc coupling interval
    integer(SHR_QTY_IN)    :: wav_cpl_offset        ! Wav coupling interval
    integer(SHR_QTY_IN)    :: rof_cpl_offset        ! Runoff coupling interval
    integer(SHR_QTY_IN)    :: esp_cpl_offset        ! Esp coupling interval
    logical                 :: end_restart           ! Write restart at end of run
    logical                :: esp_run_on_pause

   namelist /grid_nml/ grid_file, grid_format, grid_type, &
         gridcpl_file, kcatbound, kmt_file

   namelist /domain_nml/ nprocs, &
                         processor_shape,   &
                         distribution_type, &
                         distribution_wght, &
                         ew_boundary_type,  &
                         ns_boundary_type,  &
                         maskhalo_dyn,      &
                         maskhalo_remap,    &
                         maskhalo_bound

   namelist /seq_timemgr_inparm/  calendar, curr_ymd, curr_tod,  &
         stop_option, stop_n, stop_ymd, stop_tod,        &
         restart_option, restart_n, restart_ymd,         &
         history_option, history_n, history_ymd,         &
         histavg_option, histavg_n, histavg_ymd,         &
         pause_option, pause_n, pause_component_list,   &
         barrier_option, barrier_n, barrier_ymd,         &
         tprof_option, tprof_n, tprof_ymd,               &
         start_ymd, start_tod, ref_ymd, ref_tod,         &
         atm_cpl_dt, ocn_cpl_dt, ice_cpl_dt, lnd_cpl_dt, &
         atm_cpl_offset, lnd_cpl_offset, ocn_cpl_offset, &
         ice_cpl_offset, glc_cpl_dt, glc_avg_period, glc_cpl_offset,     &
         wav_cpl_dt, wav_cpl_offset, esp_cpl_dt, esp_cpl_offset,     &
         esp_run_on_pause, rof_cpl_dt, rof_cpl_offset, end_restart

!======================================================================
contains
!======================================================================


subroutine initialize_module

integer :: iunit, io

! Read calendar information 
! In 'restart' mode, this is primarily the calendar type and 'stop'
! information. The time attributes of the restart file override 
! the namelist time information. 

! FIXME : Real observations are always GREGORIAN dates ...
! but stomping on that here gets in the way of running
! a perfect_model experiment for pre-1601 AD cases.

! STOMP if ( allow_leapyear ) then
   call set_calendar_type('gregorian')
! STOMP else
! STOMP    call set_calendar_type('noleap')
! STOMP endif

! Make sure we have a cice restart file (for grid dims)
if ( .not. file_exist(ic_filename) ) then
   msgstring = 'dart_cice_mod: '//trim(ic_filename)//' not found'
   call error_handler(E_ERR,'initialize_module', &
          msgstring, source, revision, revdate)
endif

! Read CICE grid information (for grid filenames)
call find_namelist_in_file('cice_in', 'grid_nml', iunit)
read(iunit, nml = grid_nml, iostat = io)
call check_namelist_read(iunit, io, 'grid_nml')

! Read CICE grid information (for grid filenames)
call find_namelist_in_file('cice_in', 'domain_nml', iunit)
read(iunit, nml = domain_nml, iostat = io)
call check_namelist_read(iunit, io, 'domain_nml')

! Read CESM "driver" restart information 
call find_namelist_in_file('drv_in', 'seq_timemgr_inparm', iunit)
read(iunit, nml = seq_timemgr_inparm, iostat = io)
call check_namelist_read(iunit, io, 'seq_timemgr_inparm')

module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

end subroutine initialize_module

!------------------------------------------------------------------

subroutine get_ncat_dim(Ncat)

!
! Read the ncat size from the restart netcdf file.

integer, intent(out) :: Ncat   ! Number of categories in ice-thick dist

integer :: grid_id, dimid, nc_rc

if ( .not. module_initialized ) call initialize_module

call nc_check(nf90_open(trim(ic_filename), nf90_nowrite, grid_id), &
         'get_ncat_dim','open '//trim(ic_filename))

! ncat : get dimid for 'ncat' and then get value
nc_rc = nf90_inq_dimid(grid_id, 'ncat', dimid)
if (nc_rc /= nf90_noerr) then
   nc_rc = nf90_inq_dimid(grid_id, 'Ncat', dimid)
   if (nc_rc /= nf90_noerr) then
      msgstring = "unable to find either 'ncat' or 'Ncat' in file "//trim(ic_filename)
      call error_handler(E_ERR, 'get_horiz_grid_dims', msgstring, &
                         source,revision,revdate) 
   endif
endif

call nc_check(nf90_inquire_dimension(grid_id, dimid, len=Ncat), &
         'get_ncat_dim','inquire_dimension ni '//trim(ic_filename))

! tidy up

call nc_check(nf90_close(grid_id), &
         'get_ncat_dim','close '//trim(ic_filename) )

end subroutine get_ncat_dim

!------------------------------------------------------------------

subroutine get_horiz_grid_dims(Nx, Ny)

!
! Read the lon, lat grid size from the restart netcdf file.
! The actual grid file is a binary file with no header information.
!
! The file name comes from module storage ... namelist.

integer, intent(out) :: Nx   ! Number of Longitudes
integer, intent(out) :: Ny   ! Number of Latitudes

integer :: grid_id, dimid, nc_rc

if ( .not. module_initialized ) call initialize_module

call nc_check(nf90_open(trim(ic_filename), nf90_nowrite, grid_id), &
         'get_horiz_grid_dims','open '//trim(ic_filename))

! Longitudes : get dimid for 'ni' or 'nlon', and then get value
nc_rc = nf90_inq_dimid(grid_id, 'ni', dimid)
if (nc_rc /= nf90_noerr) then
   nc_rc = nf90_inq_dimid(grid_id, 'nlon', dimid)
   if (nc_rc /= nf90_noerr) then
      msgstring = "unable to find either 'ni' or 'nlon' in file "//trim(ic_filename)
      call error_handler(E_ERR, 'get_horiz_grid_dims', msgstring, &
                         source,revision,revdate) 
   endif
endif

call nc_check(nf90_inquire_dimension(grid_id, dimid, len=Nx), &
         'get_horiz_grid_dims','inquire_dimension ni '//trim(ic_filename))

! Latitudes : get dimid for 'nj' or 'nlat' ... and then get value
nc_rc = nf90_inq_dimid(grid_id, 'nj', dimid)
if (nc_rc /= nf90_noerr) then
   nc_rc = nf90_inq_dimid(grid_id, 'nlat', dimid)
   if (nc_rc /= nf90_noerr) then
      msgstring = "unable to find either 'nj' or 'nlat' in "//trim(ic_filename)
      call error_handler(E_ERR, 'get_horiz_grid_dims', msgstring, &
                         source,revision,revdate)
   endif
endif

call nc_check(nf90_inquire_dimension(grid_id, dimid, len=Ny), &
         'get_horiz_grid_dims','inquire_dimension ni '//trim(ic_filename))

! tidy up

call nc_check(nf90_close(grid_id), &
         'get_horiz_grid_dims','close '//trim(ic_filename) )

end subroutine get_horiz_grid_dims

!------------------------------------------------------------------

! CMB removed get_vert_grid_dim(Nz)

!------------------------------------------------------------------
   
subroutine get_cice_calendar(calstring)

! the initialize_module ensures that the cice namelists are read and 
! the DART time manager gets the cice calendar setting.
!
! Then, the DART time manager is queried to return what it knows ...
!
character(len=*), INTENT(OUT) :: calstring

if ( .not. module_initialized ) call initialize_module

call get_calendar_string(calstring)

end subroutine get_cice_calendar

!------------------------------------------------------------------

function set_model_time_step()

! the initialize_module ensures that the cice namelists are read.
! The restart times in the cice_in&restart_nml are used to define
! appropriate assimilation timesteps.
!
type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call initialize_module

! Check the 'restart_option' and 'restart_n' to determine
! when we can stop the model
! CMB not sure if nday is actually different than ndays, no matter here though
!if ( (trim(restart_option) == 'ndays') .or. (trim(restart_option) == 'nday' ) ) then
!   set_model_time_step = set_time(0, restart_n) ! (seconds, days)
!else if ( trim(restart_option) == 'nyears' ) then
   ! FIXME ... CMB I guess we ignore it and make the freq 1 day anyway?
   set_model_time_step = set_time(0, 1) ! (seconds, days)
!else
!   call error_handler(E_ERR,'set_model_time_step', &
!              'restart_option must be ndays or nday', source, revision, revdate)
!endif

end function set_model_time_step

!------------------------------------------------------------------

subroutine write_cice_namelist(model_time, adv_to_time)
!
type(time_type), INTENT(IN) :: model_time, adv_to_time
type(time_type) :: offset

integer :: iunit, secs, days

if ( .not. module_initialized ) call initialize_module

offset = adv_to_time - model_time
call get_time(offset, secs, days)

if (secs /= 0 ) then
   write(msgstring,*)'adv_to_time has seconds == ',secs,' must be zero'
   call error_handler(E_ERR,'write_cice_namelist', msgstring, source, revision, revdate)
endif

! call print_date( model_time,'write_cice_namelist:dart model date')
! call print_date(adv_to_time,'write_cice_namelist:advance_to date')
! call print_time( model_time,'write_cice_namelist:dart model time')
! call print_time(adv_to_time,'write_cice_namelist:advance_to time')
! call print_time(     offset,'write_cice_namelist:a distance of')
! write( *,'(''write_cice_namelist:TIME_MANAGER_NML   STOP_COUNT '',i10,'' days'')') days

!Convey the information to the namelist 'stop option' and 'stop count'

if ( (trim(stop_option) == 'nday') .or. (trim(stop_option) == 'ndays') ) then
   stop_n = days
else
   call error_handler(E_ERR,'write_cice_namelist', &
              'stop_option must be "ndays or nday"', source, revision, revdate)
endif

iunit = open_file('cice_in.DART',form='formatted',action='rewind')
write(iunit, nml=seq_timemgr_inparm)
write(iunit, '('' '')')
close(iunit)

end subroutine write_cice_namelist

!------------------------------------------------------------------
!> Open and read the binary grid file


  subroutine read_horiz_grid(nx, ny, ULAT, ULON, TLAT, TLON)

integer,                    intent(in)  :: nx, ny
real(r8), dimension(nx,ny), intent(out) :: ULAT, ULON, TLAT, TLON

!real(r8), dimension(nx,ny) :: &
!     HTN ,  &! length (cm) of north edge of T box
!     HTE ,  &! length (cm) of east  edge of T box
!     HUS ,  &! length (cm) of south edge of U box
!     HUW ,  &! length (cm) of west  edge of U box
!     ANGLE  ! angle

integer :: grid_unit, reclength
real(digits12), dimension(nx,ny) :: ULAT64, ULON64

if ( .not. module_initialized ) call initialize_module

! Check to see that the file exists.

if ( .not. file_exist(grid_file) ) then
   msgstring = 'cice_in:grid_file '//trim(grid_file)//' not found'
   call error_handler(E_ERR,'read_horiz_grid', &
          msgstring, source, revision, revdate)
endif

! Open it and read them in the EXPECTED order.
! Actually, we only need the first two, so I'm skipping the rest.

grid_unit = get_unit()
INQUIRE(iolength=reclength) ULAT64

open(grid_unit, file=trim(grid_file), form='unformatted', convert=conversion, &
            access='direct', recl=reclength, status='old', action='read' )
read(grid_unit, rec=1) ULAT64
read(grid_unit, rec=2) ULON64
!read(grid_unit, rec=3) HTN
!read(grid_unit, rec=4) HTE
!read(grid_unit, rec=5) HUS
!read(grid_unit, rec=6) HUW
!read(grid_unit, rec=7) ANGLE
close(grid_unit)

ULAT = ULAT64
ULON = ULON64

call calc_tpoints(nx, ny, ULAT, ULON, TLAT, TLON)

! convert from radians to degrees

ULAT = ULAT * rad2deg
ULON = ULON * rad2deg
TLAT = TLAT * rad2deg
TLON = TLON * rad2deg

! ensure [0,360) [-90,90]

where (ULON <   0.0_r8) ULON = ULON + 360.0_r8
where (ULON > 360.0_r8) ULON = ULON - 360.0_r8
where (TLON <   0.0_r8) TLON = TLON + 360.0_r8
where (TLON > 360.0_r8) TLON = TLON - 360.0_r8

where (ULAT < -90.0_r8) ULAT = -90.0_r8
where (ULAT >  90.0_r8) ULAT =  90.0_r8
where (TLAT < -90.0_r8) TLAT = -90.0_r8
where (TLAT >  90.0_r8) TLAT =  90.0_r8

end subroutine read_horiz_grid


!------------------------------------------------------------------
!> mimic POP grid.F90:calc_tpoints(), but for one big block.


subroutine calc_tpoints(nx, ny, ULAT, ULON, TLAT, TLON)

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

! Define some constants as in cice

pi     = c4*atan(c1)
pi2    = c2*pi
pih    = p5*pi
radian = 180.0_r8/pi

! initialize these arrays to 0. in the code below there is
! a column that is referenced by a where() construct before
! those values are set.  make sure that it doesn't cause a
! floating point exception from random memory bits which aren't
! valid floating point numbers.
TLAT(:,:) = c0
TLON(:,:) = c0

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
   write(msgstring,'(''cice_in&domain_nml:ew_boundary_type '',a,'' unknown.'')') &
                                    trim(ew_boundary_type)
   call error_handler(E_ERR,'calc_tpoints',msgstring,source,revision,revdate)
endif

end subroutine calc_tpoints


!------------------------------------------------------------------
!> Open and read the binary topography file


subroutine read_topography(nx, ny, KMT, KMU)

integer,                   intent(in)  :: nx, ny
integer, dimension(nx,ny), intent(out) :: KMT, KMU

integer  :: i, j, topo_unit, reclength

if ( .not. module_initialized ) call initialize_module

! Check to see that the file exists.

if ( .not. file_exist(kmt_file) ) then
   msgstring = 'cice_in:kmt_file '//trim(kmt_file)//' not found'
   call error_handler(E_ERR,'read_topography', &
          msgstring, source, revision, revdate)
endif

! read the binary file

topo_unit = get_unit()
INQUIRE(iolength=reclength) KMT

open( topo_unit, file=trim(kmt_file), form='unformatted', &
            access='direct', recl=reclength, status='old', action='read', convert=conversion)
read( topo_unit, rec=1) KMT
close(topo_unit)

! the equation numbered 3.2 on page 15 of this document:
!  http://www.cesm.ucar.edu/models/cesm1.0/cice2/doc/sci/POPRefManual.pdf
! is WRONG.  (WRONG == inconsistent with the current POP source code.)
!
! for any U(i,j), the T(i,j) point with the same index values is located
! south and west. so the T points which surround any U(i,j) point are
! in fact at indices i,i+1, and j,j+1 .
!
!  NO: KMU(i,j) = min(KMT(i, j), KMT(i-1, j), KMT(i, j-1), KMT(i-1, j-1))
! YES: KMU(i,j) = min(KMT(i, j), KMT(i+1, j), KMT(i, j+1), KMT(i+1, j+1))
!
! the latter matches the POP source code, on yellowstone, lines 908 and 909 in:
!  /glade/p/cesm/releases/cesm1_2_2/models/ocn/pop2/source/grid.F90
!
! wrap around longitude boundary at i == nx.  set the topmost (last) latitude
! U row to the same value in all cases. in the shifted pole grid currently in
! use all these points are on land and so are 0.  in the original unshifted
! lat/lon grid these last row U points are above the final T row and are believed
! to be unused.  for completeness we set all values in the last U row to the
! minimum of the all T row values immediately below it, for all longitudes.

do j=1,ny-1
   do i=1,nx-1
      KMU(i,j) = min(KMT(i, j), KMT(i+1, j), KMT(i, j+1), KMT(i+1, j+1))
   enddo
   KMU(nx,j) = min(KMT(nx, j), KMT(1, j), KMT(nx, j+1), KMT(1, j+1))
enddo
KMU(:,ny) = minval(KMT(:,ny))

end subroutine read_topography


!------------------------------------------------------------------
! CMB removed read_vert_grid
!------------------------------------------------------------------

subroutine get_cice_restart_filename( filename )
character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call initialize_module

filename   = trim(ic_filename)

end subroutine get_cice_restart_filename


!------------------------------------------------------------------


subroutine set_binary_file_conversion(convertstring)

character(len=*), intent(in) :: convertstring

if ( .not. module_initialized ) call initialize_module

conversion = convertstring

end subroutine set_binary_file_conversion


!------------------------------------------------------------------


end module dart_cice_mod


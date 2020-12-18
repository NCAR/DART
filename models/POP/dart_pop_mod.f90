! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module dart_pop_mod

use        types_mod, only : r4, r8, rad2deg, PI, SECPERDAY, MISSING_R8, digits12
use time_manager_mod, only : time_type, get_date, set_date, get_time, set_time, &
                             set_calendar_type, get_calendar_string, &
                             print_date, print_time, operator(==), operator(-)
use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler, &
                             find_namelist_in_file, check_namelist_read, &
                             E_ERR, E_WARN, E_MSG, find_textfile_dims, &
                             logfileunit

use  netcdf_utilities_mod, only : nc_check

use typesizes
use netcdf

implicit none
private

public :: get_pop_calendar, set_model_time_step, &
          get_horiz_grid_dims, get_vert_grid_dim, &
          read_horiz_grid, read_topography, read_vert_grid, &
          write_pop_namelist, get_pop_restart_filename, &
          set_binary_file_conversion, read_mean_dynamic_topography

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2, string3

logical, save :: module_initialized = .false.

character(len=256) :: ic_filename      = 'pop.r.nc'
!character(len=256) :: restart_filename = 'dart_pop_mod_restart_filename_not_set'

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.
logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

! if the binary grid files are big_endian and you are running on
! a little_endian machine, set this string so the convert will
! happen correctly.  options are:  'native', 'big_endian', 'little_endian'

character(len=64) :: conversion = 'native'    ! no conversion

!------------------------------------------------------------------
! The POP time manager namelist variables
!------------------------------------------------------------------

character(len=100) :: accel_file ! length consistent with POP
character(len= 64) :: stop_option, runid, dt_option, time_mix_opt
character(len=  1) :: date_separator
logical  :: impcor, laccel, allow_leapyear
real(r8) :: dtuxcel, dt_count
integer  :: iyear0, imonth0, iday0, ihour0, iminute0, isecond0
integer  :: stop_count, fit_freq, time_mix_freq
real(r8) :: robert_alpha, robert_nu


namelist /time_manager_nml/ runid, time_mix_opt, time_mix_freq, &
    impcor, laccel, accel_file, dtuxcel, iyear0, imonth0, &
    iday0, ihour0, iminute0, isecond0, dt_option, dt_count, &
    stop_option, stop_count, date_separator, allow_leapyear, fit_freq, &
    robert_alpha, robert_nu


!------------------------------------------------------------------
! The POP I/O namelist variables
!------------------------------------------------------------------

character(len=100) :: log_filename, pointer_filename
logical :: lredirect_stdout, luse_pointer_files
logical :: luse_nf_64bit_offset
integer :: num_iotasks

namelist /io_nml/ num_iotasks, lredirect_stdout, log_filename, &
    luse_pointer_files, pointer_filename, luse_nf_64bit_offset

!------------------------------------------------------------------
! The POP restart manager namelist variables
!------------------------------------------------------------------

character(len=100) :: restart_outfile ! length consistent with POP
character(len= 64) :: restart_freq_opt, restart_start_opt, restart_fmt
logical :: leven_odd_on, pressure_correction
integer :: restart_freq, restart_start, even_odd_freq

namelist /restart_nml/ restart_freq_opt, restart_freq, &
     restart_start_opt, restart_start, restart_outfile, &
    restart_fmt, leven_odd_on, even_odd_freq, pressure_correction

!------------------------------------------------------------------
! The POP initial temperature and salinity namelist
!------------------------------------------------------------------

character(len=100) :: init_ts_file ! length consistent with POP
character(len=100) :: init_ts_outfile
character(len= 64) :: init_ts_option, init_ts_suboption
character(len= 64) :: init_ts_file_fmt, init_ts_outfile_fmt
real(r8) :: init_ts_perturb
namelist /init_ts_nml/ init_ts_option, init_ts_suboption, &
                       init_ts_file, init_ts_file_fmt, &
                       init_ts_outfile, init_ts_outfile_fmt, &
                       init_ts_perturb

!------------------------------------------------------------------
! The POP domain namelist
!------------------------------------------------------------------

character(len= 64) :: clinic_distribution_type, tropic_distribution_type
character(len= 64) :: ew_boundary_type, ns_boundary_type
integer :: nprocs_clinic, nprocs_tropic
logical :: profile_barrier

namelist /domain_nml/ clinic_distribution_type, nprocs_clinic, &
                      tropic_distribution_type, nprocs_tropic, &
                      ew_boundary_type, ns_boundary_type, &
                      profile_barrier

!------------------------------------------------------------------
! The POP grid info namelist
!------------------------------------------------------------------
!
! POP grid information comes in several files:
!   horizontal grid lat/lons in one,
!   topography (lowest valid vert level) in another, and
!   the vertical grid spacing in a third.
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
! The vert grid file is ascii, with 3 columns/line:
!    cell thickness(in cm)   cell center(in m)   cell bottom(in m)
!
!------------------------------------------------------------------

character(len=100) :: horiz_grid_file, vert_grid_file, &
                      topography_file, topography_outfile, &
                      bathymetry_file, region_info_file, &
                      bottom_cell_file, region_mask_file
character(len= 64) :: horiz_grid_opt, sfc_layer_opt, vert_grid_opt, &
                      topography_opt
logical :: partial_bottom_cells, topo_smooth, flat_bottom, lremove_points, &
           l1ddyn
integer :: kmt_kmin, n_topo_smooth

namelist /grid_nml/ horiz_grid_opt, horiz_grid_file, sfc_layer_opt, &
    vert_grid_opt, vert_grid_file, topography_opt, kmt_kmin, &
    topography_file, topography_outfile, bathymetry_file, &
    partial_bottom_cells, bottom_cell_file, n_topo_smooth, &
    region_mask_file, topo_smooth, flat_bottom, lremove_points, &
    region_info_file, l1ddyn

!======================================================================
contains
!======================================================================


subroutine initialize_module

integer :: iunit, io

! Read POP calendar information
! In 'restart' mode, this is primarily the calendar type and 'stop'
! information. The time attributes of the restart file override
! the namelist time information.

call find_namelist_in_file('pop_in', 'time_manager_nml', iunit)
read(iunit, nml = time_manager_nml, iostat = io)
call check_namelist_read(iunit, io, 'time_manager_nml')

! FIXME : Real observations are always GREGORIAN dates ...
! but stomping on that here gets in the way of running
! a perfect_model experiment for pre-1601 AD cases.

! STOMP if ( allow_leapyear ) then
   call set_calendar_type('gregorian')
! STOMP else
! STOMP    call set_calendar_type('noleap')
! STOMP endif

! Read POP I/O information (for restart file ... grid dimensions)
! Read POP initial information (for input/restart filename)

call find_namelist_in_file('pop_in', 'io_nml', iunit)
read(iunit, nml = io_nml, iostat = io)
call check_namelist_read(iunit, io, 'io_nml')

call find_namelist_in_file('pop_in', 'init_ts_nml', iunit)
read(iunit, nml = init_ts_nml, iostat = io)
call check_namelist_read(iunit, io, 'init_ts_nml')

! Make sure we have a pop restart file (for grid dims)
if ( .not. file_exist(ic_filename) ) then
   string1 = 'pop_in:init_ts_file '//trim(ic_filename)//' not found'
   call error_handler(E_ERR,'initialize_module', &
          string1, source, revision, revdate)
endif

! Read POP restart information (for model timestepping)
call find_namelist_in_file('pop_in', 'restart_nml', iunit)
read(iunit, nml = restart_nml, iostat = io)
call check_namelist_read(iunit, io, 'restart_nml')

! Read POP domain information (for lon wrapping or not)
call find_namelist_in_file('pop_in', 'domain_nml', iunit)
read(iunit, nml = domain_nml, iostat = io)
call check_namelist_read(iunit, io, 'domain_nml')

! Read POP grid information (for grid filenames)
call find_namelist_in_file('pop_in', 'grid_nml', iunit)
read(iunit, nml = grid_nml, iostat = io)
call check_namelist_read(iunit, io, 'grid_nml')

module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

end subroutine initialize_module


!------------------------------------------------------------------
!> Read the lon, lat grid size from the restart netcdf file.
!> The actual grid file is a binary file with no header information.
!>
!> The file name comes from module storage ... namelist.


subroutine get_horiz_grid_dims(Nx, Ny)

integer, intent(out) :: Nx   ! Number of Longitudes
integer, intent(out) :: Ny   ! Number of Latitudes

integer :: grid_id, dimid, nc_rc

if ( .not. module_initialized ) call initialize_module

! get the ball rolling ...

call nc_check(nf90_open(trim(ic_filename), nf90_nowrite, grid_id), &
         'get_horiz_grid_dims','open '//trim(ic_filename))

! Longitudes : get dimid for 'i' or 'nlon', and then get value
nc_rc = nf90_inq_dimid(grid_id, 'i', dimid)
if (nc_rc /= nf90_noerr) then
   nc_rc = nf90_inq_dimid(grid_id, 'nlon', dimid)
   if (nc_rc /= nf90_noerr) then
      string1 = "unable to find either 'i' or 'nlon' in file "//trim(ic_filename)
      call error_handler(E_ERR, 'get_horiz_grid_dims', string1, &
                         source,revision,revdate)
   endif
endif

call nc_check(nf90_inquire_dimension(grid_id, dimid, len=Nx), &
         'get_horiz_grid_dims','inquire_dimension i '//trim(ic_filename))

! Latitudes : get dimid for 'j' or 'nlat' ... and then get value
nc_rc = nf90_inq_dimid(grid_id, 'j', dimid)
if (nc_rc /= nf90_noerr) then
   nc_rc = nf90_inq_dimid(grid_id, 'nlat', dimid)
   if (nc_rc /= nf90_noerr) then
      string1 = "unable to find either 'j' or 'nlat' in "//trim(ic_filename)
      call error_handler(E_ERR, 'get_horiz_grid_dims', string1, &
                         source,revision,revdate)
   endif
endif

call nc_check(nf90_inquire_dimension(grid_id, dimid, len=Ny), &
         'get_horiz_grid_dims','inquire_dimension i '//trim(ic_filename))

! tidy up

call nc_check(nf90_close(grid_id), &
         'get_horiz_grid_dims','close '//trim(ic_filename) )

end subroutine get_horiz_grid_dims


!------------------------------------------------------------------
!> count the number of lines in the ascii file to figure out max
!> number of vert blocks.


subroutine get_vert_grid_dim(Nz)

integer, intent(out) :: Nz

integer :: linelen ! disposable

if ( .not. module_initialized ) call initialize_module

call find_textfile_dims(vert_grid_file, Nz, linelen)

end subroutine get_vert_grid_dim


!------------------------------------------------------------------
!> the initialize_module ensures that the pop namelists are read and
!> the DART time manager gets the pop calendar setting.
!>
!> Then, the DART time manager is queried to return what it knows ...


subroutine get_pop_calendar(calstring)


character(len=*), INTENT(OUT) :: calstring

if ( .not. module_initialized ) call initialize_module

call get_calendar_string(calstring)

end subroutine get_pop_calendar


!------------------------------------------------------------------
!> the initialize_module ensures that the pop namelists are read.
!> The restart times in the pop_in&restart_nml are used to define
!> appropriate assimilation timesteps.


function set_model_time_step(seconds,days)

integer, intent(in) :: seconds ! input.nml assimilation_period_seconds
integer, intent(in) :: days    ! input.nml assimilation_period_days
type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call initialize_module

! Check to see if the input.nml:&model_nml:assimilation_period_[seconds,days]
! are specifying the nominal model output frequency/assimilation interval.

if (seconds <= 0 .and. days <= 0) then

   ! Check the 'restart_freq_opt' and 'restart_freq' to determine
   ! when we can stop the model

   if ( trim(restart_freq_opt) == 'nday' ) then
      set_model_time_step = set_time(0, restart_freq) ! (seconds, days)
   else if ( trim(restart_freq_opt) == 'nyear' ) then
      ! FIXME ... CCSM_POP uses a bogus value for this
      set_model_time_step = set_time(0, 1) ! (seconds, days)
      write(string1,*)'WARNING: POP namelist variable "restart_freq_opt" read as "nyear"'
      write(string2,*)'CESM uses bogus value - using default value of 1 day.'
      call error_handler(E_MSG,'set_model_time_step', string1, &
                 source, revision, revdate, text2=string2)
   else
      write(string1,*)'POP namelist variable "restart_freq_opt" must be "nday" or "nyear"'
      write(string2,*)'read as "'//trim(restart_freq_opt)//'"'
      call error_handler(E_ERR,'set_model_time_step', string1, &
                 source, revision, revdate, text2=string2)
   endif

elseif (seconds < 0 .or. days < 0) then
   write(string1,*)'Unable to determine the assimilation interval.'
   write(string2,*)'input.nml:&model_nml:assimilation_period_seconds must be >= 0; is ',seconds
   write(string3,*)'input.nml:&model_nml:assimilation_period_days    must be >= 0; is ',days
   call error_handler(E_ERR,'set_model_time_step', string1, &
              source, revision, revdate, text2=string2, text3=string3)
else
   set_model_time_step = set_time(seconds, days) ! (seconds, days)
endif

end function set_model_time_step


!------------------------------------------------------------------


subroutine write_pop_namelist(model_time, adv_to_time)

type(time_type), INTENT(IN) :: model_time, adv_to_time
type(time_type) :: offset

integer :: iunit, secs, days

if ( .not. module_initialized ) call initialize_module

offset = adv_to_time - model_time
call get_time(offset, secs, days)

if (secs /= 0 ) then
   write(string1,*)'adv_to_time has seconds == ',secs,' must be zero'
   call error_handler(E_ERR,'write_pop_namelist', string1, source, revision, revdate)
endif

! call print_date( model_time,'write_pop_namelist:dart model date')
! call print_date(adv_to_time,'write_pop_namelist:advance_to date')
! call print_time( model_time,'write_pop_namelist:dart model time')
! call print_time(adv_to_time,'write_pop_namelist:advance_to time')
! call print_time(     offset,'write_pop_namelist:a distance of')
! write( *,'(''write_pop_namelist:TIME_MANAGER_NML   STOP_COUNT '',i10,'' days'')') days

!Convey the information to the namelist 'stop option' and 'stop count'

if ( trim(stop_option) == 'nday' ) then
   stop_count = days
elseif ( trim(stop_option) == 'nyear' ) then
   if (days > 365) then
      stop_count = days/365   ! relying on integer arithmetic
   else
      ! CESM totally ignores this value
      write(string1,*)'POP time_manager_nml:stop_option,stop_count are ',trim(stop_option),stop_count
      write(string2,*)'DART wants to advance ',days,' "days"'
      write(string3,*)'Unable to reconcile; using original stop_option,stop_count.'
      call error_handler(E_WARN,'write_pop_namelist', string1, &
                 source, revision, revdate, text2=string2)
      continue
   endif
else
   call error_handler(E_ERR,'write_pop_namelist', &
              'stop_option must be "nday" or "nyear"', source, revision, revdate)
endif

iunit = open_file('pop_in.DART',form='formatted',action='rewind')
write(iunit, nml=time_manager_nml)
write(iunit, '('' '')')
close(iunit)

end subroutine write_pop_namelist


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

if ( .not. file_exist(horiz_grid_file) ) then
   string1 = 'pop_in:horiz_grid_file '//trim(horiz_grid_file)//' not found'
   call error_handler(E_ERR,'read_horiz_grid', &
          string1, source, revision, revdate)
endif

! Open it and read them in the EXPECTED order.
! Actually, we only need the first two, so I'm skipping the rest.

grid_unit = get_unit()
INQUIRE(iolength=reclength) ULAT64

open(grid_unit, file=trim(horiz_grid_file), form='unformatted', convert=conversion, &
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

! Define some constants as in pop

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
   write(string1,'(''pop_in&domain_nml:ew_boundary_type '',a,'' unknown.'')') &
                                    trim(ew_boundary_type)
   call error_handler(E_ERR,'calc_tpoints',string1,source,revision,revdate)
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

if ( .not. file_exist(topography_file) ) then
   string1 = 'pop_in:topography_file '//trim(topography_file)//' not found'
   call error_handler(E_ERR,'read_topography', &
          string1, source, revision, revdate)
endif

! read the binary file

topo_unit = get_unit()
INQUIRE(iolength=reclength) KMT

open( topo_unit, file=trim(topography_file), form='unformatted', convert=conversion, &
            access='direct', recl=reclength, status='old', action='read' )
read( topo_unit, rec=1) KMT
close(topo_unit)

! the equation numbered 3.2 on page 15 of this document:
!  http://www.cesm.ucar.edu/models/cesm1.0/pop2/doc/sci/POPRefManual.pdf
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
!> Open and read the ASCII vertical grid information


subroutine read_vert_grid(nz, ZC, ZG)

! The vert grid file is in ascii, with either 3 columns/line
!    cell thickness(in cm)   cell center(in m)   cell bottom(in m)
! or it can contain 2 columns/line
!    cell thickness(in cm)   1.0
! in which case we compute the cell center (ZC) and cell bottom (ZG)
! and ignore the second column.

integer,  intent(in)  :: nz
real(r8), intent(out) :: ZC(nz), ZG(nz)

integer  :: iunit, i, ios
real(r8) :: depth

logical :: three_columns
character(len=256) :: line

real(r8), parameter :: centimeters_to_meters = 0.01_r8

if ( .not. module_initialized ) call initialize_module

! Check to see that the file exists.

if ( .not. file_exist(vert_grid_file) ) then
   string1 = 'pop_in:vert_grid_file '//trim(vert_grid_file)//' not found'
   call error_handler(E_ERR,'read_vert_grid', &
          string1, source, revision, revdate)
endif

! read the ASCII file

iunit = open_file(trim(vert_grid_file), action = 'read')

! determine the number of columns
read(iunit,'(A)') line

read(line,*,iostat=ios) depth, ZC(1), ZG(1)

if(ios == 0) then
   three_columns = .true.
else
   three_columns = .false.

   ! read depth and calculate center and bottom of cells
   read(line,*,iostat=ios) depth

   ZC(1) = depth*centimeters_to_meters*0.5_r8
   ZG(1) = depth*centimeters_to_meters
endif

do i=2, nz

   if(three_columns) then
      read(iunit,*,iostat=ios) depth, ZC(i), ZG(i)
   else
      read(iunit,*,iostat=ios) depth

      ZC(i) = ZG(i-1) + depth*centimeters_to_meters*0.5_r8
      ZG(i) = ZG(i-1) + depth*centimeters_to_meters
   endif

   if ( ios /= 0 ) then ! error
      write(string1,*)'error reading depths, line ',i
      call error_handler(E_ERR,'read_vert_grid',string1,source,revision,revdate)
   endif

enddo

end subroutine read_vert_grid


!------------------------------------------------------------------
!> Open and read the mean dynamic sea surface topography
!> There is an assumed name and the shape of the variable must
!> match the POP grid being used. The hope is that that actual
!> locations MUST MATCH EXACTLY the POP SSH grid.

subroutine read_mean_dynamic_topography(fname, mdt)
character(len=*), intent(in)  :: fname
real(r8),         intent(out) :: mdt(:,:)

integer  :: ncid, VarID, io, xtype, ii
integer  :: numdims, dimIDs(NF90_MAX_DIMS), dimlen
character(len=128) :: unitsstring
real(r4) :: rmiss
real(r8) :: dmiss, dmin, dmax

if ( .not. module_initialized ) call initialize_module

mdt = MISSING_R8
dmiss = MISSING_R8

! Check to see that the netCDF file exists.

if ( .not. file_exist(fname) ) then
   string1 = 'Mean dynamic sea surface topography file not found.'
   string2 = 'Looking for filename "'//trim(fname)//'"'
   string3 = 'This filename is specified by input.nml:model_nml:mdt_reference_file_name'
   call error_handler(E_ERR,'read_mean_dynamic_topography', &
          string1, source, revision, revdate, text2=string2, text3=string3)
endif

io = nf90_open(trim(fname), NF90_NOWRITE, ncid)
call nc_check(io, 'read_mean_dynamic_topography','open "'//trim(fname)//'"')

io = nf90_inq_varid(ncid, 'mdt', VarID)
call nc_check(io, 'read_mean_dynamic_topography', &
                  '"mdt" not found but is required in "'//trim(fname)//'"')

! check variable shape against assumed shape

io = nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=numdims, xtype=xtype)
call nc_check(io, 'read_mean_dynamic_topography', &
                  'inquire_variable "mdt" from "'//trim(fname)//'"')

do ii = 1,numdims
   write(string1,*)'inquire_dimension ',ii,'for "mdt"'
   io = nf90_inquire_dimension(ncid,dimIDs(ii),len=dimlen)
   call nc_check(io, 'read_mean_dynamic_topography', string1)

   if (dimlen /= size(mdt,ii)) then
      write(string1,*)'mdt dimension mismatch'
      write(string2,*)'mdt dimension ',ii,' in file is ',dimlen
      write(string3,*)'mdt dimension ',ii,' in code is ',size(mdt,ii)
      call error_handler(E_ERR,'read_mean_dynamic_topography', &
          string1, source, revision, revdate, text2=string2, text3=string3)
   endif

enddo

io = nf90_get_var(ncid, VarID, mdt)
call nc_check(io, 'read_mean_dynamic_topography', &
                  'get_var "mdt" from "'//trim(fname)//'"')

! Replace _FillValue with something
!>@todo CHECK ... does it xtype of the variable matter when getting a _FillValue attribute. 
!> If it does not, these lines could collapse into something cleaner.

if (xtype == NF90_REAL) then
   io = nf90_get_att(ncid, VarId, '_FillValue', rmiss)
   if (io == NF90_NOERR) then
      dmiss=rmiss
      where (mdt == dmiss) mdt = MISSING_R8
   endif
elseif (xtype == NF90_DOUBLE) then
   io = nf90_get_att(ncid, VarId, '_FillValue', dmiss)
   if (io == NF90_NOERR) then
      where (mdt == dmiss) mdt = MISSING_R8
   endif
else
   call error_handler(E_ERR,'read_mean_dynamic_topography', &
        'unsupported variable type for "mdt"', source, revision, revdate, &
        text2 = 'must be "float" or "double"')
endif

! The observations are in meters, the mean dynamic topography should be in meters
! and internally in the POP model_mod, the CGS units are converted to SI.

io = nf90_get_att(ncid, VarId, 'units', unitsstring)
call nc_check(io, 'read_mean_dynamic_topography', &
                  'get_att "units" for "mdt" from "'//trim(fname)//'"')

if (unitsstring == 'centimeter' .or. unitsstring == 'cm') then
   where(mdt /= MISSING_R8) mdt = mdt/100.0_r8
elseif (unitsstring == 'm' .or. unitsstring == 'meters') then
   continue
else
   call error_handler(E_ERR,'read_mean_dynamic_topography', &
        'unsupported units for "mdt"', source, revision, revdate, &
        text2 = 'must be "centimeter", "cm", "meter", or "m"')
endif

call nc_check(nf90_close(ncid), &
              'read_mean_dynamic_topography','close ' // trim(fname) )

! Just something to do a basic check.

dmin = minval(mdt, mdt /= MISSING_R8) 
dmax = maxval(mdt, mdt /= MISSING_R8) 
write(string1,*)'..  mdt sizes are ',size(mdt,1),size(mdt,2)
write(string2,*)'_FillValue was ',dmiss,' original units "'//trim(unitsstring)//'"'
write(string3,*)'min,max (meters) ',dmin, dmax
call error_handler(E_MSG,'read_mean_dynamic_topography', &
     string1, text2=string2, text3=string3)

end subroutine read_mean_dynamic_topography


!------------------------------------------------------------------


subroutine get_pop_restart_filename( filename )

character(len=*), intent(out) :: filename

if ( .not. module_initialized ) call initialize_module

filename   = trim(ic_filename)

end subroutine get_pop_restart_filename


!------------------------------------------------------------------


subroutine set_binary_file_conversion(convertstring)

character(len=*), intent(in) :: convertstring

if ( .not. module_initialized ) call initialize_module

conversion = convertstring

end subroutine set_binary_file_conversion


!------------------------------------------------------------------


end module dart_pop_mod


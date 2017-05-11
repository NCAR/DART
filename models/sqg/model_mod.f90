! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

!========================================================================
! Assimilation interface for:
! Uniform-PV two-surface QG+1 model in spectral form (Hakim 2000)
!========================================================================

!================= m o d u l e   i n f o r m a t i o n ==================

use        types_mod, only : r8, MISSING_R8
use time_manager_mod, only : time_type, set_time, get_time, set_calendar_type, &
                             print_time, operator(<), operator(+)
use     location_mod, only : location_type,      get_close_maxdist_init, &
                             get_close_obs_init, get_close_obs, set_location, &
                             get_location, set_location_missing
use    utilities_mod, only : register_module, error_handler, nc_check, &
                             E_ERR, E_MSG, logfileunit, get_unit, close_file, &
                             dump_unit_attributes, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read
use     obs_kind_mod, only : QTY_POTENTIAL_TEMPERATURE
use     spectral_mod
use          sqg_mod, only : diffusion, init, init_jet, terrain, invert, advect, &
                             tadv, xy_to_sp, sp_to_xy, d_setup, ft_2d, norm

!========================================================================

implicit none

private

public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          end_model,              &
          static_init_model,      &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

! public subroutines / functions
public :: sqg_to_dart,            &
          dart_to_sqg,            &
          get_model_static_data

! public types
public :: model_static

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!---------------------------------------------------------------
! Namelist with default values

logical  :: output_state_vector         = .false.
real(r8) :: channel_center              = 45.0_r8
real(r8) :: channel_width               = 40.0_r8
logical  :: debug                       = .false.
integer  :: assimilation_period_days    = 0
integer  :: assimilation_period_seconds = 60*60*6

namelist /model_nml/ channel_center, channel_width, &
                     assimilation_period_days, assimilation_period_seconds, &
                     debug

real,     dimension(mmax,nmax), parameter :: Rblank = 0.0
complex,  dimension(mmax,nmax), parameter :: Cblank = 0.0

integer, parameter :: model_size = 2 * (2*kmax) * (2*lmax)

! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
    
type(time_type)     :: time_step, model_stride
character(len=129)  :: msgstring1,msgstring2

type model_static
   logical                               :: top, bot
   real                                  :: dco, lam
   complex,  allocatable, dimension(:,:) :: thbB, thbT
   real,     allocatable, dimension(:,:) :: thbyB, thbyT, ulinB, ulinT
   real,     allocatable, dimension(:,:) :: hx, hy, hu, hv
   real(r8), allocatable, dimension(:)   :: lons, lats, levs
end type model_static

type(model_static) :: sqg_static

!========================================================================

contains

!========================================================================

subroutine static_init_model()
!-----------------------------
! Initializes class data for surface quasigeostrophy model
! (all the stuff that needs to be done once.)

integer  :: iunit, io

integer  :: i, lev_index, lat_index, lon_index
real(r8) :: lon, lat, lev

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! This is where you would read a namelist, for example.
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type('Gregorian')  ! TJH addition

model_stride = set_model_stride()

! evenly spaced longitudes, latitudes and two levels
allocate( sqg_static%lons(2*kmax) )
allocate( sqg_static%lats(2*lmax) )
allocate( sqg_static%levs(2     ) )
do i = 1, 2*kmax
   sqg_static%lons(i) = 360.0*(i-1.0)/(2*kmax)
enddo
! channel is centered about "channel_center" and has a width of "channel_width"
do i = 1, 2*lmax
   sqg_static%lats(i) = channel_center + channel_width * (real(i,r8)/(2*lmax) - 0.5)
enddo
do i = 1, 2
   sqg_static%levs(i) = i
enddo

! allocate storage for locations
allocate( state_loc(model_size) )
do i = 1, model_size
   lev_index =  (i-1) / ((2*lmax)*(2*kmax)) + 1
   lat_index = ((i-1) - ((lev_index-1)*(2*lmax)*(2*kmax))) / (2*kmax) + 1
   lon_index =  (i-1) - ((lev_index-1)*(2*lmax)*(2*kmax)) - ((lat_index-1)*(2*kmax)) + 1

   lev = sqg_static%levs(lev_index)
   lat = sqg_static%lats(lat_index)
   lon = sqg_static%lons(lon_index)

   ! With the location module ... you specify that the 
   ! vertical coordinate is a 'level' by 'which_vert' == 1
   state_loc(i) = set_location(lon, lat, lev, 1)
enddo

! get dimensional model time-step,
! follow non-dimensionalization conversion for sQG (Mahajan & Hakim 2011, Table 1)
! dTs =  (Ls/Us) * dt ; Ls = 1000 km = 1e6 m ; Us = 30 m/s
time_step = set_time(int((1e6/30)*dt), 0)

! initialize derivative operators in the model
call d_setup()

! initialize diffusion coefficient
call diffusion(sqg_static%dco)

! advection flags
if     ( (model .eq. 0) .or. (model .eq. 1) ) then  ! 2D or 2sQG
   sqg_static%top = .true.; sqg_static%bot = .true.
elseif (  model .eq. 2 )                      then  ! tropo sQG
   sqg_static%top = .true.; sqg_static%bot = .false.
elseif ( (model .eq. 3) .or. (model .eq. 0) ) then  ! surface sQG
   sqg_static%top = .false.; sqg_static%bot = .true.
elseif (  model .eq. 4 )                      then  ! tropo HsQG
   sqg_static%top = .true.; sqg_static%bot = .false.
endif
!if (maxval(abs(thxyT)) .lt. 1.e-5) sqg_static%top = .false.
!if (maxval(abs(thxyB)) .lt. 1.e-5) sqg_static%bot = .false.

! initialize Hoskins-West jet
allocate(sqg_static%thbB(2*kmax,2*lmax)) ; allocate(sqg_static%thbT(2*kmax,2*lmax))
allocate(sqg_static%thbyB(mmax,nmax))    ; allocate(sqg_static%thbyT(mmax,nmax))
allocate(sqg_static%ulinB(mmax,nmax))    ; allocate(sqg_static%ulinT(mmax,nmax))
if ( hw ) then

   call init_jet(sqg_static%thbB,  sqg_static%thbT,  &
                 sqg_static%thbyB, sqg_static%thbyT, &
                 sqg_static%ulinB, sqg_static%ulinT, &
                 sqg_static%lam)

   ! barotropic wind (Ross = 0!)
   sqg_static%ulinB = sqg_static%ulinB + 0.0
   sqg_static%ulinT = sqg_static%ulinT - 0.0 * H * sqg_static%lam

else

   sqg_static%thbB  = 0.0 ; sqg_static%thbT  = 0.0
   sqg_static%thbyB = 0.0 ; sqg_static%thbyT = 0.0
   sqg_static%ulinB = 0.0 ; sqg_static%ulinT = 0.0
   sqg_static%lam   = 0.0

endif

! initialize terrain
allocate(sqg_static%hx(mmax,nmax)) ; allocate(sqg_static%hy(mmax,nmax))
allocate(sqg_static%hu(mmax,nmax)) ; allocate(sqg_static%hv(mmax,nmax))
if ( (iterr) .and. (Ross .eq. 0.0) ) then 
   call terrain(sqg_static%hx, sqg_static%hy, &
                sqg_static%hu, sqg_static%hv)
else
   sqg_static%hx = 0.0 ; sqg_static%hy = 0.0
   sqg_static%hu = 0.0 ; sqg_static%hv = 0.0
endif

end subroutine static_init_model

!========================================================================

subroutine init_conditions(x)
!----------------------------
! read initial conditions from a file,
! must only contain perturbation theta

real(r8), intent(out) :: x(:)

integer                                :: j
complex, allocatable, dimension(:,:)   :: thspB,  thspT
real,    allocatable, dimension(:,:)   :: thxyB,  thxyT
real,    allocatable, dimension(:,:,:) :: theta

! initialize the theta fields:
allocate(thxyB(2*kmax,2*lmax)) ; allocate(thxyT(2*kmax,2*lmax))
call init('sqgInput.nc',thxyB,thxyT)

! first add base-state jet:
if ( hw ) then

   allocate(thspB(2*kmax,2*lmax)) ; allocate(thspT(2*kmax,2*lmax))

   ! map into spectral space at the same resolution:
   call xy_to_sp(cmplx(thxyB,0.0),thspB,2*kmax,2*lmax,kmax,lmax)
   call xy_to_sp(cmplx(thxyT,0.0),thspT,2*kmax,2*lmax,kmax,lmax)

   thspB = thspB + sqg_static%thbB
   thspT = thspT + sqg_static%thbT

   ! map into grid-point space space at the same resolution:
   call sp_to_xy(thspB,thxyB,kmax,lmax,2*kmax,2*lmax)
   call sp_to_xy(thspT,thxyT,kmax,lmax,2*kmax,2*lmax)

   deallocate(thspB) ; deallocate(thspT)

endif

! second add linear shear:
allocate(theta(2*kmax,2*lmax,2))
do j = 1, 2*lmax
   theta(:,j,1) = thxyB(:,j) - sqg_static%lam * real(j-1) * YL/real(2*lmax) 
   theta(:,j,2) = thxyT(:,j) - sqg_static%lam * real(j-1) * YL/real(2*lmax) 
enddo

! wrap model variables into DART state vector:
call sqg_to_dart(theta,x)

deallocate(thxyB) ; deallocate(thxyT)
deallocate(theta)

end subroutine init_conditions

!========================================================================

subroutine adv_1step(x, current_time)
!-------------------------------------------------
! Does a time advance for sQG model with state vector as
! input and output.
! This interface advances from current_time to target_time
! and should be used with async = -2
! see modifications in obs_model_mod.f90

real(r8),        intent(inout)        :: x(:)
type(time_type), intent(in)           :: current_time

real(r8) :: cxB,cyB,cxT,cyT
integer  :: j
integer  :: itime, cdays, cseconds, tdays, tseconds
type(time_type) :: ctime, target_time

real,    allocatable, dimension(:,:,:) :: thxy
real,    allocatable, dimension(:,:)   :: thxyB, thxyT
real,    allocatable, dimension(:,:)   :: laplacian
real,    allocatable, dimension(:,:)   :: thxB, thxT, thyB, thyT
real,    allocatable, dimension(:,:)   :: uB, uT, vB, vT
real,    allocatable, dimension(:,:)   :: tthB, tthT

complex, allocatable, dimension(:,:)   :: thspB, thspT, thspB1, thspT1, thspB2, thspT2
complex, allocatable, dimension(:,:)   :: sB, sBold
complex, allocatable, dimension(:,:)   :: tthspB, tthspT

logical :: first

! allocate necessary space up-front:
allocate(thxy(2*kmax,2*lmax,2))
allocate(thxyB( 2*kmax,2*lmax)) ; allocate(thxyT( 2*kmax,2*lmax))

allocate(laplacian(mmax,nmax))
allocate(thxB(mmax,nmax)) ; allocate(thxT(mmax,nmax))
allocate(thyB(mmax,nmax)) ; allocate(thyT(mmax,nmax))
allocate(uB(  mmax,nmax)) ; allocate(uT(  mmax,nmax))
allocate(vB(  mmax,nmax)) ; allocate(vT(  mmax,nmax))
allocate(tthB(mmax,nmax)) ; allocate(tthT(mmax,nmax))

allocate(thspB( 2*kmax,2*lmax)) ; allocate(thspT( 2*kmax,2*lmax))
allocate(thspB1(2*kmax,2*lmax)) ; allocate(thspT1(2*kmax,2*lmax))
allocate(thspB2(2*kmax,2*lmax)) ; allocate(thspT2(2*kmax,2*lmax))

allocate(sB(2*kmax,2*lmax)) ; allocate(sBold(2*kmax,2*lmax))

allocate(tthspB(mmax,nmax)) ; allocate(tthspT(mmax,nmax))

! unwrap DART state vector into model variables:
call dart_to_sqg(x, thxy)

! first remove the linear shear:
do j = 1, 2*lmax
   thxyB(:,j) = thxy(:,j,1) + sqg_static%lam * real(j-1) * YL/real(2*lmax) 
   thxyT(:,j) = thxy(:,j,2) + sqg_static%lam * real(j-1) * YL/real(2*lmax) 
enddo

! map into spectral space at the same resolution:
call xy_to_sp(cmplx(thxyB,0.0), thspB, 2*kmax,2*lmax, kmax, lmax)
call xy_to_sp(cmplx(thxyT,0.0), thspT, 2*kmax,2*lmax, kmax, lmax)

! second remove the HW base-state:
if ( hw ) then
   thspB = thspB - sqg_static%thbB
   thspT = thspT - sqg_static%thbT
endif

! initialize certain variables for first time-step
first  = .true.
sB     = 0.0
thspB1 = 0.0; thspB2 = 0.0
thspT1 = 0.0; thspT2 = 0.0

! advance from current_time 
itime = 1
ctime = current_time
target_time = current_time + model_stride
do while ( ctime < target_time )

   ! save old stream-function for Ekman calculation:
   sBold = sB

   ! invert theta for stream-function; compute gradients for advection:
   call invert(thspB, thspT, thxB, thxT, thyB, thyT, vB, vT, uB, uT, &
               sqg_static%thbB, sqg_static%thbT, sqg_static%thbyB, sqg_static%thbyT, &
               sqg_static%ulinB, sqg_static%ulinT, &
               first, sqg_static%bot, sqg_static%top, sqg_static%lam, &
               sB, sBold, laplacian)

   ! option to compute potential enstrophy norm and growth rate:
   if ( inorm ) call norm(thspB, thspT, itime)

   ! spectral advection:
   if ( sqg_static%bot ) &
      call advect(uB, vB, thxB, thyB, &
                  sqg_static%thbyB, sqg_static%hx, sqg_static%hy, sqg_static%ulinB, &
                  tthB, sqg_static%lam, laplacian)
   if ( sqg_static%top ) &
      call advect(uT + sqg_static%hu, vT + sqg_static%hv, thxT, thyT, &
                  sqg_static%thbyT, Rblank, Rblank, sqg_static%ulinT, &
                  tthT, sqg_static%lam, Rblank)

   ! compute Courant numbers and print to stdout if debugging:
   if ( mod((itime-1),10) .eq. 0 ) then
      cxB = maxval(abs(uB + sqg_static%ulinB)) * dt / (XL/real(2*kmax))
      cyB = maxval(abs(vB)) * dt / (YL/real(2*lmax))
      cxT = maxval(abs(uT + sqg_static%ulinT)) * dt / (XL/real(2*kmax))
      cyT = maxval(abs(vT)) * dt / (YL/real(2*lmax))
      if ( debug ) write(*,'(A23,F10.3,4F8.3))') 'time,cxB,cyB,cxT,cyT = ', &
                                      real(itime-1)*dt,cxB,cyB,cxT,cyT
   endif

   ! FFT back to spectral space:
   if ( sqg_static%bot ) then
      tthspB = cmplx(tthB,0.0)
      call ft_2d(tthspB,mmax,nmax,-1)
   endif
   if ( sqg_static%top ) then
      tthspT = cmplx(tthT,0.0)
      call ft_2d(tthspT,mmax,nmax,-1)
   endif

   ! advance one time-step with explicit (hyper-) diffusion:
   if ( sqg_static%bot ) call tadv(thspB, tthspB, thspB1, thspB2, sqg_static%dco, first)
   if ( sqg_static%top ) call tadv(thspT, tthspT, thspT1, thspT2, sqg_static%dco, first)

   ! zero out k=K, l=L modes on the bottom boundary:
   thspB(kmax,:) = 0.0
   thspB(lmax,:) = 0.0
       
   itime = itime + 1
   ctime = ctime + time_step

   first = .false.

enddo

! first add the HW base-state:
if ( hw ) then
   thspB = thspB + sqg_static%thbB
   thspT = thspT + sqg_static%thbT
endif

! map into grid-point space at the same resolution:
call sp_to_xy(thspB, thxyB, kmax, lmax, 2*kmax, 2*lmax)
call sp_to_xy(thspT, thxyT, kmax, lmax, 2*kmax, 2*lmax)

! second add the linear shear:
do j = 1, 2*lmax
   thxy(:,j,1) = thxyB(:,j) - sqg_static%lam * real(j-1) * YL/real(2*lmax) 
   thxy(:,j,2) = thxyT(:,j) - sqg_static%lam * real(j-1) * YL/real(2*lmax) 
enddo

! wrap model variables into DART state vector:
call sqg_to_dart(thxy, x)

! deallocate space allocated up-front:
deallocate(thxy)
deallocate(thxyB ) ; deallocate(thxyT )

deallocate(laplacian)
deallocate(thxB) ; deallocate(thxT)
deallocate(thyB) ; deallocate(thyT)
deallocate(uB)   ; deallocate(uT)
deallocate(vB)   ; deallocate(vT)
deallocate(tthB) ; deallocate(tthT)

deallocate(thspB ) ; deallocate(thspT )
deallocate(thspB1) ; deallocate(thspT1)
deallocate(thspB2) ; deallocate(thspT2)

deallocate(sB) ; deallocate(sBold)

deallocate(tthspB) ; deallocate(tthspT)

end subroutine adv_1step

!========================================================================

function get_model_size()
!------------------------
! Returns the size of the model as an integer.

integer :: get_model_size

get_model_size = model_size

end function get_model_size

!========================================================================

subroutine init_time(time)
!-------------------------
! For now, returns the initialization time as 0

type(time_type), intent(out) :: time

time = set_time(0,0)

end subroutine init_time

!========================================================================

subroutine model_interpolate(x, location, itype, obs_val, istatus)
!-----------------------------------------------------------------
! Interpolates from state vector x to the location and returns obs_val.
! istatus = 0 suggests interpolation went ok, 1 means something went wrong.
! Argument itype specifies the type of field (for instance potential temperature in this model)

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

real(r8) :: lon, lat, lev, lon_lat_lev(3)
real(r8) :: bot_lon, top_lon, delta_lon, bot_lat, top_lat, delta_lat
real(r8) :: temp_lon, lon_fract, lat_fract, val(2,2), a(2)
integer  :: lon_below, lon_above, lat_below, lat_above, level, i

! Assume all interpolations okay for now
istatus = 0

lon_lat_lev = get_location(location)
lon = lon_lat_lev(1)
lat = lon_lat_lev(2)
lev = lon_lat_lev(3); level = int(lev)

! Level is obvious for now, but make sure that it is within valid range
if ( (level < 1) .or. (level > 2) ) then
   istatus = 1
   obs_val = MISSING_R8
   write(msgstring1,*)'level ',level,' must be 1 or 2'
   call error_handler(E_ERR,'model_mod:model_interpolate', msgstring1, source, revision, revdate)
endif

! Get globally defined lat / lon grid specs
bot_lon   = sqg_static%lons(1)
top_lon   = sqg_static%lons(2*kmax)
bot_lat   = sqg_static%lats(1)
top_lat   = sqg_static%lats(2*lmax)
delta_lon = sqg_static%lons(2) - sqg_static%lons(1)
delta_lat = sqg_static%lats(2) - sqg_static%lats(1)

! Compute bracketing lon indices
if (lon >= bot_lon .and. lon <= top_lon ) then
   lon_below = int((lon - bot_lon) / delta_lon) + 1
   lon_above = lon_below + 1
   lon_fract = (lon - ((lon_below - 1) * delta_lon + bot_lon)) / delta_lon
else
   ! At wraparound point
   lon_below = 2*kmax
   lon_above = 1
   if(lon < bot_lon) then
      temp_lon = lon + 360.0_r8
   else
      temp_lon = lon
   endif
   lon_fract = (temp_lon - top_lon) / delta_lon
endif

! Next, compute bracketing lat indices
if (lat >= bot_lat .and. lat <= top_lat) then
   lat_below = int((lat - bot_lat) / delta_lat) + 1
   lat_above = lat_below + 1
   lat_fract = (lat - ((lat_below - 1) * delta_lat + bot_lat)) / delta_lat
else
   ! Outside the channel, return with an error status
   istatus = 1
   obs_val = MISSING_R8
   return
endif

! Get values at the 4 surrounding points
val(1,1) = get_val(x, lon_below, lat_below, level)
val(1,2) = get_val(x, lon_below, lat_above, level)
val(2,1) = get_val(x, lon_above, lat_below, level)
val(2,2) = get_val(x, lon_above, lat_above, level)

! Do the weighted average for interpolation
if ( debug ) write(*,*) 'fracts ', lon_fract, lat_fract
do i = 1, 2
   a(i) = lon_fract * val(2, i) + (1.0_r8 - lon_fract) * val(1, i)
end do

obs_val = lat_fract * a(2) + (1.0_r8 - lat_fract) * a(1)

end subroutine model_interpolate

!========================================================================

function get_val(x, lon_index, lat_index, level)
!-----------------------------------------------
! Given the state vector and location, returns the value at that location

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: lon_index, lat_index, level
real(r8)              :: get_val

integer :: indx

indx = (level-1)*(2*lmax)*(2*kmax) + (lat_index-1)*(2*kmax) + lon_index

!! should not be possible now; but this error check can be commented back in.
!! (it is out for performance reasons, but if you get any strange values, this
!! is a good first check to re-enable.)
!if (indx < 1 .or. indx > size(x)) then
!   write(msgstring1,*)'index ',indx,' not between 1 and ', size(x), ' (should not be possible)'
!   call error_handler(E_ERR,'model_mod:get_val', msgstring1, source, revision, revdate)
!endif

get_val = x(indx)

end function get_val

!========================================================================

function get_model_time_step()
!-----------------------------
! Returns the time step of a single model advance, 
! NOT the dynamical timestep of the model.

type(time_type) :: get_model_time_step

get_model_time_step = model_stride

end function get_model_time_step

!========================================================================

subroutine get_state_meta_data(index_in, location, var_type)
!-----------------------------------------------------------
! Given an integer index into the state vector structure, returns the
! associated location.

integer,             intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: var_type

integer  :: lev_index, lat_index, lon_index
real(r8) :: lon, lat, lev
integer  :: location_index

! avoid out-of-range queries
if ( index_in > model_size ) then
   write(msgstring1,*)'index_in ',index_in,' must be between 1 and ', model_size
   call error_handler(E_ERR,'model_mod:get_state_meta_data', msgstring1, source, revision, revdate)
endif

lev_index =  (index_in-1) / ((2*lmax)*(2*kmax)) + 1
lat_index = ((index_in-1) - ((lev_index-1)*(2*lmax)*(2*kmax))) / (2*kmax) + 1
lon_index =  (index_in-1) - ((lev_index-1)*(2*lmax)*(2*kmax)) - ((lat_index-1)*(2*kmax)) + 1

lev = sqg_static%levs(lev_index)
lat = sqg_static%lats(lat_index)
lon = sqg_static%lons(lon_index)

! With the threed_sphere location module ... you specify that the 
! vertical coordinate is a 'level' by 'which_vert' == 1
location = set_location(lon, lat, lev, 1)

! Alternately, use state_loc that is defined in static_init_model()
!location = state_loc(index_in)

if (present(var_type)) var_type = QTY_POTENTIAL_TEMPERATURE

end subroutine get_state_meta_data

!========================================================================

subroutine end_model()
!---------------------
! Does any shutdown and clean-up needed for model. 

! good practice ... deallocate stuff from static_init_model
deallocate(sqg_static%lons) ; deallocate(sqg_static%lats) ; deallocate(sqg_static%levs)

deallocate(sqg_static%thbB ) ; deallocate(sqg_static%thbT )
deallocate(sqg_static%thbyB) ; deallocate(sqg_static%thbyT)
deallocate(sqg_static%ulinB) ; deallocate(sqg_static%ulinT)

deallocate(sqg_static%hx) ; deallocate(sqg_static%hy)
deallocate(sqg_static%hu) ; deallocate(sqg_static%hv)

deallocate(state_loc)

end subroutine end_model

!========================================================================

function nc_write_model_atts( ncFileID ) result (ierr)
!-----------------------------------------------------
! Writes the model-specific attributes to a netCDF file.

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

integer :: latDimID, lonDimID, levDimID, StateVarDimID, MemberDimID, TimeDimID
integer :: latVarID, lonVarID, levVarID, StateVarVarID, StateVarID

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1
character(len=NF90_MAX_NAME) :: filename

integer :: i

ierr = -1   ! assume things go poorly

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

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                     "nc_write_model_atts", "inquire, " // trim(filename))
call nc_check(nf90_redef(ncFileID), "nc_write_model_atts", "redef, " // trim(filename))

!--------------------------------------------------------------------
! Determine ID's from stuff already in the netCDF file
! make sure time is unlimited dimid
!--------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID), &
                            "nc_write_model_atts", "inq_dimid copy, " // trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid= TimeDimID), &
                            "nc_write_model_atts", "inq_dimid time," // trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(msgstring1,*)"Time Dimension ID ",TimeDimID, &
                     " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", msgstring1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                           len=model_size, dimid=StateVarDimID), &
                           "nc_write_model_atts", "def_dim state, " // trim(filename))

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date" ,str1), &
                          "nc_write_model_atts", "put_att creation_date")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source"  ,source), &
                          "nc_write_model_atts", "put_att model_source")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision), &
                          "nc_write_model_atts", "put_att model_revision")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate" ,revdate), &
                          "nc_write_model_atts", "put_att model_revdate")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","sqg"), &
                          "nc_write_model_atts", "put_att model")

!-------------------------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and one needs
! to know the geometry, the output variables, etc. 
! Both ways are implemented here
!-------------------------------------------------------------------------------

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create a variable for the state vector
   !----------------------------------------------------------------------------

  ! Define the state vector coordinate variable and some attributes.
   call nc_check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=NF90_INT, &
                              dimids=StateVarDimID, varid=StateVarID), &
                             "nc_write_model_atts", "def_var StateVariable")
   call nc_check(nf90_put_att(ncFileID, StateVarID,"long_name","State Variable ID"), &
                             "nc_write_model_atts", "put_att StateVariable long_name")
   call nc_check(nf90_put_att(ncFileID, StateVarID, "units",     "indexical"), &
                             "nc_write_model_atts", "put_att StateVariable units")
   call nc_check(nf90_put_att(ncFileID, StateVarID, "valid_range", (/ 1, model_size /)), &
                             "nc_write_model_atts", "put_att StateVariable valid_range")

   ! Define the actual (3D) state vector, which gets filled as time goes on ... 
   call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=NF90_REAL, &
                 dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
                 varid=StateVarVarID), "nc_write_model_atts", "def_var state")
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "model state or fcopy"), &
                             "nc_write_model_atts", "put_att state long_name")

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncFileID),"nc_write_model_atts", "state_vector enddef")

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarID, (/ (i,i=1,model_size) /)), &
                                    "nc_write_model_atts", "put_var state")

else

   !--------------------------------------------------------------------
   ! Define the new dimensions IDs
   !--------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncFileID, name="lat", len = 2*lmax, dimid = latDimID), &
                 "nc_write_model_atts",'def_dim lat, '// trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name="lon", len = 2*kmax, dimid = lonDimID), &
                 "nc_write_model_atts",'def_dim lon, '// trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name="lev", len = 2,      dimid = levDimID), &
                 "nc_write_model_atts",'def_dim lev, '// trim(filename))

   !--------------------------------------------------------------------
   ! Create the (empty) Variables and the Attributes
   !--------------------------------------------------------------------

   call nc_check(nf90_def_var(ncFileID, name="lon", xtype=nf90_double, dimids=lonDimID, varid=lonVarID), &
                 "nc_write_model_atts",'def_var lon, '//trim(filename) )

   call nc_check(nf90_put_att(ncFileID, lonVarID, "long_name" , "longitude"), &
                 "nc_write_model_atts",'put_att lon:long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, lonVarID, "cartesian_axis", "X"), &
                 "nc_write_model_atts",'put_att lon:cartesian_axis, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, lonVarID, "units", "degrees_east"), &
                 "nc_write_model_atts",'put_att lon:units, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, lonVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)), &
                 "nc_write_model_atts",'put_att lon:valid_range, '//trim(filename))


   call nc_check(nf90_def_var(ncFileID, name="lat", xtype=nf90_double, dimids=latDimID, varid=latVarID), &
                 "nc_write_model_atts",'def_var lat, '//trim(filename) )

   call nc_check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude"), &
                 "nc_write_model_atts",'put_att lat:long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, latVarID, "cartesian_axis", "Y"), &
                 "nc_write_model_atts",'put_att lat:cartesian_axis, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north"), &
                 "nc_write_model_atts",'put_att lat:units, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, latVarID, "valid_range", (/  sqg_static%lats(1), sqg_static%lats(2*lmax) /)), &
                 "nc_write_model_atts",'put_att lat:valid_range, '//trim(filename))


   call nc_check(nf90_def_var(ncFileID, name="lev", xtype=nf90_int, dimids=levDimID, varid=levVarID), &
                 "nc_write_model_atts",'def_var lev'//trim(filename) )

   call nc_check(nf90_put_att(ncFileID, levVarID, "long_name", "level"), &
                 "nc_write_model_atts",'put_att lev:long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, levVarID, "cartesian_axis", "Z"), &
                 "nc_write_model_atts",'put_att lev:cartesian_axis, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, levVarID, "units", " "), &
                 "nc_write_model_atts",'put_att lev:units, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, levVarID, "valid_range", (/  1, 2 /)), &
                 "nc_write_model_atts",'put_att lev:valid_range, '//trim(filename))

   
   call nc_check(nf90_def_var(ncid=ncFileID, name="theta", xtype=nf90_real, &
             dimids = (/ lonDimID, latDimID, levDimID, MemberDimID, unlimitedDimID /), &
             varid  = StateVarVarID),"nc_write_model_atts",'def_var theta, '//trim(filename))

   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "potential temperature"), &
                 "nc_write_model_atts",'put_att theta:long_name, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "units", "K"), &
                 "nc_write_model_atts",'put_att theta:units, '//trim(filename)) 
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "_FillValue", NF90_FILL_REAL), &
                 "nc_write_model_atts", 'put_att theta:FillValue, '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "missing_value", NF90_FILL_REAL), &
                 "nc_write_model_atts",'put_att theta:missing_value, '//trim(filename))

   !--------------------------------------------------------------------
   ! End the definition mode
   !--------------------------------------------------------------------

   call nc_check(nf90_enddef(ncfileID), "nc_write_model_atts", "prognostic enddef")
   
   !--------------------------------------------------------------------
   ! Fill the variables
   !--------------------------------------------------------------------

   call nc_check(nf90_put_var(ncFileID, lonVarID, sqg_static%lons(:) ), &
                 "nc_write_model_atts",'put_var lons, '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, latVarID, sqg_static%lats(:) ), &
                 "nc_write_model_atts",'put_var lats, '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, levVarID, (/ (i,i=1,2) /)    ), &
                 "nc_write_model_atts",'put_var lev, '//trim(filename))

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID),"nc_write_model_atts", "sync, " // trim(filename))

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts

!========================================================================

function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!-------------------------------------------------------------------------------------
! Writes the model variables to a netCDF file.

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID, StateVarVarID
real, allocatable, dimension(:,:,:) :: StateVar
character(len=128) :: filename

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
!-------------------------------------------------------------------------------

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                          "nc_write_model_vars", "inquire, " // trim(filename))

if ( output_state_vector ) then

   call nc_check(nf90_inq_varid(ncFileID, "state", StateVarVarID), &
                               "nc_write_model_vars", "inq_varid state, "  // trim(filename))
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, statevec,  &
                              start=(/ 1, copyindex, timeindex /)), &
                             "nc_write_model_vars", "put_var state, " // trim(filename))                   

else

   !--------------------------------------------------------------------
   ! unpack the state vector into prognostic variables
   !--------------------------------------------------------------------

   allocate(StateVar(2*kmax,2*lmax,2))

   call dart_to_sqg(statevec, StateVar)

   call nc_check(nf90_inq_varid(ncFileID,  "theta",  StateVarVarID), &
          "nc_write_model_vars", 'inq_varid theta, '//trim(filename))
   call nc_check(nf90_put_var( ncFileID,  StateVarVarId, StateVar(:, :, :), &
                               start=(/ 1, 1, 1, copyindex, timeindex /) ), &
                  "nc_write_model_vars", 'put_var theta, '//trim(filename))

   deallocate(StateVar)

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), "nc_write_model_vars", "sync, " // trim(filename))

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars

!========================================================================

subroutine pert_model_state(state, pert_state, interf_provided)
!--------------------------------------------------------------
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

pert_state      = state
interf_provided = .false.

end subroutine pert_model_state

!========================================================================

subroutine ens_mean_for_model(ens_mean)
!--------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

end subroutine ens_mean_for_model

!========================================================================

subroutine sqg_to_dart(theta,statevec)
!-------------------------------------
! Takes sQG variables and wraps them into a DART state vector

real,     dimension(:,:,:), intent(in)  :: theta
real(r8), dimension(:),     intent(out) :: statevec

if ( product(shape(theta)) /= product(shape(statevec)) ) then
   write(msgstring1,*)'size mismatch : product(shape(theta) /= product(shape(statevec))'
   write(msgstring2,*)'size mismatch : ', product(shape(theta)), ' /= ', product(shape(statevec))
   call error_handler(E_ERR,'sqg_to_dart',msgstring1,source,revision,revdate,text2=msgstring2)
else
   statevec = reshape(theta,shape(statevec))
endif

end subroutine sqg_to_dart

!========================================================================

subroutine dart_to_sqg(statevec,theta)
!-------------------------------------
! Takes a DART state vector and unwraps it into sQG variables

real(r8), dimension(:),     intent(in)  :: statevec
real,     dimension(:,:,:), intent(out) :: theta

if ( product(shape(statevec)) /= product(shape(theta)) ) then
   write(msgstring1,*)'size mismatch : product(shape(statevec)) /= product(shape(theta))'
   write(msgstring2,*)'size mismatch : ', product(shape(statevec)), ' /= ', product(shape(theta))
   call error_handler(E_ERR,'dart_to_sqg',msgstring1,source,revision,revdate,text2=msgstring2)
else
   theta = reshape(statevec,shape(theta))
endif

end subroutine dart_to_sqg

!========================================================================

function get_model_static_data()
!-------------------------------
! Returns static info for sQG

type(model_static) :: get_model_static_data

get_model_static_data = sqg_static

end function get_model_static_data

!========================================================================

function set_model_stride()
!------------------------------------------------------------------
! Defines the minimum amount of time to advance the model in one 'go'.
! This is NOT the dynamical timestep of the model. It is usually a
! MULTIPLE of the dynamical timestep - since most models stop that way.
!
! If we can advance the model for 6hour chunks, for example - 
!
! Also : All observations +/- half this timestep are assimilated.
! In essence, this defines the minimum window used for assimilation.

type(time_type) :: set_model_stride

! Check the user input
if ((assimilation_period_seconds < 0) .or. (assimilation_period_days < 0)) then
   write(msgstring1,*)'model_nml:assimilation_period_[seconds,days] must both be positive.'
   write(msgstring2,*)'they are : ', assimilation_period_seconds, assimilation_period_days
   call error_handler(E_ERR,'set_model_stride',msgstring1,source,revision,revdate,text2=msgstring2)
elseif ((assimilation_period_seconds == 0) .and. (assimilation_period_days == 0)) then
   write(msgstring1,*)'at least one of model_nml:assimilation_period_[seconds,days] must be positive.'
   write(msgstring2,*)'they are : ', assimilation_period_seconds, assimilation_period_days
   call error_handler(E_ERR,'set_model_stride',msgstring1,source,revision,revdate,text2=msgstring2)
else
   ! FIXME ... ensure that stride is a multiple of 'time_step'
   set_model_stride = set_time(assimilation_period_seconds, assimilation_period_days)
endif

end function set_model_stride

!========================================================================
! End of model_mod
!========================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

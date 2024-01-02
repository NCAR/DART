! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module obs_mod

! Currently uses NAG for random numbers.

use        types_mod, only : r8
use        model_mod, only : lat_max, num_lon, lon, lat, get_model_size
use     nag_wrap_mod, only : g05ddf_wrap
use    obs_tools_mod, only : conv_state_to_obs, obs_def_type, def_single_obs
!!!use transforms_mod
use loc_and_dist_mod, only : loc_type, set_loc, get_loc
use    utilities_mod, only : file_exist, open_file, close_file,            &
                             register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, find_namelist_in_file,           &
                             check_namelist_read, do_output

implicit none
private

public :: num_obs, obs_var, take_obs, ens_ics, state_to_obs, &
       init_obs, take_single_obs, get_close_state, obs_loc

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer :: num_obs = 0

! Following is to allow initialization of obs_def_type
logical :: obs_init = .false.

type (obs_def_type), allocatable :: obs_def(:)

! Array of structure for observation locations; static with time in this version
type(loc_type), allocatable :: obs_loc(:)

! Storage for the observational variance
real(r8), allocatable :: obs_variance(:)

!=======================================================================

!---- namelist with default values
! Set a cut-off for lon and lat for close obs search

real(r8) :: close_lat_window = 10.0_r8
real(r8) :: close_lon_window = 10.0_r8

namelist /obs_nml/ close_lat_window, close_lon_window

!--- module name and version number

character(len = 18), parameter :: module_name = 'obs:spectral barot'
character(len = 12), parameter :: vers_num = '10/04/2000'

!=======================================================================

contains



  subroutine init_obs
!=======================================================================
! subroutine init_obs
!
! Initializes the description a linear observations operator. For each
! observation, a list of the state variables on which it depends and the
! coefficients for each state variable is passed to def_single_obs
! which establishes appropriate data structures.

implicit none

integer  :: i, j, lon_ind_lo, lon_ind_hi, lat_ind_lo, lat_ind_hi 
integer  :: state_index(4), iunit, io
real(r8) :: olon, olat, coef(4), temp_lon, frac_lon, frac_lat
real(r8) :: temp_frac, numerator, denominator

character(len=129) :: err_string, nml_string

call register_module(source,revision,revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_nml", iunit)
read(iunit, nml = obs_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_nml")

! Record the namelist values used for the run ...
if (do_output()) then
   write(nmlfileunit, nml=obs_nml)
   write(     *     , nml=obs_nml)
endif

! Initialization for identity observations

obs_init = .true.
write(*, *) 'close_lat_window, close_lon_window ', close_lat_window, close_lon_window

! Read in the forced_barot_obs_def file

iunit = open_file('forced_barot_obs_def', action = 'read')
read(iunit, *) num_obs
allocate(obs_def(num_obs), obs_loc(num_obs), obs_variance(num_obs))

do i = 1, num_obs
   read(iunit, *) olon, olat, obs_variance(i)
   if( olon <   0.0_r8 .or. olon > 360.0_r8 .or. &
       olat < -90.0_r8 .or. olat >  90.0_r8) then
      write(*, *) 'illegal obs location read in init_obs for forced_barot'
      write(*, *) 'lon lat ', olon, olat
      stop
   end if

   call set_loc(obs_loc(i), olon, olat)
end do
call close_file(iunit)

! Four point interpolation for points wherever

do i = 1, num_obs

   call get_loc(obs_loc(i), olon, olat)

   ! Do four point bilinear interpolation
   ! Begin by finding lon bounding points
   ! First case, between first and last lon on wraparound

   if(olon < lon(1) .or. olon > lon(num_lon)) then
      lon_ind_lo  = num_lon
      lon_ind_hi  = 1
      numerator   = abs(olon - lon(num_lon))
      denominator = abs(lon(1) - lon(num_lon))
      if(  numerator > 180.0_r8) numerator   = numerator   - 180.0_r8
      if(denominator > 180.0_r8) denominator = denominator - 180.0_r8
      frac_lon = numerator / denominator
   else 
      do j = 1, num_lon - 1
         if(olon >= lon(j) .and. olon <= lon(j + 1)) then
            lon_ind_lo = j
            lon_ind_hi = j + 1
            frac_lon = (olon - lon(j)) / (lon(j + 1) - lon(j))
            goto 21
         endif
      end do
      write(*,*)'ERROR(init_obs): fell off end of longitude search'
      write(*,*)'olon is ', olon
      stop
   endif

21 continue

   ! Next find out where lat of obs is; die on extrapolation for simplicity

   if(olat < lat(1) .or. olat > lat(lat_max)) then
      write(*,*)'latitude of obs must not be outside of model grid'
      write(*,*)'This could be fixed but its some work', lat(1), lat(lat_max)
      stop
   else
      do j = 1, lat_max - 1
         if(olat >= lat(j) .and. olat <= lat(j + 1)) then
            lat_ind_lo = j
            lat_ind_hi = j + 1
            frac_lat = (olat - lat(j)) / (lat(j + 1) - lat(j))
            goto 31
         endif
      end do
      write(*,*)'ERROR(init_obs): fell off end of latitude search'
      stop
   endif

! Found the bounding points, need to turn this into indices

31 state_index(1) = lat_ind_lo + (lon_ind_lo - 1) * lat_max
   state_index(2) = lat_ind_lo + (lon_ind_hi - 1) * lat_max
   state_index(3) = lat_ind_hi + (lon_ind_lo - 1) * lat_max
   state_index(4) = lat_ind_hi + (lon_ind_hi - 1) * lat_max

   ! Temporary test with coefficients set to 1; need to get weighting

   coef(1) = (1.0 - frac_lon) * (1.0 - frac_lat)
   coef(2) = frac_lon * (1.0 - frac_lat)
   coef(3) = (1.0 - frac_lon) * frac_lat
   coef(4) = frac_lon * frac_lat

   call def_single_obs(4, state_index(1:4), coef(1:4), obs_def(i))

end do

end subroutine init_obs



  function obs_var()
!=======================================================================
! function obs_var()
!
! Defines the observational error variance. Eventually may need to be
! generalized to covariance.

implicit none

real(r8) :: obs_var(num_obs)


obs_var = obs_variance

end function obs_var



  function take_obs(x)
!=======================================================================
! function take_obs(x)
!
! Given a model state, x, returns observations for assimilation.
! For perfect model, take_obs is just state_to_obs

implicit none

real(r8), intent(in) :: x(:)
real(r8)             :: take_obs(num_obs)

! Important to initialize obs structure before taking obs

if(.not. obs_init) call init_obs

take_obs = conv_state_to_obs(x, obs_def, num_obs)

end function take_obs



  function take_single_obs(x, index)
!========================================================================
! function take_single_obs(x, index)
!
! Given a model state, x, returns observations for assimilation.
! For perfect model, take_obs is just state_to_obs

implicit none

real(r8), intent(in) :: x(:)
integer,  intent(in) :: index
real(r8)             :: take_single_obs

real(r8) :: take(1)

! Important to initialize obs structure before taking obs

if(.not. obs_init) call init_obs

take = conv_state_to_obs(x, obs_def(index:index), 1)
take_single_obs = take(1)

end function take_single_obs



  function state_to_obs(x)
!=======================================================================
! function state_to_obs(x)
!
! Given a model state, returns the associated 'observations' by applying
! the observational operator to the state.

implicit none

real(r8), intent(in) :: x(:)
real(r8)             :: state_to_obs(num_obs)

! Important to initialize obs structure before taking obs

if(.not. obs_init) call init_obs

state_to_obs = conv_state_to_obs(x, obs_def, num_obs)

end function state_to_obs



  subroutine ens_ics(x, as)
!=======================================================================
! subroutine ens_ics(x, as)
!
! Get initial state for ensemble assimilation

implicit none

real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: as(:, :)

integer :: i, j

!WARNING: MANY CHANGES

do i = 1, size(x)
   do j = 1, size(as, 2)
!      as(i, j) = x(i) + 1e5 * g05ddf_wrap(dble(0.0), dble(1.0))
       as(i, j) = x(i) + 1e6 * g05ddf_wrap(dble(0.0), dble(1.0))
   end do
end do

end subroutine ens_ics



  subroutine get_close_state(obs_num, list, max_list, num)
!========================================================================
! subroutine get_close_state(obs_num, list, max_list, num)

implicit none

integer, intent(in)  :: obs_num, max_list
integer, intent(out) :: list(max_list), num

integer  :: i, j, low_lat_ind, hi_lat_ind, lon_ind, low_lon_ind, hi_lon_ind
real(r8) :: olon, olat, low_lat, hi_lat, low_lon, hi_lon
real(r8) :: rlow_lon_ind, rhi_lon_ind, lon_width

! Important to initialize obs structure before taking obs

if(.not. obs_init) call init_obs

! Get lon and lat for this observation

call get_loc(obs_loc(obs_num), olon, olat)

! NOTE THAT THIS IS NOT ALLOWING OBS TO INFLUENCE OVER POLE!!!
! Get lat window for state

low_lat = olat - close_lat_window
hi_lat  = olat + close_lat_window

! Get corresponding indices; search is a slow way, may need to speed up???

do i = 1, lat_max
   if(low_lat < lat(i)) then
      low_lat_ind = i
      goto 11
   end if
end do

low_lat_ind = lat_max

11 do i = low_lat_ind, lat_max
   if(hi_lat < lat(i)) then
      hi_lat_ind = i - 1
      goto 21
   endif
end do
hi_lat_ind = lat_max

21 continue

! Next loop through to get longitudes for each latitude row

lon_width = close_lon_window / cos(3.14159_r8 * olat / 180.0_r8)
low_lon   = olon - lon_width

if(low_lon < 0.0_r8) low_lon = low_lon + 360.0_r8

rlow_lon_ind = (low_lon / 360.0_r8) * num_lon

! Low lon index, add 1 since index 1 is at 0 degrees

low_lon_ind = ceiling(rlow_lon_ind) + 1
hi_lon      = low_lon + 2.0_r8 * lon_width
rhi_lon_ind = (hi_lon / 360.0_r8) * num_lon
hi_lon_ind  = floor(rhi_lon_ind) + 1

num = 0
do j = low_lat_ind, hi_lat_ind
   do lon_ind = low_lon_ind, hi_lon_ind
      i = lon_ind
      if(i > num_lon) i = i - num_lon 
      num = num + 1
      if(num > max_list) then
         write(*, *) 'num_close_state > max_list in get_close_state'
         stop
      end if
      list(num) = j + (i - 1) * lat_max 
   end do
end do   

end subroutine get_close_state

!===========================================================================
! End of obs_mod.f90
!===========================================================================

end module obs_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

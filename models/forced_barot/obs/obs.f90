module obs_mod

! Currently uses NAG for random numbers.

use model_mod, only : lat_max, num_lon, lon, lat, get_model_size
use nag_wrap_mod, only : g05ddf_wrap
use obs_tools_mod, only : conv_state_to_obs, obs_def_type, def_single_obs
!!!use transforms_mod
use loc_and_dist_mod, only : loc_type, set_loc, get_loc
use utilities_mod, only : file_exist, open_file, check_nml_error, &
   print_version_number, close_file

private
public num_obs, obs_var, take_obs, ens_ics, state_to_obs, &
   init_obs, take_single_obs, get_close_state, obs_loc

integer :: num_obs = 0

! Following is to allow initialization of obs_def_type
logical :: obs_init = .false.

type (obs_def_type), allocatable :: obs_def(:)

! Array of structure for observation locations; static with time in this version
type(loc_type), allocatable :: obs_loc(:)

! Storage for the observational variance
double precision, allocatable :: obs_variance(:)

!=======================================================================

!---- namelist with default values
! Set a cut-off for lon and lat for close obs search
double precision :: close_lat_window = 10.0
double precision :: close_lon_window = 10.0

namelist /obs_nml/ close_lat_window, close_lon_window

!--- module name and version number
character(len = 18), parameter :: module_name = 'obs:spectral barot'
character(len = 12), parameter :: vers_num = '10/04/2000'

!=======================================================================

contains

!=======================================================================

subroutine init_obs

! Initializes the description a linear observations operator. For each
! observation, a list of the state variables on which it depends and the
! coefficients for each state variable is passed to def_single_obs
! which establishes appropriate data structures.

implicit none

integer :: i, j, unit_num, lon_ind_lo, lon_ind_hi, lat_ind_lo, lat_ind_hi 
integer :: state_index(4), unit, ierr, io
double precision :: olon, olat, coef(4), temp_lon, frac_lon, frac_lat
double precision :: temp_frac, numerator, denominator

! Read namelist for run time control
if(file_exist('input.nml')) then
   unit = open_file(file = 'input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(unit, nml = obs_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'obs_nml')
   enddo
 11 continue
   call close_file(unit)
endif

! Write the namelist to a log file
unit = open_file(file = 'logfile.out', action = 'append')
call print_version_number(unit, module_name, vers_num)
write(unit, nml = obs_nml)
call close_file(unit)


! Initialization for identity observations
obs_init = .true.
write(*, *) 'close_lat_window, close_lon_window ', close_lat_window, close_lon_window

! Read in the forced_barot_obs_def file
unit_num = open_file(file = 'forced_barot_obs_def', action = 'read')
read(unit_num, *) num_obs
allocate(obs_def(num_obs), obs_loc(num_obs), obs_variance(num_obs))
do i = 1, num_obs
   read(unit_num, *) olon, olat, obs_variance(i)
   if(olon < 0.0 .or. olon .gt. 360.0 .or. olat < -90.0 .or. olat > 90.0) then
      write(*, *) 'illegal obs location read in init_obs for forced_barot'
      write(*, *) 'lon lat ', olon, olat
      stop
   end if

   call set_loc(obs_loc(i), olon, olat)
end do
call close_file(unit = unit_num)


! Four point interpolation for points wherever
do i = 1, num_obs
   call get_loc(obs_loc(i), olon, olat)
! Do four point bilinear interpolation
! Begin by finding lon bounding points
! First case, between first an last lon on wraparound
   if(olon < lon(1) .or. olon > lon(num_lon)) then
      lon_ind_lo = num_lon
      lon_ind_hi = 1
      numerator = dabs(olon - lon(num_lon))
      if(numerator > 180.0) numerator = numerator - 180.0
      denominator = dabs(lon(1) - lon(num_lon))
      if(denominator > 180.0) denominator = denominator - 180.0
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
      write(*, *) 'error: fell off end of longitude search in init_obs'
      write(*, *) 'olon is ', olon
      stop
   endif

21 continue
! Next find out where lat of obs is; die on extrapolation for simplicity
   if(olat < lat(1) .or. olat > lat(lat_max)) then
      write(*, *) 'latitude of obs must not be outside of model grid'
      write(*, *) 'This could be fixed but its some work', lat(1), lat(lat_max)
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
      write(*, *) 'error: fell off end of latitude search in init_obs'
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

!=======================================================================

function obs_var()

implicit none

double precision :: obs_var(num_obs)

! Defines the observational error variance. Eventually may need to be
! generalized to covariance.

obs_var = obs_variance

end function obs_var

!=======================================================================

function take_obs(x)

implicit none

double precision :: take_obs(num_obs)
double precision, intent(in) :: x(:)

! Important to initialize obs structure before taking obs
if(.not. obs_init) call init_obs

! Given a model state, x, returns observations for assimilation.
! For perfect model, take_obs is just state_to_obs
take_obs = conv_state_to_obs(x, obs_def, num_obs)

end function take_obs

!========================================================================
function take_single_obs(x, index)

implicit none

double precision :: take_single_obs
double precision, intent(in) :: x(:)
integer, intent(in) :: index

double precision :: take(1)

! Important to initialize obs structure before taking obs
if(.not. obs_init) call init_obs

! Given a model state, x, returns observations for assimilation.
! For perfect model, take_obs is just state_to_obs
take = conv_state_to_obs(x, obs_def(index:index), 1)
take_single_obs = take(1)

end function take_single_obs

!=======================================================================

function state_to_obs(x)

implicit none

double precision, intent(in) :: x(:)
double precision :: state_to_obs(num_obs)

! Given a model state, returns the associated 'observations' by applying
! the observational operator to the state.

! Important to initialize obs structure before taking obs
if(.not. obs_init) call init_obs

state_to_obs = conv_state_to_obs(x, obs_def, num_obs)

end function state_to_obs

!=======================================================================


subroutine ens_ics(x, as)

implicit none

!  Get initial state for ensemble assimilation

double precision, intent(in) :: x(:)
double precision, intent(out) :: as(:, :)
integer :: i, j


!WARNING: MANY CHANGES
do i = 1, size(x)
   do j = 1, size(as, 2)
!       as(i, j) = x(i) + 1e5 * g05ddf_wrap(dble(0.0), dble(1.0))
       as(i, j) = x(i) + 1e6 * g05ddf_wrap(dble(0.0), dble(1.0))
   end do
end do

end subroutine ens_ics

!========================================================================

subroutine get_close_state(obs_num, list, max_list, num)

implicit none

integer, intent(in) :: obs_num, max_list
integer, intent(out) :: list(max_list), num

integer :: i, j, low_lat_ind, hi_lat_ind, lon_ind, low_lon_ind, hi_lon_ind
double precision :: olon, olat, low_lat, hi_lat, low_lon, hi_lon
double precision :: rlow_lon_ind, rhi_lon_ind, lon_width

! Important to initialize obs structure before taking obs
if(.not. obs_init) call init_obs

! Get lon and lat for this observation
call get_loc(obs_loc(obs_num), olon, olat)

! NOTE THAT THIS IS NOT ALLOWING OBS TO INFLUENCE OVER POLE!!!
! Get lat window for state
low_lat = olat - close_lat_window
hi_lat = olat + close_lat_window
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
lon_width = close_lon_window / cos(2. * 3.14159 * olat / 360.0)
low_lon = olon - lon_width
if(low_lon < 0.0) low_lon = low_lon + 360.0
rlow_lon_ind = (low_lon / 360.0) * num_lon
! Low lon index, add 1 since index 1 is at 0 degrees
low_lon_ind = ceiling(rlow_lon_ind) + 1
hi_lon = low_lon + 2. * lon_width
rhi_lon_ind = (hi_lon / 360.0) * num_lon
hi_lon_ind = floor(rhi_lon_ind) + 1

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

end module obs_mod

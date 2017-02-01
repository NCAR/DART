! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module obs_mod

! This version of the spectral barotropic obs model reads in standard
! psi field data from some gdas files. Just reads in one field each time
! it's called.

use model_mod, only : lat_max, num_lon, location_type, dp_to_grid, lon, lat, &
   num_fourier, num_spherical, barot_to_dp
use random_sequence_mod, only : init_random, random_gaussian

private
public :: num_obs, obs_var, take_obs, ens_ics, obs_location, state_to_obs

! Added to cheat on close_obs for spectral models; 12/1/99
public :: num_x_obs, num_y_obs, ob2_lon, ob2_lat

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer, parameter :: num_x_obs = 40, num_y_obs = 32
integer, parameter :: num_obs = num_x_obs * num_y_obs

! Global storage for obs locations
double precision obs_lon(num_obs), obs_lat(num_obs)
double precision ob2_lon(num_x_obs, num_y_obs), ob2_lat(num_x_obs, num_y_obs)

contains

!-----------------------------------------------------------------------

function obs_var()

implicit none

double precision :: obs_var(num_obs)
integer :: i

! Standard for most oi tests, 26 March, 1999
!obs_var = (1e5)**2


obs_var = (1e6)**2

! Following used for first successful OI tests, 24 Jan, 99
!obs_var = (1e3)**2

!obs_var = 0.0

!do i = 1, num_obs, 10
!   obs_var(i) = (5e5)**2
!end do

end function obs_var

!-----------------------------------------------------------------------

function obs_location()

implicit none

type(location_type) :: obs_location(num_obs)

integer :: i, j, num

! Get locations of equally space obs
num = 0
do i = 1, num_x_obs
   do j = 1, num_y_obs
      num = num + 1
      obs_lon(num) = (i - 1) * (360.0 / num_x_obs) + 360.0 / (2. * num_obs) 
      obs_lat(num) = (j - 1)*(180.0 / num_y_obs) + 180.0/(num_y_obs * 2.) - 90.0
! Can fix singular inversion problems by moving observations away from the pole
!      obs_lat(num) = -80.0 + (j - 1) * (160.0 / num_y_obs) + 160.0 / (num_y_obs * 2.)
      ob2_lon(i, j) = obs_lon(num)
      ob2_lat(i, j) = obs_lat(num)
      obs_location(num)%lon = obs_lon(num)
      obs_location(num)%lat = obs_lat(num)
   end do
end do

end function obs_location

!-----------------------------------------------------------------------

function take_obs(x)

! For the observed data input, the current state x is just ignored.
! For observed data, x is returned as the state space state

implicit none

double precision :: take_obs(num_obs)
double precision, intent(inout) :: x(:)
integer :: i, j, m, n, num, x_lo, x_hi, y_lo, y_hi
real :: rlat(lat_max), rlon(num_lon)
double precision :: x_grid(num_lon, lat_max)
double precision :: x_frac, y_frac

complex, dimension(0:num_fourier, 0:num_spherical) :: psisp

write(*, *) 'Working in realistic take_obs with T21 observations'
!if(1 == 1) then
!   write(*, *) 'Need to take out added noise after take_obs in main program'
!   write(*, *) 'to use this version of observed real_obs'
!   stop
!endif

! Let's read in one of the old format files to act as an initial condition
! FIRST WARNING: USING FIXED FILE NUMBER; SHOULD USE UTILITIES
! First zero the model precision spectral (must be greater than T21),
! Then read in the observed T21 state
psisp = 0.0

! Unit 81 is opened in init_conditions for now; TOO UGLY TO KEEP
!open(unit = 81, file = '/home/jla/gdas/nov_to_mar')

! WARNING: Remember that old model uses TOTAL wavenumber
psisp = 0.0
do n = 0, 21
   do m = 0, n
 31   format(1x, 2(e10.4, 1x))
      read(81, 31) psisp(m, n - m)
   end do
end do
write(*, *) 'Observed 4, 4 is ', psisp(4, 4)

! Next need to convert this to current model resolution single dimension state
x_grid = dble(dp_to_grid(barot_to_dp(psisp)))

! RETURN X as FULL STATE SPACE STATE
x = barot_to_dp(psisp)

! Loop through the set of observed points and do simple linear interpolation
! from model grid

do i = 1, num_obs
! Regular grid in longitude
   x_lo = int(obs_lon(i) / (360.0 / num_lon)) + 1
   if(x_lo /= num_lon) then
      x_hi = x_lo + 1
      x_frac = (obs_lon(i) - lon(x_lo)) / (lon(x_hi) - lon(x_lo))
   else
      x_hi = 1
      x_frac = (obs_lon(i) - lon(x_lo)) / (360.0 - lon(x_lo))
   end if
! Interpolate into latitudes (Gaussian); for now, just cheat for lats
! above or below start of grid to save programming junk
   if(obs_lat(i) < lat(1)) then
      y_lo = 1
      y_hi = 2
      y_frac = 0.0
   else if(obs_lat(i) > lat(lat_max)) then
      y_lo = lat_max - 1
      y_hi = lat_max
      y_frac = 1.0
   else
      do j = 1, lat_max - 1
         if(obs_lat(i) > lat(j) .and. obs_lat(i) < lat(j + 1)) then
            y_lo = j
            y_hi = j + 1
            y_frac = (obs_lat(i) - lat(y_lo)) / (lat(y_hi) - lat(y_lo))
         end if
      end do
   endif
!  Now do interpolation
   take_obs(i) = x_frac * y_frac * x_grid(x_hi, y_hi) + &
      (1. - x_frac) * (1. - y_frac) * x_grid(x_lo, y_lo) + &
      (1. - x_frac) * y_frac * x_grid(x_lo, y_hi) + &
      x_frac * (1. - y_frac) * x_grid(x_hi, y_lo)
! 11 format(1x, 6(e10.4, 2x))
end do

end function take_obs

!-----------------------------------------------------------------------

function state_to_obs(x)

! This is function h that maps state variables to obs

implicit none

double precision :: state_to_obs(num_obs)
double precision, intent(in) :: x(:)
integer :: i, j, m, n, num, x_lo, x_hi, y_lo, y_hi
real :: rlat(lat_max), rlon(num_lon)
double precision :: x_grid(num_lon, lat_max)
double precision :: x_frac, y_frac

! Get the gridded real data
x_grid = dble(dp_to_grid(x))

! Loop through the set of observed points and do simple linear interpolation
! from model grid

do i = 1, num_obs
! Regular grid in longitude
   x_lo = int(obs_lon(i) / (360.0 / num_lon)) + 1
   if(x_lo /= num_lon) then
      x_hi = x_lo + 1
      x_frac = (obs_lon(i) - lon(x_lo)) / (lon(x_hi) - lon(x_lo))
   else
      x_hi = 1
      x_frac = (obs_lon(i) - lon(x_lo)) / (360.0 - lon(x_lo))
   end if
! Interpolate into latitudes (Gaussian); for now, just cheat for lats
! above or below start of grid to save programming junk
   if(obs_lat(i) < lat(1)) then
      y_lo = 1
      y_hi = 2
      y_frac = 0.0
   else if(obs_lat(i) > lat(lat_max)) then
      y_lo = lat_max - 1
      y_hi = lat_max
      y_frac = 1.0
   else
      do j = 1, lat_max - 1
         if(obs_lat(i) > lat(j) .and. obs_lat(i) < lat(j + 1)) then
            y_lo = j
            y_hi = j + 1
            y_frac = (obs_lat(i) - lat(y_lo)) / (lat(y_hi) - lat(y_lo))
         end if
      end do
   endif
!  Now do interpolation
   state_to_obs(i) = x_frac * y_frac * x_grid(x_hi, y_hi) + &
      (1. - x_frac) * (1. - y_frac) * x_grid(x_lo, y_lo) + &
      (1. - x_frac) * y_frac * x_grid(x_lo, y_hi) + &
      x_frac * (1. - y_frac) * x_grid(x_hi, y_lo)
! 11 format(1x, 6(e10.4, 2x))
end do

end function state_to_obs

!-----------------------------------------------------------------------

subroutine ens_ics(x, as)

implicit none

!  Get initial state for ensemble assimilation

double precision, intent(in) :: x(:)
double precision, intent(out) :: as(:, :)
integer :: i, j
type(random_type), save :: r
logical, save :: first = .true.

if (first) then
   call init_random(r)
   first = .false.
endif

do i = 1, size(x)
   do j = 1, size(as, 2)
       as(i, j) = x(i) + 2e5 * random_gaussian(r, 0.0, 1.0)
!      write(*, *) 'as ', i, j, x(i), as(i, j)
   end do
end do

end subroutine ens_ics

!------------------------------------------------------------------------

end module obs_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

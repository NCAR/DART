! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

use netcdf
use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time
use     location_mod, only : location_type, get_dist, set_location, get_location, &
                             LocationDims, LocationName, LocationLName, &
                             get_close_maxdist_init, get_close_obs_init, get_close_obs

use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, do_output, &
                             nmlfileunit, find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term

use     obs_kind_mod, only : QTY_U_WIND_COMPONENT, &
                             QTY_V_WIND_COMPONENT, &
                             QTY_TEMPERATURE, &
                             QTY_SURFACE_PRESSURE, &
                             QTY_VERTICAL_VELOCITY
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

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
          get_close_maxdist_init, get_close_obs_init, get_close_obs, ens_mean_for_model


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Define the random number sequence variable
type(random_seq_type)   :: random_seq

! Basic model parameters controlled by nameslist; have defaults

!---------------------------------------------------------------
! Namelist with default values
!
! model_size             length of the data assimilation control
!                          vector, and is = naz*nrad*nzed*ntype
! naz                    number of gridpoints in the azimuthal direction
! nrad                   number of gridpoints in the radial direction
! nzed                   number of gridpoints in the z direction
! ntype                  number of variables at each grid point that
!                          are impacted by the data assimilation
! daz                    the gridpoint increment in the azimuthal
!                          direction (degrees)
! drad                   the gridpoint increment in the radial
!                          direction (meters)
! dzed                   the gridpoint increment in the z 
!                          direction (meters)
! inner_rad              the radius of the inner cylinder of the
!                          annulus (meters)
! outer_rad              the radius of the outer cylinder of the
!                           annulus (meters)
! depth                  the overall depth of the fluid (meters)
! delta_t                the model integration time step (in seconds)
! time_step_days         the number of days in an integration time step
! time_step_seconds      the number of seconds in an integration time step.
!                        this is not real time, this is time scales by the
!                        tank's rotation rate.  f=0.5 mean a "day" is 4sec.
!                        equate 4sec to 86400sec and determine the dart
!                        time step

integer  :: model_size        = 539400
integer  :: naz               = 120 
integer  :: nrad              = 31
integer  :: nzed              = 29
integer  :: ntype             = 5
real(r8) :: daz               = 3.00_r8
real(r8) :: drad              = 0.01_r8
real(r8) :: dzed              = 0.005_r8
real(r8) :: inner_rad         = 0.08_r8
real(r8) :: outer_rad         = 0.3_r8
real(r8) :: depth             = 0.14_r8
real(r8) :: delta_t           = 0.1_r8
integer  :: time_step_days    = 0
integer  :: time_step_seconds = 2160

namelist /model_nml/ model_size, naz, nrad, nzed, ntype, daz, drad, dzed, &
    inner_rad, outer_rad, depth, delta_t, time_step_days, time_step_seconds
!----------------------------------------------------------------


! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
integer, allocatable             :: state_kind(:)

! Define a type for mapping vector state to 3d state
! time_type is implemented as seconds and days to allow for larger intervals
type grid_type
   private
   integer :: az
   integer :: rad
   integer :: zed
end type grid_type

! Define an array that maps the vector index to the 3d-index

integer, allocatable :: the_grid(:,:,:)

! Define the location of grid variables
real(r8), allocatable :: azimuthal(:), radial(:), height(:)

! Define the time step
type(time_type) :: time_step

! Define whether you want to output the entire control vector, or
! output each prognostic variable separately
logical :: output_state_vector = .false.

! Define a type for holding the MITgcm annulus prognostic variables
type model_type
   private
   real(r8), pointer :: vars_3d(:, :, :, :)
end type model_type

! Define the variable for holding the MITgcm annulus prog vars
type(model_type) :: var

! Define a character array that contains the prognostic variable names
character (len=8), dimension(5) :: flds = &
          (/'U       ','V       ','W       ','T       ','P       '/)

contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
! Initializes information about the grid and about the variables for
! the MITgcm in annulus configuration. 

real(r8) :: x_loc
integer  :: i, iunit, ierr, io
 
real(r8) ::  az,  rad,  zed
integer  :: iaz, irad, ized
integer  :: icount, itype
character(len=129) :: err_string, nml_string

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist to the logfile
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Create space for the model prognostic variables
allocate(var%vars_3d(naz, nrad, nzed, ntype))

! Create space for the grid
allocate(the_grid(naz*ntype,nrad*ntype,nzed*ntype))

!-------------------------------------------------------
! Build a map between the vector state and the 3d state
! that takes the variable type into account
!-------------------------------------------------------
icount = 0
do ized = 1, nzed
   do irad = 1, nrad
      do iaz = 1, naz
         do itype = 1, ntype
            icount = icount + 1
            the_grid(iaz  + (itype - 1)*naz,  &
                     irad + (itype - 1)*nrad, &
                     ized + (itype - 1)*nzed) = icount
         end do
      end do
   end do
end do 

! Create space for grid variables
allocate(azimuthal(naz), radial(nrad), height(nzed))

!---------------------------
! Define the grid locations
!---------------------------

! First the azimuthal angle
do iaz = 1, naz
 azimuthal(iaz) = (iaz - 1)*daz
enddo

! Next the radial direction
do irad = 1, nrad
 radial(irad) = (irad - 1)*drad
end do

! And last the depth (height)
do ized = 1, nzed
 height(ized) = - (ized - 1)*dzed
enddo

! Create space for locations and kinds
allocate(state_loc(model_size))
allocate(state_kind(model_size))

!--------------------------------------------------
! Define the locations of the model state variables.  
!--------------------------------------------------

! For ntype=5, the five state variables altered by the data
! assimilation are U, V, W, T, and P, respectively.
! Because of the C-grid weirdness, I've set up explicit
! loops for each of these variables.
!
! Note that because the MITgcm is an ocean model, zed is
! zero at the top of the fluid and becomes increasingly
! negative as one descends into the fluid.  The azimuthal
! angle increases in an anti-clockwise direction (as viewed
! from the top of the fluid), and the radial direction 
! increases as one moves out from the center of the annulus.

! First the U variable
icount = 0
do iaz = 1, naz
   az = (iaz - 1)*daz
   do irad = 1, nrad
      rad = (irad - 1)*drad + drad/2.00_r8
      do ized = 1, nzed
         zed = - (ized - 1)*dzed - dzed/2.00_r8
         icount = icount + 1
         state_loc(icount)  = set_location(az,rad,zed,3)
         state_kind(icount) = QTY_U_WIND_COMPONENT 
      end do
   end do
end do

! Next the V variable
do iaz = 1, naz
   az = (iaz - 1)*daz + daz/2.00_r8
   do irad = 1, nrad
      rad = (irad - 1)*drad 
      do ized = 1, nzed
         zed = - (ized - 1)*dzed - dzed/2.00_r8
         icount = icount + 1
         state_loc(icount)  = set_location(az,rad,zed,3)
         state_kind(icount) = QTY_V_WIND_COMPONENT
      end do
   end do
end do

! Now the W variable
do iaz = 1, naz
   az = (iaz - 1)*daz + daz/2.00_r8
   do irad = 1, nrad
      rad = (irad - 1)*drad + drad/2.00_r8 
      do ized = 1, nzed
         zed = - (ized - 1)*dzed 
         icount = icount + 1
         state_loc(icount)  = set_location(az,rad,zed,3)
         state_kind(icount) = QTY_VERTICAL_VELOCITY
      end do
   end do
end do

! And finally the T and P variables
do iaz = 1, naz
   az = (iaz - 1)*daz + daz/2.00_r8
   do irad = 1, nrad
      rad = (irad - 1)*drad + drad/2.00_r8
      do ized = 1, nzed
         zed = - (ized - 1)*dzed - dzed/2.00_r8
         icount = icount + 1
         state_loc(icount)  = set_location(az,rad,zed,3)
         state_kind(icount) = QTY_TEMPERATURE
      end do
   end do
end do
do iaz = 1, naz
   az = (iaz - 1)*daz + daz/2.00_r8
   do irad = 1, nrad
      rad = (irad - 1)*drad + drad/2.00_r8
      do ized = 1, nzed
         zed = - (ized - 1)*dzed - dzed/2.00_r8
         icount = icount + 1
         state_loc(icount)  = set_location(az,rad,zed,3)
         state_kind(icount) = QTY_SURFACE_PRESSURE
      end do
   end do
end do

! The time_step in terms of a time type must also be initialized. 
! Note that the model timestep is typically fractions of a second
! and that this is not a timestep supported by DART.  To address
! this problem, the model timestep is scaled by the tank's rotation 
! rate (one tank revolution is defined as one "day"), and the 
! associated time step is determined.  For example, f0=0.5 means
! \Omega = 0.25, or 1 day is 4 seconds.  For a dt of 0.1, a model
! timestep is 1/40th of a day.  There ar 86400 seconds in an earth
! day, so time_step_seconds becomes 1/40th of 86400 = 2160.

time_step = set_time(time_step_seconds, time_step_days)

! Initialize a repeatable random sequence for perturbations
call init_random_seq(random_seq)

end subroutine static_init_model



subroutine comp_dt(z, dt)
!------------------------------------------------------------------
! subroutine comp_dt(z, dt)
! 
! Computes the time tendency of the model for async=0
!
! Not currently relevant to the MITgcm annulus configuration

real(r8), intent( in)        ::  z(:)
real(r8), intent(out)        :: dt(:)

end subroutine comp_dt



subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Initial conditions for the model.  Only necessary for async=0

real(r8), intent(inout) :: x(:)

end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does single time step advance for async=0

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

integer :: get_model_size

get_model_size = model_size

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
!
! Gets the initial time for a state from the model. 
!
! Not clear when this subroutine is relevant

type(time_type), intent(out) :: time

! For now, just set to 0
time = set_time(0, 0)

end subroutine init_time



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Interpolates a value of the state vector x to the specified
! location. 
!
! Will ignore this for the moment, sticking with identity obs
! to get things rolling.
   
real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus
 
end subroutine model_interpolate




function get_model_time_step()
!------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, 
! returns the associated location. All the heavy lifting is 
! done in the static_init_model subroutine

integer,             intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: var_type

! state_loc is specified in static_init_model
location = state_loc(index_in)

! add var type information here based on value of index_in.
if (present(var_type)) var_type = state_kind(index_in)

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Nothing for 
! MITgcm for now.


end subroutine end_model



subroutine model_get_close_states(o_loc, radius, inum, indices, dist, x)
!------------------------------------------------------------------
!
! Returns the number of state variables that are within a given
! radius (the units for the radius depend upon the location_mod
! module being used by the model) of an observation at location
! o_loc.  The indices in the long state vectors as well as the
! distance between each close state variable and the observation
! are also returned, provided there is sufficient storage 
! available for them in the arrays indices and dist.  
   
type(location_type), intent(in) :: o_loc
real(r8), intent(in)            :: radius  
integer, intent(out)            :: inum, indices(:)
real(r8), intent(out)           :: dist(:)
real(r8), intent(in)            :: x(:)

! the local stuff
type(location_type)             :: loc2
real(r8)                        :: loc_vect(3)
real(r8)                        :: dradcheck, dcheck
integer                         :: io_loc(3), icheck, increment
integer                         :: iaz, irad, ized, itype, icount
integer                         :: iradcheck, iazcheck
integer                         :: i, j, k, ijim, juli

! Brute force default
inum = -1

!! put the location into a vector form instead of a location type
!loc_vect = get_location(o_loc)
!
!!-----------------------------------------------------------
!! define io_loc to contain the integer location of the grid
!! point nearest the observation in o_loc.  use bisection to
!! speed up finding that grid location
!!-----------------------------------------------------------
!
!! bisection to zero in on the closest azimuthal grid point
!icheck = naz/2
!do i = 1, 3
! increment = icheck/2
! if (loc_vect(1) >= azimuthal(icheck)) then
!  icheck = icheck + increment
! else
!  icheck = icheck - increment
! end if
!end do
!i = icheck - increment
!do while (loc_vect(1) >= azimuthal(i))
! i = i + 1
!end do
!io_loc(1) = i - 1
!
!! bisection to zero in on the closest radial grid point
!icheck = nrad/2
!do i = 1, 2
! increment = icheck/2
! if (loc_vect(2) >= radial(icheck)) then
!  icheck = icheck + increment
! else
!  icheck = icheck - increment
! end if
!end do
!i = icheck - increment
!do while (loc_vect(2) >= radial(i))
! i = i + 1
!end do
!io_loc(2) = i - 1
!
!! bisection to zero in on the closest depth grid point
!icheck = nzed/2
!do i = 1, 2
! increment = icheck/2
! if (loc_vect(3) <= height(icheck)) then
!  icheck = icheck + increment
! else
!  icheck = icheck - increment
! end if
!end do
!i = icheck 
!do while (loc_vect(3) <= height(i))
! i = i + 1
!end do
!io_loc(3) = i - 1
!
!!---------------------------------------------------------------
!! io_loc now containes the grid index location nearest (or almost
!! nearest) the observation location.
!!---------------------------------------------------------------
!
!!---------------------------------------------------------------
!! record the point nearest the o_loc for each prognostic variable.
!! note that since no localization in the vertical, am looping over
!! all nzed (QQQ is this the right thing to do?QQQ)
!!---------------------------------------------------------------
!inum = 0
!do k = 1, nzed
!   do i = 1, ntype
!      inum          = inum + 1
!      iaz           = io_loc(1) + (i - 1)*(model_size/ntype)
!      irad          = io_loc(2) + (i - 1)*(model_size/ntype)
!      ized          = k         + (i - 1)*(model_size/ntype)
!!      ized          = io_loc(3) + (i - 1)*(model_size/ntype)
!      indices(inum) = the_grid(iaz, irad, ized)
!      dist(inum)    = get_dist(o_loc, state_loc(indices(inum)))
!   end do
!end do
!
!!-----------------------------------------------------------
!! loop in the positive radial direction
!!-----------------------------------------------------------
!iradcheck = io_loc(2) - 1
!dradcheck = dist(1)
!do while (dradcheck <= radius)
!   iradcheck = iradcheck + 1
!
!!-----------------------------------------------------------
!! move outward in the azimuthal direction to identify
!! points that are close (remember, things are cyclic in the
!! azimuthal direction!)
!!-----------------------------------------------------------
!   dcheck    = dradcheck
!   iazcheck  = io_loc(1)
!   icount    = 0
!   do while (dcheck <= radius)
!      iazcheck = iazcheck + 1
!      if (iazcheck > naz) iazcheck = 1
!      dcheck = get_dist(o_loc, state_loc(the_grid(iazcheck, io_loc(2), io_loc(3))))
!      if (dcheck <= radius) then
!         icount = icount + 1
!         do ijim = 1, ntype
!            do k = 1, nzed
!               inum          = inum + 1 
!               iaz           = iazcheck  + (ijim - 1)*(model_size/ntype)
!               irad          = iradcheck + (ijim - 1)*(model_size/ntype)
!               ized          = k         + (ijim - 1)*(model_size/ntype)
!               indices(inum) = the_grid(iaz, irad, ized)
!               dist(inum)    = dcheck
!            end do
!         end do
!      end if
!   end do
!
!!-----------------------------------------------------------
!! having recorded the number of steps we took in the positive
!! azimuthal direction, go that many steps in the negative
!! azimuthal direction
!!-----------------------------------------------------------
!   do juli = io_loc(1) - icount, io_loc(1)
!      iazcheck = juli
!      if (iazcheck < 1) iazcheck = naz + iazcheck
!      do ijim = 1, ntype
!         do k = 1, nzed
!            inum          = inum + 1
!            iaz           = iazcheck  + (ijim - 1)*(model_size/ntype)
!            irad          = iradcheck + (ijim - 1)*(model_size/ntype)
!            ized          = k         + (ijim - 1)*(model_size/ntype)
!            indices(inum) = the_grid(iaz, irad, ized)
!            dist(inum)    = get_dist(o_loc, state_loc(indices(inum)))
!         end do
!      end do
!   end do
!
!!-----------------------------------------------------------
!! check to see if moving another step in the radial direction is
!! already too far away (or outside the outer boundary). if not, 
!! record the location
!!-----------------------------------------------------------
!   iaz  = io_loc(1)
!   irad = iradcheck + 1
!   ized = io_loc(3)
!   if (irad > nrad) then
!      dradcheck = 100000.00_r8
!   else
!      dradcheck = get_dist(o_loc, state_loc(the_grid(iaz, irad, ized)))
!   end if
!   if (dradcheck < radius) then
!      do k = 1, nzed
!         do i = 1, ntype
!            inum          = inum + 1
!            iaz           = io_loc(1)     + (i - 1)*(model_size/ntype)
!            irad          = iradcheck + 1 + (i - 1)*(model_size/ntype)
!            ized          = k             + (i - 1)*(model_size/ntype)
!            indices(inum) = the_grid(iaz, irad, ized)
!            dist(inum)    = dradcheck
!         end do
!      end do
!   end if
!
!end do
!
!!-----------------------------------------------------------
!! loop in the negative radial direction
!!-----------------------------------------------------------
!iradcheck = io_loc(2)
!dradcheck = dist(1)
!do while (dradcheck <= radius)
!
!   iradcheck = iradcheck - 1
!
!!-----------------------------------------------------------
!! move outward in the azimuthal direction to identify
!! points that are close (remember, things are cyclic in the
!! azimuthal direction!)
!!-----------------------------------------------------------
!   dcheck    = dist(1)
!   iazcheck  = io_loc(1)
!   icount    = 0
!   do while (dcheck <= radius)
!      iazcheck = iazcheck + 1
!      if (iazcheck > naz) iazcheck = 1
!      dcheck = get_dist(o_loc, state_loc(the_grid(iazcheck, io_loc(2), io_loc(3))))
!      if (dcheck <= radius) then
!         icount = icount + 1
!         do ijim = 1, ntype
!            do k = 1, nzed
!               inum          = inum + 1 
!               iaz           = iazcheck  + (ijim - 1)*(model_size/ntype)
!               irad          = iradcheck + (ijim - 1)*(model_size/ntype)
!               ized          = k         + (ijim - 1)*(model_size/ntype)
!               indices(inum) = the_grid(iaz, irad, ized)
!               dist(inum)    = dcheck
!            end do
!         end do
!      end if
!   end do
!
!!-----------------------------------------------------------
!! having recorded the number of steps we took in the positive
!! azimuthal direction, go that many steps in the negative
!! azimuthal direction
!!-----------------------------------------------------------
!   do juli = io_loc(1) - icount, io_loc(1)
!      iazcheck = juli
!      if (iazcheck < 1) iazcheck = naz + iazcheck
!      do ijim = 1, ntype
!         do k = 1, nzed
!            inum          = inum + 1
!            iaz           = iazcheck  + (ijim - 1)*(model_size/ntype)
!            irad          = iradcheck + (ijim - 1)*(model_size/ntype)
!            ized          = k         + (ijim - 1)*(model_size/ntype)
!            indices(inum) = the_grid(iaz, irad, ized)
!            dist(inum)    = get_dist(o_loc, state_loc(indices(inum)))
!         end do
!      end do
!   end do
!
!!-----------------------------------------------------------   
!! check to see if moving another step in the radial direction is
!! already too far away (or inside the inner boundary).  if not, 
!! record the location.
!!-----------------------------------------------------------   
!   iaz  = io_loc(1)
!   irad = iradcheck - 1
!   ized = io_loc(3)
!! QQQ need to put that inner radius thing into the namelist QQQ
!   if (irad < 9) then
!      dradcheck = 100000.00_r8
!   else
!      dradcheck = get_dist(o_loc, state_loc(the_grid(iaz, irad, ized)))
!   end if
!   if (dradcheck < radius) then
!      do k = 1, nzed
!         do i = 1, ntype
!            inum          = inum + 1
!            iaz           = io_loc(1)     + (i - 1)*(model_size/ntype)
!            irad          = iradcheck + 1 + (i - 1)*(model_size/ntype)
!            ized          = k             + (i - 1)*(model_size/ntype)
!            indices(inum) = the_grid(iaz, irad, ized)
!            dist(inum)    = dradcheck
!         end do
!      end do
!   end if
!
!end do
!
!write(6,*) ' '
!write(6,*) 'Close state num, vector index, distance'
!do i = 1, inum
! write(6,*) i, indices(i), dist(i)
!Xend do

end subroutine model_get_close_states


subroutine vector_to_prog_var(x, var)
!=======================================================================
! subroutine vector_to_prog_var(x, var)
!
! This subroutine takes in the control vector, x, and splits it
! into 3d model prognostic variable fields
   
real(r8), intent(in)          :: x(:)
type(model_type), intent(out) :: var

integer :: i, j, k, it, icount
character(len=129) :: errstring

! set counter for state
icount = 1

! loop over the number of variable types and put them into
! the var variables
do it = 1, ntype
   do k = 1, nzed
      do j = 1, nrad
         do i = 1, naz
            var%vars_3d(i,j,k,it)=x(icount)
            icount = icount + 1
         enddo
      end do
   end do
end do

end subroutine vector_to_prog_var



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH Jan 24 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_04 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!--------------------------------------------------------------------
! General netCDF variables
!--------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: azDimID, radDimID, zedDimID 
integer :: azVarID, radVarID, zedVarID
integer :: UVarID, VVarID, WVarID, TVarID, PVarID

!--------------------------------------------------------------------
! netCDF variables for Location
!--------------------------------------------------------------------

integer :: LocationVarID
integer :: StateVarDimID, StateVarVarID
integer :: StateVarID, MemberDimID, TimeDimID

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

character(len=129)    :: errstring
character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1,str2

integer             :: i, itype
type(location_type) :: lctn
ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_Redef(ncFileID))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies
!-------------------------------------------------------------------------------

call check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID))
call check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID))

if ( TimeDimID /= unlimitedDimId ) then
  write(errstring,*)'Time dimension ID ',TimeDimID,'must match Unlimited Dimension ID ',unlimitedDimId
  call error_handler(E_ERR,'nc_write_model_atts', errstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-------------------------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                        len=model_size, dimid = StateVarDimID))

!--------------------------------------------------------------------
! Write Global Attributes
!--------------------------------------------------------------------
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source", source ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision", revision ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate", revdate ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "MITgcm_annulus"))

!-------------------------------------------------------------------------------
! Define the new dimensions IDs
!-------------------------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="az",     len = naz,  dimid = azDimID))
call check(nf90_def_dim(ncid=ncFileID, name="rad",    len = nrad, dimid = radDimID))
call check(nf90_def_dim(ncid=ncFileID, name="depth",  len = nzed, dimid = zedDimID))

!-------------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and their attributes 
!-------------------------------------------------------------------------------

! Grid azimuthal angle
call check(nf90_def_var(ncFileID, name="az", xtype=nf90_double, dimids=azDimID, varid=azVarID) )
call check(nf90_put_att(ncFileID, azVarID, "long_name", "azimuthal angle"))
call check(nf90_put_att(ncFileID, azVarID, "units", "degrees_east"))
call check(nf90_put_att(ncFileID, azVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)))
                  
! Grid radius
call check(nf90_def_var(ncFileID, name="rad", xtype=nf90_double, dimids=radDimID, varid=radVarID) )
call check(nf90_put_att(ncFileID, radVarID, "long_name", "radius"))   
call check(nf90_put_att(ncFileID, radVarID, "units", "meters"))
call check(nf90_put_att(ncFileID, radVarID, "valid_range", ">0"))

! Grid depth (height)
call check(nf90_def_var(ncFileID, name="zed", xtype=nf90_double, dimids=zedDimID, varid=zedVarID) )
call check(nf90_put_att(ncFileID, zedVarID,  "long_name", "depth"))
call check(nf90_put_att(ncFileID, zedVarID,  "units", "meters"))
call check(nf90_put_att(ncFileID, zedVarID,  "negative", "down"))
call check(nf90_put_att(ncFileID, zedVarID,  "valid_range", "<0"))

!--------------------------------------------------------------------
! Output either the "state vector" variables -OR- the "prognostic" 
! variables.  The "state vector" choice simply dumps the control 
! vector while the "prognostic" choice creates separate files for
! U, V, W, T, and P.
!--------------------------------------------------------------------

if ( output_state_vector ) then

   !-----------------------------------------------------------------
   ! Create the (empty) variables and attributes for the control 
   ! vector output.
   !-----------------------------------------------------------------

   ! Define the state vector coordinate variable
   call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID))
   call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"))
   call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical") )
   call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)))

   ! Define the actual state vector
   call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), varid=StateVarID))
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"))
   call check(nf90_put_att(ncFileID, StateVarId, "vector_to_prog_var", "MITgcm annulus"))
  
   ! Leave define mode so we can fill
   call check(nf90_enddef(ncfileID))

   ! Fill the state variable coordinate variable
   call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ))

else

   !-----------------------------------------------------------------
   ! Create the (empty) variables and attributes for the prognostic
   ! variable output (U, V, W, T, and P).
   !-----------------------------------------------------------------

   ! U ... itype == 1 of ntype
   itype = 1
   call check(nf90_def_var(ncid=ncFileID, name=trim(flds(itype)), xtype=nf90_real, &
          dimids = (/ azDimID, radDimID, zedDimID, MemberDimID, unlimitedDimID /), & 
          varid  = UVarID))
   call check(nf90_put_att(ncFileID, UVarID, "long_name", "Azimuthal Wind"))
   call check(nf90_put_att(ncFileID, UVarID, "units", "m/s"))
   
   ! V ... itype == 2 of ntype
   itype = 2
   call check(nf90_def_var(ncid=ncFileID, name=trim(flds(itype)), xtype=nf90_real, &
          dimids = (/ azDimID, radDimID, zedDimID, MemberDimID, unlimitedDimID /), &
          varid  = VVarID))
   call check(nf90_put_att(ncFileID, VVarID, "long_name", "Radial Wind"))
   call check(nf90_put_att(ncFileID, VVarID, "units", "m/s"))

   ! W ... itype == 3 of ntype
   itype = 3
   call check(nf90_def_var(ncid=ncFileID, name=trim(flds(itype)), xtype=nf90_real, &
          dimids = (/ azDimID, radDimID, zedDimID, MemberDimID, unlimitedDimID /), &
          varid  = VVarID))
   call check(nf90_put_att(ncFileID, VVarID, "long_name", "Vertical Wind"))
   call check(nf90_put_att(ncFileID, VVarID, "units", "m/s"))

   ! T ... itype == 4 of ntype
   itype = 4
   call check(nf90_def_var(ncid=ncFileID, name=trim(flds(itype)), xtype=nf90_real, &
          dimids = (/ azDimID, radDimID, zedDimID, MemberDimID, unlimitedDimID /), &
          varid  = TVarID))
   call check(nf90_put_att(ncFileID, TVarID, "long_name", "Temperature"))
   call check(nf90_put_att(ncFileID, TVarID, "units", "C"))

   ! P ... itype == 5 of ntype
   itype = 5
   call check(nf90_def_var(ncid=ncFileID, name=trim(flds(itype)), xtype=nf90_real, &
          dimids = (/ azDimID, radDimID, zedDimID, MemberDimID, unlimitedDimID /), &
          varid  = VVarID))
   call check(nf90_put_att(ncFileID, VVarID, "long_name", "Pressure"))
   call check(nf90_put_att(ncFileID, VVarID, "units", "Pa/1000"))
   
   ! Leave define mode so we can fill
   call check(nf90_enddef(ncfileID))
   
endif

! Fill the coordinate variables
call check(nf90_put_var(ncFileID, azVarID,  azimuthal ))
call check(nf90_put_var(ncFileID, radVarID, radial ))
call check(nf90_put_var(ncFileID, zedVarID, height ))

!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------
call check(nf90_sync(ncFileID))

write (*,*)'Model attributes written, netCDF file synched ...'

contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
      integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_atts',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check
end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!
! Writes the model-specific variables to a netCDF file
! TJH 25 June 2003
!
! The MITgcm annulus uses the C-grid, and the current configuration
! is passing five prognostic variables to DART.
!
! The routine "vector_to_prog_var" pulls the prognostic variables from
! the control vector and places them into 3d variables that are 
! output to netcdf.
   
use typeSizes
use netcdf
          
integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function
integer                            :: ifld, itype, ncfldid
!-------------------------------------------------------------------------------
! General netCDF variables
!-------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID

!-------------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------------

ierr = 0                      ! assume normal termination

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

! check to see if dumping control vector or prognostic variables
if ( output_state_vector ) then
           
   ! check that file for state variable is open, then fill it up
   call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)))

else

   ! convert control vector to prognostic variable form
   call vector_to_prog_var(statevec, var)


   ! loop over prognostic variables
   do ifld = 1, ntype

      call check(NF90_inq_varid(ncFileID, trim(flds(ifld)), ncfldid))
      call check(nf90_put_var( ncFileID, ncfldid, var%vars_3d(:,:,:,ifld), &
           start=(/ 1,1,1,copyindex,timeindex /), count=(/naz,nrad,nzed,1,1/) ))

   end do 

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
   
write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
write (*,*)'netCDF file is synched ...'

contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
   integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_vars',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check
end function nc_write_model_vars



subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
! subroutine pert_model_state(state, pert_state, interf_provided)
!
! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

integer :: i

interf_provided = .true.

do i = 1, model_size
   pert_state(i) = random_gaussian(random_seq, state(i), 0.0002_r8)
end do


end subroutine pert_model_state




subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

end subroutine ens_mean_for_model



!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

module model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

use types_mod
use location_mod, only : location_type, get_dist, set_location, get_location, &
                         nc_write_location
use time_manager_mod

private

public static_init_model, init_conditions, get_model_size, adv_1step,  &
   init_time, model_interpolate, get_model_time_step, get_state_meta_data, end_model, &
   init_model


 integer,  parameter :: model_size =   40
 real(r8), parameter ::    forcing = 8.00_r8
 real(r8), parameter ::    delta_t = 0.05_r8   ! Original timestep 
!real(r8), parameter ::    delta_t = 0.005_r8  ! timestep for assim experiments

! Define the location of the state variables in module storage
type(location_type) :: state_loc(model_size)
type(time_type)     :: time_step

contains

!======================================================================



subroutine static_init_model()
!----------------------------------------------------------------------
! subroutine static_init_model()
!
! Initializes class data for this model. For now, simply outputs the
! identity info, sets the location of the state variables, and initializes
! the time type for the time stepping (is this general enough for time???)

implicit none
real(r8) :: x_loc
integer :: i
character(len=128) :: source,revision,revdate

! let CVS fill strings ... DO NOT EDIT ...

source   = "$Source$"
revision = "$Revision$"
revdate  = "$Date$"

! Ultimately,  change output to diagnostic output block ...

write(*,*)'model attributes:'
write(*,*)'   ',trim(adjustl(source))
write(*,*)'   ',trim(adjustl(revision))
write(*,*)'   ',trim(adjustl(revdate))

! Define the locations of the model state variables
do i = 1, model_size
   x_loc = (i - 1.0) / model_size
   state_loc(i) =  set_location(x_loc)
end do


! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.
time_step = set_time(3600, 0)


end subroutine static_init_model



subroutine init_model_instance()
!---------------------------------------------------------------------
! subroutine init_model_instance
!
! Initializes instance dependent state for model. Null for L96.

end subroutine init_model_instance




subroutine comp_dt(x, dt)
!----------------------------------------------------------------------
! subroutine comp_dt(x, dt)
! 
! Computes the time tendency of the lorenz 1996 model given current state

implicit none

real(r8), intent( in) :: x(:)
real(r8), intent(out) :: dt(:)

integer :: j, jp1, jm1, jm2

do j = 1, model_size
   jp1 = j + 1
   if(jp1 > model_size) jp1 = 1
   jm2 = j - 2
   if(jm2 < 1) jm2 = model_size + jm2
   jm1 = j - 1
   if(jm1 < 1) jm1 = model_size
   
   dt(j) = (x(jp1) - x(jm2)) * x(jm1) - x(j) + forcing
end do

end subroutine comp_dt



subroutine adv_1step(x, time)
!----------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does single time step advance for lorenz 96 model
! using four-step rk time step
! The Time argument is needed for compatibility with more complex models
! that need to know the time to compute their time tendency and is not
! used in L96. Is there a better way to do this in F90 than to just hang
! this argument out everywhere?

implicit none

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

real(r8), dimension(size(x)) :: x1, x2, x3, x4, dx, inter
integer :: i

!  Compute the first intermediate step

call comp_dt(x, dx)
x1    = delta_t * dx
inter = x + x1 / 2.0_r8

!  Compute the second intermediate step

call comp_dt(inter, dx)
x2    = delta_t * dx
inter = x + x2 / 2.0_r8

!  Compute the third intermediate step

call comp_dt(inter, dx)
x3    = delta_t * dx
inter = x + x3

!  Compute fourth intermediate step

call comp_dt(inter, dx)
x4 = delta_t * dx

!  Compute new value for x

x = x + x1/6.0_r8 + x2/3.0_r8 + x3/3.0_r8 + x4/6.0_r8

end subroutine adv_1step




subroutine init_conditions(x)
!----------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Initial conditions for lorenz 96
! It is assumed that this is called before any other routines in this
! module. Should probably make that more formal and perhaps enforce for
! more comprehensive models.

implicit none

real(r8), intent(out) :: x(:)

x    = forcing
x(1) = 1.001_r8 * forcing

end subroutine init_conditions



subroutine init_time(time)
!----------------------------------------------------------------------
!
! Gets the initial time for a state from the model. Where should this info
! come from in the most general case?

implicit none

type(time_type), intent(out) :: time

! For now, just set to 0
time = set_time(0, 0)

end subroutine init_time



function model_interpolate(x, location, type)
!---------------------------------------------------------------------
!
! Interpolates from state vector x to the location. It's not particularly
! happy dumping all of this straight into the model. Eventually some
! concept of a grid underlying models but above locations is going to
! be more general. May want to wait on external infrastructure projects
! for this?

! Argument type is not used here because there is only one type of variable.
! Type is needed to allow swap consistency with more complex models.

implicit none

real(r8) :: model_interpolate
real(r8), intent(in) :: x(:)
type(location_type), intent(in) :: location
integer, intent(in) :: type

integer :: lower_index, upper_index
real(r8) :: loc, fraction

! Convert location to real
loc = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
loc = model_size * loc

lower_index = int(loc) + 1
upper_index = lower_index + 1
if(lower_index > model_size) lower_index = lower_index - model_size
if(upper_index > model_size) upper_index = upper_index - model_size

fraction = loc - int(loc)
model_interpolate = (1.0_r8 - fraction) * x(lower_index) + fraction * x(upper_index)

end function model_interpolate




function get_model_size()
!---------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

integer :: get_model_size

get_model_size = model_size

end function get_model_size



function get_model_time_step()
!------------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



subroutine get_state_meta_data(index, location)
!---------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?

implicit none

integer, intent(in) :: index
type(location_type), intent(out) :: location

location = state_loc(index)

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Nothing for L96 for now.


end subroutine end_model



subroutine nc_write_locations(ncFileID)
!----------------------------------------------------------------------------
!
! Defines a oned location variable and the attributes in the netCDF file.
! The dimension of the variable is set by passing in the dimids ... 
! in other models, dimids will be an array ...
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
!

use typeSizes
use netcdf
implicit none

integer, intent(inout)       :: ncFileID

integer :: paramDimID, LocationVarID, Nlocations
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: i,ierr

! make sure ncFileID refers to an open netCDF file, and then put into define mode.

call check(NF90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(NF90_Redef(ncFileID))

! find the DimensionID(s), create variable and atts
! then leave define mode

call check(NF90_inq_dimid(ncid=ncFileID, name="StateVariable", dimid = ParamDimID))
call check(NF90_def_var(ncFileID, name="loc1d", xtype=nf90_double, &
                                      dimids=ParamDimID, varid=LocationVarID) )
call check(NF90_put_att(ncFileID, LocationVarID, "long_name", "one-dimensional location"))
call check(NF90_put_att(ncFileID, LocationVarID, "units", "nondimensional"))
call check(NF90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)))
call check(NF90_enddef(ncfileID))   ! Leave define mode

! Write all the locations
! Question:  do we want to use get_model_size ...

call check(NF90_Inquire_Dimension(ncid=ncFileID, dimid = ParamDimID, len=Nlocations))

if ( Nlocations /= model_size ) then  
   write(*, *) 'Error: nc_write_locations : model_size is ',model_size
   write(*, *) '       but netCDF thinks we have ',Nlocations,' state variables.'
   stop
endif

do i = 1, Nlocations
!  call get_state_meta_data(i, state_loc)
   call nc_write_location(ncFileID, LocationVarID, state_loc(i), start=i)
enddo

call check(nf90_sync(ncFileID))      ! LEAVES FILE OPEN -- tjh

contains

  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(NF90_strerror(status))
      print *, 'in model_mod:nc_write_locations'
      stop
    end if
  end subroutine check

end subroutine nc_write_locations


!
!===================================================================
! End of model_mod
!===================================================================
!
end module model_mod

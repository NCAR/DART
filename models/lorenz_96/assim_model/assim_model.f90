module assim_model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

! NEED TO ADD ON ONLY CLAUSES
use types_mod
use location_mod
use time_manager_mod
use utilities_mod

private

public initialize_assim_model, init_diag_output, get_model_size, get_closest_state_time_to, &
   get_initial_condition, get_state_meta_data, get_close_states, get_num_close_states, &
   get_model_time, get_model_state_vector, copy_assim_model, advance_state, interpolate, &
   set_model_time, set_model_state_vector, write_state_restart, read_state_restart, &
   output_diagnostics, end_assim_model, assim_model_type

integer,  parameter :: model_size =   40
real(r8), parameter ::    forcing = 8.00_r8
real(r8), parameter ::    delta_t = 0.05_r8

! Eventually need to be very careful to implement this to avoid state vector copies which
! will be excruciatingly costly (storage at least) in big models.
type assim_model_type
   real(r8) :: state_vector(model_size)
   type(time_type) :: time
end type assim_model_type

logical :: output_init = .FALSE.

! Define output indices for diagnostics

integer :: diag_output_index(9)

! Define the location of the state variables in module storage

type(location_type) :: state_loc(model_size)

type(time_type) :: time_step

contains

!======================================================================



subroutine initialize_assim_model()
!----------------------------------------------------------------------
! subroutine initialize_assim_model()
!
! Initializes class data for the Lorenz-96 model. So far, this simply 
! is initializing the position of the state variables as location types.

implicit none
real(r8) :: x_loc
integer :: i

! Define the locations of the state variables;
do i = 1, model_size
   x_loc = (i - 1.0) / model_size
   state_loc(i) =  set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.

time_step = set_time(3600, 0)

end subroutine initialize_assim_model




function init_diag_output(file_name, global_meta_data, &
   copies_of_field_per_time, meta_data_per_copy)
!---------------------------------------------------------------------
!
! Initializes a diagnostic output file. Should be NetCDF shortly but
! for now just opens file and dumps stuff. A file_id is returned which
! is simply implemented as an integer unit number for now. 

implicit none

integer :: init_diag_output
character, intent(in) :: file_name(:)
character, intent(in) :: global_meta_data(:)
integer, intent(in) :: copies_of_field_per_time
character, intent(in) :: meta_data_per_copy(:, :)

integer :: i

init_diag_output = open_file(file_name)
write(init_diag_output) global_meta_data
do i = 1, copies_of_field_per_time
   write(init_diag_output, *) i, meta_data_per_copy(:, i)
end do

end function init_diag_output




function get_model_size()
!---------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

implicit none

integer :: get_model_size

get_model_size = model_size

end function get_model_size



function get_closest_state_time_to(assim_model, time)
!----------------------------------------------------------------------
!
! Returns the time closest to the given time that the model can reach
! with its state. For L96, the time step is fixed at ???

implicit none

type(time_type) :: get_closest_state_time_to
type(assim_model_type), intent(in) :: assim_model
type(time_type), intent(in) :: time

type(time_type) :: model_time, delta_time

! CAREFUL WITH FLOATING POINT DIVISION AND COMPARISONS

model_time = assim_model%time
if(model_time > time) then
   write(*, *) 'Error in get_closest_state_to_time: model time > time'
   stop
endif

delta_time = time - model_time

get_closest_state_time_to = (delta_time / time_step) * time_step + model_time

end function get_closest_state_time_to



function get_initial_condition()
!----------------------------------------------------------------------
! function get_initial_condition()
!
! Initial conditions for lorenz 96. This returns an initial assim_model_type
! which includes both a state vector and a time. Design of exactly where this 
! stuff should come from is still evolving (6 April, 2002) but for now can 
! start at time offset 0 with the initial state pulled from the original paper.
! Need to carefully coordinate this with the times for observations.

implicit none

type(assim_model_type) :: get_initial_condition

integer  :: i
real(r8) :: x_loc

! Set the state vector
do i = 1, model_size
   get_initial_condition%state_vector = forcing
end do
get_initial_condition%state_vector(1) = 1.001_r8 * forcing

! Set the time
! POTENTIAL PROBLEM: TIME MANAGER DOES NOT USE THE R8 FORMALISM: HOW TO HANDLE THIS?
get_initial_condition%time = set_time(0, 0)

end function get_initial_condition



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




subroutine get_close_states(location, radius, number, indices)
!---------------------------------------------------------------------
! subroutine get_close_states(location, radius, number, indices)
!
! Returns a list of indices for model state vector points that are
! within distance radius of the location. Might want to add an option
! to return the distances, too.

implicit none

type(location_type), intent(in) :: location
real(r8), intent(in) :: radius
integer, intent(out) :: number, indices(:)

integer :: index, i

! For large models this will have to be VERY efficient; here can just search
index = 0
do i = 1, model_size
   if(get_dist(location, state_loc(i)) < radius) then
      index = index + 1 
      if(index <= size(indices)) indices(index) = i
   end if
end do

if(index > size(indices)) then
   number = -1 * index
else 
   number = index
end if
      
end subroutine get_close_states




function get_num_close_states(location, radius)
!-----------------------------------------------------------------------
!
! Returns number of state vector points located within distance radius
! of the location.

integer :: get_num_close_states
type(location_type), intent(in) :: location
real(r8), intent(in) :: radius

integer :: i

! For large models this will have to be VERY efficient; here can just search
number = 0
do i = 1, model_size
   if(get_dist(location, state_loc(i)) < radius) number = number + 1
end do

end function get_num_close_states



function get_model_time(assim_model)
!-----------------------------------------------------------------------
!
! Returns the time component of a assim_model extended state.

type(time_type) :: get_model_time
type(assim_model_type), intent(in) :: assim_model

get_model_time = assim_model%time

end function get_model_time



function get_model_state_vector(assim_model)
!-----------------------------------------------------------------------
!
! Returns the state vector  component of a assim_model extended state.

real(r8) :: get_model_state_vector(model_size)
type(assim_model_type), intent(in) :: assim_model

get_model_state_vector = assim_model%state_vector

end function get_model_state_vector




function copy_assim_model(assim_model)
!-------------------------------------------------------------------------
!
! Does a copy of assim_model, should be overloaded to =. Still need to be
! very careful about trying to limit copies of the potentially huge state
! vectors for big models. 

implicit none

type(assim_model_type) :: copy_assim_model
type(assim_model_type), intent(in) :: assim_model
! Null stub for now

end function copy_assim_model





function advance_state(assim_model, target_time)
!-----------------------------------------------------------------------
!
! Advances the model extended state until time is equal (within roundoff?)
! of the target_time. For L96 this is relatively straightforward with 
! fixed time steps, etc.

implicit none

type(assim_model_type) :: advance_state
type(assim_model_type), intent(in) :: assim_model
type(time_type), intent(in) :: target_time

type(time_type) :: model_time

! NEED TO BE CAREFUL ABOUT FLOATING POINT TESTS: Being sloppy here

model_time = get_model_time(assim_model)
! Check for time error; use error handler when available
if(model_time > target_time) then
   write(*, *) 'Error in advance_state, target_time before model_time'
   stop
endif

advance_state = assim_model

do while(model_time < target_time)
   advance_state = adv_1step(assim_model)
   model_time = get_model_time(advance_state)
end do

end function advance_state



function interpolate(x, location)
!---------------------------------------------------------------------
!
! Interpolates from state vector x to the location. Might want to overload
! this to allow the first argument to be an assim_model_type, too. Need to
! decide if state_vectors or assim_model_types are coming down the obs
! side of the calling tree. Might help to have the time along?

implicit none

real(r8) :: interpolate
real(r8), intent(in) :: x(:)
type(location_type), intent(in) :: location

integer :: lower_index, upper_index
real(r8) :: loc, fraction

! Convert location to real
loc = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
loc = model_size * loc

lower_index = int(loc)
upper_index = lower_index + 1
if(upper_index > model_size) upper_index = 1
if(lower_index == 0) lower_index = model_size

fraction = loc - int(loc)
interpolate = (1.0_r8 - fraction) * x(lower_index) + fraction * x(upper_index)

end function interpolate



subroutine set_model_time(assim_model, time)
!-----------------------------------------------------------------------
!
! Sets the time in an assim_model type

implicit none

type(assim_model_type), intent(inout) :: assim_model
type(time_type), intent(in) :: time

assim_model%time = time

end subroutine set_model_time



subroutine set_model_state_vector(assim_model, state)
!-----------------------------------------------------------------------
!
! Sets the state vector part of an assim_model_type

implicit none

type(assim_model_type), intent(inout) :: assim_model
real(r8), intent(in) :: state(:)

! Check the size for now
if(size(state) /= model_size) then
   write(*, *) 'Error: Input state vector is wrong size in set_model_state_vector'
   stop
endif

assim_model%state_vector = state

end subroutine set_model_state_vector



subroutine write_state_restart(assim_model, file)
!----------------------------------------------------------------------
!
! Write a restart file given a model extended state and a unit number 
! opened to the restart file. (Need to reconsider what is passed to 
! identify file or if file can even be opened within this routine).

implicit none

type (assim_model_type), intent(in) :: assim_model
integer, intent(in) :: file

integer :: days, seconds

! This needs to be done more carefully, consider this an extended stub
call get_time(assim_model%time, days, seconds)
write(file, *) seconds, days

! Write the state vector
write(file, *) assim_model%state_vector

end subroutine write_state_restart



function read_state_restart(file)
!----------------------------------------------------------------------
!
! Read a restart file given a unit number (see write_state_restart)

implicit none

type(assim_model_type) :: read_state_restart
integer, intent(in) :: file

integer :: seconds, days

! WARNING : TIME TYPE IS CURRENTLY DOING REGULAR REALS, NOT R8s
read(file, *) seconds, days
read_state_restart%time = set_time(seconds, days)

! Read the state vector
read(file, *) read_state_restart%state_vector

end function read_state_restart



subroutine output_diagnostics(file_id, state_vector, time, copy_index)
!-------------------------------------------------------------------
!
! Outputs a copy of the state vector to the file (currently just an
! integer unit number), the time, and an optional index saying which
! copy of the metadata this state is associated with.

integer, intent(in) :: file_id
real(r8), intent(in) :: state_vector(:)
type(time_type), intent(in) :: time
integer, optional, intent(in) :: copy_index

! Stub for now


end subroutine output_diagnostics




subroutine end_assim_model()
!--------------------------------------------------------------------
!
! Closes down assim_model; nothing to do for L96

end subroutine end_assim_model




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



function adv_1step(state)
!----------------------------------------------------------------------
! function adv_1step(state)
!
! Does single time step advance for lorenz 96 model
! using four-step rk time step

implicit none

type(assim_model_type) :: adv_1step
type(assim_model_type), intent(in) :: state

real(r8) :: x(:)
real(r8), dimension(model_size) :: x1, x2, x3, x4, dx, inter
type(time_type) :: time
integer :: i

! Update the time: WARNING: THERE'S SOME NATURAL TIME INTERVAL FOR L96: NEED 
! TO FIND IT (CHECK WITH SHREE KHARE). For now do 1 hour (3600 secs).
time = time + set_time(3600, 0)

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

! Load x and the updated time back into output
call set_model_time(adv_1step, time)
call set_model_state_vector(adv_1step, x)

end function adv_1step



!
!===================================================================
! End of assim_model_mod
!===================================================================
!
end module assim_model_mod

module assim_model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

! NEED TO ADD ON ONLY CLAUSES
use location_mod, only : location_type, write_location, set_location, get_dist, &
   read_location, get_location
! I've had a problem with putting in the only for time_manager on the pgf90 compiler (JLA).
use time_manager_mod
use utilities_mod, only : get_unit
use types_mod

private

public static_init_assim_model, init_diag_output, get_model_size, get_closest_state_time_to, &
   get_initial_condition, get_state_meta_data, get_close_states, get_num_close_states, &
   get_model_time, get_model_state_vector, copy_assim_model, advance_state, interpolate, &
   set_model_time, set_model_state_vector, write_state_restart, read_state_restart, &
   output_diagnostics, end_assim_model, assim_model_type, init_diag_input, input_diagnostics, &
   get_diag_input_copy_meta_data, init_assim_model

integer,  parameter :: model_size =   40
real(r8), parameter ::    forcing = 8.00_r8
real(r8), parameter ::    delta_t = 0.05_r8

! Eventually need to be very careful to implement this to avoid state vector copies which
! will be excruciatingly costly (storage at least) in big models.
type assim_model_type
   private
!!!   real(r8) :: state_vector(model_size)
! Need to change to pointer to do more efficient pointer work in assim
   real(r8), pointer :: state_vector(:)
   type(time_type) :: time
end type assim_model_type

! Define the location of the state variables in module storage
type(location_type) :: state_loc(model_size)

type(time_type) :: time_step

contains

!======================================================================



subroutine static_init_assim_model()
!----------------------------------------------------------------------
! subroutine static_init_assim_model()
!
! Initializes class data for the Lorenz-96 model. So far, this simply 
! is initializing the position of the state variables as location types.

implicit none


character(len=128) :: source,revision,revdate
real(r8) :: x_loc
integer :: i

! Change output to diagnostic output block ... 

source   = "$Source$"
revision = "$Revision$"
revdate  = "$Date$"

! Change output to diagnostic output block ... 

write(*,*)'assim_model attributes:'
write(*,*)'   ',source
write(*,*)'   ',revision
write(*,*)'   ',revdate

! Define the locations of the state variables;

do i = 1, model_size
   x_loc = (i - 1.0) / model_size
   state_loc(i) =  set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.


time_step = set_time(3600, 0)

end subroutine static_init_assim_model



subroutine init_assim_model(state)
!----------------------------------------------------------------------
!
! Allocates storage for an instance of an assim_model_type. With this
! implementation, need to be VERY careful about assigment and maintaining
! permanent storage locations. Need to revisit the best way to do 
! assim_model_copy below.

implicit none

type(assim_model_type), intent(inout) :: state

allocate(state%state_vector(model_size))

end subroutine init_assim_model




function init_diag_output(file_name, global_meta_data, &
   copies_of_field_per_time, meta_data_per_copy)
!---------------------------------------------------------------------
!
! Initializes a diagnostic output file. Should be NetCDF shortly but
! for now just opens file and dumps stuff. A file_id is returned which
! is simply implemented as an integer unit number for now. 

! NOTE: Almost all of this output plus many other things can go in a 
! general model independent tool kit to avoid the burden of creating
! assim_model stuff.

implicit none

integer :: init_diag_output
character(len = *), intent(in) :: file_name, global_meta_data
integer, intent(in) :: copies_of_field_per_time
character(len = *), intent(in) :: meta_data_per_copy(copies_of_field_per_time)

integer :: i

init_diag_output = get_unit()
open(unit = init_diag_output, file = file_name)
!!!init_diag_output = open_file(file_name)
write(init_diag_output, *) global_meta_data

! Write the model size
write(init_diag_output, *) model_size

! Write number of copies of field per time plus the meta data per copy
write(init_diag_output, *) copies_of_field_per_time
do i = 1, copies_of_field_per_time
   write(init_diag_output, *) i, meta_data_per_copy(i)
end do

! Will need other metadata, too; Could be as simple as writing locations
write(init_diag_output, *) 'locat'
do i = 1, model_size
   call write_location(init_diag_output, state_loc(i))
end do

end function init_diag_output



function init_diag_input(file_name, global_meta_data, model_size, copies_of_field_per_time)
!--------------------------------------------------------------------------
!
! Initializes a model state diagnostic file for input. A file id is
! returned which for now is just an integer unit number.

implicit none

integer :: init_diag_input
character(len = *), intent(in) :: file_name
character(len = *), intent(out) ::  global_meta_data
integer, intent(out) :: model_size, copies_of_field_per_time

integer :: i

init_diag_input = get_unit()
open(unit = init_diag_input, file = file_name)
read(init_diag_input, *) global_meta_data

! Read the model size
read(init_diag_input, *) model_size

! Read the number of copies of field per time
read(init_diag_input, *) copies_of_field_per_time

end function init_diag_input



subroutine get_diag_input_copy_meta_data(file_id, model_size_out, num_copies, &
   location, meta_data_per_copy)
!-------------------------------------------------------------------------
!
! Returns the meta data associated with each copy of data in
! a diagnostic input file. Should be called immediately after 
! function init_diag_input.

implicit none

integer, intent(in) :: file_id, model_size_out, num_copies
type(location_type), intent(out) :: location(model_size_out)
character(len = *) :: meta_data_per_copy(num_copies)

character(len=129) :: header
integer :: i, j

! Should have space checks, etc here
! Read the meta data associated with each copy
do i = 1, num_copies
   read(file_id, *) j, meta_data_per_copy(i)
end do

! Will need other metadata, too; Could be as simple as writing locations
read(file_id, *) header
if(header /= 'locat') then
   write(*, *) 'Error: get_diag_input_copy_meta_data expected to read "locat"'
   stop
endif

! Read in the locations
do i = 1, model_size_out
   location(i) =  read_location(file_id)
end do

end subroutine get_diag_input_copy_meta_data





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

type(assim_model_type), intent(in) :: assim_model
type(time_type), intent(in) :: time
type(time_type) :: get_closest_state_time_to

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



subroutine get_initial_condition(x)
!----------------------------------------------------------------------
! function get_initial_condition()
!
! Initial conditions for lorenz 96. This returns an initial assim_model_type
! which includes both a state vector and a time. Design of exactly where this 
! stuff should come from is still evolving (6 April, 2002) but for now can 
! start at time offset 0 with the initial state pulled from the original paper.
! Need to carefully coordinate this with the times for observations.

implicit none

type(assim_model_type), intent(inout) :: x

integer  :: i
real(r8) :: x_loc

! Set the state vector
do i = 1, model_size
   x%state_vector(i) = forcing
end do
x%state_vector(1) = 1.001_r8 * forcing

! Set the time
x%time = set_time(0, 0)

end subroutine get_initial_condition



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




subroutine get_close_states(location, radius, number, indices, dist)
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
real(r8), intent(out) :: dist(:)

integer :: index, i
real(r8) :: this_dist

! For large models this will have to be VERY efficient; here can just search
index = 0
do i = 1, model_size
   this_dist = get_dist(location, state_loc(i))
   if(this_dist < radius) then
      index = index + 1 
      if(index <= size(indices)) indices(index) = i
      if(index <= size(dist)) dist(index) = this_dist
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
get_num_close_states = 0
do i = 1, model_size

! INTERESTING NOTE: Because of floating point round-off in comps
! this can give a 'variable' number of num close for certain obs, should fix
   if(get_dist(location, state_loc(i)) < radius) get_num_close_states= get_num_close_states + 1
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




subroutine copy_assim_model(model_out, model_in)
!-------------------------------------------------------------------------
!
! Does a copy of assim_model, should be overloaded to =? Still need to be
! very careful about trying to limit copies of the potentially huge state
! vectors for big models.  Interaction with pointer storage?

implicit none

type(assim_model_type), intent(out) :: model_out
type(assim_model_type), intent(in) :: model_in

integer :: i

! Need to make sure to copy the actual storage and not just the pointer (verify)
model_out%time = model_in%time

do i = 1, model_size
   model_out%state_vector(i) = model_in%state_vector(i)
end do

end subroutine copy_assim_model





subroutine advance_state(assim_model, target_time)
!-----------------------------------------------------------------------
!
! Advances the model extended state until time is equal (within roundoff?)
! of the target_time. For L96 this is relatively straightforward with 
! fixed time steps, etc.

implicit none

type(assim_model_type), intent(inout) :: assim_model
type(time_type), intent(in) :: target_time

type(time_type) :: model_time

! NEED TO BE CAREFUL ABOUT FLOATING POINT TESTS: Being sloppy here

model_time = get_model_time(assim_model)
! Check for time error; use error handler when available
if(model_time > target_time) then
   write(*, *) 'Error in advance_state, target_time before model_time'
   stop
endif

do while(model_time < target_time)
   call adv_1step(assim_model)
   model_time = get_model_time(assim_model)
end do

end subroutine advance_state



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



subroutine output_diagnostics(file_id, state, copy_index)
!-------------------------------------------------------------------
!
! Outputs a copy of the state vector to the file (currently just an
! integer unit number), the time, and an optional index saying which
! copy of the metadata this state is associated with.
! Need to make a much better coordinated facility for doing this,
! providing buffering, insuring that ordering is appropriate for
! time (and copy), etc. For now, this just writes what it receives
! with a header stating time and copy_index.

integer, intent(in) :: file_id
type(assim_model_type), intent(in) :: state
integer, optional, intent(in) :: copy_index

! Write the time and copy_index
call write_time(file_id, state%time)
write(file_id, *) 'fcopy '
write(file_id, *) copy_index

! Write the data, unformatted for now
write(file_id, *) state%state_vector

end subroutine output_diagnostics



subroutine input_diagnostics(file_id, state, copy_index)
!------------------------------------------------------------------
!
! Reads in diagnostic state output from file_id for copy_index
! copy. Need to make this all more rigorously enforced.

implicit none

integer, intent(in) :: file_id
! MAYBE SHOULDN'T use assim model type here, but just state and time ?
type(assim_model_type), intent(out) :: state
integer, intent(out) :: copy_index

character*5 :: header

! Read in the time
state%time = read_time(file_id)

! Read in the copy index
read(file_id, *) header
if(header /= 'fcopy')  then
   write(*, *) 'Error: expected "copy" in input_diagnostics'
   stop
endif

read(file_id, *) copy_index

! Read in the state vector
read(file_id, *) state%state_vector

end subroutine input_diagnostics




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



subroutine adv_1step(state)
!----------------------------------------------------------------------
! subroutine adv_1step(state)
!
! Does single time step advance for lorenz 96 model
! using four-step rk time step

implicit none

type(assim_model_type), intent(inout) :: state

real(r8), dimension(model_size) :: x, x1, x2, x3, x4, dx, inter
type(time_type) :: time
integer :: i

! Get the state vector from the state; More copying?
x = get_model_state_vector(state)

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

! Update the time
time = get_model_time(state)
time = time + time_step

! Load x and the updated time back into output
call set_model_time(state, time)
call set_model_state_vector(state, x)

end subroutine adv_1step



!
!===================================================================
! End of assim_model_mod
!===================================================================
!
end module assim_model_mod

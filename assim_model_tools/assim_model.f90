module assim_model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

! This module is used to wrap around the basic portions of existing dynamical models to
! add capabilities needed by the standard assimilation methods.

! NEED TO ADD ON ONLY CLAUSES
use location_mod, only : location_type, write_location, set_location, get_dist, &
   read_location, get_location
! I've had a problem with putting in the only for time_manager on the pgf90 compiler (JLA).
use time_manager_mod
use utilities_mod, only : get_unit
use types_mod
use model_mod, only : get_model_size, static_init_model, get_state_meta_data, &
   get_model_time_step, model_interpolate, init_conditions, init_time, adv_1step, &
   end_model

private

public static_init_assim_model, init_diag_output, get_model_size, get_closest_state_time_to, &
   get_initial_condition, get_state_meta_data, get_close_states, get_num_close_states, &
   get_model_time, get_model_state_vector, copy_assim_model, advance_state, interpolate, &
   set_model_time, set_model_state_vector, write_state_restart, read_state_restart, &
   output_diagnostics, end_assim_model, assim_model_type, init_diag_input, input_diagnostics, &
   get_diag_input_copy_meta_data, init_assim_model, get_state_vector_ptr


! Eventually need to be very careful to implement this to avoid state vector copies which
! will be excruciatingly costly (storage at least) in big models.
type assim_model_type
!   private
   real(r8), pointer :: state_vector(:)
   type(time_type) :: time
end type assim_model_type

! Permanent class storage for model_size
integer :: model_size


type(time_type) :: time_step

contains

!======================================================================


subroutine init_assim_model(state)
!----------------------------------------------------------------------
!
! Allocates storage for an instance of an assim_model_type. With this
! implementation, need to be VERY careful about assigment and maintaining
! permanent storage locations. Need to revisit the best way to do 
! assim_model_copy below.

implicit none

type(assim_model_type), intent(inout) :: state

! Get the model_size from the model
model_size = get_model_size()

allocate(state%state_vector(model_size))

end subroutine init_assim_model




subroutine static_init_assim_model()
!----------------------------------------------------------------------
! subroutine static_init_assim_model()
!
! Initializes class data for the assim_model. Also calls the static
! initialization for the underlying model. So far, this simply 
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

! Call the underlying model's static initialization
call static_init_model()

end subroutine static_init_assim_model




function init_diag_output(file_name, global_meta_data, &
   copies_of_field_per_time, meta_data_per_copy)
!---------------------------------------------------------------------
!
! Initializes a diagnostic output file. Should be NetCDF shortly but
! for now just opens file and dumps stuff. A file_id is returned which
! is simply implemented as an integer unit number for now. 

implicit none

integer :: init_diag_output
character(len = *), intent(in) :: file_name, global_meta_data
integer, intent(in) :: copies_of_field_per_time
character(len = *), intent(in) :: meta_data_per_copy(copies_of_field_per_time)

integer :: i, model_size
type(location_type) :: state_loc

init_diag_output = get_unit()
open(unit = init_diag_output, file = file_name)
!!!init_diag_output = open_file(file_name)
write(init_diag_output, *) global_meta_data

! Write the model size
model_size = get_model_size()
write(init_diag_output, *) model_size

! Write number of copies of field per time plus the meta data per copy
write(init_diag_output, *) copies_of_field_per_time
do i = 1, copies_of_field_per_time
   write(init_diag_output, *) i, meta_data_per_copy(i)
end do

! Will need other metadata, too; Could be as simple as writing locations
write(init_diag_output, *) 'locat'
do i = 1, model_size
   call get_state_meta_data(i, state_loc)
   call write_location(init_diag_output, state_loc)
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




  function get_closest_state_time_to(assim_model, time)
!----------------------------------------------------------------------
!
! Returns the time closest to the given time that the model can reach
! with its state. Initial implementation just assumes fixed timestep.
! Need to describe potentially more general time-stepping capabilities
! from the underlying model in the long run.

implicit none

type(assim_model_type), intent(in) :: assim_model
type(time_type), intent(in) :: time
type(time_type) :: get_closest_state_time_to

type(time_type) :: model_time, delta_time, time_step

! CAREFUL WITH FLOATING POINT DIVISION AND COMPARISONS

! Get the model time step capabilities
time_step = get_model_time_step()

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
! Initial conditions . This returns an initial assim_model_type
! which includes both a state vector and a time. Design of exactly where this 
! stuff should come from is still evolving (12 July, 2002) but for now can 
! start at time offset 0 with the initial state .
! Need to carefully coordinate this with the times for observations.

implicit none

type(assim_model_type), intent(inout) :: x

call init_conditions(x%state_vector)

call init_time(x%time)

end subroutine get_initial_condition




subroutine get_close_states(location, radius, number, indices, dist)
!---------------------------------------------------------------------
! subroutine get_close_states(location, radius, number, indices)
!
! Returns a list of indices for model state vector points that are
! within distance radius of the location. Might want to add an option
! to return the distances, too. This is written in a model independent
! form at present, hence it is in assim_model_mod. HOWEVER, for
! efficiency in large models, this will have to be model specific at
! some point. At that time, need a way to test to see if this 
! generic form should be over loaded (how to do this in F90 ) by 
! some model specific method.

implicit none

type(location_type), intent(in) :: location
real(r8), intent(in) :: radius
integer, intent(out) :: number, indices(:)
real(r8), intent(out) :: dist(:)

type(location_type) :: state_loc
integer :: index, i
real(r8) :: this_dist

! For large models this will have to be VERY efficient; here can just search
index = 0
model_size = get_model_size
do i = 1, model_size
   call get_state_meta_data(i, state_loc)
   this_dist = get_dist(location, state_loc)
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

implicit none

integer :: get_num_close_states
type(location_type), intent(in) :: location
real(r8), intent(in) :: radius

type(location_type) :: state_loc
integer :: i

! For large models this will have to be VERY efficient; here can just search
get_num_close_states = 0
do i = 1, model_size
   call get_state_meta_data(i, state_loc)
! INTERESTING NOTE: Because of floating point round-off in comps
! this can give a 'variable' number of num close for certain obs, should fix
   if(get_dist(location, state_loc) < radius) get_num_close_states= get_num_close_states + 1
end do

end function get_num_close_states



function get_model_time(assim_model)
!-----------------------------------------------------------------------
!
! Returns the time component of a assim_model extended state.

implicit none

type(time_type) :: get_model_time
type(assim_model_type), intent(in) :: assim_model

get_model_time = assim_model%time

end function get_model_time



function get_state_vector_ptr(assim_model)
!------------------------------------------------------------------------
!
! Returns a pointer directly into the assim_model state vector storage.

implicit none

real(r8), pointer :: get_state_vector_ptr(:)
type(assim_model_type), intent(in) :: assim_model

get_state_vector_ptr => assim_model%state_vector

end function get_state_vector_ptr





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

do i = 1, get_model_size()
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

type(time_type) :: model_time, time_step

! NEED TO BE CAREFUL ABOUT FLOATING POINT TESTS: Being sloppy here

model_time = get_model_time(assim_model)
! Check for time error; use error handler when available
if(model_time > target_time) then
   write(*, *) 'Error in advance_state, target_time before model_time'
   stop
endif

! At some point probably need to push the determination of the time back
! into the model itself and out of assim_model
time_step = get_model_time_step()
do while(model_time < target_time)
   call adv_1step(assim_model%state_vector)
   model_time = model_time + time_step
end do

! Set the time to updated value
assim_model%time = model_time

end subroutine advance_state



function interpolate(x, location)
!---------------------------------------------------------------------
!
! Interpolates from the state vector in an assim_model_type to the
! location. Will need to be generalized for more complex state vector
! types. It might be better to be passing an assim_model_type with
! the associated time through here, but that requires changing the
! entire observation side of the class tree. Reconsider this at a 
! later date (JLA, 15 July, 2002).

implicit none

real(r8) :: interpolate
real(r8), intent(in) :: x(:)
type(location_type), intent(in) :: location

interpolate = model_interpolate(x, location)

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
if(size(state) /= get_model_size()) then
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

implicit none

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
type(assim_model_type), intent(inout) :: state
integer, intent(out) :: copy_index

character*5 :: header

! Read in the time
state%time = read_time(file_id)

! Read in the copy index
read(file_id, *) header
if(header /= 'fcopy')  then
   write(*, *) 'Error: expected "copy" in input_diagnostics'
   write(*, *) 'header read was: ', header
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

implicit none

call end_model()

end subroutine end_assim_model



function get_model_state_vector(assim_model)
!--------------------------------------------------------------------
!
! Returns the state vector component of an assim_model extended state.

real(r8) :: get_model_state_vector(model_size)
type(assim_model_type), intent(in) :: assim_model

get_model_state_vector = assim_model%state_vector

end function get_model_state_vector


!
!===================================================================
! End of assim_model_mod
!===================================================================
!
end module assim_model_mod

! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module assim_model_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
! This module is used to wrap around the basic portions of existing dynamical models to
! add capabilities needed by the standard assimilation methods.

use    types_mod, only : r8
use location_mod, only : location_type, get_dist, write_location, read_location, &
                         LocationDims, LocationName, LocationLName
! I've had a problem with putting in the only for time_manager on the pgf90 compiler (JLA).
use time_manager_mod, only : time_type, get_time, read_time, write_time, &
                             nc_write_calendar_atts, nc_get_tindex, &
                             operator(<), operator(>), operator(+), operator(-), &
                             operator(/), operator(*), operator(==), operator(/=) 
use utilities_mod, only : get_unit, file_exist, open_file, check_nml_error, close_file, &
                          register_module, error_handler, E_ERR, E_MSG, logfileunit
use     model_mod, only : get_model_size, static_init_model, get_state_meta_data, &
            get_model_time_step, model_interpolate, init_conditions, init_time, adv_1step, &
            end_model, model_get_close_states, nc_write_model_atts, nc_write_model_vars, &
            pert_model_state

implicit none
private

public :: static_init_assim_model, init_diag_output, get_model_size, get_closest_state_time_to, &
   get_initial_condition, get_state_meta_data, get_close_states, get_num_close_states, &
   get_model_time, get_model_state_vector, copy_assim_model, advance_state, interpolate, &
   set_model_time, set_model_state_vector, write_state_restart, read_state_restart, &
   output_diagnostics, end_assim_model, assim_model_type, init_diag_input, input_diagnostics, &
   get_diag_input_copy_meta_data, init_assim_model, get_state_vector_ptr, binary_restart_files, &
   finalize_diag_output, aoutput_diagnostics, aread_state_restart, aget_closest_state_time_to, &
   awrite_state_restart, Aadvance_state, pert_model_state

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"


! Eventually need to be very careful to implement this to avoid state vector copies which
! will be excruciatingly costly (storage at least) in big models. 
type assim_model_type
   private
   real(r8), pointer :: state_vector(:)
   type(time_type) :: time
   integer :: model_size       ! TJH request
   integer :: copyID           ! TJH request
! would like to include character string to indicate which netCDF variable --
! replace "state" in output_diagnostics ...
end type assim_model_type

! Permanent class storage for model_size
integer :: model_size


type(time_type) :: time_step

!-------------------------------------------------------------
! Namelist with default values
! binary_restart_files  == .true.  -> use unformatted file format. 
!                                     Full precision, faster, smaller,
!                                     but not as portable.
! binary_restart_files  == .false.  -> use ascii file format. 
!                                     Portable, but loses precision,
!                                     slower, and larger.

logical  :: binary_restart_files = .true.

namelist /assim_model_nml/ binary_restart_files
!-------------------------------------------------------------

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
state%model_size = model_size

end subroutine init_assim_model




subroutine static_init_assim_model()
!----------------------------------------------------------------------
! subroutine static_init_assim_model()
!
! Initializes class data for the assim_model. Also calls the static
! initialization for the underlying model. So far, this simply 
! is initializing the position of the state variables as location types.

implicit none

integer :: i, iunit, ierr, io

! First thing to do is echo info to logfile ... 

call register_module(source, revision, revdate)

! Read the namelist input
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = assim_model_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'assim_model_nml')
   enddo
 11 continue
   call close_file(iunit)
endif

! Record the namelist values used for the run ... 
call error_handler(E_MSG,'static_init_assim_model','assim_model namelist values: ',' ',' ',' ')
write(logfileunit, nml=assim_model_nml)

! Call the underlying model's static initialization
call static_init_model()

end subroutine static_init_assim_model



function init_diag_output(FileName, global_meta_data, &
                  copies_of_field_per_time, meta_data_per_copy) result(ncFileID)
!--------------------------------------------------------------------------------
!
! Typical sequence:
! NF90_OPEN             ! create netCDF dataset: enter define mode
!    NF90_def_dim       ! define dimenstions: from name and length
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset
!
! Time is a funny beast ... 
! Many packages decode the time:units attribute to convert the offset to a calendar
! date/time format. Using an offset simplifies many operations, but is not the
! way we like to see stuff plotted. The "approved" calendars are:
! gregorian or standard 
!      Mixed Gregorian/Julian calendar as defined by Udunits. This is the default. 
!  noleap   Modern calendar without leap years, i.e., all years are 365 days long. 
!  360_day  All years are 360 days divided into 30 day months. 
!  julian   Julian calendar. 
!  none     No calendar. 
!
! location is another one ...
!

use typeSizes
use netcdf
implicit none

character(len=*), intent(in) :: FileName, global_meta_data
integer,          intent(in) :: copies_of_field_per_time
character(len=*), intent(in) :: meta_data_per_copy(copies_of_field_per_time)
integer                      :: ncFileID

integer             :: i, metadata_length
type(location_type) :: state_loc

integer ::   MemberDimID,   MemberVarID     ! for each "copy" or ensemble member
integer ::     TimeDimID,     TimeVarID
integer :: LocationDimID, LocationVarID
integer :: MetadataDimID, MetadataVarID

if(.not. byteSizesOK()) then
    call error_handler(E_ERR,'init_diag_output', &
   'Compiler does not support required kinds of variables.',source,revision,revdate) 
end if

metadata_length = LEN(meta_data_per_copy(1))

! Create the file
!!!call check(nf90_create(path = trim(FileName)//".nc", cmode = nf90_clobber, ncid = ncFileID))
call check(nf90_create(path = trim(FileName)//".nc", cmode = nf90_share, ncid = ncFileID))

! Define the dimensions
call check(nf90_def_dim(ncid=ncFileID, &
             name="metadatalength", len = metadata_length,        dimid = metadataDimID))

call check(nf90_def_dim(ncid=ncFileID, &
             name="locationrank",   len = LocationDims,           dimid = LocationDimID))

call check(nf90_def_dim(ncid=ncFileID, &
             name="copy",           len=copies_of_field_per_time, dimid = MemberDimID))

call check(nf90_def_dim(ncid=ncFileID, &
             name="time",           len = nf90_unlimited,         dimid = TimeDimID))

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "title", global_meta_data))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "assim_model_source", source ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "assim_model_revision", revision ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "assim_model_revdate", revdate ))

!-------------------------------------------------------------------------------
! Create variables and attributes.
! The locations are part of the model (some models have multiple grids).
! They are written by model_mod:nc_write_model_atts
!-------------------------------------------------------------------------------

!    Copy ID
call check(nf90_def_var(ncid=ncFileID, name="copy", xtype=nf90_int, dimids=MemberDimID, &
                                                                    varid=MemberVarID))
call check(nf90_put_att(ncFileID, MemberVarID, "long_name", "ensemble member or copy"))
call check(nf90_put_att(ncFileID, MemberVarID, "units",     "nondimensional") )
call check(nf90_put_att(ncFileID, MemberVarID, "valid_range", (/ 1, copies_of_field_per_time /)))


!    Metadata for each Copy
call check(nf90_def_var(ncid=ncFileID,name="CopyMetaData", xtype=nf90_char,    &
                        dimids = (/ metadataDimID, MemberDimID /),  varid=metadataVarID))
call check(nf90_put_att(ncFileID, metadataVarID, "long_name",       &
                        "Metadata for each copy/member"))

!    Time -- the unlimited dimension
call check(nf90_def_var(ncFileID, name="time", xtype=nf90_double, dimids=TimeDimID, &
                                                                  varid =TimeVarID) )
i = nc_write_calendar_atts(ncFileID, TimeVarID)     ! comes from time_manager_mod
if ( i < 0 ) then
   print *,'assim_model_mod:nc_write_calendar_atts  bombed ', i
else if ( i > 0 ) then
   print *,'assim_model_mod:nc_write_calendar_atts  bombed ', i
endif

!-------------------------------------------------------------------------------
! Leave define mode so we can fill
!-------------------------------------------------------------------------------
call check(nf90_enddef(ncfileID))

!-------------------------------------------------------------------------------
! Fill the coordinate variables.
! The time variable is filled as time progresses.
!-------------------------------------------------------------------------------

call check(nf90_put_var(ncFileID, MemberVarID,   (/ (i,i=1,copies_of_field_per_time) /) ))
call check(nf90_put_var(ncFileID, metadataVarID, meta_data_per_copy ))

!-------------------------------------------------------------------------------
! sync to disk, but leave open
!-------------------------------------------------------------------------------

call check(nf90_sync(ncFileID))

!-------------------------------------------------------------------------------
! Define the model-specific components
!-------------------------------------------------------------------------------

i =  nc_write_model_atts( ncFileID )
if ( i < 0 ) then
   print *,'assim_model_mod:nc_write_model_atts  bombed ', i
else if ( i > 0 ) then
   print *,'assim_model_mod:nc_write_model_atts  bombed ', i
endif

!-------------------------------------------------------------------------------
call check(nf90_sync(ncFileID))               ! sync to disk, but leave open
!-------------------------------------------------------------------------------

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(status)
    integer, intent ( in) :: status
    if(status /= nf90_noerr) call error_handler(E_ERR,'init_diag_output', &
                                  trim(nf90_strerror(status)), source, revision, revdate)
  end subroutine check  

end function init_diag_output



function finalize_diag_output(ncFileID) result(ierr)
!--------------------------------------------------------------------------------
!

use netcdf
implicit none

integer, intent(in) :: ncFileID
integer             :: ierr

ierr = NF90_close(ncFileID)

end function finalize_diag_output



function init_diag_input(file_name, global_meta_data, model_size, copies_of_field_per_time)
!--------------------------------------------------------------------------
!
! Initializes a model state diagnostic file for input. A file id is
! returned which for now is just an integer unit number.

implicit none

integer :: init_diag_input
character(len = *), intent(in)  :: file_name
character(len = *), intent(out) :: global_meta_data
integer,            intent(out) :: model_size, copies_of_field_per_time

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

character(len=129) :: header, errstring
integer :: i, j

! Should have space checks, etc here
! Read the meta data associated with each copy
do i = 1, num_copies
   read(file_id, *) j, meta_data_per_copy(i)
end do

! Will need other metadata, too; Could be as simple as writing locations
read(file_id, *) header
if(header /= 'locat') then
   write(errstring,*)'expected to read "locat" got ',trim(adjustl(header))
   call error_handler(E_ERR,'get_diag_input_copy_meta_data', &
        errstring, source, revision, revdate)
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

type(time_type) :: model_time

model_time = assim_model%time

get_closest_state_time_to = aget_closest_state_time_to(model_time, time)

end function get_closest_state_time_to



function aget_closest_state_time_to(model_time, time)
!----------------------------------------------------------------------
!
! Returns the time closest to the given time that the model can reach
! with its state. Initial implementation just assumes fixed timestep.
! Need to describe potentially more general time-stepping capabilities
! from the underlying model in the long run.

implicit none

type(time_type), intent(in) :: model_time, time
type(time_type) :: aget_closest_state_time_to

type(time_type) :: delta_time, time_step

character(len=129) :: errstring
integer :: is1,is2,id1,id2

! Get the model time step capabilities
time_step = get_model_time_step()

if(model_time > time) then
   call get_time(model_time,is1,id1)
   call get_time(time,is2,id2)
   write(errstring, *)'model time (',is1,id1,') > time (',is2,id2,')'
   call error_handler(E_ERR,'aget_closest_state_time_to', errstring, source, revision, revdate)
endif

delta_time = time - model_time

aget_closest_state_time_to = (delta_time / time_step) * time_step + model_time

end function aget_closest_state_time_to



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

call aget_initial_condition(x%time, x%state_vector)

end subroutine get_initial_condition



subroutine aget_initial_condition(time, x)
!----------------------------------------------------------------------
! function get_initial_condition()
!
! Initial conditions . This returns an initial assim_model_type
! which includes both a state vector and a time. Design of exactly where this 
! stuff should come from is still evolving (12 July, 2002) but for now can 
! start at time offset 0 with the initial state .
! Need to carefully coordinate this with the times for observations.

implicit none

type(time_type), intent(out) :: time
real(r8), intent(inout) :: x(:)

call init_conditions(x)

call init_time(time)

end subroutine aget_initial_condition




subroutine get_close_states(location, radius, numinds, indices, dist)
!---------------------------------------------------------------------
! subroutine get_close_states(location, radius, numinds, indices)
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
real(r8),            intent(in) :: radius
integer,             intent(out) :: numinds, indices(:)
real(r8),            intent(out) :: dist(:)

type(location_type) :: state_loc
integer :: indx, i
real(r8) :: this_dist

! If model provides a working get_close_states, use it; otherwise search
! Direct use of model dependent stuff, needs to be automated (F90 can't do this
call model_get_close_states(location, radius, numinds, indices, dist)

! If numinds returns as -1, not implemented
if(numinds == -1) then
   indx = 0
   model_size = get_model_size()
   do i = 1, model_size
      call get_state_meta_data(i, state_loc)
      this_dist = get_dist(location, state_loc)
      if(this_dist < radius) then
         indx = indx + 1
         if(indx <= size(indices)) indices(indx) = i
         if(indx <= size(dist)) dist(indx) = this_dist
      end if
   end do
   numinds = indx
endif

! If size has overflowed, indicate this with negative size return
if(numinds > size(indices) .or. numinds > size(dist)) then
   numinds = -1 * numinds
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
real(r8),            intent(in) :: radius

type(location_type) :: state_loc
integer             :: i, indices(1)
real(r8)            :: dist(1)


! call direct model get close with storage that is too 
! small and get size from this
! model_get_close_states returns -1 if it is not implemented
call model_get_close_states(location, radius, get_num_close_states, indices, dist)

if(get_num_close_states == -1) then
   ! Do exhaustive search
   get_num_close_states = 0
   do i = 1, model_size
      call get_state_meta_data(i, state_loc)
      ! INTERESTING NOTE: Because of floating point round-off in comps
      ! this can give a 'variable' number of num close for certain obs, should fix
      if(get_dist(location, state_loc) < radius) get_num_close_states= get_num_close_states + 1
   end do

endif
   
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
type(assim_model_type), intent(in)  :: model_in

integer :: i

! Need to make sure to copy the actual storage and not just the pointer (verify)
model_out%time       = model_in%time
model_out%model_size = model_in%model_size

do i = 1, model_in%model_size
   model_out%state_vector(i) = model_in%state_vector(i)
end do

end subroutine copy_assim_model





subroutine advance_state(assim_model, num, target_time, asynch)
!-----------------------------------------------------------------------
!
! Advances the model extended state until time is equal (within roundoff?)
! of the target_time. For L96 this is relatively straightforward with 
! fixed time steps, etc.
! WARNING: Revision for ifc efficiency results in a copy of ensemble storage
! here if ensembles are advance using assim_model structure. May want to 
! duplicate code below for efficiency at some point if anyone's going to
! use this.

implicit none

integer,                intent(in)    :: num
type(assim_model_type), intent(inout) :: assim_model(num)
type(time_type),        intent(in)    :: target_time
integer,                intent(in)    :: asynch

type(time_type) :: model_time(num)
real(r8) :: model_state(num, size(assim_model(1)%state_vector))
integer  :: i

! Copy the times and the states to array storage
do i = 1, num
   model_time(i)     = assim_model(i)%time
   model_state(i, :) = assim_model(i)%state_vector
end do

call Aadvance_state(model_time, model_state, num, target_time, asynch)

! Now put the times and states that are updated back into their storage
do i = 1, num
   assim_model(i)%time         = model_time(i)
   assim_model(i)%state_vector = model_state(i, :)
end do

end subroutine advance_state




subroutine Aadvance_state(model_time, model_state, num, target_time, asynch)
!-----------------------------------------------------------------------
!
! Advances the model extended state until time is equal (within roundoff?)
! of the target_time. For L96 this is relatively straightforward with 
! fixed time steps, etc.

implicit none

integer,         intent(in)    :: num
type(time_type), intent(inout) :: model_time(num)
real(r8),        intent(inout) :: model_state(:, :)
type(time_type), intent(in)    :: target_time
integer,         intent(in)    :: asynch

type(time_type) :: time_step

integer :: seconds, days, i, control_unit, ic_file_unit, ud_file_unit

character(len = 26), dimension(num) :: ic_file_name, ud_file_name 
character(len = 128) :: input_string
character(len = 129) :: errstring
integer :: is1,is2,id1,id2

! If none of these needs advancing just return
do i = 1, num
   if(model_time(i) /= target_time) goto 10
end do
return

! Loop through each model state and advance
10 do i = 1, num

   ! Check for time error; use error handler when available
   ! TJH -- cannot write time_type to stdout -- they have private
   !        components and the compiler balks. time_manager_mod
   !        provides a routine to return the seconds and days
   !        from a time type. 

   if(model_time(i) > target_time) then
      call get_time(model_time(i),is1,id1)
      call get_time(target_time,is2,id2)
      write(errstring,*)'target time ',is2,id2,' is before model_time ',is1,id1
      call error_handler(E_ERR,'Aadvance_state', errstring, source, revision, revdate)
   endif

! Two blocks here for now, one for single executable, one for asynch multiple execs

!------------- Block for single executable ----------------------------
! At some point probably need to push the determination of the time back
! into the model itself and out of assim_model

   if(asynch == 0) then

      time_step = get_model_time_step()
      call get_time(time_step, seconds, days)

      do while(model_time(i) < target_time)

         call adv_1step(model_state(i, :), model_time(i))
         model_time(i) = model_time(i) + time_step
         call get_time(model_time(i), seconds, days)

      end do

   !-------------- End single executable block ------------------------
   else 
   !-------------- Block for multiple asynch executables --------------

   ! Loop to write out state for each member to separate file, preface with target time
      if(i < 10) then
         write(ic_file_name(i), 11) 'assim_model_state_ic', i
         write(ud_file_name(i), 11) 'assim_model_state_ud', i
      else if(i < 100) then
         write(ic_file_name(i), 21) 'assim_model_state_ic', i
         write(ud_file_name(i), 21) 'assim_model_state_ud', i
      else if(i < 1000) then
         write(ic_file_name(i), 31) 'assim_model_state_ic', i
         write(ud_file_name(i), 31) 'assim_model_state_ud', i
      else if(i < 10000) then
         write(ic_file_name(i), 41) 'assim_model_state_ic', i
         write(ud_file_name(i), 41) 'assim_model_state_ud', i
      else 
         write(errstring,*)'Trying to use ',num,' model states -- too many.'
         call error_handler(E_MSG,'Aadvance_state',errstring,source,revision,revdate)
         call error_handler(E_ERR,'Aadvance_state','Use less than 10000 model states.',source,revision,revdate)
      endif

 11   format(a21, i1)
 21   format(a21, i2)
 31   format(a21, i3)
 41   format(a21, i4)
      write(*, *) 'ic and ud files ', i, ic_file_name(i), ud_file_name(i)
      write(logfileunit, *) 'ic and ud files ', i, ic_file_name(i), ud_file_name(i)

      ! Output the destination time followed by assim_model_state 
      ! Write the time to which to advance
      ! Write the assim model extended state   

      ic_file_unit = get_unit()
      if ( binary_restart_files ) then
            open(unit = ic_file_unit, file = ic_file_name(i),      form = 'unformatted')
            call write_time(ic_file_unit, target_time,                    'unformatted')
            call awrite_state_restart(model_time(i), model_state(i, :), ic_file_unit, 'unformatted')
      else
            open(unit = ic_file_unit, file = ic_file_name(i))
            call write_time(ic_file_unit, target_time)
            call awrite_state_restart(model_time(i), model_state(i, :), ic_file_unit)
      endif
      close(ic_file_unit)

   endif

   !-------------- End of multiple async executables block ------------

end do


! Also need synchronization block at the end for the asynch

if(asynch /= 0) then

   ! Write out the file names to a control file

   control_unit = get_unit()
   open(unit = control_unit, file = 'filter_control')
   write(control_unit, *) num
   do i = 1, num
      write(control_unit, '(a26)' ) ic_file_name(i)
      write(control_unit, '(a26)' ) ud_file_name(i)
   end do
   close(control_unit)

   if(asynch == 1) then

      ! We sleep 1 second here for the following reason:
      ! The first time async_filter.csh is executed,
      ! it removes the file async_may_go that may be left over.
      ! With small models, perfect_model_obs and filter can come
      ! to this point before the script async_filter.csh has time to
      ! remove the file async_may_go.

      call system('sleep 1')

      ! Create the file async_may_go to allow the async model integrations
      control_unit = get_unit()
      open(unit = control_unit, file = 'async_may_go')
      call write_time(control_unit, target_time)
      close(control_unit)

! Suspend on a read from standard in for integer value
      do
         read(*, '(a80)') input_string
         if(trim(input_string) == 'All_done:Please_proceed') exit
! Following line can allow diagnostic pass through of output
         write(*, *) input_string
      end do

      ! This sleep is necessary to insure that all files coming from
      ! remote systems are completely written.

      call system('sleep 1')

   elseif(asynch == 2) then

      call system('./advance_ens.csh ; sleep 1')

   else

      write(errstring,*)'input.nml - async is ',asynch,' must be 0, 1, or 2' 
      call error_handler(E_ERR,'Aadvance_state', errstring, source, revision, revdate)

   endif

   write(*, *) 'got clearance to proceed in Aadvance_state'

   ! All should be done, read in the states and proceed
   do i = 1, num
      ud_file_unit = get_unit()
      if ( binary_restart_files ) then
         open(unit = ud_file_unit, file = ud_file_name(i),     form = 'unformatted')
         call aread_state_restart(model_time(i), model_state(i, :), ud_file_unit, 'unformatted')
      else
         open(unit = ud_file_unit, file = ud_file_name(i))
         call aread_state_restart(model_time(i), model_state(i, :), ud_file_unit)
      endif
      close(ud_file_unit)
   end do

end if

end subroutine Aadvance_state



function interpolate(x, location, loctype)
!---------------------------------------------------------------------
!
! Interpolates from the state vector in an assim_model_type to the
! location. Will need to be generalized for more complex state vector
! types. It might be better to be passing an assim_model_type with
! the associated time through here, but that requires changing the
! entire observation side of the class tree. Reconsider this at a 
! later date (JLA, 15 July, 2002). loctype for now is an integer that
! specifies what sort of variable from the model should be interpolated.

implicit none

real(r8) :: interpolate
real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: loctype

interpolate = model_interpolate(x, location, loctype)

end function interpolate



subroutine set_model_time(assim_model, time)
!-----------------------------------------------------------------------
!
! Sets the time in an assim_model type

implicit none

type(assim_model_type), intent(inout) :: assim_model
type(time_type),        intent(in)    :: time

assim_model%time = time

end subroutine set_model_time



subroutine set_model_state_vector(assim_model, state)
!-----------------------------------------------------------------------
!
! Sets the state vector part of an assim_model_type

implicit none

type(assim_model_type), intent(inout) :: assim_model
real(r8),               intent(in)    :: state(:)

character(len=129) :: errstring

! Check the size for now
if(size(state) /= get_model_size()) then
   write(errstring,*)'state vector has length ',size(state), &
                     ' model size (',get_model_size(),') does not match.'
   call error_handler(E_ERR,'set_model_state_vector', errstring, source, revision, revdate)
endif

assim_model%state_vector = state

end subroutine set_model_state_vector



subroutine write_state_restart(assim_model, funit, fform)
!----------------------------------------------------------------------
!
! Write a restart file given a model extended state and a unit number 
! opened to the restart file. (Need to reconsider what is passed to 
! identify file or if file can even be opened within this routine).

implicit none

type (assim_model_type), intent(in)           :: assim_model
integer,                 intent(in)           :: funit
character(len=*),        intent(in), optional :: fform

if(present(fform)) then
   call awrite_state_restart(assim_model%time, assim_model%state_vector, funit, fform)
else
   call awrite_state_restart(assim_model%time, assim_model%state_vector, funit)
endif

end subroutine write_state_restart




subroutine awrite_state_restart(model_time, model_state, funit, fform)
!----------------------------------------------------------------------
!
! Write a restart file given a model extended state and a unit number 
! opened to the restart file. (Need to reconsider what is passed to 
! identify file or if file can even be opened within this routine).

implicit none

type(time_type), intent(in)                   :: model_time
real(r8), intent(in)                          :: model_state(:)
integer,                 intent(in)           :: funit
character(len=*),        intent(in), optional :: fform
integer           :: days, seconds
character(len=32) :: fileformat

fileformat = "ascii"
if (present(fform)) fileformat = trim(adjustl(fform))


! This needs to be done more carefully, consider this an extended stub
! Write time first
! Write the state vector


SELECT CASE (fileformat)
   CASE ("unf","UNF","unformatted","UNFORMATTED")
      call write_time(funit, model_time, "unformatted")
      write(funit) model_state
   CASE DEFAULT
      call write_time(funit, model_time)
      write(funit, *) model_state
END SELECT  

end subroutine awrite_state_restart



subroutine read_state_restart(assim_model, funit, fform)
!----------------------------------------------------------------------
!
! Read a restart file given a unit number (see write_state_restart)

implicit none

type(assim_model_type), intent(out)          :: assim_model
integer,                intent(in)           :: funit
character(len=*),       intent(in), optional :: fform

if(present(fform)) then
   call aread_state_restart(assim_model%time, assim_model%state_vector, funit, fform)
else
   call aread_state_restart(assim_model%time, assim_model%state_vector, funit)
endif

end subroutine read_state_restart




subroutine aread_state_restart(model_time, model_state, funit, fform)
!----------------------------------------------------------------------
!
! Read a restart file given a unit number (see write_state_restart)

implicit none

type(time_type), intent(out)                 :: model_time
real(r8), intent(out)                        :: model_state(:)
integer,                intent(in)           :: funit
character(len=*),       intent(in), optional :: fform

integer           :: seconds, days
character(len=32) :: fileformat


print *,'assim_model_mod:aread_state_restart ... reading from unit',funit


fileformat = "ascii"
if (present(fform)) fileformat = trim(adjustl(fform))

! Read the time
! Read the state vector

SELECT CASE (fileformat)
   CASE ("unf","UNF","unformatted","UNFORMATTED")
      model_time = read_time(funit, form = "unformatted")
      read(funit) model_state
   CASE DEFAULT
      model_time = read_time(funit)
      read(funit, *) model_state
END SELECT  

end subroutine aread_state_restart



subroutine output_diagnostics(ncFileID, state, copy_index)
!-------------------------------------------------------------------
! Outputs the "state" to the supplied netCDF file. 
!
! the time, and an optional index saying which
! copy of the metadata this state is associated with.
!
! ncFileID       the netCDF file identifier
! state          the copy of the state vector
! copy_index     which copy of the state vector (ensemble member ID)
!
! TJH 28 Aug 2002 original netCDF implementation 
! TJH  7 Feb 2003 [created time_manager_mod:nc_get_tindex] 
!     substantially modified to handle time in a much better manner
! TJH 24 Jun 2003 made model_mod do all the netCDF writing.
!                 Still need an error handler for nc_write_model_vars
!

implicit none

integer,                intent(in) :: ncFileID
type(assim_model_type), intent(in) :: state
integer, optional,      intent(in) :: copy_index

if(present(copy_index)) then
   call aoutput_diagnostics(ncFileID, state%time, state%state_vector, copy_index)
else
   call aoutput_diagnostics(ncFileID, state%time, state%state_vector)
endif

end subroutine output_diagnostics




subroutine aoutput_diagnostics(ncFileID, model_time, model_state, copy_index)
!-------------------------------------------------------------------
! Outputs the "state" to the supplied netCDF file. 
!
! the time, and an optional index saying which
! copy of the metadata this state is associated with.
!
! ncFileID       the netCDF file identifier
! model_time     the time associated with the state vector
! model_state    the copy of the state vector
! copy_index     which copy of the state vector (ensemble member ID)
!
! TJH 28 Aug 2002 original netCDF implementation 
! TJH  7 Feb 2003 [created time_manager_mod:nc_get_tindex] 
!     substantially modified to handle time in a much better manner
! TJH 24 Jun 2003 made model_mod do all the netCDF writing.
!                 Still need an error handler for nc_write_model_vars
!      

use typeSizes
use netcdf
implicit none

integer,           intent(in) :: ncFileID
type(time_type),   intent(in) :: model_time
real(r8),          intent(in) :: model_state(:)
integer, optional, intent(in) :: copy_index

integer :: i,ierr, timeindex, copyindex

character(len=129) :: errstring
integer :: is1,id1

if (.not. present(copy_index) ) then     ! we are dependent on the fact
   copyindex = 1                         ! there is a copyindex == 1
else                                     ! if the optional argument is
   copyindex = copy_index                ! not specified, we'd better
endif                                    ! have a backup plan

timeindex = nc_get_tindex(ncFileID, model_time)
if ( timeindex < 0 ) then
   call get_time(model_time,is1,id1)
   write(errstring,*)'model time ',is1,id1,' not in netcdf file ',ncFileID 
   call error_handler(E_ERR,'aoutput_diagnostics', errstring, source, revision, revdate)
endif

! model_mod:nc_write_model_vars knows nothing about assim_model_types,
! so we must pass the components.

i = nc_write_model_vars(ncFileID, model_state, copyindex, timeindex) 

end subroutine aoutput_diagnostics




subroutine input_diagnostics(file_id, state, copy_index)
!------------------------------------------------------------------
!
! Reads in diagnostic state output from file_id for copy_index
! copy. Need to make this all more rigorously enforced.

implicit none

integer,                intent(in)    :: file_id
! MAYBE SHOULDN'T use assim model type here, but just state and time ?
type(assim_model_type), intent(inout) :: state
integer,                intent(out)   :: copy_index

call ainput_diagnostics(file_id, state%time, state%state_vector, copy_index)

end subroutine input_diagnostics



subroutine ainput_diagnostics(file_id, model_time, model_state, copy_index)
!------------------------------------------------------------------
!
! Reads in diagnostic state output from file_id for copy_index
! copy. Need to make this all more rigorously enforced.

implicit none

integer,         intent(in)    :: file_id
type(time_type), intent(inout) :: model_time
real(r8),        intent(inout) :: model_state(:)
integer,         intent(out)   :: copy_index

character(len=5)   :: header
character(len=129) :: errstring

! Read in the time
model_time = read_time(file_id)

! Read in the copy index
read(file_id, *) header
if(header /= 'fcopy')  then
   write(errstring,*)'expected "copy", got ',header
   call error_handler(E_ERR,'ainput_diagnostics', errstring, source, revision, revdate)
endif

read(file_id, *) copy_index

! Read in the state vector
read(file_id, *) model_state

end subroutine ainput_diagnostics




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

! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> This module is used to wrap around the basic portions of existing dynamical models to
!> add capabilities needed by the standard assimilation methods.
module assim_model_mod

use    types_mod, only : r8, digits12
use location_mod, only : location_type, read_location, LocationDims
use time_manager_mod, only : time_type, get_time, read_time, write_time,           &
                             THIRTY_DAY_MONTHS, JULIAN, GREGORIAN, NOLEAP,         &
                             operator(<), operator(>), operator(+), operator(-),   &
                             operator(/), operator(*), operator(==), operator(/=), &
                             get_calendar_type
use utilities_mod, only : get_unit, close_file, register_module, error_handler,    &
                          E_ERR, E_WARN, E_MSG, E_DBG, nmlfileunit,                &
                          dump_unit_attributes, find_namelist_in_file,             &
                          check_namelist_read, nc_check, do_nml_file, do_nml_term, &
                          find_textfile_dims, file_to_text, set_output,            &
                          ascii_file_format, set_output
use     model_mod, only : get_model_size, static_init_model, get_state_meta_data_distrib,  &
                          get_model_time_step, init_conditions,                            &
                          init_time, adv_1step, end_model, nc_write_model_atts,            &
                          nc_write_model_vars, pert_model_state,                           &
                          get_close_maxdist_init, get_close_obs_init,                      &
                          model_interpolate_distrib,                                       &
                          get_close_obs_distrib, pert_model_copies

use ensemble_manager_mod, only : ensemble_type

implicit none
private

public :: static_init_assim_model, init_diag_output, get_model_size,                       &
          get_closest_state_time_to, get_initial_condition, get_state_meta_data_distrib,   &
          get_model_time, get_model_state_vector, copy_assim_model,                        &
          set_model_time, set_model_state_vector,  &
          output_diagnostics, end_assim_model, assim_model_type, init_diag_input,          &
          input_diagnostics, get_diag_input_copy_meta_data, init_assim_model,              &
          finalize_diag_output, aoutput_diagnostics,                                       &
           aget_closest_state_time_to,            &
          pert_model_state, netcdf_file_type, nc_append_time, nc_write_calendar_atts,      &
          nc_get_tindex, get_model_time_step,       &
         adv_1step, aget_initial_condition, get_close_maxdist_init,        &
          get_close_obs_init, interpolate_distrib,                                         &
          get_close_obs_distrib, pert_model_copies

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


! Type to keep model state and time together
type assim_model_type
   private
   real(r8), pointer :: state_vector(:)
   type(time_type) :: time
   integer                           :: model_size
   integer                           :: copyID
! would like to include character string to indicate which netCDF variable --
! replace "state" in output_diagnostics ...
end type assim_model_type

!----------------------------------------------------------------
! output (netcdf) file descriptor
! basically, we want to keep a local mirror of the unlimited dimension
! coordinate variable (i.e. time) because dynamically querying it
! causes unacceptable performance degradation over "long" integrations.

type netcdf_file_type
   integer :: ncid                       ! the "unit" -- sorta
   integer :: Ntimes                     ! the current working length
   integer :: NtimesMAX                  ! the allocated length.
   real(digits12),  pointer :: rtimes(:) ! times -- as a 64bit (at least) real
   type(time_type), pointer :: times(:)  ! times -- as the models use
   character(len=80)        :: fname     ! filename ...
end type netcdf_file_type

! Permanent class storage for model_size
integer :: model_size

! Ensure init code is called exactly once
logical :: module_initialized = .false.


! Global storage for error string output
character(len = 129)  :: msgstring

!-------------------------------------------------------------
!> @todo check the performance of large file support.
!> Can we have large file support as the default?
logical  :: netCDF_large_file_support  = .false.

namelist /assim_model_nml/ netCDF_large_file_support
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

integer :: iunit, io

! only execute this code once, even if called multiple times.
if (module_initialized) return

! First thing to do is echo info to logfile ... 
call register_module(source, revision, revdate)
module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "assim_model_nml", iunit)
read(iunit, nml = assim_model_nml, iostat = io)
call check_namelist_read(iunit, io, "assim_model_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=assim_model_nml)
if (do_nml_term()) write(     *     , nml=assim_model_nml)



! Call the underlying model's static initialization
call static_init_model()

end subroutine static_init_assim_model



function init_diag_output(FileName, global_meta_data, &
                  copies_of_field_per_time, meta_data_per_copy, lagID) result(ncFileID)
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
integer, OPTIONAL,intent(in) :: lagID
type(netcdf_file_type)       :: ncFileID

integer :: i, metadata_length, nlines, linelen, createmode

integer ::   MemberDimID,   MemberVarID     ! for each "copy" or ensemble member
integer ::     TimeDimID,     TimeVarID
integer :: LocationDimID
integer :: MetadataDimID, MetadataVarID
integer ::   nlinesDimID,  linelenDimID, nmlVarID

character(len=129), allocatable, dimension(:) :: textblock

if ( .not. module_initialized ) call static_init_assim_model()

if(.not. byteSizesOK()) then
    call error_handler(E_ERR,'init_diag_output', &
   'Compiler does not support required kinds of variables.',source,revision,revdate) 
end if

metadata_length = LEN(meta_data_per_copy(1))

if ( netCDF_large_file_support ) then
   createmode = NF90_64BIT_OFFSET
else
   createmode = NF90_SHARE
endif

! Create the file
ncFileID%fname = trim(adjustl(FileName))//".nc"
call nc_check(nf90_create(path = trim(ncFileID%fname), cmode = createmode, ncid = ncFileID%ncid), &
              'init_diag_output', 'create '//trim(ncFileID%fname))

write(msgstring,*)trim(ncFileID%fname), ' is ncFileID ',ncFileID%ncid
call error_handler(E_MSG,'init_diag_output',msgstring,source,revision,revdate)

! Define the dimensions
call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
              name="metadatalength", len = metadata_length,        dimid = metadataDimID), &
              'init_diag_output', 'def_dim metadatalength '//trim(ncFileID%fname))

call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
              name="locationrank",   len = LocationDims,           dimid = LocationDimID), &
              'init_diag_output', 'def_dim locationrank '//trim(ncFileID%fname))

call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
              name="copy",           len=copies_of_field_per_time, dimid = MemberDimID), &
              'init_diag_output', 'def_dim copy '//trim(ncFileID%fname))

call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
              name="time",           len = nf90_unlimited,         dimid = TimeDimID), &
              'init_diag_output', 'def_dim time '//trim(ncFileID%fname))

!-------------------------------------------------------------------------------
! Find dimensions of namelist file ... will save it as a variable.
!-------------------------------------------------------------------------------

! All DART programs require input.nml, so it is unlikely this can fail, but
! still check and in this case, error out if not found.
call find_textfile_dims("input.nml", nlines, linelen)
if (nlines <= 0 .or. linelen <= 0) then
   call error_handler(E_MSG,'init_diag_output', &
                      'cannot open/read input.nml to save in diagnostic file', &
                      source,revision,revdate)
endif

allocate(textblock(nlines))
textblock = ''

call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
              name="NMLlinelen", len = LEN(textblock(1)), dimid = linelenDimID), &
              'init_diag_output', 'def_dim NMLlinelen '//trim(ncFileID%fname))

call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
              name="NMLnlines", len = nlines, dimid = nlinesDimID), &
              'init_diag_output', 'def_dim NMLnlines '//trim(ncFileID%fname))

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "title", global_meta_data), &
              'init_diag_output', 'put_att title '//trim(ncFileID%fname))
call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "assim_model_source", source ), &
              'init_diag_output', 'put_att assim_model_source '//trim(ncFileID%fname))
call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "assim_model_revision", revision ), &
              'init_diag_output', 'put_att assim_model_revision '//trim(ncFileID%fname))
call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "assim_model_revdate", revdate ), &
              'init_diag_output', 'put_att assim_model_revdate '//trim(ncFileID%fname))

if (present(lagID)) then
   call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "lag", lagID ), &
                 'init_diag_output', 'put_att lag '//trim(ncFileID%fname))

   write(*,*)'init_diag_output detected Lag is present'

endif 

!-------------------------------------------------------------------------------
! Create variables and attributes.
! The locations are part of the model (some models have multiple grids).
! They are written by model_mod:nc_write_model_atts
!-------------------------------------------------------------------------------

!    Copy ID
call nc_check(nf90_def_var(ncid=ncFileID%ncid, name="copy", xtype=nf90_int, &
              dimids=MemberDimID, varid=MemberVarID), 'init_diag_output', 'def_var copy')
call nc_check(nf90_put_att(ncFileID%ncid, MemberVarID, "long_name", "ensemble member or copy"), &
              'init_diag_output', 'long_name')
call nc_check(nf90_put_att(ncFileID%ncid, MemberVarID, "units",     "nondimensional"), &
              'init_diag_output', 'units')
call nc_check(nf90_put_att(ncFileID%ncid, MemberVarID, "valid_range", &
              (/ 1, copies_of_field_per_time /)), 'init_diag_output', 'put_att valid_range')


!    Metadata for each Copy
call nc_check(nf90_def_var(ncid=ncFileID%ncid,name="CopyMetaData", xtype=nf90_char,    &
              dimids = (/ metadataDimID, MemberDimID /),  varid=metadataVarID), &
              'init_diag_output', 'def_var CopyMetaData')
call nc_check(nf90_put_att(ncFileID%ncid, metadataVarID, "long_name",       &
              "Metadata for each copy/member"), 'init_diag_output', 'put_att long_name')

!    input namelist 
call nc_check(nf90_def_var(ncid=ncFileID%ncid,name="inputnml", xtype=nf90_char,    &
              dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
              'init_diag_output', 'def_var inputnml')
call nc_check(nf90_put_att(ncFileID%ncid, nmlVarID, "long_name",       &
              "input.nml contents"), 'init_diag_output', 'put_att input.nml')

!    Time -- the unlimited dimension
call nc_check(nf90_def_var(ncFileID%ncid, name="time", xtype=nf90_double, dimids=TimeDimID, &
              varid =TimeVarID), 'init_diag_output', 'def_var time' )
i = nc_write_calendar_atts(ncFileID, TimeVarID)     ! comes from time_manager_mod
if ( i /= 0 ) then
   write(msgstring, *)'nc_write_calendar_atts  bombed with error ', i
   call error_handler(E_MSG,'init_diag_output',msgstring,source,revision,revdate)
endif

! Create the time "mirror" with a static length. There is another routine
! to increase it if need be. For now, just pick something.
ncFileID%Ntimes    = 0
ncFileID%NtimesMAX = 1000
allocate(ncFileID%rtimes(ncFileID%NtimesMAX), ncFileID%times(ncFileID%NtimesMAX) )

!-------------------------------------------------------------------------------
! Leave define mode so we can fill
!-------------------------------------------------------------------------------
call nc_check(nf90_enddef(ncFileID%ncid), 'init_diag_output', 'enddef '//trim(ncFileID%fname))

!-------------------------------------------------------------------------------
! Fill the coordinate variables.
! Write the input namelist as a netCDF variable.
! The time variable is filled as time progresses.
!-------------------------------------------------------------------------------

call nc_check(nf90_put_var(ncFileID%ncid, MemberVarID, (/ (i,i=1,copies_of_field_per_time) /) ), &
              'init_diag_output', 'put_var MemberVarID')
call nc_check(nf90_put_var(ncFileID%ncid, metadataVarID, meta_data_per_copy ), &
              'init_diag_output', 'put_var metadataVarID')
 
call file_to_text("input.nml", textblock)

call nc_check(nf90_put_var(ncFileID%ncid, nmlVarID, textblock ), &
              'init_diag_output', 'put_var nmlVarID')

deallocate(textblock)

!-------------------------------------------------------------------------------
! sync to disk, but leave open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID%ncid), 'init_diag_output', 'sync '//trim(ncFileID%fname))               
!-------------------------------------------------------------------------------
! Define the model-specific components
!-------------------------------------------------------------------------------

i =  nc_write_model_atts( ncFileID%ncid )
if ( i /= 0 ) then
   write(msgstring, *)'nc_write_model_atts  bombed with error ', i
   call error_handler(E_MSG,'init_diag_output',msgstring,source,revision,revdate)
endif

!-------------------------------------------------------------------------------
! sync again, but still leave open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID%ncid), 'init_diag_output', 'sync '//trim(ncFileID%fname))               
!-------------------------------------------------------------------------------

end function init_diag_output



function finalize_diag_output(ncFileID) result(ierr)
!--------------------------------------------------------------------------------
!

use netcdf
implicit none

type(netcdf_file_type), intent(inout) :: ncFileID
integer             :: ierr

ierr = NF90_close(ncFileID%ncid)

ncFileID%fname     = "notinuse"
ncFileID%ncid      = -1
ncFileID%Ntimes    = -1
ncFileID%NtimesMax = -1
if(associated(ncFileID%rtimes)) deallocate(ncFileID%rtimes, ncFileID%times )

end function finalize_diag_output



function init_diag_input(file_name, global_meta_data, model_size, copies_of_field_per_time)
!--------------------------------------------------------------------------
!
! Initializes a model state diagnostic file for input. A file id is
! returned which for now is just an integer unit number.

implicit none

integer :: init_diag_input, io
character(len = *), intent(in)  :: file_name
character(len = *), intent(out) :: global_meta_data
integer,            intent(out) :: model_size, copies_of_field_per_time

if ( .not. module_initialized ) call static_init_assim_model()

init_diag_input = get_unit()
open(unit = init_diag_input, file = file_name, action = 'read', iostat = io)
if (io /= 0) then
   write(msgstring,*) 'unable to open diag input file ', trim(file_name), ' for reading'
   call error_handler(E_ERR,'init_diag_input',msgstring,source,revision,revdate)
endif

! Read meta data
read(init_diag_input, *, iostat = io) global_meta_data
if (io /= 0) then
   write(msgstring,*) 'unable to read expected character string from diag input file ', &
                       trim(file_name), ' for global_meta_data'
   call error_handler(E_ERR,'init_diag_input',msgstring,source,revision,revdate)
endif

! Read the model size
read(init_diag_input, *, iostat = io) model_size
if (io /= 0) then
   write(msgstring,*) 'unable to read expected integer from diag input file ', &
                       trim(file_name), ' for model_size'
   call error_handler(E_ERR,'init_diag_input',msgstring,source,revision,revdate)
endif

! Read the number of copies of field per time
read(init_diag_input, *, iostat = io) copies_of_field_per_time
if (io /= 0) then
   write(msgstring,*) 'unable to read expected integer from diag input file ', &
                       trim(file_name), ' for copies_of_field_per_time'
   call error_handler(E_ERR,'init_diag_input',msgstring,source,revision,revdate)
endif

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
integer :: i, j, io

! Should have space checks, etc here
! Read the meta data associated with each copy
do i = 1, num_copies
   read(file_id, *, iostat = io) j, meta_data_per_copy(i)
   if (io /= 0) then
      write(msgstring,*) 'error reading metadata for copy ', i, ' from diag file'
      call error_handler(E_ERR,'get_diag_input_copy_meta_data', &
                         msgstring,source,revision,revdate)
   endif
end do

! Will need other metadata, too; Could be as simple as writing locations
read(file_id, *, iostat = io) header
if (io /= 0) then
   write(msgstring,*) 'error reading header from diag file'
   call error_handler(E_ERR,'get_diag_input_copy_meta_data', &
                      msgstring,source,revision,revdate)
endif
if(header /= 'locat') then
   write(msgstring,*)'expected to read "locat" got ',trim(adjustl(header))
   call error_handler(E_ERR,'get_diag_input_copy_meta_data', &
                      msgstring, source, revision, revdate)
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

type(time_type)                    :: get_closest_state_time_to
type(assim_model_type), intent(in) :: assim_model
type(time_type), intent(in) :: time

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

type(time_type) :: aget_closest_state_time_to
type(time_type), intent(in) :: model_time, time

type(time_type) :: time_step

! Get the model time step capabilities
time_step = get_model_time_step()

if(model_time > time) then
   ! If model_time is past start of obs window, don't advance it
   aget_closest_state_time_to = model_time
   return
endif

aget_closest_state_time_to = model_time

do while((time_step + 2*aget_closest_state_time_to) < 2*time)
   aget_closest_state_time_to = aget_closest_state_time_to + time_step
enddo

end function aget_closest_state_time_to



subroutine get_initial_condition(x)
!----------------------------------------------------------------------
! function get_initial_condition()
!
! Initial conditions. This returns an initial assim_model_type
! which includes both a state vector and a time. Design of exactly where this 
! stuff should come from is still evolving (12 July, 2002) but for now can 
! start at time offset 0 with the initial state.
! Need to carefully coordinate this with the times for observations.

implicit none

type(assim_model_type), intent(inout) :: x

call aget_initial_condition(x%time, x%state_vector)

end subroutine get_initial_condition



subroutine aget_initial_condition(time, x)
!----------------------------------------------------------------------
! function get_initial_condition()
!
! Initial conditions. This returns an initial state vector and a time
! for use in an assim_model_type.  Design of exactly where this 
! stuff should come from is still evolving (12 July, 2002) but for now can 
! start at time offset 0 with the initial state.
! Need to carefully coordinate this with the times for observations.

implicit none

type(time_type), intent(out) :: time
real(r8),        intent(out) :: x(:)

call init_conditions(x)

call init_time(time)

end subroutine aget_initial_condition



function get_model_time(assim_model)
!-----------------------------------------------------------------------
!
! Returns the time component of a assim_model extended state.

implicit none

type(time_type)                    :: get_model_time
type(assim_model_type), intent(in) :: assim_model

get_model_time = assim_model%time

end function get_model_time





subroutine copy_assim_model(model_out, model_in)
!-------------------------------------------------------------------------
!
! Does a copy of assim_model, should be overloaded to =? Still need to be
! very careful about trying to limit copies of the potentially huge state
! vectors for big models.  Interaction with pointer storage?

implicit none

type(assim_model_type), intent(inout) :: model_out
type(assim_model_type), intent(in)    :: model_in

integer :: i

! Need to make sure to copy the actual storage and not just the pointer (verify)
model_out%time       = model_in%time
model_out%model_size = model_in%model_size

do i = 1, model_in%model_size
   model_out%state_vector(i) = model_in%state_vector(i)
end do

end subroutine copy_assim_model

!> Pass through routine to model interpolate
subroutine interpolate_distrib(location, loctype, istatus, expected_obs, state_ens_handle)
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

type(location_type),   intent(in)    :: location
integer,               intent(in)    :: loctype
integer,               intent(out)   :: istatus(:)
type(ensemble_type),   intent(in)    :: state_ens_handle
real(r8),              intent(out)   :: expected_obs(:)

istatus = 0

call model_interpolate_distrib(state_ens_handle, location, loctype, istatus, expected_obs)

end subroutine interpolate_distrib


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

! Check the size for now
if(size(state) /= get_model_size()) then
   write(msgstring,*)'state vector has length ',size(state), &
                     ' model size (',get_model_size(),') does not match.'
   call error_handler(E_ERR,'set_model_state_vector', msgstring, source, revision, revdate)
endif

assim_model%state_vector = state

end subroutine set_model_state_vector

!-------------------------------------------------------------------


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
! Note -- ncFileId may be modified -- the time mirror needs to
! track the state of the netCDF file. This must be "inout".

implicit none

type(netcdf_file_type), intent(inout) :: ncFileID
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
! Note -- ncFileId may be modified -- the time mirror needs to
! track the state of the netCDF file. This must be "inout".

use typeSizes
use netcdf
implicit none

type(netcdf_file_type), intent(inout) :: ncFileID
type(time_type),   intent(in) :: model_time
real(r8),          intent(in) :: model_state(:)
integer, optional, intent(in) :: copy_index

integer :: i, timeindex, copyindex
integer :: is1,id1

if (.not. present(copy_index) ) then     ! we are dependent on the fact
   copyindex = 1                         ! there is a copyindex == 1
else                                     ! if the optional argument is
   copyindex = copy_index                ! not specified, we'd better
endif                                    ! have a backup plan

timeindex = nc_get_tindex(ncFileID, model_time)
if ( timeindex < 0 ) then
   call get_time(model_time,is1,id1)
   write(msgstring,*)'model time (d,s)',id1,is1,' not in ',ncFileID%fname
   write(msgstring,'(''model time (d,s) ('',i8,i5,'') is index '',i6, '' in ncFileID '',i10)') &
          id1,is1,timeindex,ncFileID%ncid
   call error_handler(E_ERR,'aoutput_diagnostics', msgstring, source, revision, revdate)
endif

   call get_time(model_time,is1,id1)
   write(msgstring,'(''model time (d,s) ('',i8,i5,'') is index '',i6, '' in ncFileID '',i10)') &
          id1,is1,timeindex,ncFileID%ncid
   call error_handler(E_DBG,'aoutput_diagnostics', msgstring, source, revision, revdate)

! model_mod:nc_write_model_vars knows nothing about assim_model_types,
! so we must pass the components.
! No need to do this anymore
i = nc_write_model_vars(ncFileID%ncid, model_state, copyindex, timeindex)

end subroutine aoutput_diagnostics




subroutine input_diagnostics(file_id, state, copy_index)
!------------------------------------------------------------------
!
! Reads in diagnostic state output from file_id for copy_index
! copy. Need to make this all more rigorously enforced.

implicit none

integer,                intent(in)    :: file_id
type(assim_model_type), intent(inout) :: state
integer,                intent(out)   :: copy_index

! MAYBE SHOULDN'T use assim model type here, but just state and time ?

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

! Read in the time
model_time = read_time(file_id)

! Read in the copy index
read(file_id, *) header
if(header /= 'fcopy')  then
   write(msgstring,*)'expected "copy", got ',header
   call error_handler(E_ERR,'ainput_diagnostics', msgstring, source, revision, revdate)
endif

read(file_id, *) copy_index

! Read in the state vector
read(file_id, *) model_state

end subroutine ainput_diagnostics




subroutine end_assim_model()
!--------------------------------------------------------------------
!
! Closes down assim_model. For now, only thing to do is tell model to end.

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




function nc_append_time(ncFileID, time) result(lngth)
!------------------------------------------------------------------------
! The current time is appended to the "time" coordinate variable.
! The new length of the "time" variable is returned.
! 
! This REQUIRES that "time" is a coordinate variable AND it is the
! unlimited dimension. If not ... bad things happen.
!
! TJH Wed Aug 28 15:40:25 MDT 2002

use typeSizes
use netcdf
implicit none

type(netcdf_file_type), intent(inout) :: ncFileID
type(time_type), intent(in) :: time
integer                     :: lngth

integer  :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer  :: TimeVarID
integer  :: secs, days, ncid
real(digits12) :: realtime         ! gets promoted to nf90_double ...

character(len=NF90_MAX_NAME)          :: varname
integer                               :: xtype, ndims, nAtts
integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
character(len=129)                    :: msgstring

type(time_type), allocatable, dimension(:) :: temptime   ! only to reallocate mirror
real(digits12),  allocatable, dimension(:) :: tempRtime  ! only to reallocate mirror

lngth = -1 ! assume a bad termination

ncid = ncFileID%ncid

call nc_check(NF90_Inquire(ncid, nDimensions, nVariables, nAttributes, unlimitedDimID), &
              'nc_append_time', 'inquire '//ncFileID%fname)
call nc_check(NF90_Inq_Varid(ncid, "time", TimeVarID), 'nc_append_time', 'inq_varid time')
call nc_check(NF90_Inquire_Variable(ncid, TimeVarID, varname, xtype, ndims, dimids, nAtts), &
             'nc_append_time', 'inquire_variable time')

if ( ndims /= 1 ) call error_handler(E_ERR,'nc_append_time', &
           '"time" expected to be rank-1',source,revision,revdate)

if ( dimids(1) /= unlimitedDimID ) call error_handler(E_ERR,'nc_append_time', &
           'unlimited dimension expected to be slowest-moving',source,revision,revdate)

! make sure the mirror and the netcdf file are in sync
call nc_check(NF90_Inquire_Dimension(ncid, unlimitedDimID, varname, lngth ), &
           'nc_append_time', 'inquire_dimension unlimited')

if (lngth /= ncFileId%Ntimes) then
   write(msgstring,*)'netCDF file has length ',lngth,' /= mirror has length of ',ncFileId%Ntimes
   call error_handler(E_ERR,'nc_append_time', &
           'time mirror and netcdf file time dimension out-of-sync', &
           source,revision,revdate,text2=msgstring)
endif

! make sure the time mirror can handle another entry.
if ( lngth == ncFileID%NtimesMAX ) then   

   write(msgstring,*)'doubling mirror length of ',lngth,' of ',ncFileID%fname
   call error_handler(E_DBG,'nc_append_time',msgstring,source,revision,revdate)

   allocate(temptime(ncFileID%NtimesMAX), tempRtime(ncFileID%NtimesMAX)) 
   temptime   = ncFileID%times            ! preserve
   tempRtime = ncFileID%rtimes            ! preserve

   deallocate(ncFileID%times, ncFileID%rtimes)

   ncFileID%NtimesMAX = 2 * ncFileID%NtimesMAX  ! double length of exising arrays

   allocate(ncFileID%times(ncFileID%NtimesMAX), ncFileID%rtimes(ncFileID%NtimesMAX) )

   ncFileID%times(1:lngth)  = temptime    ! reinstate
   ncFileID%rtimes(1:lngth) = tempRtime   ! reinstate

   deallocate(temptime, tempRtime)

endif

call get_time(time, secs, days)         ! get time components to append
realtime = days + secs/86400.0_digits12 ! time base is "days since ..."
lngth           = lngth + 1             ! index of new time 
ncFileID%Ntimes = lngth                 ! new working length of time mirror

call nc_check(nf90_put_var(ncid, TimeVarID, realtime, start=(/ lngth /) ), &
           'nc_append_time', 'put_var time')

ncFileID%times( lngth) = time
ncFileID%rtimes(lngth) = realtime

write(msgstring,*)'ncFileID (',ncid,') : ',trim(adjustl(varname)), &
         ' (should be "time") has length ',lngth, ' appending t= ',realtime
call error_handler(E_DBG,'nc_append_time',msgstring,source,revision,revdate)

end function nc_append_time



function nc_get_tindex(ncFileID, statetime) result(timeindex)
!------------------------------------------------------------------------
! 
! We need to compare the time of the current assim_model to the 
! netcdf time coordinate variable (the unlimited dimension).
! If they are the same, no problem ...
! If it is earlier, we need to find the right index and insert ...
! If it is the "future", we need to add another one ...
! If it is in the past but does not match any we have, we're in trouble.
! The new length of the "time" variable is returned.
! 
! This REQUIRES that "time" is a coordinate variable AND it is the
! unlimited dimension. If not ... bad things happen.
!
! TJH  7 Feb 2003
!
! Revision by TJH 24 Nov 2003:
! A new array "times" has been added to mirror the times that are stored
! in the netcdf time coordinate variable. While somewhat unpleasant, it
! is SUBSTANTIALLY faster than reading the netcdf time variable at every
! turn -- which caused a geometric or exponential increase in overall 
! netcdf I/O. (i.e. this was really bad)
!
! The time mirror is maintained as a time_type, so the comparison with
! the state time uses the operators for the time_type. The netCDF file,
! however, has time units of a different convention. The times are
! converted only when appending to the time coordinate variable.    
!
! Revision by TJH 4 June 2004:
! Implementing a "file type" for output that contains a unique time
! mirror for each file.

use typeSizes
use netcdf

implicit none

type(netcdf_file_type), intent(inout) :: ncFileID
type(time_type), intent(in) :: statetime
integer                     :: timeindex

integer  :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeVarID
integer  :: xtype, ndims, nAtts, nTlen
integer  :: secs, days, ncid, i

character(len=NF90_MAX_NAME)          :: varname
integer, dimension(NF90_MAX_VAR_DIMS) :: dimids

timeindex = -1  ! assume bad things are going to happen

ncid = ncFileID%ncid

! Make sure we're looking at the most current version of the netCDF file.
! Get the length of the (unlimited) Time Dimension 
! If there is no length -- simply append a time to the dimension and return ...
! Else   get the existing times ["days since ..."] and convert to time_type 
!        if the statetime < earliest netcdf time ... we're in trouble
!        if the statetime does not match any netcdf time ... we're in trouble
!        if the statetime > last netcdf time ... append a time ... 

call nc_check(NF90_Sync(ncid), 'nc_get_tindex', 'sync '//trim(ncFileID%fname))    
call nc_check(NF90_Inquire(ncid, nDimensions, nVariables, nAttributes, unlimitedDimID), &
              'nc_get_tindex', 'inquire '//trim(ncFileID%fname))
call nc_check(NF90_Inq_Varid(ncid, "time", TimeVarID), &
              'nc_get_tindex', 'inq_varid time '//trim(ncFileID%fname))
call nc_check(NF90_Inquire_Variable(ncid, TimeVarID, varname, xtype, ndims, dimids, nAtts), &
              'nc_get_tindex', 'inquire_variable time '//trim(ncFileID%fname))
call nc_check(NF90_Inquire_Dimension(ncid, unlimitedDimID, varname, nTlen), &
              'nc_get_tindex', 'inquire_dimension unlimited '//trim(ncFileID%fname))

! Sanity check all cases first.

if ( ndims /= 1 ) then
   write(msgstring,*)'"time" expected to be rank-1' 
   call error_handler(E_WARN,'nc_get_tindex',msgstring,source,revision,revdate)
   timeindex = timeindex -   1
endif
if ( dimids(1) /= unlimitedDimID ) then
   write(msgstring,*)'"time" must be the unlimited dimension'
   call error_handler(E_WARN,'nc_get_tindex',msgstring,source,revision,revdate)
   timeindex = timeindex -  10
endif
if ( timeindex < -1 ) then
   write(msgstring,*)'trouble deep ... can go no farther. Stopping.'
   call error_handler(E_ERR,'nc_get_tindex',msgstring,source,revision,revdate)
endif

! convert statetime to time base of "days since ..."
call get_time(statetime, secs, days)

if (ncFileID%Ntimes < 1) then          ! First attempt at writing a state ...

   write(msgstring,*)'current unlimited  dimension length',nTlen, &
                     'for ncFileID ',trim(ncFileID%fname)
   call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)
   write(msgstring,*)'current time array dimension length',ncFileID%Ntimes
   call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)

   nTlen = nc_append_time(ncFileID, statetime)

   write(msgstring,*)'Initial time array dimension length',ncFileID%Ntimes
   call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)

endif



TimeLoop : do i = 1,ncFileId%Ntimes

   if ( statetime == ncFileID%times(i) ) then
      timeindex = i
      exit TimeLoop
   endif

enddo TimeLoop



if ( timeindex <= 0 ) then   ! There was no match. Either the model
                             ! time precedes the earliest file time - or - 
                             ! model time is somewhere in the middle  - or - 
                             ! model time needs to be appended.

   if (statetime < ncFileID%times(1) ) then

      call error_handler(E_MSG,'nc_get_tindex', &
              'Model time precedes earliest netCDF time.', source,revision,revdate)

      write(msgstring,*)'          model time (days, seconds) ',days,secs
      call error_handler(E_MSG,'nc_get_tindex',msgstring,source,revision,revdate)

      call get_time(ncFileID%times(1),secs,days)
      write(msgstring,*)'earliest netCDF time (days, seconds) ',days,secs
      call error_handler(E_MSG,'nc_get_tindex',msgstring,source,revision,revdate)

      call error_handler(E_ERR,'nc_get_tindex', &
              'Model time precedes earliest netCDF time.', source,revision,revdate)
      timeindex = -2

   else if ( statetime < ncFileID%times(ncFileID%Ntimes) ) then  

      ! It is somewhere in the middle without actually matching an existing time.
      ! This is very bad.

      write(msgstring,*)'model time does not match any netCDF time.'
      call error_handler(E_MSG,'nc_get_tindex',msgstring,source,revision,revdate)
      write(msgstring,*)'model time (days, seconds) is ',days,secs
      call error_handler(E_MSG,'nc_get_tindex',msgstring,source,revision,revdate)

      BadLoop : do i = 1,ncFileId%Ntimes   ! just find times to print before exiting

         if ( ncFileId%times(i) > statetime ) then
            call get_time(ncFileID%times(i-1),secs,days)
            write(msgstring,*)'preceding netCDF time (days, seconds) ',days,secs
            call error_handler(E_MSG,'nc_get_tindex',msgstring,source,revision,revdate)

            call get_time(ncFileID%times(i),secs,days)
            write(msgstring,*)'subsequent netCDF time (days, seconds) ',days,secs
            call error_handler(E_ERR,'nc_get_tindex',msgstring,source,revision,revdate)
            timeindex = -3
            exit BadLoop
         endif

      enddo BadLoop

   else ! we must need to append ... 

      timeindex = nc_append_time(ncFileID, statetime)

      write(msgstring,'(''appending model time (d,s) ('',i8,i5,'') as index '',i6, '' in ncFileID '',i10)') &
          days,secs,timeindex,ncid
      call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)

   endif
   
endif

end function nc_get_tindex



function nc_write_calendar_atts(ncFileID, TimeVarID) result(ierr)
!------------------------------------------------------------------------
!
! Need this to follow conventions for netCDF output files.

use typeSizes
use netcdf

implicit none

type(netcdf_file_type), intent(in) :: ncFileID
integer,                intent(in) :: TimeVarID
integer                            :: ierr

!integer  :: unlimitedDimID, length
integer  :: ncid
!character(len=NF90_MAX_NAME) :: varname

ierr = 0

ncid = ncFileID%ncid

!call check(NF90_Sync(ncid))    
!call check(NF90_Inquire_Dimension(ncid, unlimitedDimID, varname, length))
!
!if ( TimeVarID /= unlimitedDimID ) then
!   call error_handler(E_ERR,'nc_write_calendar_atts',&
!      'unlimited dimension is not time', source,revision,revdate)
!endif

call nc_check(nf90_put_att(ncid, TimeVarID, "long_name", "time"), &
              'nc_write_calendar_atts', 'put_att long_name '//trim(ncFileID%fname))
call nc_check(nf90_put_att(ncid, TimeVarID, "axis", "T"), &
              'nc_write_calendar_atts', 'put_att axis '//trim(ncFileID%fname))
call nc_check(nf90_put_att(ncid, TimeVarID, "cartesian_axis", "T"), &
              'nc_write_calendar_atts', 'put_att cartesian_axis '//trim(ncFileID%fname))

select case( get_calendar_type() )
case(THIRTY_DAY_MONTHS)
!  call get_date_thirty(time, year, month, day, hour, minute, second)
case(GREGORIAN)
   call nc_check(nf90_put_att(ncid, TimeVarID, "calendar", "gregorian" ), &
              'nc_write_calendar_atts', 'put_att calendar '//trim(ncFileID%fname))
   call nc_check(nf90_put_att(ncid, TimeVarID, "units", "days since 1601-01-01 00:00:00"), &
              'nc_write_calendar_atts', 'put_att units '//trim(ncFileID%fname))
case(JULIAN)
   call nc_check(nf90_put_att(ncid, TimeVarID, "calendar", "julian" ), &
              'nc_write_calendar_atts', 'put_att calendar '//trim(ncFileID%fname))
case(NOLEAP)
   call nc_check(nf90_put_att(ncid, TimeVarID, "calendar", "no_leap" ), &
              'nc_write_calendar_atts', 'put_att calendar '//trim(ncFileID%fname))
case default
   call nc_check(nf90_put_att(ncid, TimeVarID, "calendar", "no calendar" ), &
              'nc_write_calendar_atts', 'put_att calendar '//trim(ncFileID%fname))
   call nc_check(nf90_put_att(ncid, TimeVarID, "units", "days since 0000-00-00 00:00:00"), &
              'nc_write_calendar_atts', 'put_att units '//trim(ncFileID%fname))
end select

end function nc_write_calendar_atts


!
!===================================================================
! End of assim_model_mod
!===================================================================
!
end module assim_model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

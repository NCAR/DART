! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_def_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Contains the basic parts of a module for defining and evaluating observation
! definitions. Can evaluate identity observations as is. The DART preprocess
! program is used to add in extra observation kinds at the indicated spots in
! the code.

use        types_mod, only : r8, missing_i, missing_r8, RAD2DEG
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, logfileunit
use     location_mod, only : location_type, read_location, write_location, &
                             interactive_location, set_location_missing
use time_manager_mod, only : time_type, read_time, write_time, set_time, set_time_missing, &
                             interactive_time
use  assim_model_mod, only : get_state_meta_data, interpolate
use     obs_kind_mod, only : assimilate_this_obs_kind, evaluate_this_obs_kind, max_obs_kinds, &
                             get_obs_kind_name, map_def_index, get_kind_from_menu
use     obs_kind_mod, only : KIND_RAW_STATE_VARIABLE, KIND_U_WIND_COMPONENT, &
                             KIND_V_WIND_COMPONENT, KIND_SURFACE_PRESSURE, &
                             KIND_TEMPERATURE, KIND_SPECIFIC_HUMIDITY, KIND_PRESSURE, &
                             KIND_VERTICAL_VELOCITY, KIND_RAINWATER_MIXING_RATIO, &
                             KIND_DEW_POINT_TEMPERATURE, KIND_DENSITY, KIND_VELOCITY, &
                             KIND_1D_INTEGRAL, KIND_RADAR_REFLECTIVITY

! DART PREPROCESS USE FOR OBS_KIND_MOD INSERTED HERE

! Additional observation kinds may require use statements to support their operation.
! These additional use statements are added below this line by the DART preprocess program.
! Need to include use statements for whatever observation kind detail modules are being used

! DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE INSERTED HERE

implicit none
private

interface assignment(=)
   module procedure copy_obs_def
end interface

public :: init_obs_def, get_obs_def_key, get_obs_def_location, get_obs_kind, get_obs_def_time, &
          get_obs_def_error_variance, set_obs_def_location, set_obs_def_kind, set_obs_def_time, &
          set_obs_def_error_variance, set_obs_def_key, interactive_obs_def, write_obs_def, &
          read_obs_def, obs_def_type, get_expected_obs_from_def, destroy_obs_def, &
          copy_obs_def, assignment(=), get_obs_name

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type obs_def_type
! In revision, obs_kind module is responsible for taking care of identity obs kinds, too
   private
   type(location_type)   :: location      ! center of mass, so to speak
   integer               :: kind
   type(time_type)       :: time
   real(r8)              :: error_variance
   integer               :: key             ! Used by specialized observation types
end type obs_def_type

logical, save :: module_initialized = .false.

contains

!----------------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

! Note that there is no namelist for this module now that obs_kind has been revised

end subroutine initialize_module


!----------------------------------------------------------------------------

subroutine init_obs_def(obs_def, location, kind, time, error_variance)
! Need to add additional component arguments as optionals as needed

! Constructor for an obs_def

type(obs_def_type), intent(out) :: obs_def
type(location_type), intent(in) :: location
integer,             intent(in) :: kind
type(time_type),     intent(in) :: time
real(r8),            intent(in) :: error_variance

if ( .not. module_initialized ) call initialize_module

obs_def%location = location
obs_def%kind = kind
obs_def%time = time
obs_def%error_variance = error_variance
! No key assigned for standard observation defs
obs_def%key = -1

end subroutine init_obs_def

!---------------------------------------------------------------------

subroutine copy_obs_def(obs_def1, obs_def2)

! Copy function to be overloaded with '='

type(obs_def_type), intent(out) :: obs_def1
type(obs_def_type), intent(in) :: obs_def2

if ( .not. module_initialized ) call initialize_module

obs_def1%location       = obs_def2%location
obs_def1%kind           = obs_def2%kind
obs_def1%time           = obs_def2%time
obs_def1%error_variance = obs_def2%error_variance
obs_def1%key            = obs_def2%key
!deallocate(obs_def1%platform_qc)
!allocate(obs_def1%platform_qc(size(obs_def2%platform_qc))
! Should this be pointer assignment or regular
!obs_def1%platform_qc >= or == obs_def2%platform_qc
!obs_def1%aperture = obs_def2%aperture

end subroutine copy_obs_def

!----------------------------------------------------------------------------

function get_obs_def_key(obs_def)

type(obs_def_type), intent(in) :: obs_def
integer                        :: get_obs_def_key

if ( .not. module_initialized ) call initialize_module

get_obs_def_key = obs_def%key

end function get_obs_def_key

!----------------------------------------------------------------------------

function get_obs_def_error_variance(obs_def)

type(obs_def_type), intent(in) :: obs_def
real(r8)                       :: get_obs_def_error_variance

if ( .not. module_initialized ) call initialize_module

get_obs_def_error_variance = obs_def%error_variance

end function get_obs_def_error_variance

!----------------------------------------------------------------------------

function get_obs_def_location(obs_def)

! Returns observation location.

type(location_type)            :: get_obs_def_location
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_obs_def_location = obs_def%location

end function get_obs_def_location

!----------------------------------------------------------------------------

function get_obs_kind(obs_def)

! Returns observation kind

integer                        :: get_obs_kind
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_obs_kind = obs_def%kind

end function get_obs_kind

!----------------------------------------------------------------------------

function get_obs_def_time(obs_def)

! Returns observation time

type(time_type)                :: get_obs_def_time
type(obs_def_type), intent(in) :: obs_def

if ( .not. module_initialized ) call initialize_module

get_obs_def_time = obs_def%time

end function get_obs_def_time

!----------------------------------------------------------------------------

function get_obs_name(obs_kind_ind)

! Returns observation name

integer, intent(in) :: obs_kind_ind
character(len = 32) :: get_obs_name

if ( .not. module_initialized ) call initialize_module

get_obs_name = get_obs_kind_name(obs_kind_ind)

end function get_obs_name

!----------------------------------------------------------------------------

subroutine set_obs_def_location(obs_def, location)

! Sets the location of an obs_def

type(obs_def_type), intent(inout) :: obs_def
type(location_type),   intent(in) :: location

if ( .not. module_initialized ) call initialize_module

obs_def%location = location

end subroutine set_obs_def_location

!----------------------------------------------------------------------------

subroutine set_obs_def_error_variance(obs_def, error_variance)

! Sets the error variance of an obs_def

type(obs_def_type), intent(inout) :: obs_def
real(r8), intent(in) :: error_variance

if ( .not. module_initialized ) call initialize_module

obs_def%error_variance = error_variance

end subroutine set_obs_def_error_variance

!----------------------------------------------------------------------------

subroutine set_obs_def_key(obs_def, key)

! Sets the key of an obs_def

type(obs_def_type), intent(inout) :: obs_def
integer,            intent(in)    :: key

if ( .not. module_initialized ) call initialize_module

obs_def%key = key

end subroutine set_obs_def_key

!----------------------------------------------------------------------------

subroutine set_obs_def_kind(obs_def, kind)

! Sets the kind of an obs_def

type(obs_def_type), intent(inout) :: obs_def
integer,               intent(in) :: kind

if ( .not. module_initialized ) call initialize_module

obs_def%kind = kind

end subroutine set_obs_def_kind

!----------------------------------------------------------------------------

subroutine set_obs_def_time(obs_def, time)

! Sets the time of an obs_def

type(obs_def_type), intent(inout) :: obs_def
type(time_type), intent(in) :: time

if ( .not. module_initialized ) call initialize_module

obs_def%time = time

end subroutine set_obs_def_time

!----------------------------------------------------------------------------

subroutine get_expected_obs_from_def(key, obs_def, obs_kind_ind, state, obs_val, &
   istatus, assimilate_this_ob, evaluate_this_ob)

! Compute forward operator for a particular obs_def
integer, intent(in) :: key
type(obs_def_type), intent(in) :: obs_def
integer, intent(in) :: obs_kind_ind
real(r8), intent(in) :: state(:)
real(r8), intent(out) :: obs_val
integer, intent(out) :: istatus
logical, intent(out) :: assimilate_this_ob, evaluate_this_ob

type(location_type) :: location

! Load up the assimilate and evaluate status for this observation kind
assimilate_this_ob = assimilate_this_obs_kind(obs_kind_ind)
evaluate_this_ob = evaluate_this_obs_kind(obs_kind_ind)

! If not being assimilated or evaluated return with missing_r8 and istatus 0
if(assimilate_this_ob .or. evaluate_this_ob) then
   location = get_obs_def_location(obs_def)
   ! Compute the forward operator;
   select case(obs_kind_ind)

      ! CASE statements and algorithms for specific observation kinds are
      ! inserted here by the DART preprocess program.
      ! DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF INSERTED HERE

      ! If the observation kind is not available, it is an error. The DART preprocess
      ! program should provide code for all available kinds.
      case DEFAULT
         call error_handler(E_ERR, 'get_expected_obs_from_def', &
            'Attempt to evaluate or assimilate undefined obs_kind type.', source, revision, revdate)
   end select
else
   ! Not computing forward operator for this kind
   obs_val = missing_r8
   istatus = 0
endif

end subroutine get_expected_obs_from_def



  subroutine read_obs_def(ifile, obs_def, key, obs_val, fform)
!----------------------------------------------------------------------------
! subroutine read_obs_def(ifile, obs_def, key, obs_val, fform)
!
! ifile
! obs_def
! key
! obs_val    needed if you want to perform operations based on value 
! fform
!
! Reads an obs_def from file which is just an integer unit number in the
! current preliminary implementation.

integer,                    intent(in)    :: ifile
type(obs_def_type),         intent(inout) :: obs_def
integer,                    intent(in)    :: key
real(r8),                   intent(inout) :: obs_val
character(len=*), optional, intent(in)    :: fform

character(len=5)  :: header
character(len=32) :: fileformat
integer           :: o_index

if ( .not. module_initialized ) call initialize_module

! this default should probably be 'FORMATTED' to be
! consistent with the alternative of 'UNFORMATTED'. 
! (in the code below anything that is not a variant of
! unformatted is the default case so for now it does not matter.)
fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Begin by reading five character ascii header, then location, kind, error variance, index

! Need to add additional error checks on read
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      continue
   CASE DEFAULT
      read(ifile, 11) header
11    Format(a5)
      if(header /= 'obdef') then
   call error_handler(E_ERR,'read_obs_def', 'Expected header "obdef" in input file', &
                      source, revision, revdate)
      endif
END SELECT

! Read the location, kind, time and error variance
obs_def%location = read_location(ifile, fileformat)
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      ! Need to map to get index associated with this integer in obs_kind
      read(ifile) o_index
      obs_def%kind = map_def_index(o_index)
   CASE DEFAULT
      read(ifile, '(a5)' ) header
      if(header /= 'kind ') then
         call error_handler(E_ERR,'read_kind', 'Expected kind header "kind " in input file', &
                            source, revision, revdate)
      endif
      ! Need to map to get index associated with this integer in obs_kind
      read(ifile, *) o_index
      ! Negative value is identity obs, doesn't need mapped
      if(o_index < 0) then
         obs_def%kind = o_index
      else
         ! Positive value must use mapping to get to proper index in obs_kind
         obs_def%kind = map_def_index(o_index)
      endif
END SELECT

! This kind may have its own module that needs to read more
select case(obs_def%kind)
   ! More complicated kinds may require reading additional information from an observation
   ! sequence file. Case code to do this is inserted here by the DART preprocess program.

! DART PREPROCESS READ_OBS_DEF INSERTED HERE

! A negative value means identity observations, just move along
   case (:-1)
      continue

   case DEFAULT
      call error_handler(E_ERR, 'read_obs_def', &
         'Attempt to read for undefined obs_kind type.', source, revision, revdate)
end select

! Read the time for the observation
obs_def%time = read_time(ifile, fileformat)

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) obs_def%error_variance
   CASE DEFAULT
      read(ifile, *) obs_def%error_variance
END SELECT

end subroutine read_obs_def

!----------------------------------------------------------------------------

subroutine write_obs_def(ifile, obs_def, key, fform)

! Writes an obs_def to file.

integer,                    intent(in) :: ifile
type(obs_def_type),         intent(in) :: obs_def
integer,                    intent(in) :: key
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

! this default should probably be 'FORMATTED' to be
! consistent with the alternative of 'UNFORMATTED'. 
! (in the code below anything that is not a variant of
! unformatted is the default case so for now it does not matter.)
fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Write the 5 character identifier for verbose formatted output
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      continue
   CASE DEFAULT
      write(ifile, 11)
11    format('obdef')
END SELECT

! Write out the location, kind and error variance
call write_location(ifile, obs_def%location, fileformat)
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) obs_def%kind
   CASE DEFAULT
      write(ifile, '(''kind'')' )
      write(ifile, *) obs_def%kind
END SELECT

! This kind may have its own module that needs to write more
select case(obs_def%kind)
   ! More complicated kinds may require writing additional information from an observation
   ! sequence file. Case code to do this is inserted here by the DART preprocess program.

   ! DART PREPROCESS WRITE_OBS_DEF INSERTED HERE

   ! A negative value means identity observations, just move along
   case (:-1)
      continue

   case DEFAULT
      call error_handler(E_ERR, 'write_obs_def', &
         'Attempt to write for undefined obs_kind type.', source, revision, revdate)
end select

call write_time(ifile, obs_def%time, fileformat)

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) obs_def%error_variance
   CASE DEFAULT
      write(ifile, *) obs_def%error_variance
END SELECT

end subroutine write_obs_def


subroutine interactive_obs_def(obs_def, key)
!---------------------------------------------------------------------------
!
! Allows interactive creation of an observation

type(obs_def_type), intent(inout) :: obs_def
integer,               intent(in) :: key

if ( .not. module_initialized ) call initialize_module

! Get the observation kind WANT A STRING OPTION, TOO?
obs_def%kind = get_kind_from_menu()

! Input any special stuff for this kind
select case(obs_def%kind)
   ! More complicated kinds may require inputting additional information to define an
   ! observation. Case code to do this is inserted here by the DART preprocess program.

   ! DART PREPROCESS INTERACTIVE_OBS_DEF INSERTED HERE

   ! A negative value means identity observations, just move along
   case (:-1)
      continue
   case DEFAULT
      call error_handler(E_ERR, 'interactive_obs_def', &
         'Attempt to interactively create undefined obs_kind type.', source, revision, revdate)
end select

! If the kind is an identity observation, don't need to call location
! Get location from state meta_data
if(obs_def%kind < 0) then
! Get the location of this from model
   call get_state_meta_data(-1 * obs_def%kind, obs_def%location)
else! Get the location
   call interactive_location(obs_def%location)
endif

! Get the time
call interactive_time(obs_def%time)

write(*, *) 'Input error variance for this observation definition '
read(*, *) obs_def%error_variance

! TJH -- might want to do some sort of error checking (i.e. for positive values)

end subroutine interactive_obs_def

!----------------------------------------------------------------

subroutine destroy_obs_def(obs_def)
! TECHNICALLY NEED TO CALL DESTRUCTORS FOR ALL SUBCOMPONENTS, NO ALLOCATED STORAGE YET

type(obs_def_type), intent(inout) :: obs_def

if ( .not. module_initialized ) call initialize_module

call set_obs_def_location(       obs_def, set_location_missing() )
obs_def%kind = missing_i
call set_obs_def_time(           obs_def, set_time_missing() )
call set_obs_def_error_variance( obs_def, missing_r8)

end subroutine destroy_obs_def

!---------------------------------------------------------------------------

end module obs_def_mod

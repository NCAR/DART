! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_kind_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

! This module is used in conjunction with the preprocessor and a set of obs_kind
! modules to construct a final fortran 90 obs_kind_mod for use with a DART
! assimilation.

use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             logfileunit, find_namelist_in_file, check_namelist_read

implicit none
private

public :: max_obs_kinds, get_obs_kind_name, assimilate_this_obs_kind, &
          evaluate_this_obs_kind, get_obs_kind_var_type, get_obs_kind_index, &
          write_obs_kind, read_obs_kind, get_kind_from_menu, map_def_index

! Public access to the raw variable types
public :: KIND_RAW_STATE_VARIABLE, KIND_U_WIND_COMPONENT, &
          KIND_V_WIND_COMPONENT, KIND_SURFACE_PRESSURE, &
          KIND_TEMPERATURE, KIND_SPECIFIC_HUMIDITY, KIND_PRESSURE, &
          KIND_VERTICAL_VELOCITY, KIND_RAINWATER_MIXING_RATIO, &
          KIND_DEW_POINT_TEMPERATURE, KIND_DENSITY, KIND_VELOCITY, &
          KIND_1D_INTEGRAL, KIND_RADAR_REFLECTIVITY, &
          KIND_GRAUPEL_MIXING_RATIO, KIND_SNOW_MIXING_RATIO, &
          KIND_GPSRO, &
          KIND_CLOUD_LIQUID_WATER, KIND_CLOUD_ICE, &
          KIND_CONDENSATIONAL_HEATING, KIND_VAPOR_MIXING_RATIO, &
          KIND_ICE_NUMBER_CONCENTRATION, KIND_GEOPOTENTIAL_HEIGHT, &
          KIND_VORTEX_LON, KIND_VORTEX_LAT, KIND_VORTEX_PMIN, KIND_VORTEX_WMAX

! Public access to the observation types is provided here
! This is constructed by the preprocessor
! DART PREPROCESS PUBLIC INSERTED HERE

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

logical, save :: module_initialized = .false.

! A list of raw variable types that must be used by observation kinds but can also
! be used by models to define the type of their variables for get_state_meta_data
! and for forward operator evaluation by models. The value of ??? is reserved
! for raw state variable types so that low order models don't have to use this
! module.

integer, parameter :: KIND_RAW_STATE_VARIABLE          = 0, &
                      KIND_U_WIND_COMPONENT            = 1, &
                      KIND_V_WIND_COMPONENT            = 2, &
                      KIND_SURFACE_PRESSURE            = 3, &
                      KIND_TEMPERATURE                 = 4, &
                      KIND_SPECIFIC_HUMIDITY           = 5, &
                      KIND_PRESSURE                    = 6, &
                      KIND_VERTICAL_VELOCITY           = 7, &
                      KIND_RAINWATER_MIXING_RATIO      = 8, &
                      KIND_DEW_POINT_TEMPERATURE       = 9, &
                      KIND_DENSITY                     = 10, &
                      KIND_VELOCITY                    = 11, &
                      KIND_RADAR_REFLECTIVITY          = 12, &
                      KIND_1D_INTEGRAL                 = 13, &
                      KIND_GRAUPEL_MIXING_RATIO        = 14, &
                      KIND_SNOW_MIXING_RATIO           = 15, &
                      KIND_GPSRO                       = 16, &
                      KIND_CLOUD_LIQUID_WATER          = 17, &
                      KIND_CLOUD_ICE                   = 18, &
                      KIND_CONDENSATIONAL_HEATING      = 19, &
                      KIND_VAPOR_MIXING_RATIO          = 20, &
                      KIND_ICE_NUMBER_CONCENTRATION    = 21, &
                      KIND_GEOPOTENTIAL_HEIGHT         = 22

integer, parameter :: KIND_VORTEX_LON                  = 81, &
                      KIND_VORTEX_LAT                  = 82, &
                      KIND_VORTEX_PMIN                 = 83, &
                      KIND_VORTEX_WMAX                 = 84

! Index values associated with each observation kind string are defined
! by the preprocessor and inserted here. The total number of obs_kinds
! is also set in parameter max_obs_kinds.
! DART PREPROCESS INTEGER DECLARATION INSERTED HERE

integer :: num_def_obs_kinds = 0
integer :: num_kind_assimilate, num_kind_evaluate

! Map from values of kind in obs_def to the fixed values in the list above.
! Initially, these are undefined and have values -1.
! For the first index 1, the value is the index in the input obs_sequence file.
! The first index 2 is the value of the corresponding index in this kind module.
integer :: map(2, max_obs_kinds) = -1

! An observation kind type links together all the information required.
! An integer index that is also associated with the parameter above,
! A character string that has the same string as the parameter above,
! an integer that indicates what kind of variable type this is (for
! instance a radiosonde u wind component is a u wind component, but so
! is a 10 meter u wind component), and two logicals that indicate
! whether observations of this kind are to be assimilated, evaluated,
! or neither. Name lengths are limited to 32 characters by compiler
! restrictions on the length of parameter identifiers.
type obs_kind_type
   integer              :: index
   character(len = 32)  :: name
   integer              :: var_type
   logical              :: assimilate
   logical              :: evaluate
end type obs_kind_type

! An obs_kind_type is defined by the preprocessor to store the association
! between obs_kinds, associated integer identifiers, the underlying type,
! and whether observations of this type should be assimilate or evaluated.
type(obs_kind_type) :: obs_kind_info(max_obs_kinds)

! Namelist array to turn on any requested observation types
character(len = 129) :: assimilate_these_obs_types(max_obs_kinds) = 'null'
character(len = 129) :: evaluate_these_obs_types(max_obs_kinds) = 'null'

namelist /obs_kind_nml/ assimilate_these_obs_types, evaluate_these_obs_types

contains

!----------------------------------------------------------------------------

  subroutine initialize_module

integer :: iunit, io, i, j
character(len = 169) :: err_string

call register_module(source, revision, revdate)
module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_kind_nml", iunit)
read(iunit, nml = obs_kind_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_kind_nml")


! DART PREPROCESS OBS_KIND_INFO INSERTED HERE


call error_handler(E_MSG,'initialize_module','obs_kind_nml values are',' ',' ',' ')
write(logfileunit, *) 'assimilate_these_obs_types'
write(*, *)
write(*, *) '-------------- ASSIMILATE_THESE_OBS_TYPES --------------'
num_kind_assimilate = 0
do i = 1, max_obs_kinds
   if(assimilate_these_obs_types(i) == 'null') goto 22
   write(logfileunit, *) trim(assimilate_these_obs_types(i))
   write(*, *) trim(assimilate_these_obs_types(i))
   num_kind_assimilate = i
end do
22 write(logfileunit, *) 'evaluate_these_obs_types'
write(*, *) '-------------- EVALUATE_THESE_OBS_TYPES --------------'
num_kind_evaluate = 0
do i = 1, max_obs_kinds
   if(evaluate_these_obs_types(i) == 'null') goto 33
   write(logfileunit, *) trim(evaluate_these_obs_types(i))
   write(*, *) trim(evaluate_these_obs_types(i))
   num_kind_evaluate = i
end do
33 continue
write(*, *) '------------------------------------------------------'
write(*, *)

! Figure out which kinds are being used, look for errors
! Start by loading up kinds to assimilate
if (num_kind_assimilate > 0) then
   do i = 1, num_kind_assimilate
      ! Search for the matching string
      do j = 1, max_obs_kinds
         if(assimilate_these_obs_types(i) == obs_kind_info(j)%name) then
            obs_kind_info(j)%assimilate = .true.
            goto 44
         endif
      end do
      ! Falling off the end is an error
      write(err_string, *) trim(assimilate_these_obs_types(i)), &
         ' from obs_kind_nml is not a legal observation kind'
      call error_handler(E_ERR, 'initialize_module', err_string, source, revision, revdate)
      44 continue
   end do
endif

! Now look for kinds to evaluate
if (num_kind_evaluate > 0) then
   do i = 1, num_kind_evaluate
      ! Search for the matching string
      do j = 1, max_obs_kinds
         if(evaluate_these_obs_types(i) == obs_kind_info(j)%name) then
            obs_kind_info(j)%evaluate = .true.
            goto 55
         endif
      end do
      ! Falling off the end is an error
      write(err_string, *) trim(evaluate_these_obs_types(i)), &
         ' from obs_kind_nml is not a legal observation kind'
      call error_handler(E_ERR, 'initialize_module', err_string, source, revision, revdate)
      55 continue
   end do
endif

! Make it an error to ask to assimilate AND evaluate the same obs kind
do i = 1, max_obs_kinds
   if(obs_kind_info(i)%evaluate .and. obs_kind_info(i)%assimilate) then
      write(err_string, *) 'Illegal to evaluate and assimilate same kind ', trim(obs_kind_info(i)%name)
      call error_handler(E_ERR, 'initialize_module', err_string, source, revision, revdate)
   endif
end do

end subroutine initialize_module

!---------------------------------------------------------------------------

function map_def_index(obs_def_index)

! Argument is the index from the obs_def; needs to be mapped to the appropriate
! Integer storage index
integer, intent(in) :: obs_def_index
integer             :: map_def_index

character(len = 169) :: err_string
integer :: i

if ( .not. module_initialized ) call initialize_module

! Need to search through the first column of map to find this obs_def_index value
! Then return the index into table in this module from corresponding row in 
! second column.
do i = 1, max_obs_kinds
   if(map(1, i) == obs_def_index) then
      map_def_index = map(2, i)
      return
   endif
end do

! Error, didn't find this obs_def_index in the map
write(err_string, *) 'Could not find obs_def_index ', obs_def_index, ' in obs_kind map'
call error_handler(E_ERR, 'map_def_index', err_string, source, revision, revdate)

end function map_def_index

!----------------------------------------------------------------------------

function get_obs_kind_name(obs_kind_ind)

! Returns observation name

integer, intent(in) :: obs_kind_ind
character(len = 32) :: get_obs_kind_name

if ( .not. module_initialized ) call initialize_module

get_obs_kind_name = obs_kind_info(obs_kind_ind)%name

end function get_obs_kind_name

!----------------------------------------------------------------------------

function get_obs_kind_index(obs_kind_name)

! Returns the integer index corresponding to an observation type string name
! Returns a -1 if this string is not in list

character(len = 32), intent(in) ::obs_kind_name
integer                         :: get_obs_kind_index

integer :: i

if ( .not. module_initialized ) call initialize_module

do i = 1, max_obs_kinds
   if(trim(adjustl(obs_kind_name)) == trim(adjustl(obs_kind_info(i)%name))) then
      get_obs_kind_index = i
      return
   endif
end do

get_obs_kind_index = -1

end function get_obs_kind_index

!----------------------------------------------------------------------------

function assimilate_this_obs_kind(obs_kind_ind)

! Returns true if this obs_kind is being assimilated

logical             :: assimilate_this_obs_kind
integer, intent(in) :: obs_kind_ind

if ( .not. module_initialized ) call initialize_module

assimilate_this_obs_kind = obs_kind_info(obs_kind_ind)%assimilate

end function assimilate_this_obs_kind

!----------------------------------------------------------------------------

function evaluate_this_obs_kind(obs_kind_ind)

! Returns true if this obs_kind is being assimilated

logical             :: evaluate_this_obs_kind
integer, intent(in) :: obs_kind_ind

if ( .not. module_initialized ) call initialize_module

evaluate_this_obs_kind = obs_kind_info(obs_kind_ind)%evaluate

end function evaluate_this_obs_kind

!----------------------------------------------------------------------------

function get_obs_kind_var_type(obs_kind_ind)

! Returns underlying variable type of this observation

integer, intent(in) :: obs_kind_ind
integer             :: get_obs_kind_var_type

if ( .not. module_initialized ) call initialize_module

get_obs_kind_var_type = obs_kind_info(obs_kind_ind)%var_type

end function get_obs_kind_var_type

!----------------------------------------------------------------------------

subroutine write_obs_kind(ifile, fform)

! Writes out the observation kind strings and corresponding integer
! indices as a header for an obs_sequence file.

integer,                    intent(in) :: ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat
integer :: i

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Write the 5 character identifier for verbose formatted output
SELECT CASE (fileformat)
   ! This header needs to be written for formatted OR unformatted
   ! If it's not present, it means to use the default old definitions
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) 'obs_kind_definitions'
   CASE DEFAULT
      write(ifile, 11)
11    format('obs_kind_definitions')
END SELECT

! Loop through the list to write out the integer indices and strings
! For all the defined observation types
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      ! Write the number of defined kinds, then the list
      write(ifile) max_obs_kinds
      do i = 1, max_obs_kinds
         write(ifile) obs_kind_info(i)%index, obs_kind_info(i)%name
      end do
   CASE DEFAULT
      write(ifile, *) max_obs_kinds
      do i = 1, max_obs_kinds
         write(ifile, *) obs_kind_info(i)%index, obs_kind_info(i)%name
      end do
END SELECT

end subroutine write_obs_kind

!----------------------------------------------------------------------------


subroutine read_obs_kind(ifile, pre_I_format, fform)

! Reads the observation kind strings and corresponding integer
! indices as a header for an obs_sequence file. If this isn't
! present, need to revert to default mapping for backwards
! compatibility.

integer,                    intent(in) :: ifile
logical,                    intent(in) :: pre_I_format
character(len=*), intent(in), optional :: fform

character(len=20)  :: header
character(len=32)  :: fileformat, o_name
character(len=129) :: msg_string
integer            :: i, num_def_kinds, o_index, list_index

if ( .not. module_initialized ) call initialize_module

! If this is old format, there's no obs_kind header to read
! Still need to initialize input kind map to use the order in
! the obs_kind file. It's users responsibility to make sure
! that this order is consistent with what the obs_sequence
! file thinks.
if(pre_I_format) then
   do i = 1, max_obs_kinds
      map(1, i) = i; map(2, i) = i
   end do
   return
endif

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Read the 5 character identifier for verbose formatted output
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      ! Need to look for header string
      read(ifile) header
      if(header /= 'obs_kind_definitions') then
         call error_handler(E_ERR, 'read_obs_kind', 'Didnt find obs_kind_definitions string', &
            source, revision, revdate)
      endif
   CASE DEFAULT
      read(ifile, 11) header
11    format(a20)
      if(header /= 'obs_kind_definitions') then
         call error_handler(E_ERR, 'read_obs_kind', 'Didnt find obs_kind_definitions string', &
            source, revision, revdate)
      endif
END SELECT

! Loop through the list to read the integer indices and strings
! For all the defined observation types
! Set up the map from kinds in the obs_sequence file to those
! in the data structure in this module.
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) num_def_kinds
      do i = 1, num_def_kinds
         read(ifile) o_index, o_name
         ! What is the integer associated with this o_name in this module?
         list_index = get_obs_kind_index(o_name)
         ! Check for error
         if(list_index == -1) then
            write(msg_string, *) 'didnt find observation kind ', o_name, ' in obs_kind_mod list'
            call error_handler(E_ERR, 'read_obs_kind', msg_string, &
               source, revision, revdate)
         endif
         map(1, i) = o_index
         map(2, i) = list_index
      end do
   CASE DEFAULT
      read(ifile, *) num_def_kinds
      do i = 1, num_def_kinds
         read(ifile, *) o_index, o_name
         ! What is the integer associated with this o_name in this module?
         list_index = get_obs_kind_index(o_name)
         ! Check for error
         if(list_index == -1) then
            write(msg_string, *) 'didnt find observation kind ', o_name, ' in obs_kind_mod list'
            call error_handler(E_ERR, 'read_obs_kind', msg_string, &
               source, revision, revdate)
         endif
         map(1, i) = o_index
         map(2, i) = list_index
      end do
END SELECT

end subroutine read_obs_kind

!----------------------------------------------------------------------------

function get_kind_from_menu()

integer :: get_kind_from_menu

integer :: i, ierr
character(len=32) :: in

if ( .not. module_initialized ) call initialize_module

! Should only do kinds that have been selected by preprocessor, so those
! are ones that are being evaluated or assimilated.
21 write(*, *) '     ', 'Input -1 * state variable index for identity observations'
write(*, *) '     ', 'OR input the name of the observation kind from table below:'
write(*, *) '     ', 'OR input the integer index, BUT see documentation...'
do i = 1, max_obs_kinds
   if(assimilate_this_obs_kind(i) .or. evaluate_this_obs_kind(i)) &
      write(*, *) '     ',  obs_kind_info(i)%index, trim(obs_kind_info(i)%name)
end do

! Read the input as a string, convert to integers as appropriate 
read(*, 11) in
11 format(A)

! If string is a positive or negative number, convert it to integer
read(in, *, IOSTAT = ierr) get_kind_from_menu
if(ierr /= 0) then
   get_kind_from_menu = get_obs_kind_index(in)
   ! If string isn't found, need to reprompt
   if(get_kind_from_menu == -1) then
      write(*, *) trim(in) // ' is not a supported kind: Please try again.'
      goto 21
   endif
else
   ! Make sure that number entered isn't 0 or too larg
   if(get_kind_from_menu == 0 .or. get_kind_from_menu > max_obs_kinds) then
      write(*, *) get_kind_from_menu, 'is not a legal entry: Please try again.'
      goto 21
   endif
endif

! Make sure 

end function get_kind_from_menu

!----------------------------------------------------------------------------

end module obs_kind_mod
